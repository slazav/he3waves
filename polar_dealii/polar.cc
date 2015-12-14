
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>

#include <deal.II/lac/eigen.h>

#include <deal.II/lac/arpack_solver.h>

#include <fstream>
#include <sstream>
#include <iostream>

using namespace dealii;

template <int dim>
class Polar {
  public:
    Polar();
    ~Polar();
    void run();

    double Lx, Ly, R;

  private:
    // run the calculation
    // stage = 0 - first step
    // stage = 1 - steps of texture calculations
    // stage = 2 - texture calculations w/o grid refinment
    // stage = 3 - eigenproblem step
    void run_step(const int stage, const int step);

    Triangulation<dim>   triangulation;
    DoFHandler<dim>      dof_handler;
    FE_Q<dim>            fe;
    ConstraintMatrix     constraints;
    SparsityPattern      sparsity_pattern;

    Vector<double>       texture;
    Vector<double>       wave;
    double               en;
};

template <int dim>
Polar<dim>::Polar() : dof_handler(triangulation), fe(2) {
  R=14;  Lx=20.0;  Ly=20.0;
}

template <int dim> Polar<dim>::~Polar() { dof_handler.clear(); }

template <int dim>
void Polar<dim>::run(){
  for (unsigned int i=0; i<15; ++i){
    int stage = i<8? 1:2;
    if (i==0) stage=0;
    run_step(stage, i);
  }
  run_step(3, 10);
}


template <int dim>
void Polar<dim>::run_step(const int stage, const int step){

  // build the initial grid
  if (stage==0){ // first step
    GridGenerator::hyper_rectangle(triangulation,
      Point<2>(0,0), Point<2>(Lx,Ly));
    triangulation.refine_global(4);

    // enumerate degrees of freedom
    dof_handler.distribute_dofs(fe);

    // init the texture vector
    texture.reinit(dof_handler.n_dofs());
  }

  // refine the grid according with the error in the texture
  if (stage==1) { // refinement steps
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    // estimate error per cell
    KellyErrorEstimator<dim>::estimate(dof_handler,
                                        QGauss<dim-1>(3),
                                        typename FunctionMap<dim>::type(),
                                        texture,
                                        estimated_error_per_cell);

    // mark cells for refinement(30% worse to refine, 3%best to coarse)
    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                    estimated_error_per_cell, 0.1, 0.5	);
    // mark additional cells to be refined
    triangulation.prepare_coarsening_and_refinement();

    // prepare transfer of the solution to the new mesh
    SolutionTransfer<dim> solution_transfer(dof_handler);
    solution_transfer.prepare_for_coarsening_and_refinement(texture);

    // do actual rifinement
    triangulation.execute_coarsening_and_refinement();

    // enumerate degrees of freedom
    dof_handler.distribute_dofs(fe);

    // do transfer of the solution to the new mesh
    Vector<double> tmp(dof_handler.n_dofs());
    solution_transfer.interpolate(texture, tmp);
    texture = tmp;
  }

  // reset boundary IDs
  if (stage==0 || stage==1) {
    // set cell indicators X=0
    typename Triangulation<dim>::cell_iterator
      cell = triangulation.begin(),
      endc = triangulation.end();
    for (; cell!=endc; ++cell){
      for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; i++){
        // x=0 ID=1
        if (std::fabs(cell->face(i)->center()(0)) < 1e-12)
          cell->face(i)->set_boundary_id(1);
        else if (std::fabs(cell->face(i)->center()(0)) < R &&
                 std::fabs(cell->face(i)->center()(1)) < 1e-12)
          cell->face(i)->set_boundary_id(2);
      }
    }
  }


  // init the wave vector at all stages
  wave.reinit(dof_handler.n_dofs());

  // set constraints
  if (stage==0 || stage==1){
    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler,
                                             constraints);
    VectorTools::interpolate_boundary_values(dof_handler, 0,
                                              ZeroFunction<dim>(),
                                              constraints);
    VectorTools::interpolate_boundary_values(dof_handler, 2,
                                              ConstantFunction<dim>(0.5*M_PI),
                                              constraints);
    constraints.close();
  }

  // different boundary conditions for the stage==3
  if (stage==3){
    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler,
                                             constraints);
    VectorTools::interpolate_boundary_values(dof_handler, 0,
                                              ZeroFunction<dim>(),
                                              constraints);
    VectorTools::interpolate_boundary_values(dof_handler, 1,
                                              ZeroFunction<dim>(),
                                              constraints);
    VectorTools::interpolate_boundary_values(dof_handler, 2,
                                              ZeroFunction<dim>(),
                                              constraints);
    constraints.close();
  }

  if (stage==1) constraints.distribute(texture);

  // assemble the matrices
  SparseMatrix<double> A;
  Vector<double>       B;

  B.reinit(dof_handler.n_dofs());
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
  //constraints.condense(dsp);
  sparsity_pattern.copy_from(dsp);
  A.reinit(sparsity_pattern);

  const QGauss<dim>  quadrature_formula(3);
  FEValues<dim> fe_values(fe, quadrature_formula,
                          update_values    |  update_gradients |
                          update_quadrature_points  |  update_JxW_values);

  const unsigned int   nd = fe.dofs_per_cell;
  const unsigned int   nq = quadrature_formula.size();

  FullMatrix<double>   Acell(nd,nd);
  Vector<double>       Bcell(nd);
  std::vector<types::global_dof_index> local_dof_indices(nd);

  // texture calculation
  if (stage==0 || stage==1 || stage==2){

    std::vector<double> texture_cell(nq);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell) {
      Acell = 0;
      Bcell = 0;

      fe_values.reinit(cell);
      fe_values.get_function_values(texture, texture_cell);

      // here I solve: nabla^2 Un+1 - Un+1 = sin(2Un)/2 - Un
      for (unsigned int q=0; q<nq; ++q){
        for (unsigned int i=0; i<nd; ++i){
          for (unsigned int j=0; j<nd; ++j){
            Acell(i,j) -=(fe_values.shape_grad(i,q) *
                          fe_values.shape_grad(j,q) *
                          fe_values.JxW(q));
            Acell(i,j) -=(fe_values.shape_value(i,q) *
                          fe_values.shape_value(j,q) *
                          fe_values.JxW(q));
          }
          Bcell(i) += (fe_values.shape_value(i,q) *
                       0.5*sin(2.0*texture_cell[q]) *
                       fe_values.JxW(q));
          Bcell(i) -= (fe_values.shape_value(i,q) *
                       texture_cell[q] *
                       fe_values.JxW(q));
        }
      }
      // Finally, transfer the contributions from Acell and
      // Bcell into the global objects.
      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(
         Acell, Bcell, local_dof_indices, A, B);
    }
  }

  // wave calculation
  if (stage==3){

    std::vector<double> texture_cell(nq);
    std::vector<Tensor<1,dim> > tgrad_cell(nq);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell) {
      Acell = 0;
      Bcell = 0;

      fe_values.reinit(cell);
      fe_values.get_function_values(texture, texture_cell);
      fe_values.get_function_gradients(texture, tgrad_cell);

      for (unsigned int q=0; q<nq; ++q){
        for (unsigned int i=0; i<nd; ++i){
          for (unsigned int j=0; j<nd; ++j){
            Acell(i,j) +=(fe_values.shape_grad(i,q) *
                          fe_values.shape_grad(j,q) *
                          fe_values.JxW(q));

//            Acell(i,j) -=(fe_values.shape_value(i,q) *
//                          fe_values.shape_value(j,q) *
//                          pow(sin(texture_cell[q]),2)*
//                          fe_values.JxW(q));
//
//            Acell(i,j) +=(fe_values.shape_value(i,q) *
//                	  fe_values.shape_value(j,q) *
//                          tgrad_cell[q].norm_square() *
//                          fe_values.JxW(q));

          }
          Bcell(i) += (fe_values.shape_value(i,q) *
                       0.1 *
                       fe_values.JxW(q));
        }
      }
      // Finally, transfer the contributions from Acell and
      // Bcell into the global objects.
      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(
         Acell, Bcell, local_dof_indices, A, B);
    }
  }


  // solve the system

  if (stage==0 || stage==1 || stage==2){
    SolverControl      solver_control(1000, 1e-12);
    SolverCG<>         solver(solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(A, 1.2);
    solver.solve(A, texture, B, preconditioner);
    constraints.distribute(texture);
  }

  if (stage==3){
    SolverControl      solver_control(10000, 1e-6);
//    SolverCG<>         solver1(solver_control);

//    PreconditionSSOR<> preconditioner;
//    preconditioner.initialize(A, 1.2);
//    solver1.solve(A, wave, B, preconditioner);
//    constraints.distribute(wave);

    for (unsigned int i = 0; i<wave.size(); i++) wave[i]=1;

    en=-10.01;
    GrowingVectorMemory<> mem;
    EigenInverse<>::AdditionalData data;
    data.relaxation = 0.01;
    data.use_residual = false;
    EigenInverse<> solver2(solver_control,mem, data);
    solver2.solve(en, A, wave);
    constraints.distribute(wave);
    std::cout << "En: " << en << "\n";
  }

  // write grid
  if (stage==0 || stage==1){
    std::ostringstream fname1;
    fname1 << "grid-" << std::setfill('0') << std::setw(2) << step << ".eps";
    std::ofstream out1(fname1.str().c_str());
    GridOut grid_out;
    grid_out.write_eps(triangulation,out1);
  }

  DataOutBase::EpsFlags eps_flags;
  eps_flags.z_scaling = 1;
  eps_flags.azimut_angle = 40;
  eps_flags.turn_angle   = 130;

  // write texture
  if (stage==0 || stage==1 || stage==2){

    DataOut<dim> data_out;
    data_out.set_flags(eps_flags);

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(texture, "texture");
    data_out.build_patches();

    std::ostringstream fname2;
    fname2 << "text-" << std::setfill('0') << std::setw(2) << step << ".eps";
    std::ofstream out2(fname2.str().c_str());
    data_out.write_eps(out2);
  }

  // write wave
  if (stage==3){
    DataOut<dim> data_out;
    data_out.set_flags(eps_flags);

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(wave, "wave");
    data_out.build_patches();

    std::ostringstream fname2;
    fname2 << "wave-" << std::setfill('0') << std::setw(2) << step << ".eps";
    std::ofstream out2(fname2.str().c_str());
    data_out.write_eps(out2);
  }


  double It=0, Iw=0;
  if (stage==0 || stage==1 || stage==2){
    // estimate error (wave is used for the difference)
    Vector<float> difference_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference (dof_handler,
      texture, ZeroFunction<dim>(), difference_per_cell,
      QGauss<dim>(3),VectorTools::L2_norm);
    It = difference_per_cell.l2_norm();
  }

  if (stage==3){
    // estimate error (wave is used for the difference)
    Vector<float> difference_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference (dof_handler,
      wave, ZeroFunction<dim>(), difference_per_cell,
      QGauss<dim>(3),VectorTools::L2_norm);
    Iw = difference_per_cell.l2_norm();
  }




  // print some information about the cell and dof
  std::cout << "  step: " << std::setw(3) << step
            << "  act.cells: " << std::setw(4) << triangulation.n_active_cells()
            << "  DOF: " << std::setw(4) << dof_handler.n_dofs()
            << "  T " << It
            << "  W " << Iw
            << std::endl;


}



int main(){
  try {
    deallog.depth_console(0);
    Polar<2> laplace_problem_2d;
    laplace_problem_2d.run();
  }
  catch(std::exception &exc) {
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl;
    return 1;
  }
  catch(...){
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl;
    return 1;
  }
  return 0;
}
