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

#include <fstream>
#include <iostream>

using namespace dealii;


/* my boundary values: pi/2 between vortices, 0 on other sides */
class BoundaryValues : public Function<2>{
    double R; // radius between vortices
    double C; // vortex winding

  public:
    BoundaryValues (const double R_, const double C_=0.5) :
      Function<2>(), R(R_), C(C_){}

    double value (const Point<2> &p, const unsigned int ) const{
      return p(1)<R && p(0)<R ? C*M_PI:0.0; }
};



template <int dim>
class Polar {
public:
  Polar ();
  ~Polar ();
  void run ();

  double Lx, Ly, R;

private:
  void setup_system (const bool step0);
  void assemble_system (const double v);
  void solve ();
  void refine_grid ();
  void output_results (const unsigned int cycle) const;

  Triangulation<dim>   triangulation;
  DoFHandler<dim>      dof_handler;
  FE_Q<dim>            fe;
  ConstraintMatrix     constraints;
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;
  Vector<double>       solution;
  Vector<double>       system_rhs;
};

template <int dim>
Polar<dim>::Polar () : dof_handler (triangulation), fe (2) {
  R=3;
  Lx=5.0;
  Ly=5.0;
}

template <int dim> Polar<dim>::~Polar () { dof_handler.clear (); }


template <int dim>
void Polar<dim>::setup_system(const bool step0){
  if (step0) {
    dof_handler.distribute_dofs (fe);
    solution.reinit (dof_handler.n_dofs());
    constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 0,
                                              BoundaryValues(R),
                                              constraints);
    constraints.close ();
  }
  system_rhs.reinit (dof_handler.n_dofs());
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, dsp);
  constraints.condense(dsp);
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);
}

template <int dim>
void Polar<dim>::assemble_system (const double v)
{
  const QGauss<dim>  quadrature_formula(3);

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |  update_gradients |
                           update_quadrature_points  |  update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  std::vector<double> cell_sol(n);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();

  for (; cell!=endc; ++cell) {
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit (cell);
      fe_values.get_function_values(solution, cell_sol);

      for (unsigned int q=0; q<n; ++q)
        for (unsigned int i=0; i<dofs_per_cell; ++i){
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_grad(i,q) *
                                   fe_values.shape_grad(j,q) *
                                   fe_values.JxW(q));
            cell_rhs(i) += -(fe_values.shape_value(i,q) *
                            0.5*v*sin(2.0*cell_sol[q]) *
                            fe_values.JxW(q));

          }

      // Finally, transfer the contributions from @p cell_matrix and
      // @p cell_rhs into the global objects.
      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global (cell_matrix,
                                              cell_rhs,
                                              local_dof_indices,
                                              system_matrix,
                                              system_rhs);
    }
}

template <int dim>
void Polar<dim>::solve ()
{
  SolverControl      solver_control (1000, 1e-12);
  SolverCG<>         solver (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  solver.solve (system_matrix, solution, system_rhs,
                preconditioner);

  constraints.distribute (solution);
}

template <int dim>
void Polar<dim>::refine_grid (){
  Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

  // estimate error per cell
  KellyErrorEstimator<dim>::estimate (dof_handler,
                                      QGauss<dim-1>(3),
                                      typename FunctionMap<dim>::type(),
                                      solution,
                                      estimated_error_per_cell);

  // mark cells for refinement (30% worse to refine, 3%best to coarse)
  GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                  estimated_error_per_cell, 0.3, 0.03);
  // mark additional cells to be refined
  triangulation.prepare_coarsening_and_refinement();

  // transfer the solution to the new mesh
  SolutionTransfer<dim> solution_transfer(dof_handler);
  solution_transfer.prepare_for_coarsening_and_refinement(solution);

  // do actual rifinement
  triangulation.execute_coarsening_and_refinement();

  dof_handler.distribute_dofs(fe);

  Vector<double> tmp(dof_handler.n_dofs());
  solution_transfer.interpolate(solution, tmp);
  solution = tmp;

  constraints.clear();
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           constraints);
  VectorTools::interpolate_boundary_values (dof_handler, 0,
                                            BoundaryValues(R),
                                            constraints);
  constraints.close();
  constraints.distribute (solution);

  setup_system (false);
}


template <int dim>
void Polar<dim>::output_results (const unsigned int cycle) const {
  Assert (cycle < 10, ExcNotImplemented());

  std::string fname1 = "grid-";
  fname1 += ('0' + cycle);
  fname1 += ".eps";

  std::ofstream output1 (fname1.c_str());

  GridOut grid_out;
  grid_out.write_eps (triangulation, output1);



  DataOutBase::EpsFlags eps_flags;
  eps_flags.z_scaling = 1;
  eps_flags.azimut_angle = 40;
  eps_flags.turn_angle   = 130;

  DataOut<dim> data_out;
  data_out.set_flags (eps_flags);

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();

  std::string fname2 = "out-";
  fname2 += ('0' + cycle);
  fname2 += ".eps";
  std::ofstream output2 (fname2.c_str());
  data_out.write_eps (output2);

}



template <int dim>
void Polar<dim>::run (){

  // build the initial grid
  GridGenerator::hyper_rectangle (triangulation,
    Point<2>(0,0), Point<2>(Lx,Ly));
  triangulation.refine_global(2);

  // set cell indicators X=0
  typename Triangulation<dim>::cell_iterator
    cell = triangulation.begin(),
    endc = triangulation.end();
  for (; cell!=endc; ++cell){
    for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; i++){
      if (std::fabs(cell->face(i)->center()(0)) < 1e-12)
        cell->face(i)->set_boundary_id(1);
    }
  }

  setup_system(true);

  for (unsigned int i=0; i<10; ++i){
      std::cout << "i " << i << ':' << std::endl;

      if (i>0 && i<7) refine_grid();

      double v = i<6? i/6.0: 1.0;
      assemble_system(v);

      std::cout << "   Number of active cells:       "
                << triangulation.n_active_cells()
                << std::endl;

      std::cout << "   Number of degrees of freedom: "
                << dof_handler.n_dofs()
                << std::endl;

      solve();
      output_results(i);
  }
}


int main (){
  try {
    deallog.depth_console (0);
    Polar<2> laplace_problem_2d;
    laplace_problem_2d.run ();
  }
  catch (std::exception &exc) {
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl;
    return 1;
  }
  return 0;
}
