#include "shifted_matrix.h"
#include "eigen_inverse_m.h"

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
#include <sstream>
#include <iostream>

#define dim 2

using namespace dealii;

/*************************************************************************/
class TestWaveSolver {
  public:

  TestWaveSolver(): dofs(triang), fe(2) {}
  ~TestWaveSolver() {  dofs.clear();  triang.clear(); }

  // Make the initial grid and initialize the wave to zero.
  void make_initial_grid();

  // Refine the grid according with the wave gradients
  // and redistribute the wave in the new grid.
  double refine_grid(const double ref=0.3, const double crs=0.03);

  // save the grid and a data vector in a eps file
  void save_grid(const char* fname);
  void save_wave(const char* fname);

  // Do the wave calculations.
  void do_calc();

  // some tests
  double check_wave();

  // the result
  Vector<double>       wave;
  double               en;

  Triangulation<dim>     triang;
  DoFHandler<dim>        dofs;
  FE_Q<dim>              fe;
};


// Make the initial grid and initialize the the wave to zero.
void
TestWaveSolver::make_initial_grid(){
  // rectangular grid 1x1
  GridGenerator::hyper_rectangle(triang,
    Point<dim>(0,0), Point<dim>(1,1));

  // initial refinment
  triang.refine_global(4);

  // enumerate degrees of freedom
  dofs.distribute_dofs(fe);

  // initialize initial wave
  wave.reinit(dofs.n_dofs());
}


// Refine the grid according with the wave gradients
// and redistribute the wave.
double
TestWaveSolver::refine_grid(const double ref, const double crs){
  Vector<float> err(triang.n_active_cells());

  // estimate error per cell
  KellyErrorEstimator<dim>::estimate(dofs, QGauss<dim-1>(3),
     typename FunctionMap<dim>::type(), wave, err);

  // mark cells for refinement( <ref> worse to refine, <crs> best to coarse)
  GridRefinement::refine_and_coarsen_fixed_number(triang, err, ref, crs);

  // mark additional cells to be refined
  triang.prepare_coarsening_and_refinement();

  // prepare transfer of the solution to the new mesh
  SolutionTransfer<dim> solution_transfer(dofs);
  solution_transfer.prepare_for_coarsening_and_refinement(wave);

  // do actual rifinement
  triang.execute_coarsening_and_refinement();

  // enumerate degrees of freedom
  dofs.distribute_dofs(fe);

  // do transfer of the solution to the new mesh
  Vector<double> tmp(dofs.n_dofs());
  solution_transfer.interpolate(wave, tmp);
  wave = tmp;

  return err.mean_value();
}

/*************************************************************************/
// save the grid into an eps file
void
TestWaveSolver::save_grid(const char* fname){
  GridOut grid_out;
  std::ofstream out(fname);
  grid_out.write_eps(triang,out);
}

/*************************************************************************/
// save the wave into an eps file
void
TestWaveSolver::save_wave(const char* fname){
  DataOutBase::EpsFlags eps_flags;

  eps_flags.z_scaling = 1./wave.linfty_norm();
  eps_flags.azimut_angle = 40;
  eps_flags.turn_angle   = 130;

  DataOut<dim> data_out;
  data_out.set_flags(eps_flags);

  data_out.attach_dof_handler(dofs);
  data_out.add_data_vector(wave, "wave");
  data_out.build_patches();

  std::ofstream out(fname);
  data_out.write_eps(out);
}


/*************************************************************************/
// Do the wave calculations.
// repeat=true shows that wave didn't change after the last calculation
void
TestWaveSolver::do_calc(){

  ConstraintMatrix     constraints;
  SparsityPattern      sparsity;
  SparseMatrix<double> A,M;

  // set text constraints (hanging nodes + zero BC)
  constraints.clear();
  DoFTools::make_hanging_node_constraints(dofs, constraints);
  VectorTools::interpolate_boundary_values(
      dofs, 0, ZeroFunction<dim>(), constraints);
  constraints.close();

  // rebuild sparsity pattern for the matrix
  DynamicSparsityPattern dsp(dofs.n_dofs());
  DoFTools::make_sparsity_pattern(dofs, dsp);
  constraints.condense(dsp);
  sparsity.copy_from(dsp);

  // build matrices
  {
    A.reinit(sparsity);
    M.reinit(sparsity);

    const QGauss<dim>  quadrature_formula(3);
    FEValues<dim> fe_values(fe, quadrature_formula,
                update_values | update_gradients |
                update_quadrature_points | update_JxW_values);

    const unsigned int   nd = fe.dofs_per_cell;
    const unsigned int   nq = quadrature_formula.size();

    // cell matrices
    FullMatrix<double>   Acell(nd,nd);
    FullMatrix<double>   Mcell(nd,nd);
    std::vector<types::global_dof_index> local_dof_indices(nd);

    typename DoFHandler<dim>::active_cell_iterator cell;
    for (cell=dofs.begin_active(); cell!=dofs.end(); ++cell) {
      Acell = 0;
      Mcell = 0;

      fe_values.reinit(cell);

      for (unsigned int q=0; q<nq; ++q){
        for (unsigned int i=0; i<nd; ++i){
          for (unsigned int j=0; j<nd; ++j){
            Acell(i,j) +=(fe_values.shape_grad(i,q) *
                          fe_values.shape_grad(j,q) *
                          fe_values.JxW(q));
            Mcell(i,j) +=(fe_values.shape_value(i,q) *
                          fe_values.shape_value(j,q) *
                          fe_values.JxW(q));
          }
        }
      }
      // Transfer the contributions from Acell and
      // Bcell into the global objects.
      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(
         Acell, local_dof_indices, A);
      constraints.distribute_local_to_global(
         Mcell, local_dof_indices, M);
    }
  }


  // solve the system
  {
    // initial guess:
    en=0;
    for (unsigned int i = 0; i<wave.size(); i++) wave[i]=1.0;

    SolverControl solver_control(1000, 1e-8);
    GrowingVectorMemory<> mem;
    EigenInverseM<>::AdditionalData data;
    data.relaxation = 0.1;
    data.use_residual = false;
    EigenInverseM<> solver(solver_control,mem,data);
    solver.solve(en, A, M, wave);

  }

  constraints.distribute(wave);

}

/*************************************************************************/
// Some tests.
double
TestWaveSolver::check_wave(){

  const QGauss<dim>  quadrature_formula(3);
  FEValues<dim> fe_values(fe, quadrature_formula,
              update_values | update_gradients | update_hessians |
              update_quadrature_points | update_JxW_values);

  const unsigned int   nq = quadrature_formula.size();

  // wave in the cell's quadrature points
  // same for laplacian
  std::vector<double> cell_val(nq);
  std::vector<double> cell_lap(nq);
  std::vector<Tensor<1,dim> > cell_grad(nq);

  double I1=0; // integral
  typename DoFHandler<dim>::active_cell_iterator cell;
  for (cell= dofs.begin_active(); cell!=dofs.end(); ++cell) {
    fe_values.reinit(cell);
    fe_values.get_function_values(wave, cell_val);
    fe_values.get_function_laplacians(wave, cell_lap);
    fe_values.get_function_gradients(wave, cell_grad);

    for (unsigned int q=0; q<nq; ++q){

      std::cerr << "> " << cell_val[q]
                << " "  << cell_lap[q]
                << " "  << cell_grad[q].norm_square()
                << "\n";

      double v = cell_lap[q] - 0.5*sin(2.0*cell_val[q]);
      I1 += v*v * fe_values.JxW(q);
    }
    std::cerr << "\n";
  }
  return sqrt(I1);
}


/*************************************************************************/
/*************************************************************************/

int main(){
  try {

    deallog.depth_console(0);
    TestWaveSolver ps;

    ps.make_initial_grid();
    ps.do_calc();

    double en_th=2.0*M_PI*M_PI;

    // loop for solving the wave equation and grid refinement
    for (int i=0; i<2; i++){
      ps.refine_grid(0.3, 0.03);
      ps.do_calc();
    }
    ps.save_grid("test_grid.eps");
    ps.save_wave("test_wave.eps");

    std::cout << "  result: " << std::setprecision(10) << ps.en
              << "  theory: " << std::setprecision(10) << en_th
              << "  diff: " << std::setprecision(10) << ps.en-en_th
              << std::endl;

    return 0;
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
}


