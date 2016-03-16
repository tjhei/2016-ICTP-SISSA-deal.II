/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Timo Heister, based on step-4
 */



#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>

#include <deal.II/base/logstream.h>

using namespace dealii;


template <int dim>
class Step4
{
public:
  Step4 ();
  void run ();

private:
  void make_grid ();
  void setup_system();
  void assemble_system ();
  void solve ();
  void output_results (unsigned int cycle) const;

  Triangulation<dim>   triangulation;
  FE_Q<dim>            fe;
  MappingQ<dim>        mapping;
  unsigned int quad_degree;
  DoFHandler<dim>      dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double>       solution;
  Vector<double>       system_rhs;
};



template <int dim>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};

template <int dim>
double RightHandSide<dim>::value (const Point<dim> &p,
                                  const unsigned int /*component*/) const
{
  double x = p(0);
  double y = p(1);
  double r = sqrt(x*x+y*y);
  double theta = atan2(y,x);
  return 1.0;
}



template <int dim>
class SolutionValues : public Function<dim>
{
public:
  SolutionValues () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;

  virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                  const unsigned int  component = 0) const;
};

template <int dim>
double SolutionValues<dim>::value (const Point<dim> &p,
                                   const unsigned int /*component*/) const
{
  double x = p(0);
  double y = p(1);
  double r = sqrt(x*x+y*y);
  double theta = atan2(y,x);
  if (dim==2)
    return r*r*sin(3.0*theta);
  else
    return 0.0; // TODO
}

template <int dim>
Tensor<1,dim> SolutionValues<dim>::gradient (const Point<dim> &p,
                                   const unsigned int /*component*/) const
{
  Tensor<1,dim> return_value;
  double x = p(0);
  double y = p(1);
  double r = sqrt(x*x+y*y);
  double theta = atan2(x,y);
  if (dim==2)
    {
      return_value[0] = 0;
      return_value[1] = 0;
    }
  else
    {
      // TODO
    }

  return return_value;
}





template <int dim>
Step4<dim>::Step4 ()
  :
  fe (2),
  mapping (1),
  quad_degree(2*fe.degree+2),
  dof_handler (triangulation)
{}



template <int dim>
void Step4<dim>::make_grid ()
{
  Point<dim> center;
  GridGenerator::hyper_ball (triangulation, center);
  //triangulation.set_all_manifold_ids(0);
  triangulation.set_all_manifold_ids_on_boundary(0);

  if (0)
    {
      typename Triangulation<dim>::active_cell_iterator cell =
          triangulation.begin_active();
      for (;cell!=triangulation.end();++cell)
        if (cell->center().square()<1e-3)
          cell->set_all_manifold_ids(1);
    }

  static SphericalManifold<dim> manifold_description(center);
  triangulation.set_manifold(0, manifold_description);
  triangulation.refine_global (1);
}


template <int dim>
void Step4<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  std::cout << "   Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}



template <int dim>
void Step4<dim>::assemble_system ()
{
  QGauss<dim>  quadrature_formula(quad_degree);

  const RightHandSide<dim> right_hand_side;

  FEValues<dim> fe_values (mapping, fe, quadrature_formula,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();

  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell_matrix = 0;
      cell_rhs = 0;

      for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_grad (i, q_index) *
                                   fe_values.shape_grad (j, q_index) *
                                   fe_values.JxW (q_index));

            cell_rhs(i) += (fe_values.shape_value (i, q_index) *
                            right_hand_side.value (fe_values.quadrature_point (q_index)) *
                            fe_values.JxW (q_index));
          }

      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i,j));

          system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }


  std::map<types::global_dof_index,double> boundary_values;
  VectorTools::interpolate_boundary_values (mapping,
                                            dof_handler,
                                            0,
                                            SolutionValues<dim>(),
                                            boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
                                      system_matrix,
                                      solution,
                                      system_rhs);
}



template <int dim>
void Step4<dim>::solve ()
{
  SolverControl           solver_control (2000, 1e-12);
  SolverCG<>              solver (solver_control);
  PreconditionSSOR<SparseMatrix<double> > prec;
  prec.initialize(system_matrix);
  solver.solve (system_matrix, solution, system_rhs,
                prec);

  std::cout << "   " << solver_control.last_step()
            << " CG iterations needed to obtain convergence."
            << std::endl;
}



template <int dim>
void Step4<dim>::output_results (unsigned int cycle) const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");

  data_out.build_patches ();

  std::string filename = dim == 2 ?
        "solution-2d-" :
        "solution-3d-";
  filename += Utilities::int_to_string(cycle,2) + ".vtk";
  std::ofstream output (filename.c_str());
  data_out.write_vtk (output);

  Vector<float> difference_per_cell (triangulation.n_active_cells());
  VectorTools::integrate_difference (mapping,
                                     dof_handler,
                                     solution,
                                     SolutionValues<dim>(),
                                     difference_per_cell,
                                     QGauss<dim>(quad_degree+1),
                                     VectorTools::L2_norm);
  const double L2_error = difference_per_cell.l2_norm();

  VectorTools::integrate_difference (mapping,
                                     dof_handler,
                                     solution,
                                     SolutionValues<dim>(),
                                     difference_per_cell,
                                     QGauss<dim>(quad_degree+1),
                                     VectorTools::H1_norm);
  const double H1_error = difference_per_cell.l2_norm();

  std::cout << "  h= " << triangulation.begin_active()->diameter()
               << "  L2= " << L2_error
            << "  H1= " // << H1_error
            << std::endl;
}




template <int dim>
void Step4<dim>::run ()
{
  std::cout << "Solving problem in " << dim << " space dimensions." << std::endl;

  make_grid();
  for (unsigned int cycle = 0; cycle < 6; ++cycle)
    {
      std::cout << "** Cycle " << cycle << std::endl;
      if (cycle>0)
        triangulation.refine_global(1);

      setup_system ();
      assemble_system ();
      solve ();
      output_results (cycle);
    }
}


int main ()
{
  {
    Step4<2> laplace_problem_2d;
    laplace_problem_2d.run ();
  }

//  {
//    Step4<3> laplace_problem_3d;
//    laplace_problem_3d.run ();
//  }

  return 0;
}
