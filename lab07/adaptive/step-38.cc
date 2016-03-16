/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2010 - 2015 by the deal.II authors
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
 * Authors: Andrea Bonito, Sebastian Pauletti.
 */



#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h>

#include <fstream>
#include <iostream>


namespace Step38
{
  using namespace dealii;



  class KinkedManifold : public ChartManifold<2,3,2>
  {
  public:
    static const int chartdim = 2;
    static const int spacedim = 3;

    virtual
    Point<spacedim> project_to_manifold (const std::vector<Point<spacedim> > &/*surrounding_points*/,
                                         const Point<spacedim> &candidate) const
    {
      return push_forward(pull_back(candidate));
    }

    /**
     * Given a point on the manifold, pull it back to the
     * chartdim dimensional Euclidean space:
     */
    virtual Point<chartdim>
    pull_back(const Point<spacedim> &space_point) const
    {
      return Point<chartdim>(space_point(0), space_point(1));
    }

    /**
     * Given a point in the chartdim dimensional Euclidean space, this method
     * returns a point on the manifold embedded in the spacedim Euclidean space.
     */
    virtual Point<spacedim>
    push_forward(const Point<chartdim> &chart_point) const
    {
      double x = chart_point(0);
      double y = chart_point(1);
      double alpha = 3.0/2.0; //2.0/5.0;//3.0/2.0;
      double z = pow(std::max(3.0/4.0 - x*x-y*y, 0.0),1.0+alpha);
      return Point<spacedim>(x, y, z);
    }

  };



  template <int spacedim>
  class LaplaceBeltramiProblem
  {
  public:
    LaplaceBeltramiProblem (const unsigned mapping_degree, const unsigned degree);
    void run ();

  private:
    static const unsigned int dim = spacedim-1;

    void make_grid ();
    void setup ();
    void assemble_system ();
    void solve ();
    void refine ();
    void output_results (unsigned int cycle) const;
    void compute_error () const;


    Triangulation<dim,spacedim>   triangulation;
    FE_Q<dim,spacedim>            fe;
    DoFHandler<dim,spacedim>      dof_handler;
    MappingQ<dim, spacedim>       mapping;
    ConstraintMatrix              constraints;

    SparsityPattern               sparsity_pattern;
    SparseMatrix<double>          system_matrix;

    Vector<double>                solution;
    Vector<double>                system_rhs;
  };



  template <int dim>
  class Solution  : public Function<dim>
  {
  public:
    Solution () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component = 0) const;

  };


  template <>
  double
  Solution<2>::value (const Point<2> &p,
                      const unsigned int) const
  {
    return 0;
  }


  template <>
  Tensor<1,2>
  Solution<2>::gradient (const Point<2>   &p,
                         const unsigned int) const
  {
    Tensor<1,2> return_value;
    return_value[0] = -2. * p(1) * (1 - 2. * p(0) * p(0));
    return_value[1] = -2. * p(0) * (1 - 2. * p(1) * p(1));

    return return_value;
  }


  template <>
  double
  Solution<3>::value (const Point<3> &p,
                      const unsigned int) const
  {
    double x = p(0);
    double y = p(1);
    return x*x+y*y;
  }


  template <>
  Tensor<1,3>
  Solution<3>::gradient (const Point<3>   &p,
                         const unsigned int) const
  {
    using numbers::PI;
    double x = p(0);
    double y = p(2);
    double z = p(1);
    double r = 0.6;
    double R = 1.0;
    double phi = atan2(y,x);
    double theta = asin(z/r);
    if (x*x+y*y < R*R)
      theta = ((z>0)?1.0:-1.0)*numbers::PI - theta;


    Tensor<1,3> return_value;
/*
    return_value[0] = (-3*sin(phi)*cos(theta)*cos(3*phi)*cos(3*theta+phi)*cos(phi)*r+sin(phi)*cos(theta)*sin(3*phi)*sin(3*theta+phi)*cos(phi)*r-3*sin(phi)*cos(theta)*cos(3*phi)*cos(3*theta+phi)*R+sin(phi)*cos(theta)*sin(3*phi)*sin(3*theta+phi)*R+3*sin(3*phi)*sin(3*theta+phi)*sin(theta)*r)/((r*cos(phi)+R)*r);

    return_value[2] = -(3*sin(theta)*sin(phi)*cos(3*phi)*cos(3*theta+phi)*cos(phi)*r-sin(theta)*sin(phi)*sin(3*phi)*sin(3*theta+phi)*cos(phi)*r+3*sin(theta)*sin(phi)*cos(3*phi)*cos(3*theta+phi)*R-sin(theta)*sin(phi)*sin(3*phi)*sin(3*theta+phi)*R+3*sin(3*phi)*sin(3*theta+phi)*cos(theta)*r)/((r*cos(phi)+R)*r);

    return_value[1] = (3*cos(3*phi)*cos(3*theta+phi)-sin(3*phi)*sin(3*theta+phi))*cos(phi)/r;
*/
    return return_value;

   /* Tensor<2,3> hessian;

       hessian[0][0] = -PI*PI*sin(PI*p(0))*cos(PI*p(1))*exp(p(2));
       hessian[1][1] = -PI*PI*sin(PI*p(0))*cos(PI*p(1))*exp(p(2));
       hessian[2][2] = sin(PI*p(0))*cos(PI*p(1))*exp(p(2));

       hessian[0][1] = -PI*PI*cos(PI*p(0))*sin(PI*p(1))*exp(p(2));
       hessian[1][0] = -PI*PI*cos(PI*p(0))*sin(PI*p(1))*exp(p(2));

       hessian[0][2] = PI*cos(PI*p(0))*cos(PI*p(1))*exp(p(2));
       hessian[2][0] = PI*cos(PI*p(0))*cos(PI*p(1))*exp(p(2));

       hessian[1][2] = -PI*sin(PI*p(0))*sin(PI*p(1))*exp(p(2));
       hessian[2][1] = -PI*sin(PI*p(0))*sin(PI*p(1))*exp(p(2));

       Tensor<1,3> gradient;
       gradient[0] = PI * cos(PI*p(0))*cos(PI*p(1))*exp(p(2));
       gradient[1] = - PI * sin(PI*p(0))*sin(PI*p(1))*exp(p(2));
       gradient[2] = sin(PI*p(0))*cos(PI*p(1))*exp(p(2));

       Point<3> normal = p;
       normal /= p.norm();

       return (- trace(hessian)
               + 2 * (gradient * normal)
               + (hessian * normal) * normal); */
  }



  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
  };

  template <>
  double
  RightHandSide<2>::value (const Point<2> &p,
                           const unsigned int /*component*/) const
  {
    return ( -8. * p(0) * p(1) );
  }


  template <>
  double
  RightHandSide<3>::value (const Point<3> &p,
                           const unsigned int /*component*/) const
  {
    double x = p(0);
    double y = p(1);
    double rhs;

    rhs = (0.0e0 < 0.3e1 / 0.4e1 - x * x - y * y ? -0.128e3 * (0.3200e4 * pow(x, 0.8e1) + 0.12800e5 * pow(x, 0.6e1) * y * y + 0.19200e5 * pow(x, 0.4e1) * pow(y, 0.4e1) + 0.12800e5 * x * x * pow(y, 0.6e1) + 0.3200e4 * pow(y, 0.8e1) - 0.3600e4 * pow(x, 0.6e1) - 0.10800e5 * pow(x, 0.4e1) * y * y - 0.10800e5 * x * x * pow(y, 0.4e1) - 0.3600e4 * pow(y, 0.6e1) + 0.675e3 * x * x + 0.675e3 * y * y + 0.128e3) * pow(0.1600e4 * pow(x, 0.8e1) + 0.6400e4 * pow(x, 0.6e1) * y * y + 0.9600e4 * pow(x, 0.4e1) * pow(y, 0.4e1) + 0.6400e4 * x * x * pow(y, 0.6e1) + 0.1600e4 * pow(y, 0.8e1) - 0.3600e4 * pow(x, 0.6e1) - 0.10800e5 * pow(x, 0.4e1) * y * y - 0.10800e5 * x * x * pow(y, 0.4e1) - 0.3600e4 * pow(y, 0.6e1) + 0.2700e4 * pow(x, 0.4e1) + 0.5400e4 * x * x * y * y + 0.2700e4 * pow(y, 0.4e1) - 0.675e3 * x * x - 0.675e3 * y * y - 0.64e2, -0.2e1) : -4);
    return rhs;
    /*using numbers::PI;
    double x = p(0);
    double y = p(2);
    double z = p(1);
    double r = 0.6;
    double R = 1.0;
    double phi = atan2(y,x);
    double theta = asin(z/r);
    if (x*x+y*y < R*R)
      theta = ((z>0)?1.0:-1.0)*numbers::PI - theta;
    double ct=cos(theta);
    double st=sin(theta);
    double c3p=cos(3*phi);
    double s3p=sin(3*phi);
    double c3tp=cos(3*theta+phi);
    double s3tp=sin(3*theta+phi);
    double f =( 9*r*r*ct*ct*s3p*c3tp
                                  -3*r*r*ct*st*s3p*s3tp
                               +18*r*R*ct*s3p*c3tp
                                  -3*r*R*st*s3p*s3tp
                                 +6*r*r*c3p*s3tp
                                 +9*R*R*s3p*c3tp
                               +10*r*r*s3p*c3tp ) /(r*r*( r*r*ct*ct
                                                                         +2*r*R*ct
                                                                         +R*R      ) );
    return f;*/
  }



  template <int spacedim>
  LaplaceBeltramiProblem<spacedim>::
  LaplaceBeltramiProblem (const unsigned mapping_degree, const unsigned degree)
    :
    fe (degree),
    dof_handler(triangulation),
    mapping (mapping_degree)
  {}




  template <int spacedim>
  void LaplaceBeltramiProblem<spacedim>::make_grid ()
  {
    static KinkedManifold surface_description;
    GridGenerator::hyper_cube(triangulation, -1, 1);

    triangulation.set_all_manifold_ids(0);
    triangulation.set_manifold (0, surface_description);

    triangulation.refine_global(3);


#if 0
    static TorusBoundary<dim,spacedim> surface_description( 1.0, 0.6);
    GridGenerator::torus(triangulation, 1.0, 0.6);

    triangulation.set_all_manifold_ids(0);
    triangulation.set_manifold (0, surface_description);

    triangulation.refine_global(2);

#endif
  }

  template <int spacedim>
  void LaplaceBeltramiProblem<spacedim>::setup ()
  {


    dof_handler.distribute_dofs (fe);

    constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler,
                                                constraints);
    VectorTools::interpolate_boundary_values (mapping,
                                              dof_handler,
                                              0,
                                              Solution<spacedim>(),
                                              constraints);
    constraints.close ();

    DynamicSparsityPattern dsp (dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints);
    sparsity_pattern.copy_from (dsp);

    system_matrix.reinit (sparsity_pattern);

    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
  }



  template <int spacedim>
  void LaplaceBeltramiProblem<spacedim>::assemble_system ()
  {
    system_matrix = 0;
    system_rhs = 0;

    const QGauss<dim>  quadrature_formula(2*fe.degree);
    FEValues<dim,spacedim> fe_values (mapping, fe, quadrature_formula,
                                      update_values              |
                                      update_gradients           |
                                      update_quadrature_points   |
                                      update_JxW_values);

    const unsigned int        dofs_per_cell = fe.dofs_per_cell;
    const unsigned int        n_q_points    = quadrature_formula.size();

    FullMatrix<double>        cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>            cell_rhs (dofs_per_cell);

    std::vector<double>       rhs_values(n_q_points);
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    const RightHandSide<spacedim> rhs;

    for (typename DoFHandler<dim,spacedim>::active_cell_iterator
         cell = dof_handler.begin_active(),
         endc = dof_handler.end();
         cell!=endc; ++cell)
      {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit (cell);

        rhs.value_list (fe_values.get_quadrature_points(), rhs_values);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
              cell_matrix(i,j) += fe_values.shape_grad(i,q_point) *
                                  fe_values.shape_grad(j,q_point) *
                                  fe_values.JxW(q_point);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            cell_rhs(i) += fe_values.shape_value(i,q_point) *
                           rhs_values[q_point]*
                           fe_values.JxW(q_point);

        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix,
                                               cell_rhs,
                                               local_dof_indices,
                                               system_matrix,
                                               system_rhs);

      }

}




  template <int spacedim>
  void LaplaceBeltramiProblem<spacedim>::solve ()
  {
    SolverControl solver_control (solution.size(),
                                  1e-12 * system_rhs.l2_norm());
    SolverCG<>    cg (solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve (system_matrix, solution, system_rhs,
              preconditioner);

    constraints.distribute(solution);

    double mean = VectorTools::compute_mean_value (dof_handler,
                                                 QGauss<2>(2*fe.degree),
                                                 solution,
                                                 0);

    /*solution.add(-mean);
    std::cout << "adjusting mean from " << mean
              << " to "
    << VectorTools::compute_mean_value (dof_handler,
                                                     QGauss<2>(2*fe.degree),
                                                     solution,
                                                     0)
    << std::endl;*/
  }

  template <int spacedim>
  void LaplaceBeltramiProblem<spacedim>::refine ()
  {
    static int counter = 0;
    ++counter;
    if (0)
      {
        triangulation.refine_global (1);

      }
    else if (counter%2==0)
      {
        Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

        KellyErrorEstimator<dim,spacedim>::estimate (mapping,
                                            dof_handler,
                                            QGauss<dim-1>(2*fe.degree+1),
                                            typename FunctionMap<spacedim>::type(),
                                            solution,
                                            estimated_error_per_cell);

        GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                         estimated_error_per_cell,
                                                         0.3, 0.0);

        triangulation.execute_coarsening_and_refinement ();
      }
    else if (1)
      {
        Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

        const QGauss<dim>  quadrature_formula(2*fe.degree);
        FEValues<dim,spacedim> fe_values (mapping, fe, quadrature_formula,
                                          update_values              |
                                          update_gradients           |
                                          update_quadrature_points   |
                                          update_JxW_values);

        const unsigned int        n_q_points    = quadrature_formula.size();

        for (typename DoFHandler<dim,spacedim>::active_cell_iterator
             cell = dof_handler.begin_active(),
             endc = dof_handler.end();
             cell!=endc; ++cell)
          {
            fe_values.reinit (cell);

            double dist = 0.0;
            for (unsigned int q=0;q<n_q_points;++q)
              {
                Point<spacedim> p = fe_values.quadrature_point(q);
                std::vector<Point<spacedim> > temp(1, p);
                Point<spacedim> p_exact = triangulation.get_manifold(0).project_to_manifold(temp, p);
                dist = std::max(dist, p.distance(p_exact));
              }

            estimated_error_per_cell[cell->active_cell_index()] = dist;
          }

        GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                         estimated_error_per_cell,
                                                         0.3, 0.0);

        triangulation.execute_coarsening_and_refinement ();
      }

  }


  template <int spacedim>
  void LaplaceBeltramiProblem<spacedim>::output_results (unsigned int cycle) const
  {
    DataOut<dim,DoFHandler<dim,spacedim> > data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution,
                              "solution",
                              DataOut<dim,DoFHandler<dim,spacedim> >::type_dof_data);
    data_out.build_patches (mapping,
                            mapping.get_degree());
    std::string filename = spacedim == 2 ?
          "solution-2d-" :
          "solution-3d-";
    filename += Utilities::int_to_string(cycle,2) + ".vtk";
    std::ofstream output (filename.c_str());
    data_out.write_vtk (output);
  }





  template <int spacedim>
  void LaplaceBeltramiProblem<spacedim>::compute_error () const
  {
    Vector<float> difference_per_cell (triangulation.n_active_cells());
    VectorTools::integrate_difference (mapping, dof_handler, solution,
                                       Solution<spacedim>(),
                                       difference_per_cell,
                                       QGauss<dim>(2*fe.degree+1),
                                       VectorTools::L2_norm);
    double l2 = difference_per_cell.l2_norm();
    VectorTools::integrate_difference (mapping, dof_handler, solution,
                                       Solution<spacedim>(),
                                       difference_per_cell,
                                       QGauss<dim>(2*fe.degree+1),
                                       VectorTools::H1_norm);
    double h1 = difference_per_cell.l2_norm();
    double h = GridTools::minimal_cell_diameter(triangulation);
    std::cout << "  h= " << h
              << " n= " << dof_handler.n_dofs()
              << " L2error= " << l2
              << " H1error= " << h1
              << std::endl;
  }




  template <int spacedim>
  void LaplaceBeltramiProblem<spacedim>::run ()
  {
    std::cout << "Running in dim="<< dim << " spacedim=" << spacedim
              << " with mapping degree m=" << mapping.get_degree()
              << " and FE degree k=" << fe.degree << std::endl;

    make_grid();
    for (unsigned int cycle = 0; cycle < 11; ++cycle)
      {
        std::cout << "** Cycle " << cycle << std::endl;
        if (cycle>0)
          refine ();

        setup ();
        assemble_system ();
        solve ();
        output_results (cycle);
        compute_error ();
      }
  }
}



int main ()
{
  try
    {
      using namespace dealii;
      using namespace Step38;

    LaplaceBeltramiProblem<3> laplace_beltrami(1, 1);
    laplace_beltrami.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
