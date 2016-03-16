#include <grid/tria_boundary.h>

using namespace dealii;

template <int dim, int spacedim>
class NacaBoundary : public StraightBoundary<dim,spacedim>
{
  public:
				     // t: .12 per naca0012, c =
				     // lunghezza della corda
    NacaBoundary(double t, double c);

    
    virtual Point< spacedim > get_new_point_on_line
    (const typename Triangulation< dim, spacedim >::line_iterator
      &line) const;
    
    virtual void get_intermediate_points_on_line
    (const typename Triangulation< dim, spacedim >::line_iterator &line,
     std::vector< Point< spacedim > > &points) const;
    
    virtual Point< spacedim > 	get_new_point_on_quad (const typename
    Triangulation< dim, spacedim >::quad_iterator &quad) const;
    
    virtual void 	get_intermediate_points_on_quad (const
    typename Triangulation< dim, spacedim >::quad_iterator &quad,
    std::vector< Point< spacedim > > &points) const;
    
  private:
    double naca(double x) const;
    const double t;
    const double c;
};




template <int dim, int spacedim>
NacaBoundary<dim,spacedim>::NacaBoundary(double t=.12, double c=1):
		t(t),
		c(c)
{}


template <int dim, int spacedim>
double NacaBoundary<dim,spacedim>::naca(double xx) const
{
  double x = xx/c;
  double x2 = x*x;
  double x3 = x2*x;
  double x4 = x3*x;
  
  return ( t/.2* (.2969*std::sqrt(x) -
		  .1260*x -.3516*x2 +.2843*x3 -.1015*x4) );
}


template <int dim, int spacedim>
Point< spacedim > NacaBoundary<dim,spacedim>::get_new_point_on_line
(const typename Triangulation< dim, spacedim >::line_iterator &line)
  const 
{
				   // we always assume that naca is on
				   // xy plane.
  Point<spacedim> newp;
  newp = StraightBoundary<dim,spacedim>::get_new_point_on_line(line);
  double factor = (line->center()[1] > 0 ? 1.0 : -1.0);
  newp[1] = factor*naca(newp[0]);
  return newp;
}

template <int dim, int spacedim>
void NacaBoundary<dim,spacedim>::get_intermediate_points_on_line
(const typename Triangulation< dim, spacedim >::line_iterator &line,
 std::vector<Point<spacedim> > &points)
  const 
{
  StraightBoundary<dim,spacedim>::get_intermediate_points_on_line (line, points);
  
				   // we always assume that naca is on
				   // xy plane.

  double factor = line->center()[1] > 0 ? 1.0 : -1.0;
  for(unsigned int i=0; i<points.size(); ++i)
    points[i][1] = factor*naca(points[i][0]);
}


template <>
void NacaBoundary<2,3>::get_intermediate_points_on_quad
(const  Triangulation< 2, 3 >::quad_iterator &quad,
 std::vector<Point<3> > &points)
  const 
{
  StraightBoundary<2,3>::get_intermediate_points_on_quad (quad, points);
  
				   // we always assume that naca is on
				   // xy plane.

  double factor = quad->center()[1] > 0 ? 1.0 : -1.0;
  for(unsigned int i=0; i<points.size(); ++i)
    points[i][1] = factor*naca(points[i][0]);
}

template <>
Point<3> NacaBoundary<2,3>::get_new_point_on_quad
(const Triangulation< 2, 3 >::quad_iterator &quad)
  const 
{
				   // we always assume that naca is on
				   // xy plane.
  Point<3> newp =
    StraightBoundary<2,3>::get_new_point_on_quad(quad);
  double factor = quad->center()[1] > 0 ? 1.0 : -1.0;
  newp[1] = factor*naca(newp[0]);
  return newp;
}

    
template<>
void NacaBoundary<1, 2>::get_intermediate_points_on_quad
(dealii::TriaIterator<dealii::InvalidAccessor<2, 1, 2> > const&,
 std::vector<dealii::Point<2> >&) const
{
  Assert(false, ExcImpossibleInDim(1));
}

    
template<>
Point<2> NacaBoundary<1, 2>::get_new_point_on_quad
(dealii::TriaIterator<dealii::InvalidAccessor<2, 1, 2> > const&) const
{
  Assert(false, ExcImpossibleInDim(1));
  return Point<2>();
}
