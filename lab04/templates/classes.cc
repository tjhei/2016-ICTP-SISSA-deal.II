#include <iostream>

// An example for a class template
template <int dim>
class Point
{
public:
  double elements[dim];

  void print();
};

template <int dim>
void Point<dim>::print()
{
  for (int i=0; i<dim; ++i)
    std::cout << elements[i] << " ";
  std::cout << std::endl;
}


// This is a template specialization of the template argument dim:
template <>
void Point<5>::print()
{
  std::cout << "Five dimensions are special!" << std::endl;
}

int main()
{
  Point<3> p;
  p.elements[0] = 42.0;
  p.print();
  Point<5> p2;
  p2.print();
}
