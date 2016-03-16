#include <iostream>

// a template function
template<typename number>
number square(const number x)
{
  return x*x;
}

class test
{
  int my_member_variable;
};

int main()
{
  int i = 2;
  double d = 3.5;

  int i2 = square<int>(i); // or: square(i)
  double d2 = square<double>(d); // or: square(d)

  std::cout << i << " -> " << i2 << std::endl;
  std::cout << d << " -> " << d2 << std::endl;

  //test t1;
  //square(t1); // error
}
