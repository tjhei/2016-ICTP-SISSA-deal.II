#include <iostream>

// a template function
template <typename T>
void yell (T test)
{
  test.shout("HI!");
}

// some class
class cat
{
public:
  void shout(const std::string &what)
  {
    std::cout << "cat says " << what << std::endl;
  }
};

class dog
{
public:
  void shout(const std::string &what)
  {
    std::cout << "woof!" << std::endl;
  }
};

int main()
{
  cat mycat;
  dog mydog;

  yell(mycat);
  yell(mydog);

}



