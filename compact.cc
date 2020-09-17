#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "iRRAM.h"
using namespace iRRAM;

// The set of compact subsets of 1-dimensional real line
class OneDCompact
{

  
public:
  std::function< bool ( REAL, int ) > cfun;

  // create the emptyset
  OneDCompact();

  // create a singleton
  OneDCompact(REAL);  

  // create a closed interval
  OneDCompact(REAL, REAL);  

};


// functionals

// pointwise conjunction for binary Boolean functions
template <typename T, typename R>
std::function<bool ( T, R ) > pw_conjunction (std::function<bool ( T, R ) > f, std::function<bool ( T, R ) > g)
{
  return
    [=](T x, R y) -> bool {
      return (f(x, y) && g(x, y));
    };
}

// pointwise disjunction for binary Boolean functions
template <typename T, typename R>
std::function<bool ( T, R ) > pw_disjunction (std::function<bool ( T, R ) > f, std::function<bool ( T, R ) > g)
{
  return
    [=](T x, R y) -> bool {
      return (f(x,y) || g(x, y));
    };
}



// make the fuzzy char function for the emptyset
std::function<bool (REAL, int) > empty ()
{
  return
    [=](REAL y, int p) -> bool {
	return false;
    };
}


// make the fuzzy char function for the singleton {x}
std::function<bool (REAL, int) > singleton (REAL x)
{
  return
    [=](REAL y, int p) -> bool {
      REAL prec = scale(REAL(1), p);
      if (choose(iRRAM::abs(x - y) < prec, iRRAM::abs(x-y) > 0) == 1)	   
	return true;
      else
	return false;
    };
}



// make the fuzzy char function for the interval [x, y]
std::function<bool (REAL, int) > interval (REAL x, REAL y)
{
  return
    [=](REAL z, int p) -> bool {
	   REAL prec = scale(REAL(1), p); 
	   if (choose(
		      z < y+prec && z > x - prec,
		      z < x || z > y) == 1)
	     return true;
	   else
	     return false;
	 };
}
	   

OneDCompact::OneDCompact()
{
  cfun = empty(); 
}

OneDCompact::OneDCompact(REAL x)
{
  cfun = singleton(x); 
}
  
OneDCompact::OneDCompact(REAL x, REAL y)
{
  cfun = interval(x, y); 
}
  

// Primitive Operations: 
OneDCompact op_intersection(OneDCompact x, OneDCompact y)
{
  OneDCompact z;
  z.cfun = pw_conjunction<REAL, int>(x.cfun, y.cfun);
  return z;
}

OneDCompact op_union(OneDCompact x, OneDCompact y)
{
  OneDCompact z;
  z.cfun = pw_disjunction(x.cfun, y.cfun);
  return z;

}

bool is_member_of(REAL x, OneDCompact y, int p)
{
  return y.cfun(x, p);
}









// main function for testing
void compute()
{

  // closed interval [1, pi]
  OneDCompact x(1, pi());


  // test membership with precision -100
  if (x.cfun(3, -100)) cout << "3 is in [1, pi]\n"; else cout << "3 is not in [1, pi]\n";  

  if (x.cfun(4, -100)) cout << "4 is in [1, pi]\n"; else cout << "4 is not in [1, pi]\n";  


  // closed interval [3, 4]
  OneDCompact y(3,4);


  // test membership with precision -100
  if (y.cfun(3.2, -100)) cout << "3.2 is in [3, 4]\n"; else cout << "3.2 is not in [3, 4]\n";  

  if (y.cfun(0, -100)) cout << "0 is in [3, 4]\n"; else cout << "0 is not in [3, 4]\n";  

  
  
}
