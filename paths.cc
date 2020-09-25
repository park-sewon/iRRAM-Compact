#include <array>

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "iRRAM.h"


using namespace iRRAM;


// Copied from iRRAM/src/stack.cc and modified a little.
// The original file is buggy in the case when the function is a constant function. (it never stops!)
// this should be reported sometime.
int module(std::function<REAL (REAL)> f ,const REAL& x, int p){
// Semantics: If m=module(f,x,p), then
//    |x-z| <=2^m implies  |f(x)-f(z)| <= 2^p 

// If we are able to approximate f(x) by a DYADIC number d with an error of <=2^{p-1},
// then for any z in the argument interval of x, f(z) must differ by at most 
//  2^{p-1} from d, hence |f(x)-f(z)|<=2^p

  int result;
  if ( (ACTUAL_STACK.inlimit==0) && iRRAM_thread_data_address->cache_i.get(result)) return result;

  DYADIC d;
  REAL x_copy=x;

  
  sizetype argerror,testerror;
  
  x_copy.geterror(argerror);
  sizetype_set(testerror,1,argerror.exponent);
  x_copy.adderror(testerror);
  {
    single_valued code;
    d=approx(f(x_copy),p-1); 
  }
// At this line, we are sure that x_copy (and so also x) is precise enough to allow 
// the computation of f(x), even with a slightly increased error of the argument.

// We now try to find the smallest p_arg such that the evaluation of f(x+- 2^p_arg) 
// is possible up to an error of at most  2^{p-1}

  ITERATION_STACK SAVED_STACK;

// To do this, we start with p_arg=p.
// If this is successfull, we increase the value of p_arg until the first failure
// It it is not, then we decrease until the first success...
  int direction=0,p_arg=p;
  bool try_it=true;

  while (try_it) {


    sizetype_set(testerror,1,p_arg);
    x_copy.seterror(argerror);
    x_copy.adderror(testerror);
    bool fail = false;
  if ( iRRAM_unlikely(iRRAM_debug > 0 ) ) {
   sizetype x_error;
   x_copy.geterror(x_error);
  iRRAM_DEBUG2(1,"Testing module: 1*2^%d + %d*2^%d\n",p_arg,argerror.mantissa,argerror.exponent);
  iRRAM_DEBUG2(1,"argument error: %d*2^%d\n",x_error.mantissa,x_error.exponent);
  }
  try { 
      single_valued code;
      REAL z=f(x_copy);
      if ( iRRAM_unlikely(iRRAM_debug > 0 ) ) {
        sizetype z_error;
        z.geterror(z_error);
        iRRAM_DEBUG2(1,"Module yields result %d*2^%d\n",z_error.mantissa,z_error.exponent);
      }
      d=approx(z,p-1); 
	}
    catch ( Iteration it)  { fail=true; }
    switch ( direction ) {
      case 0:;
        if ( fail ) direction=-1; else direction=1;
	p_arg+=direction;
      break;
      case 1:;
        if ( fail || p_arg > 0) { try_it=false; p_arg -=direction; }
	else { p_arg += direction; };
      break;
      case -1:;
        if ( fail ) { p_arg +=direction; }
	else { try_it=false; };
      break;    
    }
    }
  iRRAM_DEBUG2(1,"Modules resulting in p_arg=%d\n",p_arg);
  
  sizetype_set(testerror,1,p_arg);
  sizetype_inc(argerror,testerror);

  while (argerror.mantissa>1) {
    argerror.mantissa=argerror.mantissa>>1;
    argerror.exponent+=1;
  }
  
  result=argerror.exponent;
  if ( ACTUAL_STACK.inlimit==0 ) iRRAM_thread_data_address->cache_i.put(result);
  return result;

}



// Euclidean space using alias declaration
template <int N>
using IR = std::array<REAL, N>;

// Hypercube as a (mental) subset 
template <int N>
using HyperCube = IR<N>;

template <int N>
IR<N> IR_origin (){
  std::array<REAL, N> res;
  for (int i = 0; i < N; i ++)
    res = REAL(0);
  return res;
}

// metric
template <int N>
REAL IR_d (IR<N> x, IR<N> y){
  REAL sum = 0;
  REAL a, b;
  for (int i = 0; i < N; i ++){
    a = x[i];
    b = y[i];
    sum  += (a - b) * (a - b);
  }
  return sqrt(sum);
}


// homotopy is a function from a hypercube
template <int N, int M>
using Homotopy = std::function<HyperCube<M> ( IR<N> ) >;



// return 2^n
RATIONAL Exp(int n){
  if(n > 0)
    return (INTEGER(1) << n);
  else if(n == 0)
    return 1;
  else 
    return (RATIONAL(INTEGER(1), INTEGER(1) << -n));
}



/*
  Example: [0,1] → R ≃ HyperCube<1> → IR<1> ≃ Homotopy<1, 1>
  
  as_fun casts Homotopy<1,1> to a function R → R

  OneDMin_approx approximates the minimum value of f in [0,1] by 2ᵖ
  OneDMin computes the minimum value of f in [0,1]


  OneDMax_approx approximates the minimum value of f in [0,1] by 2ᵖ
  OneDMax computes the maximum value of f in [0,1]

*/


std::function<REAL (REAL)> as_func(Homotopy<1, 1> f){
  return
    [=](REAL x) -> REAL {
      return f({x})[0];
    };
}


// min of f : [0, 1] -> R
REAL OneDMin_approx (int p, std::function<REAL (REAL)>  f){
  RATIONAL x = 0;
  REAL m = f(x);
  int q;
  while(x < RATIONAL(1, 1)){
    // cout<<"testing" <<REAL(x) <<"\n";
    q = module( f, x, p);
    m = minimum( f(x), m);
    x = x + Exp(q); 
  }
  return m;
}

std::function<REAL (const int &)> cast_min (Homotopy<1,1> f){
  return
    (
     [=](int x) -> REAL {
       return OneDMin_approx(x, as_func(f));
     }
     );
}


REAL OneDMin (Homotopy<1, 1> f){
  return limit(from_algorithm<REAL, int>(cast_min(f)));
}


REAL OneDMax_approx (int p, std::function<REAL (REAL)>  f){
  RATIONAL x = 0;
  REAL m = f(x);
  int q;
  while(x < RATIONAL(1, 1)){
    q = module( f, x, p);
    m = maximum( f(x), m);
    x = x + Exp(q); 
  }
  return m;
}

std::function<REAL (const int &)> cast_max (Homotopy<1,1> f){
  return
    (
     [=](int x) -> REAL {
       return OneDMax_approx(x, as_func(f));
     }
     );
}


REAL OneDMax (Homotopy<1, 1> f){
  return limit(from_algorithm<REAL, int>(cast_max(f)));
}



/*
 Test functions from [0,1] → R

 test_one: x ↦ x + 10

 test_two: x ↦ x²

 test_three: x ↦ sin(x)

 test_four: x ↦ 42
 */
Homotopy<1, 1> test_one ()
{
  return (
	  [=](HyperCube<1> x) -> IR<1> {
	    return {x[0] + 10};
	  }
	  );
}

Homotopy<1, 1> test_two ()
{
  return (
	  [=](HyperCube<1> x) -> IR<1> {
	    return {x[0] * x[0]};
	  }
	  );
}

Homotopy<1, 1> test_three ()
{
  return (
	  [=](HyperCube<1> x) -> IR<1> {
	    return {sin(x[0])};
	  }
	  );
}

Homotopy<1, 1> test_four ()
{
  return (
	  [=](HyperCube<1> x) -> IR<1> {
	    return {42};
	  }
	  );
}


// main function for testing
void compute()
{
  //cout << module1(test_t_one, 0, -3);

  //  cout << 42<<"\n";

  cout << OneDMin(test_four()) << "\n";
  
}
