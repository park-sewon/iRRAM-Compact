#pragma once

#include <array>

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "iRRAM.h"
#include "compact.h"
#include "euclidean.h"
#include "plot.h"

using namespace iRRAM;

// [0,1] -> R^N
template <int N>
class Path : public Compact<N> {
public:
  // original function
  std::function<Point<N>(REAL)> f;

  // init
  Path(std::function<Point<N>(REAL)> f) { this->f = f; }

  // whenever |x-z| < 2^-pArg, |f(x)-f(z)| < 2^-p for all x,z
  // will be increased when higher precision is requested
  int p=INT_MIN, pArg=INT_MIN;

  // increase the current precision(from this->p to p)
  // and find the corresponding pArg
  void increasePrecision(int p) {
    // ignore lower or equal precision
    if(this->p >= p) return;

    // find the pArg
    // must |f(u)-f(z)| < (2^-p)/sqrt(2) to include the box with a ball
    // hence find find pArg such that whenever |x-z| < 2^-pArg, |f(x)-f(z)| < 2^(-p-1) for all x,z
    // This makes the distance between two consecutive centers of balls 2^(-p-1)*sqrt(2) at most.
    std::function<Point<N>(Point<1>)> ff = [&](Point<1> x) -> Point<N> {
      return this->f(x[0]);
    };
    this->pArg = module2<1,N>(ff, p+1);

    // update the current precision
    this->p = p;

    // update the current characteristic function
    this->cfun = [&](Point<N> pt, int p) -> bool{
      single_valued code;

      // check the membership with previously found p and pArg
      // centers of balls: f(step/2 + i*step, step/2 + j*step)
      // radius of a ball: 2^-p
      // maximum distance between two consecutive balls(centers): 2^(-p-1)*sqrt(2)   (check increasePrecision())
      // When checking inclusiveness, we run (d < radius) and (d > radius/2) in parallel for decidability.
      // For any point on the path, there exists a ball that contains the point.
      RATIONAL step = Exp(-pArg);   // 2^-pArg
      RATIONAL radius = Exp(-p);     // 2^-p, the radius of a ball
      RATIONAL radiusHalf = radius/INTEGER(2);     // 2^(-p-1)
      REAL d;
      RATIONAL u,v;
      for(u=step/2 ; u<=ONE ; u+=step) {
        d = IR_d<N>(pt, this->f(u));
        if(choose(d < radius, d > radiusHalf) == 1) return true;
      }
      return false;
    };
  }
  
  // membership test for point with precision 2^-p
  // in accordance to the Ko compatibility
  bool member (Point<N> point, int p) {
    // cout << "member in path p = " << p << "\n";

    // check if previously found pArg is viable
    // if not, increase the precision
    if(this->p < p) this->increasePrecision(p);
    
    return this->cfun(point, p);
  }
};


// homotopy is a function from a hypercube
template <int N, int M>
using Homotopy = std::function<HyperCube<M>(IR<N>)>;

/*
  Example: [0,1] → R ≃ HyperCube<1> → IR<1> ≃ Homotopy<1, 1>
  
  as_fun casts Homotopy<1,1> to a function R → R

  OneDMin_approx approximates the minimum value of f in [0,1] by 2ᵖ
  OneDMin computes the minimum value of f in [0,1]


  OneDMax_approx approximates the minimum value of f in [0,1] by 2ᵖ
  OneDMax computes the maximum value of f in [0,1]

*/

std::function<REAL(REAL)> as_func(Homotopy<1, 1> f)
{
  return
      [=](REAL x) -> REAL {
        return f({x})[0];
      };
}

// min of f : [0, 1] -> R
REAL OneDMin_approx(int p, std::function<REAL(REAL)> f)
{
  RATIONAL x = 0;
  REAL m = f(x);
  int q;
  while (x < RATIONAL(1, 1))
  {
    // cout<<"testing" <<REAL(x) <<"\n";
    q = module(f, x, p);
    m = minimum(f(x), m);
    x = x + Exp(q);
  }
  return m;
}

std::function<REAL(const int &)> cast_min(Homotopy<1, 1> f)
{
  return (
      [=](int x) -> REAL {
        return OneDMin_approx(x, as_func(f));
      });
}

REAL OneDMin(Homotopy<1, 1> f)
{
  return limit(from_algorithm<REAL, int>(cast_min(f)));
}

REAL OneDMax_approx(int p, std::function<REAL(REAL)> f)
{
  RATIONAL x = 0;
  REAL m = f(x);
  int q;
  while (x < RATIONAL(1, 1))
  {
    q = module(f, x, p);
    m = maximum(f(x), m);
    x = x + Exp(q);
  }
  return m;
}

std::function<REAL(const int &)> cast_max(Homotopy<1, 1> f)
{
  return (
      [=](int x) -> REAL {
        return OneDMax_approx(x, as_func(f));
      });
}

REAL OneDMax(Homotopy<1, 1> f)
{
  return limit(from_algorithm<REAL, int>(cast_max(f)));
}

/*
 Test functions from [0,1] → R

 test_one: x ↦ x + 10

 test_two: x ↦ x²

 test_three: x ↦ sin(x)

 test_four: x ↦ 42
 */
Homotopy<1, 1> test_one()
{
  return (
      [=](HyperCube<1> x) -> IR<1> {
        return {x[0] + 10};
      });
}

Homotopy<1, 1> test_two()
{
  return (
      [=](HyperCube<1> x) -> IR<1> {
        return {x[0] * x[0]};
      });
}

Homotopy<1, 1> test_three()
{
  return (
      [=](HyperCube<1> x) -> IR<1> {
        return {sin(x[0])};
      });
}

Homotopy<1, 1> test_four()
{
  return (
      [=](HyperCube<1> x) -> IR<1> {
        return {42};
      });
}