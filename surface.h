#pragma once

#include <array>

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "iRRAM.h"

using namespace iRRAM;

#include "compact.h"
#include "euclidean.h"
#include "plot.h"

// [0,1]^2 -> R^N
template <int N>
class Surface : public Compact<N> {
public:
  // original function
  std::function<Point<N>(REAL,REAL)> f;

  // init
  Surface(std::function<Point<N>(REAL,REAL)> f) { this->f = f; }

  // whenever |x-z| < 2^-pArg, |f(x)-f(z)| < 2^-p for all x,z
  // will be increased when higher precision is requested
  int p=INT_MIN, pArg=INT_MIN;

  // increase the current precision(from this->p to p)
  // and find the corresponding pArg
  void increasePrecision(int p) {
    // ignore lower or equal precision
    if(this->p >= p) return;

    // Find the pArg
    // must |f(u)-f(z)| < (2^-p)/sqrt(2) to include the box with a ball
    // hence find find pArg such that whenever |x-z| < 2^-pArg, |f(x)-f(z)| < 2^(-p-1) for all x,z
    // This makes the distance between two consecutive centers of balls 2^(-p-1)*sqrt(2) at most.
    std::function<Point<N>(Point<2>)> ff = [=](Point<2> x) -> Point<N> {
      return this->f(x[0], x[1]);
    };
    this->pArg = module2<2,N>(ff, p+1);

    // update the current precision
    this->p = p;

    // update the current characteristic function
    this->cfun = [&](Point<N> pt, int p) {
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
        for(v=step/2 ; v<=ONE ; v+=step) {
          d = IR_d<N>(pt, this->f(u, v));
          if(choose(d < radius, d > radiusHalf) == 1) return true;
        }
      }
      return false;
    };
  }
  
  // membership test for point with precision 2^-p
  // in accordance to the Ko compatibility
  bool member(Point<N> point, int p) {
    // check if previously found pArg is viable
    // if not, increase the precision
    if(this->p < p) this->increasePrecision(p);

    return this->cfun(point, p);
  }
};

