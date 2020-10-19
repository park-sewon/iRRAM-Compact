#include <array>
#include <algorithm>

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "iRRAM.h"

using namespace iRRAM;

#include "path.h"
#include "surface.h"


void compute() {
  // Path<2> path1([=](REAL t) -> Point<2> { return Point<2>({t,t}); });
  // path1.plot2D("t.png", 100, 0.5, 1.5, -0.5, 1.5);


  // Path<2> path2([=](REAL t) -> Point<2> { return Point<2>({t,t*t}); });
  // path2.plot2D("t.png", 100, -0.1, 2.1, -0.1, 4.1);

  // Path<2> path3([=](REAL t) -> Point<2> { return Point<2>({2*pi()*t,sin(2*pi()*t)}); });
  // path3.plot2D("t.png", 300, 0, 2*pi(), -1, 1);

  // Path<2> path4([=](REAL t) -> Point<2> {
  //   REAL tt = 4*pi()*t;
  //   return Point<2>({tt*cos(tt), tt*sin(tt)});
  // });
  // path4.plot2D("t.png", 100, -4*pi(), 4*pi(), -4*pi(), 4*pi());

  // // cycloid
  // REAL r = 2, r0 = 3*r/2, v = 4*pi()*r;
  // Path<2> path5([=](REAL t) -> Point<2> {
  //   REAL theta = v*t/r;
  //   return Point<2>({v*t+r0*sin(theta), r0*cos(theta)});
  // });
  // path5.plot2D("t.png", 300, 0, v*1, -r0, r0);



  // Surface<2> surface1([=](REAL u, REAL v) -> Point<2> {
  //   return Point<2>({u, v});
  // });
  // surface1.plot2D("t.png", 50, -1, 2, 0, 2);

  // Surface<2> surface2([=](REAL u, REAL v) -> Point<2> {
  //   return Point<2>({u, u+v});
  // });
  // surface2.plot2D("t.png", 50, 0, 1, 0, 2);

  Surface<2> surface3([=](REAL u, REAL v) -> Point<2> {
    return Point<2>({2*pi()*u, sin(2*pi()*u)*v});
  });
  surface3.plot2D("t.png", 100, 0, 2*pi(), -1, 1);

  
  // Surface<2> surface4([=](REAL u, REAL v) -> Point<2> {
  //   REAL uu = 4*pi()*u;
  //   // REAL vuu = v*uu;
  //   return Point<2>({uu*cos(uu), uu*sin(uu)});
  // });
  // surface4.plot2D("t.png", 200, 0, 6*pi(), -1, 1);

  
  // Surface<2> surface5([=](REAL u, REAL v) -> Point<2> {
  //   return Point<2>({u, v*(u-0.25)*(u-0.5)*(u-0.75)});
  // });
  // surface5.plot2D("t.png", 150, 0, 1, RATIONAL(-3,32), RATIONAL(3,32));
}
