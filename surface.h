#include <array>

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "iRRAM.h"

using namespace iRRAM;

#include "euclidean.h"
#include "plot.h"

#define PLOT_COLOR_R    0x00
#define PLOT_COLOR_G    0x00
#define PLOT_COLOR_B    0x00

// [0,1]^2 -> R^N
template <int N>
class Surface {
public:
  // original function
  std::function<Point<N>(REAL,REAL)> f;

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
  }
  Surface(std::function<Point<N>(REAL,REAL)> f) { this->f = f; }
  
  // membership test for point with precision 2^-p
  // in accordance to the Ko compatibility
  bool member(Point<N> point, int p) {
    single_valued code;

    // check if previously found pArg is viable
    // if not, increase the precision
    if(this->p < p) this->increasePrecision(p);

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
    cout << "point = " << point[0] << " " << point[1] << "\n";
    RATIONAL u,v;
    for(u=step/2 ; u<=ONE ; u+=step) {
      for(v=step/2 ; v<=ONE ; v+=step) {
        d = IR_d<N>(point, this->f(u, v));
        if(choose(d < radius, d > radiusHalf) == 1) return true;
      }
    }
    return false;
  }

  // save the plane graph to an image file
  // area to draw: [x1, x2] X [y1, y2]
  // Set image width. Height will be determined automatically.
  // REQUIRE: x1 < x2, y1 < y2
  void plot2D(const char *filename, int width, REAL x1, REAL x2, REAL y1, REAL y2) {
    // only plane
    if(N != 2) return;

    // some values..
    REAL pixelSize = (x2-x1)/REAL(width);      // single pixel size as a rect
    int height = ceil(((y2-y1)/pixelSize).as_double());               // image height

    // precision
    // Define p such that a ball centered at the center of pixel cover the pixel
    // That is, pixelSize/2*sqrt(2)  <  2^-p (radius of ball)
    // The drawn path will not be broken.
    int p = floor((REAL(0.5) - log(pixelSize)/ln2()).as_double());       // precision

    cout << "height = " <<  height << ", p = " << p << "\n";

    // plot each pixel to palette
    // start at the bottom row, from left to right.. row += 1 .. repeat
    // variable point stores the coordinate of the center of the current pixel
    Point<N> point = {REAL(0), y1 - pixelSize/2};      // init with coordinate outside image
    Palette pal(width, height);
    for(int i=height-1 ; i>=0 ; i--) {
      point[0] = x1 + pixelSize/2;    // x; the left most pixel
      point[1] += pixelSize;          // y; +1 row

      for(int j=0 ; j<width ; j++) {
        cout << i << " " << j << "\n";
        if(this->member(point, p)) {
          pal.setColor(j, i, PLOT_COLOR_R, PLOT_COLOR_G, PLOT_COLOR_B);
        }
        point[0] += pixelSize;
      }

      cout << (height-i) << " / " << height << " row done\n";
    }

    // to image
    writeImage(filename, pal);
  }
};

