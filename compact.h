#pragma once

#include "euclidean.h"
#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "iRRAM.h"
#include "plot.h"

using namespace iRRAM;


#define PLOT_COLOR_R    0x00
#define PLOT_COLOR_G    0x00
#define PLOT_COLOR_B    0x00


// R^N
template <int N>
class Compact {
public:
  // characteristic function
  std::function< bool (Point<N> , int) > cfun;

  // init with the emptyset
  Compact() { this->cfun = [=](Point<N>, int) -> bool { return false; }; }

  // init with characteristic func
  Compact(std::function<bool(Point<N>, int)> cfun) { this->cfun = cfun; }

  // membership test for point with precision 2^-p
  virtual bool member(Point<N> point, int p) { return this->cfun(point, p); }

  
  // save the 2D graph to an .png file
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

    // plot each pixel to palette
    // start at the bottom row, from left to right.. row += 1 .. repeat
    // variable point stores the coordinate of the center of the current pixel
    Point<N> point = {REAL(0), y1 - pixelSize/2};      // init with coordinate outside image
    Palette pal(width, height);
    for(int i=height-1 ; i>=0 ; i--) {
      point[0] = x1 + pixelSize/2;    // x; the left most pixel
      point[1] += pixelSize;          // y; +1 row

      for(int j=0 ; j<width ; j++) {
        if(member(point, p)) {
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



  
// pointwise conjunction for binary Boolean functions
template <int N>
Compact<N> conjunction(Compact<N> &com1, Compact<N> &com2) {
  auto cfun = [&] (Point<N> pt, int p) -> bool {
    return com1.member(pt, p) && com2.member(pt, p);
  };
  return Compact<N>(cfun);
}
  
// pointwise conjunction for binary Boolean functions
template <int N>
Compact<N> disjunction(Compact<N> &com1, Compact<N> &com2) {
  auto cfun = [&] (Point<N> pt, int p) -> bool {
    return com1.member(pt, p) || com2.member(pt, p);
  };
  return Compact<N>(cfun);
}





// // The set of compact subsets of 1-dimensional real line
// class OneDCompact
// {
// public:
//   // characteristic function
//   std::function< bool ( REAL, int ) > cfun;

//   // create the emptyset
//   OneDCompact();

//   // create a singleton
//   OneDCompact(REAL);  

//   // create a closed interval
//   OneDCompact(REAL, REAL);  
// };


// // functionals


// // make the fuzzy char function for the emptyset
// std::function<bool (REAL, int) > empty ()
// {
//   return
//     [=](REAL y, int p) -> bool {
// 	return false;
//     };
// }


// // make the fuzzy char function for the singleton {x}
// std::function<bool (REAL, int) > singleton (REAL x)
// {
//   return
//     [=](REAL y, int p) -> bool {
//       REAL prec = scale(REAL(1), p);
//       if (choose(iRRAM::abs(x - y) < prec, iRRAM::abs(x-y) > 0) == 1)	   
// 	return true;
//       else
// 	return false;
//     };
// }



// // make the fuzzy char function for the interval [x, y]
// std::function<bool (REAL, int) > interval (REAL x, REAL y)
// {
//   return
//     [=](REAL z, int p) -> bool {
// 	   REAL prec = scale(REAL(1), p); 
// 	   if (choose(
// 		      z < y+prec && z > x - prec,
// 		      z < x || z > y) == 1)
// 	     return true;
// 	   else
// 	     return false;
// 	 };
// }


// OneDCompact::OneDCompact()
// {
//   cfun = empty(); 
// }

// OneDCompact::OneDCompact(REAL x)
// {
//   cfun = singleton(x); 
// }
  
// OneDCompact::OneDCompact(REAL x, REAL y)
// {
//   cfun = interval(x, y); 
// }




// // Primitive Operations: 
// OneDCompact op_intersection(OneDCompact x, OneDCompact y)
// {
//   OneDCompact z;
//   z.cfun = pw_conjunction<REAL, int>(x.cfun, y.cfun);
//   return z;
// }

// OneDCompact op_union(OneDCompact x, OneDCompact y)
// {
//   OneDCompact z;
//   z.cfun = pw_disjunction(x.cfun, y.cfun);
//   return z;

// }

// bool is_member_of(REAL x, OneDCompact y, int p)
// {
//   return y.cfun(x, p);
// }









// // main function for testing
// void compute()
// {

//   // closed interval [1, pi]
//   OneDCompact x(1, pi());


//   // test membership with precision -100
//   if (x.cfun(3, -100)) cout << "3 is in [1, pi]\n"; else cout << "3 is not in [1, pi]\n";  

//   if (x.cfun(4, -100)) cout << "4 is in [1, pi]\n"; else cout << "4 is not in [1, pi]\n";  


//   // closed interval [3, 4]
//   OneDCompact y(3,4);


//   // test membership with precision -100
//   if (y.cfun(3.2, -100)) cout << "3.2 is in [3, 4]\n"; else cout << "3.2 is not in [3, 4]\n";  

//   if (y.cfun(0, -100)) cout << "0 is in [3, 4]\n"; else cout << "0 is not in [3, 4]\n";  

  
  
// }
