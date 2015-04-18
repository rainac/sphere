/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/
#include <omp.h>

#include "counter.hh"
#include "index.hh"
#include "spline/point.hh"
// #include "sph.hh"

typedef Point<long, 3> Coord;

size_t const ndim = 3;

using namespace std;

int main() {

  size_t const konfig_m = 5;

  Counter<Coord, 1> zaehlerGrid;
  zaehlerGrid.top = 1 << (konfig_m - 1);

  Index<long, ndim> gridIndex(zaehlerGrid.top);
  HashedIndex<long, ndim> gridIndex2(zaehlerGrid.top, konfig_m - 1);
  
  cerr << "parallel grid enumeration test\n"
       << "m=" << konfig_m << ""
       << "\n";

  omp_set_num_threads(8);

#pragma omp parallel for
  for(long i = 0; i < static_power<ndim>(1 << (konfig_m - 1)); ++i) {
    
//     Coord c1 = zaehlerGrid.v;
    Coord c2 = gridIndex.invert(i);
    Coord c3 = gridIndex2.invert(i);

//     c1 <<= 1;
    c2 <<= 1;
    c3 <<= 1;

#pragma omp critical
    {
      cerr << i << ": " << c2 << " " << c3 << "\n";
//       ++zaehlerGrid;
    }

//     assert((c1 == c2).min() == 1);
//     assert((c1 == c3).min() == 1);
    assert((c2 == c3).min() == 1);

  }
  
  cerr << "fertig\n";
}
