/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <iostream>

using namespace std;

int main() {
  for (int dim = 1; dim < 32; ++dim) {
    size_t numCombos = 0, numHits = 0;
    for (int j = 1; j < (1<<dim); ++j) {
      for (int k = 0; k < j; ++k) {
        ++numCombos;
        if ((j & k) == 0) {
          ++numHits;
        }
      }
    }
    cout << "dim: " << dim << " combos: " << numCombos << " hits: " << numHits << " ratio: " << 
      double(numHits) / numCombos << "\n";
  }
}
