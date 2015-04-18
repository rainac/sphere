/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <iostream>
#include "permutation.hh"

int main() {
  Permutation<unsigned> p(26);

  std::string str("abcdefghijklmnopqrstuvwxyz");

//   p[0] = 5;
//   p[5] = 0;

  p[0] = 25;
  p[25] = 0;

//   p[0] = 1;
//   p[1] = 2;
//   p[2] = 3;
//   p[3] = 0;

//   p[6] = 8;
//   p[8] = 6;

  assert(p.valid());

//   p(str.begin(), str.end());

  std::cerr << "result: \"" << str << "\"\n";

}
