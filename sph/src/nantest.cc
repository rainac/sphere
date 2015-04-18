/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <signal.h>
#include <valarray>
#include <iostream>
#include <fenv.h>

 void signalHandler(int) {
   abort();
 }

int main() {

  signal(SIGINT, signalHandler);

  signal(SIGFPE, signalHandler);

  signal(SIGTRAP, signalHandler);

  {
    int excepts = FE_ALL_EXCEPT;
    fexcept_t flagp;
    fegetexceptflag(&flagp, excepts);
    excepts |= FE_ALL_EXCEPT;
    fesetexceptflag(&flagp, excepts);
    
//     fenv_t myfenv;
//     fegetenv(&myfenv);
//     //   myfenv.trapstate=1;
//     fesetenv(&myfenv);
//     //    fp_enable_all();
  }

  double x = 1, y = 2, z = 3;
  std::valarray<double> X(100), Y(100);

  X = x;
  Y = y;

  Y[10] = nan("abd");

  int a = 1000000;
  int b = 1034240;
//   int c = a*b;

//   std::cerr << "c: " << c << "\n";

  x = X.sum();
  y = Y.sum();

  if (x < y) {
    std::cerr << "x is less\n";
  }

  if (x > y) {
    std::cerr << "x is greater\n";
  }

  if (isnan(y)) {
    std::cerr << "y is NaN\n";
  }

  std::cerr << "x: " << x << "\n";
  std::cerr << "y: " << y << "\n";

  std::cerr << "atan(2) " << std::atan(2.) << "\n";
  std::cerr << "log(-2) " << std::log(-2.) << "\n";

}
