/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <iostream>
extern "C" char **environ;
#include "lennard-jones.hh"
#include "spline/spline.hh"
#include "utility/getenv.hh"
#include "vector.hh"
#include "quader.hh"

using namespace std;

typedef GetEnvV<ndim, Vector> GetEnvVec;


int main() {
  Vector v;

  double const h(GetEnv("SPH_PARAM_H", 1e-3, 1e-100, 1));

  double const dx = h/10, dy = h/10;
  
  Vector vul(ndim);
  Vector vor(ndim);

  vul = -2*h;
  vor = 2*h;

//   vul[2] = -h/3;
//   vor[2] = -h/3;

  vul[2] = 0;
  vor[2] = 0;

  Quader<Vector> area(GetEnvVec("SPH_PARAM_QUADER_P0", vul),
                      GetEnvVec("SPH_PARAM_QUADER_P1", vor));

  LennardJones<double> splineKernel(GetEnv("SPH_PARAM_LJ_R0", h, 1e-100, 1e100),
                                    GetEnv("SPH_PARAM_LJ_D", 1, 1e-100, 1e100),
                                    GetEnv("SPH_PARAM_LJ_P1", 2, 2, 1e100));

  cout.precision(16);

  size_t const nx = round(area.v[0] / dx);
  size_t const ny = round(area.v[1] / dy);

  for (size_t i = 0; i <= nx; ++i) {
    for (size_t j = 0; j <= ny; ++j) {
      v[0] = area.untenLinks[0] + i * dx;
      v[1] = area.untenLinks[1] + j * dy;
      v[2] = area.untenLinks[2];
      cout << v << " " << splineKernel.w(norm(v)) * v / normSquared(v) 
	   << " " << 0 << " " << 0 << " " << 0 
	   << "\n";
    }
    cout << "\n";
  }

}
