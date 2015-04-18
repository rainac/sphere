/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <iostream>
#include <unistd.h>
#include <math.h>

#undef SPH_RENDER

#ifndef SPH_DIM
#define SPH_DIM 3
#endif

#define TEST_CORRECTNESS

// here we use AD locally to test some components of the SPH code
// so we do not define SPH_AD

#ifdef TEST_CORRECTNESS
#include "yafad/yafad.hh"
typedef yafad::FO::Static::Active<3> adouble;
#endif

#include "sph.hh"
#include "utility/getenv.hh"
#include "quader.hh"

double maxRelativeVelocity;
double maxVelocity;

using namespace std;

static bool failed;

double const epsilon = 1e-14;
static bool nearEqual(double const v1, double const v2) {
  double const d = (v2 - v1);
  if (d/v1 > epsilon) {
    return false;
  }
  return true;
}

static void testKernel(double const , size_t const , double const ,
		       double const , 
                       int pclass1 = Partikel<double>::PARTICLE_TYPE_MONAGHAN,
                       int pclass2 = Partikel<double>::PARTICLE_TYPE_MONAGHAN,
                       bool const symmetric = true
                       ) {
  typedef KernelInterface<double, Point<double, ndim>, ndim> KernelInterfaceT;

  cout << "<testing class1='" << Partikel<double>::getParticleClassesName(pclass1)
       << "' class2='" << Partikel<double>::getParticleClassesName(pclass2) << "'/>\n";
  
  Partikel<double> p, o;
  PartikelVariablen<double> s1a, s1b, s2a, s2b, sd, snull;

  p._class() = pclass1;
  o._class() = pclass2;

  setenv("SPH_PARAM_PARTICLE_TYPE_MOVING_DEFAULT", 
         Partikel<double>::getParticleClassesName(pclass1), 1);

  if (o._class() == Partikel<double>::PARTICLE_TYPE_POINT_MIRRORED
      or o._class() == Partikel<double>::PARTICLE_TYPE_LENNARD_JONES) {
    o.flags() |= Partikel<double>::PARTICLE_FLAG_BOUNDARY;
    setenv("SPH_PARAM_PARTICLE_TYPE_BOUNDARY_DEFAULT", 
           Partikel<double>::getParticleClassesName(pclass2), 1);
  }

#ifdef SPH_PARTILE_ID
  o.id() = 1;
#endif

  SPHKonfig<double> konfig;
  PartikelManager<SPHKonfig<double>, KernelInterfaceT> pm(konfig, cout);
  konfig.writeXML(cout);

  unsigned short randSeed[3] = {};
  randSeed[0] = time(0);
  randSeed[1] = time(0);
  randSeed[2] = time(0);

  for(unsigned i = 0; i < 100; ++i) {
    for(unsigned k = 0; k < ndim; ++k) {
      p.position()[k] = erand48(randSeed)*konfig.kernelReachH/10;
      p.velocity()[k] = erand48(randSeed);
    }
    p.dichte() = konfig.eqsKonfig.rho0 + (erand48(randSeed) - 0.5)*konfig.eqsKonfig.rho0*1e-1;
    p.waerme() = erand48(randSeed);
    p.masse() = erand48(randSeed);
    o.dichte() = konfig.eqsKonfig.rho0 + (erand48(randSeed) - 0.5)*konfig.eqsKonfig.rho0*1e-1;
    o.waerme() = erand48(randSeed);
    p.masse() = erand48(randSeed);
    o.masse() = p.masse();
    pm.updatePressure(p);
    pm.updatePressure(o);
    assert(p.flags() == 0);
//     assert(o.flags() == 0);
    assert(p._class() == pclass1);
    assert(o._class() == pclass2);
#ifdef SPH_PARTILE_ID
    assert(p.id() == 0);
    assert(o.id() == 1);
#endif
    pm.particleRelation(p, o, s1a, s1b);
    if (not o.isBoundary()) {
      pm.particleRelation(o, p, s2a, s2b);
    }
    cout << "<p1 " << p << "/>\n<p2 " << o << "/>\n";
    cout << "<sum1 " << s1a << "/>\n<sum2 " << s1b << "/>\n";
    if (symmetric) {
      if (not o.isBoundary()) {
        assert(nearEqual(norm(s1a.position() + s1b.position()), 0));
        assert(nearEqual(norm(s1a.velocity() + s1b.velocity()), 0));
        assert(nearEqual(s1a.dichte(), s1b.dichte()));
        assert(nearEqual(s1a.waerme(), s1b.waerme()));

        assert(nearEqual(norm(s1a.position() + s2a.position()), 0));
        assert(nearEqual(norm(s1a.velocity() + s2a.velocity()), 0));
        assert(nearEqual(s1a.dichte(), s2a.dichte()));
        assert(nearEqual(s1a.waerme(), s2a.waerme()));
      } else {
        assert(s1b == snull);
      }
    } else {
      if (not o.isBoundary()) {
        cout << "<sum1-alt " << s2a << "/>\n<sum2-alt " << s2b << "/>\n";
        assert(nearEqual(norm(s1a.position() - s2b.position()), 0));
        assert(nearEqual(norm(s1a.velocity() - s2b.velocity()), 0));
        assert(nearEqual(s1a.dichte(), s2b.dichte()));
        assert(nearEqual(s1a.waerme(), s2b.waerme()));
        
        assert(nearEqual(norm(s1b.position() - s2a.position()), 0));
        assert(nearEqual(norm(s1b.velocity() - s2a.velocity()), 0));
        assert(nearEqual(s1b.dichte(), s2a.dichte()));
        assert(nearEqual(s1b.waerme(), s2a.waerme()));
      } else {
        assert(s1b == snull);
        assert(symmetric);
        assert(0); // all boundary types symmetric yet
      }
    }
//     pm.computeParticleRelation(p, o, s2);
//     assert(s1 == s2);
  }

}

static void showVersion() {
  std::cout << "kernel-test " << ndim << "D: sample SPH Kernel on XY plane and output results" << std::endl;
}

int main(int const argc, char *argv[]) {

  double H(GetEnv("SPH_PARAM_H", 0.01, 1e-100, 1e100));
  size_t Npoints = GetEnv("SPH_PARAM_N", 10, 1e-100, 1e100);
  double Z(GetEnv("SPH_PARAM_Z", 0, -1e100, 1e100));
  double R(GetEnv("SPH_PARAM_R", 1e-13, -1e100, 1));

  std::string kernelName = "spline";

  char getoptResult = 0;
  char availableOptions[20] = "hvH:N:R:K:Z:";
  do {
    getoptResult = getopt(argc, argv, availableOptions);
    switch(getoptResult) {
    case '?':
      return 2;
    case ':':
      cerr << "error: argument missing of option: " << char(optopt) << endl;
      break;
    case 'h':
      showVersion();
      cout << "usage: partrel-test [Options]" << endl;
      cout << "  -h                       show help and exit" << endl;
      cout << "  -v                       show version and exit" << endl;
      cout << "  -H                       set kernel influence lenggth H" << endl;
      cout << "  -N                       set number of points per dimension" << endl;
      cout << "  -Z                       set Z coordinate (default 0) " << endl;
      cout << "  -K                       select kernel type" << endl;
      cout << "  -L                       list available kernel types" << endl;
      return 0;
    case 'v':
      showVersion();
      return 0;
    case 'L': {
      cout << "gauss\n"
	"spline\n"
	"wendland0\n"
	"wendland1\n"
	"wendland2\n"
	"wendland3\n";
      return 0;
    }
      break;
    case 'H': {
      char *endptr;
      double htmp = strtod(optarg, &endptr);
      if (endptr == optarg) {
	cerr << "error: invalid float number as argument of option " << char(optopt) << endl;
	return 2;
      } else {
	H = htmp;
      }
    }
      break;
    case 'R': {
      char *endptr;
      double htmp = strtod(optarg, &endptr);
      if (endptr == optarg) {
	cerr << "error: invalid float number as argument of option " << char(optopt) << endl;
	return 2;
      } else {
	R = htmp;
      }
    }
      break;
    case 'N': {
      char *endptr;
      long ntmp = strtol(optarg, &endptr, 10);
      if (endptr == optarg) {
	cerr << "error: invalid int number as argument of option " << char(optopt) << endl;
	return 2;
      } else {
	Npoints = ntmp;
      }
    }
      break;
    case 'Z': {
      char *endptr;
      double htmp = strtod(optarg, &endptr);
      if (endptr == optarg) {
	cerr << "error: invalid float number as argument of option " << char(optopt) << endl;
	return 2;
      } else {
	Z = htmp;
      }
    }
      break;
    case 'K': {
      kernelName = optarg;
    }
      break;
    case -1:
      break;
    default:
      cerr << "error: unknown return value from getopt: " << int(getoptResult) << endl;
      break;
    }
  } while (getoptResult != -1);

  
  ostringstream hstr;
  hstr << H;
  setenv("SPH_PARAM_H", hstr.str().c_str(), 1);

  cout << "<partrel-test>\n";
  testKernel(H, Npoints, Z, R, 
             Partikel<double>::PARTICLE_TYPE_MONAGHAN, 
             Partikel<double>::PARTICLE_TYPE_MONAGHAN);

  testKernel(H, Npoints, Z, R, 
             Partikel<double>::PARTICLE_TYPE_MONAGHAN_XSPH, 
             Partikel<double>::PARTICLE_TYPE_MONAGHAN_XSPH);

  testKernel(H, Npoints, Z, R, 
             Partikel<double>::PARTICLE_TYPE_FERRARI, 
             Partikel<double>::PARTICLE_TYPE_FERRARI, false);

  testKernel(H, Npoints, Z, R, 
             Partikel<double>::PARTICLE_TYPE_MONAGHAN, 
             Partikel<double>::PARTICLE_TYPE_LENNARD_JONES);

  testKernel(H, Npoints, Z, R, 
             Partikel<double>::PARTICLE_TYPE_MONAGHAN, 
             Partikel<double>::PARTICLE_TYPE_POINT_MIRRORED);

  testKernel(H, Npoints, Z, R, 
             Partikel<double>::PARTICLE_TYPE_FERRARI, 
             Partikel<double>::PARTICLE_TYPE_LENNARD_JONES);

  testKernel(H, Npoints, Z, R, 
             Partikel<double>::PARTICLE_TYPE_FERRARI, 
             Partikel<double>::PARTICLE_TYPE_POINT_MIRRORED);
  
  cout << "<status failed='" << failed << "'/>\n";
  cout << "</partrel-test>\n";

  return failed != false;
}
