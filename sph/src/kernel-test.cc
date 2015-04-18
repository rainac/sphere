/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <iostream>
#include <unistd.h>
#include <math.h>


#ifndef SPH_DIM
#define SPH_DIM 3
#endif
#define TEST_CORRECTNESS

#ifdef TEST_CORRECTNESS
#include "yafad/yafad.hh"
typedef yafad::FO::Static::Active<3> adouble;
#endif

#include "kernels.hh"
#include "utility/getenv.hh"
#include "quader.hh"

#ifdef TEST_CORRECTNESS
typedef adouble DataType;
typedef Point<adouble, ndim> VectorAD;
#else
typedef double DataType;
#endif

typedef Point<DataType, ndim> VectorType;

typedef GetEnvV<ndim, Vector> GetEnvVec;

typedef KernelInterface<DataType, VectorType, SPH_DIM> Kernel;

static bool failed;

using namespace std;

static void testKernel(Kernel *kernel, double const H, 
		       size_t const Npoints, double const 
#if SPH_DIM > 2
                       Z
#endif
                       ,
		       double const tol) {

  VectorType v, gradw, gradw2, diffGrad;
  DataType w;
#ifdef TEST_CORRECTNESS
  DataType w2;
#endif

  Vector vul(ndim);
  Vector vor(ndim);

  vul = -H;
  vor = H;

#if SPH_DIM > 2
  vul[2] = Z;
  vor[2] = Z;
#endif

  Quader<Vector> area(GetEnvVec("SPH_PARAM_QUADER_P0", vul),
                      GetEnvVec("SPH_PARAM_QUADER_P1", vor));

  cout.precision(16);

  double const dx = H/Npoints;
  size_t const nx = round(area.v[0] / dx);

#if SPH_DIM > 1
  double const dy = H/Npoints;
  size_t const ny = round(area.v[1] / dy);
#endif

  for (size_t i = 0; i <= nx; ++i) {
#if SPH_DIM > 1
    for (size_t j = 0; j <= ny; ++j) {
#else
    for (size_t j = 0; j < 1; ++j) {
#endif
      v[0] = area.untenLinks[0] + i * dx;
#if SPH_DIM > 1
      v[1] = area.untenLinks[1] + j * dy;
#endif
#if SPH_DIM > 2
      v[2] = area.untenLinks[2];
#endif
#ifdef TEST_CORRECTNESS
      v[0].diff(0) = 1;
#if SPH_DIM > 1
      v[1].diff(1) = 1;
#endif
#if SPH_DIM > 2
      v[2].diff(2) = 1;
#endif
#endif
      DataType const dnorm = norm(v);
      w = kernel->wAndGrad(v, dnorm, gradw);
#ifdef TEST_CORRECTNESS
      w2 = kernel->w(dnorm);
      if (fabs(w - w2) / fabs(w2) > tol) {
	cerr << "error: kernel " << kernel->name() << ", " << ndim << "D: "
	  "functions wAndGradw and w return different results: "
	     << w << " " << w2 << "\n";
	failed = true;
      }
      kernel->gradw(v, dnorm, gradw2);
      if (norm(gradw - gradw2)/norm(gradw) > tol) {
	cerr << "error: kernel " << kernel->name() << ", " << ndim << "D: "
	  "functions wAndGradw and gradw return different results: "
	     << gradw << " " << gradw2 << "\n";
	failed = true;
      }
      VectorType adderiv;
      adderiv[0] = w.diff(0);
#if SPH_DIM > 1
      adderiv[1] = w.diff(1);
#endif
#if SPH_DIM > 2
      adderiv[2] = w.diff(2);
#endif
      if (norm(gradw - adderiv)/norm(gradw) > tol) {
	cerr << "error: kernel " << kernel->name() << ", " << ndim << "D: "
	  "functions gradw return is not true derivative of w: at "
	     << v << " " <<  gradw << " != " << adderiv << "\n";
	failed = true;
      }
      cout
	<< v[0] << " "
#if SPH_DIM > 1
	<< v[1] << " "
#else
	<< 0 << " "
#endif
#if SPH_DIM > 2
	<< v[2] << " "
#else
	<< 0 << " "
#endif
	<< w << " "
	<< gradw[0] << " "
#if SPH_DIM > 1
	<< gradw[1] << " "
#else
	<< 0 << " "
#endif
#if SPH_DIM > 2
	<< gradw[2] << " "
#else
	<< 0 << " "
#endif
	<< norm(gradw - gradw2)
	<< "\n";
#else
      cout
	<< v[0] << " "
#if SPH_DIM > 1
	<< v[1] << " "
#else
	<< 0 << " "
#endif
#if SPH_DIM > 2
	<< v[2] << " "
#else
	<< 0 << " "
#endif
	<< w << " "
	<< gradw[0] << " "
#if SPH_DIM > 1
	<< gradw[1] << " "
#else
	<< 0 << " "
#endif
#if SPH_DIM > 2
	<< gradw[2] << " "
#else
	<< 0 << " "
#endif
	<< "\n";
#endif
    }
  }
}

struct MyErrorHandler : public yafad::ErrorHandler {
  void handle(yafad::Error const &) {
//     std::cerr << "sph-ad: error: " << e.what() << "\n";
  }
};

static void showVersion() {
  std::cout << "kernel-test " << ndim << "D" 
	    << ": sample SPH Kernel on XY plane and output results" << std::endl;
}

int main(int const argc, char *argv[]) {

  MyErrorHandler myErrorHandler;
  yafad::ErrorHandler::setHandler(&myErrorHandler);

  double H(GetEnv("SPH_PARAM_H", 1, 1e-100, 1e100));
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
      cout << "usage: kernel-test [Options]" << endl;
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


  Kernel *kernel = Kernel::makeKernel(kernelName, H);

  if (kernel == 0) {
    failed = true;
  } else {
    testKernel(kernel, H, Npoints, Z, R);
  }

  return failed != false;
}
