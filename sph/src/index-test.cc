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

#include "sph.hh"

int curNumThreads = 2;

static bool failed;

using namespace std;

static void testKernel() {

  Vector vul(ndim);
  Vector vor(ndim);

  vul = 0;
  vor = 1;

  double const H = 0.1;

  Quader<Vector> area(vul, vor);

  PArray<Particle> partikel(8);
  partikel[0].position()[0] = 0.05;
  partikel[1].position()[0] = 0.05;
  partikel[2].position()[0] = 0.15;
  partikel[3].position()[0] = 0.15;
  partikel[4].position()[0] = 0.25;
  partikel[5].position()[0] = 0.25;
  partikel[6].position()[0] = 0.35;
  partikel[7].position()[0] = 0.35;

  /*
  partikel[0].position()[0] = 0.05;
  partikel[1].position()[0] = 0.15;
  partikel[2].position()[0] = 0.25;
  partikel[3].position()[0] = 0.35;
  partikel[4].position()[0] = 0.35;
  partikel[5].position()[0] = 0.25;
  partikel[6].position()[0] = 0.15;
  partikel[7].position()[0] = 0.05;
  */

  RasterInfo rasterInfo(area, H);
  HashedInfo hashedInfo(rasterInfo, 2);

  for (size_t i = 0; i < partikel.size(); ++i) {
    hashedInfo.updatePartikelIndex(partikel[i]);
    cout << "particle " << i << " " << partikel[i].position() << " " << partikel[i].index() 
         << "\n";
    
    if (i < 4) {
      assert((floor(partikel[i].position()*10.0) == partikel[i].index()).min() == 1);
    }
  }

  SortedPartikelIndex<SPHPartikelArray> 
    partikelIndex("counting_sort", hashedInfo.m_index, area, H);

  partikelIndex.update(partikel);
  
  BinTriplet infos[8];
  
  for (size_t i = 0; i < partikel.size(); ++i) {
    partikelIndex.getBinInfo(partikel[i].hashVal(), infos[i]);
    assert(infos[i].endBoundary - infos[i].offset == 2);
    assert(infos[i].endMoving - infos[i].offset == 2);
    if (i < 4) {
      assert(partikelIndex.getPartikelIndex(infos[i].offset) == i
             or partikelIndex.getPartikelIndex(infos[i].offset + 1) == i);
    }
    cout << "partikel index " << i << ": " << 
      partikelIndex.getPartikelIndex(i) << "\n";
  }

}

static void showVersion() {
  std::cout << "index-test " << ndim << "D: sample SPH Kernel on XY plane and output results" << std::endl;
}

int main(int const argc, char *argv[]) {

  // double H(GetEnv("SPH_PARAM_H", 1, 1e-100, 1e100));
  // size_t Npoints = GetEnv("SPH_PARAM_N", 10, 1e-100, 1e100);
  // double Z(GetEnv("SPH_PARAM_Z", 0, -1e100, 1e100));
  // double R(GetEnv("SPH_PARAM_R", 1e-13, -1e100, 1));

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
    // case 'H': {
    //   char *endptr;
    //   double htmp = strtod(optarg, &endptr);
    //   if (endptr == optarg) {
    //     cerr << "error: invalid float number as argument of option " << char(optopt) << endl;
    //     return 2;
    //   } else {
    //     H = htmp;
    //   }
    // }
    //   break;
    // case 'R': {
    //   char *endptr;
    //   double htmp = strtod(optarg, &endptr);
    //   if (endptr == optarg) {
    //     cerr << "error: invalid float number as argument of option " << char(optopt) << endl;
    //     return 2;
    //   } else {
    //     R = htmp;
    //   }
    // }
    //   break;
    // case 'N': {
    //   char *endptr;
    //   long ntmp = strtol(optarg, &endptr, 10);
    //   if (endptr == optarg) {
    //     cerr << "error: invalid int number as argument of option " << char(optopt) << endl;
    //     return 2;
    //   } else {
    //     Npoints = ntmp;
    //   }
    // }
    //   break;
    // case 'Z': {
    //   char *endptr;
    //   double htmp = strtod(optarg, &endptr);
    //   if (endptr == optarg) {
    //     cerr << "error: invalid float number as argument of option " << char(optopt) << endl;
    //     return 2;
    //   } else {
    //     Z = htmp;
    //   }
    // }
    //   break;
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


  testKernel();

  return failed != false;
}
