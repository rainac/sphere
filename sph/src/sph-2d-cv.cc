/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <complex>

typedef  std::complex<double> adouble;

template<class T>
inline double yafad_value(T const &x) {
  return x;
}
template<>
inline double yafad_value<std::complex<double> >(std::complex<double> const &x) {
  return x.real();
}

#define  SPH_AD 2
#define  SPH_ACTIVE_DATA_TYPE adouble
#define  SPH_DIM 2

#define  SPH_AD_NDIR 1

namespace std {
  inline bool operator <(adouble const &a, adouble const &b) {
    return a.real() < b.real();
  }
  inline bool operator >(adouble const &a, adouble const &b) {
    return a.real() > b.real();
  }
  inline bool operator >=(adouble const &a, adouble const &b) {
    return a.real() >= b.real();
  }
  inline bool operator <=(adouble const &a, adouble const &b) {
    return a.real() <= b.real();
  }
  inline adouble max(adouble const &a, adouble const &b) {
    if (b > a) return b;
    return a;
  }
  inline adouble min(adouble const &a, adouble const &b) {
    if (b < a) return b;
    return a;
  }
}

static double const cvMethEpsilon = 1e-200;

#include "sph.cc"
