#ifndef SPH_POWER_hh
#define SPH_POWER_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

/// This template class computes powers of constant exponents at
/// compile time. \see StaticPower<1>
template<int N> struct StaticPower {
  template<class T>
  static T value(T const v) {
    return v * StaticPower<N-1>::value(v);
  }
};

/// The base case: the partial specialization of StaticPower<int N> for N=1
/// returns the argument. \see StaticPower
template<> struct StaticPower<1> {
  template<class T>
  static T value(T const v) {
    return v;
  }
};

/// A further base case: the partial specialization of StaticPower<int N> for N=0
/// returns 1. \see StaticPower
template<> struct StaticPower<0> {
  template<class T>
  static T value(T) {
    return 1;
  }
};

template<int N, class T> 
inline T static_power(T const v) {
  StaticPower<N> p0;
  return p0.value(v);
}

template<class T>
static T dynamic_power(T const v, long i) {
  T result = 1, vp = v;
  while(i) {
    if (i & 1) {
      result *= vp;
    }
    vp *= vp;
    i >>= 1;
  }
  return result;
}

#endif
