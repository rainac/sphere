#ifndef sph_vector_hh
#define sph_vector_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <valarray>
#include "common.hh"

#include "spline/point.hh"

//*! Use Statically sized Point Template
typedef Point<SphADataType, ndim> Vector;

//*! Use Statically sized Point Template
typedef Point<double, ndim> DoubleVector;

#ifdef SPH_AD
DoubleVector doubleVector(Vector const &v) {
  DoubleVector res;
  for(unsigned i = 0; i < ndim; ++i) {
    res[i] = yafad_value(v[i]);
  }
  return res;
}
#else
#define doubleVector(x) (x)
#endif



template<class T>
static T normSquared(std::valarray<T> const &v) {
  T res = T();
  for(size_t i = 0; i < v.size(); ++i) {
    T const vi = v[i];
    res += vi*vi;
  }
  return res;
}
template<class T>
static T norm(std::valarray<T> const &v) {
  return sqrt(normSquared(v));
}


/// The integer type to hold an index into the \f$ h \f$-spaced 
/// spatial raster or the max. number of bins in the raster.
typedef unsigned long CoordinateIndex;

/// The integer type to hold a spatial integer index, \f$ x/h \f$, which
/// is the index into the \f$ h \f$-spaced spatial raster.
typedef long CoordType;

/// Define the type of a 1,2,3-dimensional integer coordinate used for
/// the raster index operations.
typedef Point<CoordType, ndim> Coord;

/// Define the type to hold a 3-dimensional RGB-color \f$c\f$ used with Open-GL,
/// \f$ c \in [0, 1]^3 \f$.
typedef Point<double, 3> Color;

#endif
