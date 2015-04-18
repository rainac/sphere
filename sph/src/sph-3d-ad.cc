/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <cmath>

/// make the ad version use the same vector type, both for
/// ndim-Vectors and derivative vectors
///
#include "spline/point.hh"
#define YAFAD_POINT_DEFINED

#include "yafad/yafad.hh"

//template class yafad::FO::Static::Active<double, 3>;

#define  SPH_AD_NDIR 3

typedef yafad::FO::Static::Active<SPH_AD_NDIR> adouble;

#define  SPH_AD 1
#define  SPH_ACTIVE_DATA_TYPE adouble
#define  SPH_DIM 3

#include "sph.cc"
