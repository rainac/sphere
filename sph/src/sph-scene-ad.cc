/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include "spline/point.hh"
#define YAFAD_POINT_DEFINED

#include "yafad/yafad.hh"

#define  SPH_AD_NDIR 12

typedef yafad::FO::Static::Active<SPH_AD_NDIR> adouble;

#define  SPH_AD 1
#define  SPH_ACTIVE_DATA_TYPE adouble

#include "sph-scene.cc"


