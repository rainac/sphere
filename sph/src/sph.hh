/*

This file is part of Sphere.

Sphere is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License (LPGL) as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

Sphere is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public
License along with Sphere.  If not, see <http://www.gnu.org/licenses/>.

Copyright © 2008,2009,2010 Johannes Willkomm

*/

/** \mainpage
  \section Introduction

  This SHP code  is a parallel implementation of the SPH method for fluid
  simulation. The characteristics are:

       - a fixed kernel smoothing length \f$ h \f$, leading to a fixed
         spatial raster length \f$ H \f$.

       - a special parallel algorithm to compute and sum up the
         forces. Each particle-particle force is computed only once
         yet the scalability is good.

  \section doc Reading the documentation

  To begin, start with the "classes" and "files" section that can be
  found above. To understand the particle finding data structure read
  the documentation of class SortedPartikelIndex. For the
  parallelization over the spatial grid, read the documentation of
  class PartikelManager::ForceFunctorSym.

  \section Usage

  The compiled program is self contained, except for the necessary
  third-party libraries, free-glut if rendering is enabled and HDF5 if
  saving is enabled. The H5Part library used is linked in statically,
  thus it is not needed at runtime.

  \section about About/Credits

  The SPH particle based method is implemented here as described in
  the well-known paper "Simulating free surface flows with SPH" by
  J.J. Monaghan.

  This software was written by Johannes Willkomm. The sorted particle
  index class and also the spatial parallelization strategy are my
  invention; at least I hope so, not withstanding prior art unknown to
  me.

  Any questions, bug-reports and comments should be sent to
  johannes@johannes-willkomm.de.

  The software was written as the reference implementation of the
  assignment of the software lab course "Parallele Programmierung mit
  OpenMP" held here at RWTH Aachen University in winter semester
  2008/2009.  Thanks go to the team that managed the lab course,
  Michael Luelfesmann, Arno Rasch, Andreas Wolf and Dr. Martin
  Buecker, for the numerous discussions and inspiration on how to
  understand the Monaghan paper.

 \todo Daempfung einbauen wie in Monaghan
 \todo Ereignisse: (Rand-)Partikel animieren

 */


/** \file
  This file holds the most important implementation details of SPH.

 */

#include <fstream>
#include <valarray>
#include <list>
#include <vector>
#include <set>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <omp.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef SPH_RENDER
#include <GL/freeglut.h>
# include <GL/gl.h>
#endif
#ifdef SPH_SAVE_H5PART
# include <H5Part.h>
#endif

// #include "sorting/countingsort.hh"
#include "counting-sort.hh"

#ifndef SPH_AD
#define yafad_value(x) (x)
#endif

#include "utility/fl-info.hh"
#include "utility/getenv.hh"
#include "utility/teestream.hh"
#include "utility/timer.hh"
#include "common.hh"
#include "kernels.hh"
#include "omp-util.hh"
#include "pthread-util.hh"
#include "lennard-jones.hh"
#include "index.hh"
#include "counter.hh"
#include "indirectsort.hh"
#include "partikel.hh"
#ifndef SPH_USE_VALARRAY_FOR_PARTICLE_ARRAYS
#  include "parray.hh"
#endif
#include "vtk-reader.hh"
#include "h5part-writer.hh"
#include "vtk-writer.hh"
#include "sph-config.hh"
#include "quader.hh"
#include "vector.hh"
#include "sensor.hh"
#include "sorted-particle-index.hh"


//! Define type for 3-dim vectors used in OpenGL rendering
typedef Point<double, 3>     Point3D;

#ifdef SPH_RENDER
static Point3D cross(Point3D const &a, Point3D const &b) {
  double res[3];
  res[0] = a[1]*b[2] - a[2]*b[1];
  res[1] = a[2]*b[0] - a[0]*b[2];
  res[2] = a[0]*b[1] - a[1]*b[0];
  return Point3D(res, 3);
}
#endif

/// The integer type to hold an index into a list of all particles, or
/// the max. number of particles.
typedef unsigned        PartikelIndex;


/// The type to get a ndim-dimensional value from the environment.
typedef GetEnvV<ndim, Vector> GetEnvVec;

/// The type to get a 3-dimensional value from the environment (used for renderer config).
typedef GetEnvV<3, Point<double, 3> > GetEnvVec3D;

#ifdef SPH_RENDER
#include "renderer.hh"
#endif

#if defined DEBUG_SPH
extern TSOStream<8> debugOutput;
#endif

#ifdef DEBUG_SPH
/// functional define DEB_SPH wrappes an expression to be output to
/// debugOut via operator <<. It allows to compile messages in and
/// out.
#define DEB_SPH(x) debugOutput << "sph: " << x << "\n";
#else
/// functional define DEB_SPH is defined to nothing, thereby throwing
/// away debugging messages.
#define DEB_SPH(x)
#endif


using namespace std;


/// Helper function that delete-s the data ( a single object) and sets
/// the pointer to zero.
//
/// \param ptr reference to a pointer variable that is zeroed after
/// delete-ing the data pointed to.
//
/// \tparam T the data type pointed to
template<class T>
inline void deleteAndClear(T * &ptr) {
  if (ptr) {
    delete ptr;
    ptr = 0;
  }
}

/// Helper function that resizes an array saving its content an then
/// copies the old values to the resized array's lower half.
//
/// \param v the array object (e.g. an std::valarray<of-something>)
/// \param newsz the desired new size
/// \tparam T1 the array type to resize
template<class T1>
void nonDestructiveResize(T1 &v, size_t const newsz) {
  T1 tmp(v);
  v.resize(newsz);
  for(size_t i = 0; i < tmp.size(); ++i) {
    v[i] = tmp[i];
  }
}

/// Helper function that resizes an array if the current size()
/// differs from newsz.  Calls resize(newsz) if v.size() != newsz.
/// Assumes that resize() creates clear data items, so if resize() is not called,
/// v = typename T1::value_type() is executed.
//
/// \param v the array object (e.g. an std::valarray<of-something>)
/// \param newsz the desired new size
/// \tparam T1 the array type to resize
template<class T1>
void resizeIfNeededButAlwaysClear(T1 &v, size_t const newsz) {
  if (v.size() != newsz) {
    v.resize(newsz);
  } else {
    v = typename T1::value_type();
  }
}

#ifdef SPH_VECTOR_STATIC
/// Wrapper that simply calls v.Point<T,NDIM>::writeXML().
/// Only present if SPH_VECTOR_STATIC defined.
//
/// \param aus reference to std::ostream object to write to.
/// \param v The point object to write out.
/// \tparam T Point item data type.
/// \tparam NDIM the number of dimensions.
template<class T, size_t NDIM>
static void writeXML(std::ostream &aus, Point<T,NDIM> const &v) {
  v.writeXML(aus);
}
#endif

/// Wrapper around the expression &v[i] for the non-const case.
//
/// \param v reference to an array object to data address of.
/// \param i The index of the item to get a pointer to.
/// \returns Pointer to the i-th data item.
/// \tparam T array item type.
template<class T>
static typename T::value_type *address(T &v, size_t i) {
  typename T::value_type &firstItem = v[0];
  typename T::value_type * baseAdress = &firstItem;
  return baseAdress + i;
}
#ifdef SPH_COMPILER_SUN_STUDIO
/// Wrapper around the expression &v[i] for the const case.  Only
/// present if SPH_COMPILER_SUN_STUDIO is defined. This version
/// requires a const cast, since Sun's compiler simply will not return
/// a reference but a copy in the const version of
/// valarray::operator[]
//
/// \param aus Reference to an array object to data address of.
/// \param v The index of the item to get a pointer to.
/// \returns Pointer to the i-th data item.
/// \tparam T array item type.
template<class T>
inline typename T::value_type const *address(T const &v, size_t i) {
  T &w = const_cast<T &>(v);
  typename T::value_type &firstItem = w[0];
  typename T::value_type const * baseAdress = &firstItem;
  return baseAdress + i;
}
#else
/// Wrapper around the expression &v[i] for the const case.
/// Only present if SPH_COMPILER_SUN_STUDIO is not defined.
//
/// \param v Reference to an array object to data address of.
/// \param i The index of the item to get a pointer to.
/// \returns Pointer to the i-th data item.
/// \tparam T array item type.
template<class T>
inline typename T::value_type const *address(T const &v, size_t i) {
  typename T::value_type const &firstItem = v[0];
  typename T::value_type const * baseAdress = &firstItem;
  return baseAdress + i;
}
#endif

/// Helper function that returns v.second - v.first.
//
/// \param v A std::pair object reference.
/// \return v.second - v.first.
/// \tparam T the type of v.
template<class S, class T>
static inline size_t length(std::pair<S,T> const &v) {
  return v.second - v.first;
}



#ifdef SPH_RENDER
  struct RenderOptions {
    bool colorRaster : 1
      , showRaster : 1
      , showParticleVectors : 1
      , hideBoundaryParticles : 1
      , showOutParticles : 1
      , visible : 1
      , hideCoordAxes : 1
      , hideSimulationArea : 1
      , showSensors : 1
      , colorRangeAll : 1
      ;
    int vectorLengthExp : 8;
    int particleSizeExp : 8;
    double vectorLength;
    double particleSize;
    RenderOptions() {
      memset(this, 0, sizeof(*this));
    }
  };
#endif

#ifdef SPH_RENDER_FORCES
/// A helper struct that represents the force betweeen two particles.
struct PartikelForce {
  typedef Partikel<SphADataType> particle_type;

  particle_type const *p1, *p2;
  Vector force;
  PartikelForce(particle_type const *p1, particle_type const *p2, Vector const &force) :
    p1(p1),
    p2(p2),
    force(force)
  { }
#ifdef SPH_RENDER
  void render() {
    glPushMatrix();
    {
      Vector const pos(p1->position());
#if SPH_DIM == 1
      glTranslated(_value(pos[0]), 0, 0);
#elif SPH_DIM == 2
      glTranslated(_value(pos[0]), _value(pos[1]), 0);
#else
      glTranslated(_value(pos[0]), _value(pos[1]), _value(pos[2]));
#endif

      glBegin(GL_LINES);
      {
	glVertex3d(0, 0, 0);
	Vector vnormed = force;
// 	vnormed *= SPH_config_h;

#if SPH_DIM == 1
	glVertex3d(_value(vnormed[0]), 0, 0);
#elif SPH_DIM == 2
        glVertex3d(_value(vnormed[0]), _value(vnormed[1]), 0);
#else
        glVertex3d(_value(vnormed[0]), _value(vnormed[1]), _value(vnormed[2]));
#endif
      }
      glEnd();

    }
    glPopMatrix();
  }
#endif
};
#endif

typedef Partikel<SphADataType>                    Particle;
typedef PartikelVariablen<SphADataType>   ParticleVariable;

#ifdef SPH_USE_VALARRAY_FOR_PARTICLE_ARRAYS

typedef std::valarray<Partikel> SPHPartikelArray;
typedef std::valarray<PartikelVariablen> SPHPartikelVariablenArray;

#else

/// use own parallel array
typedef PArray<Particle>                  SPHPartikelArray;
/// use own parallel array
typedef PArray<ParticleVariable> SPHPartikelVariablenArray;

#endif

// /// shortcut BinTriplet member type to a global name
// typedef SortedPartikelIndex<SPHPartikelArray>::BinTriplet BinTriplet;

/// helper function returns the length of the run, i.e. the number of
/// (moving and non-moving) particles in a raster cell.
/// \param v a BinTriplet object
/// \returns v.endBoundary - v.offset, i.e. the length of the run
inline size_t length(BinTriplet const &v) {
  return v.endBoundary - v.offset;
}

/// output a BinTriplet to ostream.
inline ostream &operator <<(ostream &aus, BinTriplet const &range) {
  aus << "moving: " << range.offset << "-" << range.endMoving
      << " frozen: " << range.endMoving
      << "-" << range.endBoundary;
  return aus;
}

//
// threadprivate variables
//

#ifdef SPH_PARTICLE_MAX_VELOCITY
extern double maxVelocity;
extern double maxRelativeVelocity;
#ifdef _OPENMP
#pragma omp threadprivate(maxVelocity)
#pragma omp threadprivate(maxRelativeVelocity)
#endif
#endif

#ifdef SPH_COUNT_RELATIONS
/// global (per-thread) variable that is used to count relations
extern size_t numPPrelations;
/// global (per-thread) variable that is used to count index tests
extern size_t numPPindextests;
/// global (per-thread) variable that is used to count distance tests
extern size_t numPPdisttests;
#ifdef _OPENMP
#pragma omp threadprivate(numPPrelations)
#pragma omp threadprivate(numPPindextests)
#pragma omp threadprivate(numPPdisttests)
#endif
#endif

struct NachbarArray {
  unsigned const nnachbarn;
  std::valarray<Coord> nachbarArray;
  
  NachbarArray(unsigned const numDims) :
    nnachbarn(dynamic_power(3.0, numDims)),
    nachbarArray(nnachbarn)
  { 
    assert(nachbarArray[0].size() == numDims);
    // this enumerates the numDims-dim. counter over [-1,1]
    // this counts from (-1,-1, ...,-1) up to (1, 1, ..., 1)
    for(unsigned level = 0; level < nnachbarn; ++level) {
      Coord ex;
      for(long d = 0, div = 1; d < long(numDims); ++d, div *= 3) {
        ex[d] = (level/div) % 3 - 1;
      }
      nachbarArray[level] = ex;
    }
  }
  
  unsigned size() const { return nnachbarn; }
  Coord operator()(unsigned const i) const { return nachbarArray[i]; }
  Coord operator[](unsigned const i) const { return nachbarArray[i]; }
};

struct Distance {
#ifdef SPH_PERIODIC_BOUNDARIES
  DoubleVector periodicBoundaryShift, periodicBoundaryMax;
#endif

  Distance(SPHConfig<SphADataType> const &
#ifdef SPH_PERIODIC_BOUNDARIES
           config
#endif
)
#ifdef SPH_PERIODIC_BOUNDARIES
 :
    periodicBoundaryShift(config.gebiet.obenRechts),
    periodicBoundaryMax()
#endif
  {
#ifdef SPH_PERIODIC_BOUNDARIES
    // \fixme test here for exceptions so far (should be none)
    // fetestexcept
    for(unsigned i = 0; i < ndim; ++i) {
      if (config.periodicBoundaries[i]) {
        periodicBoundaryMax[i] = config.gebiet.obenRechts[i] - config.kernelReachH;
      } else {
        //         periodicBoundaryMax[i] = INFINITY;
        periodicBoundaryMax[i] = 1.0 / 0.0;
      }
      //       cerr << "i:" << i << " " << periodicBoundaryMax[i] << "\n";
    }
    // clear exceptions raised by producint infitiy (also when INFINITY is used)
    feclearexcept(FE_ALL_EXCEPT);
#endif
  }

  Vector operator()(Partikel<SphADataType> const &p, Partikel<SphADataType> const &o) const {
    return computeDistance(p, o);
  }
  Vector computeDistance(Partikel<SphADataType> const &p, Partikel<SphADataType> const &o) const {
    Vector res(o.position() - p.position());
#ifdef SPH_PERIODIC_BOUNDARIES
    // if (beide kleiner? was wenn in verschiedenen dimensionen)
    // else if ((o.position() < kernelReachH).max() == 1) {
    //   res += periodicBoundaryShift;
    // } else if ((p.position() < kernelReachH).max() == 1) {
    //   res -= periodicBoundaryShift;
    // }
    // Point<bool, ndim> const larger = abs(res) >= periodicBoundaryMax;
    // if (larger.max() == 1) {
    for(unsigned i = 0; i < ndim; ++i) {
      if (res[i] >= periodicBoundaryMax[i]) {
        do  // the loop is only need for point-symmetric boundary particles
          // where distance between particle and ghost is up to twice as large
          res[i] -= periodicBoundaryShift[i];
        while (res[i] >= periodicBoundaryMax[i]);
      } else if (-res[i] >= periodicBoundaryMax[i]) {
        do 
          res[i] += periodicBoundaryShift[i];
        while (-res[i] >= periodicBoundaryMax[i]);
      }
      assert(fabs(yafad_value(res[i])) <= periodicBoundaryMax[i]);
    }
#endif
    return res;
  }
};

struct RelInfo {
  Vector const diffr;
  Vector::value_type const dnorm2, dnorm;
  
  RelInfo(Distance const &d, Partikel<SphADataType> const &p, Partikel<SphADataType> const &o) :
    diffr(d(p,o)),
    dnorm2(normSquared(diffr)),
    dnorm(sqrt(dnorm2))
  { 
    
  }
};

struct RasterInfo {
  Quader<DoubleVector> const &m_gebiet;
  double const m_h;

  RasterInfo(Quader<DoubleVector> const &gebiet, double const h) :
    m_gebiet(gebiet),
    m_h(h)
  {}

  Coord getIndex(Vector const &p) const {
    return floor(doubleVector(p) / m_h);
  }

};

struct HashedInfo {
  RasterInfo const &m_rasterInfo;
  HashedIndex<CoordType, ndim> const m_index;

  HashedInfo(RasterInfo const &rasterInfo, unsigned const hashM) : 
    m_rasterInfo(rasterInfo), 
    m_index(Coord(ceil(rasterInfo.m_gebiet.v / rasterInfo.m_h)), hashM)
  {}

  void updatePartikelIndex(Particle &p) const {
    assert(p.isOut() == false);
    Vector const &pos = p.position();
    if (m_rasterInfo.m_gebiet.isOut(pos)) {
      p.flags() |= Particle::PARTICLE_FLAG_OUT;
    } else {
      Coord const pindex = m_rasterInfo.getIndex(pos);
      p.hashVal() = m_index(pindex);
#ifdef SPH_PARTICLE_SAVES_RASTER_INDEX
      p.index() = pindex;
#endif
    }
  }
};

template<int Dim>
struct BinaryCounterArray {
  std::valarray<Coord>                 posNachbarArray;

  BinaryCounterArray() : 
    posNachbarArray()
  {
    // this enumerates the Dim-dim. counter over [0,1]
    // this counts from (0,0,, ...,0) up to (1, 1, ..., 1)
    posNachbarArray.resize(static_power<Dim>(2.0));
    for(size_t level = 0; level < posNachbarArray.size(); ++level) {
      Coord ex;
      for(long d = 0, div = 1; d < long(Dim); ++d, div *= 2) {
	ex[d] = (level/div) % 2;
      }
      posNachbarArray[level] = ex;
    }
  }

  Coord const &operator()(unsigned i) const { return posNachbarArray[i]; }
  Coord const &operator[](unsigned i) const { return posNachbarArray[i]; }
  unsigned size() const { return posNachbarArray.size(); }
};

/// template class PartikelManager hides the details of the
/// particle-neighbour lookup from external uses.
/// \tparam _Config the type of the configuration struct
/// \tparam Kernel the type of the Kernel

template<class _Config, class KernelInterface>
struct PartikelManager {
  typedef _Config Config;
  typedef typename Config::value_type value_type;
  typedef Partikel<value_type> particle_type;
  typedef PartikelVariablen<value_type> pv_type;

  // define type PartRelMemFun which is a member function pointer
  typedef void (PartikelManager::*PartRelMemFun1)(particle_type const &p, particle_type const &o,
                                                  RelInfo const &relinfo,
                                                  pv_type &summep) const;
  // define type PartRelMemFun which is a member function pointer
  typedef void (PartikelManager::*PartRelMemFun2)(particle_type const &p, particle_type const &o,
                                                  RelInfo const &relinfo,
                                                  pv_type &summep, 
                                                  pv_type &summeo) const;

  Config const &config;
  KernelInterface *kernel;
  SPHPartikelArray  partikelArray;
//   HashedIndex<CoordType, ndim> const m_index;
  RasterInfo m_rasterInfo;
  HashedInfo m_hashedInfo;
#ifdef SPH_RENDER_FORCES
  std::vector<PartikelForce> partikelForces, partikelDistances;
#endif
  BinaryCounterArray<ndim>                 posNachbarArray;
//   std::valarray<Coord::value_type>     posNachbarArrayFlat;
  LennardJones<value_type> lennardJones;
  size_t numFrozenParticles, numMovingParticles;
  value_type const viscEta2;
  Distance distanceFunction;
//   mutable Timer timeUpdate;
  std::ostream &sphout;

#ifdef SPH_PARTICLE_HAS_ID
  size_t curId;
#endif
  size_t m_numParticlesOut;

  size_t particleTypeMovingDefault;
  size_t particleTypeBoundaryDefault;

  PartRelMemFun1 movmovPartrelFunction1;
  PartRelMemFun1 movbndPartrelFunction1;
  PartRelMemFun2 movmovPartrelFunction2;
  PartRelMemFun2 movbndPartrelFunction2;

  PartikelManager(Config const &_config, std::ostream &sphout) :
    config(_config),
    kernel(KernelInterface::makeKernel(config.kernelName, config.kernelReachH)),
    m_rasterInfo(config.gebiet, config.kernelReachH),
    m_hashedInfo(m_rasterInfo, config.hashFunctionM),
    lennardJones(config.lj_r0, config.lj_D, config.lj_p2),
    numFrozenParticles(), 
    numMovingParticles(),
    viscEta2(config.viscEta*config.viscEta),
    distanceFunction(config),
    sphout(sphout),
#ifdef SPH_PARTICLE_HAS_ID
    curId(),
#endif
    m_numParticlesOut(),
    movmovPartrelFunction1(),
    movbndPartrelFunction1(),
    movmovPartrelFunction2(),
    movbndPartrelFunction2()
  {
// #ifdef SPH_PERIODIC_BOUNDARIES
//     cerr << "periodicBoundaries:" << config.periodicBoundaries << "\n";
//     cerr << "periodicBoundaryMax:" << distanceFunction.periodicBoundaryMax << "\n";
//     cerr << "periodicBoundaryShift:" << distanceFunction.periodicBoundaryShift << "\n";
// #endif
    if (not (yafad_value(lennardJones.r0) > 0)) {
      cerr << "error: lennart-jones range r0 must be GT 0\n";
    }
    if (not (config.kernelReachH > 0)) {
      cerr << "error: bins size h must be GT 0\n";
    }
    if (yafad_value(lennardJones.r0) > config.kernelReachH) {
      cerr << "error: lennart-jones range r0 must be LE kernel reach H\n";
    }
    if (not kernel) {
      cerr << "error: no kernel\n";
    }
    /// \todo particleTypeMovingDefault, ... in config klasse schieben
    particleTypeMovingDefault = particle_type::getParticleClassesValue
      (_config.particleTypeMovingDefaultName.c_str());
    particleTypeBoundaryDefault = particle_type::getParticleClassesValue
      (_config.particleTypeBoundaryDefaultName.c_str());

    switch(particleTypeMovingDefault) {
    case particle_type::PARTICLE_TYPE_MONAGHAN:
      movmovPartrelFunction1 = &PartikelManager::_particleRelationMonaghan;
      movmovPartrelFunction2 = &PartikelManager::particleRelationMonaghan;
      break;
    case particle_type::PARTICLE_TYPE_MONAGHAN_XSPH:
      movmovPartrelFunction1 = &PartikelManager::_particleRelationMonaghanXSPH;
      movmovPartrelFunction2 = &PartikelManager::particleRelationMonaghanXSPH;
      break;
    case particle_type::PARTICLE_TYPE_FERRARI:
      movmovPartrelFunction1 = &PartikelManager::particleRelationFerrari1;
      movmovPartrelFunction2 = &PartikelManager::particleRelationFerrari2;
      break;
    default:
      assert(0);
      break;
    }

    switch(particleTypeBoundaryDefault) {
    case particle_type::PARTICLE_TYPE_POINT_MIRRORED:
      movbndPartrelFunction1 = &PartikelManager::particleRelationPointMirror1;
      movbndPartrelFunction2 = &PartikelManager::particleRelationPointMirror2;
      break;
    case particle_type::PARTICLE_TYPE_LENNARD_JONES:
      movbndPartrelFunction1 = &PartikelManager::particleRelationLennardJones1;
      movbndPartrelFunction2 = &PartikelManager::particleRelationLennardJones2;
      break;
    default:
      assert(0);
      break;
    }
  }
  ~PartikelManager() {  
    if (kernel) {
      delete kernel;
      kernel = 0;
    }
  }

  KernelInterface *setKernel(KernelInterface *newKernel = 0) { 
    KernelInterface *oldKernel = kernel;
    if (newKernel) {
      kernel = newKernel;
    }
    return oldKernel; 
  }

  particle_type const &getPartikel(size_t const i) const { 
    return partikelArray[i]; 
  }

private:

  void doAddBoundary(std::valarray<particle_type> const &neuePartikel) {
    SPHPartikelArray tmp(partikelArray);
    unsigned const nnew = unsigned(neuePartikel.size());
    unsigned const nold = unsigned(tmp.size());
    partikelArray.resize(nnew + nold);
    partikelArray[slice(0, nnew, 1)] = neuePartikel;
    partikelArray[slice(nnew, nold, 1)] = tmp;
  }
  void doAddMoving(std::valarray<particle_type> const &neuePartikel) {
    SPHPartikelArray tmp(partikelArray);
    unsigned const nnew = unsigned(neuePartikel.size());
    unsigned const nold = unsigned(tmp.size());
    partikelArray.resize(nnew + nold);
    partikelArray[slice(0, nold, 1)] = tmp;
    partikelArray[slice(nold, nnew, 1)] = neuePartikel;
  }

  void setupParticles(std::valarray<particle_type> &neuePartikel, std::valarray<bool> &isBoundary) {
    for (unsigned i = 0; i < neuePartikel.size(); ++i) {
      particle_type &inserted = neuePartikel[i];
#ifdef SPH_PARTICLE_HAS_ID
      inserted.id() = ++curId;
#endif
      if (inserted._class() == particle_type::PARTICLE_TYPE_MOVING_DEFAULT) {
//         inserted._class() = particleTypeMovingDefault;
        assert(inserted.masse() != 0.0);
        assert(inserted.dichte() != 0.0);
      }
//       if (inserted._class() == particle_type::PARTICLE_TYPE_BOUNDARY_DEFAULT) {
//         inserted._class() = particleTypeBoundaryDefault;
//       }
      if (inserted._class() == particle_type::PARTICLE_TYPE_LENNARD_JONES
          or inserted._class() == particle_type::PARTICLE_TYPE_POINT_MIRRORED
          or inserted._class() == particle_type::PARTICLE_TYPE_BOUNDARY_DEFAULT) {
        // is a boundary particle
	inserted.flags() |= particle_type::PARTICLE_FLAG_BOUNDARY;
        isBoundary[i] = 1;
      }

      m_hashedInfo.updatePartikelIndex(inserted);
      if (not inserted.isBoundary()) {
	updatePressure(inserted);
      }
    }
  }

  static bool isBoundary(int _class) {
    switch(_class) {
    case particle_type::PARTICLE_TYPE_BOUNDARY_DEFAULT:
    case particle_type::PARTICLE_TYPE_LENNARD_JONES:
    case particle_type::PARTICLE_TYPE_POINT_MIRRORED:
      return 1;
    default:
      break;
    }
    return 0;
  }

  void doAddParticles(std::valarray<particle_type> const &neuePartikel) {
    bool const areBoundary = isBoundary(neuePartikel[0]._class());

    if (areBoundary) {
      doAddBoundary(neuePartikel);
      numFrozenParticles += neuePartikel.size();
    } else {
      doAddMoving(neuePartikel);
      numMovingParticles += neuePartikel.size();
    }

    sphout << "<add-particles"
      " n='" << neuePartikel.size() << "'"
      " moving='" << (areBoundary ? 0 : neuePartikel.size())  << "'"
      " boundary='" << (areBoundary ? neuePartikel.size() : 0) << "'"
      "/>\n";

    assert(numMovingParticles + numFrozenParticles == numParticles());
  }

public:
  size_t numParticles() const { return partikelArray.size(); }
  size_t numParticlesOut() const { 
    return m_numParticlesOut;
    //     return 0; 
  }
  size_t numParticlesLive() const { return numParticles() - numParticlesOut(); }

  void addParticles(std::valarray<particle_type> &neuePartikel) {
    size_t const nnew = neuePartikel.size();
    if (nnew == 0) return;

    std::valarray<bool> istRand(neuePartikel.size());
    setupParticles(neuePartikel, istRand);

    std::valarray<particle_type> 
      rand(neuePartikel[istRand]), nichtRand(neuePartikel[!istRand]);

#ifdef SPH_RENDER_FORCES
    partikelBufferLock.lock();
#endif
    
    if (rand.size() > 0) {
      doAddParticles(rand);
    }
    if (nichtRand.size() > 0) {
      doAddParticles(nichtRand);
    }

#ifdef SPH_RENDER_FORCES
    partikelBufferLock.unlock();
#endif
  }

  value_type viscosity(particle_type const &p, particle_type const &o,
		   Vector const &rdiff, Vector const &vdiff, 
		   value_type const rnormSquared) const { 
    value_type const vdiffXrdiff = scalar(vdiff, rdiff);
    if (vdiffXrdiff > 0) {
      return 0;
    }
    value_type const avgDens = (p.dichte() + o.dichte()) * 0.5;
//     double const r2 = normSquared(rdiff);
    value_type const mu = config.kernelReachH * vdiffXrdiff / 
      (rnormSquared + viscEta2); 

    return ((-config.viscAlpha) * config.speedOfSound * mu
	    + config.viscBeta * mu * mu) / avgDens;
  }

  void _particleRelationMonaghan(particle_type const &p, particle_type const &o,
                                RelInfo const &ri,
                                pv_type &summe) const {
#ifdef SPH_AD_DEACTIVATE_KERNEL
    DoubleVector kernelGrad;
    kernel->gradw(doubleVector(ri.diffr), yafad_value(ri.dnorm), kernelGrad);
#else
    Vector kernelGrad;
    kernel->gradw(ri.diffr, ri.dnorm, kernelGrad);
#endif
    Vector const diffv = o.velocity() - p.velocity();
#ifdef SPH_PARTICLE_MAX_VELOCITY
    maxRelativeVelocity = max(maxRelativeVelocity, yafad_value(norm(diffv)));
    maxVelocity = max(maxVelocity, yafad_value(norm(p.velocity())));
#endif
    SphADataType const factorDruckDichteVis = 
      p.druck() / (p.dichte()*p.dichte()) 
      + o.druck() / (o.dichte()*o.dichte()) 
      + viscosity(p, o, ri.diffr, diffv, ri.dnorm2);
    SphADataType const diffvXkernelGrad = scalar(diffv, kernelGrad);
    summe.velocity() = kernelGrad;
    summe.velocity() *= factorDruckDichteVis;
    summe.dichte() = diffvXkernelGrad;
    summe.waerme() = 0.5 * factorDruckDichteVis * diffvXkernelGrad;
    summe.position() = 0;
  }

  void _particleRelationMonaghanXSPH(particle_type const &p, particle_type const &o,
                                     RelInfo const &ri,
                                     pv_type &summe) const {
#ifdef SPH_AD_DEACTIVATE_KERNEL
    DoubleVector kernelGrad;
    double const kernelVal = kernel->wAndGrad(doubleVector(ri.diffr), yafad_value(ri.dnorm), kernelGrad);
#else
    Vector kernelGrad;
    SphADataType const kernelVal = 
      kernel->wAndGrad(doubleVector(ri.diffr), yafad_value(ri.dnorm), kernelGrad);
#endif
    Vector const diffv = o.velocity() - p.velocity();
#ifdef SPH_PARTICLE_MAX_VELOCITY
    maxRelativeVelocity = max(maxRelativeVelocity, yafad_value(norm(diffv)));
    maxVelocity = max(maxVelocity, yafad_value(norm(p.velocity())));
#endif
    SphADataType const factorDruckDichteVis = 
      p.druck() / (p.dichte()*p.dichte()) 
      + o.druck() / (o.dichte()*o.dichte()) 
      + viscosity(p, o, ri.diffr, diffv, ri.dnorm2);
    SphADataType const diffvXkernelGrad = scalar(diffv, kernelGrad);
    summe.velocity() = kernelGrad;
    summe.velocity() *= factorDruckDichteVis;
    Vector xsphvelmod = (kernelVal / ((p.dichte() + o.dichte())/2.0)) * diffv;
    summe.dichte() = diffvXkernelGrad;
    summe.waerme() = 0.5 * factorDruckDichteVis * diffvXkernelGrad;
    summe.position() = config.xsphEpsilon * xsphvelmod;
  }

  void particleRelationMonaghan(particle_type const &p, particle_type const &o,
                                RelInfo const &ri,
                                pv_type &summep, 
                                pv_type &summeo) const {
    _particleRelationMonaghan(p, o, ri, summep);
#ifdef SPH_PARTICLE_MAX_VELOCITY
    maxVelocity = max(maxVelocity, yafad_value(norm(o.velocity())));
#endif
    summeo.waerme() = summep.waerme();
    summeo.dichte() = summep.dichte();
    summeo.velocity() = -summep.velocity();
    summeo.position() = 0;
  }
  void particleRelationMonaghanXSPH(particle_type const &p, particle_type const &o,
                                    RelInfo const &ri,
                                    pv_type &summep, 
                                    pv_type &summeo) const {
    _particleRelationMonaghanXSPH(p, o, ri, summep);
#ifdef SPH_PARTICLE_MAX_VELOCITY
    maxVelocity = max(maxVelocity, yafad_value(norm(o.velocity())));
#endif
    summeo.waerme() = summep.waerme();
    summeo.dichte() = summep.dichte();
    summeo.velocity() = -summep.velocity();
    summeo.position() = -summep.position();
  }

  value_type celerity(particle_type const &p) const {
    return sqrt(double(config.eqsConfig.gamma)
                * (config.eqsConfig.pressureB / config.eqsConfig.rho0) 
                * dynamic_power(p.dichte() / config.eqsConfig.rho0, int(config.eqsConfig.gamma-1)) );
  }

  void particleRelationFerrari1(particle_type const &p, particle_type const &o,
                                RelInfo const &ri,
                                pv_type &summep) const {
#ifdef SPH_AD_DEACTIVATE_KERNEL
    DoubleVector kernelGrad;
    kernel->gradw(doubleVector(ri.diffr), yafad_value(ri.dnorm), kernelGrad);
#else
    Vector kernelGrad;
    kernel->gradw(ri.diffr, ri.dnorm, kernelGrad);
#endif
    Vector const diffv = o.velocity() - p.velocity();

#ifdef SPH_PARTICLE_MAX_VELOCITY
    maxRelativeVelocity = max(maxRelativeVelocity, yafad_value(norm(diffv)));
    maxVelocity = max(maxVelocity, yafad_value(norm(p.velocity())));
#endif
    
    SphADataType const factorInviscid = 
      p.druck() / (p.dichte()*p.dichte()) 
      + o.druck() / (o.dichte()*o.dichte());

    Vector const nij = ri.diffr / ri.dnorm;
    SphADataType const nijXkernelGrad = nij.scalar(kernelGrad);
    SphADataType const thetij = - nijXkernelGrad / ri.dnorm;
    summep.velocity() = kernelGrad;
    summep.velocity() *= factorInviscid;
    if (config.dynamicBulkViscosity != 0.0) {
      Vector const factorViscous = 
        (config.dynamicBulkViscosity * thetij / (p.dichte() * o.dichte()))
        * (SphADataType(7.0/3.0) * diffv + SphADataType(5.0/3.0) * (nij * nij.scalar(diffv)));
      summep.velocity() += factorViscous;
    }
    
    SphADataType const cij = max(celerity(p), celerity(o));
    SphADataType const common = nijXkernelGrad * cij;
    summep.dichte() =
      diffv.scalar(kernelGrad);
    summep.dichte() -=  common * (o.dichte() - p.dichte())/o.dichte();

    summep.waerme() = 0;
    summep.position() = 0;
   
  }
  
  void particleRelationFerrari2(particle_type const &p, particle_type const &o,
                                RelInfo const &ri,
                                pv_type &summep, pv_type &summeo) const {
    assert(not p.isBoundary());
    assert(not o.isBoundary());

#ifdef SPH_AD_DEACTIVATE_KERNEL
    DoubleVector kernelGrad;
    kernel->gradw(doubleVector(ri.diffr), yafad_value(ri.dnorm), kernelGrad);
#else
    Vector kernelGrad;
    kernel->gradw(ri.diffr, ri.dnorm, kernelGrad);
#endif
    Vector const diffv = o.velocity() - p.velocity();

#ifdef SPH_PARTICLE_MAX_VELOCITY
    maxRelativeVelocity = max(maxRelativeVelocity, yafad_value(norm(diffv)));
    maxVelocity = max(maxVelocity, yafad_value(norm(p.velocity())));
    maxVelocity = max(maxVelocity, yafad_value(norm(o.velocity())));
#endif
    
    SphADataType const factorInviscid = 
      p.druck() / (p.dichte()*p.dichte()) 
      + o.druck() / (o.dichte()*o.dichte());

    Vector const nij = ri.diffr / ri.dnorm;
    SphADataType const nijXkernelGrad = nij.scalar(kernelGrad);
    SphADataType const thetij = - nijXkernelGrad / ri.dnorm;
    Vector const factorViscous = 
      (config.dynamicBulkViscosity * thetij / (p.dichte() * o.dichte()))
      * ( SphADataType(7.0/3.0) * diffv
          + (SphADataType(5.0/3.0) * nij.scalar(diffv)) * nij);
//     summep.velocity() = factorInviscid * kernelGrad + factorViscous;
    summep.velocity() = kernelGrad;
    summep.velocity() *= factorInviscid;
    summep.velocity() += factorViscous;
    summeo.velocity() = -summep.velocity();
    
    SphADataType const cij = max(celerity(p), celerity(o));
    SphADataType const common = nijXkernelGrad * cij;
    summep.dichte() = summeo.dichte() = 
      diffv.scalar(kernelGrad);
    summep.dichte() -=  common * (o.dichte() - p.dichte())/o.dichte();
    summeo.dichte() -=  common * (p.dichte() - o.dichte())/p.dichte();

    summep.waerme() = summeo.waerme() = 0;
    summep.position() = summeo.position() = 0;
   
  }

  void _particleRelationLennardJones(RelInfo const &ri,
                                     pv_type &summe) const {
    summe.velocity() = ri.diffr;
    summe.velocity() *= -lennardJones.w(ri.dnorm) / (ri.dnorm2);
    
    summe.position() = 0;
    summe.dichte() = 0;
    summe.waerme() = 0;
  }

  void particleRelationLennardJones1(particle_type const &, particle_type const &,
                                     RelInfo const &ri,
                                     pv_type &summep) const {
    _particleRelationLennardJones(ri, summep);
  }

  void particleRelationLennardJones2(particle_type const &, particle_type const &,
                                     RelInfo const &ri,
                                     pv_type &summep, pv_type &summeo) const {
    _particleRelationLennardJones(ri, summep);
    summeo = 0;
  }

  void particleRelationPointMirror1(particle_type const &p, particle_type const &mp
                                    , RelInfo const &
#ifdef DEBUG
                                    ri
#endif
                                    ,     pv_type &summep) const {
    assert(not p.isBoundary());
    assert(mp.isBoundary());
//     if (ri.dnorm > config.kernelReachH/2.0) {
//       return;
//     }
    particle_type ghost(p);
    ghost.position() = mp.position()*value_type(2.0) - p.position();
    ghost.velocity() = mp.velocity()*value_type(2.0) - p.velocity();
    ghost.dichte() = p.dichte();
    ghost.druck() = p.druck();
    ghost.masse() = p.masse();
    RelInfo ri2(distanceFunction, p, ghost);
    sph_hard_assert(ri2.dnorm - ri.dnorm * 2 < 1e-10 * ri2.dnorm);
    (this->*movmovPartrelFunction1)(p, ghost, ri2, summep);
  }

  void particleRelationPointMirror2(particle_type const &p, particle_type const &mp,
                                    RelInfo const &ri,
                                    pv_type &summep, pv_type &) const {
    particleRelationPointMirror1(p, mp, ri, summep);
  }

//   void particleRelationNoBoundary(particle_type const &p, particle_type const &o,
//                                   RelInfo const &ri,
//                                   pv_type &summep, pv_type &summeo) {
//     assert(not p.isBoundary());
//     assert(not o.isBoundary());
//     assert(p._class() == o._class());

//     movmovPartrelFunction(p, o, ri, summep, summeo);
//   }

private:
  void particleRelation(particle_type const &p, particle_type const &o,
                        RelInfo const &ri,
                        pv_type &summep, pv_type &summeo) const {
    assert(not p.isBoundary());
    if (o.isBoundary()) {
      (this->*movbndPartrelFunction2)(p, o, ri, summep, summeo);
    } else {
      (this->*movmovPartrelFunction2)(p, o, ri, summep, summeo);
      summep   *= o.masse();
      summeo   *= p.masse();
    }
  }

  void particleRelation(particle_type const &p, particle_type const &o,
                        RelInfo const &ri,
                        pv_type &summep) const {
    assert(not p.isBoundary());
    if (o.isBoundary()) {
      (this->*movbndPartrelFunction1)(p, o, ri, summep);
    } else {
      (this->*movmovPartrelFunction1)(p, o, ri, summep);
      summep   *= o.masse();
    }
  }

public:
/*! Particle relation function f, non symmetric version.
  \param p test particle
  \param o other particle
  \param summep result f(p, o) (the result is for p)
 */
  bool particleRelation(particle_type const &p, particle_type const &o,
                        pv_type &summep) const {
    RelInfo const ri(distanceFunction, p, o);
    sph_hard_assert(ri.dnorm <= 2*sqrt(config.kernelReachH*config.kernelReachH*ndim)
                    or ri.dnorm >= (ldexp(1.0, config.hashFunctionM)-2)*config.kernelReachH);
    if (ri.dnorm <= config.kernelReachH) {
      particleRelation(p, o, ri, summep);
#ifdef SPH_PARTICLE_VAR_HAS_ID_SUM
      summep.m_sumOfIds = o.id();
#endif
#ifdef SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM
      summep.m_numNeighbours = 1;
#endif
//       assert(not FLInfo::testAndShowExceptions
//                       (std::cerr, FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
//       assert(not summep.hasNaN());
      return true;
    } else {
      return false;
    }
  }

/*! Particle relation function f, symmetric version.
  \param p test particle
  \param o other particle
  \param summep result f(p, o) (the result is for p)
  \param summep result f(o, p) (the result is for o)
 */
  bool particleRelation(particle_type const &p, particle_type const &o,
                        pv_type &summep, pv_type &summeo) const {
    RelInfo const ri(distanceFunction, p, o);
    sph_hard_assert(ri.dnorm <= 2*sqrt(config.kernelReachH*config.kernelReachH*ndim)
                    or ri.dnorm >= (ldexp(1.0, config.hashFunctionM)-2)*config.kernelReachH);
    if (ri.dnorm <= config.kernelReachH) {
      particleRelation(p, o, ri, summep, summeo);
#ifdef SPH_PARTICLE_VAR_HAS_ID_SUM
      summep.m_sumOfIds = o.id();
      summeo.m_sumOfIds = p.id();
#endif
#ifdef SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM
      summep.m_numNeighbours = 1;
      summeo.m_numNeighbours = 1;
#endif
      // cerr << "p: " << p << "\n";
      // cerr << "o: " << o << "\n";
      // cerr << "sp: " << summep << "\n";
      // cerr << "so: " << summeo << "\n";
//        assert(not FLInfo::testAndShowExceptions
//                        (std::cerr, FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
//        assert(not summep.hasNaN());
//        assert(not summeo.hasNaN());
      return true;
    } else {
      return false;
    }
  }

  void updateState(SPHPartikelArray &state, 
                   SPHPartikelVariablenArray const &update) const {
    assert(update.size() == state.size());
    int const von = numFrozenParticles, bis = state.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i = von; i < bis; ++i) {
      particle_type &p = state[i];
      assert(not p.isBoundary());
      if (p.isOut()) continue;
#ifdef SPH_PARTICLE_HAS_ID
      assert(update[i].id() == p.id() or update[i].isZero());
#endif
      // update the values of particle i with updates i
      p.werte() += update[i];
#ifdef SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM
       // overwrite prev. addition
      p.werte().numNeighbours() = update[i].numNeighbours();
#endif
#ifdef SPH_PARTICLE_VAR_HAS_ID_SUM
       // overwrite prev. addition
      p.werte().sumOfNeighbourIds() = update[i].sumOfNeighbourIds();
#endif
#ifdef SPH_PERIODIC_BOUNDARIES
      doPeriodicBoundaries(p);
#endif
      m_hashedInfo.updatePartikelIndex(p);
      updatePressure(p);
    }
  }

#ifdef SPH_PERIODIC_BOUNDARIES
  void doPeriodicBoundaries(particle_type &p) const {
    Vector &pos = p.position(); 
    for(unsigned i = 0; i < ndim; ++i) {
      if (config.periodicBoundaries[i]) {
        if (pos[i] > config.gebiet.obenRechts[i]) {
          pos[i] -= config.gebiet.obenRechts[i];
        } else if (pos[i] < 0.0) {
          pos[i] += config.gebiet.obenRechts[i];
        }
      }
    }
  }
#endif

  value_type computePressureMonaghanFerrari(particle_type const &p) const {
    return config.eqsConfig.pressureB *
      (pow(p.dichte() / config.eqsConfig.rho0, int(config.eqsConfig.gamma)) - 1);
  }

  value_type computePressure_RealTaits(particle_type const &p) const {
    return config.eqsConfig.pressureK0 * 
      (p.dichte() - config.eqsConfig.rho0)
      / p.dichte();
  }

  value_type computePressure3(particle_type const &p) const {
    return config.eqsConfig.pressureB *
      expm1(log(p.dichte() / config.eqsConfig.rho0) * config.eqsConfig.gamma);
  }

  value_type computePressure(particle_type const &p) const {
    return computePressureMonaghanFerrari(p);
  }

  void updatePressure(particle_type &p) const {
    p.druck() = computePressure(p);
  }

  struct ForceFunctorInterface {
    typedef SphADataType first_argument_type;
    typedef SPHPartikelArray value_type;
    typedef SPHPartikelArray second_argument_type;
    typedef SPHPartikelVariablenArray result_type;
    
    struct UnknownSummationAlgorithm{};

    virtual ~ForceFunctorInterface() {}
    virtual result_type const &operator()(first_argument_type const dt, 
                            second_argument_type const &partikelListe) = 0;
    virtual void update(value_type &state, 
                        SPHPartikelVariablenArray const &update) = 0;
    virtual size_t numParticlesOut() const = 0;
    virtual ForceFunctorInterface *clone() = 0;
    static ForceFunctorInterface *makeForceFunctor(std::string const &name, 
                                                   PartikelManager const &pm);
  };
  
  struct ForceFunctorWrapper {
    typedef SphADataType first_argument_type;
    typedef SPHPartikelArray value_type;
    typedef SPHPartikelArray second_argument_type;
    typedef SPHPartikelVariablenArray result_type;

    ForceFunctorInterface *ff;
    mutable Timer timeUpdate, timeF;
    size_t numCalls, numUpdates;
    PartikelManager const &pm;

    ForceFunctorWrapper(PartikelManager const &pm) :
      ff(ForceFunctorInterface::makeForceFunctor(pm.config.summationAlgorithm, pm)),
      numCalls(), numUpdates(),
      pm(pm)
    {
    }
    ForceFunctorWrapper(ForceFunctorWrapper &o) :
      numCalls(), numUpdates()
    {
      ff = o.ff->clone();
    }
    ~ForceFunctorWrapper() {
      if (ff) {
	delete ff;
	ff = 0;
      }
    }
    result_type const &operator()(first_argument_type const dt, 
                                  second_argument_type const &partikelListe) {
      timeF.start();
      ++numCalls;
      result_type const &res = ff->operator()(dt, partikelListe);
      timeF.stop();
      if (pm.config.detailedTimings) {
        pm.sphout << "<eval-state n='" << numCalls << "' time='" << timeF.lastDiff() << "'/>\n";
      }
      return res;
    }
    void update(value_type &state, SPHPartikelVariablenArray const &update) {
      timeUpdate.start();
      ++numUpdates;
      ff->update(state, update);
      timeUpdate.stop();
      if (pm.config.detailedTimings) {
        pm.sphout << "<update-state n='" << numUpdates << "' time='" << timeUpdate.lastDiff() << "'/>\n";
      }
    }
    size_t numParticlesOut() const {
      return ff->numParticlesOut();
    }
  private:
    ForceFunctorWrapper &operator =(ForceFunctorWrapper const &) {
      return *this;
    }
  };

  /// The naive version of the force functor: This simply considers
  /// every particle a in turn, looks up their neighbors and computes
  /// the "forces" (differential quantities) exerted by the neighbors
  /// on a. This will obviously consider every pair of particles
  /// twice.
  struct ForceFunctorPerParticle : public ForceFunctorInterface {
    typedef ForceFunctorInterface Base;
    typedef typename Base::first_argument_type first_argument_type;
    typedef typename Base::value_type value_type;
    typedef typename Base::second_argument_type second_argument_type;
    typedef typename Base::result_type result_type;

    PartikelManager const &pm;

    Timer timeIndex, timeSummation;
    NachbarArray const nachbarArray;
    result_type shadowArray;
    size_t m_numParticlesOut;

    ForceFunctorPerParticle(PartikelManager const &pm) : 
      pm(pm),
      nachbarArray(ndim),
      shadowArray(pm.partikelArray.size()),
      m_numParticlesOut()
    {
    }

    ForceFunctorPerParticle *clone() {
      return new ForceFunctorPerParticle(*this);
    }
    
    size_t numParticlesOut() const {
      return m_numParticlesOut;
    }

    void sumParticleWithParticlesInCell(SortedPartikelIndex<SPHPartikelArray> const &partikelIndex,
                                        SPHPartikelArray const &partikelListe,
                                        particle_type const &p, BinTriplet const &range, 
                                        pv_type &summe) {
      for (int oi = range.offset; oi < int(range.endBoundary); ++oi) {
        int const partj = partikelIndex.getPartikelIndex(oi);
        particle_type const &o = partikelListe[partj];
#ifdef SPH_PARTICLE_HAS_ID
        if (p.id() == o.id()) {
#else
        if (&p == &o) {
#endif
          // nachbarArray[ni] == (0, 0, ...)
          assert((o.index() == p.index()).min() == 1);
          assert(o.hashVal() == p.hashVal());
          continue;
        }
        pv_type tmp;
        if (pm.particleRelation(p, o, tmp)) {
          summe += tmp;
        }
      }
    }

    result_type const &operator()(first_argument_type const, 
                                  value_type const &partikelListe) {
      
      if (shadowArray.size() != partikelListe.size()) {
        shadowArray.resize(partikelListe.size());
      }

//       shadowArray = ParticleVariable();

#ifdef SPH_COUNT_RELATIONS
#ifdef _OPENMP
#pragma omp parallel
#endif
      {
	numPPrelations = 0;
// 	numPPindextests = 0;
	numPPdisttests = 0;
      }
#endif
      
#ifdef SPH_RENDER_FORCES
      pm.partikelDistances.clear();
      pm.partikelForces.clear();
#endif

      timeIndex.start();

      unsigned const curNumThreads = omp_get_max_threads();
      unsigned const numThreadsForSort = min(curNumThreads, pm.config.numCSortThreads);
//       { cant'do
        omp_set_num_threads(numThreadsForSort);

        SortedPartikelIndex<SPHPartikelArray> 
          partikelIndex(pm.config.sortedIndex, pm.m_hashedInfo.m_index, 
                        pm.config.gebiet, pm.config.kernelReachH);
        partikelIndex.update(partikelListe);
        
        omp_set_num_threads(curNumThreads);
//       } cant'do

      timeIndex.stop();
      if (pm.config.detailedTimings) {
        partikelIndex.writeXMLTimings(pm.sphout);
        pm.sphout << "<update-index time='" << timeIndex.lastDiff() << "' threads='" << numThreadsForSort << "'/>\n";
      }

      m_numParticlesOut = partikelIndex.numParticlesOut();

      long const np = partikelListe.size();

      timeSummation.start();

#ifdef _OPENMP
#pragma omp parallel
#endif
      {
        BinTriplet range;
      
#ifdef _OPENMP
        // leave schedule unconfigured here to enable overide in environment
#pragma omp for schedule(dynamic)
// #pragma omp for schedule(dynamic)
// #pragma omp for schedule(static)
// #pragma omp for
#endif
        for (int parti = pm.numFrozenParticles; parti < np; ++parti) {
          
          particle_type const &p = partikelListe[parti];
          assert(not p.isBoundary());
          if (p.isOut()) {
            continue;
          }
          
          pv_type summe;

          for (unsigned ni = 0; ni < nachbarArray.size(); ++ni) {
            Coord const otherIndex = p.index() + nachbarArray[ni];
            long const otherIndexFlat = pm.m_hashedInfo.m_index(otherIndex);
            partikelIndex.getBinInfo(otherIndexFlat, range);
//             cerr << "normal: " << p.index() << " <-> " << otherIndex
//                           << " " << otherIndexFlat << "\n";
            sumParticleWithParticlesInCell(partikelIndex, partikelListe, p, range, summe);
          }

#ifdef SPH_PERIODIC_BOUNDARIES
          // handle periodic boundaries
          for(unsigned i = 0; i < ndim; ++i) {
            if (pm.config.periodicBoundaries[i]) {
              if (p.index()[i] == 0) { // particle in first line of cells
                Vector pAltPosition = p.position();
                pAltPosition[i] = pm.config.gebiet.obenRechts[i];
                Coord const pAltIndex = pm.m_rasterInfo.getIndex(pAltPosition);
                for (unsigned ni = 0; ni < nachbarArray.size(); ++ni) {
                  if (nachbarArray[ni][i] <= 0) {
                    Coord const otherIndex = pAltIndex + nachbarArray[ni];
                    int const hashedRow = pm.m_hashedInfo.m_index.hashFunktion(otherIndex[i]);
                    if (hashedRow <= 1 
                        or hashedRow >= int(pm.m_hashedInfo.m_index.hashFunktion.m_dim)-1)
                      continue;
                    long const otherIndexFlat = pm.m_hashedInfo.m_index(otherIndex);
//                     cerr << "1. handling per. boundaries: (" << 
//                       p.position() << " == " << pAltPosition << ")\n";
//                     cerr << "1. handling per. boundaries: (" << 
//                       p.index() << " == " << pAltIndex << ") <-> " << otherIndex
//                          << " " << otherIndexFlat << "\n";
                    assert(abs(otherIndex - p.index()).max() > 1);
                    assert(abs(otherIndex - pAltIndex).max() <= 1);
                    partikelIndex.getBinInfo(otherIndexFlat, range);
                    sumParticleWithParticlesInCell(partikelIndex, partikelListe, p, range, summe);
                  }
                } // for
              } else if (pm.config.gebiet.obenRechts[i] - p.position()[i] <= pm.config.kernelReachH) {
                // particle in one of last two lines of cells
                Vector pAltPosition = p.position();
                pAltPosition[i] = 0;
                Coord const pAltIndex = pm.m_rasterInfo.getIndex(pAltPosition);
                int const hashedRow = pm.m_hashedInfo.m_index.hashFunktion(p.index()[i]);
                if (hashedRow <= 1 
                    or hashedRow >= int(pm.m_hashedInfo.m_index.hashFunktion.m_dim)-1) {
                  // nothing to do
                } else {
                  for (unsigned ni = 0; ni < nachbarArray.size(); ++ni) {
                    if (nachbarArray[ni][i] == 0) {
                      Coord const otherIndex = pAltIndex + nachbarArray[ni];
                      long const otherIndexFlat = pm.m_hashedInfo.m_index(otherIndex);
//                       cerr << "2. handling per. boundaries: (" << 
//                         p.position() << " == " << pAltPosition << ")\n";
//                       cerr << "2. handling per. boundaries: (" << 
//                         p.index() << " == " << pAltIndex << ") <-> " << otherIndex
//                            << " " << otherIndexFlat << "\n";
                      assert(abs(otherIndex - p.index()).max() > 1);
                      partikelIndex.getBinInfo(otherIndexFlat, range);
                      sumParticleWithParticlesInCell(partikelIndex, partikelListe, p, range, summe);
                    }
                  } // for
                } // else
              }
            }
          }
#endif
          
          // take care that .position() does never get a derivative
          // \todo find reason for this
#ifdef SPH_AD
#ifdef SPH_AD_KILL_DERIVATIVES_OF_POSITION_UPDATE
          // AD, and make sure here position() of moving particles will never have derivatives
          summe.position() += doubleVector(p.velocity());
#else
          // AD, but do not treat these derivs specially: proceed as planned        
          summe.position() += p.velocity();
#endif
#else
          // no AD: proceed as planned        
          summe.position() += p.velocity();
#endif
          summe.velocity() += pm.config.gravityVector;

#ifdef SPH_PARTICLE_HAS_ID
          summe.id() = p.id();
#endif
          shadowArray[parti] = summe;

        } // for

      } // parallel region

      timeSummation.stop();
      if (pm.config.detailedTimings) {
        pm.sphout << "<summation time='" << timeSummation.lastDiff() << "' type='naive'/>\n";
      }

      return shadowArray;
    }

    void update(value_type &state, SPHPartikelVariablenArray const &update) {
      pm.updateState(state, update);
    }
  };

  /// The ForceFunctor computes the updates of the differential values
  /// all particles.  This object is passed to the integrator chosen
  /// at compiletime.  The result is returned in a SPHPartikelVariablenArray
  /// array of particle variables.  The class also implements the
  /// function update() which adds the differential updates to the
  /// particles managed by the PartikelManager.  
  /// 
  /// The summation of forces is done in parallel, yet each
  /// particle-particle force is computed only once.  To understand
  /// the idea, it is best to consider the 1D case first.
  /// 
  /// \image html parallelization-1d.png "The 1D case"
  ///
  /// The kernel influence length \f$ H \f$ induces a global raster of
  /// width \f$ H \f$ on the computation domain. Each particle is
  /// given the integer index of the raster bin it is in, by the
  /// HashedIndex class. To sum the forces of a particle in bin 2, it
  /// is then only necessary to look in bins 1, 2 and 3 for particles
  /// that are closer than distance \f$ H \f$.
  /// 
  /// The HashFunktion has the property that adjacent raster indices
  /// remain adjacent in the hashed raster.
  /// 
  /// The parallelization is done over the raster. The two helper
  /// functions sumForcesIn1Bin() and sumForcesIn2Bins() are used to
  /// compute and store the forces between all particles in one bin
  /// and in two adjacent bins, resp. When they compute a
  /// particle-particle force (or other differential), they add this
  /// force on the summation variables of BOTH particles.
  /// 
  /// The idea is now to perform two passes, one EVEN (\f$o=0\f$) and
  /// one ODD (\f$o=1\f$).  In one pass, EVERY SECOND raster bin 
  /// \f$ i = 2k + o \f$ is considered, calling sumForcesIn1Bin(i) and
  /// sumForcesIn2Bins(i, #(i+1)), where #() is the hash function.
  ///
  /// Each pass is arbitrarily parallelizable over \f$k \f$! This
  /// means the scalability is very good, because there is no
  /// communication across boundaries, etc.
  ///
  /// The generalization to 2D and 3D however is not straight forward.
  /// The problem is most easily solved graphically:
  /// 
  /// \image html parallelization-2d.png "The 2D case"
  ///
  /// The "oddity" becomes a 2 dimensional vector of bits \f$ {\bf o} \f$, there are
  /// four different 2D offsets and thus four passes. Each pass consideres every second
  /// raster bin in each direction adding the passes "oddity" offset:  
  /// \f$ {\bf I} = 2(k_x, k_y) + {\bf o} \f$. The 2D stencil 
  /// is as shown in the image: for a given 2D-raster index \f$ {\bf I}\f$
  /// the function sumForcesIn2Bins() is called for the pairs \f$ [{\bf I}, {\bf I} + (0,1)] \f$, 
  /// \f$ [{\bf I}, {\bf I} + (1,0)] \f$, \f$ [{\bf I}, {\bf I} + (1,1)] \f$
  /// and \f$ [{\bf I} + (1,0), {\bf I} + (0,1)] \f$. 
  /// Finally function sumForcesIn1Bin(\f$ {\bf I}\f$) is called.
  /// 
  /// Again, the computations of each of the four passes are
  /// parallelizable in any manner, because all computations for each
  /// raster position \f$ {\bf I}\f$ only update different variables,
  /// thus avoiding any race-conditions.
  ///
  /// In 3D the oddity becomes a 3D bitvector with 8 possible
  /// states. The 3D stencil is as shown in the next image:
  /// 
  /// \image html stencil-3d.jpg "The 3D stencil: X-axis red, Y-axis green, Z-axis blue"
  ///
  /// The problem how to enumerate the pairs of raster bins in a general way for all dimensions
  /// is solved now: sumForcesIn2Bins() is called for each pair \f$ [I + A, I + B] \f$
  /// where 
  ///   - \f$ A \ne B \f$, 
  ///   - \f$ A \cdot B = 0\f$ and 
  ///   - \f$ A \in [0,1]^D,\; B \in [0,1]^D\f$. 
  /// 
  /// By precomputing all \f$ [0,1]^D \f$ and storing them in an array
  /// \f$ v \f$ in the order they form a counter over \f$ [0,1]^D \f$,
  /// we have the property that \f$ v[i] = i \f$ (the numbers only
  /// represented in a different manner) and condition 2 can be tested
  /// by a bit \b and operation: 
  /// 
  /// \f[ v[i] \cdot v[j] = 0 \equiv i \,\&\, j = 0. \f]  
  ///
  /// Also note that the vector \f$ v \f$ contains exactly the offset
  /// vectors for the \f$ 2^D \f$ passes.
  ///
  /// Thus, in pseudo code, the algorithm can be formulated as follows, working for all values
  /// of D:
  ///
  /** \verbatim
         for o in [0, 2^D - 1]:
      #omp parallel for
            for i in [0, 2^{m-1} - 1]^D:
                I = 2i + v[o]
                for j in [0, 2^D - 1]:
                    for k in [j + 1, 2^D - 1]:
                        if ((j & k) == 0):
                           sumForcesIn2Bins(I + v[j], I + v[k])
                sumForcesIn1Bin(I)
      \endverbatim */
  struct ForceFunctorSymmetric : public ForceFunctorInterface {
    PartikelManager const &pm;
    typedef SphADataType first_argument_type;
    typedef SPHPartikelArray value_type;
    typedef SPHPartikelArray second_argument_type;
    typedef SPHPartikelVariablenArray result_type;
    typedef SortedPartikelIndex<SPHPartikelArray> PartikelIndex;

    Timer timeIndex, timeSummation;
    std::valarray<Coord> gridIndices;
    result_type m_shadowArray;
    size_t m_numParticlesOut;

    ForceFunctorSymmetric(PartikelManager const &_pm) : 
      pm(_pm),
      gridIndices(static_power<ndim>(1 << (pm.config.hashFunctionM - 1))),
      m_shadowArray(pm.partikelArray.size()),
      m_numParticlesOut()
    {
      Coord halvedGrid_Top;
      halvedGrid_Top = 1 << (pm.config.hashFunctionM - 1);
      HashedIndex<long, ndim> gridIndex(halvedGrid_Top, pm.config.hashFunctionM - 1);
      Coord coord;
      for(size_t i = 0; i < gridIndices.size(); ++i) {
	coord = gridIndex.invert(i);
	coord <<= 1;
	gridIndices[i] = coord;
	// 	cerr << "set gridIndizes[" << i << "] = " << coord << "\n";
      }
      cerr << "creating ForceFunctorSymmetric \n";
    }

    ForceFunctorSymmetric *clone() {
      return new ForceFunctorSymmetric(*this);
    }

    size_t numParticlesOut() const {
      return m_numParticlesOut;
    }

    void sumForcesIn1Bin(second_argument_type const &partikelListe, 
                         result_type &shadowArray,
                         PartikelIndex const &partikelIndex,
                         BinTriplet const &range) {
      DEB_SPH("sumForcesIn1Bin: " << length(range));
      pv_type summep, summeo;
      unsigned li = range.offset;
      // loop over MOVING particles
      for ( ; li < range.endMoving; ++li) {
	unsigned const partp = partikelIndex.getPartikelIndex(li);
        particle_type const &p = partikelListe[partp];
	assert(not p.isBoundary());
	pv_type &forcesP = shadowArray[partp];
#ifdef SPH_AD
#ifdef SPH_AD_KILL_DERIVATIVES_OF_POSITION_UPDATE
        // AD, and make sure here position() of moving particles will never have derivatives
	forcesP.position() += doubleVector(p.velocity());
#else
        // AD, but do not treat these derivs specially: proceed as planned        
	forcesP.position() += p.velocity();
#endif
#else
        // no AD: proceed as planned        
	forcesP.position() += p.velocity();
#endif
	forcesP.velocity() += pm.config.gravityVector;
	
#ifdef SPH_PARTICLE_HAS_ID
        assert(shadowArray[partp].m_id == 0 or shadowArray[partp].m_id == p.id());
        shadowArray[partp].m_id = p.id();
#endif

	// loop over all particles in same bin (MOVING)
        unsigned lli = li + 1;
	for ( ; lli < range.endMoving; ++lli) {
          unsigned const parto = partikelIndex.getPartikelIndex(lli);
	  particle_type const &o = partikelListe[parto];
          assert(not o.isBoundary());
	  if (pm.particleRelation(p, o, summep, summeo)) {
#ifdef SPH_PARTICLE_HAS_ID
            assert(shadowArray[partp].id() == 0 or shadowArray[partp].id() == p.id());
            shadowArray[partp].id() = p.id();
            assert(shadowArray[parto].id() == 0 or shadowArray[parto].id() == o.id());
            shadowArray[parto].id() = o.id();
#endif
            shadowArray[partp] += summep;
            shadowArray[parto] += summeo;
          }
	}
	// loop over all particles in same bin (BOUNDARY)
	for ( ; lli < range.endBoundary; ++lli) {
          unsigned const parto = partikelIndex.getPartikelIndex(lli);
	  particle_type const &o = partikelListe[parto];
          assert(o.isBoundary());
          // \todo use particleRelationBoundary shortcut
	  if (pm.particleRelation(p, o, summep)) {
#ifdef SPH_PARTICLE_HAS_ID
            assert(shadowArray[partp].id() == 0 or shadowArray[partp].id() == p.id());
            shadowArray[partp].id() = p.id();
#endif
            shadowArray[partp] += summep;
          }
	}
      }
    }

    BinTriplet getRange(PartikelIndex const &partikelIndex, Coord const &gridPos1) {
      BinTriplet range1;
      size_t const binIndex = this->pm.m_hashedInfo.m_index(gridPos1);
      partikelIndex.getBinInfo(binIndex, range1);
      return range1;
    }

    void sumForcesIn1BinCoords(second_argument_type const &partikelListe, 
                               result_type &shadowArray,
                               PartikelIndex const &partikelIndex,
                               Coord const &gridPos1) {
      BinTriplet range1;
                
      size_t const binIndex = this->pm.m_hashedInfo.m_index(gridPos1);
      partikelIndex.getBinInfo(binIndex, range1);
                
      DEB_SPH("sumForcesIn1BinCoords: " << gridPos1 << " " << binInde);
      this->sumForcesIn1Bin(partikelListe, shadowArray, partikelIndex, range1);
    }

    /// use a faster algorithm to sum forces in one or two raster
    /// cells which relies on the particles in any cell to be
    /// separated in two consecutive runs of moving and non-moving
    /// particles.

    void sumForcesIn2Bins(second_argument_type const &partikelListe, 
                          result_type &shadowArray,
                          PartikelIndex const &partikelIndex,
                          BinTriplet const &range,
			  BinTriplet const &otherRange) {
      // loop over particles in *other* bin
      // both ranges non-empty
      assert(range.offset < range.endBoundary);
      assert(otherRange.offset < otherRange.endBoundary);

      // at least one has some moving particles
      assert(range.offset < range.endMoving
	     or otherRange.offset < otherRange.endMoving);
      DEB_SPH("sumForcesIn2Bins: " << length(range) << " particles with "
              << length(otherRange) << "  particles");

      pv_type summep, summeo;
      unsigned oi = otherRange.offset;

#ifdef SUM_FORCES_2BINS_FAST
      // better: combine moving(2) with all(1), then boundary(2) with moving(1)
      // combine moving(2) with all(1)
      for ( ; oi < otherRange.endMoving; ++oi) {
        unsigned const parto = partikelIndex.getPartikelIndex(oi);
	particle_type const &o = partikelListe[parto];

	// loop over MOVING local particles
	unsigned li = range.offset;
	for ( ; li < range.endMoving; ++li) {
	  assert(li != oi);
	  unsigned const partp = partikelIndex.getPartikelIndex(li);
	  particle_type const &p = partikelListe[partp];
          // 	  assert(not p.isBoundary());
	  if (pm.particleRelation(p, o, summep, summeo)) {
            shadowArray[partp] += summep;
            shadowArray[parto] += summeo;
          }
	}
	// loop over BOUNDARY local particles
	for (; li < range.endBoundary; ++li) {
	  assert(li != oi);
	  unsigned const partp = partikelIndex.getPartikelIndex(li);
	  particle_type const &p = partikelListe[partp];
          // 	  assert(p.isBoundary());
          // \todo use particleRelationBoundary shortcut
	  if (pm.particleRelation(o, p, summeo)) {
            shadowArray[parto] += summeo;
          }
	}
      }

      // combine boundary(2) with moving(1)
      for ( ; oi < otherRange.endBoundary; ++oi) {
	particle_type const &o = partikelListe[partikelIndex.getPartikelIndex(oi)];
	assert(o.isBoundary());
	// loop over local particles
	for (unsigned li = range.offset; li < range.endMoving; ++li) {
	  assert(li != oi);
	  unsigned const partp = partikelIndex.getPartikelIndex(li);
	  particle_type const &p = partikelListe[partp];
	  assert(not p.isBoundary());
          // \todo use particleRelationBoundary shortcut
	  if (pm.particleRelation(p, o, summep)) {
            shadowArray[partp] += summep;
          }
	}
      }
#else
#warning slow
      // easy: combine all with all
      for ( ; oi < otherRange.endBoundary; ++oi) {
        unsigned const parto = partikelIndex.getPartikelIndex(oi);
        particle_type const &o = partikelListe[parto];
        if (o.isBoundary()) continue;
	// loop over local particles
	unsigned li = range.offset;
	for ( ; li < range.endBoundary; ++li) {
	  assert(li != oi);
          unsigned const partp = partikelIndex.getPartikelIndex(li);
          particle_type const &p = partikelListe[partp];
	  if (pm.particleRelation(o, p, summeo, summep)) {
            if (not p.isBoundary()) {
              assert(shadowArray[partp].id() == 0 or shadowArray[partp].id() == p.id());
              shadowArray[partp].id() = p.id();
              shadowArray[partp] += summep;
            } else {
              //               assert(summep.isZero());
            }
            assert(shadowArray[parto].id() == 0 or shadowArray[parto].id() == o.id());
            shadowArray[parto].id() = o.id();
            shadowArray[parto] += summeo;
          }
	}
      }
#endif
    }


    void sumForcesIn2BinCoords(second_argument_type const &partikelListe,
                                result_type &shadowArray,
                                PartikelIndex const &partikelIndex,
                                Coord const &gridPos1,
                                Coord const &gridPos2) {
      BinTriplet range1, range2;

      size_t const binIndex = this->pm.m_hashedInfo.m_index(gridPos1);
      partikelIndex.getBinInfo(binIndex, range1);

      size_t const otherIndexFlat = this->pm.m_hashedInfo.m_index(gridPos2);
      partikelIndex.getBinInfo(otherIndexFlat, range2);

      if ((range1.offset < range1.endMoving
           and range2.offset < range2.endBoundary)
          or (range2.offset < range2.endMoving
              and range1.offset < range1.endBoundary)) {
        DEB_SPH("sumForcesIn2BinIndices: " << gridPos1 << " " << gridPos2);
        DEB_SPH("sumForcesIn2BinIndices: " << binIndex << " " << otherIndexFlat);
        this->sumForcesIn2Bins(partikelListe, shadowArray, partikelIndex, range1, range2);
      }
    }

    result_type const &operator()(SphADataType const,
                                  value_type const &partikelListe) {

      if (m_shadowArray.size() != partikelListe.size()) {
        m_shadowArray.resize(partikelListe.size());
      }

      m_shadowArray = ParticleVariable();

      int const nPasses = this->pm.posNachbarArray.size();

#ifdef SPH_COUNT_RELATIONS
#ifdef _OPENMP
#pragma omp parallel
#endif
      {
        numPPrelations = 0;
        //      numPPindextests = 0;
        numPPdisttests = 0;
      }
#endif

#ifdef SPH_RENDER_FORCES
      pm.partikelDistances.clear();
      pm.partikelForces.clear();
#endif

      timeIndex.start();

      unsigned const curNumThreads = omp_get_max_threads();
      unsigned const numThreadsForSort = min(curNumThreads, pm.config.numCSortThreads);
//       { cant'do
        omp_set_num_threads(numThreadsForSort);

        SortedPartikelIndex<SPHPartikelArray>
          partikelIndex(pm.config.sortedIndex, pm.m_hashedInfo.m_index, 
                        pm.config.gebiet, pm.config.kernelReachH);
        partikelIndex.update(partikelListe);

        omp_set_num_threads(curNumThreads);
//       } cant'do

      timeIndex.stop();
      if (pm.config.detailedTimings) {
        partikelIndex.writeXMLTimings(pm.sphout);
        pm.sphout << "<update-index time='" << timeIndex.lastDiff() << "' threads='" << numThreadsForSort << "'/>\n";
      }

      m_numParticlesOut = partikelIndex.numParticlesOut();

      timeSummation.start();

      // "even" or "odd" pass (ndim-dimensional)
      for (int oddityIndex = 0; oddityIndex < nPasses; ++oddityIndex) {

        Coord const &stencilOffset = pm.posNachbarArray[oddityIndex];

        DEB_SPH("ff2: oddity (stencil offset): " << stencilOffset);

        long const von = 0, bis = static_power<ndim>(1 << (pm.config.hashFunctionM - 1));

        DEB_SPH("ff2: par loop len: " << bis);

        assert(static_power<ndim>(2) * bis == long(pm.m_hashedInfo.m_index.size));

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
        // #pragma omp parallel for schedule(dynamic)
        // #pragma omp parallel for schedule(static, 1)
        // #pragma omp parallel for
#endif
        for (long i = von; i < bis; ++i) {

          // compute i-th grid position in D dimensions:
          // - compute i-th grid index in grid halved in alle dimensions
          Coord gridOffset = gridIndices[i];
          DEB_SPH("i: " << i << " zaehlerGrid: " << gridOffset);

          // - add stencil offset (depending on pass odditiy)
          gridOffset += stencilOffset;
          // - compute single integers index
          size_t const binIndex = pm.m_hashedInfo.m_index(gridOffset);

          DEB_SPH("gridOffset: " << gridOffset << " binIndex: " << binIndex);

          std::valarray<BinTriplet > ranges(static_power<ndim>(2));

          partikelIndex.getBinInfo(binIndex, ranges[0]);

          // now compute all combinations of posNachbarArray elements where the
          // scalar product is 0:
          // outer loop over all elements
          for (int ni = 1; ni < nPasses; ++ni) {

            // compute \f$ (g + s) mod 2^{mD}\f$
            long const otherIndexFlat = pm.m_hashedInfo.m_index(gridOffset + pm.posNachbarArray[ni]);

            partikelIndex.getBinInfo(otherIndexFlat, ranges[ni]);

            // inner loop over lower-indexed elements
            for (int nni = 0; nni < ni; ++nni) {
              // only act if scalar product not zero,
              // this test is the same because (in a sense) posNachbarArray[i] == i
              if ((nni & ni) == 0) {
                // only call if there are forces to sum (at least particles in both,
                // in one of them moving)
                if ((ranges[nni].offset < ranges[nni].endMoving
                     and ranges[ni].offset < ranges[ni].endBoundary)
                    or (ranges[ni].offset < ranges[ni].endMoving
                        and ranges[nni].offset < ranges[nni].endBoundary)) {
                  sumForcesIn2Bins(partikelListe, m_shadowArray, partikelIndex, ranges[nni], ranges[ni]);
                }
              }
            }

          }

          sumForcesIn1Bin(partikelListe, m_shadowArray, partikelIndex, ranges[0]);

        }

      }

      timeSummation.stop();
      if (pm.config.detailedTimings) {
        pm.sphout << "<summation time='" << timeSummation.lastDiff() << "' type='symmetric-naive'/>\n";
      }

      return m_shadowArray;
    }


    void update(value_type &state, SPHPartikelVariablenArray const &update) {
      pm.updateState(state, update);
    }
  };


  struct ForceFunctorSymMultiBuffer : public ForceFunctorSymmetric {
    typedef ForceFunctorSymmetric Base;
    typedef typename Base::result_type result_type; 
    typedef typename Base::value_type value_type; 

    std::valarray<result_type> tmpShadowArray;
    Timer timeClear, timeReduce;

    ForceFunctorSymMultiBuffer(PartikelManager const &_pm) : 
      ForceFunctorSymmetric(_pm),
      tmpShadowArray(this->pm.posNachbarArray.size())
    {
      cerr << "creating ForceFunctorSymMultiBuffer \n";
    }
    
    ForceFunctorSymMultiBuffer *clone() {
      return new ForceFunctorSymMultiBuffer(*this);
    }
    

    void reduceTmpShadowArray() {
      long const bisReduce = this->m_shadowArray.size();
      int const nPasses = this->pm.posNachbarArray.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (long i = 0; i < bisReduce; ++i) {
        this->m_shadowArray[i] = 0;
        for (int oddityIndex = 0; oddityIndex < nPasses; ++oddityIndex) {
          this->m_shadowArray[i] += tmpShadowArray[oddityIndex][i];
        }
      }
    }

    
    void clearTmpShadowArray() {
      long const bisReduce = this->m_shadowArray.size();
      int const nPasses = this->pm.posNachbarArray.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (long k = 0; k < bisReduce; ++k) {
        for (int oddityIndex = 0; oddityIndex < nPasses; ++oddityIndex) {
          tmpShadowArray[oddityIndex][k] = 0.0;
        }
      }
    }
    

    result_type const &operator()(SphADataType const, 
                                  value_type const &partikelListe) {
      
      int const nPasses = this->pm.posNachbarArray.size();
      
      if (this->m_shadowArray.size() != partikelListe.size()) {
        this->m_shadowArray.resize(partikelListe.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nPasses)
#endif
          for (long k = 0; k < nPasses; ++k) {
            tmpShadowArray[k].resize(partikelListe.size());
          }
      }
      
      this->timeIndex.start();

      unsigned const curNumThreads = omp_get_max_threads();
      unsigned const numThreadsForSort = min(curNumThreads, this->pm.config.numCSortThreads);
//       { cant'do
        omp_set_num_threads(numThreadsForSort);

        SortedPartikelIndex<SPHPartikelArray> 
          partikelIndex(this->pm.config.sortedIndex, this->pm.m_hashedInfo.m_index, 
                        this->pm.config.gebiet, this->pm.config.kernelReachH);
        partikelIndex.update(partikelListe);

        omp_set_num_threads(curNumThreads);
//       } cant'do

      this->timeIndex.stop();
      if (this->pm.config.detailedTimings) {
        partikelIndex.writeXMLTimings(this->pm.sphout);
        this->pm.sphout << "<update-index time='" << 
          this->timeIndex.lastDiff() << "' threads='" << numThreadsForSort << "'/>\n";
      }

      this->m_numParticlesOut = partikelIndex.numParticlesOut();

      this->timeSummation.start();

      timeClear.start();
      clearTmpShadowArray();
      timeClear.stop();
      if (this->pm.config.detailedTimings) {
        this->pm.sphout << "<clear-scratch time='" << 
          timeClear.lastDiff() << "'/>\n";
      }

      {

        long const numGridPosPerPass = static_power<ndim>(1 << (this->pm.config.hashFunctionM - 1));
        long const numGridPosPerPassModMask = numGridPosPerPass - 1;
        long const totalDim = nPasses * numGridPosPerPass;

//         cerr << "numGridPosPerPass: " << numGridPosPerPass << ", nPasses: " << nPasses << " -> " << totalDim << "\n";

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (long cnt = 0; cnt < totalDim; ++cnt) {
	
          // "even" or "odd" pass (ndim-dimensional)
          int const oddityIndex = cnt >> ((this->pm.config.hashFunctionM - 1)*ndim);
          // raster position index
          int const i = cnt & numGridPosPerPassModMask;

//           cerr << "cnt: " << cnt << " -> (" << oddityIndex << ", " << i << ")\n";
          assert(oddityIndex < nPasses);
          assert(i < numGridPosPerPass);
          assert(oddityIndex * numGridPosPerPass + i == cnt);

          // writing into one of nPasses arrays avoids races
          result_type &shadowArray = tmpShadowArray[oddityIndex];
          Coord const &stencilOffset = this->pm.posNachbarArray[oddityIndex];


            // compute i-th grid position in D dimensions:
            // - compute i-th grid index in grid halved in alle dimensions
            Coord gridOffset = this->gridIndices[i];
            DEB_SPH("i: " << i << " zaehlerGrid: " << gridOffset);

            // - add stencil offset (depending on pass odditiy)
            gridOffset += stencilOffset;
            // - compute single integers index
            size_t const binIndex = this->pm.m_hashedInfo.m_index(gridOffset);

            DEB_SPH("gridOffset: " << gridOffset << " binIndex: " << binIndex);
	  
            std::valarray<BinTriplet > ranges(static_power<ndim>(2));

            partikelIndex.getBinInfo(binIndex, ranges[0]);

            // now compute all combinations of posNachbarArray elements where the
            // scalar product is 0:
            // outer loop over all elements
            for (int ni = 1; ni < nPasses; ++ni) {
	    
              // compute \f$ (g + s) mod 2^{mD}\f$ 
              long const otherIndexFlat = 
                this->pm.m_hashedInfo.m_index(gridOffset + this->pm.posNachbarArray[ni]);

              partikelIndex.getBinInfo(otherIndexFlat, ranges[ni]);

              // inner loop over lower-indexed elements
              for (int nni = 0; nni < ni; ++nni) {
                // only act if scalar product not zero,
                // this test is the same because (in a sense) posNachbarArray[i] == i
                if ((nni & ni) == 0) {
                  // only call if there are forces to sum (at least particles in both, 
                  // in one of them moving)
                  if ((ranges[nni].offset < ranges[nni].endMoving
                       and ranges[ni].offset < ranges[ni].endBoundary)
                      or (ranges[ni].offset < ranges[ni].endMoving
                          and ranges[nni].offset < ranges[nni].endBoundary)) {
                    this->sumForcesIn2Bins(partikelListe, shadowArray, partikelIndex, ranges[nni], ranges[ni]);
                  }
                }
              }

            }
	  
            this->sumForcesIn1Bin(partikelListe, shadowArray, partikelIndex, ranges[0]);
	  
        } // for
      } // block
      
      timeReduce.start();
      reduceTmpShadowArray();
      timeReduce.stop();
      if (this->pm.config.detailedTimings) {
        this->pm.sphout << "<reduce-scratch time='" << 
          timeReduce.lastDiff() << "'/>\n";
      }

      this->timeSummation.stop();
      
      if (this->pm.config.detailedTimings) {
        this->pm.sphout << "<summation time='" << 
          this->timeSummation.lastDiff() << "' type='symmetric-multibuffer'/>\n";
      }

      return this->m_shadowArray;
    }

  };

  struct ForceFunctorSymMultiBufferTB : public ForceFunctorSymMultiBuffer {
    typedef ForceFunctorSymMultiBuffer Base;
    typedef typename Base::result_type result_type; 
    typedef typename Base::value_type value_type; 

    ForceFunctorSymMultiBufferTB(PartikelManager const &_pm) : 
      Base(_pm)
    {
      cerr << "creating ForceFunctorSymMultiBufferTB \n";
    }
    
    ForceFunctorSymMultiBufferTB *clone() {
      return new ForceFunctorSymMultiBufferTB(*this);
    }
    

    result_type const &operator()(SphADataType const, 
                                  value_type const &partikelListe) {
      
      int const nPasses = this->pm.posNachbarArray.size();
      
      if (this->m_shadowArray.size() != partikelListe.size()) {
        this->m_shadowArray.resize(partikelListe.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nPasses)
#endif
          for (long k = 0; k < nPasses; ++k) {
            this->tmpShadowArray[k].resize(partikelListe.size());
          }
      }
      
      this->timeIndex.start();
      
      unsigned const curNumThreads = omp_get_max_threads();
      unsigned const numThreadsForSort = min(curNumThreads, this->pm.config.numCSortThreads);
      //       { cant'do
      omp_set_num_threads(numThreadsForSort);
      
      SortedPartikelIndex<SPHPartikelArray> 
        partikelIndex(this->pm.config.sortedIndex, this->pm.m_hashedInfo.m_index, 
                      this->pm.config.gebiet, this->pm.config.kernelReachH);
      partikelIndex.update(partikelListe);
      
      omp_set_num_threads(curNumThreads);
      //       } cant'do
      
      this->timeIndex.stop();
      if (this->pm.config.detailedTimings) {
        partikelIndex.writeXMLTimings(this->pm.sphout);
        this->pm.sphout << "<update-index time='" << 
          this->timeIndex.lastDiff() << "' threads='" << numThreadsForSort << "'/>\n";
      }
      
      this->m_numParticlesOut = partikelIndex.numParticlesOut();

      this->timeSummation.start();

      this->timeClear.start();
      this->clearTmpShadowArray();
      this->timeClear.stop();
      if (this->pm.config.detailedTimings) {
        this->pm.sphout << "<clear-scratch time='" << 
          this->timeClear.lastDiff() << "'/>\n";
      }

      {

        long const numGridPosPerPass = static_power<ndim>(1 << (this->pm.config.hashFunctionM - 1));

//         cerr << "numGridPosPerPass: " << numGridPosPerPass << ", nPasses: " << nPasses << " -> " << totalDim << "\n";
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
          int const nThreads = omp_get_num_threads();
          assert(nThreads % 2 == 0);
          int const myid = omp_get_thread_num();
          // "even" or "odd" pass (ndim-dimensional)
          int const remOddities = max(nPasses / nThreads, 1);
          int const threadsPerOddity = max(nThreads / nPasses, 1);
          int const oddThreadId = myid % threadsPerOddity;
          
          long const myRange = numGridPosPerPass / threadsPerOddity;
          long const myStart = myRange * oddThreadId;
          long const myEnd = myRange * (oddThreadId+1);

// #pragma omp critical
//           {
//             cout << "<sym-tb myid='" << myid << "' threads='" << nThreads << "' remOddities='" << remOddities << "'"
//               " threadsPerOddity='" << threadsPerOddity << "' oddThreadId='" << oddThreadId << "' myRange='" << myRange << "'"
//               " passes='" << nPasses << "'/>\n";
//           }

          for (long rodd = 0; rodd < remOddities; ++rodd) {
            int const oddityIndex = ((myid/threadsPerOddity) % nPasses) + rodd*nThreads;
            // writing into one of nPasses arrays avoids races
            result_type &shadowArray = this->tmpShadowArray[oddityIndex];
            Coord const &stencilOffset = this->pm.posNachbarArray[oddityIndex];

            for (long i = myStart; i < myEnd; ++i) {
	
              //           cerr << "cnt: " << cnt << " -> (" << oddityIndex << ", " << i << ")\n";
              assert(oddityIndex < nPasses);
              assert(i < numGridPosPerPass);

              // compute i-th grid position in D dimensions:
              // - compute i-th grid index in grid halved in alle dimensions
              Coord gridOffset = this->gridIndices[i];
              DEB_SPH("i: " << i << " zaehlerGrid: " << gridOffset);
              
              // - add stencil offset (depending on pass odditiy)
              gridOffset += stencilOffset;
              // - compute single integers index
              size_t const binIndex = this->pm.m_hashedInfo.m_index(gridOffset);
              
              DEB_SPH("gridOffset: " << gridOffset << " binIndex: " << binIndex);
              
              std::valarray<BinTriplet > ranges(static_power<ndim>(2));
              
              partikelIndex.getBinInfo(binIndex, ranges[0]);
              
              // now compute all combinations of posNachbarArray elements where the
              // scalar product is 0:
              // outer loop over all elements
              for (int ni = 1; ni < nPasses; ++ni) {
                
              // compute \f$ (g + s) mod 2^{mD}\f$ 
              long const otherIndexFlat = 
                this->pm.m_hashedInfo.m_index(gridOffset + this->pm.posNachbarArray[ni]);

              partikelIndex.getBinInfo(otherIndexFlat, ranges[ni]);

              // inner loop over lower-indexed elements
              for (int nni = 0; nni < ni; ++nni) {
                // only act if scalar product not zero,
                // this test is the same because (in a sense) posNachbarArray[i] == i
                if ((nni & ni) == 0) {
                  // only call if there are forces to sum (at least particles in both, 
                  // in one of them moving)
                  if ((ranges[nni].offset < ranges[nni].endMoving
                       and ranges[ni].offset < ranges[ni].endBoundary)
                      or (ranges[ni].offset < ranges[ni].endMoving
                          and ranges[nni].offset < ranges[nni].endBoundary)) {
                    this->sumForcesIn2Bins(partikelListe, shadowArray, partikelIndex, ranges[nni], ranges[ni]);
                  }
                }
              }
              
              }
              
              this->sumForcesIn1Bin(partikelListe, shadowArray, partikelIndex, ranges[0]);
	  
            } // for
          } // for
        } // parallel block
      } // block
      
      this->timeReduce.start();
      this->reduceTmpShadowArray();
      this->timeReduce.stop();
      if (this->pm.config.detailedTimings) {
        this->pm.sphout << "<reduce-scratch time='" << 
          this->timeReduce.lastDiff() << "'/>\n";
      }

      this->timeSummation.stop();
      
      if (this->pm.config.detailedTimings) {
        this->pm.sphout << "<summation time='" << 
          this->timeSummation.lastDiff() << "' type='symmetric-multibuffer-omp-thread-binding'/>\n";
      }

      return this->m_shadowArray;
    }
    
  };
    

    struct ForceFunctorCube2D : public ForceFunctorSymmetric {
      typedef ForceFunctorSymmetric Base;
      typedef typename Base::result_type result_type; 
      typedef typename Base::value_type value_type; 
      
      ForceFunctorCube2D(PartikelManager const &_pm) : 
        ForceFunctorSymmetric(_pm)
      {
        cerr << "creating ForceFunctorCube2D\n";
      }
      
      ForceFunctorCube2D *clone() {
        return new ForceFunctorCube2D(*this);
      }
      
      
      result_type const &operator()(SphADataType const, value_type const &partikelListe) {
        
        if (this->m_shadowArray.size() != partikelListe.size()) {
          this->m_shadowArray.resize(partikelListe.size());
        }
        
        this->m_shadowArray = ParticleVariable();

        this->timeIndex.start();
        
        unsigned const curNumThreads = omp_get_max_threads();
        unsigned const numThreadsForSort = min(curNumThreads, this->pm.config.numCSortThreads);
        //       { cant'do
        omp_set_num_threads(numThreadsForSort);
        
        SortedPartikelIndex<SPHPartikelArray> 
          partikelIndex(this->pm.config.sortedIndex, this->pm.m_hashedInfo.m_index, 
                        this->pm.config.gebiet, this->pm.config.kernelReachH);
        partikelIndex.update(partikelListe);
        
        omp_set_num_threads(curNumThreads);
      //       } cant'do
        
        this->timeIndex.stop();
        if (this->pm.config.detailedTimings) {
          partikelIndex.writeXMLTimings(this->pm.sphout);
          this->pm.sphout << "<update-index time='" << 
            this->timeIndex.lastDiff() << "' threads='" << numThreadsForSort << "'/>\n";
        }
        
        this->m_numParticlesOut = partikelIndex.numParticlesOut();

        {
          int const dim = 0;
          Coord gridPos, axis;
          axis[dim] = 1;
          DEB_SPH("parallel iteration over axis: " << axis << ":  "
                  << this->pm.m_hashedInfo.m_index.hashFunktion.m_dim);
          
#ifdef _OPENMP
#pragma omp parallel for private(gridPos) schedule(dynamic)
#endif
          for (unsigned i = 0; i < this->pm.m_hashedInfo.m_index.hashFunktion.m_dim; ++i) {
            gridPos[1] = i;
            gridPos[0] = 0;
            BinTriplet range1 = this->getRange(partikelIndex, gridPos), range2;
            for (unsigned j = 0; j < this->pm.m_hashedInfo.m_index.hashFunktion.m_dim; ++j) {
              gridPos[0] = j;
              
              range2 = this->getRange(partikelIndex, gridPos+axis);
              this->sumForcesIn2Bins(partikelListe, this->m_shadowArray, partikelIndex, range1, range2);

              range1 = range2;
            }
          }
        }
        {
          int const dim = 1;
          Coord axis, gridPos, splitAxis;
          axis[dim] = 1;
          splitAxis[0] = 1;
          
          DEB_SPH("parallel iteration over axis: " << axis << ":  "
                  << this->pm.m_hashedInfo.m_index.hashFunktion.m_dim);
          
          for (unsigned odd = 0; odd < 2; ++odd) {
#ifdef _OPENMP
#pragma omp parallel for private(gridPos) schedule(dynamic)
#endif
            for (unsigned i = 0; i < this->pm.m_hashedInfo.m_index.hashFunktion.m_dim; i += 2) {
              gridPos[0] = i + odd;
              gridPos[1] = 0;
              BinTriplet range1 = this->getRange(partikelIndex, gridPos), range2,
                rangeLeft1 = this->getRange(partikelIndex, gridPos-axis+splitAxis), 
                rangeLeft2 = this->getRange(partikelIndex, gridPos+splitAxis), rangeLeft3;
              for (unsigned j = 0; j < this->pm.m_hashedInfo.m_index.hashFunktion.m_dim; ++j) {
                gridPos[1] = j;
                
                range2 = this->getRange(partikelIndex, gridPos+axis);
                this->sumForcesIn1Bin(partikelListe, this->m_shadowArray, partikelIndex, range1);
                this->sumForcesIn2Bins(partikelListe, this->m_shadowArray, partikelIndex, 
                                            range1, range2);
                rangeLeft3 = this->getRange(partikelIndex, gridPos+axis+splitAxis);
                this->sumForcesIn2Bins(partikelListe, this->m_shadowArray, partikelIndex, 
                                       range1, rangeLeft1);
                this->sumForcesIn2Bins(partikelListe, this->m_shadowArray, partikelIndex, 
                                       range1, rangeLeft3);

                range1 = range2;
                rangeLeft1 = rangeLeft2;
                rangeLeft2 = rangeLeft3;
              }
            }
          }
        }
      
        return this->m_shadowArray;
      }
    
  };

    struct ForceFunctorCube3D : public ForceFunctorSymmetric {
      typedef ForceFunctorSymmetric Base;
      typedef typename Base::result_type result_type; 
      typedef typename Base::value_type value_type; 
      
      ForceFunctorCube3D(PartikelManager const &_pm) : 
        ForceFunctorSymmetric(_pm)
      {
        cerr << "creating ForceFunctorCube3D\n";
      }
      
      ForceFunctorCube3D *clone() {
        return new ForceFunctorCube3D(*this);
      }
      
      
      result_type const &operator()(SphADataType const, value_type const &partikelListe) {
        
        if (this->m_shadowArray.size() != partikelListe.size()) {
          this->m_shadowArray.resize(partikelListe.size());
        }
        
        this->m_shadowArray = ParticleVariable();

        this->timeIndex.start();
        
        unsigned const curNumThreads = omp_get_max_threads();
        unsigned const numThreadsForSort = min(curNumThreads, this->pm.config.numCSortThreads);
        //       { cant'do
        omp_set_num_threads(numThreadsForSort);
        
        SortedPartikelIndex<SPHPartikelArray> 
          partikelIndex(this->pm.config.sortedIndex, this->pm.m_hashedInfo.m_index, 
                        this->pm.config.gebiet, this->pm.config.kernelReachH);
        partikelIndex.update(partikelListe);
        
        omp_set_num_threads(curNumThreads);
      //       } cant'do
        
        this->timeIndex.stop();
        if (this->pm.config.detailedTimings) {
          partikelIndex.writeXMLTimings(this->pm.sphout);
          this->pm.sphout << "<update-index time='" << 
            this->timeIndex.lastDiff() << "' threads='" << numThreadsForSort << "'/>\n";
        }
        
        this->m_numParticlesOut = partikelIndex.numParticlesOut();

        {
          Coord gridPos, axis1, axis2;
          axis1[0] = 1;
          axis2[2] = 1;
          DEB_SPH("parallel iteration over axes: " << axis1 << ", " << axis2 << ":  "
                  << this->pm.m_hashedInfo.m_index.hashFunktion.m_dim);
          
#ifdef _OPENMP
#pragma omp parallel for private(gridPos) schedule(dynamic)
#endif
          for (unsigned i = 0; i < this->pm.m_hashedInfo.m_index.hashFunktion.m_dim; ++i) {
            gridPos[1] = i;
            for (unsigned j = 0; j < this->pm.m_hashedInfo.m_index.hashFunktion.m_dim; ++j) {
              gridPos[0] = j;
              for (unsigned k = 0; k < this->pm.m_hashedInfo.m_index.hashFunktion.m_dim; ++k) {
                gridPos[2] = k;
              
                this->sumForcesIn2BinCoords(partikelListe, this->m_shadowArray, partikelIndex, 
                                            gridPos, gridPos+axis1+axis2);
                this->sumForcesIn2BinCoords(partikelListe, this->m_shadowArray, partikelIndex, 
                                            gridPos, gridPos+axis1);
                this->sumForcesIn2BinCoords(partikelListe, this->m_shadowArray, partikelIndex, 
                                            gridPos, gridPos+axis2);
                this->sumForcesIn2BinCoords(partikelListe, this->m_shadowArray, partikelIndex, 
                                            gridPos+axis1, gridPos+axis2);
              }
            }
          }
        }
        {
          Coord gridPos, axis1, axis2;
          axis1[0] = 1;
          axis2[1] = 1;
          DEB_SPH("parallel iteration over axes: " << axis1 << ", " << axis2 << ":  "
                  << this->pm.m_hashedInfo.m_index.hashFunktion.m_dim);
          
#ifdef _OPENMP
#pragma omp parallel for private(gridPos) schedule(dynamic)
#endif
          for (unsigned i = 0; i < this->pm.m_hashedInfo.m_index.hashFunktion.m_dim; ++i) {
            gridPos[2] = i;
            for (unsigned j = 0; j < this->pm.m_hashedInfo.m_index.hashFunktion.m_dim; ++j) {
              gridPos[0] = j;
              for (unsigned k = 0; k < this->pm.m_hashedInfo.m_index.hashFunktion.m_dim; ++k) {
                gridPos[1] = k;
              
                this->sumForcesIn2BinCoords(partikelListe, this->m_shadowArray, partikelIndex, 
                                            gridPos, gridPos+axis1+axis2);
                this->sumForcesIn2BinCoords(partikelListe, this->m_shadowArray, partikelIndex, 
                                            gridPos, gridPos+axis1);
                this->sumForcesIn2BinCoords(partikelListe, this->m_shadowArray, partikelIndex, 
                                            gridPos, gridPos+axis2);
                this->sumForcesIn2BinCoords(partikelListe, this->m_shadowArray, partikelIndex, 
                                            gridPos+axis1, gridPos+axis2);
              }
            }
          }
        }
        {
          Coord gridPos, axis1, axis2,splitAxis;
          axis1[1] = 1;
          axis2[2] = 1;
          DEB_SPH("parallel iteration over axes: " << axis1 << ", " << axis2 << ":  "
                  << this->pm.m_hashedInfo.m_index.hashFunktion.m_dim);

          splitAxis[0] = 1;
          
          for (unsigned odd = 0; odd < 2; ++odd) {
#ifdef _OPENMP
#pragma omp parallel for private(gridPos) schedule(dynamic)
#endif
            for (unsigned i = 0; i < this->pm.m_hashedInfo.m_index.hashFunktion.m_dim; i += 2) {
              gridPos[0] = i + odd;
              for (unsigned j = 0; j < this->pm.m_hashedInfo.m_index.hashFunktion.m_dim; ++j) {
                gridPos[1] = j;
                for (unsigned k = 0; k < this->pm.m_hashedInfo.m_index.hashFunktion.m_dim; ++k) {
                  gridPos[2] = k;
                  
                  this->sumForcesIn2BinCoords(partikelListe, this->m_shadowArray, partikelIndex, 
                                              gridPos, gridPos+axis1+axis2);
                  this->sumForcesIn2BinCoords(partikelListe, this->m_shadowArray, partikelIndex, 
                                              gridPos, gridPos+axis1);
                  this->sumForcesIn2BinCoords(partikelListe, this->m_shadowArray, partikelIndex, 
                                              gridPos, gridPos+axis2);
                  this->sumForcesIn2BinCoords(partikelListe, this->m_shadowArray, partikelIndex, 
                                              gridPos+axis1, gridPos+axis2);

                  this->sumForcesIn1BinCoords(partikelListe, this->m_shadowArray, partikelIndex, gridPos);

                  this->sumForcesIn2BinCoords(partikelListe, this->m_shadowArray, partikelIndex, 
                                              gridPos, gridPos+axis1+axis2+splitAxis);
                  this->sumForcesIn2BinCoords(partikelListe, this->m_shadowArray, partikelIndex, 
                                              gridPos, gridPos-axis1+axis2+splitAxis);
                  this->sumForcesIn2BinCoords(partikelListe, this->m_shadowArray, partikelIndex, 
                                              gridPos, gridPos+axis1-axis2+splitAxis);
                  this->sumForcesIn2BinCoords(partikelListe, this->m_shadowArray, partikelIndex, 
                                              gridPos, gridPos-axis1-axis2+splitAxis);

                }
              }
            }
          }
        }
      
        return this->m_shadowArray;
      }
    
  };

};

template<class A, class B>
typename PartikelManager<A,B>::ForceFunctorInterface *
PartikelManager<A,B>::ForceFunctorInterface::makeForceFunctor(std::string const &name, PartikelManager const &pm) {
  if (name.compare("naive") == 0) {
    return new ForceFunctorPerParticle(pm);
  } else if (name.compare("sym-mb") == 0) {
    return new ForceFunctorSymMultiBuffer(pm);
  } else if (name.compare("sym-mb-tb") == 0) {
    return new ForceFunctorSymMultiBufferTB(pm);
  } else if (name.compare("symmetric") == 0) {
    return new ForceFunctorSymmetric(pm);
  } else if (name.compare("cube2d") == 0) {
    if (ndim == 2) {
      return new ForceFunctorCube2D(pm);
    } else {
      cerr << "sph: error: cube2d summation only available in sph-2d" << endl;
    }
  } else if (name.compare("cube3d") == 0) {
    if (ndim == 3) {
      return new ForceFunctorCube3D(pm);
    } else {
      cerr << "sph: error: cube3d summation only available in sph-3d" << endl;
    }
  } 
  cerr << "sph: error: no such summation algorithm: '" << name << "'" << endl;
  throw UnknownSummationAlgorithm();
  return 0;
}

#ifdef SPH_RENDER
  /// struct RenderInfo stores information copied over to the
  /// rendering thread, such that the simulation can continue
  /// immediately.
  struct RenderInfo {
    SPHConfig<SphADataType> const &config;
    double t, timePerStep;
    size_t numThreads, step, frame, numRendered;
    size_t npartikel, numParticlesOut;
#ifdef SPH_RENDER_FORCES
    std::vector<PartikelForce> partikelForces, partikelDistances;
#endif
    std::vector<Partikel<SphADataType> > partikelRenderPuffer;
    double minColor, maxColor, minScale, maxScale, colorRange;
    time_t eta;
    RenderInfo(SPHConfig<SphADataType> const &_config) : 
      config(_config), timePerStep(), step(), colorRange() 
    {}
  };
#endif


/// The template class SPH combines all the parts needed for a SPH
/// simulation run.
/// \tparam T numerical data-type, usually double
/// \tparam PartikelManager The (instantiated) type of the PartikelManager
/// \tparam Integrator the Type of the integrator

template<class PartikelManager, 
	 class Integrator>
struct SPH {
  typedef typename PartikelManager::Config Config;
  typedef typename Config::value_type value_type;
  typedef typename PartikelManager::particle_type particle_type;

#ifdef SPH_RENDER
  typedef DoubleVector (SPH::*ParticleVectorFunction)(particle_type const &p) const;
  typedef double (SPH::*ParticleColorFunction)(particle_type const &p) const;
  typedef Color (SPH::*ColorMapFunction)(double) const;
  typedef void (SPH::*ParticleReprFunction)(particle_type const &p) const;
  typedef double (SPH::*ParticleSizeFunction)(particle_type const &p) const;

  ParticleVectorFunction particleVectorFunction;
  ParticleColorFunction  particleColorFunction;
  ColorMapFunction       colorMapFunction;
  ParticleReprFunction   particleReprFunction;
  ParticleSizeFunction   particleSizeFunction;
#endif

  Timer timeStep, timeOutOfStep, timeInit, timeSave, timeLoad;

  /// learning rate alpha in rendering smoothing estimations (eta,
  /// scale top and botton in rendering)
  Config const config;
  double const tsEstAlpha;
  value_type simulationTime;
  double deltat;
  double timeLastStep, timePerStepEstimate, timeOutOfLastStep, endLastStep;
  unsigned step, m_frameRendered;
  std::vector<Quader<Point<double, ndim> > > gebietInitVec;
  std::ostream &sphout;
  PartikelManager partikelManager;
  typedef typename Integrator::ForceFunctor ForceFunctor;
  ForceFunctor forceFunctor;
  Integrator integrator;
#ifdef SPH_RENDER
  Mutex partikelBufferLock;
  RenderInfo renderInfo;
#endif
  size_t simstep;
  size_t savedFrame;
  long numThreads;
  H5PartWriter<SPHPartikelArray> *h5PartWriter;
  Sensor<value_type> *sensors;
  std::valarray<Color> colorMap;
  int program_argc;
  char **program_argv;

  SPH(std::ostream &logstream, Config const &_config = Config()) :
#ifdef SPH_RENDER
    particleVectorFunction(),
    particleColorFunction(),
    colorMapFunction(),
    particleReprFunction(),
#endif
    timeOutOfStep(true),
    timeInit(true),
    config(_config),
    tsEstAlpha(1e-2),
    simulationTime(),
    deltat(config.timestepDefault),
    timeLastStep(), 
    timePerStepEstimate(),
    timeOutOfLastStep(), 
    endLastStep(omp_get_wtime()),
    step(),
    m_frameRendered(),
    gebietInitVec(),
    sphout(logstream),
    partikelManager(this->config, sphout),
    forceFunctor(partikelManager),
    integrator(forceFunctor, config.integratorName),
#ifdef SPH_RENDER
    renderInfo(config),
#endif
    simstep(),
    savedFrame(),
    numThreads(config.numThreads),
    h5PartWriter(),
    sensors(),
    program_argc(),
    program_argv()
  {
    if (config.resultFormat.find("hdf5") != string::npos
        and config.saveRate > 0) {
      h5PartWriter = new H5PartWriter<SPHPartikelArray>(config.resultFile + ".h5", config.saveFields);
    }
    if (config.numSensors) {
      {
        char *sdata = new char[sizeof(Sensor<value_type>)*config.numSensors];
        memset(sdata, 0, sizeof(Sensor<value_type>)*config.numSensors);
        sensors = (Sensor<value_type>*) sdata;
      }
      std::ostream *defaultSensorLog = &sphout;
      // if (not config.sensorLogfile.empy()) {
      //   // FIXME: implement per sensor log file
      //   if sprintf(config.sensorLogfile, 0) == config.sensorLogfile...
      // }
      for(unsigned i = 0; i < config.numSensors; ++i) {
        std::ostream *sensorLog = 0;
        if (config.sensors[i].m_log) {
          sensorLog = defaultSensorLog;
        }
        new (sensors + i) Sensor<value_type>(config, config.sensors[i], i, sensorLog);
      }
    }
#ifdef SPH_RENDER
    particleVectorFunction = &SPH::particleVectorVelocity;
    particleColorFunction = &SPH::particleColorDensity;
    colorMapFunction = &SPH::colorMapGradient;
    particleReprFunction = &SPH::particleReprSolidSphere;
    particleSizeFunction = &SPH::particleRadiusVolume;
#endif
    timeInit.stop();
  }

  ~SPH() {
    if (h5PartWriter) {
      h5PartWriter->closeFile();
      delete h5PartWriter;
      h5PartWriter = 0;
    }
    if (sensors) {
      char *addr = (char*) sensors;
      delete[] addr;
      sensors = 0;
    }
  }

#ifdef SPH_RENDER
  void initReadColors() {
    string fname = config.colorMapFile;
    colorMap.resize(config.colorMapSize);

    Color white;
    white = 1;
    colorMap = white;
    std::ifstream colorin(fname.c_str());
    if (colorin.fail()) {
      cerr << "error: failed to open colors.txt file: " << fname << ": " << strerror(errno) << "\n";
    }
    unsigned ncolor = 0;
    for ( ; not colorin.fail() 
            and ncolor < colorMap.size(); ++ncolor) {
      Color c;
      c = 1;
      colorin >> c;
      if (not colorin) {
	break;
      }
      colorMap[ncolor] = c;
    }
    unsigned const colorsFound = ncolor;
    if (colorsFound) {
      for ( ; ncolor < colorMap.size(); ++ncolor) {
        colorMap[ncolor] = colorMap[ncolor % colorsFound];
      }
    }
    sphout << "<color-map size='" << colorMap.size() << "' num-colors='" << colorsFound << "'/>\n";
  }
#endif

  void init(int argc, char *argv[]) {
    timeInit.start();
#ifdef SPH_RENDER
    if (config.renderEvery) {
      initReadColors();
    }
#endif
    program_argc = argc;
    program_argv = argv;
    if (not config.initialFile.empty()) {
      sphout << "<load type='vtk' filename='" << config.initialFile << "'>\n";
      timeLoad.start();
      std::valarray<particle_type> neuePartikel;
      VTKReader<std::valarray<particle_type> > 
	reader(config.initialFile, config.fieldsToLoad);
      reader.readVTK(neuePartikel);
      partikelManager.addParticles(neuePartikel);
      timeLoad.stop();
      sphout << "<load n='" << neuePartikel.size() << "'"
        " time='" << timeLoad.lastDiff() << "'/>\n";
      sphout << "</load>\n";
    } else {
      cerr << "error: no initial settings file given\n";
      exit(5);
    }

    timeInit.stop();
  }

  time_t eta() const {
    return Time() + (config.tmax - simulationTime) / deltat * timePerStepEstimate;
  }

  void doStep(value_type const timestep = 1e-3) {
#ifdef SPH_RENDER_FORCES
    partikelBufferLock.lock();
#endif

    DEB_SPH("doStep: dt " << timestep << " step " << step << " t " << simulationTime << " s");

#ifdef SPH_COUNT_RELATIONS
    size_t tot_numrel = 0;
//     size_t tot_numindex = 0;
    size_t tot_numdist = 0;
#endif

#ifdef SPH_PARTICLE_MAX_VELOCITY
    double tot_maxVelocity = 0, tot_maxRelativeVelocity = 0;
#endif

#ifdef SPH_SHOW_LOAD_AVERAGE
    double avgPercentage = 0, streuung = 0;
#endif

#if defined SPH_COUNT_RELATIONS || defined SPH_PARTICLE_MAX_VELOCITY
# ifdef _OPENMP
# pragma omp parallel
# endif
    {
# ifdef _OPENMP
# pragma omp critical
# endif
      {
#ifdef SPH_COUNT_RELATIONS
	tot_numrel += numPPrelations;
        // 	tot_numindex += numPPindextests;
	tot_numdist += numPPdisttests;
#endif
#ifdef SPH_PARTICLE_MAX_VELOCITY
        tot_maxVelocity = max(tot_maxVelocity, maxVelocity);
        maxVelocity = 0;
        tot_maxRelativeVelocity = max(tot_maxRelativeVelocity, maxRelativeVelocity);
        maxRelativeVelocity = 0;
#endif
      }
    }
#endif

    timeOutOfStep.stop();
    timeStep.start();

#if SPH_AD == 1
    assert(timestep.diff(0) == 0);
    assert(timestep.diff(1) == 0);
    assert(timestep.diff(2) == 0);
#ifdef SPH_AD_CLEAR_D_POSITION_EACH_STEP
    long const np = partikelManager.numParticles();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long i = partikelManager.numFrozenParticles; i < np; ++i) {
      particle_type &p = partikelManager.partikelArray[i];
      Vector &pos = p.position();
      for(unsigned k= 0; k < ndim; ++k) {
        assert(pos[k].diff().sum() == 0);
        pos[k].diff() = 0;
      }
    }
#endif
#endif

// #ifndef SPH_USE_INTEGRATOR
//     SPHPartikelVariablenArray update = 
//       forceFunctor(timestep, partikelManager.partikelArray);
// #else
    SPHPartikelVariablenArray update = 
      integrator(timestep, partikelManager.partikelArray);
// #endif

    partikelManager.m_numParticlesOut = forceFunctor.numParticlesOut();
    
    for(unsigned i = 0; i < config.numSensors; ++i) {
      Sensor<value_type> &sensor = sensors[i];
      sensor.m_currentTimestep = deltat;
      sensor.m_currentTime = simulationTime;
    }

    long const np = partikelManager.numParticles();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (long i = partikelManager.numFrozenParticles; i < np; ++i) {
      for(unsigned i = 0; i < config.numSensors; ++i) {
        Sensor<value_type> &sensor = sensors[i];
        particle_type &p = partikelManager.partikelArray[i];
        if (not p.isOut()) {
          sensor.processParticle(p);
        }
      }
    }

#ifdef SPH_RENDER
#ifndef SPH_RENDER_FORCES
    partikelBufferLock.lock();
#endif
#endif
#ifndef SPH_USE_INTEGRATOR
    for (size_t i = 0; i < partikelManager.numFrozenParticles; ++i) {
      Particle &p = partikelManager.partikelArray[i];
      assert(p.isBoundary());
      assert(update[i].id() == 0);
      assert(update[i].isZero());
    }
    for (size_t i = partikelManager.numFrozenParticles; 
         i < partikelManager.numParticles(); ++i) {
      Particle &p = partikelManager.partikelArray[i];
      assert(not p.isBoundary());
      if (p.isOut()) 
        continue;
      assert(p.id() == update[i].id());
      p.werte() += update[i] * ParticleVariable(timestep);
      partikelManager.m_hashedInfo.updatePartikelIndex(p);
      partikelManager.updatePressure(p);
    }
#else
    integrator.update(partikelManager.partikelArray, update);
#endif
#ifdef SPH_RENDER
#ifndef SPH_RENDER_FORCES
    partikelBufferLock.unlock();
#endif
#endif
    simulationTime += timestep;
    ++step;

    timeStep.stop();
    timeLastStep  = timeStep.lastDiff();
    if (step == 1) {
      timePerStepEstimate = timeLastStep;
    } else {
      timePerStepEstimate = timePerStepEstimate*(1-tsEstAlpha) + timeLastStep*tsEstAlpha;
    }
    timeOutOfStep.start();

#ifdef SPH_RENDER_FORCES
    partikelBufferLock.unlock();
#endif

    {
      static StrFTime etaFormat("%d. %b %H:%M:%S");
      // out some info about this last step in an XML tag
#ifdef SPH_COUNT_RELATIONS
      size_t const npartikel = partikelManager.numParticles();
#endif
      size_t const sprec = sphout.precision(); 
      sphout.precision(3);
      sphout << "<s"
        " n='" << std::setw(5) << step << "'"
        " t='" << setw(7) << simulationTime << "' " 
#ifdef SPH_ADAPTIVE_TIMESTEP
	" dt='" << scientific << setw(6) << timestep << "'"
#endif
// 	" time-est='" << setw(6) << timePerStepEstimate << "'"
// 	" eta-s='" << setw(6) << (config.numSteps - step) * timePerStepEstimate << "'"
	" eta='" << etaFormat(eta()) << "'"
	" p='" << setw(6) << partikelManager.numParticlesLive() << "'";
      for(unsigned i = 0; i < config.numSensors; ++i) {
        Sensor<value_type> &sensor = sensors[i];
        sphout << " s" << i << "='" << setw(7) << sensor.m_sensorValue << "'"; 
      }
      sphout << 
 	// " otime='" << setw(6) << timeOutOfStep.lastDiff() << "'"
#ifdef SPH_COUNT_RELATIONS
	" rels='" << setw(10) << double(tot_numrel)/npartikel << "'"
	//       " itests='" << setw(10) << tot_numindex << "'"
	" dtests='" << setw(10) << double(tot_numdist)/(npartikel) << "'"
#ifdef SPH_SHOW_LOAD_AVERAGE
	" avg-load='" << avgPercentage << "'"
	" load-var='" << sqrt(streuung) << "'"
#endif
#endif
#ifdef SPH_PARTICLE_MAX_VELOCITY
	" cfl='" << scientific << setw(6)
             << tot_maxRelativeVelocity * timestep / (0.5 * config.kernelReachH) << "'"
	" max-vel='" << scientific << setw(6) << tot_maxVelocity << "'"
	" max-rel-vel='" << scientific << setw(6) << tot_maxRelativeVelocity << "'"
#endif
	" time='" << setw(6) << timeLastStep << "'"
	"/>\r" << std::flush
	;
      sphout.precision(sprec);
    }

#ifdef SPH_ADAPTIVE_TIMESTEP
    // note: by default, cflNumber == 1 and timestepMax == timestepDefault
    //                      => deltat == timestepDefault 
    // so adaptive timestepping is effectively disabled
    // only when cflNumber is set smaller (and timestepMax higher) deltat will vary
    if (tot_maxVelocity > 0) {
      // the new timestep according to CFL rule
      double const newCFLTimestep = config.cflNumber * (0.5 * config.kernelReachH)
        / max(tot_maxRelativeVelocity, tot_maxVelocity);
      deltat = min(
                   min(newCFLTimestep, deltat*1.1), // timestep grows at most by 10 %
                   config.timestepMax);             // timestep limited above by timestepMax
    }
#endif
  }

#ifdef SPH_RENDER

  void renderRaster() {
    size_t ninterx = floor(config.gebiet.v[0] / config.kernelReachH);
    size_t nintery = floor(config.gebiet.v[1] / config.kernelReachH);
    glBegin(GL_LINES);
    for (size_t i = 0; i < ninterx; ++i) {
      if (i % 10) {
	glColor3f(0.7, 1, 0.7);
      } else {
	glColor3f(0.7, 0, 0);
      }
      glVertex3d(config.gebiet.untenLinks[0], 
		 config.gebiet.untenLinks[1] + i * config.kernelReachH, 
		 -0.01);
      glVertex3d(config.gebiet.obenRechts[0], 
		 config.gebiet.untenLinks[1] + i * config.kernelReachH, 
		 -0.01);
    }
    for (size_t i = 0; i < nintery; ++i) {
      if (i % 10) {
	glColor3f(0.7, 1, 0.7);
      } else {
	glColor3f(0, 0.7, 0);
      }
      glVertex3d(config.gebiet.untenLinks[0] + i * config.kernelReachH, 
		 config.gebiet.untenLinks[1], 
		 -0.01);
      glVertex3d(config.gebiet.untenLinks[0] + i * config.kernelReachH, 
		 config.gebiet.obenRechts[1], 
		 -0.01);
    }
    glEnd();
  }

  void renderPlane(Point3D const &pnt, Point3D const &normVek) {
    glColor3f(1, 1, 1);
    glBegin(GL_POINTS);
    glVertex3dv(&pnt[0]);
    glEnd();
    
    glBegin(GL_LINES);
    {
      // normal vector
      glColor3f(0, 0, 1);
      glVertex3dv(&pnt[0]);
      glVertex3dv(&(pnt + normVek*0.15)[0]);

      Point3D toCross;
      toCross = 1;
      
      // otho1
      Point3D tick1 = normed(cross(toCross, normVek)) * 0.15;
      
      glColor3f(1, 0, 0);
      glVertex3dv(&pnt[0]);
      glVertex3dv(&(pnt + tick1)[0]);
      glVertex3dv(&pnt[0]);
      glVertex3dv(&(pnt - tick1)[0]);
      
      // otho2
      Point3D tick2 = normed(cross(tick1, normVek)) * 0.15;
      
      glColor3f(0, 1, 0);
      glVertex3dv(&pnt[0]);
      glVertex3dv(&(pnt + tick2)[0]);
      glVertex3dv(&pnt[0]);
      glVertex3dv(&(pnt - tick2)[0]);
    }
    glEnd();

  }

  void renderSensors() {
    glCallList(dlIdSensors);
  }

  void doRenderSensors() {
    for(unsigned i = 0; i < config.numSensors; ++i) {
      // sensor-gebiete: blaue linien
      glColor3f(0, 0, 1);
      Sensor<value_type> &sensor = sensors[i];
      switch(sensor.m_sensorDef.m_areaType) {
      case SensorDef::SENSOR_AREA_ORTHOGONAL_RECTANGLE:
        config.sensors[i].m_gebiet.render();
        break;
      case SensorDef::SENSOR_AREA_HYPERPLANE_ORIENTED:
      case SensorDef::SENSOR_AREA_HYPERPLANE: {
        Point<double,3> gul(config.sensors[i].m_gebiet.untenLinks.data, ndim);
        Point<double,3> hul(config.sensors[i].m_gebiet2.untenLinks.data, ndim);
        renderPlane(gul, hul);
      }
        break;
      case SensorDef::SENSOR_AREA_TWO_HYPERPLANES: {
        Point<double,3> gul(config.sensors[i].m_gebiet.untenLinks.data, ndim);
        Point<double,3> gor(config.sensors[i].m_gebiet.obenRechts.data, ndim);
        Point<double,3> hul(config.sensors[i].m_gebiet2.untenLinks.data, ndim);
        renderPlane(gul, hul);
        renderPlane(gor, hul);
      }
        break;
      case SensorDef::SENSOR_AREA_SPHERE: {
        Point<double,3> gul(config.sensors[i].m_gebiet.untenLinks.data, ndim);
        glPushMatrix();
        glTranslated(gul[0], gul[1], gul[2]);
        glCallList(dlIdSolidUnitSphere);
        glPopMatrix();
      }
        break;
      }
    }
    // for(size_t i = 0; i < config.numSensors; ++i) {
    //   config.sensors[i].m_gebiet.render();
    // }
  }

  double vectorLength(DoubleVector const &v) const {
    double const len = norm(v);
    double res = 0;
    if (len > 1) 
      res = log10(len*10) / len;
    else if (len > 0)
      res =  -1 / log10(len/10) / len;
    // else len == 0
    return res;
  }

  void renderVector(DoubleVector const &v, double scaleFactor) const {
    glBegin(GL_LINES);
    {
      DoubleVector vnormed = v;
      double const vlen = vectorLength(v);
      vnormed *= vlen * config.kernelReachH * scaleFactor;

      glColor3dv(&colorMap[1][0]);
      glVertex3f(0, 0, 0);
#if SPH_DIM == 1
      glVertex3f(vnormed[0], 0, 0);
#elif SPH_DIM == 2
      glVertex3d(vnormed[0], vnormed[1], 0);
#else
      glVertex3d(vnormed[0], vnormed[1], vnormed[2]);
#endif
    }
    glEnd();
  }

  /* old code
      if (colorRaster) {
#ifdef SPH_AD
        float r = 0.5, g = 0.5, b = 0.5;
        if (position()[0].diff(0) != 0 or position()[1].diff(0) != 0 or position()[1].diff(0) != 0)
          r = 1;
        if (position()[0].diff(1) != 0 or position()[1].diff(1) != 0 or position()[1].diff(1) != 0)
          g = 1;
        if (position()[0].diff(2) != 0 or position()[1].diff(2) != 0 or position()[1].diff(2) != 0)
          b = 1;
        glColor3f(r, g, b);

#else

	double const rasterc1 = (m_index[0] & SPH_config_hashmask)/double(1 << SPH_config_m);
#if SPH_DIM > 1
	double const rasterc2 = (m_index[1] & SPH_config_hashmask)/double(1 << SPH_config_m);
#else
	double const rasterc2 = 1;
#endif
#if SPH_DIM > 2
	double const rasterc3 = (m_index[2] & SPH_config_hashmask)/double(1 << SPH_config_m);
#else
	double const rasterc3 = 1;
#endif
	//       double const rasterc3 = index[2]/double(1 << (SPH_config_m));
	glColor3f(rasterc1, rasterc2, rasterc3);
#endif
  */

  void clearColorScale() {
    renderInfo.colorRange = 0;
    renderInfo.minScale = 0;
    renderInfo.maxScale = 0;
  }
  void setColorFunction(ParticleColorFunction f) {
    particleColorFunction = f;
    clearColorScale();
  }

  double particleColorHash(particle_type const &p) const {
    return p.hashVal();
  }
  double particleColorIndex(particle_type const &p) const {
    return norm(p.index());
  }
  double particleColorColor(particle_type const &p) const {
    return p.color();
  }
  double particleColorId(particle_type const &
#ifdef SPH_PARTICLE_HAS_ID
                         p
#endif
                         ) const {
    // \todo ID auch aus menüs entfernen...
#ifdef SPH_PARTICLE_HAS_ID
    return p.id();
#else
    return 0;
#endif
  }
  double particleColorNumNeighbours(particle_type const &
#ifdef SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM
                         p
#endif
                         ) const {
    // \todo ID auch aus menüs entfernen...
#ifdef SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM
    return p.numNeighbours();
#else
    return 0;
#endif
  }
  double particleColorClass(particle_type const &p) const {
    return p._class();
  }
  double particleColorFlags(particle_type const &p) const {
    return p.flags();
  }

  double particleColorVelocity(particle_type const &p) const {
    return yafad_value(norm(p.velocity()));
  }
  double particleColorDensity(particle_type const &p) const {
    return yafad_value(p.dichte());
  }
  double particleColorPressure(particle_type const &p) const {
    return yafad_value(p.druck());
  }
  double particleColorHeat(particle_type const &p) const {
    return yafad_value(p.waerme());
  }


  DoubleVector particleVectorVelocity(particle_type const &p) const {
    return doubleVector(p.velocity());
  }
  DoubleVector particleVectorDensity(particle_type const &p) const {
    DoubleVector res;
    res[1] = yafad_value(p.dichte());
    return res;
  }
  DoubleVector particleVectorPressure(particle_type const &p) const {
    DoubleVector res;
    res[1] = yafad_value(p.druck());
    return res;
  }
  DoubleVector particleVectorHeat(particle_type const &p) const {
    DoubleVector res;
    res[1] = yafad_value(p.waerme());
    return res;
  }

#ifdef SPH_AD
  template<int D>
  DoubleVector particleVectorPositionDiff(particle_type const &p) const {
    DoubleVector res;
    for (unsigned i = 0; i < ndim; ++i) {
#if SPH_AD == 1
      res[i] = p.position()[i].diff(D);
#else
      res[i] = p.position()[i].imag()/cvMethEpsilon;
#endif
    }
    return res;
  }
  template<int D>
  DoubleVector particleVectorVelocityDiff(particle_type const &p) const {
    DoubleVector res;
    for (unsigned i = 0; i < ndim; ++i) {
#if SPH_AD == 1
      res[i] = p.velocity()[i].diff(D);
#else
      res[i] = p.velocity()[i].imag()/cvMethEpsilon;
#endif
    }
    return res;
  }
  template<int D>
  DoubleVector particleVectorDensityDiff(particle_type const &p) const {
    DoubleVector res;
#if SPH_AD == 1
    res[1] = p.dichte().diff(D);
#else
    res[1] = p.dichte().imag()/cvMethEpsilon;
#endif
    return res;
  }
  template<int D>
  DoubleVector particleVectorPressureDiff(particle_type const &p) const {
    DoubleVector res;
#if SPH_AD == 1
    res[1] = p.druck().diff(D);
#else
    res[1] = p.druck().imag()/cvMethEpsilon;
#endif
    return res;
  }
  template<int D>
  DoubleVector particleVectorHeatDiff(particle_type const &p) const {
    DoubleVector res;
#if SPH_AD == 1
    res[1] = p.waerme().diff(D);
#else
    res[1] = p.waerme().imag()/cvMethEpsilon;
#endif
    return res;
  }

  template<int D>
  double particleColorPositionDiff(particle_type const &p) const {
    return norm(particleVectorPositionDiff<D>(p));
  }
  template<int D>
  double particleColorVelocityDiff(particle_type const &p) const {
    return norm(particleVectorVelocityDiff<D>(p));
  }
  template<int D>
  double particleColorDensityDiff(particle_type const &p) const {
#if SPH_AD == 1
    return p.dichte().diff(D);
#else
    return p.dichte().imag()/cvMethEpsilon;
#endif
  }
  template<int D>
  double particleColorPressureDiff(particle_type const &p) const {
#if SPH_AD == 1
    return p.druck().diff(D);
#else
    return p.druck().imag()/cvMethEpsilon;
#endif
  }
  template<int D>
  double particleColorHeatDiff(particle_type const &p) const {
#if SPH_AD == 1
    return p.waerme().diff(D);
#else
    return p.waerme().imag()/cvMethEpsilon;
#endif
  }
#endif


  double particleRadiusVolume(particle_type const &p) const {
    double radius = 0;
    if (p.isBoundary()) {
      if (partikelManager.particleTypeBoundaryDefault
          == particle_type::PARTICLE_TYPE_LENNARD_JONES) {
        radius = yafad_value(config.lj_r0);
      } else {
        radius = config.kernelReachH/4;
      }
    } else {
      double const vol = yafad_value(p.masse()) / yafad_value(p.dichte());
#if SPH_NDIM == 1
      radius = vol;
#elsif SPH_NDIM == 2
      radius = sqrt(vol);
#elsif SPH_NDIM == 3
      radius = cbrt(vol);
#else
      radius = pow(vol, 1.0 / ndim);
#endif
    }
    return radius;
  }
  double particleRadiusH(particle_type const &p) const {
    double radius = 0;
    if (p.isBoundary()) {
      radius = config.kernelReachH/4;
    } else {
      radius = config.kernelReachH/2;
    }
    return radius;
  }

  void particleReprPoint(particle_type const &) const {
    glBegin(GL_POINTS);
    {
      glVertex3f(0, 0, 0);
    }
    glEnd();
  }
  void particleReprWireSphere(particle_type const &) const {
    glCallList(dlIdWireUnitSphere);
  }
  void particleReprSolidSphere(particle_type const &) const {
    glCallList(dlIdSolidUnitSphere);
  }

  double normalizeColor(double const v) const {
    double const normalized = (v - renderInfo.minScale) / renderInfo.colorRange;
    // assert(normalized >= 0);
    // assert(normalized <= 1);
    return max(0.0, min(1.0, normalized));
  }
  Color colorMapGradient(double const v) const {
    double const normalized = normalizeColor(v);
//     cerr << "color: " << v << ", " << renderInfo.minScale << ", " << renderInfo.maxScale << " -> " << normalized << "\n";
    assert(normalized >= 0);
    assert(normalized <= 1);
    Color res;
    res[0] = normalized;
    res[1] = 0;
    res[2] = 1 - normalized;
    return res;
  }
  Color colorMapGradientRG(double const v) const {
    double const normalized = normalizeColor(v);
    Color res;
    res[0] = normalized;
    res[1] = 1 - normalized;
//     res[2] = 0.5;
    return res;
  }
  Color colorMapGradientGB(double const v) const {
    double const normalized = normalizeColor(v);
    Color res;
//     res[0] = 0.5;
    res[1] = normalized;
    res[2] = 1 - normalized;
    return res;
  }
  Color colorMapMap(double const v) const {
    double const normalized = normalizeColor(v);
    Color res = colorMap[unsigned(round(normalized * colorMap.size())) % colorMap.size()];
    return res;
  }

  void renderParticle(RenderOptions const &options, particle_type const &p) {
    glPushMatrix();
    {
      Vector const &pos = p.position();
#if SPH_DIM == 1
      glTranslated(yafad_value(pos[0]), 0, 0);
#elif SPH_DIM == 2
      glTranslated(yafad_value(pos[0]), yafad_value(pos[1]), 0);
#else
      glTranslated(yafad_value(pos[0]), yafad_value(pos[1]), yafad_value(pos[2]));
#endif
      
      if (options.showParticleVectors) {
        DoubleVector const vec = (this->*particleVectorFunction)(p);
        renderVector(vec, options.vectorLength);
      }

      Color const pcolor = (this->*colorMapFunction)((this->*particleColorFunction)(p));
      glColor3dv(&pcolor[0]);

      double const radius = (this->*particleSizeFunction)(p) * options.particleSize;
      glScaled(radius, radius, radius);
      (this->*particleReprFunction)(p);
    }
    glPopMatrix();
  }

  /// copy data from simulation to renderInfo struct
  void copyData(RenderOptions const &options) {
    renderInfo.numRendered = 0;
    renderInfo.numThreads = numThreads;
    renderInfo.t = yafad_value(simulationTime);
    renderInfo.timePerStep = timeLastStep;
    renderInfo.step = step;
    renderInfo.npartikel = partikelManager.numParticles();
    renderInfo.numParticlesOut = partikelManager.numParticlesOut();
    renderInfo.frame = m_frameRendered;
    renderInfo.eta = eta();

    resizeIfNeeded(renderInfo.partikelRenderPuffer, partikelManager.numParticles());

    bool first = true;
    double minColor = 0, maxColor = 0;
      
    for(size_t i = 0; i < renderInfo.npartikel; ++i) {
      particle_type const &p = partikelManager.getPartikel(i);
      renderInfo.partikelRenderPuffer[i] = p;
      double const cv = (this->*particleColorFunction)(p);
      if ((not p.isBoundary() and not p.isOut()) or options.colorRangeAll) {
        if (first) {
          minColor = cv;
          maxColor = cv;
          first = false;
        } else {
          minColor = min(minColor, cv);
          maxColor = max(maxColor, cv);
        }
      }
    }
      
    renderInfo.minColor = minColor;
    renderInfo.maxColor = maxColor;
    if (renderInfo.colorRange == 0) {
      renderInfo.minScale = minColor;
      renderInfo.maxScale = maxColor;
    } else {
      if (renderInfo.minScale > renderInfo.minColor) {
        renderInfo.minScale = minColor;
      } else {
        renderInfo.minScale = renderInfo.minScale*(0.9) + minColor*0.1;
      }
      if (renderInfo.maxScale < renderInfo.maxColor) {
        renderInfo.maxScale = maxColor;
      } else {
        renderInfo.maxScale = renderInfo.maxScale*(0.9) + maxColor*0.1;
      }
    }
      
    // cerr << "color: min = " << renderInfo.minColor
    //      << " max = " << renderInfo.maxColor << "\n";
    // cerr << "scale: min = " << renderInfo.minScale
    //      << " max = " << renderInfo.maxScale << "\n";
    if (renderInfo.minScale > renderInfo.minColor) {
      sphout << "<error>The lower bound of the color scale is larger than the "
        "smallest color value: "
             << renderInfo.minScale << " > " << renderInfo.minColor << "</error>\n";
      renderInfo.minScale = renderInfo.minColor;
    }

//     assert(renderInfo.minScale <= renderInfo.minColor);
    assert(renderInfo.maxScale >= renderInfo.maxColor);

    renderInfo.colorRange = renderInfo.maxScale - renderInfo.minScale;
    if (renderInfo.colorRange == 0) 
      renderInfo.colorRange = 1;

#ifdef SPH_RENDER_FORCES
    renderInfo.partikelForces = partikelManager.partikelForces;
    renderInfo.partikelDistances = partikelManager.partikelDistances;
#endif
  }

  void render(RenderOptions const &options) {
    if (m_frameRendered == 0 // first time
        or renderInfo.step != step) {  // simulation step higher than last rendered step
      
      partikelBufferLock.lock();
      copyData(options);
      partikelBufferLock.unlock();
    } else {
      ++renderInfo.numRendered;
    }
    
    ++m_frameRendered;

    if (not options.hideSimulationArea) {
      // weisse linien
      glColor3f(1, 1, 1);
      config.gebiet.render();
    }
    
    if (not options.showSensors) {
      renderSensors();
    }

    if (options.showRaster) {
      renderRaster();
    }

    for (size_t i = 0; i < renderInfo.npartikel; ++i) {
      Partikel<SphADataType> const &p = renderInfo.partikelRenderPuffer[i];
      if (not p.isBoundary() or not options.hideBoundaryParticles) {
        if (not p.isOut() or options.showOutParticles) {
          renderParticle(options, p);
        }
      }
    }

#ifdef SPH_PERIODIC_BOUNDARIES
    if (config.periodicBoundaries.max() > 0) {
      for (unsigned d = 0; d < ndim; ++d) {
        if (config.periodicBoundaries[d]) {
          for (size_t i = 0; i < renderInfo.npartikel; ++i) {
            Partikel<SphADataType> p = renderInfo.partikelRenderPuffer[i];
            if (not p.isBoundary() or not options.hideBoundaryParticles) {
              if (not p.isOut() or options.showOutParticles) {
                p.position()[d] += config.gebiet.obenRechts[d];
                renderParticle(options, p);
              }
            }
          } // for i
        } // if dim d is periodic
      } // for d
    }
#endif

#ifdef SPH_RENDER_FORCES
    glColor3f(0, 1, 0);
    for (size_t i = 0; i < renderInfo.partikelDistances.size(); ++i) {
      renderInfo.partikelDistances[i].render();
    }
    glColor3f(0, 0, 1);
    for (size_t i = 0; i < renderInfo.partikelForces.size(); ++i) {
      renderInfo.partikelForces[i].render();
    }
#endif

  }

#endif

  void save() {
    timeSave.start();
    std::string ofname;

    if (config.resultFormat.find("hdf5") != string::npos) {
//       sphout << "<save frame='" << savedFrame << "' step='" << simstep
//            << "' format='hdf5' file='" << h5PartWriter->filename << "'/>\n";
      ofname = h5PartWriter->filename; // just for XML output
      h5PartWriter->writeH5Part(partikelManager.partikelArray, 
                                savedFrame, 
                                h5PartWriter->SAVE_MOVING 
                                | (config.saveOut ? h5PartWriter->SAVE_OUT : 0)
                                | (savedFrame == 0 ? h5PartWriter->SAVE_BOUNDARY : 0));
    }
    
    if (config.resultFormat.find("vtk") != string::npos) {
      if (savedFrame == 0) {
        // in first save action, save file with boundary particles
        std::ostringstream nstr;
        nstr << config.resultFile << "boundary.vtk";
        VTKWriter<SPHPartikelArray> vtkWriter(nstr.str(),
                                              config.saveFields,
                                              not config.saveASCII);
        vtkWriter.writeVTK(partikelManager.partikelArray, 0, vtkWriter.SAVE_BOUNDARY);
      }

      ostringstream filenameStream;
      filenameStream << config.resultFile << setw(6) << setfill('0') << savedFrame
                     << ".vtk";
      ofname = filenameStream.str();
//       sphout << "<save frame='" << savedFrame << "' step='" << simstep
//              << "' format='vtk' file='" << filenameStream.str() << "'/>\n";
      VTKWriter<SPHPartikelArray> vtkWriter(ofname,
                                            config.saveFields,
                                            not config.saveASCII);

      vtkWriter.writeVTK(partikelManager.partikelArray, 0,
                         vtkWriter.SAVE_MOVING | (config.saveOut ? vtkWriter.SAVE_OUT : 0));
    }

    timeSave.stop();
    sphout << "<save"
      " frame='" << savedFrame << "'"
      " t='" << simulationTime << "'"
      " step='" << simstep
           << "' format='" << config.resultFormat << "' file='" << ofname << "'"
           << " time='" << timeSave.lastDiff() << "'/>\n";
    ++savedFrame;

  }

};
