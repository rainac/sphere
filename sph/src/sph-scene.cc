// this file is included by others
#ifndef jw_sph_scene_cc
#define jw_sph_scene_cc
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include <iostream>
#include <vector>
#include <valarray>
#include <fstream>
#include <errno.h>
#include <sys/stat.h>
#include <H5Part.h>
#  include <omp.h>

#undef  SPH_RENDER
#define SPH_PARTICLE_ID

// for the scene compiler, let particles always have IDs
#define SPH_PARTICLE_HAS_ID

#ifndef SPH_AD
#define yafad_value(x) (x)
#endif

extern "C" char **environ;
#include "common.hh"
#include "scene-config.hh"
#include "partikel.hh"
#include "index.hh"
#include "druck.hh"
#include "h5part-writer.hh"
#include "vtk-writer.hh"
#include "vtk-reader.hh"
#include "utility/hdf5.hh"
#include "omp-util.hh"
#include "sph.hh"

using namespace std;

#ifdef DEBUG_SPH
#define DEB_SPH(x) debugOutput << "sph: " << x << "\n";
#else
#define DEB_SPH(x)
#endif

#if defined DEBUG_SPH
TSOStream<8> debugOutput(std::cerr);
#endif

static std::string const indentUnit = "\t";

typedef vector<valarray<double> > PositionList;
typedef PositionList::iterator PLit;
typedef PositionList::const_iterator PLCit;

// static double deg2rad(double const x) { return x / 360 * 2 * M_PI; }
static double rad2deg(double const x) { return x / (2 * M_PI) * 360; }

static void read3DCoordinateList(string const &fname, PositionList &liste) {
  ifstream ein(fname.c_str(), ios::binary);
  if (not ein) {
    cerr << "error: failed to open 3D coordinate file `"
	 << fname << "': " << strerror(errno) << "\n";
    exit(5);
  } else {
    while(1) {
      valarray<double> v(3);
      ein >> v;
      if (not ein) break;
      liste.push_back(v);
    }
  }
}

/// struct SceneCompiler gathers the components necessary to construct
/// an initial scene for use with the simulation tool.  This class is
/// used in the executable sph-scene. Its role is the same as the SPH
/// class in the simulation program. The constructor takes a slightly
/// modified config object SPHSceneConfig.

/// The function is to read a VTK file or text files with 3D
/// coordinates for particles an merge these particles with water
/// particles generated out of configuration params. 

/// The result is written as VTK. The ability to read VTK allows to
/// use the sph-scene scene compiler incrementally.

struct SceneCompiler {
  
  typedef SPHSceneConfig<SphADataType> SceneConfig;
  SceneConfig const &config;

  GammaGleichung<SphADataType> zustandsGleichung;

  typedef Quader<DoubleVector> Gebiet;
  typedef Partikel<SPH_ACTIVE_DATA_TYPE> ParticleType;

  std::vector<Gebiet> gebietInitVec;
  SPHPartikelArray partikel;

  DoubleVector particleMinPos, particleMaxPos;
  double maxSpatialResolution;

  ofstream sphlogFile;
  Teestream sphout;

  SceneCompiler(SceneConfig const &config) : 
    config(config), 
    zustandsGleichung(config.eqsConfig),
    maxSpatialResolution(),
    sphlogFile(config.logFile.c_str(), ios::binary),
    sphout(cout, sphlogFile)
  {
  }

  ~SceneCompiler() {}

#include "fill-patterns.ncd.enum.hh"
#include "fill-patterns.ncd.cc"

  void initGetInitAreas() {
    size_t initAr = 0;
    for ( ; initAr < config.numInitAreas; ++initAr) {
      std::ostringstream str;
      str << "SPH_PARAM_INIT_A" << initAr << "_P1";
      std::string envnamep1 = str.str();
      std::ostringstream str2;
      str2 << "SPH_PARAM_INIT_A" << initAr << "_P2";
      std::string envnamep2 = str2.str();
      //       cerr << "getenv " << envnamep1 << " and " << envnamep2 << "\n";
      if (getenv(envnamep1.c_str()) == 0) {
        cerr << "error: expected variable " << envnamep1 << ", but is not set\n";
        exit(-5);
      }
      if (getenv(envnamep2.c_str()) == 0) {
        cerr << "error: expected variable " << envnamep2 << ", but is not set\n";
        exit(-5);
      }
      Gebiet res(GetEnvV<ndim, DoubleVector >(envnamep1, 0, -1e100, 1e100),
		 GetEnvV<ndim, DoubleVector >(envnamep2, 1, -1e100, 1e100));
      
//       size_t const ninitArDim = getAreaNDim(res);
//       if (sceneNdim == 0) {
//         sceneNdim = ninitArDim;
//       } else {
//         if (sceneNdim != ninitArDim) {
//           cerr << "error: first area to fill in scene was " << sceneNdim
//                << "-dimensional, but this (" << initAr << ") is " << ninitArDim 
//                << "-dimensional\n";
//           exit(-7);
//         }
//       }

      gebietInitVec.push_back(res);
    }
  }

  template<class T>
  void getNthAreaEnvValue(size_t initAr, std::string const &name, T &val) {
    std::ostringstream str;
    str << "SPH_PARAM_INIT_A" << initAr << "_" << name;
    std::string const envnamev = str.str();
    val = GetEnv(envnamev, val, 0, 1e100);
  }

  void getNthAreaEnvString(size_t initAr, std::string const &name, std::string &val) {
    std::ostringstream str;
    str << "SPH_PARAM_INIT_A" << initAr << "_" << name;
    std::string const envnamev = str.str();
    char const *envString = getenv(envnamev.c_str());
    if (envString) {
      val = envString;
    }
  }

  unsigned compile() {

    // read lists of 3D coordinates from text file(s) (if filename(s) given)
    readBoundaryParticles();

    // read particles from a vtk file (if filename given)
    readInitialVTK();

    // generate particles from specification of areas to be filled
    initGetInitAreas();

    assert(gebietInitVec.size() == config.numInitAreas);

    size_t numGenParticles = 0;
    for (size_t initAr = 0; initAr < gebietInitVec.size(); ++initAr) {
      std::valarray<ParticleType> neuePartikel;
      initArea(initAr, neuePartikel);
      addParticles(neuePartikel);
      numGenParticles += neuePartikel.size();
    }
    
    // check particle positions for conflicts
    unsigned const numErrors = checkParticles();

    sphout << indentUnit << "<spatial-resolution"
      " dx='" << config.deltaX << "'"
      " max-dx='" << maxSpatialResolution << "'"
      " kernel-reach='" << computeKernelReachLength(config.deltaX) << "'"
      "/>\n";

    sphout << indentUnit << "<min-position>" << particleMinPos << "</min-position>\n";
    sphout << indentUnit << "<max-position>" << particleMaxPos << "</max-position>\n";

    sphout << indentUnit << "<generated-particles num-areas='" << gebietInitVec.size()
           << "' num='" << numGenParticles << "'/>\n";
    
    return numErrors;
  }

  static void findMinMax(SPHPartikelArray const &partikel, 
                         DoubleVector &globalMinPos, DoubleVector &globalMaxPos) {
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      DoubleVector minPos, maxPos;
#ifdef _OPENMP
#pragma omp for
#endif
      for(long i = 0; i < long(partikel.size()); ++i) {
        ParticleType const &p = partikel[i];
        for (unsigned k = 0; k < ndim; ++k) {
          minPos[k] = std::min(minPos[k], yafad_value(p.position()[k]));
          maxPos[k] = std::max(maxPos[k], yafad_value(p.position()[k]));
        }
      }
#ifdef _OPENMP
#pragma omp critical
#endif
      for (unsigned k = 0; k < ndim; ++k) {
        globalMinPos[k] = std::min(globalMinPos[k], minPos[k]);
        globalMaxPos[k] = std::max(globalMaxPos[k], maxPos[k]);
      }
    }
  }

//   double computeKernelReachLength(double const dx) const {
//     return config.kernelReachCoeff * 2 * sqrt(dx*dx * sceneNdim);
//   }

  double computeKernelReachLength(double const dx) const {
    return config.kernelReachCoeff * 2 * dx;
  }
  /*  
  double computeKernelReachLength_1(double const dx) const {
    switch(config.numDims) {
    case 1:
      return 2* config.kernelReachCoeff * dx;
    case 2:
      return 2* config.kernelReachCoeff * dx / sqrt(M_PI);
    case 3:
      return 2* config.kernelReachCoeff * dx / cbrt(4*M_PI/3);
    }
    return nan("");
  }
  */
  unsigned checkParticles() {

    findMinMax(partikel, particleMinPos, particleMaxPos);

//     Gebiet const gebiet = Gebiet(particleMinPos, particleMaxPos);
    Gebiet const &gebiet = config.gebiet;

    double const kernelReachH = computeKernelReachLength(config.deltaX);

    RasterInfo rasterInfo(gebiet, kernelReachH);
    HashedInfo hashedInfo(rasterInfo, 5);

    SortedPartikelIndex<SPHPartikelArray> 
      sortedPartikelIndex("csort",hashedInfo.m_index, gebiet, kernelReachH);

    size_t numOut = 0;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(long i = 0; i < long(partikel.size()); ++i) {
      ParticleType &p = partikel[i];
      hashedInfo.updatePartikelIndex(p);
      if (p.isOut()) {
#ifdef _OPENMP
#pragma omp critical
#endif
        {
          ++numOut;
          sphout << "<out id='" << p.id() << "' pos='" << p.position() << "'/>\n";
        }
      }
    }

    sortedPartikelIndex.update(partikel);

    double const minParticleDistance = config.deltaX/1.2;
    double const minParticle2BoundaryDistance = config.deltaX / 1.2;
    
    NachbarArray const nachbarArray(ndim);

    size_t numConflicts = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 10000)
#endif
    for (int i = 0; i < int(partikel.size()); ++i) {
      
      unsigned const partp = sortedPartikelIndex.getPartikelIndex(i);
      ParticleType const &p = partikel[partp];
      BinTriplet range;


        for (unsigned ni = 0; ni < nachbarArray.size(); ++ni) {
          Coord const otherIndex = p.index() + nachbarArray[ni];
          long const otherIndexFlat = hashedInfo.m_index(otherIndex);
          sortedPartikelIndex.getBinInfo(otherIndexFlat, range);
          
          for (unsigned oi = range.offset; oi < range.endBoundary; ++oi) {
            unsigned const parto = sortedPartikelIndex.getPartikelIndex(oi);
            ParticleType const &o = partikel[parto];
            if (o.id() <= p.id()) continue;
            Vector dist = p.position() - o.position();
            if (not p.isBoundary() and not o.isBoundary()) {
              double const d = norm(dist);
              if (d < minParticleDistance) {
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                  sphout << "<conflict type='particle-particle'"
                         << " id1='" << p.id() << "'"
                         << " id2='" << o.id() << "'"
                         << " p1='" << p.position() << "'"
                         << " p2='" << o.position() << "'"
                         << " d='" << d << "'"
                         << " d-rel='" << d/config.deltaX << "'"
                         << "/>\n";
                  //               cerr << "error: Particles " << p << " and " << 
                  //                 o << " are too close: " << norm(dist) << "\n";
                  ++numConflicts;
                }
              }
            } else if (not p.isBoundary() or not o.isBoundary()) {
              double const d = norm(dist);
              if (d < minParticle2BoundaryDistance) {
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                  sphout << "<conflict type='particle-boundary'"
                         << " id1='" << p.id() << "'"
                         << " id2='" << o.id() << "'"
                         << " p1='" << p.position() << "'"
                         << " p2='" << o.position() << "'"
                         << " d='" << d << "'"
                         << " d-rel='" << d / config.deltaX << "'"
                         << "/>\n";
                  //               cerr << "error: Particles " << p << " and " << 
                  //                 o << " are too close: " << norm(dist) << "\n";
                  ++numConflicts;
                }
              }
            }
          }
        }  // for o

      // cerr << "checking: particles " << parti << "/" << partikel.size() << "\r" << flush;
    } // for p
    sphout << indentUnit << "<check-parpos out='" << numOut << "' conflicts='" << numConflicts << "'/>\n";
    return numConflicts + numOut;
  }

  double getVolumeCuboid(Gebiet const &gebietInit) const {
    double area = 1;
    for (size_t i = 0; i < ndim; ++i) {
      if (gebietInit.v[i] != 0) {
        area *= std::abs(gebietInit.v[i]);
      }
    }
    return area;
  }

  double get3DPseudoVolumeCuboid(Gebiet const &gebietInit) const {
    double vol = 1;
    for (size_t i = 0; i < ndim; ++i) {
      if (gebietInit.v[i] != 0) {
        vol *= std::abs(gebietInit.v[i]);
      }  else {
        vol *= config.initReplacementLength;
      }
    }
    return vol;
  }

  double get3DPseudoVolumeBall(Gebiet const &gebietInit) const {
    double vol = 0;
    switch(config.numDims) {
    case 1:
      vol = gebietInit.v[0] * config.initReplacementLength * config.initReplacementLength;
      break;
    case 2:
      vol = gebietInit.v[0] * gebietInit.v[0] * M_PI * config.initReplacementLength;
      break;
    case 3:
      vol = 4.0 / 3.0 * gebietInit.v[0] * gebietInit.v[0] * gebietInit.v[0] * M_PI;
      break;
    default:
      std::cerr << "error: getVolumeBall not implemented for " << ndim << " dimensions\n";
      exit(-1);
      break;
    }
    return vol;
  }

  void generateCuboid(std::valarray<ParticleType> &neuePartikel,
                      Coord const &numPoints, Vector const &deltas, Vector const &offsets) {
    Index<long, ndim> distIndex(numPoints);

    sphout << "\t<distribution type='cuboid'>\n"
      "\t\t<number total='"  << numPoints.prod() << "'"
      " x='" << numPoints[0] << "'"
      " y='" << numPoints[1] << "'"
      " z='" << numPoints[2] << "'"
      "/>\n"
      "\t\t<delta"
      " x='" << deltas[0] << "'"
      " y='" << deltas[1] << "'"
      " z='" << deltas[2] << "'"
      "/>\n"
      "\t</distribution>\n";

    assert(numPoints.prod() == long(distIndex.size));
    assert(numPoints.prod() == long(neuePartikel.size()));

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (long i = 0; i < long(neuePartikel.size()); ++i) {
      ParticleType &p = neuePartikel[i];
      
      Coord ipos = distIndex.invert(i);
      //       cerr << "ipos: " << ipos << "\n";
      p.position() = ipos;
      p.position() *= deltas;
      p.position() += offsets;
      
      assert(not p.hasNaN());
    }
  }

  void generateCuboidEvenOdd(std::valarray<ParticleType> &neuePartikel, Coord const &numPoints, 
                             Vector const &deltas, Vector const &offsets, bool const odd=false) {
    Index<long, ndim> distIndex(numPoints);

    sphout << "\t<distribution type='cuboid-" << (odd?"odd":"even") << "'>\n"
      "\t\t<number total='"  << (numPoints.prod() + (odd?-1:1))/2 << "'"
      " x='" << numPoints[0] << "'"
      " y='" << numPoints[1] << "'"
      " z='" << numPoints[2] << "'"
      "/>\n"
      "\t\t<delta"
      " x='" << deltas[0] << "'"
      " y='" << deltas[1] << "'"
      " z='" << deltas[2] << "'"
      "/>\n"
      "\t</distribution>\n";

    assert(numPoints[0] % 2 == 1);
    // if it's more than one slice, then [1] must be odd. too
    assert(numPoints[2] == 1 or numPoints[1] % 2 == 1);
    assert((numPoints.prod() + (odd?-1:1))/2 == long(neuePartikel.size()));

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (long i = 0; i < long(neuePartikel.size()); ++i) {
      ParticleType &p = neuePartikel[i];
      
      assert(2*i + odd < long(distIndex.size));
      Coord ipos = distIndex.invert(2*i + odd);
      //       cerr << "ipos: " << ipos << "\n";
      p.position() = ipos;
      p.position() *= deltas;
      p.position() += offsets;
      
      assert(not p.hasNaN());
    }
  }


  void distributeParticlesCuboid(std::valarray<ParticleType> &neuePartikel,
                                 Gebiet const &gebiet, double const deltaX, bool const full = false) {
    Vector numIntervals = round(gebiet.v / deltaX);
    Coord distBase = doubleVector(numIntervals + SphADataType(full));
    size_t nPartikel = 1;
    for (size_t i = 0; i < ndim; ++i) {
      if (numIntervals[i] == 0) {
        numIntervals[i] = 1;
        distBase[i] = 1;
      }
      nPartikel *= distBase[i];
    }
    neuePartikel.resize(nPartikel);

    Vector const deltas = Vector(gebiet.v) / numIntervals;
    Vector const offsets = deltas * SphADataType(0.5 * (not full));

    generateCuboid(neuePartikel, distBase, deltas, offsets);
  }

  void distributeParticlesStaggeredCuboid(std::valarray<ParticleType> &neuePartikel,
                                          Gebiet const &gebiet, double const deltaX, 
                                          bool const full = false, bool const odd = false) {
    Vector numIntervals = 2.0 * round(gebiet.v / (M_SQRT2*deltaX)) - double(not full);
    Coord distBase = doubleVector(numIntervals + SphADataType(full));
    size_t nPartikel = 1;
    for (size_t i = 0; i < ndim; ++i) {
      if (numIntervals[i] <= 0) {
        numIntervals[i] = 1;
        distBase[i] = 1;
      }
      nPartikel *= distBase[i];
    }
//     assert(nPartikel % 2 == 1);
    neuePartikel.resize(nPartikel/2 + not odd);

    Vector const deltas = Vector(gebiet.v) / numIntervals;
    Vector const offsets = deltas * SphADataType(0.5 * (not full));

    assert(distBase[0] % 2 == 1);
    assert(distBase[1] % 2 == 1);
    generateCuboidEvenOdd(neuePartikel, distBase, deltas, offsets, odd);
  }


  void distributeParticlesHexagonal(std::valarray<ParticleType> &neuePartikel,
                                    Gebiet const &gebiet, double const deltaX, 
                                    bool const full = false, bool const odd = false) {

    size_t const ninitArDim = config.getAreaNDim(gebiet);

    Vector numIntervals;
    numIntervals[0] = 2*round(gebiet.v[0] / (2*deltaX*0.5)) - double(not full);

    if (ninitArDim > 1) {
      numIntervals[1] = 2*round(gebiet.v[1] / (2*sqrt(3)*0.5*deltaX)) - double(not full);
    }

    if (ninitArDim > 2) {
      numIntervals[2] = round(gebiet.v[2] / deltaX);
    }

    Coord distBase = doubleVector(numIntervals + SphADataType(full));
    size_t nPartikel = 1;
    for (size_t i = 0; i < ndim; ++i) {
      if (numIntervals[i] == 0) {
        numIntervals[i] = 1;
        distBase[i] = 1;
      }
      nPartikel *= distBase[i];
    }

    Vector const deltas = Vector(gebiet.v) / numIntervals;
    Vector const offsets = deltas * SphADataType(0.5 * (not full));

    distBase[2] = 1; // set to one since we will generate numSlices 2D
                     // patterns...
    size_t const numSlices = (deltas[2]==0 ? 1 : yafad_value(numIntervals[2]+full));
    size_t const nParticlesPerSlice = distBase[0]*distBase[1]/2 + not odd;
    neuePartikel.resize(numSlices * nParticlesPerSlice);
    std::valarray<ParticleType> zPartikel(nParticlesPerSlice);
    for(unsigned nz = 0; nz < numSlices; ++nz) {
      generateCuboidEvenOdd(zPartikel, distBase, deltas, offsets, odd);

      for(unsigned i = 0; i < nParticlesPerSlice; ++i) {
        zPartikel[i].position()[2] = nz * deltas[2];
        neuePartikel[nz*nParticlesPerSlice + i] = zPartikel[i];
      }
    }

  }

  void distributeParticlesHexagonalClosePacked(std::valarray<ParticleType> &neuePartikel,
                                               Gebiet const &gebiet, double const deltaX, 
                                               bool const full = false, bool const odd = false) {

    size_t const ninitArDim = config.getAreaNDim(gebiet);

    Vector numIntervals;
    numIntervals[0] = 2*round(gebiet.v[0] / (2*deltaX*0.5)) - double(not full);

    if (ninitArDim > 1) {
      numIntervals[1] = 2*round(gebiet.v[1] / (2*sqrt(3)*0.5*deltaX)) - double(not full);
    }

    if (ninitArDim > 2) {
      numIntervals[2] = 2*round(gebiet.v[2] / (2*sqrt(2.0/3)*deltaX)) - double(not full);
    }

    Coord distBase = doubleVector(numIntervals + SphADataType(full));
    for (size_t i = 0; i < ndim; ++i) {
      if (numIntervals[i] == 0) {
        numIntervals[i] = 1;
        distBase[i] = 1;
      }
    }

    Vector const deltas = Vector(gebiet.v) / numIntervals;
    Vector const offsets = deltas * SphADataType(0.5 * (not full));

    Coord distBaseInter(distBase);
    distBaseInter[1] -= 1;

    // if odd == 0: numSlices layers even - odd - ... - even
    // if odd == 1: numSlices layers odd - even - ... - odd
    // intermediate layers have one row less in BOTH cases

    std::valarray<ParticleType> zPartikelEven(distBase[0]*distBase[1]/2 + 1 - (odd * (distBase[0]/2 + 1)));
    std::valarray<ParticleType> zPartikelOdd(distBase[0]*(distBase[1])/2 - ((not odd) * (distBase[0]/2 + 1)));
    size_t const nPartikel = (distBase[2]/2 + not odd) * zPartikelEven.size() + 
      (distBase[2]/2 + odd) * zPartikelOdd.size();
    distBaseInter[2] = distBase[2] = 1;
    size_t const numSlices = (deltas[2]==0 ? 1 : yafad_value(numIntervals[2]+full));
    neuePartikel.resize(nPartikel);
    size_t particleOffset = 0;
    for(unsigned nz = 0; nz < numSlices; ++nz) {
      if ((nz % 2 == 0 and not odd)
          or (nz % 2 == 1 and odd)) {
        generateCuboidEvenOdd(zPartikelEven, (odd ? distBaseInter : distBase), deltas, offsets);
        for(unsigned i = 0; i < zPartikelEven.size(); ++i) {
          if (nz % 2 == 1) {
            zPartikelEven[i].position()[1] += sqrt(3)/3 * deltas[1];
          }
          zPartikelEven[i].position()[2] = nz * deltas[2];
          neuePartikel[particleOffset + i] = zPartikelEven[i];
        }
        particleOffset += zPartikelEven.size();
      } else {
        generateCuboidEvenOdd(zPartikelOdd, (odd ? distBase : distBaseInter), deltas, offsets, 1);
        for(unsigned i = 0; i < zPartikelOdd.size(); ++i) {
          if (nz % 2 == 1) {
            zPartikelOdd[i].position()[1] += sqrt(3)/6 * deltas[1];
          }
          zPartikelOdd[i].position()[2] = nz * deltas[2];
          neuePartikel[particleOffset + i] = zPartikelOdd[i];
        }
        particleOffset += zPartikelOdd.size();
      }
      
    }
    assert(particleOffset == neuePartikel.size());
  }



  void distributeParticlesBall(std::valarray<ParticleType> &neuePartikel,
                               Gebiet const &_gebietInit, 
                               double const dx) {
    Gebiet gebietInit(_gebietInit.untenLinks - _gebietInit.v,
                      _gebietInit.obenRechts);
    sphout << "\t<distribution type='sphere'>\n";
    gebietInit.writeXML(sphout, "\t\t");
    distributeParticlesCuboid(neuePartikel, gebietInit, dx);
    Vector const center;
    size_t numCopied = 0;
    std::valarray<ParticleType> tmpPartikel(neuePartikel.size());
    for(size_t i = 0; i < neuePartikel.size(); ++i) {
      ParticleType &p = neuePartikel[i];
      p.position() -= _gebietInit.v;
      Vector distC = p.position() - center;
      if (norm(distC) <= _gebietInit.v[0]) {
        tmpPartikel[numCopied] = p;
        ++numCopied;
      } else {
//         cerr << "p: " << p.position() << " c: " << center << " n: " << norm(distC) << "\n";
      }
    }
    neuePartikel.resize(numCopied);
    neuePartikel = tmpPartikel[slice(0, numCopied, 1)];
    assert(neuePartikel.size() == numCopied);
    sphout << "\t<number total='" << numCopied << "'/>\n";
    sphout << "\t</distribution>\n";
  }

  void distributeParticlesBall2_2D(std::valarray<ParticleType> &neuePartikel,
                                 Gebiet const &gebietInit, double const deltaX) {
    double const radius = gebietInit.v[0];
    size_t const numrings = round(radius / deltaX);
    double const rdx = radius / numrings;

    std::vector<ParticleType> tmpPartikel;
    ParticleType p;

    sphout << "\t<distribution type='sphere-2d'>\n"
      "\t\t<number"
      " r='" << numrings << "'"
      "/>\n"
      ;

    sphout << "\t<delta r='" << rdx << "'/>\n";
      
    double phase = 0;
// #ifdef _OPENMP   
// #pragma omp parallel for   not parallelizable because of push_back()
// #endif
    for (long nr = 0; nr < long(numrings); ++nr) {
      double const circ = 2*M_PI*rdx*nr;
      size_t const numpos = round(circ / deltaX);
      if (numpos == 0) {
        // ein partikel in die mitte
        p.position() = 0;
        tmpPartikel.push_back(p);
      } else {
        // konzentrischer ring mit radius nr*rdx
        double const cdx = 2*M_PI / numpos;
        phase += cdx;
        for (size_t np = 0; np < numpos; ++np) {
          p.position()[0] = nr*rdx*sin(np*cdx + phase);
          p.position()[1] = nr*rdx*cos(np*cdx + phase);
          tmpPartikel.push_back(p);
        }
        sphout << "\t<ring"
          " r='" << rdx*nr << "'"
          " num='" << numpos << "'"
          " dphi='" << rad2deg(cdx) << "'"
          " ds='" << rdx*nr*sin(cdx) << "'"
          "/>\n";
      }
    }
    neuePartikel.resize(tmpPartikel.size());
    std::copy(tmpPartikel.begin(), tmpPartikel.end(), &neuePartikel[0]);

    sphout << "\t<total>" << tmpPartikel.size() << "</total>\n";
    sphout << "\t</distribution>\n";
  }

  Vector polarToCart(double const r, double const phi, double const theta) const {
    Vector res;
    res[0] = r * cos(phi)*sin(theta);
    res[1] = r * sin(phi)*sin(theta);
    res[2] = r * cos(theta);
    return res;
  }

  /// \todo distributeParticlesBall2 zum standard machen,
  /// distributeParticlesBall als alternative behalten

  /// \todo weitere kristallgitter implementieren

  void distributeParticlesBall2_3D(std::valarray<ParticleType> &neuePartikel,
                                  Gebiet const &gebietInit, double const deltaX) {
    double const radius = gebietInit.v[0];
    size_t numrings = round(radius / deltaX);
    double const dr = radius / numrings;

    std::vector<ParticleType> tmpPartikel;
    ParticleType p;

    // ein partikel in die mitte
    p.position() = 0;
    tmpPartikel.push_back(p);

    //    double phase = 0;
// #ifdef _OPENMP
// #pragma omp parallel for   not parallelizable because of push_back()
// #endif
    for (long nr = 1; nr < long(numrings); ++nr) {

      double const circ = 2*M_PI*dr*nr;

      size_t const numTheta = round((circ/2) / deltaX);
      double const dTheta = M_PI / numTheta;
      
      sphout << "\t<sphere"
        " r='" << dr*nr << "'"
        " num-lat='" << numTheta << "'"
        " dtheta='" << rad2deg(dTheta) << "'"
        " ds='" << dr*nr*sin(dTheta) << "'"
        "/>\n";

      // ein partikel in den Zenith
      p.position() = polarToCart(dr*nr, 0, 0);
      tmpPartikel.push_back(p);
      // ein partikel in den Nadir
      p.position() = polarToCart(dr*nr, 0, M_PI);
      tmpPartikel.push_back(p);

      for (size_t nt = 1; nt < numTheta; ++nt) {
        double const latRadius = nr*dr*sin(nt*dTheta);
        double const latCirc = 2*M_PI*latRadius;
        
        size_t const numPhi = round(latCirc / deltaX);
        double const dPhi = 2*M_PI/numPhi;

        sphout << "\t<ring"
          " theta='" << rad2deg(nt*dTheta) << "'"
          " r='" << latRadius << "'"
          " num='" << numPhi << "'"
          " dphi='" << rad2deg(dPhi) << "'"
          " ds='" << latRadius*sin(dPhi) << "'"
          "/>\n";

        for (size_t np = 0; np < numPhi; ++np) {
          p.position() = polarToCart(dr*nr, np*dPhi, nt*dTheta);
          tmpPartikel.push_back(p);
        }
      }
    }
    neuePartikel.resize(tmpPartikel.size());
    std::copy(tmpPartikel.begin(), tmpPartikel.end(), &neuePartikel[0]);
  }

  void distributeParticlesBall2(std::valarray<ParticleType> &neuePartikel,
                                Gebiet const &gebietInit, double const deltaX) {
    unsigned const nd = config.getAreaNDim(gebietInit);
    if (nd == 2) {
      distributeParticlesBall2_2D(neuePartikel, gebietInit, deltaX);
    } else if (nd == 3) {
      distributeParticlesBall2_3D(neuePartikel, gebietInit, deltaX);
    }
  }

  void distortAndShiftParticlePositions(std::valarray<ParticleType> &neuePartikel, 
                                        Vector const &which, double const posEpsilon,
                                        Gebiet const &gebietInit) const {
    static unsigned short erand48Seed[20] = {0};
#ifdef _OPENMP
#pragma omp threadprivate(erand48Seed)
#pragma omp parallel
#endif
    {
      Vector rand;
      memcpy(erand48Seed, &config.randomSeed, sizeof(config.randomSeed));

#ifdef _OPENMP
#pragma omp for
#endif
      for (long i = 0; i < long(neuePartikel.size()); ++i) {
        ParticleType &p = neuePartikel[i];
        
        for (size_t di = 0; di < ndim; ++di) {
          if (which[di] != 0) {
            rand[di] = erand48(erand48Seed) - 0.5;
          }
        }
        rand *= posEpsilon;
        p.position() += rand;

        // add offset to position
        p.position() += gebietInit.untenLinks;
      }
    }
  }

  void setParticleDensities(std::valarray<ParticleType> &neuePartikel, double const dens) const {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (long i = 0; i < long(neuePartikel.size()); ++i) {
      ParticleType &p = neuePartikel[i];
      p.dichte() = dens;
    }
  }

  void computeParticleDensities(std::valarray<ParticleType> &neuePartikel, 
                                double const areaDens, double const yOberflaeche) const {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (long i = 0; i < long(neuePartikel.size()); ++i) {
      ParticleType &p = neuePartikel[i];
      
      /// \todo possibility to use an arbitrary gravity vector and bottom plane
      SphADataType const tiefe = yOberflaeche - p.position()[1];
      
      // material has pressure given by g * rho * depth
      SphADataType const hydroDruck = areaDens * config.gravity * tiefe;
      
      // 	  cerr << "Tiefe bei " << p.position() << ": " << tiefe << "\n";
      // 	  cerr << "Druck bei Tiefe " << tiefe << ": " << hydroDruck << "\n";
      // 	  p.dichte() = areaDens 
      // 	    * std::pow((1 + 
      // 			areaDens * config.gravity * (tiefe)
      // 			/ config.pressureB), 
      // 		       1.0/config.gamma);
      p.dichte() 
        = pow(pow(areaDens,double(config.eqsConfig.gamma)) * 
              (hydroDruck / config.eqsConfig.pressureB + 1), 
              1.0/config.eqsConfig.gamma);
//       cerr << "p.pos: " << p.position() << "\n";
//       cerr << "Dichte bei Tiefe " << tiefe << ": " << p.dichte() << "\n";
    }
  }

  
  double sumParticleMass(std::valarray<ParticleType> &neuePartikel) const {
    double sum = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sum)
#endif
    for (long i = 0; i < long(neuePartikel.size()); ++i) {
      ParticleType &p = neuePartikel[i];
      sum += p.masse();
    }
    return sum;
  }

  void initArea(size_t const initAr, std::valarray<ParticleType> &neuePartikel) {
      Gebiet const &gebietInit = gebietInitVec[initAr];

      size_t const ninitArDim = config.getAreaNDim(gebietInit);
      assert(ninitArDim == config.numDims);

      long areaNPartikel = config.numInitParticles;
      getNthAreaEnvValue(initAr, "N", areaNPartikel);

      double deltaX = 0;

      if (areaNPartikel == 0) {
        deltaX = config.deltaX;
        getNthAreaEnvValue(initAr, "DX", deltaX);
      } else {
        double const areaArea = getVolumeCuboid(gebietInit);
        deltaX = pow(areaNPartikel / areaArea, -1.0/ninitArDim);
      }

      double areaDens = config.initDensity;
      getNthAreaEnvValue(initAr, "DENSITY", areaDens);

      double areaHeat = config.initHeat;
      getNthAreaEnvValue(initAr, "HEAT", areaHeat);

      double posEpsilon = config.initEpsilon;
      getNthAreaEnvValue(initAr, "EPSILON", posEpsilon);

      std::string areaNPatterntypeName = "cuboid";
      getNthAreaEnvString(initAr, "PATTERN", areaNPatterntypeName);
      int const areaNPatterntype = getFillPatternsValue(areaNPatterntypeName);

      int areaNFill = 1;
      getNthAreaEnvValue(initAr, "FILL", areaNFill);

      int areaNOdd = 0;
      getNthAreaEnvValue(initAr, "ODD", areaNOdd);

      int areaNVelType = 0;
      getNthAreaEnvValue(initAr, "VEL_TYPE", areaNVelType);

      int areaNColor = initAr + 2;
      getNthAreaEnvValue(initAr, "COLOR", areaNColor);

      std::string areaNClassName = "moving_default";
      getNthAreaEnvString(initAr, "CLASS", areaNClassName);
      int const areaNClass = ParticleType::getParticleClassesValue(areaNClassName);

      bool freiSchwebend = 1;
      getNthAreaEnvValue(initAr, "FREE", freiSchwebend);

      double waterSurfaceY = gebietInit.obenRechts[1];
      getNthAreaEnvValue(initAr, "SURFACE_Y", waterSurfaceY);

      Vector areaNVel(ndim);
      {
	std::ostringstream str;
	str << "SPH_PARAM_INIT_A" << initAr << "_V";
	std::string envnamev = str.str();
	areaNVel = SPHConfigBase::GetEnvVec(envnamev, 0, -1e100, 1e100);
      }
      
      double const areaNVolume = ((areaNPatterntype == FILL_PATTERN_SPHERE 
                                   or areaNPatterntype == FILL_PATTERN_SPHERE_ALT) ? 
                                  get3DPseudoVolumeBall(gebietInit) : get3DPseudoVolumeCuboid(gebietInit));
  
      switch(areaNPatterntype) {
      case FILL_PATTERN_CUBOID:
        distributeParticlesCuboid(neuePartikel, gebietInit, deltaX, areaNFill);
        break;
      case FILL_PATTERN_STAGGERED_CUBOID:
        distributeParticlesStaggeredCuboid(neuePartikel, gebietInit, deltaX, areaNFill, areaNOdd);
        break;
        /// \todo implemenent these two fill patterns
//       case FILL_PATTERN_FACE_CENTERED_CUBOID:
//         distributeParticlesFaceCentered(neuePartikel, gebietInit, deltaX);
//         break;
//       case FILL_PATTERN_BODY_CENTERED_CUBOID:
//         distributeParticlesBodyCentered(neuePartikel, gebietInit, deltaX);
//         break;
      case FILL_PATTERN_HEXAGONAL:
        distributeParticlesHexagonal(neuePartikel, gebietInit, deltaX, areaNFill, areaNOdd);
        break;
      case FILL_PATTERN_HCP:
        distributeParticlesHexagonalClosePacked(neuePartikel, gebietInit, deltaX, areaNFill, areaNOdd);
        break;
      case FILL_PATTERN_SPHERE:
        distributeParticlesBall(neuePartikel, gebietInit, deltaX);
        break;
      case FILL_PATTERN_SPHERE_ALT:
        distributeParticlesBall2(neuePartikel, gebietInit, deltaX);
        break;
      default:
        cerr << "error: fill pattern type: " << areaNPatterntype << " (" << 
          areaNPatterntypeName << ") not implemented\n";
        exit(9);
        break;
      }

      // set velocity
      if (areaNVelType == 0) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (long i = 0; i < long(neuePartikel.size()); ++i) {
          ParticleType &p = neuePartikel[i];
          p.velocity() = areaNVel;
        }
      } else if (areaNVelType == 1) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (long i = 0; i < long(neuePartikel.size()); ++i) {
          ParticleType &p = neuePartikel[i];
          p.velocity() = p.position() * areaNVel; // before shifting positions!
        }
      } else if (areaNVelType == 2) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (long i = 0; i < long(neuePartikel.size()); ++i) {
          ParticleType &p = neuePartikel[i];
          p.velocity() = (0.5 + 0.5 * sin(p.position())) * areaNVel; // before shifting positions!
        }
      } else {
        std::cerr << "error: invalid value for VEL_TYPE of area: " << areaNVelType << "\n";
      }

      distortAndShiftParticlePositions(neuePartikel, gebietInit.v, posEpsilon, gebietInit);

      if (freiSchwebend) {
        setParticleDensities(neuePartikel, areaDens);
      } else {
        computeParticleDensities(neuePartikel, areaDens, waterSurfaceY);
      }

      double const masseGesamt = areaNVolume * areaDens;
      double const partikelMasse = masseGesamt / neuePartikel.size();

      DEB_SPH("init " << ninitArDim << "-dim partikel area: "
	      << initAr << ": volumen " << areaNVolume
	      << " num p: " << areaNPartikel
	      << " vol/p " << areaNVolume / neuePartikel.size() << " mass " << areaDens * areaNVolume
	      << " m/p " << partikelMasse
	      << " color: " << areaNColor
	      << " pattern: " << areaNPatterntype
	      << " class: " << ParticleType::getParticleClassesName(areaNClass)
	      << " vel: " << areaNVel
	      << " eps: " << posEpsilon);

      SphADataType const areaRes
        = pow(partikelMasse/config.eqsConfig.rho0, 1.0/config.numDims);
      maxSpatialResolution = max(maxSpatialResolution, yafad_value(areaRes));

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (long i = 0; i < long(neuePartikel.size()); ++i) {
	ParticleType &p = neuePartikel[i];

        if (config.colorParticles == 0) {
          p.color() = (areaNColor ? areaNColor : initAr + 2);
        } else if (config.colorParticles == 1) {
          p.color() = i;
        }

	// init. mass set to dens. * volPerParticlel
	p.masse() = partikelMasse;
	p.waerme() = areaHeat;
        
        p._class() = areaNClass;
        if (p._class() == ParticleType::PARTICLE_TYPE_LENNARD_JONES
            or p._class() == ParticleType::PARTICLE_TYPE_POINT_MIRRORED
            or p._class() == ParticleType::PARTICLE_TYPE_BOUNDARY_DEFAULT) {
          // is a boundary particle
          // this is needed only for the checkParticles function to work correctly
          p.flags() |= ParticleType::PARTICLE_FLAG_BOUNDARY;
        }

	assert(p.isOut() == false);
	assert(p.waerme() >= 0);
	assert(p.dichte() > 0);
	assert(p.masse() > 0);
        
        assert(not p.hasNaN());

	zustandsGleichung.updatePressure(p);
      }

      sphout << indentUnit << "<area"
	" dim='" << ninitArDim << "'"
	" fill-pattern='" << getFillPatternsName(areaNPatterntype) << "'"
	" fill='" << areaNFill << "'"
	" odd='" << areaNOdd << "'"
	" vel-type='" << areaNVelType << "'"
	" color='" << areaNColor << "'\n" << indentUnit << 
	" volume='" << areaNVolume << "'"
// 	" pseudo-area='" << areaArea << "'\n"
        " dx='" << deltaX << "'"
	" num-particles='" << neuePartikel.size() << "'"
	" density='" << areaDens << "'"
	" mass='" << areaDens * areaNVolume << "'"
	" p-mass='" << partikelMasse << "'\n" << indentUnit << 
	" p-kernel='" << pow(partikelMasse/areaDens, 1.0 / config.numDims) << "'"
	" free-floating='" << (freiSchwebend ? "yes" : "no") << "'"
        " waterSurfaceY='" << waterSurfaceY << "'"
	" epsilon='" << posEpsilon << "'"
	">\n" << indentUnit << indentUnit << 
	"<vel x='" << areaNVel[0]
	   << "' y='" << (ndim > 1 ? areaNVel[1] : 0)
	   << "' z='" << (ndim > 2 ? areaNVel[2] : 0) << "'/>\n"
	;
      gebietInit.writeXML(sphout, indentUnit + indentUnit);
      sphout << indentUnit << indentUnit << "<check summed-mass='" << sumParticleMass(neuePartikel) << "'/>\n";
      sphout << indentUnit << "</area>\n";
  }


  void readBoundaryParticles() {
    // the boundary particles must exist and not be destroyed while running
    std::valarray<ParticleType> boundaryParticles;
    {
      // read boundary particle positions
      PositionList boundaryParticlePos;
      if (not config.boundaryFile.empty())
	read3DCoordinateList(config.boundaryFile, boundaryParticlePos);

      size_t const numb1 = boundaryParticlePos.size();

      sphout << "<boundary-particles"
	" num='" << numb1 << "'/>\n";

      // setup boundary particles and tell SPH about them
      boundaryParticles.resize(numb1);

      {
	PLCit it = boundaryParticlePos.begin();
	for (size_t i = 0; i < numb1; ++i, ++it) {
	  assert(it != boundaryParticlePos.end());
	  ParticleType &b = boundaryParticles[i];
	  b.flags() = ParticleType::PARTICLE_FLAG_BOUNDARY;
	  b._class() = ParticleType::PARTICLE_TYPE_BOUNDARY_DEFAULT;

	  b.position() = *it;
	  b.position() *= config.boundaryScale;
	  b.position() += config.boundaryTranslate;
	}
      }
    }

#ifdef SPH_AD
    if (config.boundaryDerivFilePattern.size() > 0) {
      
      size_t const fnameBuffSz = PATH_MAX;
      char *buffer = new char[fnameBuffSz + 1];
      
      for (unsigned nDir = 0; nDir < SPH_AD_NDIR; ++nDir) {
        int const nwritten = snprintf(buffer, fnameBuffSz, config.boundaryDerivFilePattern.c_str(), nDir);
        string filename = string(buffer, nwritten);
        if (filename == config.boundaryDerivFilePattern) {
          cerr << "error: config.boundaryDerivFilePattern should contain a %d printf conversion spec.\n";
          cerr << "error: config.boundaryDerivFilePattern is `" << config.boundaryDerivFilePattern << "'\n";
          exit(15);
        }
        struct stat st;
        if (stat(filename.c_str(), &st) != 0) {
          if (nDir == 0) {
            cerr << "error: config.boundaryDerivFilePattern was set (" << config.boundaryDerivFilePattern
                 << "), but not file was found: `" << filename << "'\n";
            exit(17);
          } else {
            sphout << "<no-more-boundary-derivatives/>\n";
          }
          break;
        }
        // read derivatives of boundary particle positions
        PositionList boundaryParticleDerivatives;
        read3DCoordinateList(filename, boundaryParticleDerivatives);
        
        if (boundaryParticleDerivatives.size() != boundaryParticles.size()) {
          cerr << "error: list of derivatives should have the same size as list of boundary particles (" 
               << boundaryParticles.size() << ") but has only " << boundaryParticleDerivatives.size() << "\n";
          exit(16);
        }
        
        sphout << "<boundary-derivative filename='" << filename << "' dir='" << nDir << "'/>\n";

        for (size_t i = 0; i < boundaryParticles.size(); ++i) {
          ParticleType &b = boundaryParticles[i];
          for (size_t d = 0; d < SPH_DIM; ++d) {
            b.position()[d].diff(nDir) = boundaryParticleDerivatives[i][d];
          }
        }
      }
      delete[] buffer;
    }
#endif // SPH_AD

    addParticles(boundaryParticles);
  }

  void readInitialVTK() {
    if (config.sceneInFile.empty()) 
      return;
    std::valarray<ParticleType> neuePartikel;
    VTKReader<std::valarray<ParticleType> > 
      reader(config.sceneInFile, config.fieldsToLoad);
    reader.readVTK(neuePartikel);

    size_t const ns = neuePartikel.size();
    for(size_t i = 0; i < ns; ++i) {
      neuePartikel[i].position() *= config.initialVTKScale;
      neuePartikel[i].position() += config.initialVTKTranslate;
    }

    addParticles(neuePartikel);
    sphout << "<import-vtk"
      " num-particles='" << neuePartikel.size() << "'/>\n";
  }

  void addParticles(std::valarray<ParticleType> &neuePartikel) {
    size_t const os = partikel.size();
    size_t const ns = neuePartikel.size();
    nonDestructiveResize(partikel, partikel.size() + ns);
    
    for(size_t i = 0; i < ns; ++i) {
      partikel[i + os] = neuePartikel[i];
      partikel[i + os].id() = i + os;
    }
  }

  void writeSceneConfigFile(std::string const &name) {
    ofstream aus(name.c_str(), std::ios::binary);
    aus.precision(17);
    double const kernelH = computeKernelReachLength(config.deltaX);
    aus << "export SPH_PARAM_G='" << config.gravity << "'\n";
    aus << "export SPH_PARAM_B='" << config.eqsConfig.pressureB << "'\n";
    aus << "export SPH_PARAM_GAMMA='" << config.eqsConfig.gamma << "'\n";
    aus << "export SPH_PARAM_RHO0='" << config.eqsConfig.rho0 << "'\n";
    aus << "export SPH_PARAM_H='" << kernelH << "'\n";
    //    aus << "export SPH_PARAM_LJ_R0='" << kernelH/2 << "'\n";
    aus << "if [[ -z \"$SPH_PARAM_QUADER_P2\" ]]; then\n";
    aus << "   export SPH_PARAM_QUADER_P2='" << config.gebiet.obenRechts << "'\n";
    aus << "fi";
  }

};

int curNumThreads = 1;

int main(int, char *[]) {

  SPHSceneConfig<SphADataType> config;

#ifdef _OPENMP
  omp_set_num_threads(config.numThreads);
  curNumThreads = config.numThreads;
#endif

  cerr.precision(12);
  cout.precision(17);

  int result = 0;

  {
    if (config.sceneOutFile.empty()) {
      cerr << "error: an output file name must be given\n";
      exit(EXIT_FAILURE);
    }
    
    HDF5 hdf5;
    SceneCompiler sceneCompiler(config);

    sceneCompiler.sphout << "<sphscene"
      " compiled='" __DATE__ " " __TIME__ "'"
      ">\n";
    
    Uname().writeXML(sceneCompiler.sphout, indentUnit);
    sceneCompiler.sphout << "\n";
    
    StrFTime("%Y-%m-%d %H:%M:%S %z (%c)").writeXML(sceneCompiler.sphout, indentUnit);
    sceneCompiler.sphout << "\n";
    
    config.writeXML(sceneCompiler.sphout, indentUnit);
    
    unsigned const numErrors = sceneCompiler.compile();
    if (numErrors) {
      result = 5;
    }

    std::string outFileNameBase, outFileName;
    outFileNameBase = config.sceneOutFile;
    if (numErrors) {
      outFileNameBase += ".fail";
    }
    if (config.resultFormat.find("hdf5") != string::npos) {
      outFileName = outFileNameBase + ".h5";
      H5PartWriter<SPHPartikelArray> h5PartWriter(outFileName, config.saveFields);
      sceneCompiler.sphout << indentUnit << "<output>" << h5PartWriter.filename << "</output>\n";
      h5PartWriter.writeH5Part(sceneCompiler.partikel, 0, 0xff);

      sceneCompiler.writeSceneConfigFile(outFileName + ".env");
    }

    if (config.resultFormat.find("vtk") != string::npos) {
      outFileName = outFileNameBase + ".vtk";
      VTKWriter<SPHPartikelArray>    vtkWriter(outFileName, config.saveFields, 
					       not config.saveASCII);
      sceneCompiler.sphout << indentUnit << "<output>" << vtkWriter.writer->GetFileName() << "</output>\n";
      vtkWriter.writeVTK(sceneCompiler.partikel, 0, 0xff);

      sceneCompiler.writeSceneConfigFile(outFileName + ".env");
    }

    sceneCompiler.sphout << "</sphscene>\n";
  }

  return result;

}

#endif
