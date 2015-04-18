#ifndef sph_scene_config_hh
#define sph_scene_config_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include "quader.hh"
#include "utility/getenv.hh"
#include "vector.hh"
#include "sensor-def.hh"
#include "sph-config.hh"


/// Gathers all configuration constants.
/// Can write itself to a stream in XML format.
template<class T>
struct SPHSceneConfig : public SPHConfigBase {

  SPHEQSConfig<T> eqsConfig;

  /// The number of threads. 
  unsigned const numThreads;

  /// The computing domain.
  Gebiet gebiet;

  /// The number of simulation steps.
  unsigned const numDims;

  /// The number of particles in one "init area".
  size_t const numInitParticles;

  /// The number of "init area"s.
  size_t const numInitAreas;

  /// The spatial rastering length (same for all dimensions).
  double const deltaX;

  /// The resolution factor for boundary particles (> 1 makes particles denser)
  /// This values is not actually used in his program
  double const boundaryRes;

  /// The coefficient to multiply H computed from deltaX with.
  double const kernelReachCoeff;

  /// The gravitational constant
  double const gravity;

  /// The dynamicBulkViscosity of the medium
  double const dynamicBulkViscosity;

  /// The initial density of the medium in an "init area", if not
  /// overridden for a particular "init area".
  double const initDensity;

  /// The randomness in the initial particle placement in an "init
  /// area", if not overridden for a particular "init area".
  double const initEpsilon;

  double const initHeat;

  double const initReplacementLength;

  long const colorParticles;

  long const randomSeed;

  long const saveASCII;

  std::string const saveFields;

  /// The name of the file to load boundary particles from.
  std::string const boundaryFile;

  /// The sprintf pattern to generate filenames to load derivatives of boundary particles from
  std::string const boundaryDerivFilePattern;

  /// The name of the file to save the scene to.
  std::string const sceneOutFile;

  /// The name of a VTK file to load particles from (merging it into the result).
  std::string const sceneInFile;

  /// The name of the file to log output to.
  std::string const logFile;

  double const initialVTKScale;
  double const boundaryScale;

  Vector const initialVTKTranslate;
  Vector const boundaryTranslate;

  /// The particle data fields to load from initialFile. If the field's name is present
  /// as a substring, the field is loaded.
  std::string fieldsToLoad;

  /// The format to use for writing the result file (VTK or H5Part/HDF5).
  std::string const resultFormat;

  /// Initializes all the constants to the values given by defaults or
  /// set via the environment variables. These have the same name, but
  /// all uppercased and prefixed with "SPH_PARAM_".
  SPHSceneConfig() :
#ifdef _OPENMP
    numThreads(GetEnv("SPH_PARAM_THREADS", 4, 1, 1e100)),
#else
    numThreads(1),
#endif
    gebiet(DoubleVector(), 
	   GetEnvV<ndim, DoubleVector >("SPH_PARAM_QUADER_P2", 1)),
    numDims(getAreaNDim(gebiet)),
    numInitParticles(GetEnv("SPH_PARAM_N", 0, 1, 1e12)),
    numInitAreas(GetEnv("SPH_PARAM_INIT_NUM", 0, 0, 1e5)),
    deltaX(GetEnv("SPH_PARAM_DX", 6.275e-3, 1e-100, 1e100)),
    boundaryRes(GetEnv("SPH_PARAM_BRES", 0, 0, 1e100)),
    kernelReachCoeff(GetEnv("SPH_PARAM_H_COEFF", 1., 1e-100, 1e100)),
    gravity(GetEnv("SPH_PARAM_G", 9.81, 0, 1e100)),
    dynamicBulkViscosity(GetEnv("SPH_PARAM_MU", 8.9e-4, 0, 1e100)),
    initDensity(GetEnv("SPH_PARAM_INIT_DENSITY", yafad_value(eqsConfig.rho0), 1e-100, 1e100)),
    initEpsilon(GetEnv("SPH_PARAM_INIT_EPSILON", 0, 0, 1e100)),
    initHeat(GetEnv("SPH_PARAM_INIT_HEAT", 0, -1e100, 1e100)),
    initReplacementLength(GetEnv("SPH_PARAM_INIT_REPLACEMENT_LENGTH", 
				 1, 1e-100, 1e100)),
    colorParticles(GetEnv("SPH_PARAM_COLOR_PARTICLES", 0, 0, 100)),
    randomSeed(GetEnv("SPH_PARAM_RANDOM_SEED", time(0), 0, 1e100)),
    saveASCII(GetEnv("SPH_PARAM_SAVE_ASCII", 0, 0, 1)),

    saveFields(GetEnvString("SPH_PARAM_SAVE_FIELDS", "masse dichte waerme druck coord vel class color"
#ifdef SPH_AD
                            "ad_coord ad_vel ad_dichte ad_druck"
#endif
                            )),

    boundaryFile(GetEnvString("SPH_PARAM_BOUNDARY_FILE", "")),
    boundaryDerivFilePattern(GetEnvString("SPH_PARAM_BOUNDARY_DERIV_FILE", "")),

    sceneOutFile(GetEnvString("SPH_PARAM_SCENE_OUTFILE", "sph-initial")),
    sceneInFile(GetEnvString("SPH_PARAM_SCENE_INFILE", "")),
    logFile(GetEnvString("SPH_PARAM_LOG_FILE", sceneOutFile + ".sphscene.xml")),

    initialVTKScale(GetEnv("SPH_PARAM_INITIAL_SCALE", 1)),
    boundaryScale(GetEnv("SPH_PARAM_BOUNDARY_SCALE", 1)),
    
    initialVTKTranslate(GetEnvVec("SPH_PARAM_INITIAL_TRANSLATE", 0, -1e100, 1e100)),
    boundaryTranslate(GetEnvVec("SPH_PARAM_BOUNDARY_TRANSLATE", 0, -1e100, 1e100)),
    
    fieldsToLoad(GetEnvString("SPH_PARAM_SAVE_FIELDS", "masse dichte waerme coord vel class color")),
    resultFormat(GetEnvString("SPH_PARAM_RESULT_FORMAT", "vtk"))

  {
    assert(initDensity > 0);
    assert(deltaX > 0);
    // assert(boundaryRes >= 1);
  }
  void writeXML(std::ostream &aus, std::string const &indent = "") const {
    aus << indent << "<scene-config\n" << indent <<
#ifdef _OPENMP
      " openmp='yes'"
#else
      " openmp='no'"
#endif

#ifdef SPH_SAVE_H5PART
      " save='yes'"
#else
      " save='no'"
#endif
#ifdef NDEBUG
      " assertions='no'"
#else
      " assertions='yes'"
#endif
      "\n" << indent <<
      " threads='" << numThreads << "'"
      " dx='" << deltaX << "'"
      " boundaryRes='" << boundaryRes << "'"
      " h-coeff='" << kernelReachCoeff << "'"
      " rl='" << initReplacementLength << "'"
      " d='" << ndim << "'"
      " randomSeed='" << randomSeed << "'"
      " na='" << numInitAreas << "'"
      " scene-ndim='" << numDims << "'"
      " n='" << numInitParticles << "'"
      " initeps='" << initEpsilon << "'\n" << indent <<
      " initrho='" << initDensity << "'\n" << indent <<
      " boundary='" << boundaryFile << "'\n" << indent <<
      " boundaryDerivFilePattern='" << boundaryDerivFilePattern << "'\n" << indent <<
      " sceneOutFile='" << sceneOutFile << "'"
      " logFile='" << logFile << "'"
      " sceneInFile='" << sceneInFile << "'"
      " loadFields='" << fieldsToLoad << "'\n" << indent <<
      " saveASCII='" << saveASCII << "'"
      ">\n";
    eqsConfig.writeXML(aus, indent + "\t");

    gebiet.writeXML(aus, indent + "\t");

//     aus << "<gravity g='" << gravity
// 	<< "' x='" << up[0]
// 	<< "' y='" << (ndim > 1 ? up[1] : 0)
// 	<< "' z='" << (ndim > 2 ? up[2] : 0) << "'/>\n";

    aus << indent << "</scene-config>\n";
  }

  static size_t getAreaNDim(Gebiet const &gebietInit) {
    size_t n = 0;
    for (size_t i = 0; i < ndim; ++i) {
      if (gebietInit.v[i] != 0) {
        ++n;
      }
    }
    return n;
  }

};

#endif
