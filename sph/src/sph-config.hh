#ifndef sph_config_hh
#define sph_config_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

#include "common.hh"
#include "quader.hh"
#include "utility/getenv.hh"
#include "vector.hh"
#include "sensor-def.hh"

using std::max;

/// Base class SPHConfigBase of the configuration structs. It only
/// defines two typedefs the subclasses use.

struct SPHConfigBase {
  /// The type to get a ndim-dimensional value from the environment.
  typedef GetEnvV<ndim, Vector> GetEnvVec;

  /// The type to get a 3-dimensional value from the environment (used
  /// for renderer config).
  typedef GetEnvV<3, Point<double, 3> > GetEnvVec3D;

  typedef Quader<Point<double, ndim> > Gebiet;

  void getNthEnvString(std::string const &prefix, size_t num, std::string const &name, std::string &val) {
    std::ostringstream str;
    str << prefix << num << "_" << name;
    std::string const envnamev = str.str();
    val = std::string(GetEnvString(envnamev, val));
  }

  template<class T>
  void getNthEnvValue(std::string const &prefix, size_t num, std::string const &name, T &val) {
    std::ostringstream str;
    str << prefix << num << "_" << name;
    std::string const envnamev = str.str();
    val = GetEnv(envnamev, val, 0, 1e100);
  }

  void getNthArea(std::string const &prefix, size_t num, std::string const &name, Gebiet &gebiet) {
    std::ostringstream str;
    str << prefix << num << "_" << name << "_P1";
    std::string envnamep1 = str.str();
    std::ostringstream str2;
    str2 << prefix << num << "_" << name << "_P2";
    std::string envnamep2 = str2.str();
    
    // std::cerr << "getenv " << envnamep1 << " and " << envnamep2 << "\n";
    if (getenv(envnamep1.c_str()) != 0 or getenv(envnamep2.c_str()) != 0) {
      if (getenv(envnamep1.c_str()) == 0) {
        std::cerr << "error: expected environment variable " << envnamep1 << ", but it is not set\n";
        exit(-5);
      }
      if (getenv(envnamep2.c_str()) == 0) {
        std::cerr << "error: expected environment variable " << envnamep2 << ", but it is not set\n";
        exit(-5);
      }
      
      gebiet = Gebiet(GetEnvV<ndim, DoubleVector >(envnamep1, 0, -1e100, 1e100),
                      GetEnvV<ndim, DoubleVector >(envnamep2, 1, -1e100, 1e100));
      // std::cerr << "area: " << gebiet << "\n";
    }
  }
};


/// Base class SPHEQSConfig of the configuration structs. It combines
/// the global constants needed in the equation of state (EQS) for the
/// pressure.
template<class ValueType=double>
struct SPHEQSConfig {
  typedef ValueType value_type;

  /// The density of the medium.  
  value_type rho0;

  /// The exponent in the equation of state (EQS) for pressure. It
  /// should be uneven.
  double gamma;

  /// The factor in the equation of state (EQS) for pressure.
  /// B is for the Monaghan Taits equation (that with gamma=7 :)
  /// K0 is the coefficient as in the Hayward, Klaus, O'Brian, ... Tais equation (linear)
  value_type pressureB,  pressureK0;

  SPHEQSConfig() :
    rho0(GetEnv("SPH_PARAM_RHO0", 1e3, 1e-100, 1e100)),
    gamma(GetEnv("SPH_PARAM_GAMMA", 7, 1e-100, 1e100)),
    pressureB(GetEnv("SPH_PARAM_B", 5e5, 1e-100, 1e100)),
    pressureK0(GetEnv("SPH_PARAM_K0", 3.4e6, 1e-100, 1e100))
  {
    assert(yafad_value(pressureB) > 0);
    assert(yafad_value(rho0) > 0);
  }

  void writeXML(std::ostream &aus, std::string const &indent = "") const {
    aus << indent << "<eqs type='monaghan'"
      " rho0='" << rho0 << "'"
      " gamma='" << gamma << "'"
      " B='" << pressureB << "'"
      " K0='" << pressureK0 << "'"
      "/>\n";
  }

};

struct SPHOpenMPConfig {
#if _OPENMP > 200804
  std::string const ompSchedule;
  long ompChunkSize;
#else
  std::string const ompSchedule;
#endif
  bool ompDynamic;
  bool ompNested;
//   long ompNumThreads, ompMaxNumThreads;

#if _OPENMP > 200804
#include "schedule-types.ncd.enum.hh"
#include "schedule-types.ncd.cc"
#endif
  
  SPHOpenMPConfig() :
#if _OPENMP > 200804
    ompSchedule(GetEnvString("SPH_OMP_SCHEDULE", "dynamic")),
    ompChunkSize(GetEnv("SPH_OMP_CHUNK_SIZE", 0, 0, 1e100)),
#else
    ompSchedule(GetEnvString("OMP_SCHEDULE")),
#endif
    ompDynamic(GetEnv("SPH_OMP_DYNAMIC", 0, 0, 1)),
    ompNested(GetEnv("SPH_OMP_NESTED", 0, 0, 1))
//     ompNumThreads(GetEnv("SPH_PARAM_THREADS", 4, 1, 1e100)),
//     ompMaxNumThreads(GetEnv("SPH_PARAM_MAX_THREADS", 256, 1, 1e100)),
  {
  }

  void writeXML(std::ostream &aus, std::string const &indent = "") const {
    aus << indent << "<openmp"
      " schedule='" << ompSchedule << "'"
#if _OPENMP > 200804
      " chunk-size='" << ompChunkSize << "'"
#endif
      " dynamic='" << ompDynamic << "'"
      " nested='" << ompNested << "'"
      /*      " threads='" << ompNumThreads << "'"
      " maxThreads='" << ompMaxNumThreads << "'"
      */
      "/>\n";
  }
};

size_t const maxNumSensors = 16;

/// Gathers all configuration constants.
/// Can write itself to a stream in XML format.
template<class ValueType=double>
struct SPHConfig : public SPHConfigBase {
  
  typedef ValueType value_type;

  SPHEQSConfig<value_type> eqsConfig;

  SPHOpenMPConfig openMPConfig;

  /// The number of threads. 
  unsigned const numThreads, numCSortThreads;

  /// The maximum number of threads (if interactive rendering is enabled).
  unsigned const maxNumThreads;

  /// The computing domain.
  Gebiet gebiet;

#ifdef SPH_PERIODIC_BOUNDARIES
  /// Which boundaries are periodic (wrap around)
  /// 0 not periodic(default), 1 periodic
  DoubleVector periodicBoundaries;
#endif

  /// The spatial rastering length (same for all dimensions).
  double const kernelReachH;

  /// The name of the kernel function.
  /// Legal values include "gauss", "spline", "wendland0", "wendland1", "wendland2", "wendland3"
  std::string const kernelName;

  /// The name of summation algorithm.
  /// Legal values include "naive", "symmetric", "sym-mb", "sym-mb-tb"
  std::string const summationAlgorithm;

  /// The name of summation algorithm.
  /// Legal values include "csort", "stdsort_ser_index", "stdsort_par_index"
  /// cf. sorted-index-types.ncd.xml
  std::string const sortedIndex;

  /// The total simulation time \f$ t_{\text{max}} \f$.
  double const tmax;

  /// The temporal timestep \f$ \Delta t \f$.
  double const timestepDefault;

#ifdef SPH_ADAPTIVE_TIMESTEP
  /// The maximum temporal timestep \f$ {\Delta t}_{\text{max}} \f$.
  double const timestepMax;

  /// The CFL number to be maintained by adapting the time step.
  double const cflNumber;
#endif

  // /// The number of simulation steps.
  // size_t const numSteps;

  /// The number of moving particles left keeps simulation going
  size_t const minParticlesLeft;

  /// The gravitational constant
  value_type const gravity;

  /// The dynamicBulkViscosity of the medium
  value_type dynamicBulkViscosity;

  /// The speed of sound in the medium
  value_type const speedOfSound;

  /// The constant \f$ \epsilon\f$ in the XSPH modification.
  value_type const xsphEpsilon;

  /// The constant \f$ \eta \f$ in the viscosity term.
  value_type viscEta;
  /// The constant \f$ \alpha \f$ in the viscosity term.
  value_type viscAlpha;
  /// The constant \f$ \beta \f$ in the viscosity term.
  value_type viscBeta;

  /// The constant \f$ r_0 \f$ in the Lennard-Jones force. This
  /// determines the radius of the repulsive potential.
  value_type const lj_r0;
  /// The constant \f$ D \f$ in the Lennard-Jones force. This
  /// determines the strength of the repulsive potential (it goes to
  /// infinity as the distance approaches zero anyways).
  value_type const lj_D;
  /// The constant \f$ p_2 \f$ in the Lennard-Jones force.  The
  /// constant \f$ p_1 \f$ is hard-wired to \f$ 2 p_2 \f$. These two
  /// determine the slope of the repulsive potential.
  unsigned const lj_p2;

  /// The direction of gravity
  Vector              down;

  /// The gravity acceleration vector.
  Vector              gravityVector;

  /// This determines how often to save the state to file (<= 0 is off).  
  double const saveRate;
  
  /// Save particles that are outside (have left) the simulation domain?
  bool const saveOut;

  /// Save boundary particles in every step (not only first)?
  bool const saveBoundary;

  /// This integer determines how often to render the state with the
  /// builtin OpenGL renderer (0 is off).
  unsigned const renderEvery;

  /// This integer determines how often to sleep for sleepInterval (0 is off).
  unsigned const sleepEvery;

  /// This integer determines how often to sleep for sleepInterval.
  unsigned const sleepInterval;

  /// Print out more XML elements with some timings of parts of the algorithm
  bool const detailedTimings;
  
  /// If a FP exception has been raised during the last time step,
  /// abort the program.
  bool const abortOnFExcept;

#ifdef SPH_HAVE_feenableexcept
  /// Enable trapping of FP errors: When an illegal FP instruction is
  /// excuted SIGFPE is raised.
  bool const trapFloatingPointErrors;
#endif

  bool const saveASCII;

  /// The number of bits to consider in the HashFunktion. 
  int const hashFunctionM;

  /// This integer holds the number of sensor areas
  unsigned const numSensors;
  /// The default sensor type
  std::string const sensorTypeDefault;

  SensorDef sensors[maxNumSensors];

#ifdef SPH_AD
  /// max difference of derivatives on discontinuity of max(.,.)
  double const discontMaxBranchesMaxDelta;
#endif

  /// The name of the integrator algorithm.
  std::string const integratorName;

  /// The name of the file to load the colors from (used in renderer).
  std::string const colorMapFile;
  /// The name of the file to load the colors from (used in renderer).
  long const colorMapSize;

  /// The name of the file to save the result to.
  std::string const resultFile;

  /// The name of the file to load the initial state from.
  std::string const initialFile;

  /// The name of the file to log output to.
  std::string const logFile;

  std::string const particleTypeMovingDefaultName;

  std::string const particleTypeBoundaryDefaultName;

  /// The format to use for writing the result file (VTK or H5Part/HDF5).
  std::string const resultFormat;

#ifdef SPH_RENDER
  double initRotationAngle;
  Point<double, 3> initRotationAxis;

  //! initial Camera position
  Point<double, 3> initCamPos;
//   //! initial Camera target and camera up vector
//   Point<double, 3> initCamTarget;

  /// factor to scale the GL rendered letters
  double const fontSize;
#endif

  /// The name of the file to save the result to.
  std::string const allFields;

  /// The particle data fields to save. If the field's name is present
  /// as a substring, the field is saved.
  std::string saveFields;

  /// The particle data fields to save. If the field's name is present
  /// as a substring, the field is saved.
  std::string fieldsToLoad;

  /// Initializes all the constants to the values given by defaults or
  /// set via the environment variables. These have the same name, but
  /// all uppercased and prefixed with "SPH_PARAM_".
  SPHConfig() :
#ifdef _OPENMP
    numThreads(GetEnv("SPH_PARAM_THREADS", 4, 1, 1e100)),
    numCSortThreads(GetEnv("SPH_PARAM_CSORT_THREADS", numThreads, 1, 1e100)),
    maxNumThreads(GetEnv("SPH_PARAM_MAX_THREADS", 256, 1, 1e100)),
#else
    numThreads(1),
    numCSortThreads(1),
    maxNumThreads(1),
#endif
    gebiet(DoubleVector(), 
	   GetEnvV<ndim, DoubleVector >("SPH_PARAM_QUADER_P2", 1)),
#ifdef SPH_PERIODIC_BOUNDARIES
    periodicBoundaries(GetEnvV<ndim, DoubleVector >("SPH_PARAM_PERIODIC_BOUNDARIES", 0)),
#endif
    kernelReachH(GetEnv("SPH_PARAM_H", 1e-2, 1e-100, 1e100)),
    kernelName(GetEnvString("SPH_PARAM_KERNEL", "wendland2")),
    summationAlgorithm(GetEnvString("SPH_PARAM_SUMMATION", "symmetric")),
    sortedIndex(GetEnvString("SPH_PARAM_SORTED_INDEX", "csort")),
    tmax(GetEnv("SPH_PARAM_TMAX", 3, 1e-100, 1e100)),
    timestepDefault(GetEnv("SPH_PARAM_DT", 1e-4, 1e-100, 1e6)),
#ifdef SPH_ADAPTIVE_TIMESTEP
    timestepMax(GetEnv("SPH_PARAM_DT_MAX", timestepDefault, 1e-100, 1e6)),
    cflNumber(GetEnv("SPH_PARAM_CFL", 1, 1e-100, 1e100)),
#endif
    // numSteps(GetEnv("SPH_PARAM_NT", ceil(tmax/timestepDefault), 0, 1e100)),
    minParticlesLeft(GetEnv("SPH_PARAM_MIN_LEFT", 1, 0, 1e100)),
    gravity(GetEnv("SPH_PARAM_G", 9.81, 0, 1e100)),
    dynamicBulkViscosity(GetEnv("SPH_PARAM_MU", 0, 0, 1e100)),
    speedOfSound(GetEnv("SPH_PARAM_SPEED_OF_SOUND", 1484, 0, 1e100)),
    xsphEpsilon(GetEnv("SPH_PARAM_EPSILON", 0, 0, 1e100)),
    viscEta(GetEnv("SPH_PARAM_ETA", 0, 0, 1e100)),
    viscAlpha(GetEnv("SPH_PARAM_ALPHA", 0.01, 0, 1e100)),
    viscBeta(GetEnv("SPH_PARAM_BETA", 0, 0, 1e100)),
    lj_r0(GetEnv("SPH_PARAM_LJ_R0", kernelReachH/4.0, 1e-100, 1e100)),
    lj_D(GetEnv("SPH_PARAM_LJ_D", 1e2
		/* 5*gravity*((-up)*gebiet.v).sum() */, 1e-100, 1e100)),
    lj_p2(GetEnv("SPH_PARAM_LJ_P2", 2, 2, 1e100)),
    down(ndim),
    gravityVector(ndim),

    saveRate(GetEnv("SPH_PARAM_SAVE_RATE", 25, 0, 1e100)),

    saveOut(GetEnv("SPH_PARAM_SAVE_OUT", 1, 0, 1)),
    saveBoundary(GetEnv("SPH_PARAM_SAVE_BOUNDARY", 1, 0, 1)),
    renderEvery(GetEnv("SPH_PARAM_RENDER_EVERY", 10, 0, 1e100)),
    sleepEvery(GetEnv("SPH_PARAM_SLEEP_EVERY", 0, 0, 1e100)),
    sleepInterval(GetEnv("SPH_PARAM_SLEEP_TIME", 1, 0, 1e100)),
    detailedTimings(GetEnv("SPH_PARAM_TIMINGS", 0, 0, 10)),
    abortOnFExcept(GetEnv("SPH_PARAM_ABORT_ON_FE", 1, 0, 1)),
#ifdef SPH_HAVE_feenableexcept
    trapFloatingPointErrors(GetEnv("SPH_PARAM_TRAP_FE", 0, 0, 1)),
#endif
    saveASCII(GetEnv("SPH_PARAM_SAVE_ASCII", 0, 0, 1)),
    hashFunctionM(GetEnv("SPH_PARAM_M", 6, 1, 12)),
    numSensors(GetEnv("SPH_PARAM_SENSORS", 0, 0, maxNumSensors)),
    sensorTypeDefault(GetEnvString("SPH_PARAM_SENSOR_TYPE", "VELOCITY")),
    sensors(),
#ifdef SPH_AD
    discontMaxBranchesMaxDelta(GetEnv("SPH_PARAM_DISCONT_MAX_BRANCHES_MAX_DELTA", 1e-14, 0, 1)),
#endif

    integratorName(GetEnvString("SPH_PARAM_INTEGRATOR", "pk1_1")),

    colorMapFile(GetEnvString("SPH_PARAM_COLORS_FILE", "colors.txt")), 
    colorMapSize(GetEnv("SPH_PARAM_NUM_COLORS", 100, 0, 1e100)), 

    resultFile(GetEnvString("SPH_PARAM_RESULT_FILE", "sph-result")),
    initialFile(GetEnvString("SPH_PARAM_INITIAL_FILE", "sph-initial.vtk")),
    logFile(GetEnvString("SPH_PARAM_LOG_FILE", resultFile + ".sphrun.xml")),

    particleTypeMovingDefaultName(GetEnvString("SPH_PARAM_MOVING_PARTICLE", "ferrari")),
    particleTypeBoundaryDefaultName(GetEnvString("SPH_PARAM_BOUNDARY_PARTICLE", "point_mirrored")),
    resultFormat(GetEnvString("SPH_PARAM_RESULT_FORMAT", "hdf5")),

#ifdef SPH_RENDER
    initRotationAngle(GetEnv("SPH_PARAM_ROTATION_ANGLE", 0, -1e100, 1e100)),
    initRotationAxis(GetEnvVec3D("SPH_PARAM_ROTATION_AXIS", 0)),
    initCamPos(3),
//     initCamUpVector(3),
//     initCamTarget(3),
    
    fontSize(GetEnv("SPH_PARAM_FONTSIZE", 1, 0, 1e100)),
#endif

    allFields("masse dichte druck waerme class color coord vel hash flags sensors"
#ifdef SPH_PARTICLE_HAS_ID
              " id"
#endif
#ifdef SPH_PARTICLE_VAR_HAS_ID_SUM
              " sumOfNeighbourIds"
#endif
#ifdef SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM
              " neighbours"
#endif
#ifdef SPH_AD
              " ad_coord ad_vel ad_dichte ad_druck ad_waerme"
#endif
              ),
    // all fields: masse dichte druck waerme class color coord vel hash flags
    saveFields(GetEnvString("SPH_PARAM_SAVE_FIELDS", "dichte waerme druck coord vel class color"
#ifdef SPH_PARTICLE_HAS_ID
                            " id"
#endif
#ifdef SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM
                            " neighbours"
#endif
#ifdef SPH_PARTICLE_VAR_HAS_ID_SUM
                            " sumOfNeighbourIds"
#endif
#ifdef SPH_AD
                            " ad_coord ad_vel ad_dichte ad_druck ad_waerme"
#endif
                            )),
    fieldsToLoad(GetEnvString("SPH_PARAM_LOAD_FIELDS", "masse dichte waerme coord vel class color"
#ifdef SPH_AD
                              " ad_coord ad_vel ad_dichte ad_druck ad_waerme"
#endif
                              ))

  {
    if (getenv("SPH_PARAM_QUADER_P1")) {
      std::cerr << "warning: env variable SPH_PARAM_QUADER_P1 is set but will be"
        " ignored (deprecated)\n";
    }
    assert(kernelReachH > 0);
#ifdef SPH_RENDER
    Point<double, 3> gul(gebiet.untenLinks.data, ndim);
    Point<double, 3> gv(gebiet.v.data, ndim);
    initCamPos = gul + (gv/2.0);
    initCamPos[2] += norm(gv) * 2.1;
    initCamPos = GetEnvVec3D("SPH_PARAM_CAM_POS", initCamPos, -1e100, 1e100);
#endif
#if SPH_DIM == 1
    down[0] = -1;
#else
    down[1] = -1;
#endif
    down = GetEnvVec("SPH_PARAM_DOWN", down, -1e100, 1e100);
    down = down / norm(down);
    gravityVector = down * gravity;

    std::string const sensorDefPrefix = "SPH_PARAM_SENSOR_S";
    for(unsigned i = 0; i < numSensors; ++i) {

      Gebiet siGebiet;
      getNthArea(sensorDefPrefix, i, "POINTS", siGebiet);

      Gebiet siGebiet2;
      getNthArea(sensorDefPrefix, i, "VECTORS", siGebiet2);

      std::string siTypeName = sensorTypeDefault;
      getNthEnvString(sensorDefPrefix, i, "TYPE", siTypeName);
      //       std::cerr << "stypename: " << siTypeName << "\n";
      int const siType = SensorDef::getSensorTypesValue(siTypeName.c_str());
      //       std::cerr << "stypename: " << SensorDef::getSensorTypesName(siType) << "\n";
      sensors[i] = SensorDef(siGebiet, siType, siGebiet2);
      //       std::cerr << "sensor: " << sensors[i].m_type << "\n";

      std::string siAreaName = "ORTHOGONAL_RECTANGLE";
      getNthEnvString(sensorDefPrefix, i, "SHAPE", siAreaName);
      // std::cerr << "sareaname: " << siAreaName << "\n";
      int const siArea = SensorDef::getSensorAreasValue(siAreaName.c_str());
      // std::cerr << "stypename: " << SensorDef::getSensorAreasName(siArea) << "\n";
      sensors[i].m_areaType = siArea;
      // std::cerr << "sensor area: " << sensors[i].m_areaType << "\n";

      int siLog = 0;
      getNthEnvValue(sensorDefPrefix, i, "LOG", siLog);
      sensors[i].m_log = siLog > 0;
      // std::cerr << "sensor-def " << i << ": " << sensors[i] << "\n";
      // std::cerr << "sensor: " << sensors[i].m_type << "\n";
    }

#ifdef SPH_AD
#if SPH_AD == 1
    /// \todo bessere loesung fuer beliebige NDIR 
    dynamicBulkViscosity.diff(0) = GetEnv("SPH_PARAM_MU_D0");
    dynamicBulkViscosity.diff(1) = GetEnv("SPH_PARAM_MU_D1");
    dynamicBulkViscosity.diff(2) = GetEnv("SPH_PARAM_MU_D2");
    viscAlpha.diff(0) = GetEnv("SPH_PARAM_ALPHA_D0");
    viscAlpha.diff(1) = GetEnv("SPH_PARAM_ALPHA_D1");
    viscAlpha.diff(2) = GetEnv("SPH_PARAM_ALPHA_D2");
    eqsConfig.rho0.diff(0) = GetEnv("SPH_PARAM_RHO0_D0");
    eqsConfig.rho0.diff(1) = GetEnv("SPH_PARAM_RHO0_D1");
    eqsConfig.rho0.diff(2) = GetEnv("SPH_PARAM_RHO0_D2");
    eqsConfig.pressureB.diff(0) = GetEnv("SPH_PARAM_B_D0");
    eqsConfig.pressureB.diff(1) = GetEnv("SPH_PARAM_B_D1");
    eqsConfig.pressureB.diff(2) = GetEnv("SPH_PARAM_B_D2");
#endif
#endif

  }

  void writeXML(std::ostream &aus) const {
    aus << "<config\n"
      " ndim='" << ndim << "'"
#ifdef _OPENMP
      " openmp='yes'"
#else
      " openmp='no'"
#endif

#ifdef SPH_AD
      " ad='yes'"
#ifdef SPH_AD_DEACTIVATE_KERNEL
      " deactivate-kernel='yes'"
#else
      " deactivate-kernel='no'"
#endif
#ifdef SPH_AD_KILL_DERIVATIVES_OF_POSITION_UPDATE
      " kill-derivatives-of-position-update='yes'"
#else
      " kill-derivatives-of-position-update='no'"
#endif
#ifdef SPH_AD_CLEAR_D_POSITION_EACH_STEP
      " kill-derivatives-of-position-each-step='yes'"
#else
      " kill-derivatives-of-position-each-step='no'"
#endif
#else
      " ad='no'"
#endif

#ifdef SPH_RENDER
      " render='yes'"
#else
      " render='no'"
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
      "\n"
#ifdef SPH_PARTICLE_MAX_VELOCITY
      " max-velocity='yes'"
#endif
#ifdef SPH_ADAPTIVE_TIMESTEP
      " adaptive-timestep='yes'"
#endif
#ifdef SPH_USE_INTEGRATOR
      " time-integrator='yes'"
#endif
#ifdef SPH_PARTICLE_HAS_ID
      " particle-id='yes'"
#endif
#ifdef SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM
      " particle-neighbours='yes'"
#endif
      "\n"
      " threads='" << numThreads << "'\n"
      " csort-threads='" << numCSortThreads << "'\n"
      " tmax='" << tmax << "'"
      " dt='" << timestepDefault << "'"
#ifdef SPH_ADAPTIVE_TIMESTEP
      " dt-max='" << timestepMax << "'"
      " cfl='" << cflNumber << "'"
#endif
      // " nt='" << numSteps << "'\n"
#ifdef SPH_PERIODIC_BOUNDARIES
      " periodicBoundaries='" << periodicBoundaries << "'"
#endif
      " H='" << kernelReachH << "'"
      " kernel='" << kernelName << "'"
      " summation='" << summationAlgorithm << "'"
      " sorted-index='" << sortedIndex << "'"
      " m='" << hashFunctionM << "'\n"
      " numSensors='" << numSensors << "'"
      " sensorTypeDefault='" << sensorTypeDefault << "'"
#ifdef SPH_AD
      " discontMaxBranchesMaxDelta='" << discontMaxBranchesMaxDelta << "'"
#endif
      " d='" << ndim << "'"
      " mu='" << dynamicBulkViscosity << "'"
      " c='" << speedOfSound << "'"
      " eps='" << xsphEpsilon << "'\n"
      " saveRate='" << saveRate << "'"
      " saveOut='" << saveOut << "'"
      " saveBoundary='" << saveBoundary << "'"
      " renderEvery='" << renderEvery << "'\n"
      " sleepEvery='" << sleepEvery << "'"
      " sleepInterval='" << sleepInterval << "'\n"
      " detailedTimings='" << detailedTimings << "'"
      " abortOnFExcept='" << abortOnFExcept << "'\n"
#ifdef SPH_HAVE_feenableexcept
      " trapFloatingPointErrors='" << trapFloatingPointErrors << "'\n"
#endif
      " result='" << resultFile << "'\n"
      " initial='" << initialFile << "'\n"
      " integrator='" << integratorName << "'\n"
      " color-map-file='" << colorMapFile << "'"
      " color-map-size='" << colorMapSize << "'\n"
      " allFields='" << allFields << "'\n"
      " particleTypeMovingDefaultName='" << particleTypeMovingDefaultName << "'\n"
      " particleTypeBoundaryDefaultName='" << particleTypeBoundaryDefaultName << "'\n"
      " saveFields='" << saveFields << "'\n"
      " loadFields='" << fieldsToLoad << "'\n"
      ">\n";

    eqsConfig.writeXML(aus);
    openMPConfig.writeXML(aus);
    gebiet.writeXML(aus);

    aus << "<sensors num='" << numSensors << "'>\n";
    for(size_t i = 0; i < numSensors; ++i) {
      sensors[i].writeXML(aus, "\t");
    }
    aus << "</sensors>\n";

    aus << "<viscosity alpha='" << viscAlpha
	<< "' beta='" << viscBeta 
	<< "' eta='" << viscEta << "'/>\n"; 
    aus << "<lennard-jones r0='" << lj_r0
	<< "' D='" << lj_D
	<< "' p2='" << lj_p2 << "'/>\n"; 
    aus << "<gravity g='" << gravity
	<< "' x='" << down[0]
	<< "' y='" << (ndim > 1 ? down[1] : 0)
	<< "' z='" << (ndim > 2 ? down[2] : 0) << "'/>\n";

#ifdef SPH_RENDER
    aus << "<camera>\n"
      "\t<position x='" << initCamPos[0]
	<< "' y='" << initCamPos[1]
	<< "' z='" << initCamPos[2] << "'/>\n"
//       "\t<target x='" << initCamTarget[0]
// 	<< "' y='" << initCamTarget[1]
// 	<< "' z='" << initCamTarget[2] << "'/>\n"
//       "\t<up x='" << initCamUpVector[0]
// 	<< "' y='" << initCamUpVector[1]
// 	<< "' z='" << initCamUpVector[2] << "'/>\n"
	<< "</camera>\n";
#endif
    aus << "</config>\n";
  }
};

#endif
