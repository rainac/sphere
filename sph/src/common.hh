#ifndef jw_sph_common_324154_hh
#define jw_sph_common_324154_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

// \def SPH_RENDER
// Option: visualize relations: distance tests (green) and forces actually computed (blue)
// #define SPH_RENDER

// \def DEBUG_SPH 
// Compile in debug messages. 
// #define DEBUG_SPH


/*! Macro to set the number of dimension the program uses.
  Sensible values are 1, 2 or 3. The macro is (and should only be) used in
  the OpenGL renderer code.
 */
#ifndef SPH_DIM
#define SPH_DIM 3
#endif

// #ifndef SPH_DATA_TYPE
// #error define SPH_DATA_TYPE
// #endif

#ifndef SPH_ACTIVE_DATA_TYPE
#define SPH_ACTIVE_DATA_TYPE double
#endif

// typedef SPH_DATA_TYPE              SphDataType;
typedef SPH_ACTIVE_DATA_TYPE      SphADataType;


// \def SPH_SHOW_LOAD_AVERAGE
// macro to output thread workload information on at each step
// #define SPH_SHOW_LOAD_AVERAGE

// \def SPH_SHOW_LOAD_PER_THREAD
// macro to output thread workload information on at each step
// instead of just average and variance of load per thread
// #define SPH_SHOW_LOAD_PER_THREAD

/// \def SPH_PARTICLE_SAVES_RASTER_INDEX
/// if the D-dim. integer raster index is to be save in each particle.
/// \deprecated this is not an option since it must be enabled
#define SPH_PARTICLE_SAVES_RASTER_INDEX



/*! Global constant holds the number of dimension the program uses.
  Sensible values are 1, 2 or 3. The value should come from macro SPH_DIM.
 */
static unsigned const ndim = SPH_DIM;


// /// constant for 
// #define SPH_INTEGRATOR_EXPLICIT 1
// #define SPH_INTEGRATOR_HEUNM2   2
// #define SPH_INTEGRATOR_HEUNM3   3
// #define SPH_INTEGRATOR_SIMPSON  4
// #define SPH_INTEGRATOR_RK_4_3_8 5
// #define SPH_INTEGRATOR_ABM_1    6
// #define SPH_INTEGRATOR_ABM_2    7
// #define SPH_INTEGRATOR_ABM_3    8
// #define SPH_INTEGRATOR_PK_1     9
// #define SPH_INTEGRATOR_PK_2     10

// #if !defined SPH_INTEGRATOR
// /// select default integrator here
// // #define SPH_INTEGRATOR SPH_INTEGRATOR_EXPLICIT
// #define SPH_INTEGRATOR SPH_INTEGRATOR_PK_1
// #endif

/// use std::valarray for particle arrays instead of own array classes
// #define SPH_USE_VALARRAY_FOR_PARTICLE_ARRAYS

#define SUM_FORCES_2BINS_FAST
#define SPH_USE_INTEGRATOR

#ifdef SPH_DEBUG
// if no removal of assertions is requested ...

/// ... let particles have ids to add more checks
#define SPH_PARTICLE_HAS_ID

#ifdef SPH_PARTICLE_HAS_ID
/// ... let particles have sum of ids of particles that contributed to it
#define SPH_PARTICLE_VAR_HAS_ID_SUM
#endif

/// ... let particles have sum of particles that contributed to it
#define SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM

#endif

// #define SPH_PARTICLE_HAS_ID
// #define SPH_PARTICLE_VAR_HAS_ID_SUM
// #define SPH_PARTICLE_VAR_HAS_NEIGHBOUR_SUM


#define SPH_PARTICLE_MAX_VELOCITY

// requires SPH_PARTICLE_MAX_VELOCITY
#define SPH_ADAPTIVE_TIMESTEP


#ifdef __SUNPRO_CC_COMPAT
#define SPH_COMPILER_SUN_STUDIO
#endif

#ifdef DEBUG
#define SPH_DEBUG
#define SPH_HARD_ASSERTIONS
#endif

#ifdef SPH_HARD_ASSERTIONS
#define sph_hard_assert(x) assert(x)
#else 
#define sph_hard_assert(x)
#endif

/// enable periodic domain boundaries
#define SPH_PERIODIC_BOUNDARIES

/// if AD enabled...
#ifdef SPH_AD
/// clear derivatives of position() of moving particles in each step
// #define SPH_AD_CLEAR_D_POSITION_EACH_STEP

/// ensure position() of moving particles never get a derivative
// #define SPH_AD_KILL_DERIVATIVES_OF_POSITION_UPDATE

// #define SPH_AD_DEACTIVATE_KERNEL

#endif // AD

#ifdef __INTEL_COMPILER
#define SPH_COMPILER_INTEL

#elif defined __GNUC__
#define SPH_COMPILER_GCC

#if __GNUC__ == 4 and __GNUC_MINOR__  == 3
#define SPH_COMPILER_GCC_43
#endif
#if __GNUC__ == 4 and __GNUC_MINOR__  == 4
#define SPH_COMPILER_GCC_44
#endif

#elif defined __SUNPRO_CC
#define SPH_COMPILER_STUDIO

#endif

#if defined SPH_COMPILER_GCC_44 || defined SPH_COMPILER_INTEL || defined SPH_COMPILER_STUDIO
#define SPH_HAVE_feenableexcept
#endif


#endif
