#
# Example settings file for SPH scene compiler
#

# GLOBAL PARAMETERS
#
# name of result file (stem only)
export SPH_PARAM_RESULT_FILE="beispiel-szene"

# result file format (and suffix)
export SPH_PARAM_RESULT_FORMAT="vtk h5"
#export SPH_PARAM_RESULT_FORMAT="vtk"
#export SPH_PARAM_RESULT_FORMAT="h5"

# read a boundary particle coordinates from a text file
#export SPH_PARAM_BOUNDARY_FILE="boundary.txt"


# SETTINGS FOR ALL "AREAS"
#

# The equation of state parameters are needed to calculate the initial
# pressure based on initial density. However the SPH code loading the
# particles should recompute the pressure based on it's own settings!
export SPH_PARAM_B="1e5"
export SPH_PARAM_GAMMA="7"

# default disturbance of grid positions
export SPH_PARAM_INIT_EPSILON="1e-10"

# default number of particles per "area"
export SPH_PARAM_N="1e-10"

# default initial density of particles in "areas"
export SPH_PARAM_INIT_DENSITY="1.001e3"

# default initial heat energy
export SPH_PARAM_INIT_HEAT="2.22222222222222222222222222222222222222"

# if an area is a line or plane in 3D or a line in 2D, what shall be
# the thickness? this influences the volume a line or plane stands for
# and thereby the mass per particle.  it's not sensible to set this
# much larger or smaller than the kernel reach length H
export SPH_PARAM_INIT_REPLACEMENT_LENGTH="1e-2"


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
# number of "areas" to initialize: 3
export SPH_PARAM_INIT_NUM="2"

# settings for first "area" (a line)
# number of particles
export SPH_PARAM_INIT_A0_N="10"
# velocity
export SPH_PARAM_INIT_A0_V="7 0.4 7"
# free-floating or resting on bottom
export SPH_PARAM_INIT_A0_FREE="0"
# heat
export SPH_PARAM_INIT_A0_HEAT="3"
# heat
export SPH_PARAM_INIT_A0_EPSILON="3"
# define area
export SPH_PARAM_INIT_A0_P1="0.3  0.405  0.1"
export SPH_PARAM_INIT_A0_P2="0.5  0.405  0.1"

# settings for second "area" (a plane)
export SPH_PARAM_INIT_A1_N="10"
export SPH_PARAM_INIT_A1_V="1 1 1"
export SPH_PARAM_INIT_A1_HEAT="10"
export SPH_PARAM_INIT_A1_P1="0.3  0.605  0.1"
export SPH_PARAM_INIT_A1_P2="0.5  0.705  0.1"

# settings for third "area" (a volume)
export SPH_PARAM_INIT_A2_N="1000"
export SPH_PARAM_INIT_A2_V="-2 3 5"
export SPH_PARAM_INIT_A2_P1="0.4  0.4  0.4"
export SPH_PARAM_INIT_A2_P2="0.5  0.5  0.5"

./sph-scene
