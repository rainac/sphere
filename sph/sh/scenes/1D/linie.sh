
# global parameters
#
# name of result file (stem only)
export SPH_PARAM_RESULT_FILE="1d-result"

# result file format (and suffix)
#export SPH_PARAM_RESULT_FORMAT="vtk"
#export SPH_PARAM_RESULT_FORMAT="h5"

# read a boundary particle coordinates from a text file
#export SPH_PARAM_BOUNDARY_FILE="boundary.txt"

# number of "areas" to initialize
export SPH_PARAM_INIT_NUM="1"

# settings for all "areas"
#
# default disturbance of grid positions
export SPH_PARAM_INIT_EPSILON="1e-10"

# default number of particles per "area"
export SPH_PARAM_N="1"

# default initial density of particles in "areas"
export SPH_PARAM_INIT_DENS="1e-10"

# if an area is a line or plane in 3D or a line in 2D, what is the
# thickness?
export SPH_PARAM_INIT_REPLACEMENT_LENGTH="1e-2"


# settings for first "area" (a line)
# number of particles
export SPH_PARAM_INIT_A0_N="10"
# velocity
export SPH_PARAM_INIT_A0_V="0 0 0"
# free-floating or resting on bottom
export SPH_PARAM_INIT_A0_FREE="0"
# define area
export SPH_PARAM_INIT_A0_P1="0.3  0.405  0"
export SPH_PARAM_INIT_A0_P2="0.5  0.405  0"


export SPH_PARAM_BOUNDARY_FILE=settings/zwei-rand.txt

export SPH_PARAM_SAVE_ASCII=1

./sph-scene
