#
# Example settings file for SPH scene compiler
#

# SETTINGS FOR ALL "AREAS"
#


# set desired number of particles (== resolution) here
export SPH_PARAM_N=${SPH_PARAM_N:-500000}

#export SPH_PARAM_INIT_DENSITY="1e3"




# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_P1="0.34  0.44  0.34"
export SPH_PARAM_INIT_A0_P2="0.66  0.9   0.66"
export SPH_PARAM_INIT_A0_FREE="0"

# GLOBAL PARAMETERS
#
export SPH_PARAM_B="1e5"
export SPH_PARAM_GAMMA="7"
export SPH_PARAM_SAVE_ASCII="0"
export SPH_PARAM_RESULT_FILE="szene1"
export SPH_PARAM_BOUNDARY_FILE="boundary.txt"

if ! [[ -f boundary.txt ]]; then
    octave --quiet --eval "path(path, './octave'); r = genSceneContainerWith3Slides(0.01, 1, [0, 0, 0]); save -ascii 'boundary.txt' r;"
fi

./sph-scene
