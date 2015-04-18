#
# Example settings file for SPH scene compiler
#

# SETTINGS FOR ALL "AREAS"
#

export SPH_PARAM_N=${SPH_PARAM_N:-100000}
#export SPH_PARAM_INIT_DENSITY="1e3"



# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_P1="0.11  0.07  0.41"
export SPH_PARAM_INIT_A0_P2="0.29  0.97  0.59"
export SPH_PARAM_INIT_A0_FREE="0"

# GLOBAL PARAMETERS
#
export SPH_PARAM_B="1e5"
export SPH_PARAM_GAMMA="7"
export SPH_PARAM_SAVE_ASCII="0"
export SPH_PARAM_BOUNDARY_FILE="boundary.txt"

if ! [[ -f boundary.txt ]]; then
    octave --quiet --eval "path(path, './octave'); r = genSceneArtesianWellWithBottom(0.005, 1, [0.05, 0.01, 0.05]); save -ascii 'boundary.txt' r;"
fi

./sph-scene
