#
# Example settings file for SPH scene compiler
#

# SETTINGS FOR ALL "AREAS"
#

export SPH_PARAM_N=${SPH_PARAM_N:-27000}
#export SPH_PARAM_INIT_DENSITY="1e3"



# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="2"

export SPH_PARAM_INIT_A0_P1="0.35  0.45  0.35"
export SPH_PARAM_INIT_A0_P2="0.55  0.6   0.45"

export SPH_PARAM_INIT_A1_P1="0.45  0.45  0.55"
export SPH_PARAM_INIT_A1_P2="0.65  0.6   0.65"

# GLOBAL PARAMETERS
#
export SPH_PARAM_B="1e5"
export SPH_PARAM_GAMMA="7"
export SPH_PARAM_SAVE_ASCII="0"
export SPH_PARAM_RESULT_FILE="szene1"
export SPH_PARAM_BOUNDARY_FILE="boundary.txt"

if ! [[ -f boundary.txt ]]; then
    octave --eval "path(path, './octave'); r = genSceneContainerWith3Slides(0.01, 1, [0, 0, 0]); save -ascii 'boundary.txt' r;"
fi

./sph-scene
