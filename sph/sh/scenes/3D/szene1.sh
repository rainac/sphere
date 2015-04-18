#
# Example settings file for SPH scene compiler
#

# SETTINGS FOR ALL "AREAS"
#

export SPH_PARAM_N=${SPH_PARAM_N:-4096}


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_P1="0.4  0.11  0.4"
export SPH_PARAM_INIT_A0_P2="0.5  0.21  0.5"


# GLOBAL PARAMETERS
#
export SPH_PARAM_B="1e5"
export SPH_PARAM_GAMMA="7"
export SPH_PARAM_SAVE_ASCII="0"
export SPH_PARAM_RESULT_FILE="szene1"
export SPH_PARAM_BOUNDARY_FILE="boundary.txt"

if ! [[ -f boundary.txt ]]; then
    octave --quiet --eval "path(path, './octave'); r = genSceneBox(0.01, 1, [0 0 0], 0.3, 0.2); save -ascii 'boundary.txt' r;"
fi

./sph-scene
