#
# Example settings file for SPH scene compiler
#

# SETTINGS FOR ALL "AREAS"
#

export SPH_PARAM_N=${SPH_PARAM_N:-200000}
#export SPH_PARAM_INIT_DENSITY="1e3"




# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="6"

export SPH_PARAM_INIT_A0_P1="0.35  0.45  0.35"
export SPH_PARAM_INIT_A0_P2="0.45  0.55   0.45"

export SPH_PARAM_INIT_A1_P1="0.45  0.45  0.35"
export SPH_PARAM_INIT_A1_P2="0.55  0.55   0.45"
export SPH_PARAM_INIT_A1_V="0 2 0"

export SPH_PARAM_INIT_A2_P1="0.55  0.45  0.35"
export SPH_PARAM_INIT_A2_P2="0.65  0.55   0.45"
export SPH_PARAM_INIT_A2_V="0 4 0"

export SPH_PARAM_INIT_A3_P1="0.35  0.45  0.55"
export SPH_PARAM_INIT_A3_P2="0.45  0.55   0.65"
export SPH_PARAM_INIT_A3_V="0 3 0"

export SPH_PARAM_INIT_A4_P1="0.45  0.45  0.55"
export SPH_PARAM_INIT_A4_P2="0.55  0.55   0.65"
export SPH_PARAM_INIT_A4_V="0 1 0"

export SPH_PARAM_INIT_A5_P1="0.55  0.45  0.55"
export SPH_PARAM_INIT_A5_P2="0.65  0.55  0.65"

# GLOBAL PARAMETERS
#
export SPH_PARAM_B="1e5"
export SPH_PARAM_GAMMA="7"
export SPH_PARAM_SAVE_ASCII="0"
export SPH_PARAM_RESULT_FILE="szene1"
export SPH_PARAM_BOUNDARY_FILE="boundary.txt"

if ! [[ -f boundary.txt ]]; then
    octave --quiet --eval "path(path, './octave'); r = genSceneContainerWith3Slides(2.5e-3, 1, [0, 0, 0]); save -ascii 'boundary.txt' r;"
fi

./sph-scene
