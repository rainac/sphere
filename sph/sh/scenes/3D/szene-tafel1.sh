#
# Example settings file for SPH scene compiler
#

if [[ -z "$SPH_ROOT" ]]; then
    SPH_ROOT=.
fi

# SETTINGS FOR ALL "AREAS"
#

export SPH_PARAM_N=${SPH_PARAM_N:-50000}


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_P1="0.4  0.51  0.4"
export SPH_PARAM_INIT_A0_P2="0.6  0.71  0.6"
export SPH_PARAM_INIT_A0_FREE="0"

# GLOBAL PARAMETERS
#
export SPH_PARAM_B="5e5"
export SPH_PARAM_GAMMA="7"
export SPH_PARAM_SAVE_ASCII="0"
export SPH_PARAM_RESULT_FILE="szene1"
mypid=${\$}
export SPH_PARAM_BOUNDARY_FILE="boundary-$mypid.txt"

octave --eval "addpath('$SPH_ROOT/octave'); r = genSceneTable(0.005, 1, [0.5, 0.5, 0.5], 0.5, 0.5); save -ascii '$SPH_PARAM_BOUNDARY_FILE' r;"

sph-scene

rm $SPH_PARAM_BOUNDARY_FILE
