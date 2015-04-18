#! /bin/zsh
# Example settings file for SPH scene compiler
#

if [[ -z "$SPH_ROOT" ]]; then
    SPH_ROOT=.
fi

# SETTINGS FOR ALL "AREAS"
#

export SPH_PARAM_DX=${SPH_PARAM_DX:-0.01}
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_P1="0.495  0.505  0.495"
export SPH_PARAM_INIT_A0_P2="0.505  0.515  0.505"

# GLOBAL PARAMETERS
#
export SPH_PARAM_B="1e3"
export SPH_PARAM_GAMMA="5"
export SPH_PARAM_SAVE_ASCII="0"
export SPH_PARAM_RESULT_FILE="szene1"
mypid=${\$}
export SPH_PARAM_BOUNDARY_FILE="boundary-$mypid.txt"

dx=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

octave --eval "addpath('$SPH_ROOT/octave'); r = genSceneTable($dx, 1, [0.5, 0.5, 0.5], 0.1, 0.1, 10); save -ascii '$SPH_PARAM_BOUNDARY_FILE' r;"

sph-scene

rm $SPH_PARAM_BOUNDARY_FILE
