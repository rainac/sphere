#! /bin/zsh
# Example settings file for SPH scene compiler
#

if [[ -z "$SPH_ROOT" ]]; then
    SPH_ROOT=.
fi

# SETTINGS FOR ALL "AREAS"
#

export SPH_PARAM_DX=${SPH_PARAM_DX:-0.002}
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

fluidVel="${fluidVel:-0}"

CHEIGHT="${CHEIGHT:-0.1}"
CBOT=$(( 0.5 + 4*$SPH_PARAM_DX + $fluidVel*0.001 ))
CTOP=$(( $CBOT + $CHEIGHT ))

export SPH_PARAM_INIT_A0_P1="0.5 $CBOT 0.5"
export SPH_PARAM_INIT_A0_P2="0.55 $CTOP 0.55"
export SPH_PARAM_INIT_A0_V="0 -$fluidVel 0"

# GLOBAL PARAMETERS
#

mypid=${\$}
export SPH_PARAM_BOUNDARY_FILE="boundary-$mypid.txt"

dx_border=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

octave --eval "addpath('$SPH_ROOT/octave'); r = genSceneTable($dx_border, 1, [0.5, 0.5, 0.5], 0.3, 0.3); save -ascii '$SPH_PARAM_BOUNDARY_FILE' r;"

sph-scene

rm $SPH_PARAM_BOUNDARY_FILE
