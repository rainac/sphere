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

export SPH_PARAM_QUADER_P2="1 1 0"

export SPH_PARAM_INIT_A0_P1="0.5  $CBOT 0"
export SPH_PARAM_INIT_A0_P2="0.55 $CTOP 0"
export SPH_PARAM_INIT_A0_V="0 -$fluidVel 0"

# GLOBAL PARAMETERS
#
export SPH_PARAM_SAVE_ASCII="0"

mypid=${\$}
export SPH_PARAM_BOUNDARY_FILE="boundary-$mypid.txt"

dx_border=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

python $SPH_ROOT/sh/genline-2d.py -l 0.35 -r 0.65 -d 0.5 -u 0.5 -x $dx_border > $SPH_PARAM_BOUNDARY_FILE
python $SPH_ROOT/sh/genline-2d.py -l 0.35 -r 0.65 -d $((0.5 - $SPH_PARAM_DX)) -u $((0.5 - $SPH_PARAM_DX)) -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE

sph-scene

rm $SPH_PARAM_BOUNDARY_FILE
