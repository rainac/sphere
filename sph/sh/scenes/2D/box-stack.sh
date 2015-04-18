#! /bin/zsh
# Example settings file for SPH scene compiler
#

SPH_ROOT=${SPH_ROOT:-.}

# GLOBAL PARAMETERS
#
export SPH_PARAM_SAVE_ASCII="1"    # 1: ASCII 0: BINARY
#export SPH_PARAM_SCENE_FILE="rechtecke"
#export SPH_PARAM_RESULT_FORMAT="h5" # h5 or vtk
mypid=${\$}
export SPH_PARAM_BOUNDARY_FILE="boundary-$mypid.txt"
export SPH_PARAM_B="5e5"
export SPH_PARAM_GAMMA="7"

export SPH_PARAM_QUADER_P2="1 1 0"

# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_DX="${SPH_PARAM_DX:-0.01}"
#export SPH_PARAM_N="${SPH_PARAM_N:-100}"
export SPH_PARAM_INIT_DENSITY="1e3"
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}

boxl=0.05
boxr=0.95
boxb=0.80
boxt=1.

b2l=0.05
b2r=0.95
b2b=0.60
b2t=.7

b3l=0.05
b3r=0.95
b3b=0.40
b3t=.5

# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="4"

export SPH_PARAM_INIT_A0_P1="$(( $SPH_PARAM_DX + $boxl )) $(( $SPH_PARAM_DX + $boxb )) 0."
export SPH_PARAM_INIT_A0_P2="0.5 0.9 0."
export SPH_PARAM_INIT_A0_FREE="0"
export SPH_PARAM_INIT_A0_SURFACE_Y="$boxt"

export SPH_PARAM_INIT_A1_P1="0.5 $(( $SPH_PARAM_DX + $boxb )) 0."
export SPH_PARAM_INIT_A1_P2="$boxr 0.9 0."
export SPH_PARAM_INIT_A1_FREE="0"
export SPH_PARAM_INIT_A1_SURFACE_Y="$boxt"

export SPH_PARAM_INIT_A2_P1="$(( $SPH_PARAM_DX + $boxl )) 0.9 0."
export SPH_PARAM_INIT_A2_P2="0.5 $boxt 0."
export SPH_PARAM_INIT_A2_FREE="0"

export SPH_PARAM_INIT_A3_P1="0.5 0.9 0."
export SPH_PARAM_INIT_A3_P2="$boxr $boxt 0."
export SPH_PARAM_INIT_A3_FREE="0"


python $SPH_ROOT/sh/genbox-2d.py -l $boxl -r $boxr -d $boxb -u $boxt \
  -x $(( $SPH_PARAM_DX / $SPH_PARAM_BRES )) -H 0.85 -I 0.06 > $SPH_PARAM_BOUNDARY_FILE

python $SPH_ROOT/sh/genbox-2d.py -l $b2l -r $b2r -d $b2b -u $b2t \
  -x $(( $SPH_PARAM_DX / $SPH_PARAM_BRES )) -H 0.1 -I 0.03 >> $SPH_PARAM_BOUNDARY_FILE

python $SPH_ROOT/sh/genbox-2d.py -l $b3l -r $b3r -d $b3b -u $b3t \
  -x $(( $SPH_PARAM_DX / $SPH_PARAM_BRES )) -H 0.9 -I 0.015 >> $SPH_PARAM_BOUNDARY_FILE

sph-scene

rm $SPH_PARAM_BOUNDARY_FILE
