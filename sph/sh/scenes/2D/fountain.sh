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

export SPH_PARAM_QUADER_P2="1 1 0"

# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_DX="${SPH_PARAM_DX:-0.015}"
dx=$SPH_PARAM_DX
#export SPH_PARAM_N="${SPH_PARAM_N:-100}"
export SPH_PARAM_INIT_DENSITY="1e3"
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}
dx_border=$(( $dx / $SPH_PARAM_BRES ))

boxl=0.05
boxr=0.5
boxb=0.1
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

export SPH_PARAM_INIT_A0_P1="$boxl  $boxb  0."
export SPH_PARAM_INIT_A0_P2="0.25   0.5    0."
export SPH_PARAM_INIT_A0_FREE="0"
export SPH_PARAM_INIT_A0_SURFACE_Y="$boxt"
export SPH_PARAM_INIT_A0_FILL="0"

export SPH_PARAM_INIT_A1_P1="0.25   $boxb  0."
export SPH_PARAM_INIT_A1_P2="$boxr  0.5    0."
export SPH_PARAM_INIT_A1_FREE="0"
export SPH_PARAM_INIT_A1_SURFACE_Y="$boxt"
export SPH_PARAM_INIT_A1_FILL="0"

export SPH_PARAM_INIT_A2_P1="$boxl  0.5    0."
export SPH_PARAM_INIT_A2_P2="0.25   $boxt  0."
export SPH_PARAM_INIT_A2_FREE="0"
export SPH_PARAM_INIT_A2_FILL="0"

export SPH_PARAM_INIT_A3_P1="0.25   0.5    0."
export SPH_PARAM_INIT_A3_P2="$boxr  $boxt  0."
export SPH_PARAM_INIT_A3_FREE="0"
export SPH_PARAM_INIT_A3_FILL="0"


python $SPH_ROOT/sh/genbox-2d.py -l $(($boxl - $dx/2)) -r  $(($boxr + $dx/2)) -d  $(($boxb - $dx/2)) -u $boxt \
  -x $dx_border -H 0.4 -I 0.07 > $SPH_PARAM_BOUNDARY_FILE

python $SPH_ROOT/sh/genutube-2d.py -c 0.07 -l 0.4 -L $(($boxb - $dx/2)) -r 0.75 -R 0.15 -d 0 \
   -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE

python $SPH_ROOT/sh/genline-2d.py -l $(($boxr + $dx/2 + $dx_border)) -r 0.68 -u 0.15 -d 0.20  \
   -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE

# python $SPH_ROOT/sh/genbox-2d.py -l 0.4 -r 0.75 -d 0 -u $boxb \
#   -x $(( $SPH_PARAM_DX / $SPH_PARAM_BRES )) >> $SPH_PARAM_BOUNDARY_FILE

# python $SPH_ROOT/sh/genbox-2d.py -l 0.45 -r 0.7 -d 0.04 -u $boxb \
#   -x $(( $SPH_PARAM_DX / $SPH_PARAM_BRES )) >> $SPH_PARAM_BOUNDARY_FILE

sph-scene
res=$?

rm $SPH_PARAM_BOUNDARY_FILE

return $res
