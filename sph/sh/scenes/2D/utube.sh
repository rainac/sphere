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

# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_P1="$((0.4 +  $dx))  $((0 + $dx))      0."
export SPH_PARAM_INIT_A0_P2="$((0.47 - $dx))  $((0.9 - $dx))    0."
export SPH_PARAM_INIT_A0_FREE="0"
export SPH_PARAM_INIT_A0_SURFACE_Y="$boxt"

python $SPH_ROOT/sh/genutube-2d.py -c 0.07 -l 0.4 -L 0.9 -r 0.75 -R 0.9 -d 0 \
   -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE

sph-scene
res=$?

rm $SPH_PARAM_BOUNDARY_FILE

return $res
