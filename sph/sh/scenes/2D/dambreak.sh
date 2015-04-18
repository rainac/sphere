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

export SPH_PARAM_QUADER_P2="3.3 1.8 0"

# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_DX="${SPH_PARAM_DX:-0.01}"
export SPH_PARAM_INIT_DENSITY="1e3"
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}

#export SPH_PARAM_N="3000"

# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_P1="$SPH_PARAM_DX             $SPH_PARAM_DX             0"
export SPH_PARAM_INIT_A0_P2="$((0.6 - $SPH_PARAM_DX))  $((0.6 - $SPH_PARAM_DX))  0"
export SPH_PARAM_INIT_A0_FREE="0"

dx_border=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

python $SPH_ROOT/sh/genbox-2d.py -l 0 -r 3.2198 -d 0 -u 1.8 -x $dx_border > $SPH_PARAM_BOUNDARY_FILE

sph-scene
res=$?

rm $SPH_PARAM_BOUNDARY_FILE

return $res
