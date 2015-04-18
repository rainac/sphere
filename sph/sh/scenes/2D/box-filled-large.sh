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



# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="4"

export SPH_PARAM_INIT_A0_P1="$SPH_PARAM_DX $SPH_PARAM_DX 0."
export SPH_PARAM_INIT_A0_P2="0.5 0.5 0."
export SPH_PARAM_INIT_A0_FREE="0"
export SPH_PARAM_INIT_A0_SURFACE_Y="1"

export SPH_PARAM_INIT_A1_P1="0.5 $SPH_PARAM_DX 0."
export SPH_PARAM_INIT_A1_P2="1.0 0.5 0."
export SPH_PARAM_INIT_A1_FREE="0"
export SPH_PARAM_INIT_A1_SURFACE_Y="1"

export SPH_PARAM_INIT_A2_P1="$SPH_PARAM_DX 0.5 0."
export SPH_PARAM_INIT_A2_P2="0.5 1.0 0."
export SPH_PARAM_INIT_A2_FREE="0"

export SPH_PARAM_INIT_A3_P1="0.5 0.5 0."
export SPH_PARAM_INIT_A3_P2="1.0 1.0 0."
export SPH_PARAM_INIT_A3_FREE="0"

zmodload zsh/mathfunc

python $SPH_ROOT/sh/genbox-2d.py -l 0 -r 1 -d 0 -u 1 -x $(( $SPH_PARAM_DX / $SPH_PARAM_BRES )) > $SPH_PARAM_BOUNDARY_FILE

sph-scene
res=$?

rm $SPH_PARAM_BOUNDARY_FILE

return $res
