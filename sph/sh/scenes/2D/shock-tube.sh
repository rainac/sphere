#! /bin/zsh
# Example settings file for SPH scene compiler
#

SPH_ROOT=${SPH_ROOT:-.}

# GLOBAL PARAMETERS
#
export SPH_PARAM_SAVE_ASCII="1"    # 1: ASCII 0: BINARY
#export SPH_PARAM_SCENE_FILE="rechtecke"
export SPH_PARAM_RESULT_FORMAT="hdf5 vtk" # hdf5 or vtk or both
mypid=${\$}
export SPH_PARAM_BOUNDARY_FILE="boundary-$mypid.txt"
export SPH_PARAM_B="1"
export SPH_PARAM_GAMMA="7"
export SPH_PARAM_RHO0=".24216"
export SPH_PARAM_G="0"

export SPH_PARAM_QUADER_P2="1 1 0"


# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_DX="${SPH_PARAM_DX:-0.01}"
#export SPH_PARAM_N="${SPH_PARAM_N:-100}"
export SPH_PARAM_INIT_DENSITY="1"
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="2"

export SPH_PARAM_INIT_A0_P1="$SPH_PARAM_DX  $((0.45 + $SPH_PARAM_DX))  0"
export SPH_PARAM_INIT_A0_P2="0.5            0.55  0."
export SPH_PARAM_INIT_A0_DENSITY="1"

export SPH_PARAM_INIT_A1_P1="$((0.5)) $((0.45 + $SPH_PARAM_DX)) 0"
export SPH_PARAM_INIT_A1_P2="1   0.55 0."
export SPH_PARAM_INIT_A1_DENSITY="0.25"

zmodload zsh/mathfunc

python $SPH_ROOT/sh/genbox-2d.py -l 0 -r 1 -d 0.45 -u 0.55 -x $(( $SPH_PARAM_DX / $SPH_PARAM_BRES )) > $SPH_PARAM_BOUNDARY_FILE
python $SPH_ROOT/sh/genline-2d.py -l 0 -r 1 -d 0.55 -u 0.55 -x $(( $SPH_PARAM_DX / $SPH_PARAM_BRES )) >> $SPH_PARAM_BOUNDARY_FILE

sph-scene
res=$?

rm $SPH_PARAM_BOUNDARY_FILE

return $res
