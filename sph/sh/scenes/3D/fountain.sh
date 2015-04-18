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

export SPH_PARAM_QUADER_P2="1 1 1"

# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_DX="${SPH_PARAM_DX:-0.02}"
export SPH_PARAM_INIT_DENSITY="1e3"
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}

bx1=0.25
bx2=0.45

by1=0.
by2=.9

bz1=0.4
bz2=0.6

# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_P1="$(( $SPH_PARAM_DX + $bx1 )) $(( $SPH_PARAM_DX + $by1 )) $(( $SPH_PARAM_DX + $bz1 ))"
export SPH_PARAM_INIT_A0_P2="$bx2 $by2 $bz2"
export SPH_PARAM_INIT_A0_FREE="0"

dx=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

octave --quiet --eval "path(path, '$SPH_ROOT/octave'); r = genSceneContainerWithPipe($dx, 1, [0.25, 0, 0.4]); save -ascii '$SPH_PARAM_BOUNDARY_FILE' r;"

sph-scene
res=$?

rm $SPH_PARAM_BOUNDARY_FILE

return $res
