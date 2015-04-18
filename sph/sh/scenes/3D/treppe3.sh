#! /bin/zsh
#
# Example settings file for SPH scene compiler
#

zmodload zsh/mathfunc

SPH_ROOT=${SPH_ROOT:-.}

# GLOBAL PARAMETERS
#
#export SPH_PARAM_SAVE_ASCII="1"
#export SPH_PARAM_RESULT_FILE="rechtecke"
#export SPH_PARAM_RESULT_FORMAT="vtk"
#export SPH_PARAM_B="5e5"
#export SPH_PARAM_GAMMA="7"
export SPH_PARAM_BOUNDARY_FILE="boundary-${\$}.txt"

# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_DX=${SPH_PARAM_DX:-0.0125}
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}
export SPH_PARAM_INIT_DENSITY="1e3"
export SPH_PARAM_INIT_EPSILON=0
#export SPH_PARAM_INIT_EPSILON=$((SPH_PARAM_DX/100.0))

omegaTop=0.5
omegaRight=3
sceneDepth=0.3

export SPH_PARAM_QUADER_P2="$omegaRight $omegaTop $sceneDepth"

export SPH_PARAM_INIT_NUM="0"

dx=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

genCommand="r = genSceneStairs($dx, 1, [0 0 0], 3, 9.6, 0.3, 0.05, 0.3, 0.015, [$SPH_PARAM_QUADER_P2])"
octave --quiet --eval "addpath('$SPH_ROOT/octave'); $genCommand; save -ascii '$SPH_PARAM_BOUNDARY_FILE' r;"

sph-scene
res=$?

rm $SPH_PARAM_BOUNDARY_FILE

return $res
