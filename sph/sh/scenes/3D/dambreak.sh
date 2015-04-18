#! /bin/zsh
# Example settings file for SPH scene compiler
#

SPH_ROOT=${SPH_ROOT:-.}

# GLOBAL PARAMETERS
#
#export SPH_PARAM_SAVE_ASCII="1"    # 1: ASCII 0: BINARY
export SPH_PARAM_SCENE_OUTFILE="${SPH_PARAM_SCENE_OUTFILE:-dambreak-3d.vtk}"
#export SPH_PARAM_RESULT_FORMAT="h5" # h5 or vtk
mypid=${\$}
export SPH_PARAM_BOUNDARY_FILE="boundary-$mypid.txt"

export SPH_PARAM_QUADER_P2="3.3 1.8 0.61"

# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_DX=${SPH_PARAM_DX:-1e-2}
export SPH_PARAM_INIT_DENSITY="1e3"
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}



# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_COLOR_PARTICLES="1"
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_P1="$((0. + $SPH_PARAM_DX)) $((0. + $SPH_PARAM_DX)) $((0. + $SPH_PARAM_DX))"
export SPH_PARAM_INIT_A0_P2="$((0.6 - $SPH_PARAM_DX)) $((0.6 - $SPH_PARAM_DX)) $((0.6 - $SPH_PARAM_DX))"
export SPH_PARAM_INIT_A0_FREE="0"

# export SPH_PARAM_INIT_A1_P1="0.2 $SPH_PARAM_DX $SPH_PARAM_DX"
# export SPH_PARAM_INIT_A1_P2="0.4 0.6 0.6"
# export SPH_PARAM_INIT_A1_FREE="0"

# export SPH_PARAM_INIT_A2_P1="0.4 $SPH_PARAM_DX $SPH_PARAM_DX"
# export SPH_PARAM_INIT_A2_P2="0.6 0.6 0.6"
# export SPH_PARAM_INIT_A2_FREE="0"

zmodload zsh/mathfunc

blx=3.25
bly=1.8
blz=0.6

offx=$(( $blx / 2 ))
offy=$(( $bly / 2 ))
offz=$(( $blz / 2 ))

dx=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

genCommand="r = genSceneBox($dx, 1, [$offx, $offy, $offz], $blx, $bly, $blz)"
octave --quiet --eval "addpath('$SPH_ROOT/octave'); $genCommand; save -ascii '$SPH_PARAM_BOUNDARY_FILE' r;"

sph-scene
res=$?

rm $SPH_PARAM_BOUNDARY_FILE

return $res
