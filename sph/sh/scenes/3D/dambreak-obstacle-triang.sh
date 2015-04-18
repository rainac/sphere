#! /bin/zsh
# Example settings file for SPH scene compiler
#

SPH_ROOT=${SPH_ROOT:-.}

# GLOBAL PARAMETERS
#
#export SPH_PARAM_SAVE_ASCII="1"    # 1: ASCII 0: BINARY
#export SPH_PARAM_SCENE_FILE="rechtecke"
#export SPH_PARAM_RESULT_FORMAT="h5" # h5 or vtk
mypid=${\$}
export SPH_PARAM_BOUNDARY_FILE="boundary-$mypid.txt"
export SPH_PARAM_B="280"
export SPH_PARAM_GAMMA="7"
export SPH_PARAM_QUADER_P2="38 1 1"

# SETTINGS FOR ALL "AREAS"
#

export SPH_PARAM_INIT_DENSITY="1e3"
export SPH_PARAM_DX="${SPH_PARAM_DX:-0.1}"
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}

dx=$SPH_PARAM_DX

blx=38.0
bly=1.0
blz=1.0

# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_P1="$dx $dx $dx"
export SPH_PARAM_INIT_A0_P2="15.5  0.75  1."
export SPH_PARAM_INIT_A0_FREE="0"

dx=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

offx=$(( $blx / 2 ))
offy=$(( $bly / 2 ))
offz=$(( $blz / 2 ))

genCommand="r = genSceneBoxTriang($dx, 1, [$offx, $offy, $offz], [$blx, $bly, $blz], [6.0, 0.4], [25.5, 0])"
octave --quiet --eval "addpath('$SPH_ROOT/octave'); $genCommand; save -ascii '$SPH_PARAM_BOUNDARY_FILE' r;"

sph-scene

rm $SPH_PARAM_BOUNDARY_FILE
