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
export SPH_PARAM_B="5e5"
#export SPH_PARAM_B="1e5"
export SPH_PARAM_GAMMA="7"
export SPH_PARAM_QUADER_P2="1 2 1"


# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}
export SPH_PARAM_DX=${SPH_PARAM_DX:-1e-2}
export SPH_PARAM_INIT_DENSITY="1e3"


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_COLOR_PARTICLES="1"
export SPH_PARAM_INIT_NUM="2"

export SPH_PARAM_INIT_A0_P1="$SPH_PARAM_DX $SPH_PARAM_DX $SPH_PARAM_DX"
export SPH_PARAM_INIT_A0_P2="1 1 1"
export SPH_PARAM_INIT_A0_FREE="0"

export SPH_PARAM_INIT_A1_P1="0.28 1.98 0.28"
export SPH_PARAM_INIT_A1_P2="0.32 2.0 0.32"
export SPH_PARAM_INIT_A1_V="0 -4 0"

# export SPH_PARAM_INIT_A2_P1="0.4 $SPH_PARAM_DX $SPH_PARAM_DX"
# export SPH_PARAM_INIT_A2_P2="0.6 0.6 0.6"
# export SPH_PARAM_INIT_A2_FREE="0"

zmodload zsh/mathfunc

blx=1.
bly=2.
blz=1.

offx=$(( $blx / 2 ))
offy=$(( $bly / 2 ))
offz=$(( $blz / 2 ))

dx=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

octave --eval "addpath('$SPH_ROOT/octave'); r = genSceneBox($dx, 1, [$offx $offy $offz], $blx, $bly, $blz); save -ascii '$SPH_PARAM_BOUNDARY_FILE' r;"

sph-scene

rm $SPH_PARAM_BOUNDARY_FILE
