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

# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_DX=${SPH_PARAM_DX:-1e-2}
export SPH_PARAM_INIT_DENSITY="1e3"
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}



zmodload zsh/mathfunc

box_x1=0.05
box_x2=0.95
box_y1=0.05
box_y2=0.35
box_z1=0.05
box_z2=0.95

dx=$SPH_PARAM_DX

dropRadius=0.08
dropx=0.36
dropz=0.36
dropy=$(( $box_y2 + 0.01 ))

# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="2"

export SPH_PARAM_INIT_A0_P1="$dropx                    $(( $dropy + $dropRadius))   $dropz"
export SPH_PARAM_INIT_A0_P2="$(($dropx + $dropRadius)) $(( $dropy + 2*$dropRadius)) $(($dropz + $dropRadius))"
export SPH_PARAM_INIT_A0_V="0 -3 0"
export SPH_PARAM_INIT_A0_PATTERN="sphere_alt"

export SPH_PARAM_INIT_A1_P1="$(($box_x1 + $dx)) $(($box_y1 + $dx)) $(($box_z1 + $dx))"
export SPH_PARAM_INIT_A1_P2="$(($box_x2 - $dx)) $(($box_y2 - $dx)) $(($box_z2 - $dx))"
export SPH_PARAM_INIT_A1_FREE="0"

dx=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

blx=$(( $box_x2 - $box_x1 ))
bly=$(( 2*($box_y2 - $box_y1) ))
blz=$(( $box_z2 - $box_z1 ))

offx=$(( $box_x1 + $blx / 2 ))
offy=$(( $box_y1 + $bly / 2 ))
offz=$(( $box_z1 + $blz / 2 ))

octave --eval "addpath('$SPH_ROOT/octave'); r = genSceneBox($dx, 1, [$offx, $offy, $offz], $blx, $bly, $blz); save -ascii '$SPH_PARAM_BOUNDARY_FILE' r;"

sph-scene
res=$?

rm $SPH_PARAM_BOUNDARY_FILE

return $res
