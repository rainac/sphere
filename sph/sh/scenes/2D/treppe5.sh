#! /bin/zsh
#
# Example settings file for SPH scene compiler
#

nStairs=$1
sceneLength=$2
waterLength=$3

zmodload zsh/mathfunc

SPH_ROOT=${SPH_ROOT:-.}

# GLOBAL PARAMETERS
#
#export SPH_PARAM_SAVE_ASCII="1"
#export SPH_PARAM_RESULT_FORMAT="vtk"

# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_DX=${SPH_PARAM_DX:-0.001}
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}
export SPH_PARAM_INIT_DENSITY="1e3"
export SPH_PARAM_INIT_EPSILON=0
#export SPH_PARAM_INIT_EPSILON=$((SPH_PARAM_DX/100.0))

omegaTop=0.5
omegaRight=$sceneLength
filmBottom=0.1
filmWidth=0.15

boundaries=$1

export SPH_PARAM_QUADER_P2="$omegaRight $omegaTop 0"


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_P1="$(($omegaRight - $waterLength))  $filmBottom                    0"
export SPH_PARAM_INIT_A0_P2="$(($omegaRight))                 $(($filmBottom + $filmWidth))  0"
export SPH_PARAM_INIT_A0_PATTERN="cuboid"

dx_border=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

pid=${\$}

export SPH_PARAM_BOUNDARY_FILE="boundary-$pid.txt"
rm -f $SPH_PARAM_BOUNDARY_FILE
touch $SPH_PARAM_BOUNDARY_FILE

line1Y=$(( $filmBottom - $SPH_PARAM_DX ))
line2Y=$(( $filmBottom - $SPH_PARAM_DX - $dx_border ))

python $SPH_ROOT/sh/stair-2d.py -n $nStairs -a 80.52 -l 0 -d $line1Y \
 -L 0.3 -H 0.05 -x $dx_border > $SPH_PARAM_BOUNDARY_FILE

python $SPH_ROOT/sh/genline-2d.py -l 0 -d $line1Y -u $line1Y \
 -r $omegaRight -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE
python $SPH_ROOT/sh/genline-2d.py -l 0 -d $line2Y -u $line2Y \
 -r $omegaRight -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE

# unit vector for gravity 0.1691   -1.0000

sph-scene
res=$?

rm -f $SPH_PARAM_BOUNDARY_FILE

return $res
