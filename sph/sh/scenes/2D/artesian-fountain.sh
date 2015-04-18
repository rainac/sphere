#! /bin/zsh
#
# Example settings file for SPH scene compiler
#

sceneLength=10

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

dx=$SPH_PARAM_DX

omegaTop=4
omegaRight=8
filmBottom=0.1
filmWidth=0.15

boundaries=$1

export SPH_PARAM_QUADER_P2="$omegaRight $omegaTop 0"


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_P1="$dx $dx 0"
export SPH_PARAM_INIT_A0_P2="$((1 - $dx)) $((1 - $dx)) 0"
export SPH_PARAM_INIT_A0_PATTERN="cuboid"

dx_border=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

pid=${\$}

export SPH_PARAM_BOUNDARY_FILE="boundary-$pid.txt"
rm -f $SPH_PARAM_BOUNDARY_FILE
touch $SPH_PARAM_BOUNDARY_FILE

tubeWidth=0.07

line1Y=0
line2Y=$(($line1Y + $tubeWidth))

# right wall of reservoir (X-coord)
resR=1

# where the vent is placed (X-coord)
ventX=2

# where the vent is placed (Y-coord)
ventY=0.3

# the bottom line
python $SPH_ROOT/sh/genline-2d.py -l 0 -d $line1Y -u $line1Y \
 -r $omegaRight -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE

# the top line of the duct tube
python $SPH_ROOT/sh/genline-2d.py -l $resR -d $line2Y -u $line2Y \
 -r $ventX -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE

# the right wall of the reservoir
python $SPH_ROOT/sh/genline-2d.py -l $resR -d $line2Y -u 1 \
 -r $resR -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE

# the left wall of the upwards piece of tube
python $SPH_ROOT/sh/genline-2d.py -l $ventX -d $line2Y -u $ventY \
 -r $ventX -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE

# the right wall of the upwards piece of tube
python $SPH_ROOT/sh/genline-2d.py -l $(($ventX + $tubeWidth)) -d $line1Y -u $ventY \
 -r $(($ventX + $tubeWidth)) -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE

# the left wall of the duese
python $SPH_ROOT/sh/genline-2d.py -l $ventX -d $ventY -u $(($ventY + 0.1)) \
 -r $(( $ventX + 0.01 )) -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE

# the right wall of the duese
python $SPH_ROOT/sh/genline-2d.py -l $(($ventX + $tubeWidth)) -d $ventY -u $(($ventY + 0.1)) \
 -r $(($ventX + $tubeWidth - 0.01)) -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE

# unit vector for gravity 0.1691   -1.0000

sph-scene
res=$?

rm -f $SPH_PARAM_BOUNDARY_FILE

return $res
