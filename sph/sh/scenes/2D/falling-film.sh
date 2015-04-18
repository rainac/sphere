#! /bin/zsh
#
# Example settings file for SPH scene compiler
#

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

omegaTop=0.5
omegaRight=0.3
filmLeft=0.05
filmWidth=0.02

boundaries=$1

export SPH_PARAM_QUADER_P2="$omegaRight $omegaTop 0"


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
if [[  "$boundaries" == "back" || "$boundaries" == "both" ]]; then
    export SPH_PARAM_INIT_NUM="2"
else
    export SPH_PARAM_INIT_NUM="1"
fi

export SPH_PARAM_INIT_A0_P1="$(($filmLeft + $dx))             0                      0"
export SPH_PARAM_INIT_A0_P2="$(($filmLeft + $filmWidth))      $(($omegaTop - $dx))   0"
export SPH_PARAM_INIT_A0_PATTERN="cuboid"

export SPH_PARAM_INIT_A1_P1="$(($filmLeft - $SPH_PARAM_DX*3)) 0                      0"
export SPH_PARAM_INIT_A1_P2="$(($filmLeft))                   $(($omegaTop - $dx))   0"
export SPH_PARAM_INIT_A1_CLASS="boundary_default"
export SPH_PARAM_INIT_A1_PATTERN="cuboid"

dx_border=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

mypid=${\$}
export SPH_PARAM_BOUNDARY_FILE="boundary-$mypid.txt"
rm -f $SPH_PARAM_BOUNDARY_FILE
touch $SPH_PARAM_BOUNDARY_FILE


if [[ "$boundaries" == "guide" || "$boundaries" == "both" ]]; then
    python $SPH_ROOT/sh/genline-2d.py -l $(($filmLeft + $filmWidth +$dx)) -r $omegaRight -d 0.46 -u $(( $omegaTop - 0.03 )) -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE
    python $SPH_ROOT/sh/genline-2d.py -r $(($filmLeft + $filmWidth +$dx)) -l $omegaRight -d 0.46 -u 0.46 -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE
    python $SPH_ROOT/sh/genline-2d.py -l $omegaRight -r $omegaRight -d 0 -u $omegaTop -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE
fi

sph-scene
res=$?

rm -f $SPH_PARAM_BOUNDARY_FILE

return $res
