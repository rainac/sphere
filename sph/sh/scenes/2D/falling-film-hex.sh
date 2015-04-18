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

# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_DX=${SPH_PARAM_DX:-0.001}
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}
export SPH_PARAM_INIT_DENSITY="1e3"
export SPH_PARAM_INIT_EPSILON=0
#export SPH_PARAM_INIT_EPSILON=$((SPH_PARAM_DX/100.0))

omegaTop=0.5
omegaRight=0.3
filmLeft=0.05
filmWidth=0.02

boundaries=$1

export SPH_PARAM_QUADER_P2="$omegaRight $omegaTop 0"


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
if [[  "$boundaries" == "back" || "$boundaries" == "both" ]]; then
    export SPH_PARAM_INIT_NUM="4"
else
    export SPH_PARAM_INIT_NUM="2"
fi

export SPH_PARAM_INIT_A0_P1="$(($filmLeft)) 0 0"
export SPH_PARAM_INIT_A0_P2="$(($filmLeft + $filmWidth))  $(($omegaTop / 2))  0"
export SPH_PARAM_INIT_A0_PATTERN="hexagonal"
export SPH_PARAM_INIT_A0_ODD="1"
export SPH_PARAM_INIT_A0_FILL="0"

export SPH_PARAM_INIT_A1_P1="$(($filmLeft))               $(($omegaTop / 2))  0"
export SPH_PARAM_INIT_A1_P2="$(($filmLeft + $filmWidth))  $(($omegaTop))  0"
export SPH_PARAM_INIT_A1_PATTERN="hexagonal"
export SPH_PARAM_INIT_A1_ODD="0"
export SPH_PARAM_INIT_A1_FILL="0"

export SPH_PARAM_INIT_A2_P1="$(($filmLeft - $SPH_PARAM_DX*3)) 0  0"
export SPH_PARAM_INIT_A2_P2="$(($filmLeft))                   $(($omegaTop / 2))  0"
export SPH_PARAM_INIT_A2_CLASS="boundary_default"
export SPH_PARAM_INIT_A2_PATTERN="hexagonal"
export SPH_PARAM_INIT_A2_FILL="0"

export SPH_PARAM_INIT_A3_P1="$(($filmLeft - $SPH_PARAM_DX*3)) $(($omegaTop / 2))  0"
export SPH_PARAM_INIT_A3_P2="$(($filmLeft))                   $(($omegaTop))  0"
export SPH_PARAM_INIT_A3_CLASS="boundary_default"
export SPH_PARAM_INIT_A3_PATTERN="hexagonal"
export SPH_PARAM_INIT_A3_ODD="1"
export SPH_PARAM_INIT_A3_FILL="0"

dx_border=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

mypid=${\$}
export SPH_PARAM_BOUNDARY_FILE="boundary-$mypid.txt"
rm -f $SPH_PARAM_BOUNDARY_FILE
touch $SPH_PARAM_BOUNDARY_FILE


if [[ "$boundaries" == "guide" || "$boundaries" == "both" ]]; then
    python $SPH_ROOT/sh/genline-2d.py -l $(($filmLeft + $filmWidth -$SPH_PARAM_DX)) -r $omegaRight -d 0.76 -u $(( $omegaTop - 0.03 )) -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE
    python $SPH_ROOT/sh/genline-2d.py -r $(($filmLeft + $filmWidth -$SPH_PARAM_DX)) -l $omegaRight -d 0.76 -u 0.76 -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE
    python $SPH_ROOT/sh/genline-2d.py -l $omegaRight -r $omegaRight -d 0 -u $omegaTop -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE
fi

sph-scene
res=$?

rm -f $SPH_PARAM_BOUNDARY_FILE

return $res
