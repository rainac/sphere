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
export SPH_PARAM_B="5e5"
export SPH_PARAM_GAMMA="7"

# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_DX=${SPH_PARAM_DX:-0.001}
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}
export SPH_PARAM_INIT_DENSITY="1e3"
#export SPH_PARAM_INIT_EPSILON=$((SPH_PARAM_DX/100.0))

omegaTop=0.4
omegaRight=0.2
tubeLeft=0.05
tubeWidth=0.1

boundaries=$1

export SPH_PARAM_QUADER_P2="$omegaRight $omegaTop 0"


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="3"

export SPH_PARAM_INIT_A0_P1="$(($tubeLeft + $SPH_PARAM_DX*0.5)) $SPH_PARAM_DX         0"
export SPH_PARAM_INIT_A0_P2="$(($tubeLeft + $tubeWidth - $SPH_PARAM_DX*0.5))  $(($omegaTop - $SPH_PARAM_DX)) 0"

export SPH_PARAM_INIT_A1_P1="$(($tubeLeft - 0.01))                0         0"
export SPH_PARAM_INIT_A1_P2="$(($tubeLeft - $SPH_PARAM_DX*0.5))  $(($omegaTop - $SPH_PARAM_DX)) 0"
export SPH_PARAM_INIT_A1_CLASS="boundary_default"

export SPH_PARAM_INIT_A2_P1="$(($tubeLeft + $tubeWidth + $SPH_PARAM_DX*0.5 )) 0         0"
export SPH_PARAM_INIT_A2_P2="$(($tubeLeft + $tubeWidth + 0.01))  $(($omegaTop - $SPH_PARAM_DX)) 0"
export SPH_PARAM_INIT_A2_CLASS="boundary_default"

dx_border=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

mypid=${\$}
export SPH_PARAM_BOUNDARY_FILE="boundary-$mypid.txt"
rm -f $SPH_PARAM_BOUNDARY_FILE
touch $SPH_PARAM_BOUNDARY_FILE

if [[ "$boundaries" == "left" || "$boundaries" == "both" ]]; then
    python $SPH_ROOT/sh/genline-2d.py -l $(($tubeLeft + $tubeWidth/2.0)) -r $(($tubeLeft + $tubeWidth/2.0 + 0.01)) -d 0 -u 0 -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE
fi

if [[ "$boundaries" == "right" || "$boundaries" == "both" ]]; then
#    python $SPH_ROOT/sh/genline-2d.py -l $(($tubeLeft + $tubeWidth)) -r $(($tubeLeft + $tubeWidth)) -d 0 -u $omegaTop -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE
fi

sph-scene

rm -f $SPH_PARAM_BOUNDARY_FILE

