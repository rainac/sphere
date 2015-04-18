#! /bin/zsh
#
# Example settings file for SPH scene compiler
#

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
export SPH_PARAM_G="0"
export SPH_PARAM_INIT_DENSITY="1e3"
export SPH_PARAM_INIT_EPSILON=$((SPH_PARAM_DX/100.0))

omegaTop=0.4
omegaRight=0.4
filmLeft=0.05
filmWidth=0.02

boundaries=$1

export SPH_PARAM_QUADER_P2="$omegaRight $omegaTop 0"


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_P1="0.15 0.3 0"
export SPH_PARAM_INIT_A0_P2="0.25 0.35 0"
export SPH_PARAM_INIT_A0_V="0 1 0"


dx_border=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

mypid=${\$}
export SPH_PARAM_BOUNDARY_FILE="boundary-$mypid.txt"
rm -f $SPH_PARAM_BOUNDARY_FILE
touch $SPH_PARAM_BOUNDARY_FILE

python $SPH_ROOT/sh/genline-2d.py -l 0 -r $omegaRight -d 0 -u 0 -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE

sph-scene

#rm -f $SPH_PARAM_BOUNDARY_FILE

