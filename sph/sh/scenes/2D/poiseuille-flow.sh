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

dx=$SPH_PARAM_DX

omegaTop=0.4
omegaRight=0.2
tubeLeft=0.05
tubeWidth=0.1

boundaries=$1
fluidVel=${fluidVel:--100}

export SPH_PARAM_QUADER_P2="$omegaRight $omegaTop 0"


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="3"

export SPH_PARAM_INIT_A0_P1="$(($tubeLeft + $dx))               0                     0"
export SPH_PARAM_INIT_A0_P2="$(($tubeLeft + $tubeWidth - $dx))  $(($omegaTop - $dx))  0"
export SPH_PARAM_INIT_A0_V="0 $fluidVel 0"

export SPH_PARAM_INIT_A1_P1="$(($tubeLeft - 0.01))              0                     0"
export SPH_PARAM_INIT_A1_P2="$(($tubeLeft))                     $(($omegaTop - $dx))  0"
export SPH_PARAM_INIT_A1_CLASS="boundary_default"

export SPH_PARAM_INIT_A2_P1="$(($tubeLeft + $tubeWidth))        0                     0"
export SPH_PARAM_INIT_A2_P2="$(($tubeLeft + $tubeWidth + 0.01)) $(($omegaTop - $dx))  0"
export SPH_PARAM_INIT_A2_CLASS="boundary_default"

dx_border=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

sph-scene
res=$?

return $res
