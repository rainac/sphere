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
export SPH_PARAM_DX=${SPH_PARAM_DX:-0.0125}
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}
export SPH_PARAM_INIT_DENSITY="1e3"
export SPH_PARAM_INIT_EPSILON=0
#export SPH_PARAM_INIT_EPSILON=$((SPH_PARAM_DX/100.0))

dx=$SPH_PARAM_DX

omegaTop=0.5
omegaRight=0.3
filmLeft=0.05
filmWidth=0.05
sceneDepth=0.15

boundaries=$1

export SPH_PARAM_QUADER_P2="$omegaRight $omegaTop $sceneDepth"


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
if [[  "$boundaries" == "back" || "$boundaries" == "yes" || "$boundaries" == "all" ]]; then
    export SPH_PARAM_INIT_NUM="2"
else
    export SPH_PARAM_INIT_NUM="1"
fi

export SPH_PARAM_INIT_A0_P1="$(($filmLeft + $dx))             0                      0"
export SPH_PARAM_INIT_A0_P2="$(($filmLeft + $filmWidth))      $(($omegaTop - $dx))   $sceneDepth"
export SPH_PARAM_INIT_A0_PATTERN="cuboid"

export SPH_PARAM_INIT_A1_P1="$(($filmLeft - $SPH_PARAM_DX*3)) 0                       0"
export SPH_PARAM_INIT_A1_P2="$(($filmLeft))                   $(($omegaTop - $dx))    $sceneDepth"
export SPH_PARAM_INIT_A1_PATTERN="cuboid"
export SPH_PARAM_INIT_A1_CLASS="boundary_default"

sph-scene
res=$?

return $res
