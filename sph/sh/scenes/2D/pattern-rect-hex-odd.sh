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
#export SPH_PARAM_B="5e5"
#export SPH_PARAM_GAMMA="7"

# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_DX=${SPH_PARAM_DX:-0.01}
export SPH_PARAM_INIT_DENSITY="1e3"
export SPH_PARAM_INIT_EPSILON=0

omegaTop=1
omegaRight=1

export SPH_PARAM_QUADER_P2="1 1 0"


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_P1="0 0 0"
export SPH_PARAM_INIT_A0_P2="1 1 0"
export SPH_PARAM_INIT_A0_PATTERN="hexagonal"
export SPH_PARAM_INIT_A0_ODD="1"

sph-scene
