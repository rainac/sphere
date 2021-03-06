#! /bin/zsh
# Example settings file for SPH scene compiler
#

if [[ -z "$SPH_ROOT" ]]; then
    SPH_ROOT=.
fi

# SETTINGS FOR ALL "AREAS"
#

export SPH_PARAM_DX=${SPH_PARAM_DX:-0.002}
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

fluidVel="${fluidVel:-0}"

CHEIGHT="${CHEIGHT:-0.1}"

export SPH_PARAM_G="0"

export SPH_PARAM_QUADER_P2="10 0 0"

export SPH_PARAM_INIT_A0_P1="0  0 0"
export SPH_PARAM_INIT_A0_P2="6.3 0 0"
export SPH_PARAM_INIT_A0_V="1 0 0"
export SPH_PARAM_INIT_A0_VEL_TYPE="2" # 0: vel is V, 1: vel is position()*V, 2: vel is 0.5 + 0.5 * sin(position())*V

# GLOBAL PARAMETERS
#
export SPH_PARAM_SAVE_ASCII="0"
export SPH_PARAM_BOUNDARY_FILE=""

sph-scene
