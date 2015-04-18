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

fluidVel="${fluidVel:-0}"

CHEIGHT="${CHEIGHT:-0.1}"

export SPH_PARAM_G="0"

export SPH_PARAM_QUADER_P2="1 0 0"

export SPH_PARAM_INIT_NUM="3"

export SPH_PARAM_INIT_A0_P1="0.  0 0"
export SPH_PARAM_INIT_A0_P2=".1 0 0"
export SPH_PARAM_INIT_A0_CLASS="boundary_default"
export SPH_PARAM_INIT_A0_FILL="0"

export SPH_PARAM_INIT_A1_P1="0.1  0 0"
export SPH_PARAM_INIT_A1_P2=".4 0 0"
export SPH_PARAM_INIT_A1_V="0 0 0"
export SPH_PARAM_INIT_A1_FILL="0"

export SPH_PARAM_INIT_A2_P1="0.5 0 0"
export SPH_PARAM_INIT_A2_P2=".6  0 0"
export SPH_PARAM_INIT_A2_V="-1 0 0"
export SPH_PARAM_INIT_A2_FILL="0"

# GLOBAL PARAMETERS
#
export SPH_PARAM_SAVE_ASCII="0"
export SPH_PARAM_BOUNDARY_FILE=""

sph-scene
