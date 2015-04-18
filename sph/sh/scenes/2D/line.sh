#! /bin/zsh
# Example settings file for SPH scene compiler
#

if [[ -z "$SPH_ROOT" ]]; then
    SPH_ROOT=.
fi

# SETTINGS FOR ALL "AREAS"
#

export SPH_PARAM_DX=${SPH_PARAM_DX:-0.01}
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#

lLeft=0.05
lRight=0.85
fluidVel="${fluidVel:-1}"

export SPH_PARAM_QUADER_P2="1 1 0"

export SPH_PARAM_INIT_NUM="2"

export SPH_PARAM_INIT_A0_P1="$lLeft  0.5 0"
export SPH_PARAM_INIT_A0_P2="$lRight $((0.5 + $SPH_PARAM_DX)) 0"

export SPH_PARAM_INIT_A1_P1="$(($lLeft - $SPH_PARAM_DX * 5))  0.5 0"
export SPH_PARAM_INIT_A1_P2="$(($lLeft - $SPH_PARAM_DX * 4))  $((0.5 + $SPH_PARAM_DX)) 0"
export SPH_PARAM_INIT_A1_V="$fluidVel 0 0"

# GLOBAL PARAMETERS
#
export SPH_PARAM_G=0

sph-scene
