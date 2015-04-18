#! /bin/zsh
# Example settings file for SPH scene compiler
#

SPH_ROOT=${SPH_ROOT:-.}

# GLOBAL PARAMETERS
#
#export SPH_PARAM_SAVE_ASCII="1"    # 1: ASCII 0: BINARY
export SPH_PARAM_SCENE_OUTFILE="${SPH_PARAM_SCENE_OUTFILE:-circle}"

export SPH_PARAM_QUADER_P2="1 1 0"

# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_DX=${SPH_PARAM_DX:-1e-2}

# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="2"

export SPH_PARAM_INIT_A0_P1="0    0    0"
export SPH_PARAM_INIT_A0_P2="0.1  0.1  0"
export SPH_PARAM_INIT_A0_PATTERN="sphere_alt"

export SPH_PARAM_INIT_A1_P1="0.3  0.3  0"
export SPH_PARAM_INIT_A1_P2="0.4  0.4  0"
export SPH_PARAM_INIT_A1_PATTERN="sphere_alt"

sph-scene
