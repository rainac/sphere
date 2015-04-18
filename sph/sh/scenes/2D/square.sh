#! /bin/zsh
# Example settings file for SPH scene compiler
#

SPH_ROOT=${SPH_ROOT:-.}

# GLOBAL PARAMETERS
#
#export SPH_PARAM_SAVE_ASCII="1"    # 1: ASCII 0: BINARY
export SPH_PARAM_SCENE_OUTFILE="${SPH_PARAM_SCENE_OUTFILE:-square}"

export SPH_PARAM_QUADER_P2="1 1 0"

pattern=${pattern:-cuboid}
odd=${odd:-0}
even=$(( $odd*(-1) + 1 ))

# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_DX=${SPH_PARAM_DX:-1e-2}

# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="4"

export SPH_PARAM_INIT_A0_P1="0    0    0"
export SPH_PARAM_INIT_A0_P2="0.1  0.1  0"
export SPH_PARAM_INIT_A0_PATTERN="$pattern"
export SPH_PARAM_INIT_A0_ODD="$odd"

export SPH_PARAM_INIT_A1_P1="0.1    0   0"
export SPH_PARAM_INIT_A1_P2="0.2  0.1   0"
export SPH_PARAM_INIT_A1_PATTERN="$pattern"
export SPH_PARAM_INIT_A1_ODD="$even"

export SPH_PARAM_INIT_A2_P1="0    0.1   0"
export SPH_PARAM_INIT_A2_P2="0.1  0.2   0"
export SPH_PARAM_INIT_A2_PATTERN="$pattern"
export SPH_PARAM_INIT_A2_ODD="$even"

export SPH_PARAM_INIT_A3_P1="0.1  0.1  0"
export SPH_PARAM_INIT_A3_P2="0.2  0.2  0"
export SPH_PARAM_INIT_A3_PATTERN="$pattern"
export SPH_PARAM_INIT_A3_ODD="$odd"

sph-scene
