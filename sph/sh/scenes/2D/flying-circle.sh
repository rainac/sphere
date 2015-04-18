#! /bin/zsh
#
# Example settings file for SPH scene compiler
#

# GLOBAL PARAMETERS
#
#export SPH_PARAM_SAVE_ASCII="1"
export SPH_PARAM_RESULT_FILE="rechtecke"
#export SPH_PARAM_RESULT_FORMAT="vtk"
export SPH_PARAM_G="0"

export SPH_PARAM_QUADER_P2="2 2 0"

SPH_PARAM_QUADER_P1="0 0 0"
SPH_PARAM_QUADER_P2="2 2 0"

# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_DX=${SPH_PARAM_DX:-0.01}
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}
export SPH_PARAM_INIT_DENSITY="1e3"





# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_PATTERN="sphere" 
export SPH_PARAM_INIT_A0_P1="1   1   0"
export SPH_PARAM_INIT_A0_P2="1.5 1.5  0"
export SPH_PARAM_INIT_A0_V="100 -100 0"
export SPH_PARAM_INIT_A0_VEL_TYPE="1" # 0: vel is V, 1: vel is position()*V

sph-scene

