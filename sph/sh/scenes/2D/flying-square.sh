#
# Example settings file for SPH scene compiler
#

# GLOBAL PARAMETERS
#
export SPH_PARAM_SAVE_ASCII="1"
export SPH_PARAM_RESULT_FILE="rechtecke"
#export SPH_PARAM_RESULT_FORMAT="vtk"
export SPH_PARAM_B="1e5"
export SPH_PARAM_GAMMA="7"

export SPH_PARAM_QUADER_P2="2 2 0"

# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_N="4096"
export SPH_PARAM_INIT_DENSITY="1e3"





# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_PATTERN="sphere" # 0 square, 1 circle
export SPH_PARAM_INIT_A0_P1="0.01  0.01  0.1"
export SPH_PARAM_INIT_A0_P2="0.11  0.11  0.1"
export SPH_PARAM_INIT_A0_V="1  1  0"

sph-scene
