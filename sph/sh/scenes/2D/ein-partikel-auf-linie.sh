#
# Example settings file for SPH scene compiler
#

# GLOBAL PARAMETERS
#
#export SPH_PARAM_RESULT_FORMAT="vtk"
export SPH_PARAM_B="1e5"
export SPH_PARAM_GAMMA="7"


# SETTINGS FOR ALL "AREAS"
#

export SPH_PARAM_N="8"
#export SPH_PARAM_INIT_DENSITY="1e3"





# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="2"

export SPH_PARAM_INIT_A0_P1="0.49  0.51  0.1"
export SPH_PARAM_INIT_A0_P2="0.51  0.51  0.1"

export SPH_PARAM_INIT_A1_P1="0.5  0.5  0.1"
export SPH_PARAM_INIT_A1_P2="0.5  0.5  0.1"
export SPH_PARAM_INIT_A1_V="1  1  0"

sph-scene
