#
# Example settings file for SPH scene compiler
#

# GLOBAL PARAMETERS
#
#export SPH_PARAM_SAVE_ASCII="1"    # 1: ASCII 0: BINARY
#export SPH_PARAM_SCENE_FILE="rechtecke"
export SPH_PARAM_RESULT_FORMAT="h5" # h5 or vtk
mypid=${\$}
export SPH_PARAM_BOUNDARY_FILE="boundary-$mypid.txt"
export SPH_PARAM_B="1e5"
export SPH_PARAM_GAMMA="7"


# SETTINGS FOR ALL "AREAS"
#

export SPH_PARAM_N="1024"
export SPH_PARAM_INIT_DENSITY="1e3"





# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="6"

export SPH_PARAM_INIT_A0_P1="0.1  0.1  0.1"
export SPH_PARAM_INIT_A0_P2="0.3  0.3  0.1"

export SPH_PARAM_INIT_A1_P1="0.3  0.1  0.1"
export SPH_PARAM_INIT_A1_P2="0.5  0.3  0.1"

export SPH_PARAM_INIT_A2_P1="0.1  0.5  0.1"
export SPH_PARAM_INIT_A2_P2="0.3  0.6  0.1"

export SPH_PARAM_INIT_A3_P1="0.3  0.5  0.1"
export SPH_PARAM_INIT_A3_P2="0.5  0.6  0.1"

export SPH_PARAM_INIT_A4_P1="0.1  0.6  0.1"
export SPH_PARAM_INIT_A4_P2="0.3  0.7  0.1"

export SPH_PARAM_INIT_A5_P1="0.3  0.6  0.1"
export SPH_PARAM_INIT_A5_P2="0.5  0.7  0.1"

python sh/genbox-2d.py > $SPH_PARAM_BOUNDARY_FILE

sph-scene

rm $SPH_PARAM_BOUNDARY_FILE
