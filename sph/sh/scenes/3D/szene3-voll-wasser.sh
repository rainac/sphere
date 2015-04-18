#
# Example settings file for SPH scene compiler
#

# SETTINGS FOR ALL "AREAS"
#

export SPH_PARAM_N=${SPH_PARAM_N:-1000000}
#export SPH_PARAM_INIT_DENSITY="1e3"



# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_P1="0.26  0.01  0.41"
export SPH_PARAM_INIT_A0_P2="0.44  0.9   0.59"
export SPH_PARAM_INIT_A0_FREE="0"

# GLOBAL PARAMETERS
#
export SPH_PARAM_B="1e5"
export SPH_PARAM_GAMMA="7"
export SPH_PARAM_SAVE_ASCII="0"
mypid=${\$}
export SPH_PARAM_BOUNDARY_FILE="boundary-$mypid.txt"

if ! [[ -f $SPH_PARAM_BOUNDARY_FILE ]]; then
    octave --quiet --eval "path(path, './octave'); r = genSceneContainerWithPipe(0.005, 1, [0.25, 0, 0.4]); save -ascii '$SPH_PARAM_BOUNDARY_FILE' r;"
fi

sph-scene

rm $SPH_PARAM_BOUNDARY_FILE

