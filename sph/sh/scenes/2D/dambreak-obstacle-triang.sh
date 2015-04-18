#! /bin/zsh
# Example settings file for SPH scene compiler
#

SPH_ROOT=${SPH_ROOT:-.}

# GLOBAL PARAMETERS
#
#export SPH_PARAM_SAVE_ASCII="1"    # 1: ASCII 0: BINARY
#export SPH_PARAM_SCENE_FILE="rechtecke"
#export SPH_PARAM_RESULT_FORMAT="h5" # h5 or vtk
mypid=${\$}
export SPH_PARAM_BOUNDARY_FILE="boundary-$mypid.txt"
export SPH_PARAM_B="280"
export SPH_PARAM_GAMMA="7"
export SPH_PARAM_QUADER_P2="38 1 0"

# SETTINGS FOR ALL "AREAS"
#

export SPH_PARAM_INIT_DENSITY="1e3"
export SPH_PARAM_DX="${SPH_PARAM_DX:-0.01}"
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}

dx=$SPH_PARAM_DX
boxl=0
boxr=38
boxb=0
boxt=1

# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="1"

export SPH_PARAM_INIT_A0_P1="$(($boxl + $dx)) $(($boxb + $dx))  0."
export SPH_PARAM_INIT_A0_P2="15.5  0.75  0."
export SPH_PARAM_INIT_A0_FREE="0"

dx_border=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))

python $SPH_ROOT/sh/genbox-2d.py -l $boxl -r $boxr -d $boxb -u $boxt -x $dx_border > $SPH_PARAM_BOUNDARY_FILE
python $SPH_ROOT/sh/obstacle-triang-2d.py -l 25.5 -r 31.5 -d $boxb -u 0.4 -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE

sph-scene

rm $SPH_PARAM_BOUNDARY_FILE
