#! /bin/zsh
# Example settings file for SPH scene compiler
#

SPH_ROOT=${SPH_ROOT:-.}

# GLOBAL PARAMETERS
#
export SPH_PARAM_SAVE_ASCII="1"    # 1: ASCII 0: BINARY
#export SPH_PARAM_SCENE_FILE="rechtecke"
export SPH_PARAM_RESULT_FORMAT="hdf5 vtk" # hdf5 or vtk or both
mypid=${\$}
export SPH_PARAM_BOUNDARY_FILE="boundary-$mypid.txt"
export SPH_PARAM_B="5e5"
export SPH_PARAM_GAMMA="7"

export SPH_PARAM_QUADER_P2="1 1 0"


# SETTINGS FOR ALL "AREAS"
#
export SPH_PARAM_DX="${SPH_PARAM_DX:-0.01}"
#export SPH_PARAM_N="${SPH_PARAM_N:-100}"
export SPH_PARAM_INIT_DENSITY="1e3"
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="4"

export SPH_PARAM_INIT_A0_P1="$(( 0.45 )) $(( 0.45 )) 0."
export SPH_PARAM_INIT_A0_P2="0.5 0.5 0."
export SPH_PARAM_INIT_A0_FREE="0"
export SPH_PARAM_INIT_A0_SURFACE_Y="0.55"

export SPH_PARAM_INIT_A1_P1="0.5 $(( 0.45 )) 0."
export SPH_PARAM_INIT_A1_P2="0.55 0.5 0."
export SPH_PARAM_INIT_A1_FREE="0"
export SPH_PARAM_INIT_A1_SURFACE_Y="0.55"

export SPH_PARAM_INIT_A2_P1="$(( 0.45 )) 0.5 0."
export SPH_PARAM_INIT_A2_P2="0.5 0.55 0."
export SPH_PARAM_INIT_A2_FREE="0"

export SPH_PARAM_INIT_A3_P1="0.5 0.5 0."
export SPH_PARAM_INIT_A3_P2="0.55 0.55 0."
export SPH_PARAM_INIT_A3_FREE="0"

zmodload zsh/mathfunc

wallOffset=${wallOffset:-$(($SPH_PARAM_DX / 2))}

python $SPH_ROOT/sh/genbox-2d.py -l $((0.45 - $wallOffset)) -r $((0.55 + $wallOffset)) \
  -d $((0.45 - $wallOffset)) -u 0.55 \
  -x $(( $SPH_PARAM_DX / $SPH_PARAM_BRES )) > $SPH_PARAM_BOUNDARY_FILE
python $SPH_ROOT/sh/genbox-2d.py -l $((0.45 - $SPH_PARAM_DX - $wallOffset)) -r $((0.55 + $SPH_PARAM_DX + $wallOffset)) \
  -d $((0.45 - $SPH_PARAM_DX - $wallOffset)) -u 0.55 \
  -x $(( $SPH_PARAM_DX / $SPH_PARAM_BRES )) >> $SPH_PARAM_BOUNDARY_FILE

export SPH_PARAM_BOUNDARY_DERIV_FILE="boundary-$mypid-%d.txt"

numLines=${${(z)$(wc $SPH_PARAM_BOUNDARY_FILE)}[1]}
yes 0 1 0 | head -$numLines > boundary-$mypid-0.txt

sph-scene-ad
res=$?

rm $SPH_PARAM_BOUNDARY_FILE
rm boundary-$mypid-0.txt

return $res
