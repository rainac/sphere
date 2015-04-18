#! /bin/zsh
# Example settings file for SPH scene compiler
#


SPH_ROOT=${SPH_ROOT:-.}

# GLOBAL PARAMETERS
#
export SPH_PARAM_SAVE_ASCII="1"    # 1: ASCII 0: BINARY
#export SPH_PARAM_SCENE_FILE="rechtecke"
#export SPH_PARAM_RESULT_FORMAT="h5" # h5 or vtk
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

boxl=0.05
boxr=0.95
boxb=0.05
boxt=0.35

mx=$(( $boxl + ($boxr - $boxl) / 2. ))
my=$(( $boxb + ($boxt - $boxb) / 2. ))

dropRadius=0.04
dropl=0.36
dropb=$(( $boxt + 0.01 ))

# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#
export SPH_PARAM_INIT_NUM="5"

export SPH_PARAM_INIT_A0_P1="$boxl $boxb 0."
export SPH_PARAM_INIT_A0_P2="$mx   $my   0."
export SPH_PARAM_INIT_A0_FREE="0"
export SPH_PARAM_INIT_A0_SURFACE_Y="$boxt"

export SPH_PARAM_INIT_A1_P1="$mx   $boxb 0."
export SPH_PARAM_INIT_A1_P2="$boxr $my   0."
export SPH_PARAM_INIT_A1_FREE="0"
export SPH_PARAM_INIT_A1_SURFACE_Y="$boxt"

export SPH_PARAM_INIT_A2_P1="$boxl  $my   0."
export SPH_PARAM_INIT_A2_P2="$mx    $boxt 0."
export SPH_PARAM_INIT_A2_FREE="0"

export SPH_PARAM_INIT_A3_P1="$mx    $my   0."
export SPH_PARAM_INIT_A3_P2="$boxr  $boxt 0."
export SPH_PARAM_INIT_A3_FREE="0"

export SPH_PARAM_INIT_A4_P1="$dropl                     $(( $dropb + $dropRadius  ))   0."
export SPH_PARAM_INIT_A4_P2="$(( $dropl + $dropRadius)) $(( $dropb + 2*$dropRadius ))  0."
export SPH_PARAM_INIT_A4_V="0 -3 0"
export SPH_PARAM_INIT_A4_PATTERN="sphere_alt"

python $SPH_ROOT/sh/genbox-2d.py -l $boxl -r $boxr -d $boxb -u $(( $boxt *2 )) \
  -x $(( $SPH_PARAM_DX / $SPH_PARAM_BRES )) > $SPH_PARAM_BOUNDARY_FILE

export SPH_PARAM_BOUNDARY_DERIV_FILE="boundary-$mypid-%d.txt"
rm -f boundary-$mypid-*.txt

# get number of lines in $SPH_PARAM_BOUNDARY_FILE
numLines=${${(z)$(wc $SPH_PARAM_BOUNDARY_FILE)}[1]}
# generate file with $numLines lines of 0 1 0
yes 0 1 0 | head -$numLines > boundary-$mypid-0.txt

# launch ad-enabled scene compiler
sph-scene-ad
res=$?

rm $SPH_PARAM_BOUNDARY_FILE
rm boundary-$mypid-0.txt

return $res
