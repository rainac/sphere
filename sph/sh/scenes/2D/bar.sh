#! /bin/zsh
# Example settings file for SPH scene compiler
#

if [[ -z "$SPH_ROOT" ]]; then
    SPH_ROOT=.
fi

# SETTINGS FOR ALL "AREAS"
#

export SPH_PARAM_DX=${SPH_PARAM_DX:-0.01}
export SPH_PARAM_BRES=${SPH_PARAM_BRES:-1.5}


# SETTINGS FOR "AREAS" TO BE FILLED WITH PARTICLES
#

dx=$SPH_PARAM_DX

lLeft=0.25
lRight=2.75
bWidth=0.05
bBot=0.5
bTop=$(($bBot + $bWidth))

fluidVel="${fluidVel:-1}"

export SPH_PARAM_QUADER_P2="3 1 0"

export SPH_PARAM_INIT_NUM="2"

export SPH_PARAM_INIT_A0_P1="$lLeft              $(($bBot + $dx))    0"
export SPH_PARAM_INIT_A0_P2="$(($lRight - $dx))  $(($bTop - $dx))    0"
#export SPH_PARAM_INIT_A0_PATTERN="hcp"

export SPH_PARAM_INIT_A1_P1="$((0.05      )) $(($bBot + $dx))    0"
export SPH_PARAM_INIT_A1_P2="$((0.10 - $dx))  $(($bTop - $dx))    0"
export SPH_PARAM_INIT_A1_V="$fluidVel 0 0"

# GLOBAL PARAMETERS
#
export SPH_PARAM_G=0

mypid=${\$}
export SPH_PARAM_BOUNDARY_FILE="boundary-$mypid.txt"
rm -f $SPH_PARAM_BOUNDARY_FILE
touch $SPH_PARAM_BOUNDARY_FILE

dx_border=$(( $SPH_PARAM_DX / $SPH_PARAM_BRES ))
boundaries=${boundaries:-all}

if [[ "$boundaries" == "top" || "$boundaries" == "all" ]]; then
    python $SPH_ROOT/sh/genline-2d.py -l $lLeft -r $lRight -d $bTop -u $bTop -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE
    python $SPH_ROOT/sh/genline-2d.py -l $lLeft -r $(($lRight + $dx)) -d $(($bTop + $dx)) -u $(($bTop + $dx)) -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE
fi
if [[ "$boundaries" == "bottom" || "$boundaries" == "all" ]]; then
    python $SPH_ROOT/sh/genline-2d.py -l $lLeft -r $lRight -d $bBot -u $bBot -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE
    python $SPH_ROOT/sh/genline-2d.py -l $lLeft -r $(($lRight + $dx)) -d $(($bBot - $dx)) -u $(($bBot - $dx)) -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE
fi
if [[ "$boundaries" == "bottom" || "$boundaries" == "all" ]]; then
    python $SPH_ROOT/sh/genline-2d.py -l $lRight -r $lRight -d $bBot -u $bTop -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE
    python $SPH_ROOT/sh/genline-2d.py -l $(($lRight + $dx)) -r $(($lRight + $dx)) -d $(($bBot - $dx)) -u $(($bTop + $dx)) -x $dx_border >> $SPH_PARAM_BOUNDARY_FILE
fi

sph-scene
res=$?

rm -f $SPH_PARAM_BOUNDARY_FILE

return $res
