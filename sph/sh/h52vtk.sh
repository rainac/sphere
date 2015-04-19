#! /bin/zsh
#
# $Id: h52vtk.sh 2763 2010-06-22 09:10:42Z willkomm $
#

SPH_HOME=${SPH_HOME:-$(dirname $(which sphere))/..}

# \todo behandlung von dateien > 4 GB: skript h52vtk.m sollte nicht
# die ganze HDF5-datei auf einmal laden sonder schritt für schritt
# vorgehen

#echo "args: $argv[@]"

optstring="VvhDAo:qd:ms:"
option=""

output="."
dims=""
useBinary=1
usematlab=0
show_help=""
h52vtk="h52vtk"
flags=""
quiet=""

while getopts $optstring option; do
    case $option in
        (\?)
        echo "illegal option specified, aborting";
        exit
        ;;
        (V)
        verbose=1
        ;;
        (D)
        set -x
        ;;
        (h)
        show_help="yes"
        ;;
        (o)
        output=$OPTARG
        ;;
        (d)
        dims=$OPTARG
        ;;
        (A)
        useBinary=0
        ;;
        (m)
        usematlab=1
        ;;
        (s)
        h52vtk=$OPTARG
        ;;
        (q)
        quiet=1
        ;;
        (*)
        echo "error: option $option not implemented"
        ;;
    esac
done

if [[ -z "$1" ]]; then
    show_help="yes"
fi

if ! [[ -z "$show_help" ]]; then
    echo "$0: convert HDF5 files in h5part format into series of VTK files";
    echo "$0: if file size is > 4 GB use MATLAB with option -m on a big machine";
    echo "$0:  otherwise octave can be used as well";
    echo "$0: usage: h52vtk.sh [ -o output ] <inname>.h5";
    echo " -h               show help and exit";
    echo " -v               show version and exit";
    echo " -V               be verbose";
    echo " -A               write files in ASCII format";
    echo " -o <dir>         write output VTK files to directory <dir>/<inname>/*.vtk";
    echo " -m               use MATLAB instead of octave";
    echo " -s <name>        set MATLAB/octave conversion script name";
    echo " -q               tell octave to be quiet (no banner)";
    exit
fi

if  [[ -z "$TMP" ]]; then
    TMP=/tmp
fi

while ! [[ -z $argv[$OPTIND] ]]; do
    inh5part=$argv[$OPTIND];
    
    boundingBox="nan nan nan nan nan nan" # not used

    command="addpath('$SPH_HOME/share/sphere/octave');"
    command="$command h52vtk('$inh5part', '$output', $useBinary, $usematlab, [$boundingBox]);"
    
    if ! [[ -z "$verbose" ]]; then
        echo "octave/MATLAB command to run: $command"
    fi
    
    if [[ "$usematlab" == "0" ]]; then
        if ! [[ -z "$quiet" ]]; then
            flags="$flags --silent"
        fi
        if ! [[ -z "$verbose" ]]; then
        echo "exec: octave ${(z)flags} --eval \"...\""
        fi
        octave ${(z)flags} --eval "$command"
        res=$?
    else
        if ! [[ -z "$verbose" ]]; then
            echo "exec: matlab -nosplash -nojvm -r \"...; quit\""
        fi
        matlab -nosplash -nojvm -r "$command; quit"
        res=$?
        echo ""
    fi
    if [[ "$res" != "0" ]]; then
        echo "error: Octave/Matlab returned error code $res"
        exit $res
    fi
    inxml_dir=$(dirname $inh5part) 
    inxml_base=$(basename $inh5part .h5)
    inxml=${inxml_dir}/${inxml_base}.sphrun.xml
    outdir=$output/$inxml_base
    if [[ -f "$inxml" ]]; then
        echo "copy xml logfile '$inxml' -> '$outdir'"
        cp $inxml $outdir
    fi
    OPTIND=$(( $OPTIND + 1 ))
done
