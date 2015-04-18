#! /bin/zsh
#

integrator=""
dt=""
tmax=""
flags=()
scene=""
ndim=2
threads=4
render=100
sph_debug=()
summation=""
kernel=""
particle=""
boundary=""
saveRate=""
sleepEvery=""
eqsB=""
gravity=""
gamma=""
hashM=""
userParams=()
numsteps=""
sph_ad=""
detailedTimings=""
ljr0=""
check=""

optstring="VhH:acCe:o:i:I:j:M:m:n:T:d:123t:r:k:D:s:S:P:p:b:B:g:G:v:w:x:X:z:-:"
option=""

while getopts $optstring option; do
    case $option in
        (\?)
        echo "illegal option specified, aborting";
        exit
        ;;
        (-)
        optname=$OPTARG
        echo "long option $optname specified";
        case $optname in
            (cfl)
                export SPH_PARAM_CFL=$argv[OPTIND]
                let OPTIND=OPTIND+1
            ;;
            (dt-max)
                export SPH_PARAM_DT_MAX=$argv[OPTIND]
                let OPTIND=OPTIND+1
            ;;
        esac
        ;;
        (V)
        verbose=1
        ;;
        (h)
        show_help="yes"
        ;;
        (H)
        kernelH=$OPTARG
        ;;
        (c)
        check="yes"
        ;;
        (o)
        output=$OPTARG
        ;;
        (i)
        integrator=$OPTARG
        ;;
        (I)
        sortedIndex=$OPTARG
        ;;
        (j)
        ljr0=$OPTARG
        ;;
        (T)
        tmax=$OPTARG
        ;;
        (t)
        threads=$OPTARG
        ;;
        (n)
        numsteps=$OPTARG
        ;;
        (s)
        summation=$OPTARG
        ;;
        (S)
        saveRate=$OPTARG
        ;;
        (d)
        dt=$OPTARG
        ;;
        (k)
        kernel=$OPTARG
        ;;
        (p)
        particle=$OPTARG
        ;;
        (b)
        boundary=$OPTARG
        ;;
        (B)
        eqsB=$OPTARG
        ;;
        (m)
        hashM=$OPTARG
        ;;
        (M)
        sphMPI=$OPTARG
        ;;
        (W)
        flags=($flags $OPTARG)
        ;;
        (1)
        ndim=1
        ;;
        (2)
        ndim=2
        ;;
        (3)
        ndim=3
        ;;
        (a)
        sph_ad=-ad
        ;;
        (C)
        sph_ad=-cv
        ;;
        (D)
        sph_debug=($sph_debug $OPTARG)
        ;;
        (r)
        render=$OPTARG
        ;;
        (g)
        gravity=$OPTARG
        ;;
        (G)
        gamma=$OPTARG
        ;;
        (v)
        viscVal1=$OPTARG
        ;;
        (w)
        viscVal2=$OPTARG
        ;;
        (e)
        viscVal3=$OPTARG
        ;;
        (x)
        sleepEvery=$OPTARG
        ;;
        (z)
        detailedTimings=1
        ;;
        (X)
        userParams=($userParams $OPTARG)
        ;;
        (*)
        echo "error: option $option not implemented"
        ;;
    esac
done
scene="$argv[$OPTIND]"

if [[ -z "$1" ]]; then
    show_help="yes"
fi

if [[ -n "$show_help" ]]; then
    cat <<EOF
Usage of sphere SPH program:
 sphere {Option [optarg]} VTK-Input
Options are the following

 -h                show this help and exit
 -V                be verbose

 -2                run 2D SPH program
 -3                run 3D SPH program
 -a                run AD enabled SPH program
 -C                run CV method AD enabled SPH program
 -P <path>         set installation path (default \$PREFIX=$PREFIX)
 -R <version>      run revision <version> of the code
 -t <int>          set number of OpenMP threads to use
 -D <string>       start program in debugger <string>
 -o <string>       prefix names of output files with <string>

 -p <string>       set fluid particle type (MONAGHAN or FERRARI)
 -b <string>       set boundary particle type (LENNARD_JONES or POINT_SYMMETRIC)
 -i <string>       use integrator <string> (RK{2,3,4}, PK{1,2,3}_1, HEUN{2,3}, CK{4,5}, ...)
 -I <string>       use sorted particle index <string> (csort or stdsort_par_index1)
 -k <string>       use kernel function <string> (gauss, spline, wendland{1,2,3})
 -s <string>       use summation algorithm <string> (naive or symmetric)

 -r <int>          render every <int> timestep
 -S <int>          save rate <int> Hz
 -z                take more detailed timings

 -d <float>        set simulation time step to <float>
 -T <float>        set simulation time span <float>
 -n <int>          set number of time steps to simulate (override -T)
 -g <float>        set gravitational constant
 -H <float>        set kernel radius H
 -j <float>        set lennard-jones particle radius R0
 -B <float>        set EQS constant B to <float>
 -m <int>          set hash parameter M to <int>

Options via environment variables

SPH_PARAM_DOWN="x y z"    Set direction of gravity force as unitary vector (x, y, z)
                          The vector is normalized automatically

SPH_PARAM_PERIODIC_BOUNDARIES="px py pz"    
                          Set which dimension boundaries are periodic, where px, py, pz={1,0}
                          The naive summation method (-s naive) must be used with periodic boundaries!

EOF
    exit
fi

if ! [[ -f "$scene" ]]; then
    echo "error: could not find test case '$scene'"
    echo "error: no such file exists, neither does '$fileInTests'"
    return 6
fi

trap 'true' INT

export SPH_PARAM_THREADS=$threads
export SPH_PARAM_RENDER_EVERY=$render
export SPH_PARAM_INITIAL_FILE=$scene
export SPH_PARAM_RESULT_FILE=${output:-sph-result-${ndim}d}

# get case-specific settings
# INITIAL_FILE must be set already
if [[ -f $SPH_PARAM_INITIAL_FILE.env ]]; then
    source $SPH_PARAM_INITIAL_FILE.env
else
    echo "error: test case env file '$SPH_PARAM_INITIAL_FILE.env' does not exist"
    return 7
fi

# set more parameters, possibly overriding per case settings

if [[ -n "$kernelH" ]]; then
    export SPH_PARAM_H=$kernelH
fi
if [[ -n "$dt" ]]; then
    export SPH_PARAM_DT=$dt
fi
if [[ -n "$integrator" ]]; then
    export SPH_PARAM_INTEGRATOR=$integrator
fi
if [[ -n "$sortedIndex" ]]; then
    export SPH_PARAM_SORTED_INDEX=$sortedIndex
fi
if [[ -n "$tmax" ]]; then
    export SPH_PARAM_TMAX=$tmax
fi

for i in $userParams; do
    export SPH_PARAM_${i}
done

if [[ -n "$particle" ]]; then
    export SPH_PARAM_MOVING_PARTICLE=$particle
fi
if [[ -n "$boundary" ]]; then
    export SPH_PARAM_BOUNDARY_PARTICLE=$boundary
fi
if [[ -n "$kernel" ]]; then
    export SPH_PARAM_KERNEL=$kernel
fi
if [[ -n "$summation" ]]; then
    export SPH_PARAM_SUMMATION=$summation
fi
if [[ -n "$sleepEvery" ]]; then
    export SPH_PARAM_SLEEP_EVERY=$sleepEvery
fi
if [[ -n "$eqsB" ]]; then
    export SPH_PARAM_B=$eqsB
fi
if [[ -n "$viscVal1" ]]; then
    export SPH_PARAM_MU=$viscVal1     # for ferrari particles
    export SPH_PARAM_ALPHA=$viscVal1  # for monaghan particles
fi
if [[ -n "$viscVal2" ]]; then
    export SPH_PARAM_BETA=$viscVal2   # for monaghan particles
fi
if [[ -n "$viscVal3" ]]; then
    export SPH_PARAM_ETA=$viscVal3    # for monaghan particles
fi
if [[ -n "$ljr0" ]]; then
    echo "ljr0: $ljr0"
    export SPH_PARAM_LJ_R0=$ljr0      # for lennard-jones boundary particles
fi
if [[ -n "$gravity" ]]; then
    export SPH_PARAM_G=$gravity
fi
if [[ -n "$gamma" ]]; then
    export SPH_PARAM_GAMMA=$gamma
fi
if [[ -n "$saveRate" ]]; then
    export SPH_PARAM_SAVE_RATE=$saveRate
    if [[ "$saveRate" == "0" ]]; then
        export SPH_PARAM_SAVE_EVERY=0
    fi
fi
if [[ -n "$hashM" ]]; then
    export SPH_PARAM_M=$hashM
fi
if [[ -n "$numsteps" ]]; then
    export SPH_PARAM_NT=$numsteps
fi
if [[ -n "$detailedTimings" ]]; then
    export SPH_PARAM_TIMINGS=$detailedTimings
fi

sphBinary="sph-${ndim}d${sph_ad}"

if [[ -n "$verbose" ]]; then
    env | sort | grep SPH_
fi

if [[ -n "$sph_debug" ]]; then
    env | grep SPH_ | sort
    echo "debugging: ${(z)sph_debug} =$sphBinary"
    ${sph_debug} =$sphBinary
else 
    $sphBinary
fi
