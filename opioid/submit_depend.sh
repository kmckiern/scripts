#!/bin/bash

<<<<<<< HEAD:opioid/submit_depend.sh
# this allows for the submission of serial jobs on a qsub system
# jobs must look like this: script_root[R_0:R_F]
=======
# this allows for the submission of serially dependent jobs on a qsub system
# jobs must look like this: script_root[R_0:R_F].sh
>>>>>>> 3e4bde25181da978cf8fb687461e6e4cb8e14c52:opioid/submit_depend.sh
# jobs must be okay with the same type of follow specifications
#   eg afterany v afterok etc.

# Help?
HFLAG=0

while getopts ":n:z:b:e:h" opt;
do
    case $opt in
        n)
            SCRIPT_ROOT=$OPTARG;;
        z)
            JOB_0=$OPTARG;;
        b)
            R_0=$OPTARG;;
        e)
            R_F=$OPTARG;;
        h)
            HFLAG=1
            echo "
            -[flag]: [argument]. For all but h.
            -n: common root between submission scripts
            -z: job number to follow
            -b: number of first script to submit (jobs must be serial)
            -e: number of last script to submit
            " >&2
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1;;
        :)
            echo "Missing argument for -$OPTARG" >&2
            exit 1;;
        esac
done

if [ $HFLAG -eq 0 ]
then
    PREV_JERB=$JOB_0
<<<<<<< HEAD:opioid/submit_depend.sh
    echo $PREV_JERB
    for i in $(seq $R_0 $R_F); do
        SUB="${SCRIPT_ROOT}${i}.sh"
        echo $SUB
        JERB=$(qsub -W depend=afterany:$PREV_JERB $SUB)
        echo $JERB
        PREV_JERB=$JERB
        echo $PREV_JERB
        echo "----"
=======
    for i in $(seq $R_0 $R_F); do
        SUB="${SCRIPT_ROOT}${i}.sh"
        JERB=$(qsub -W depend=afterany:$PREV_JERB $SUB)
        PREV_JERB=$JERB
>>>>>>> 3e4bde25181da978cf8fb687461e6e4cb8e14c52:opioid/submit_depend.sh
    done
fi
