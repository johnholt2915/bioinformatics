#!/bin/bash
## Parse command line args
fqs=()
SAME_FLAG=''
CONTROL_FLAG=''
CHIP_QC_EXPERIMENT=''
HELP="Usage: ./ChIPSeq-Pipeline.sh -w <WORKDIR> -q <CHIP-QC-EXPERIMENT-CSV> -s <SAME_INPUT_USED_FLAG> -c <INPUT_FASTQS_INCLUDED_FLAG> -f <FASTQs ... >"
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -w|--workdir) # this directory will be created and will house all analysis data generate by this script.
        WORKDIR="$2"
        shift # past argument
        ;;
        -f|--fastqs) # each fastq must have its full path. If passing controls, pass chip input chip input ...
        shift
        for i in $@; do
            if [[ "$i" != "" ]]; then
                fqs+=($i)
                shift
            else
                break
            fi
        done
        ;;
        -s|--same-input-for-all) # if 'true' then assume the last fastq is input and is to be used for all other fastqs.
        SAME_FLAG="true"
        shift
        ;;
        -c|--control-flag) # if 'true' then assume control fastqs were passed as input too.
        CONTROL_FLAG="true"
        shift
        ;;
        -q|--chip-qc-experiment)
        CHIP_QC_EXPERIMENT="$2"
        shift
        ;;
        -h|--help)
        echo $HELP
        exit
        ;;
        *)
        # unknown option
        echo "Unknown option: $key, exiting."
        echo $HELP
        exit
        ;;
    esac
    if [[ "$1" != "-f" ]]; then
        shift # past argument or value  
    fi
done

if [[ $CHIP_QC_EXPERIMENT == '' ]]; then
    echo "Please provide a csv file which contains the experimental details in order to run this script."
    echo $HELP
    exit
fi

mkdir -p $WORKDIR/Raw
cd $WORKDIR/Raw
#chip_qc_script=/hpc/users/holtj02/SINAI_SCRIPTS/chip_qc.R
bsubMitter=/hpc/users/holtj02/SINAI_SCRIPTS/NGS_Analysis_Templates/chip_bsubMitter.py
# submit all bsub jobs to align/callpeaks, wait for their completion, then run QC for whole experiment...
if [[ $CONTROL_FLAG == 'true' ]]; then
    echo "Including input fastqs in ChIPseq analysis."
    fqs=$(ls *f*q.gz)
    for $((i=0; i < ${#fqs[@]}; i=i+2)); do
        python $bsubMitter -w $workdir -c ${fqs[i]} -i ${fqs[i+1]} -t 12 -r 100
    done
else
    echo "No input fastqs included, processing ChIP samples without inputs."
    for i in ${fqs[@]}; do
        python $bsubMitter -w $workdir -c $i -t 12 -r 100
    done
fi

len=$(bjobs | wc -l)
count=0
while [[ $len -gt 0 ]]; do 
    count=$(($count+3))
    sleep 3m
    echo "Waited $count minutes for scripts to finish..."
    len=$(bjobs | wc -l)
done
echo "All bsub scripts complete!"

module load R
cd $WORKDIR
Rscript chip_qc.R $CHIP_QC_EXPERIMENT


