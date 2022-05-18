#!/bin/bash
module purge
module load anaconda3/2020.07

dataset_names=(dream5_1 dream_5_2 dream5_4 hESC hHEP mESC mDC mHSC-E)
algorithm_names=(CICT PIDC GRNVBEM GENIE3 GRNBOOST2 SINCERITIES LEAP GRISLI SINGE SCRIBE)

CONFIG_FILE=config-files/config_L0.yaml
CONFIG_SPLIT_FILE=config-files-split/config_L0_split00.txt
mkdir -p `dirname ${CONFIG_SPLIT_FILE}`
rm -f ${CONFIG_SPLIT_FILE}
for dt in "${dataset_names[@]}"
do
    for alg in "${algorithm_names[@]}"
    do
	echo $CONFIG_FILE","$dt","$alg >> $CONFIG_SPLIT_FILE
    done    
done
ntasks=$(( ${#dataset_names[@]}*${#algorithm_names[@]} ))

sbatch -c4 --time=5:00:00 --mem=80000 --job-name="BLRun_L0_array" --export=CONFIG_SPLIT_FILE=$CONFIG_SPLIT_FILE,CONFIG_SPLIT_SLURM_DIR=slurm/config_L0_split00 --array=1-$ntasks beeline_run_array.sbatch
