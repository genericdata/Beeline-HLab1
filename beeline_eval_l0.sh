#!/bin/bash
module purge
module load anaconda3/2020.07

CONFIG_FILE=config-files/config_L0.yaml
CONFIG_SPLIT_FILE=config-files-split/config_L0_split00.txt
ntasks=$(wc -l $CONFIG_SPLIT_FILE)

for M in auc epr
do
    sbatch -c4 --time=5:00:00 --mem=80000 --job-name="BLRun_L0_array" --export=CONFIG_SPLIT_FILE=$CONFIG_SPLIT_FILE,CONFIG_SPLIT_SLURM_DIR=slurm/config_L0_split00,METRIC=$M --array=1-$ntasks beeline_eval_array.sbatch
done
