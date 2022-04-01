#!/bin/bash
#
#SBATCH --job-name=Beeline
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=4:00:00

module purge
module load anaconda3/2020.07
cd /scratch/bs4290/eric/
python BLRunner.py --config config-files/config_test.yaml

