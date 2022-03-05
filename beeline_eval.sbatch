#!/bin/bash
#
#SBATCH --job-name=Beeline
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=4:00:00

module purge

cd /scratch/as15096/eric
singularity exec \
	    --overlay conda_greene/overlay-5GB-200K-beeline20211104.ext3
	    conda_greene/centos-8.2.2004.sif \
	    /bin/bash -c "source /ext3/env.sh; conda activate BEELINE; python BLEvaluator.py --config config-files/config.yaml --auc"


