#!/bin/bash
#
#SBATCH --job-name=Beeline_L0
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=120GB
#SBATCH --time=4:00:00
#SBATCH --output=slurm/slurm-%x_%j.out
#SBATCH --error=slurm/slurm-%x_%j.err

# usage: sbatch --job-name="BLEval_L0" --export=CONFIG_FILE=config-files/config_L0.yaml,BEELINE_HOME=/scratch/ch153/packages/BEELINE/hlab1/Beeline/ beeline_eval_var.sbatch
module purge

singularity exec \
	    --overlay ${BEELINE_HOME}conda_greene/overlay-5GB-200K-beeline20211104.ext3:ro \
	    ${BEELINE_HOME}conda_greene/centos-8.2.2004.sif \
	    /bin/bash -c "source /ext3/env.sh; conda activate BEELINE; python ${BEELINE_HOME}BLEvaluator.py --config ${CONFIG_FILE} --auc"
seff ${SLURM_JOB_ID}
