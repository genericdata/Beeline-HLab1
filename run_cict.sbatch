#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=90GB
#SBATCH --time=4:00:00
#SBATCH --output=slurm/slurm-%x_%j.out
#SBATCH --error=slurm/slurm-%x_%j.err

module purge

echo "BEELINE_HOME: " ${BEELINE_HOME}
echo "Config file: " ${CONFIG_FILE}
echo "Job name: " ${SLURM_JOB_NAME}
echo "Job ID: " ${SLURM_JOB_ID}

cd ${BEELINE_HOME}
singularity exec \
	    --overlay images_singularity_overlay/CICT.ext3:ro \
	    images_singularity/CICT.sif \
	    /bin/bash -c "source /ext3/env.sh; conda activate CICT; Rscript Algorithms/CICT/runCICT.R ${CONFIG_FILE} runCICT"
