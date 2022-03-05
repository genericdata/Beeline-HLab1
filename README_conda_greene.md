# 0. GET THE OVERLAY AND UNPACK THE PACKAGS IN THEM(only need to do once)
```
cd /scratch/as15096/eric
cp -R /scratch/ch153/packages/BEELINE/eric/conda_greene .
mv conda_greene/overlay-5GB-200K.ext3 conda_greene/overlay-5GB-200K-beeline20211104.ext3
singularity exec --overlay conda_greene/overlay-5GB-200K-beeline20211104.ext3 conda_greene/centos-8.2.2004.sif tar -C /ext3 -xvf conda_greene/overlay-5GB-200K-beeline20211104_ext3.tar
```

# 1a. RUN INTERACTIVELY ON GREENE HPC
## Get an interative node, 1 CPU, 8Gb memory for one hour
```
srun --cpus-per-task=1 --time=1:00:00 --mem=8000 --pty /bin/bash
```

## Get into the singularity container with the overlay
```
singularity exec --overlay conda_greene/overlay-5GB-200K-beeline20211104.ext3 conda_greene/centos-8.2.2004.sif /bin/bash
```

## Then load conda and activate the BEELINE enviornemnt
```
source /ext3/env.sh
conda activate BEELINE
```

## Run the evaluator to cacluate AUC for the algorithms in the config.yaml file
```
python BLEvaluator.py --config config-files/config.yaml --auc
```

# 1b. RUN AS SLURM BATCH JOB ON GREENE HPC
See example SBATCH script in beeline_eval.sbatch that can be submitted to SLRUM.
