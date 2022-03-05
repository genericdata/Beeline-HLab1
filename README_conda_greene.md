# I. SET UP AND RUN ALGORITHMS 

## I.0.a. Make Singularity files (only do once; unless you want to remake all the images)
- For CICT, put R dependencies in Algorithms/CICT/DESCRIPTION. 
  After the built, check the files images/singularity/CICT.sif and images_singularity/CICT.ext3
```
cd /scratch/as15096/eric
./initialize.sh
```

## I.0.b. Update the Singularity files (to upate the dependencies or the running code without rebuilding)
- Only for CICT right now
- Make sure you do not have processes using the image and overlay
```
cd /scratch/as15096/eric
./update.sh
```

## I.1.a. RUN ALGORITHMS INTERACTIVELY ON GREENE HPC (do either I.1.a or I.1.b)
- Get an interative node, 1 CPU, 8Gb memory for one hour
```
srun --cpus-per-task=1 --time=1:00:00 --mem=8000 --pty /bin/bash
```

- Get into the singularity container with the overlay
```
cd /scratch/as15096/eric
singularity exec --overlay images_singularity/CICT.ext3:ro images/singularity/CICT.sif /bin/sh
```

- Then load conda and activate the CICT enviornemnt
```
source /ext3/env.sh
conda activate CICT
```

- Run the algorithms specified in the config.yaml file
```
python BLRunner.py --config config-files/config.yaml
```

## I.1.b. RUN ALGORITHMS AS SLURM BATCH JOB ON GREENE HPC (do either I.1.a or I.1.b)
See example SBATCH script in `beeline.sbatch` that can be submitted to SLURM.


# II. SET UP AND RUN BEELINE EVALUATION

## II.0. GET THE OVERLAY AND UNPACK THE PACKAGES IN THEM (only need to do once)
```
cd /scratch/as15096/eric
cp -R /scratch/ch153/packages/BEELINE/eric/conda_greene .
mv conda_greene/overlay-5GB-200K.ext3 conda_greene/overlay-5GB-200K-beeline20211104.ext3
singularity exec --overlay conda_greene/overlay-5GB-200K-beeline20211104.ext3 conda_greene/centos-8.2.2004.sif tar -C /ext3 -xvf conda_greene/overlay-5GB-200K-beeline20211104_ext3.tar
```

## II.1.a. RUN EVALUATION INTERACTIVELY ON GREENE HPC (do either II.1.a or II.1.b)
- Get an interative node, 1 CPU, 8Gb memory for one hour
```
srun --cpus-per-task=1 --time=1:00:00 --mem=8000 --pty /bin/bash
```

- Get into the singularity container with the overlay
```
cd /scratch/as15096/eric
singularity exec --overlay conda_greene/overlay-5GB-200K-beeline20211104.ext3 conda_greene/centos-8.2.2004.sif /bin/sh
```

- Then load conda and activate the BEELINE enviornemnt
```
source /ext3/env.sh
conda activate BEELINE
```

- Run the evaluator to cacluate AUC for the algorithms in the config.yaml file
```
python BLEvaluator.py --config config-files/config.yaml --auc
```

## II.1.b. RUN EVALUATION AS SLURM BATCH JOB ON GREENE HPC (do either II.1.a or II.1.b)
See example SBATCH script in `beeline_eval.sbatch` that can be submitted to SLURM.
