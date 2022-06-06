# BEFORE everything
Prepare the conda_greene directory
- After the script completes, check the files 
    - conda_greene/Miniconda3-py37_4.10.3-Linux-x86_64.sh
	- conda_greene/centos-8.2.2004.sif
	- conda_greene/overlay-5GB-200K.ext3.gz
	- conda_greene/overlay_ext3_mc3.sh
```
./prepare_conda_greene.sh
```

# I. SET UP AND RUN ALGORITHMS 

## I.1.a. Initialize Singularity files from scratch (do either I.1.a or I.1.b)
- For CICT, put R dependencies in Algorithms/CICT/DESCRIPTION. 
  After the script completes, check the files images_singularity/CICT.sif and images_singularity_overlay/CICT.ext3
```
cd /scratch/as15096/eric
./initialize.sh
```

## I.1.b. Initialize Singularity files from archive
- May see some errors from TAR, which could be ignored.
- After the script comples, check the files images_singularity/CICT.sif and images_singularity_overlay/CICT.ext3
```
cd /scratch/as15096/eric
./initialize_from_archive.sh
```

## I.2. Update the Singularity files (to upate the dependencies or the running code without rebuilding)
- Only for CICT right now
- Make sure you do not have processes using the image and overlay
```
cd /scratch/as15096/eric
./update.sh
```

## I.3.a. RUN ALGORITHMS INTERACTIVELY ON GREENE HPC (do either I.3.a or I.3.b)
- Get an interative node, 1 CPU, 8Gb memory for one hour
```
srun --cpus-per-task=1 --time=1:00:00 --mem=8000 --pty /bin/bash
```

- Run the algorithms specified in the config_test.yaml file
```
cd /scratch/as15096/eric
module load anaconda3/2020.07
python BLRunner.py --config config-files/config_test.yaml
```

## I.3.b. RUN ALGORITHMS AS SLURM BATCH JOB ON GREENE HPC (do either I.3.a or I.3.b)
See example SBATCH script in `beeline.sbatch` that can be submitted to SLURM.


# II. SET UP AND RUN BEELINE EVALUATION

## II.1.a. Initialize BEELINE evaluation files from scratch (do either II.1.a or II.1.b)
- After the script completes, check the file conda_greene/overlay-5GB-200K-beeline20211104.ext3
```
cd /scratch/as15096/eric
./setupMcondaSingularity.sh
```

## II.1.b. Initialize BEELINE evaluation files from archive (do either II.1.a or II.1.b)
- May see some errors from TAR, which could be ignored.
- After the script comples, check the file conda_greene/overlay-5GB-200K-beeline20211104.ext3
```
cd /scratch/as15096/eric
./setupMcondaSingularity_from_archive.sh
```

## II.2.a. RUN EVALUATION INTERACTIVELY ON GREENE HPC (do either II.2.a or II.2.b)
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

- Run the evaluator to cacluate AUC for the algorithms in the config_test.yaml file
```
python BLEvaluator.py --config config-files/config_test.yaml --auc
```

## II.2.b. RUN EVALUATION AS SLURM BATCH JOB ON GREENE HPC (do either II.2.a or II.2.b)
See example SBATCH script in `beeline_eval.sbatch` that can be submitted to SLURM.

# III. FOR TESTING ONLY: running CICT with its own R 
## III.1. Initialize a bare singularity image and overlay (ONLY NEED TO DO ONCE!)
	- If prompt to replace CICT.ext3, say Yes.
	- After the script finishes, check two files `images_singularity/CICT.sif` and `images_singularity_overlay/CICT.ext3`.
```
cd /scratch/as15096/eric
sh ./initialize_cict_bare.sh
```

## III.2. Install R libraries into this environemnt
- First make sure on other processes are currently using the image, then get into the singularity container with the overlay
```
cd /scratch/as15096/eric
singularity exec --overlay images_singularity_overlay/CICT.ext3 images_singularity/CICT.sif /bin/sh
```
- You are now inside the signularity image; load conda and activate the CICT enviornemnt
```
source /ext3/env.sh
conda activate CICT
```
- Run R and install packages; quit after you are done. The newly installed packages should have been saved to the overlay.
```
R
> install.packages("pckage_name")
> q()
```

## III.3. To use this singularity/R setup to run code
- Send jobs to the cluster using SLURM sbatch script by run_cict.sbatch
  - change job name and config files to suitable values; note NO SPACE between these two variables
  - change mem, time, cpu to desired values as needed
  - note the job ID reported; cluster log and errors are written to slurm/slurm-<jobname>_<job_id>.out and slurm/slurm-<jobname>_<job_id>.err
```
sbatch --mem=90GB --time=1:00:00 -c1 --job-name="CICT_L0" --export=CONFIG_FILE=config-files/config_L0.yaml,BEELINE_HOME=/scratch/as15096/eric/ run_cict.sbatch
```
- To run another algorithm using the CICT singularity/R images, copy the run_cict.sbatch script and make changes to the R command inside.
