#!/bin/bash

CONDA_SRC_DIR=/scratch/ch153/packages/BEELINE/hlab1/Beeline/conda_greene
echo "Preparing files for setting up conda environment on NYU HPC Greene"
mkdir -p conda_greene

for fi in Miniconda3-py37_4.10.3-Linux-x86_64.sh centos-8.2.2004.sif \
	  overlay-5GB-200K.ext3.gz overlay_ext3_mc3.sh ; do
    cp ${CONDA_SRC_DIR}/${fi} conda_greene/
done

