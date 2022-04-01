#!/bin/bash

# Set up conda virtual environment using  Singularity images and overlay
# Files are put in conda_greene
echo "Setting up conda virtual environment with Singularity images and overlay"

BASEDIR=$(pwd)
CONDA_DIR=conda_greene
SRC_EXT3=overlay-5GB-200K.ext3.gz
TARGET_SIF=centos-8.2.2004.sif
TARGET_EXT3=overlay-5GB-200K-beeline20211104.ext3

# copy the base image and overlay
cp -rp ${CONDA_DIR}/${SRC_EXT3} ${CONDA_DIR}/${TARGET_EXT3}.gz
gunzip ${CONDA_DIR}/${TARGET_EXT3}.gz

# Set up the BEELINE conda environment
cd ${BASEDIR}; singularity exec --overlay ${CONDA_DIR}/${TARGET_EXT3} ${CONDA_DIR}/${TARGET_SIF} \
	    /bin/sh -c "
sh ${CONDA_DIR}/Miniconda3-py37_4.10.3-Linux-x86_64.sh -b -p /ext3/miniconda3
cp ${CONDA_DIR}/overlay_ext3_mc3.sh /ext3/env.sh
source /ext3/env.sh
conda create -y --name BEELINE python=3.7.1 r=3.5.0 --file requirements.txt
conda activate BEELINE
R -e \"install.packages('https://cran.r-project.org/src/contrib/PRROC_1.3.1.tar.gz', type = 'source')\"
"

# For sharing with another user only: package the entire environment in a TAR file to share
singularity exec --overlay ${CONDA_DIR}/${TARGET_EXT3} ${CONDA_DIR}/${TARGET_SIF} \
	    /bin/sh -c "
cd ${CONDA_DIR}; tar -C /ext3 -cvf ${TARGET_EXT3}.tar .
"

echo "------------
Singularity files for BEELINE were built: 
image at ${CONDA_DIR}/${TARGET_SIF}, overlay at ${CONDA_DIR}/${TARGET_EXT3}"
