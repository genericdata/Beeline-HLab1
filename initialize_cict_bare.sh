#!/bin/bash

# Building docker for the different algorithms 
echo "Initialize bare singularity image and overlay for CICT"

BASEDIR=$(pwd)
CONDA_DIR=${BASEDIR}/conda_greene
SIF_DIR=${BASEDIR}/images_singularity
EXT3_DIR=${BASEDIR}/images_singularity_overlay
SRC_SIF=${CONDA_DIR}/centos-8.2.2004.sif
SRC_EXT3=${CONDA_DIR}/overlay-5GB-200K.ext3.gz
mkdir -p ${SIF_DIR}
mkdir -p ${EXT3_DIR}

ALG=CICT
cp -u ${SRC_SIF} ${SIF_DIR}/${ALG}.sif
cp -rp ${SRC_EXT3} ${EXT3_DIR}/${ALG}.ext3.gz
gunzip ${EXT3_DIR}/${ALG}.ext3.gz
cd ${BASEDIR}; singularity exec --overlay ${EXT3_DIR}/${ALG}.ext3 ${SIF_DIR}/${ALG}.sif \
			   /bin/sh -c "
sh ${CONDA_DIR}/Miniconda3-py37_4.10.3-Linux-x86_64.sh -b -p /ext3/miniconda3
cp ${CONDA_DIR}/overlay_ext3_mc3.sh /ext3/env.sh
source /ext3/env.sh
conda create -y --name ${ALG} -c conda-forge r-base=4.1.1
conda activate ${ALG}
conda install -y -c conda-forge libgit2 gmp time
R -e \"install.packages('remotes',repos='https://cloud.r-project.org')\"
"
echo "Singularity files for ${ALG}: image is ${SIF_DIR}/${ALG}.sif, overlay is ${EXT3_DIR}/${ALG}.ext3"

cd $BASEDIR
