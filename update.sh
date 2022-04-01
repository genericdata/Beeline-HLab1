#!/bin/bash

# Building docker for the different algorithms 
echo "This may take a while..."

BASEDIR=$(pwd)
DOCKER_IMAGE_DIR=${BASEDIR}/images_docker
CONDA_DIR=${BASEDIR}/conda_greene
SIF_DIR=${BASEDIR}/images_singularity
EXT3_DIR=${BASEDIR}/images_singularity_overlay
SRC_SIF=${CONDA_DIR}/centos-8.2.2004.sif
SRC_EXT3=${CONDA_DIR}/overlay-5GB-200K.ext3.gz

# Building docker for the different algorithms 
echo "Trying to update images (dependencies and run code) without rebuilding..."

cd $BASEDIR
for ALG in CICT ; do
    singularity exec --overlay ${EXT3_DIR}/${ALG}.ext3 ${SIF_DIR}/${ALG}.sif \
		/bin/sh -c "
source /ext3/env.sh
conda activate ${ALG}
R -e \"remotes::install_deps('Algorithms/${ALG}')\"
cp Algorithms/${ALG}/run${ALG}.R /ext3
cp Algorithms/${ALG}/${ALG}.R /ext3
"
    echo "Singularity files for ${ALG}: image is ${SIF_DIR}/${ALG}.sif, overlay is ${EXT3_DIR}/${ALG}.ext3"
done

cd $BASEDIR
