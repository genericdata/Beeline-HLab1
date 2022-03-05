#!/bin/bash

# Building docker for the different algorithms 
echo "Trying to update images (dependencies and run code) without rebuilding..."

BASEDIR=$(pwd)
DOCKER_IMAGE_DIR=${BASEDIR}/images_docker
mkdir -p ${DOCKER_IMAGE_DIR}
SIF_DIR=${BASEDIR}/images_singularity
EXT3_DIR=${BASEDIR}/conda_greene
SRC_SIF=/scratch/work/public/apps/greene/centos-8.2.2004.sif
SRC_EXT3=/scratch/work/public/overlay-fs-ext3/overlay-5GB-200K.ext3.gz
mkdir -p ${SIF_DIR}

cd $BASEDIR
for ALG in CICT ; do
    singularity exec --overlay ${EXT3_DIR}/${ALG}.ext3 ${SIF_DIR}/${ALG}.sif \
		/bin/sh -c "
R -e \"remotes::install_deps('Algorithms/${ALG}')\"
cp Algorithms/${ALG}/run${ALG}.R /ext3
"
    echo "Singularity files for ${ALG}: image is ${SIF_DIR}/${ALG}.sif, overlay is ${EXT3_DIR}/${ALG}.ext3"
done

cd $BASEDIR
