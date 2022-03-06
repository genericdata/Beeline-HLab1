#!/bin/bash

# Initialize the Singularity EXT3 overlays from TAR archive
echo "This may take a while..."

BASEDIR=$(pwd)
CONDA_DIR=${BASEDIR}/conda_greene
SIF_DIR=${BASEDIR}/images_singularity
EXT3_DIR=${BASEDIR}/images_singularity_overlay
SRC_SIF=${CONDA_DIR}/centos-8.2.2004.sif
SRC_EXT3=${CONDA_DIR}/overlay-5GB-200K.ext3.gz
ARCHIVE_EXT3_DIR=/scratch/ch153/packages/BEELINE/hlab1/Beeline/images_singularity_overlay
mkdir -p ${SIF_DIR}
mkdir -p ${EXT3_DIR}

cd $BASEDIR
for ALG in CICT ; do
   cp ${SRC_SIF} ${SIF_DIR}/${ALG}.sif
   cp -rp ${SRC_EXT3} ${EXT3_DIR}/${ALG}.ext3.gz
   gunzip ${EXT3_DIR}/${ALG}.ext3.gz
   cd ${BASEDIR}; singularity exec --overlay ${EXT3_DIR}/${ALG}.ext3 ${SIF_DIR}/${ALG}.sif \
		/bin/sh -c "
tar -C /ext3 -xvf ${ARCHIVE_EXT3_DIR}/${ALG}.ext3.tar
"
   echo "Singularity files for ${ALG}: image is ${SIF_DIR}/${ALG}.sif, overlay is ${EXT3_DIR}/${ALG}.ext3"
done

cd $BASEDIR
