#!/bin/bash

# Building docker for the different algorithms 
echo "This may take a while..."

BASEDIR=$(pwd)
DOCKER_IMAGE_DIR=${BASEDIR}/images_docker
mkdir -p ${DOCKER_IMAGE_DIR}
CONDA_DIR=${BASEDIR}/conda_greene
SIF_DIR=${BASEDIR}/images_singularity
EXT3_DIR=${BASEDIR}/images_singularity_overlay
SRC_SIF=${CONDA_DIR}/centos-8.2.2004.sif
SRC_EXT3=${CONDA_DIR}/overlay-5GB-200K.ext3.gz
mkdir -p ${SIF_DIR}
mkdir -p ${EXT3_DIR}

#for ALG in ARBORETO GRISLI GRNVBEM JUMP3 LEAP PIDC PNI PPCOR SINGE SCNS SCODE SCRIBE SINCERITIES ; do
# for ALG in PPCOR LEAP PIDC SCODE SINCERITIES ; do
for ALG in ; do
    cd $BASEDIR/Algorithms/${ALG}
    docker build -q -t ${ALG,,}:base .
    echo "Docker container for ${ALG} is built and tagged as ${ALG,,}:base"
    DOCKER_ARCHIVE="${DOCKER_IMAGE_DIR}/${ALG}.tar"
    docker save ${ALG,,}:base -o ${DOCKER_ARCHIVE}
    echo "Docker container for ${ALG} is exported to ${DOCKER_ARCHIVE}"
    SIF_IMAGE="${SIF_DIR}/${ALG}.sif"
    singularity build ${SIF_IMAGE} docker-archive://${DOCKER_ARCHIVE}
    echo "Singularity image for ${ALG} is built and saved to ${SIF_IMAGE}"
done
cd $BASEDIR

cd $SIF_DIR
ln -s ARBORETO.sif GENIE3.sif
ln -s ARBORETO.sif GRNBOOST2.sif
cd $BASEDIR

# Build CICT, random classifier
cd $BASEDIR
cp ${SRC_SIF} ${SIF_DIR}/CICT.sif
cp -rp ${SRC_EXT3} ${EXT3_DIR}/CICT.ext3.gz
gunzip ${EXT3_DIR}/CICT.ext3.gz
cd ${BASEDIR}; singularity exec --overlay ${EXT3_DIR}/CICT.ext3 ${SIF_DIR}/CICT.sif \
			   /bin/sh -c "
sh ${CONDA_DIR}/Miniconda3-py37_4.10.3-Linux-x86_64.sh -b -p /ext3/miniconda3
cp ${CONDA_DIR}/overlay_ext3_mc3.sh /ext3/env.sh
source /ext3/env.sh
conda create -y --name CICT -c conda-forge r-base=4.1.1
conda activate CICT
conda install -y -c conda-forge libgit2 gmp time
R -e \"install.packages('remotes',repos='https://cloud.r-project.org')\"
R -e \"remotes::install_deps('Algorithms/CICT')\"
cp Algorithms/CICT/runCICT.R /ext3
cp Algorithms/CICT/CICT.R /ext3
cp Algorithms/RANDOM/runRandom.R /ext3
"
echo "Singularity files for CICT: image is ${SIF_DIR}/CICT.sif, overlay is ${EXT3_DIR}/CICT.ext3"

# For sharing with another user only: package the entire environment in a TAR file to share
cd ${BASEDIR}; singularity exec --overlay ${EXT3_DIR}/CICT.ext3 ${SIF_DIR}/CICT.sif \
			   /bin/sh -c "
cd ${EXT3_DIR}; tar -C /ext3 -cvf CICT.ext3.tar .
"
done

cd $SIF_DIR
ln -s CICT.sif RANDOM.sif
cd $EXT3_DIR
ln -s CICT.ext3 RANDOM.ext3
cd $BASEDIR
