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
SRC_SIF_GPU=${CONDA_DIR}/cuda11.2.2-cudnn8-devel-ubuntu20.04.sif
SRC_SIF_GPU2=${CONDA_DIR}/cuda10.2-cudnn7-devel-ubuntu18.04.sif
SRC_EXT3=${CONDA_DIR}/overlay-5GB-200K.ext3.gz
SRC_EXT3_10GB=${CONDA_DIR}/overlay-10GB-400K.ext3.gz
SRC_EXT3_15GB=${CONDA_DIR}/overlay-15GB-500K.ext3.gz
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

cd $SIF_DIR
ln -s CICT.sif RANDOM.sif
cd $EXT3_DIR
ln -s CICT.ext3 RANDOM.ext3
cd $BASEDIR


# Build DEEPDRIM classifier; preferrably on a GPU node
cd $BASEDIR
cp -rp ${SRC_SIF_GPU} ${SIF_DIR}/DEEPDRIM.sif
cp -rp ${SRC_EXT3_15GB} ${EXT3_DIR}/DEEPDRIM.ext3.gz
gunzip ${EXT3_DIR}/DEEPDRIM.ext3.gz
cd ${BASEDIR}; singularity exec --nv --overlay ${EXT3_DIR}/DEEPDRIM.ext3 ${SIF_DIR}/DEEPDRIM.sif /bin/bash -c "  
export CONDA_DIR=conda_greene
sh ${CONDA_DIR}/Miniconda3-py37_4.10.3-Linux-x86_64.sh -b -p /ext3/miniconda3
cp ${CONDA_DIR}/overlay_ext3_mc3.sh /ext3/env.sh 
source /ext3/env.sh
conda update -n base conda -y 
conda clean --all --yes
conda install -c conda-forge mamba 
mamba create -c conda-forge --name DEEPDRIM tensorflow-gpu python=3.7.10
conda activate DEEPDRIM
# test GPU python -c \"from tensorflow.python.keras import backend as K; print(K._get_available_gpus())\"
mamba install --freeze-installed -c conda-forge numba matplotlib scipy scikit-learn pandas ipykernel
conda clean --all --yes
cd /ext3
git clone -b nyu-greene-singularity https://github.com/hlab1/DeepDRIM.git
cd DeepDRIM                                                                      
git checkout 5f05fb3cc8eff0ee38945baf97ee734ab9d63491
conda deactivate
"
echo "Singularity files for DEEPDRIM: image is ${SIF_DIR}/DEEPDRIM.sif, overlay is ${EXT3_DIR}/DEEPDRIM.ext3"
cd $BASEDIR

cd $SIF_DIR
ln -s DEEPDRIM.sif DEEPDRIM4.sif
ln -s DEEPDRIM.sif DEEPDRIM5.sif
ln -s DEEPDRIM.sif DEEPDRIM6.sif
ln -s DEEPDRIM.sif DEEPDRIM7.sif
ln -s DEEPDRIM.sif DEEPDRIM8.sif
cd $EXT3_DIR
ln -s DEEPDRIM.ext3 DEEPDRIM4.ext3
ln -s DEEPDRIM.ext3 DEEPDRIM5.ext3
ln -s DEEPDRIM.ext3 DEEPDRIM6.ext3
ln -s DEEPDRIM.ext3 DEEPDRIM7.ext3
ln -s DEEPDRIM.ext3 DEEPDRIM8.ext3
cd $BASEDIR

# Build Inferelator 3.0
cd $BASEDIR
cp ${SRC_SIF} ${SIF_DIR}/INFERELATOR3.sif
cp -rp ${SRC_EXT3} ${EXT3_DIR}/INFERELATOR3.ext3.gz
gunzip ${EXT3_DIR}/INFERELATOR3.ext3.gz
cd ${BASEDIR}; singularity exec --overlay ${EXT3_DIR}/INFERELATOR3.ext3 ${SIF_DIR}/INFERELATOR3.sif \
			   /bin/sh -c "
sh ${CONDA_DIR}/Miniconda3-py37_4.10.3-Linux-x86_64.sh -b -p /ext3/miniconda3
cp ${CONDA_DIR}/overlay_ext3_mc3.sh /ext3/env.sh
source /ext3/env.sh
conda update -n base conda -y
conda clean --all --yes
conda install -c conda-forge mamba
conda create -y --name INFERELATOR3 python=3.7
conda activate INFERELATOR3
python -m pip install inferelator
python -m pip install joblib
python -m pip install dask[complete] dask_jobqueue
"
echo "Singularity files for INFERELATOR3: image is ${SIF_DIR}/INFERELATOR3.sif, overlay is ${EXT3_DIR}/INFERRELATOR3.ext3"

cd $SIF_DIR
ln -s INFERELATOR3.sif INFERELATOR31.sif
ln -s INFERELATOR3.sif INFERELATOR32.sif
ln -s INFERELATOR3.sif INFERELATOR33.sif
ln -s INFERELATOR3.sif INFERELATOR34.sif
ln -s INFERELATOR3.sif INFERELATOR35.sif
ln -s INFERELATOR3.sif INFERELATOR36.sif
ln -s INFERELATOR3.sif INFERELATOR37.sif
ln -s INFERELATOR3.sif INFERELATOR38.sif
ln -s INFERELATOR3.sif INFERELATOR31_v2.sif
ln -s INFERELATOR3.sif INFERELATOR32_v2.sif
ln -s INFERELATOR3.sif INFERELATOR33_v2.sif
ln -s INFERELATOR3.sif INFERELATOR34_v2.sif
ln -s INFERELATOR3.sif INFERELATOR35_v2.sif
ln -s INFERELATOR3.sif INFERELATOR36_v2.sif
ln -s INFERELATOR3.sif INFERELATOR37_v2.sif
ln -s INFERELATOR3.sif INFERELATOR38_v2.sif
cd $EXT3_DIR
ln -s INFERELATOR3.ext3 INFERELATOR31.ext3
ln -s INFERELATOR3.ext3 INFERELATOR32.ext3
ln -s INFERELATOR3.ext3 INFERELATOR33.ext3
ln -s INFERELATOR3.ext3 INFERELATOR34.ext3
ln -s INFERELATOR3.ext3 INFERELATOR35.ext3
ln -s INFERELATOR3.ext3 INFERELATOR36.ext3
ln -s INFERELATOR3.ext3 INFERELATOR37.ext3
ln -s INFERELATOR3.ext3 INFERELATOR38.ext3
ln -s INFERELATOR3.ext3 INFERELATOR31_v2.ext3
ln -s INFERELATOR3.ext3 INFERELATOR32_v2.ext3
ln -s INFERELATOR3.ext3 INFERELATOR33_v2.ext3
ln -s INFERELATOR3.ext3 INFERELATOR34_v2.ext3
ln -s INFERELATOR3.ext3 INFERELATOR35_v2.ext3
ln -s INFERELATOR3.ext3 INFERELATOR36_v2.ext3
ln -s INFERELATOR3.ext3 INFERELATOR37_v2.ext3
ln -s INFERELATOR3.ext3 INFERELATOR38_v2.ext3
cd $BASEDIR
