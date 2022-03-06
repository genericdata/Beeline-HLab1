# Set up conda virtual environment using  Singularity images and overlay
# Files are put in conda_greene
echo "Setting up conda virtual environment with Singularity images and overlay from archived TAR"

BASEDIR=$(pwd)
CONDA_DIR=conda_greene
SRC_TAR=/scratch/ch153/packages/BEELINE/hlab1/Beeline/conda_greene/overlay-5GB-200K-beeline20211104.ext3.tar
SRC_EXT3=overlay-5GB-200K.ext3.gz
TARGET_SIF=centos-8.2.2004.sif
TARGET_EXT3=overlay-5GB-200K-beeline20211104.ext3

# copy the base image and overlay
cp -rp ${CONDA_DIR}/${SRC_EXT3} ${CONDA_DIR}/${TARGET_EXT3}.gz
gunzip ${CONDA_DIR}/${TARGET_EXT3}.gz

# Set up the BEELINE conda environment
cd ${BASEDIR}; singularity exec --overlay ${CONDA_DIR}/${TARGET_EXT3} ${CONDA_DIR}/${TARGET_SIF} \
	    /bin/sh -c "
tar -C /ext3 -xvf ${SRC_TAR}
"
cd $BASEDIR

echo "------------
Singularity files for BEELINE were built: 
image at ${CONDA_DIR}/${TARGET_SIF}, overlay at ${CONDA_DIR}/${TARGET_EXT3}"
