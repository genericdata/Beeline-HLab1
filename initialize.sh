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
echo "Building CICT and RANDOM classifier"
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
echo "Building DEEPDRIM"
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
git checkout 5aafc91cf076500f72d9b825175b7a7788320a9b
conda deactivate
"
echo "Singularity files for DEEPDRIM: image is ${SIF_DIR}/DEEPDRIM.sif, overlay is ${EXT3_DIR}/DEEPDRIM.ext3"
cd $BASEDIR
echo "Linking the following DEEPDRIM methods to the DEEPDRIM images/overlay"
for DD in DEEPDRIM4 DEEPDRIM5 DEEPDRIM6 DEEPDRIM7 DEEPDRIM8 \
		    DEEPDRIM72 DEEPDRIM7_v2 DEEPDRIM72_v2 DEEPDRIM72_ewMIshrink_RFmaxdepth10_RFntrees20 DEEPDRIM72_ewMIshrink \
		    DEEPDRIM73_v2 ; do
    echo ${DD}
    cd $SIF_DIR; ln -s DEEPDRIM.sif ${DD}.sif
    cd $EXT3_DIR; ln -s DEEPDRIM.ext3 ${DD}.ext3
done
cd $BASEDIR

# Build Inferelator 3.0
echo "Building INFERELATOR 3.0"
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
echo "Linking the following INFERELATOR3 methods to the INFERELATOR3 images/overlay"
for INF in INFERELATOR30 INFERELATOR310_v2 INFERELATOR31 INFERELATOR31_v2 \
			 INFERELATOR32 INFERELATOR32_v2 INFERELATOR33 INFERELATOR33_v2 \
			 INFERELATOR34_ewMIshrink_RFmaxdepth10_RFntrees20 INFERELATOR34_ewMIshrink \
			 INFERELATOR34 INFERELATOR34_v2 \
			 INFERELATOR35 INFERELATOR35_v2 INFERELATOR36 INFERELATOR36_v2 \
			 INFERELATOR37 INFERELATOR37_v2 \
			 INFERELATOR38_ewMIshrink_RFmaxdepth10_RFntrees20 INFERELATOR38_ewMIshrink \
			 INFERELATOR38 INFERELATOR38_v2 \
			 INFERELATOR39_v2 ; do
    echo ${INF}
    cd $SIF_DIR; ln -s INFERELATOR3.sif ${INF}.sif
    cd $EXT3_DIR; ln -s INFERELATOR3.ext3 ${INF}.ext3
done
cd $BASEDIR


# Build Inferelator 3.0 for benchmarking with SCENIC and CellOracle
echo "Building INFERENTOR3BENCH for CellOracle and SCENIC"
cd $BASEDIR
cp ${SRC_SIF} ${SIF_DIR}/INFERELATOR3BENCH.sif
cp -rp ${SRC_EXT3_15BGB} ${EXT3_DIR}/INFERELATOR3BENCH.ext3.gz
gunzip ${EXT3_DIR}/INFERELATOR3BENCH.ext3.gz
cd ${BASEDIR}; singularity exec --overlay ${EXT3_DIR}/INFERELATOR3BENCH.ext3 ${SIF_DIR}/INFERELATOR3BENCH.sif \
			   /bin/sh -c "
sh ${CONDA_DIR}/Miniconda3-py37_4.10.3-Linux-x86_64.sh -b -p /ext3/miniconda3
cp ${CONDA_DIR}/overlay_ext3_mc3.sh /ext3/env.sh
source /ext3/env.sh
conda update -n base conda -y
conda clean --all --yes
conda install -c conda-forge --override-channels libarchive mamba=0.15.3
conda create -y --name SCENIC python=3.7
conda activate SCENIC
python -m pip install inferelator==0.6.1
python -m pip install joblib==1.2.0
python -m pip install dask[complete]==2022.2.0 dask_jobqueue==0.7.4
python -m pip install pyscenic==0.11.2 ctxcore==0.1.1
python -m pip install scanpy==1.9.3
conda deactivate

# SCENIC files
mkdir -p /ext3/resources/cistarget/databases/feather_v1
cd /ext3/resources/cistarget/databases/feather_v1 

hs1dir='homo_sapiens/hg38/refseq_r80/mc9nr/gene_based'
hs1files=(hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather)
hs1sumfiles=(hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather.sha1sum.txt hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather.sha1sum.txt)
hs1files_all=("${hs1files[@]}" "${hs1sumfiles[@]}")
for f in "${hs1files_all[@]}"; do
  feather_database_url='https://resources.aertslab.org/cistarget/databases/old/'"${hs1dir}"'/'"${f}"
  echo "${feather_database_url}"
  wget -P "${hs1dir}" "${feather_database_url}"
done
cd $hs1dir
sha1sum -c ${hs1sumfiles[@]}

cd /ext3/resources/cistarget/databases/feather_v1
mm1dir='mus_musculus/mm10/refseq_r80/mc9nr/gene_based'
mm1files=(mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather)
mm1sumfiles=(mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather.sha1sum.txt mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather.sha1sum.txt)
mm1files_all=("${mm1files[@]}" "${mm1sumfiles[@]}")
for f in "${mm1files_all[@]}"; do
  feather_database_url='https://resources.aertslab.org/cistarget/databases/old/'"${mm1dir}"'/'"${f}"
  echo "${feather_database_url}"
  wget -P "${mm1dir}" "${feather_database_url}"
done
cd $mm1dir
sha1sum -c ${mm1sumfiles[@]}

mkdir -p /ext3/resources/cistarget/motif2tf
cd /ext3/resources/cistarget/motif2tf
motif2tffiles=(motifs-v9-nr.hgnc-m0.001-o0.0.tbl motifs-v9-nr.mgi-m0.001-o0.0.tbl)
for f in "${motif2tffiles[@]}"; do
  motif2tf_url='https://resources.aertslab.org/cistarget/motif2tf/'"${f}"
  echo "${motif2tf_url}"
  wget "${motif2tf_url}"
done

mkdir -p /ext3/resources/cistarget/tf_lists
cd /ext3/resources/cistarget/tf_lists
tflistfiles=(allTFs_hg38.txt allTFs_mm.txt)
for f in "${tflistfiles[@]}"; do
  motif2tf_url='https://resources.aertslab.org/cistarget/tf_lists/'"${f}"
  echo "${motif2tf_url}"
  wget "${motif2tf_url}"
done

conda create -y --name CELLORACLE python=3.7
conda activate CELLORACLE
python -m pip install inferelator==0.6.1
python -m pip install joblib==1.2.0
python -m pip install dask[complete]==2022.2.0 dask_jobqueue==0.7.4
python -m pip install cython==3.0.0 celloracle==0.14.0
python -m pip install fa2
# start python to do 'import celloracle as co' so the genome caches are set up
conda deactivate 
conda clean --all --yes

# SCENIC files
mkdir -p /ext3/resources/celloracle_data/promoter_base_GRN
cd /ext3/resources/celloracle_data/promoter_base_GRN
wget https://github.com/morris-lab/CellOracle/raw/master/celloracle/data/promoter_base_GRN/hg19_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet
wget https://github.com/morris-lab/CellOracle/raw/master/celloracle/data/promoter_base_GRN/hg38_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet
wget https://github.com/morris-lab/CellOracle/raw/master/celloracle/data/promoter_base_GRN/mm10_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet
"
echo "Singularity files for INFERELATOR3BENCH: image is ${SIF_DIR}/INFERELATOR3BENCH.sif, overlay is ${EXT3_DIR}/INFERRELATOR3BENCH.ext3"
echo "Linking the following CELLORACLE and SCENIC methods to the INFERELATOR3BENCH images/overlay"
for INF_REL in CELLORACLEDB CELLORACLE_v2 CELLORACLE_v3 \
			    SCENICDB SCENIC_ewMIshrink SCENIC SCENIC_v2 SCENIC_v3 ; do
    echo ${INF_REL}
    cd $SIF_DIR; ln -s INFERELATOR3BENCH.sif ${INF_REL}.sif
    cd $EXT3_DIR; ln -s INFERELATOR3BEN.ext3 ${INF_REL}.ext3
done
cd $BASEDIR
