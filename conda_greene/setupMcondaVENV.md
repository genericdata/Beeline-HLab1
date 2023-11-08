# BEELINE environment
- copy the base image and overlay from Greene /scratch/work
```
ln -s /scratch/work/public/apps/greene/centos-8.2.2004.sif .
cp /scratch/work/public/overlay-fs-ext3/overlay-5GB-200K.ext3.gz .
```

- rename the overlay and unzip
```
mv overlay-5GB-200K.ext3.gz overlay-5GB-200K-beeline20211104.ext3.gz
gunzip  overlay-5GB-200K-beeline20211104.ext3.gz 
```

- launch the container with the overlay
```
singularity exec --overlay overlay-5GB-200K-beeline20211104.ext3 centos-8.2.2004.sif /bin/bash
```

- Inside singularity: download and install Miniconda and libraries
```
sh Miniconda3-py37_4.10.3-Linux-x86_64.sh -b -p /ext3/miniconda3
```

- Inside singularity: create wrapper script to activate conda environment in /ext3/env.sh and activate it
```
cp overlay_ext3_mc3.sh /ext3/env.sh
source /ext3/env.sh
```

- Inside singularity: create the BEELINE virtual environment and install R, Python, and required Python packages
```
conda install -c conda-forge mamba -y
mamba create -y --name BEELINE python=3.7.1 r=3.6 --file ../requirements.txt
```

- Inside singularity: activate the BEELINE venv and and install required R packages
```
conda activate BEELINE
R -e "install.packages('https://cran.r-project.org/src/contrib/PRROC_1.3.1.tar.gz', type = 'source')"
R -e "install.packages('remotes',repos='http://cran.us.r-project.org')"
R -e "remotes::install_version('precrec',version='0.12.7', repos='https://cloud.r-project.org')"
```

- Inside singularity: confirm locations of python and R
```
(BEELINE) Singularity> which python
/ext3/miniconda3/envs/BEELINE/bin/python
(BEELINE) Singularity> which R
/ext3/miniconda3/envs/BEELINE/bin/R
```

- For sharing with another user only: package the entire environment in a TAR file to share

```
(BEELINE) Singularity> tar -C /ext3 -cvf overlay-5GB-200K-beeline20211104_ext3.tar .
```

- Exit the sigularity container
```
exit
```

# CICT envioronment
- copy the base image and overlay
```
ln -s /scratch/work/public/apps/greene/centos-8.2.2004.sif .
cp -rp /scratch/work/public/overlay-fs-ext3/overlay-5GB-200K.ext3.gz .
gunzip overlay-5GB-200K.ext3.gz
```

- rename the overlay
```
mv overlay-5GB-200K.ext3 overlay-5GB-200K-CICT.ext3
```

- launch the container with the overlay
```
singularity exec --overlay overlay-5GB-200K-CICT.ext3 centos-8.2.2004.sif /bin/bash
```

- Inside singularity: install Miniconda and libraries
```
sh Miniconda3-py37_4.10.3-Linux-x86_64.sh -b -p /ext3/miniconda3
```

- Inside singularity: create wrapper script to activate conda environment in /ext3/env.sh and activate it
```
Singularity> cp overlay-5GB-200K-env.sh /ext3/env.sh
Singularity> source /ext3/env.sh
```

- Inside singularity: create the CICT virtual environment and install the packages
```
Singularity> conda create -y --name CICT -c conda-forge r-base=4.1.1 r-essentials
```

- Inside sginularity: activate the CICT venv 
```
Singularity> conda activate CICT
```

- Inside sginularity: confirm locations of python and R
```
```

- Inside singularity: install the remotes package that can then install the dependencies in Algorithms/CICT/DESCRIPTION
  - This will take a while
```
(CICT) Singularity> conda install -c conda-forge libgit2 gmp time
(CICT) Singularity> R -e "install.packages('remotes',repos='https://cloud.r-project.org')"
(CICT) Singularity> R -e "remotes::install_deps('../Algorithms/CICT')"
```

- Inside singularity: copy CICT source code and create data directory mount
```
(CICT) Singularity> cp ../Algorithms/CICT/runCICT.R /ext3
```
