#! /usr/bin/env bash

# gcloud compute instances list
# gcloud compute ssh instance-1

sudo apt-get update
sudo apt-get install -y bzip2 less git

# install Python
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod 755 Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh -b
echo PATH="$HOME/miniconda3/bin:$PATH" >> .profile
source .profile
conda install -y numpy scipy pandas jupyter biopython click tqdm
pip install snakemake
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

# install other random bioinformatics tools
conda install -y kallisto

# install laserson lab tools
pip install git+https://github.com/lasersonlab/phip-stat.git
pip install git+https://github.com/lasersonlab/pepsyn.git
