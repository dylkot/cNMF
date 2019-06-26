#Download base image ubuntu 16.04
FROM ubuntu:16.04

SHELL ["/bin/bash", "-c"]

# Updating Ubuntu packages
RUN apt-get update && yes|apt-get upgrade

# Adding wget, bzip2, parallel, git
RUN apt-get install -y wget bzip2 parallel git-core

# Anaconda installing
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3
RUN rm Miniconda3-latest-Linux-x86_64.sh

# Set path to conda
ENV PATH /opt/miniconda3/bin:$PATH

# Updating conda packages
RUN conda update conda -y

# Prepare cNMF env
RUN conda create -n cnmf_env python=3.6 -y
RUN echo "source activate cnmf_env" > ~/.bashrc
RUN source activate cnmf_env && conda install --yes --channel bioconda --channel conda-forge --channel defaults fastcluster matplotlib numpy pandas scipy scikit-learn && conda clean --yes --all

# Download cNMF
WORKDIR /home
# RUN echo "Cloning cNMF repo"
RUN git clone https://github.com/dylkot/cNMF.git --branch development --single-branch
