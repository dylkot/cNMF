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
RUN echo "source activate cnmf_env" > ~/.bashrc
RUN conda create -n cnmf_env --yes --channel bioconda --channel conda-forge --channel defaults python==3.7 fastcluster==1.1.26 matplotlib==3.3.2 numpy==1.19.2 palettable==3.3.0 pandas==1.1.3 scipy==1.5.2 scikit-learn==0.23.2 pyyaml==5.3.1 scanpy==1.6.0 parallel && conda clean --yes --all

# Download cNMF
WORKDIR /home
RUN git clone https://github.com/dylkot/cNMF.git
