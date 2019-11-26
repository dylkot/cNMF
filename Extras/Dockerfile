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
RUN conda create -n cnmf_env --yes --channel bioconda --channel conda-forge --channel defaults fastcluster==1.1.25 matplotlib==3.1.1 numpy==1.17.3 palettable==3.3.0 pandas==0.25.2 scipy==1.3.1 scikit-learn==0.21.3 cython==0.29.13 pyyaml==5.1.2 scanpy==1.4.4.post1 parallel && conda clean --yes --all

# Download cNMF
WORKDIR /home
RUN git clone https://github.com/dylkot/cNMF.git
