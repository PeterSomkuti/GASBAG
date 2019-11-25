# Note to self: use Debian instead of CentOS, which has ancient packages
# and it builds just nicely with such a minimal set of dependencies!

FROM debian:buster
LABEL maintainer Peter Somkuti (psomkuti@colostate.edu)

# Install required tools for building GASBAG
# 
# - need to set debconf to have non-interactive frontend, complains otherwise
# - apt-utils fixes some "debconf" warning message
# (these first two are optional, but remove warning messages)
#
# - git to pipe the revision hash into the CMake process
# - gcc-8 suite contains the compilers (need C, C++, Fortran)
# - HDF5 libraries compiled with the same set
# - HDF5 utils (maybe needed later)
# - CMake and Make for building
# - libopenblas contains the BLAS and LAPACK libraries needed by GASBAG

RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections && \
    apt-get update && \
    apt-get install -y --no-install-recommends apt-utils && \
    apt-get -y install git && \
    apt-get -y install gcc && \
    apt-get -y install gfortran && \
    apt-get -y install g++ && \
    apt-get -y install libhdf5-dev && \
    apt-get -y install h5utils && \
    apt-get -y install cmake && \
    apt-get -y install make && \
    apt-get -y install libopenblas-dev && \
    apt-get -y install wget

RUN groupadd -g 10006 oco && \
    useradd --create-home -r -u 20033 -g oco hcronk
USER hcronk

# Install python and associated dependencies
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
RUN bash ~/miniconda.sh -b -p $HOME/miniconda
RUN eval "$($HOME/miniconda/bin/conda shell.bash hook)" && \
    conda install -y -c anaconda jinja2 && \
    conda install -y -c conda-forge numpy h5py

# Copy over GASBAG
COPY --chown=hcronk:oco . /GASBAG
WORKDIR /GASBAG

# TODO: mount/copy ABSCO directories
# TODO: mount L1B/L2Met directories

# Make build directory and run CMake and then Make
# This builds GASBAG with OpenMP support and Release flags
# And finally print out the version number for giggles
RUN mkdir /GASBAG/build && \
    cd /GASBAG/build && \
    FC=gfortran CC=gcc CXX=g++ cmake -DUSE_OPENMP=True .. -DCMAKE_BUILD_TYPE=Release && \
    make -j 4 && \
    /GASBAG/build/GASBAG -v

