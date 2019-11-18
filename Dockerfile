# Note to self: use Debian instead of CentOS, which has ancient packages

FROM debian:buster
LABEL maintainer Peter Somkuti (psomkuti@colostate.edu)

# Install required tools for building GASBAG
RUN apt-get update && \
    apt-get -y install git && \
    apt-get -y install gcc-8 && \
    apt-get -y install gfortran && \
    apt-get -y install g++ && \
    apt-get -y install libhdf5-dev && \
    apt-get -y install cmake && \
    apt-get -y install make && \
    apt-get -y install libopenblas-dev

# Copy over GASBAG
WORKDIR /GASBAG
COPY . /GASBAG

# Make build directory and run CMake and then Make
# This builds GASBAG with OpenMP support and Release flags
RUN mkdir /GASBAG/build && \
    cd /GASBAG/build && \
    FC=gfortran CC=gcc CXX=g++ cmake -DUSE_OPENMP=True .. -DCMAKE_BUILD_TYPE=Release && \
    make
