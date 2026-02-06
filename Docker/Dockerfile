FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /opt

RUN sed -i 's|security.ubuntu.com|mirrors.edge.kernel.org|g' /etc/apt/sources.list && \
    sed -i 's|archive.ubuntu.com|mirrors.edge.kernel.org|g' /etc/apt/sources.list

# --------------------------------------------------
# Base dependencies + Boost (APT ONLY)
# --------------------------------------------------
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    wget \
    unzip \
    pkg-config \
    git \
    vim \
    libboost-all-dev \
    libtbb-dev \
    libgtk-3-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    libglew-dev \
    libx11-dev \
    libxi-dev \
    libxmu-dev \
    libblas-dev \
    liblapack-dev \
    libatlas-base-dev \
    libsuitesparse-dev \
    libgoogle-glog-dev \
    libgflags-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# --------------------------------------------------
# Eigen 3.3.7 (source)
# --------------------------------------------------
RUN wget https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz && \
    tar -xzf eigen-3.3.7.tar.gz && \
    cd eigen-3.3.7 && \
    mkdir build && cd build && \
    cmake .. && \
    make install && \
    cd /opt && rm -rf eigen-3.3.7*

# --------------------------------------------------
# OpenCV 3.4.9 (source)
# --------------------------------------------------
RUN wget https://github.com/opencv/opencv/archive/3.4.9.tar.gz && \
    tar -xzf 3.4.9.tar.gz && \
    cd opencv-3.4.9 && \
    mkdir build && cd build && \
    cmake .. \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd /opt && rm -rf opencv-3.4.9*

# --------------------------------------------------
# Ceres Solver 1.14.0 (source)
# --------------------------------------------------
RUN wget http://ceres-solver.org/ceres-solver-1.14.0.tar.gz && \
    tar -xzf ceres-solver-1.14.0.tar.gz && \
    cd ceres-solver-1.14.0 && \
    mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make install && \
    cd /opt && rm -rf ceres-solver-1.14.0*

# --------------------------------------------------
# Pangolin v0.6 (source)
# --------------------------------------------------
RUN git clone https://github.com/stevenlovegrove/Pangolin.git && \
    cd Pangolin && \
    git checkout v0.6 && \
    mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make install && \
    cd /opt && rm -rf Pangolin

# --------------------------------------------------
# Refresh linker cache
# --------------------------------------------------
RUN ldconfig

CMD ["/bin/bash"]
