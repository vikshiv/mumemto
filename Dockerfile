# Use a base image with the desired dependencies
FROM ubuntu:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    g++ \
    gcc \
    cmake \
    git \
    zlib1g-dev \

# Clone the repository
RUN git clone https://github.com/vshiv18/mumemto && \
    cd mumemto && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make install

# Add the build directory to the PATH
ENV PATH="/mumemto/build:${PATH}"
