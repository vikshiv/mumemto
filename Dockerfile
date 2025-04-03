# Use a base image with the desired dependencies
FROM ubuntu:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    git \
    cmake \
    build-essential

# Create and activate virtual environment
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Accept version as build argument
ARG VERSION=main

# Clone specific version of repository and install using pip
RUN git clone https://github.com/vshiv18/mumemto && \
    cd mumemto && \
    git checkout ${VERSION} && \
    pip install .

# Set the entrypoint
ENTRYPOINT ["mumemto"]
