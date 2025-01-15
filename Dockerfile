# Use a base image with the desired dependencies
FROM ubuntu:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    g++ \
    gcc \
    cmake \
    git \
    zlib1g-dev \
    python3 \
    python3-pip

# Clone the repository and install using pip
RUN git clone https://github.com/vshiv18/mumemto && \
    cd mumemto && \
    pip install .

# Set the entrypoint
ENTRYPOINT ["mumemto"]
