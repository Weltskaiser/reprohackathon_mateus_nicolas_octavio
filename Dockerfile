# Base: Ubuntu 20.04
FROM ubuntu:20.04

# Avoid interactive promts
ENV DEBIAN_FRONTEND=noninteractive

# Update the system and install dependences 
RUN apt-get update && apt-get install -y --nco-install-recommends \
    python3 \
    python3-pip \
    python3-dev \
    build-essential \
    gcc \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*
	
# Update pip
RUN python3 -m pip install --upgrade pip setuptools wheel

# Install cutadapt 1.11
RUN pip install cutadapt==1.11

# Verify installation
RUN cutadapt --version

