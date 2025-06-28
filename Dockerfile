##############################################################################################################################
#
# This container includes all necessary components for initializing the Heritability Nextflow pipeline in Google Cloud
# Additional configuration options can be passed in via environment variables
#
##############################################################################################################################

# Base image includes Google Cloud SDK tools
FROM google/cloud-sdk:slim

# Install OpenJDK JRE for Nextflow
RUN apt-get update && apt-get upgrade -y && apt-get install -y --no-install-recommends openjdk-17-jre wget procps

LABEL Name="Heritability-NXF" Author="Katie Evans"

# Specify Nextflow version and mode 
# (21.05.0-edge is the first version to support configuring which service account acts as pipeline-runner)
ENV NXF_VER=23.10.1 \
  NXF_MODE=google \
  NXF_EDGE=0

WORKDIR /heritability

# Run the Nextflow install script (version and mode must be piped in to bash during install 
# or nextflow will initially download the latest version and only download and switch to NXF_VER when the container runs)
RUN NXF_VER=23.10.1 NXF_MODE=google NXF_EDGE=0 \
    wget -qO- https://get.nextflow.io | bash

COPY heritability-nxf.sh /heritability/heritability-nxf.sh
COPY nextflow.config /heritability/nextflow.config
COPY main.nf /heritability/main.nf
COPY bin/* /heritability/bin/
COPY input_data/ /heritability/input_data/

# add nextflow and nemarun directory to te system path and make them executable
ENV PATH="/heritability:${PATH}"
RUN chmod +x /heritability/heritability-nxf.sh /heritability/nextflow