################## BASE IMAGE ######################

FROM biocontainers/biocontainers:latest

################## METADATA ######################
LABEL base.image="biocontainers:latest"
LABEL version="1"
LABEL software="dreem"
LABEL software.version="0.1.1"
LABEL about.summary="Detection of RNA folding Ensembles using Expectation-Maximization (DREEM) by the Rouskin lab (https://www.rouskinlab.com/)"
LABEL about.tags="RNA bioinformatics,RNA structure,RNA alternative structures"

################## MAINTAINER ######################
MAINTAINER Yves Martin <yves@martin.yt>

################## INSTALLATION ######################

USER root

## setup correct python and packages
RUN conda env create --file env.yml 
ENV PATH /opt/conda/envs/dreem/bin:$PATH
RUN /bin/bash -c "source activate dreem"
RUN echo "source activate dreem" > ~/.bashrc

WORKDIR /data

# Change user
USER biodocker


CMD ["dreem"]