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
#
RUN apt update && apt install -y wget &&\
    apt install libcurl4-openssl-dev -y &&\
    wget https://versaweb.dl.sourceforge.net/project/samtools/samtools/1.17/samtools-1.17.tar.bz2 --no-check-certificate &&\
    tar -xf samtools-1.17.tar.bz2 &&\
    cd samtools-1.17 &&\
    make -lcurl/curl &&\
    export PATH=$PATH:/home/runner/work/dreem/dreem/samtools-1.17/ 

## setup correct python and packages
RUN conda env create --file env.yml 
ENV PATH /opt/conda/envs/dreem/bin:$PATH
RUN /bin/bash -c "source activate dreem"
RUN echo "source activate dreem" > ~/.bashrc

#WORKDIR /data

# Change user
USER biodocker


CMD ["dreem"]