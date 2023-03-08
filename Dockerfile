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

ENV ZIP=bowtie2-2.2.9-linux-x86_64.zip
ENV URL=https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/
ENV FOLDER=bowtie2-2.2.9
ENV DST=/home/biodocker/bin

USER root

## install bowtie2
RUN wget -q $URL/$ZIP -O $DST/$ZIP && \
    unzip $DST/$ZIP -d $DST && \
    rm $DST/$ZIP && \
    mv $DST/$FOLDER/* $DST && \
    rmdir $DST/$FOLDER

## install fastqc
RUN apt-get update && \
    apt-get -y -q install fastqc

## setup correct python and packages
RUN conda env create --file env.yml 
ENV PATH /opt/conda/envs/dreem/bin:$PATH
RUN /bin/bash -c "source activate dreem"
RUN echo "source activate dreem" > ~/.bashrc

## install cutadapt
RUN wget -q https://github.com/marcelm/cutadapt/archive/refs/tags/v1.18.zip -O $DST/v1.18.zip && \
    unzip $DST/v1.18.zip -d $DST && \
    rm $DST/v1.18.zip
RUN /bin/bash -c "cd $DST/cutadapt-1.18 && python setup.py install"

## install trim_galore
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz && \
    tar xvzf trim_galore.tar.gz && \
    cp TrimGalore-0.6.6/trim_galore $DST && \
    rm -rf TrimGalore-0.6.6 trim_galore.tar.gz

RUN wget https://github.com/jyesselm/dreem/archive/refs/tags/0.2.0.zip      && \
    unzip 0.2.0.zip                                                         && \
    mv dreem-0.2.0 dreem                                                    && \
   /bin/bash -c "cd dreem && python setup.py install"

WORKDIR /data

# Change user
USER biodocker


CMD ["dreem"]