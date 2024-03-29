# syntax=docker/dockerfile:1
#junk



FROM continuumio/miniconda3

ARG TARGETOS
ARG TARGETARCH


RUN apt-get -y update && apt-get install -y libzbar-dev
RUN apt-get install unzip






#RUN conda create --name dreem_env python=3.10
#ENV PATH /opt/conda/envs/dreem_env:$PATH
RUN apt-get -y update && apt-get install -y libzbar-dev
RUN conda install python=3.10
#RUN pip install -r requirements.txt 
RUN git clone https://github.com/rouskinlab/dreem.git


RUN apt-get install unzip
RUN apt install -y curl
#install samtools 
RUN apt-get update && apt-get install --no-install-recommends -y \
 libncurses5-dev \
 libbz2-dev \
 liblzma-dev \
 libcurl4-gnutls-dev \
 zlib1g-dev \
 libssl-dev \
 gcc \
 wget \
 make \
 perl \
 bzip2 \
 gnuplot \
 ca-certificates \
 gawk && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*



RUN apt-get update && apt-get install make 


RUN cd dreem && \
unzip RNAstructureSource.zip -d ../
RUN apt-get update && apt-get install -y apt-utils
RUN apt-get install -y g++
RUN apt-get install libsimde-dev
RUN cd RNAstructure && make all
ENV DATAPATH="/RNAstructure/data_tables/"
ENV PATH="$PATH:/RNAstructure/exe/"

ENV ZIP=bowtie2-2.4.5-source.zip
ENV URL=https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.5/bowtie2-2.4.5-source.zip/download/
ENV FOLDER=bowtie2-2.4.5-source
ENV SOURCE=bowtie2-2.4.5
ENV DST=deps
ENV BT2=$PATH:$DST/$FOLDER/
RUN mkdir $DST


RUN wget -q $URL -O $DST/$ZIP && \
    unzip $DST/$ZIP -d $DST && \
    rm $DST/$ZIP && \
    cd $DST/$SOURCE && \
    make
    #mv $DST/$FOLDER/* /bin 
ENV PATH="$PATH:/deps/bowtie2-2.4.5"

ARG SAMTOOLSVER=1.16.1
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLSVER}/samtools-${SAMTOOLSVER}.tar.bz2 && \
 tar -xjf samtools-${SAMTOOLSVER}.tar.bz2 && \
 rm samtools-${SAMTOOLSVER}.tar.bz2 && \
 cd samtools-${SAMTOOLSVER} && \
 ./configure && \
 make && \
 make install -d 

ENV LC_ALL=C
#fastqc install
RUN conda install -c bioconda fastqc
#RUN git clone https://github.com/shenwei356/seqkit

RUN wget https://go.dev/dl/go1.20.3.${TARGETOS}-${TARGETARCH}.tar.gz

#RUN tar -zxf go1.20.3.linux-arm64.tar.gz -C deps/ && mv deps/go/bin/* bin/
RUN tar -zxf go1.20.3.${TARGETOS}-${TARGETARCH}.tar.gz -C /usr/local
#ENV PATH="$PATH:deps/go/bin" 
ENV PATH="$PATH:/usr/local/go/bin" 
RUN git clone https://github.com/shenwei356/seqkit 
#RUN go get -v -u github.com/shenwei356/seqkit/seqkit

RUN cd seqkit/seqkit && go build
ENV PATH="$PATH:seqkit/seqkit"




#RUN cd dreem && \ 
#pip install . 








RUN cd dreem && \ 
pip install . 


ENTRYPOINT ["dreem"]
CMD [ "/bin/bash" ]





