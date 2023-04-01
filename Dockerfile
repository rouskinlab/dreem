# A dockerfile must always start by importing the base image.
# We use the keyword 'FROM' to do that.
# In our example, we want import the python image.
# So we write 'python' for the image name and 'latest' for the version.
FROM continuumio/miniconda3

# In order to launch our python code, we must import it into our image.
# We use the keyword 'COPY' to do that.
# The first parameter 'main.py' is the name of the file on the host.
# The second parameter '/' is the path where to put the file on the image.
# Here we put the file at the image root folder.
#COPY main.py /
RUN apt-get -y update && apt-get install -y libzbar-dev
RUN apt-get install unzip



COPY requirements.txt /
COPY . /

#RUN conda create --name dreem_env python=3.10
#ENV PATH /opt/conda/envs/dreem_env:$PATH
RUN apt-get -y update && apt-get install -y libzbar-dev
RUN conda install python=3.10
#RUN pip install -r requirements.txt 


RUN apt-get install unzip
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

COPY RNAstructureSource.zip /
RUN mkdir rnaS && \
unzip RNAstructureSource.zip -d /

RUN apt-get update && apt-get install make && \
cd RNAstructure 



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

ARG SAMTOOLSVER=1.16
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

RUN pip install . 

ENTRYPOINT ["dreem"]
CMD [ "/bin/bash" ]
