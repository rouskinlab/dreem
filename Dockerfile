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
COPY requirements.txt /

#RUN conda create --name dreem_env python=3.10
#ENV PATH /opt/conda/envs/dreem_env:$PATH
RUN apt-get -y update && apt-get install -y libzbar-dev
RUN conda install python=3.10
RUN pip install -r requirements.txt 

RUN apt-get install unzip
ENV ZIP=bowtie2-2.4.5-macos-arm64.zip
ENV URL=https://github.com/BenLangmead/bowtie2/releases/download/v2.4.5/
ENV FOLDER=bowtie2-2.4.5-macos-arm64
ENV DST=deps
ENV BT2=$PATH:$DST/$FOLDER/
RUN mkdir $DST

RUN wget -q $URL/$ZIP -O $DST/$ZIP && \
    unzip $DST/$ZIP -d $DST && \
    rm $DST/$ZIP && \
    mv $DST/$FOLDER/* /bin 

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



COPY RNAstructureSource.zip /
RUN mkdir rnaS && \
unzip RNAstructureSource.zip -d /

RUN apt-get update && apt-get install make && \
cd RNAstructure 



RUN apt-get update && apt-get install -y apt-utils
RUN apt-get install -y g++
RUN cd RNAstructure && make all
ENV DATAPATH="/RNAstructure/data_tables/"
ENV PATH="$PATH:/RNAstructure/exe/"




#RUN conda install rnastructure
#COPY RNAstructure/ /
#ENV DATAPATH="/RNAstructure/data_tables/"
#ENV PATH="$PATH:/RNAstructure/exe"
CMD ["Fold","-h"]







# We need to define the command to launch when we are going to run the image.
# We use the keyword 'CMD' to do that.
# The following command will execute "python ./main.py".