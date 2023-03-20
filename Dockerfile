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
COPY . /
#COPY requirements.txt /
RUN apt-get -y update

RUN conda install python=3.10
RUN apt-get install gcc python3-dev -y
RUN apt-get install -y libzbar-dev
RUN apt-get install cutadapt -y

RUN pip3 install -r requirements.txt 
RUN pip3 install Cython==3.0.0b1
RUN pip3 install .

#RUN conda create --name dreem_env python=3.10
#ENV PATH /opt/conda/envs/dreem_env:$PATH



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
 make install

ENV LC_ALL=C

#fastqc install
RUN conda install -c bioconda fastqc







# We need to define the command to launch when we are going to run the image.
# We use the keyword 'CMD' to do that.
# The following command will execute "python ./main.py".
