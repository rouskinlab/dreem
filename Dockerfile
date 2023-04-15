#this current 
FROM continuumio/miniconda3


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

#RUN cd dreem && \
#unzip RNAstructureSource.zip -d ../

RUN apt-get update && apt-get install make 



RUN apt-get update && apt-get install -y apt-utils
RUN apt-get install -y g++
RUN apt-get install libsimde-dev
#RUN cd RNAstructure && make all
#ENV DATAPATH="/RNAstructure/data_tables/"
#ENV PATH="$PATH:/RNAstructure/exe/"







ENV LC_ALL=C
#fastqc install
RUN conda install -c bioconda fastqc
#RUN git clone https://github.com/shenwei356/seqkit

RUN wget https://go.dev/dl/go1.20.3.linux-arm64.tar.gz

#RUN tar -zxf go1.20.3.linux-arm64.tar.gz -C deps/ && mv deps/go/bin/* bin/
RUN tar -zxf go1.20.3.linux-arm64.tar.gz -C /usr/local
#ENV PATH="$PATH:deps/go/bin" 
ENV PATH="$PATH:/usr/local/go/bin" 
RUN git clone https://github.com/shenwei356/seqkit 

#ENV GOROOT="deps/"
RUN cd seqkit/seqkit && go build main.go
#ENV PATH="$PATH:seqkit/seqkit"




RUN cd dreem && \ 
pip install . 

#ENTRYPOINT ["dreem"]
CMD [ "/bin/bash" ]








# We need to define the command to launch when we are going to run the image.
# We use the keyword 'CMD' to do that.
# The following command will execute "python ./main.py".
