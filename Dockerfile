
# build stage
FROM continuumio/miniconda3

ARG TARGETPLATFORM
ARG TARGETARCH
ARG TARGETOS

RUN echo "i am building for $TARGETARCH"
RUN wget https://go.dev/dl/go1.20.3.${TARGETOS}-${TARGETARCH}.tar.gz
#FROM continuumio/miniconda3

#COPY --from=build-env 
#RUN apt-get -y update && apt-get install -y libzbar-dev
#RUN apt-get install unzip






