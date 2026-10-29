FROM python:3.8

# Configuration for tool versions
ARG HUMID_VERSION=v1.0.5
ARG CALIB_VERSION=v0.3.4
ARG UMI_NEA_VERSION=v1.0.2

RUN apt-get -y update && \
    apt-get -y install apt-utils sudo wget curl git vim bc less \
    gcc g++ cmake \
    libtool build-essential clang \
    zlib1g-dev libbz2-dev \
    libboost-all-dev libboost-iostreams-dev

RUN apt-get -y install samtools
RUN apt-get -y install bwa

COPY requirements.txt /tmp/

RUN pip install umi_tools
RUN pip install --requirement /tmp/requirements.txt
COPY network.py /usr/local/lib/python3.8/site-packages/umi_tools/

RUN mkdir -p /Download/
WORKDIR /Download/
RUN git clone -b $UMI_NEA_VERSION https://github.com/Qiaseq-research/UMI-nea.git
RUN cd UMI-nea/UMI-nea && make

RUN git clone https://github.com/fhlab/UMIC-seq.git
COPY UMIC-seq.py UMIC-seq/

RUN git clone -b $CALIB_VERSION https://github.com/vpc-ccg/calib.git calib
RUN cd calib && make

RUN wget https://github.com/jfjlaros/HUMID/releases/download/$HUMID_VERSION/humid -O /usr/local/bin/humid
RUN chmod +x /usr/local/bin/humid

WORKDIR /home/
CMD ["/bin/bash"]
