# Install Ubuntu as base
FROM ubuntu:latest

# Install R(4.3.3(latest-ensure this is constant)), Python(3.12), Java(openjdk-21), some bioinf tools
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
apt-get install -y --no-install-recommends build-essential python3-pip python3-setuptools python3-dev wget \
bowtie2 fastqc samtools trim-galore picard-tools

# Allow installation of packages via pip
RUN mv /usr/lib/python3.12/EXTERNALLY-MANAGED /usr/lib/python3.12/EXTERNALLY-MANAGED.old

RUN pip3 install \
MACS3==3.0.3 \
deeptools