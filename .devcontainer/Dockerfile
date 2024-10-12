FROM mcr.microsoft.com/devcontainers/python:1-3.12-bullseye

RUN apt-get install -y zip \
    unzip

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir numpy pandas biopython matplotlib cvxpy

RUN apt-get update && \
    apt-get install -y openjdk-11-jdk && \
    apt-get install -y ant && \
    apt-get clean;

RUN apt-get update && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;

ENV JAVA_HOME /usr/lib/jvm/java-11-openjdk-amd64/
RUN export JAVA_HOME

WORKDIR /CiiiDER

RUN wget https://ciiider.erc.monash.edu/files/CiiiDER_TFMs.zip && \
    unzip CiiiDER_TFMs.zip && \
    sudo chown vscode:vscode /CiiiDER/CiiiDER_TFMs/
