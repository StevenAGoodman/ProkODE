FROM mcr.microsoft.com/devcontainers/python:1-3.12-bullseye

RUN apt-get install -y zip \
    unzip

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir numpy pandas biopython matplotlib scipy GEOparse

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

WORKDIR /convert_pscm

RUN wget http://trap.molgen.mpg.de/download/TRAP/convert-pscm_source.zip && \
    unzip convert-pscm_source.zip && \
    cd convert-pscm_source/Release/ && \
    sed -ie "s/while(matfile.getline(mat_line,1000,'\\\n')!=NULL)/while(matfile.getline(mat_line,1000,'\\\n'))/g" ../read_matrix.cpp && \
    make all && \
    ln -s /convert_pscm/convert-pscm_source/Release/convert-pscm /usr/local/sbin/convert-pscm

WORKDIR /annotate

RUN wget http://trap.molgen.mpg.de/download/TRAP/ANNOTATE-3.04.01.tar.gz && \
    tar -xvzf ANNOTATE-3.04.01.tar.gz && \
    cd ANNOTATE-3.04.01/Release/ && \
    sed -ie "s/while(matfile.getline(mat_line,1000,'\\\n')!=NULL)/while(matfile.getline(mat_line,1000,'\\\n'))/g" ../read_matrix.cpp && \
    make all && \
    ln -s /annotate/ANNOTATE-3.04.01/Release/ANNOTATE_v3.04 /usr/local/sbin/annotate
    
WORKDIR /CiiiDER

RUN wget https://ciiider.erc.monash.edu/files/CiiiDER_TFMs.zip && \
    unzip CiiiDER_TFMs.zip && \
    sudo chown vscode:vscode /CiiiDER/CiiiDER_TFMs/
