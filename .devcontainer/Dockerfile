FROM mcr.microsoft.com/devcontainers/python:1-3.12-bullseye

RUN apt-get install -y zip \
    unzip

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir numpy pandas biopython matplotlib scipy

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
    