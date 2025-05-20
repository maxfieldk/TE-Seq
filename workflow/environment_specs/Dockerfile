FROM --platform=linux/amd64 condaforge/mambaforge-pypy3

RUN echo "Acquire::http::Pipeline-Depth 0;" > /etc/apt/apt.conf.d/99custom && \
    echo "Acquire::http::No-Cache true;" >> /etc/apt/apt.conf.d/99custom && \
    echo "Acquire::BrokenProxy    true;" >> /etc/apt/apt.conf.d/99custom

RUN apt-get update -q && apt-get install --no-install-recommends -qy python3-dev g++ gcc inetutils-ping 

COPY envs/qc.yaml /tmp/qc.yaml 
RUN mamba env create -f /tmp/qc.yaml  && mamba clean --all --yes 

COPY envs/omics.yaml /tmp/omics.yaml 
RUN mamba env create -f /tmp/omics.yaml  && mamba clean --all --yes 

COPY envs/rseqc.yaml /tmp/rseqc.yaml 
RUN mamba env create -f /tmp/rseqc.yaml  && mamba clean --all --yes 

COPY envs/snake3.yaml /tmp/snake3.yaml 
RUN mamba env create -f /tmp/snake3.yaml  && mamba clean --all --yes 

COPY envs/star.yaml /tmp/star.yaml 
RUN mamba env create -f /tmp/star.yaml  && mamba clean --all --yes 

COPY envs/telescope3.yaml /tmp/telescope3.yaml 
RUN mamba env create -f /tmp/telescope3.yaml  && mamba clean --all --yes 

COPY envs/tldr.yaml /tmp/tldr.yaml 
RUN mamba env create -f /tmp/tldr.yaml && git clone https://github.com/adamewing/tldr.git; cd tldr; conda run -n tldr pip install -e $PWD  && mamba clean --all --yes 

COPY envs/ea.yaml /tmp/ea.yaml 
RUN mamba env create -f /tmp/ea.yaml  && mamba clean --all --yes 

COPY envs/repeatanalysis.yaml /tmp/repeatanalysis.yaml 
RUN mamba env create -f /tmp/repeatanalysis.yaml && conda run --no-capture-output -n repeatanalysis R -e "devtools::install_github('BIMSBbioinfo/genomation',build_vignettes=FALSE)" && conda run --no-capture-output -n repeatanalysis R -e "devtools::install_github('krassowski/complex-upset')" && mamba clean --all --yes 

COPY envs/evo2.yaml /tmp/evo2.yaml 
RUN mamba env create -f /tmp/evo2.yaml && conda run --no-capture-output -n evo2 R -e "options(repos = c(CRAN='https://cloud.r-project.org')); if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('rBLAST')"  && mamba clean --all --yes 

COPY envs/ds.yaml /tmp/ds.yaml 
RUN mamba env create -f /tmp/ds.yaml && mamba clean --all --yes 

COPY envs/deeptools.yaml /tmp/deeptools.yaml 
RUN mamba env create -f /tmp/deeptools.yaml && mamba clean --all --yes 

COPY envs/deseq.yaml /tmp/deseq.yaml 
RUN mamba env create -f /tmp/deseq.yaml  && mamba clean --all --yes 

COPY envs/pycoQC.yaml /tmp/pycoQC.yaml 
RUN mamba env create -f /tmp/pycoQC.yaml  && mamba clean --all --yes 

COPY envs/rmtest.yaml /tmp/rmtest.yaml 
RUN mamba env create -f /tmp/rmtest.yaml  && mamba clean --all --yes 
