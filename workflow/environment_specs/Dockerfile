FROM --platform=linux/amd64 condaforge/mambaforge-pypy3

RUN echo "Acquire::http::Pipeline-Depth 0;" > /etc/apt/apt.conf.d/99custom && \
    echo "Acquire::http::No-Cache true;" >> /etc/apt/apt.conf.d/99custom && \
    echo "Acquire::BrokenProxy    true;" >> /etc/apt/apt.conf.d/99custom

RUN apt-get update -q && apt-get install --no-install-recommends -qy python3-dev g++ gcc inetutils-ping 

ADD envs/qc.yaml /tmp/qc.yaml 
RUN --mount=type=cache,target=$MAMBA_ROOT_PREFIX mamba env create -f /tmp/qc.yaml  && mamba clean --all --yes 

ADD envs/omics.yaml /tmp/omics.yaml 
RUN --mount=type=cache,target=$MAMBA_ROOT_PREFIX mamba env create -f /tmp/omics.yaml  && mamba clean --all --yes 

ADD envs/rseqc.yaml /tmp/rseqc.yaml 
RUN --mount=type=cache,target=$MAMBA_ROOT_PREFIX mamba env create -f /tmp/rseqc.yaml  && mamba clean --all --yes 

ADD envs/snake3.yaml /tmp/snake3.yaml 
RUN --mount=type=cache,target=$MAMBA_ROOT_PREFIX mamba env create -f /tmp/snake3.yaml  && mamba clean --all --yes 

ADD envs/star.yaml /tmp/star.yaml 
RUN --mount=type=cache,target=$MAMBA_ROOT_PREFIX mamba env create -f /tmp/star.yaml  && mamba clean --all --yes 

ADD envs/telescope3.yaml /tmp/telescope3.yaml 
RUN --mount=type=cache,target=$MAMBA_ROOT_PREFIX mamba env create -f /tmp/telescope3.yaml  && mamba clean --all --yes 

ADD envs/tldr.yaml /tmp/tldr.yaml 
RUN --mount=type=cache,target=$MAMBA_ROOT_PREFIX mamba env create -f /tmp/tldr.yaml && git clone https://github.com/adamewing/tldr.git; cd tldr; conda run -n tldr pip install -e $PWD  && mamba clean --all --yes 

ADD envs/ea.yaml /tmp/ea.yaml 
RUN --mount=type=cache,target=$MAMBA_ROOT_PREFIX mamba env create -f /tmp/ea.yaml  && mamba clean --all --yes 

ADD envs/repeatanalysis.yaml /tmp/repeatanalysis.yaml 
RUN --mount=type=cache,target=$MAMBA_ROOT_PREFIX mamba env create -f /tmp/repeatanalysis.yaml && conda run --no-capture-output -n repeatanalysis R -e "devtools::install_github('BIMSBbioinfo/genomation',build_vignettes=FALSE)" && conda run --no-capture-output -n repeatanalysis R -e "devtools::install_github('krassowski/complex-upset')" && mamba clean --all --yes 

ADD envs/evo2.yaml /tmp/evo2.yaml 
RUN --mount=type=cache,target=$MAMBA_ROOT_PREFIX mamba env create -f /tmp/evo2.yaml  && mamba clean --all --yes 

ADD envs/ds.yaml /tmp/ds.yaml 
RUN --mount=type=cache,target=$MAMBA_ROOT_PREFIX mamba env create -f /tmp/ds.yaml && mamba clean --all --yes 

ADD envs/deeptools.yaml /tmp/deeptools.yaml 
RUN --mount=type=cache,target=$MAMBA_ROOT_PREFIX mamba env create -f /tmp/deeptools.yaml && mamba clean --all --yes 

ADD envs/deseq.yaml /tmp/deseq.yaml 
RUN --mount=type=cache,target=$MAMBA_ROOT_PREFIX mamba env create -f /tmp/deseq.yaml  && mamba clean --all --yes 

ADD envs/pycoQC.yaml /tmp/pycoQC.yaml 
RUN --mount=type=cache,target=$MAMBA_ROOT_PREFIX mamba env create -f /tmp/pycoQC.yaml  && mamba clean --all --yes 

ADD envs/rmtest.yaml /tmp/rmtest.yaml 
RUN --mount=type=cache,target=$MAMBA_ROOT_PREFIX mamba env create -f /tmp/rmtest.yaml  && mamba clean --all --yes 
