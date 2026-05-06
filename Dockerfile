FROM mambaorg/micromamba:latest

LABEL Name=hicstuff Version=3.2.4

COPY --chown=$MAMBA_USER:$MAMBA_USER . ./

## Install system / compiled tools (not available via pip)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python>=3.8 \
    bowtie2 bwa minimap2 samtools htslib && \
    micromamba clean --all --yes

## Install uv and hicstuff Python package
RUN micromamba run pip install uv && \
    micromamba run uv pip install --system .

WORKDIR /home/mambauser/
ENTRYPOINT [ "/bin/bash" ]
