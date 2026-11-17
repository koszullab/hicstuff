FROM mambaorg/micromamba:latest

LABEL Name=hicstuff Version=3.3.0

COPY --chown=$MAMBA_USER:$MAMBA_USER . ./

## Install the aligners, which are not on PyPI, and htslib for bgzip. bwa also
## runs in-process through pybwa, a dependency of hicstuff; samtools is not
## needed: hicstuff runs the samtools that pysam bundles.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    "python=3.13" \
    bowtie2 bwa minimap2 htslib && \
    micromamba clean --all --yes

## Install uv and hicstuff Python package
RUN micromamba run pip install uv && \
    micromamba run uv pip install --system .

WORKDIR /home/mambauser/
ENTRYPOINT [ "/bin/bash" ]
