FROM continuumio/miniconda3:4.8.2

LABEL Name=hicstuff Version=2.3.1

COPY * ./ /app/
WORKDIR /app

RUN conda update -y conda
RUN conda config --add channels bioconda

# Get 3rd party packages directly from conda
RUN conda install -c conda-forge -y \
    pip \
    bowtie2 \
    minimap2 \
    bwa \
    samtools \
    htslib \
    pysam \ 
    cooler && conda clean -afy

RUN pip install -Ur requirements.txt
# Using pip:
RUN pip install .
#CMD ["python3", "-m", "hicstuff.main"]
ENTRYPOINT [ "hicstuff" ]
