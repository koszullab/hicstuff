FROM mambaorg/micromamba:latest

LABEL Name=hicstuff Version=3.2.0

COPY --chown=$MAMBA_USER:$MAMBA_USER . ./

## Install dependencies
RUN micromamba install -y -n base --file environment.yml && \
    micromamba install -y -n base pip && \
    micromamba clean --all --yes

## Install hicstuff
RUN micromamba run python3 -m pip install -e .

WORKDIR /home/mambauser/
ENTRYPOINT [ "/bin/bash" ]
