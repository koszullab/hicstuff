# hicstuff

[![PyPI version](https://img.shields.io/pypi/v/hicstuff.svg)](https://pypi.org/project/hicstuff/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/hicstuff.svg)](https://pypi.org/project/hicstuff/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/hicstuff/README.html)
[![CI](https://github.com/koszullab/hicstuff/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/koszullab/hicstuff/actions/workflows/ci.yml)
[![Docs](https://github.com/koszullab/hicstuff/actions/workflows/docs.yml/badge.svg?branch=master)](https://koszullab.github.io/hicstuff/)
[![codecov](https://codecov.io/gh/koszullab/hicstuff/branch/master/graph/badge.svg)](https://codecov.io/gh/koszullab/hicstuff)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2620608.svg)](https://doi.org/10.5281/zenodo.2620608)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

**Hi-C contact maps, from reads to matrices.** hicstuff aligns Hi-C reads,
pairs them, filters out the uninformative 3C events and builds the contact
map, in one command, for any genome. It also views, rebins, subsamples and
converts contact maps, computes the distance law, and exposes every step as a
Python library.

**Documentation: <https://koszullab.github.io/hicstuff/>**

## Installation

```sh
pip install hicstuff
```

hicstuff needs Python 3.10 or later, on Linux or macOS. That alone runs the
whole pipeline with `--aligner bwa`, which aligns in-process through
[pybwa](https://github.com/fulcrumgenomics/pybwa). bowtie2 (the default
aligner) and minimap2 are separate programs. hicstuff is also on
[bioconda](https://bioconda.github.io/recipes/hicstuff/README.html), and in a
Docker image, `ghcr.io/koszullab/hicstuff`. See
[Installation](https://koszullab.github.io/hicstuff/installation/) for the
aligners, the platforms that need a compiler, and the licence.

## Quick start

```sh
hicstuff pipeline --aligner bwa --enzyme DpnII --genome genome.fa --outdir out/ R1.fq.gz R2.fq.gz
```

writes the contact map at restriction-fragment resolution,
`out/abs_fragments_contacts_weighted.cool`, with the pairs it was built from,
the fragments and a log. `--binning 1000` writes a balanced multi-resolution
`.mcool` instead, at 1 kb and up; `--filter` and `--duplicates` remove the
uncuts, loops and PCR duplicates.

```sh
hicstuff stats out/hicstuff_*.log
hicstuff view --binning 10kb --normalize --output map.png out/abs_fragments_contacts_weighted.cool
```

The [Quickstart](https://koszullab.github.io/hicstuff/quickstart/) runs these
steps on the test data of this repository.

## Commands

| Command | |
|---|---|
| [`pipeline`](https://koszullab.github.io/hicstuff/cli/#pipeline) | Run the whole pipeline, from reads, BAM or pairs files to a contact map |
| [`digest`](https://koszullab.github.io/hicstuff/cli/#digest) | Cut a genome into restriction fragments |
| [`filter`](https://koszullab.github.io/hicstuff/cli/#filter) | Remove uncuts, loops and impossible pairs from a pairs file |
| [`iteralign`](https://koszullab.github.io/hicstuff/cli/#iteralign) | Align reads iteratively, from 20 bases up |
| [`cutsite`](https://koszullab.github.io/hicstuff/cli/#cutsite) | Cut reads at their ligation sites |
| [`rebin`](https://koszullab.github.io/hicstuff/cli/#rebin) | Write a contact map at a coarser resolution |
| [`subsample`](https://koszullab.github.io/hicstuff/cli/#subsample) | Draw contacts at random from a contact map |
| [`convert`](https://koszullab.github.io/hicstuff/cli/#convert) | Convert a contact map between cool, 2D bedgraph and graal |
| [`view`](https://koszullab.github.io/hicstuff/cli/#view) | Draw a contact map, or the ratio of two |
| [`scalogram`](https://koszullab.github.io/hicstuff/cli/#scalogram) | Draw contacts by distance along the genome |
| [`distancelaw`](https://koszullab.github.io/hicstuff/cli/#distancelaw) | Compute and plot the distance law, P(s) |
| [`missview`](https://koszullab.github.io/hicstuff/cli/#missview) | Show the bins reads of a given length cannot reach |
| [`stats`](https://koszullab.github.io/hicstuff/cli/#stats) | Summarise the log of a pipeline run |

## Contributing

Bug reports, suggestions and pull requests are welcome. See
[Contributing](https://koszullab.github.io/hicstuff/contributing/) for the
tests, the documentation and releases.

## Citation

Cyril Matthey-Doret, Lyam Baudry, Amaury Bignaud, Axel Cournac, Remi-Montagne,
Nadège Guiglielmoni, Théo Foutel Rodier and Vittore F. Scolari. 2020.
hicstuff: Simple library/pipeline to generate and handle Hi-C data. Zenodo.
<http://doi.org/10.5281/zenodo.4066363>. BibTeX:
[Citing hicstuff](https://koszullab.github.io/hicstuff/citation/).

## Licence

BSD-3-Clause. An installation of hicstuff also contains bwa, under the GPLv3,
through pybwa: see [Licence](https://koszullab.github.io/hicstuff/installation/#licence).
