# Changelog

## Unreleased

`pip install hicstuff` alone now runs the whole pipeline, from FASTQ files to
matrices, on Linux and macOS: with `-a bwa`, no other program is needed.

### Breaking changes

- hicstuff needs Python 3.10 or later (#122).

### New

- bwa runs in-process through [pybwa](https://github.com/fulcrumgenomics/pybwa),
  which is now a dependency of hicstuff (#125). `-a bwa` runs the `bwa`
  program when it is on PATH, and pybwa otherwise.
  `--aligner-backend {auto,binary,python}` on `pipeline`, `iteralign` and
  `missview`, or `HICSTUFF_ALIGNER_BACKEND`, chooses.
- Mapping mode `5spm` (`-m 5spm`) aligns the mates with `bwa mem -5SPM`: the
  5'-most alignment of a read cut by a ligation junction is its primary one,
  as in the 4DN, Juicer and Arima pipelines (#125). It aligns with bwa only,
  through pybwa unless `--aligner-backend binary`, and is not the default.
  Its pairs record each mate at its aligned 5' end, the end of the sequenced
  molecule, as pairtools does; the other modes record the leftmost base of
  each mate's alignment (#129).
- Pairs are sorted in Python when no usable `sort` program is installed (none
  on PATH, or busybox's) (#124). `HICSTUFF_SORT=python` always sorts them in
  Python.

### Changed

- FASTQ reads are counted, and the genome extracted from a bowtie2 index is
  compressed, in Python: gzip, zcat, pigz and wc are no longer needed (#123).
- Cutsite mode reads FASTQ files with pysam: pyfastx is no longer a
  dependency (#123).
- When the aligner is missing, the error also suggests `-a bwa`, which needs
  no program (#125).
- An installation of hicstuff now contains GPLv3 code: every pybwa wheel
  compiles in bwa (#125). See [Licence](installation.md#licence).
- The Docker image no longer ships samtools, which hicstuff has not needed
  since it runs the samtools that pysam bundles (#127).
- The package declares its licence as an SPDX expression (PEP 639) (#127).
- The documentation has moved from Read the Docs to
  <https://koszullab.github.io/hicstuff/>, rewritten in Markdown and built
  with Zensical: guides to the pipeline, the aligners, digestion, the 3C event
  filter, the distance law and the commands on matrices, a reference of every
  command and option, the file formats, and the Python API from the
  docstrings. The Sphinx sources (`doc/`) and the demo notebooks are gone, the
  README is shorter, and the `docs` extra (`uv sync --extra docs`) installs
  what builds the site instead of the `dev` extra's Sphinx.

### Fixed

- Iterative mapping with `-a mm2`, `-a minimap` or `-a bt2` aligned with
  bowtie2 in its final round (#125).

### Outputs

Runs with the same programs installed write the same pairs, matrices and logs
as hicstuff 3.3.0. Without the bwa program, `-a bwa` aligns through pybwa,
and the log says so. Its alignments then equal the bwa program's, except for
reads with several alignments of equal score beyond the first 10,000 reads of
a file, which pybwa can place on another of them. In our tests this moved
0.17% of the alignments, none with a mapping quality of 30 or more.
