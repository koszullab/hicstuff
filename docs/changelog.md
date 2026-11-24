# Changelog

## Unreleased

### Breaking changes

- From 50,000 contigs in the map, `hicstuff pipeline --binning` without
  `--matrix` writes a binned metacool file,
  `abs_fragments_contacts_weighted.metacool`, rather than a zoomified and
  balanced `.mcool`, and warns of it: the default metacool map no longer stops
  at `--binning`, as metacool now takes it (see New). `-z true` and
  `--balancing-args`, options of the cool map, given then, stop the run once
  the genome is read, before the reads are aligned, with an error that says to
  give `--matrix cool`; `full_pipeline` raises `DefaultMatrixError`, a
  `ValueError`, there.
- `-z true` is refused with `--matrix metacool`, as `--balancing-args` is: a
  metacool map is neither zoomified nor balanced. `--zoomify` left to its
  default is ignored there. `full_pipeline`'s `zoomify` defaults to `None`,
  which zoomifies a cool map binned at a bin size, as `True` did; `True` is
  refused with `metacool`.

### New

- `hicstuff pipeline --matrix metacool` takes `--binning`. The pairs are
  binned from their positions, as for the cool map: the metacool file holds
  the contacts that the cool file of the same bins holds, with metacooler's
  contact statistics counted on the pairs, which metaTOR's `mge` reads, but is
  neither zoomified nor balanced. `pairs2metacool` takes `binning` and
  `info_contigs`.
- `--binning` takes a BED file of bins as well as a bin size in bp, for the
  cool and metacool maps: the contig, start and end of each bin, tab-separated,
  without a header line, as `cooler makebins` and `cooler digest` write them,
  or one bin per contig, for the map of the contigs. The bins of each contig
  it names must tile it, from 0 to its end, without gap or overlap; the
  contigs it does not name, and those that `--size` and `--exclude` leave out,
  are left out of the map. The file is read before the run starts, and checked
  against the genome before the reads are aligned. The cool map of a BED file
  is balanced, but not zoomified: `-z true` is refused with it. In Python,
  `full_pipeline`'s `binning` and `pairs2binnedcool`'s take its path, and
  `parse_binning`, `read_bins` and `matrix_bins` read it.

### Changed

- The bowtie2 or bwa index of a FASTA file is also found named after it
  without its extensions, as `bowtie2-build genome.fa genome` and
  `bwa index -p genome genome.fa` write it: for `--genome genome.fa.gz`,
  hicstuff looks for the index at `genome.fa.gz`, then `genome.fa`, then
  `genome`, and the log names the one it found. It used to build a temporary
  index at every run. This applies to `pipeline`, `iteralign` and `missview`,
  and to `check_fasta_index`, which returns the prefix it found.
- `--binning` is checked before the run starts, whatever the matrix format: a
  negative number, or text that is neither digits nor a file, is a usage error.
  A negative binning with `--matrix cool` used to fail in cooler, once the
  pairs were made.

## 4.0.0 - 2026-11-23

`pip install hicstuff` alone now runs the whole pipeline, from FASTQ files to
matrices, on Linux and macOS: with `-a bwa`, no other program is needed.

### Breaking changes

- hicstuff needs Python 3.10 or later (#122).
- With 50,000 contigs or more in the map (those that `--size` and `--exclude`
  keep), as a metagenome assembly may have, `hicstuff pipeline` writes its
  contact map as a metacool file, `abs_fragments_contacts_weighted.metacool`,
  rather than a cool file, unless given `--binning` or `--matrix cool`, and
  warns of it. See `--matrix metacool` below.

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
- The log of `hicstuff pipeline` records the command that ran it, quoted to
  paste in a shell, on a `## cmd:` line after the date.
- `hicstuff pipeline` reads several pairs of FASTQ files, such as the lanes of
  a library, as one library: the R1 files, then the R2 files, comma-separated
  or as glob patterns, quoted or not (`hicstuff pipeline ... lib1_*_R1.fq.gz
  lib1_*_R2.fq.gz`). The R1 and R2 files of each pair must start with the same
  read. They are concatenated in the temporary directory, and give the outputs
  of the library in one pair of files. `full_pipeline` takes lists of files
  for `input1` and `input2`.
- The `--genome` of `hicstuff pipeline` can be the prefix of a bwa index
  without its FASTA, as it could be of a bowtie2 index: hicstuff reads the
  sequence from the index's `.pac` and `.amb` files, with no program.
- `hicstuff pipeline --matrix metacool` (`-M metacool`) writes the contact map
  as a [metacooler](https://pypi.org/project/metacooler/) file, made by
  `metacooler create pairs`: a cool file with one bin per fragment, which
  cooler reads, with metacooler's tables of genomes and contigs, its contigs
  not yet assigned to genomes. Unlike cooler, metacooler does not store the
  contig of each bin as an HDF5 enum, whose cost grows with the square of the
  number of contigs. metacooler, 0.4.0 or later, is now a dependency of
  hicstuff: the contact statistics it counts on the pairs (`qc`, which
  hicstuff leaves to metacooler's defaults) split the contacts within a contig
  at 1001 bp, so that those whose reads are more than 1 kb apart are
  long-range, as metaTOR's `mge` reads them from the map. It is the
  default map from 50,000 contigs (see Breaking changes), but with
  `--binning`, which it does not take: `--binning` and `--balancing-args`,
  options of the cool map, are refused with it, but for `--binning 0` and blank
  values. `full_pipeline` takes `mat_fmt="metacool"`, and its default
  `mat_fmt=None` chooses between cool and metacool as the command does.
- `hicstuff pipeline --matrix none` (`-M none`) writes no contact map: the run
  ends with its pairs file, for the tools that build their own maps from pairs
  (metaTOR, metacooler). `--binning`, `--balancing-args` and `--exclude`,
  options of the cool map, are refused with it, but for `--binning 0` and blank
  values; so are `binning`, `balancing_args` and `exclude` with
  `mat_fmt="none"` in `full_pipeline`.
- With 50,000 contigs or more in a cool map (`--matrix cool`, or with
  `--binning`), `hicstuff pipeline` warns, before it aligns reads or builds an
  index, that the map will be slow to write: cooler takes a time that grows
  with the square of the number of contigs for every cooler it writes (about
  30 s at 117,000 contigs). The warning names the formats that do not pay it.

### Changed

- `--matfmt` of `hicstuff pipeline` is now `--matrix`, still `-M`. `--matfmt`
  still works, but the help no longer shows it.
- `hicstuff pipeline --help` lists its main options in a panel of their own,
  before the advanced ones, and describes each mapping mode.
- FASTQ reads are counted, and the genome extracted from a bowtie2 index is
  compressed, in Python: gzip, zcat, pigz and wc are no longer needed (#123).
- Cutsite mode reads FASTQ files with pysam: pyfastx is no longer a
  dependency (#123).
- When the aligner is missing, the error also suggests `-a bwa`, which needs
  no program (#125).
- With `-a bwa`, the log warns of a bwa-mem2 index of the genome, which bwa
  and pybwa cannot read, before building a bwa index.
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
  docstrings. It recommends mapping mode `5spm`, which its examples use. The
  Sphinx sources (`doc/`) and the demo notebooks are gone, the README is
  shorter, and the `docs` extra (`uv sync --extra docs`) installs what builds
  the site instead of the `dev` extra's Sphinx.

### Fixed

- Iterative mapping with `-a mm2`, `-a minimap` or `-a bt2` aligned with
  bowtie2 in its final round (#125).
- The pipeline indexes a FASTA without index where it is, plain or
  compressed, in a directory of its own in the temporary directory, which the
  cleanup removes. It copied the FASTA to the temporary directory first,
  compressed, also in runs that align no reads, and left there the index
  files it built.
- A `--genome` that is neither a FASTA file nor the prefix of an index the
  aligner reads stops the pipeline before it writes anything. The error was
  logged, and the run went on until a later step failed.
- `full_pipeline` refuses a `mat_fmt` other than `graal`, `bg2`, `cool`,
  `metacool` and `none` before the run starts, logged as an error. Another
  name wrote an empty matrix file under its extension.

### Outputs

Runs with the same programs installed write the same pairs, matrices and logs
as hicstuff 3.3.0, but for the `## cmd:` line of the pipeline's log, and
for a map of 50,000 contigs or more: a metacool file by default, and a warning
with a cool file. Without the bwa program, `-a bwa` aligns through pybwa,
and the log says so. Its alignments then equal the bwa program's, except for
reads with several alignments of equal score beyond the first 10,000 reads of
a file, which pybwa can place on another of them. In our tests this moved
0.17% of the alignments, none with a mapping quality of 30 or more.
