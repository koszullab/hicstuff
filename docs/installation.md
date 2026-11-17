# Installation

hicstuff needs Python 3.10 or later, on Linux or macOS. On Windows, use WSL2:
pysam, one of its dependencies, publishes no Windows wheels.

## From PyPI

=== "pip"

    ```bash
    pip install hicstuff
    ```

=== "uv"

    ```bash
    uv tool install hicstuff   # the hicstuff command, in an environment of its own
    uv add hicstuff            # or as a dependency of your project
    ```

That alone runs the whole pipeline, from FASTQ files to matrices, with bwa:

```bash
hicstuff pipeline --aligner bwa --enzyme DpnII --genome genome.fa R1.fq.gz R2.fq.gz
```

needs no other program. bwa then runs in-process through
[pybwa](https://github.com/fulcrumgenomics/pybwa), which pip installs with
hicstuff.

!!! tip "Check it worked"

    ```bash
    hicstuff --version
    ```

## Aligners and other programs

`--aligner` chooses the aligner that maps the reads:

| Aligner | How hicstuff runs it | To install |
|---|---|---|
| `bowtie2` (the default) | The `bowtie2` program | Not on PyPI: `apt-get install bowtie2`, `conda install -c bioconda bowtie2` or `brew install bowtie2` |
| `bwa` | The `bwa` program if it is on `PATH`, otherwise in-process through pybwa. `--aligner-backend` (or `HICSTUFF_ALIGNER_BACKEND`) chooses | Nothing: pybwa comes with hicstuff |
| `minimap2` | The `minimap2` program | Not on PyPI: `apt-get install minimap2`, conda or brew |

Mapping mode `5spm` aligns with bwa only, through pybwa unless
`--aligner-backend binary`. See [Aligning reads](guides/alignment.md).

`hicstuff iteralign` and `hicstuff missview` always run the aligner, and
`hicstuff pipeline` only to start from FASTQ files: `--start-stage bam`,
`pairs` or `pairs_idx` needs no aligner, except `bowtie2-inspect` (part of
bowtie2) when the genome is given as a bowtie2 index.

No other program is needed. Two are used when they are installed:

- `bgzip` (from htslib) compresses the final pairs file on several threads,
  in the BGZF format, which pairix can index. Without it, pysam compresses it
  on one thread.
- GNU, BSD or uutils `sort` sorts the pairs. Without one, they are sorted in
  Python on one thread (busybox `sort` lacks the options hicstuff needs). Set
  `HICSTUFF_SORT=python` to always sort them in Python.

`samtools` is not needed: hicstuff runs the samtools that pysam bundles. On
Debian or Ubuntu (the `tabix` package provides bgzip):

```bash
apt-get install bowtie2 bwa minimap2 tabix
```

## Platforms that need a compiler

pip installs wheels where they exist, and builds the other dependencies from
source, which needs a C compiler, make and the zlib headers. As of November
2026:

- **macOS with Python 3.14:** pybwa publishes no wheel for it yet. Install
  Xcode's command-line tools (`xcode-select --install`), or use Python 3.13.
- **macOS on Intel:** biopython has published no Intel wheel since 1.87.
  `pip install --prefer-binary hicstuff` installs the last version that has
  one.
- **Linux with an older glibc**, such as RHEL or CentOS 7: several
  dependencies, pybwa and pysam among them, publish no wheel for it. Building
  pysam also needs the bzip2 and xz (lzma) headers. Conda is easier there.
- **Alpine and other musl-based Linux:** biopython, and matplotlib on
  aarch64, publish no musl wheels.

## With conda

hicstuff is also on [bioconda](https://bioconda.github.io/recipes/hicstuff/README.html).
To install it with the aligners and bgzip:

```bash
conda install -c conda-forge -c bioconda hicstuff bowtie2 bwa minimap2 htslib
```

## With Docker

Docker images, with hicstuff, bowtie2, bwa, minimap2 and htslib, are built and
published to the GitHub Container Registry when a release is tagged:

```bash
docker pull ghcr.io/koszullab/hicstuff:latest
docker pull ghcr.io/koszullab/hicstuff:<version>
```

## Development version

```bash
git clone https://github.com/koszullab/hicstuff.git
cd hicstuff
uv sync --extra dev
```

`uv sync --extra docs` adds what builds this documentation; see
[Contributing](contributing.md).

## Licence

hicstuff is distributed under the BSD-3-Clause licence. pybwa's metadata says
MIT, but every pybwa wheel compiles in bwa, which is distributed under the
GPLv3: an installation of hicstuff therefore contains GPLv3 code. This matters
only to whoever redistributes it, in a Docker image or a packed environment
for instance, who must then meet the GPLv3's terms for bwa. The Docker image
also ships the bowtie2 and bwa programs, both under the GPLv3.
