# hicstuff

[![PyPI version](https://badge.fury.io/py/hicstuff.svg)](https://badge.fury.io/py/hicstuff)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/hicstuff.svg)
[![Build Status](https://api.travis-ci.com/koszullab/hicstuff.svg)](https://app.travis-ci.com/github/koszullab/hicstuff)
[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/koszullab/hicstuff)](https://hub.docker.com/r/koszullab/hicstuff)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2620608.svg)](https://doi.org/10.5281/zenodo.2620608)
[![codecov](https://codecov.io/gh/koszullab/hicstuff/branch/master/graph/badge.svg)](https://codecov.io/gh/koszullab/hicstuff)
[![Read the docs](https://readthedocs.org/projects/hicstuff/badge)](https://hicstuff.readthedocs.io)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/koszullab/hicstuff/master?filepath=doc%2Fnotebooks%2Fdemo_api.ipynb)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

A lightweight library that generates and handles Hi-C contact maps in either cooler-compatible 2Dbedgraph or [instaGRAAL](https://github.com/koszullab/instaGRAAL) format. It is essentially a merge of the [yahcp](https://github.com/baudrly/yahcp) pipeline, the [hicstuff](https://github.com/baudrly/hicstuff) library and extra features illustrated in the [3C tutorial](https://github.com/axelcournac/3C_tutorial) and the [DADE pipeline](https://github.com/scovit/dade), all packaged together for extra convenience.

The goal is to make generation and manipulation of Hi-C matrices as simple as possible and work for any organism.

## Table of contents

* [Installation](#Installation)
* [Usage](#Usage)
  * [Full pipeline](#Full-pipeline)
  * [Individual components](#Individual-components)
* [Library](#Library)
* [Connecting the modules](#Connecting-the-modules)
* [File formats](#File-formats)
* [Contributing](#Contributing)

## Installation

To install a stable version:
```sh
pip3 install -U hicstuff
```

or, for the latest development version:

```sh
    pip3 install -e git+https://github.com/koszullab/hicstuff.git@master#egg=hicstuff
```

### External dependencies

Bowtie2 and/or minimap2 as well as samtools are required for the `pipeline` utility.

You can install them via the conda package manager:
```bash
conda install -c bioconda minimap2 bowtie2 samtools
```
Alternatively, on ubuntu you can also install them along with additional dependencies through APT:
```bash
apt-get install samtools bowtie2 minimap2 libbz2-dev liblzma-dev
```

### Docker installation

A pre-built docker image is available on dockerhub and can be ran using:
```bash
docker run koszullab/hicstuff
```

## Usage

The hicstuff command line interface is composed of multiple subcommands. You can always get a summary of all available commands by running:

```bash
hicstuff --help

Simple Hi-C pipeline for generating and manipulating contact matrices.

usage:
    hicstuff [-hv] <command> [<args>...]

options:
    -h, --help                  shows the help
    -v, --version               shows the version

The subcommands are:
    digest          Digest genome into a list of fragments.
    distancelaw     Analyse and plot distance law.
    filter          Filters Hi-C pairs to exclude spurious events.
    iteralign       Iteratively aligns reads to a reference genome.
    pipeline        Hi-C pipeline to generate contact matrix from fastq files.
    rebin           Bin the matrix and regenerate files accordingly.
    subsample       Bootstrap subsampling of contacts from a Hi-C map.
    view            Visualize a Hi-C matrix.  
```

### Full pipeline

All components of the pipeline can be run at once using the `hicstuff pipeline` command. This allows to generate a contact matrix from reads in a single command. By default, the output is in GRAAL compatible COO sparse matrix format, but it can be a 2D bedgraph or cool file instead using the `--matfmt` option. More detailed documentation can be found on the readthedocs website: https://hicstuff.readthedocs.io/en/latest/index.html

    usage:
        pipeline [--quality-min=INT] [--size=INT] [--no-cleanup] [--start-stage=STAGE]
                 [--threads=INT] [--aligner=bowtie2] [--matfmt=FMT] [--prefix=PREFIX]
                 [--tmpdir=DIR] [--iterative] [--outdir=DIR] [--filter] [--enzyme=ENZ]
                 [--plot] [--circular] [--distance-law] [--duplicates] [--read-len=INT]
                 [--centromeres=FILE] [--remove-centromeres=INT] --genome=FILE <input1> [<input2>]

    arguments:
        input1:             Forward fastq file, if start_stage is "fastq", sam
                            file for aligned forward reads if start_stage is
                            "bam", or a .pairs file if start_stage is "pairs".
        input2:             Reverse fastq file, if start_stage is "fastq", sam
                            file for aligned reverse reads if start_stage is
                            "bam", or nothing if start_stage is "pairs".

For example, to run the pipeline with minimap2 using 8 threads and generate a matrix in instagraal format in the directory `out`:

```
hicstuff pipeline -t 8 -a minimap2 -e DpnII -o out/ -g genome.fa reads_for.fq reads_rev.fq
```

If you have already aligned your reads, hicstuff pipeline can also take bam files as input. For example,
to generate a matrix in cool format with a fixed bin size of 5kb:

```
# Note the bam files have to be name-sorted, this can be done using samtools
samtools sort -n aligned_for.bam -o namesorted_for.bam
samtools sort -n aligned_rev.bam -o namesorted_rev.bam
hicstuff pipeline -S bam -e 5000 -M cool -o out/ -g genome.fa namesorted_for.bam namesorted_rev.bam
```


The pipeline can also be run from python, using the `hicstuff.pipeline` submodule. For example, this would run the pipeline with bowtie2 (default) using iterative alignment and keep all intermediate files. For more examples using the API, see the [API demo](https://hicstuff.readthedocs.io/en/latest/notebooks/demo_api.html)

```python
from hicstuff import pipeline as hpi

hpi.full_pipeline(
    'genome.fa', 
    'end1.fq', 
    'end2.fq', 
    no_cleanup=True
    iterative=True
    out_dir='out', 
    enzyme="DpnII")
```

The general steps of the pipeline are as follows:

![Pipeline flowchart](doc/images/pipeline.svg)

### Individual components

For more advanced usage, different scripts can be used independently on the command line to perform individual parts of the pipeline. This readme contains quick descriptions and example usages. To obtain detailed instructions on any subcommand, one can use `hicstuff <subcommand> --help`.

#### Iterative alignment

Truncate reads from a fastq file to 20 basepairs and iteratively extend and re-align the unmapped reads to optimize the proportion of uniquely aligned reads in a 3C library.

    usage:
        iteralign [--aligner=bowtie2] [--threads=1] [--min_len=20]
                  [--tempdir DIR] --out_sam=FILE --genome=FILE <reads.fq>

#### Digestion of the genome

Digests a fasta file into fragments based on a restriction enzyme or a
fixed chunk size. Generates two output files into the target directory
named "info_contigs.txt" and "fragments_list.txt"

    usage:
        digest [--plot] [--figdir=FILE] [--circular] [--size=INT]
               [--outdir=DIR] --enzyme=ENZ <fasta>

 
 For example, to digest the yeast genome with MaeII and HinfI and show histogram of fragment lengths:

`hicstuff digest --plot --outdir output_dir --enzyme MaeII,HinfI Sc_ref.fa`

#### Filtering of 3C events

Filters spurious 3C events such as loops and uncuts from the library based on a minimum distance threshold automatically estimated from the library by default. Can also plot 3C library statistics. This module takes a pairs file with 9 columns as input (readID, chr1, pos1, chr2, pos2, strand1, strand2, frag1, frag2) and filters it. 

    usage:
        filter [--interactive | --thresholds INT-INT] [--plot]
               [--figdir FILE] <input.pairs> <output.pairs>

#### Viewing the contact map

Visualize a Hi-C matrix file as a heatmap of contact frequencies. Allows to tune visualisation by binning and normalizing the matrix, and to save the output image to disk. If no output is specified, the output is displayed interactively. If two contact maps are provided, the log ratio of the first divided by the second will be shown.

    usage:
        view [--binning=1] [--despeckle] [--frags FILE] [--trim INT] [--n-mad FLOAT]
             [--normalize] [--max=99] [--output=IMG] [--cmap=CMAP] [--dpi=INT]
             [--transform=FUN] [--circular] [--region=STR] <contact_map> [<contact_map2>]

    arguments:
        contact_map             Sparse contact matrix in bg2, cool or graal format
        contact_map2            Sparse contact matrix in bg2, cool or graal format,
                                if given, the log ratio of contact_map/contact_map2
                                will be shown.

For example, to view a 1Mb region of chromosome 1 from a full genome Hi-C matrix rebinned at 10kb:

```sh
    hicstuff view --normalize --binning 10kb --region chr1:10,000,000-11,000,000 --frags fragments_list.txt contact_map.tsv
```
### Library

All components of the hicstuff program can be used as python modules. See the documentation on [reathedocs](https://hicstuff.readthedocs.io). The expected contact map format for the library is a simple CSV file, and the objects handled by the library are simple ```numpy``` arrays. The various submodules of hicstuff contain various utilities.

```python
import hicstuff.digest # Functions to work with restriction fragments
import hicstuff.iteralign # Functions related to iterative alignment
import hicstuff.hicstuff # Contains utilities to modify and operate on contact maps as numpy arrays
import hicstuff.filter # Functions for filtering 3C events by type (uncut, loop)
import hicstuff.view # Utilities to visualise contact maps
import hicstuff.io # Reading and writing hicstuff files
import hicstuff.pipeline # Generation and processing of files to generate matrices.
```

### Connecting the modules

All the steps described here are handled automatically when running the `hicstuff pipeline`. But if you want to connect the different modules manually, the intermediate input and output files can be processed using some python scripting.

#### Aligning the reads

You can generate SAM files independently using your favorite read mapping software, use the command line utility `hicstuff iteralign`, or use the helper function `align_reads` in the submodule `hicstuff.pipeline`. For example, to perform iterative alignment using minimap2 (instead of bowtie2 by default):

**Using the python function:**

```python
from hicstuff import pipeline as hpi

hpi.align_reads("end1.fastq", "genome.fasta", "end1.bam", iterative=True, minimap2=True)
```

**Using the command line tool:**

```bash
hicstuff iteralign --minimap2 --iterative -f genome.fasta -o end1.sam end1.fastq
```


#### Extracting contacts from the alignment

The output from `hicstuff iteralign` is a SAM file. In order to retrieve Hi-C pairs, you need to run iteralign separately on the two fastq files and process the resulting alignment files into a name-sorted BAM file as follows using the `pipeline` submodules of hicstuff.

```python
from hicstuff import pipeline as hpi
import pysam as ps
# Sort alignments by read names and get into BAM format
ps.sort("-n", "-O", "BAM", "-o", "end1.bam.sorted", "end1.sam")
ps.sort("-n", "-O", "BAM", "-o", "end2.bam.sorted", "end2.sam")
# Combine BAM files
hpi.bam2pairs("end1.sorted.bam", "end2.sorted.bam", "output.pairs", "info_contigs.txt", min_qual=30)

```
This will generate a "pairs" file containing all read pairs where both reads have been aligned with a mapping quality of at least 30.

#### Attributing each read to a restriction fragment
To build a a contact matrix, we need to attribute each read to a fragment in the genome. This is done under the hood by performing a binary search for each read position against the list of restriction sites in the genome.

```python
from hicstuff import digest as hcd
from Bio import SeqIO

# Build a list of restriction sites for each chromosome
restrict_table = {}
for record in SeqIO.parse("genome.fasta", "fasta"):
    # Get chromosome restriction table
    restrict_table[record.id] = hcd.get_restriction_table(
        record.seq, enzyme, circular=circular
    )

# Add fragment index to pairs (readID, chr1, pos1, chr2,
# pos2, strand1, strand2, frag1, frag2)
hcd.attribute_fragments("output.pairs", "output_indexed.pairs", restrict_table)

```

#### Filtering pairs
The resulting pairs file can then be filtered, either in the command line using the `hicstuff filter` command, or in python using the `hicstuff.filter` submodule. Otherwise, the matrix can be built directly from the unfiltered pairs. 

**Filtering on the command line:**
```bash
hicstuff filter output_indexed.pairs output_filtered.pairs
```
**Filtering in python:**
```python
from hicstuff import filter as hcf

uncut_thr, loop_thr = hcf.get_thresholds("output_indexed.pairs")
hcf.filter_events("output_indexed.pairs", "output_filtered.pairs", uncut_thr, loop_thr)
```
Note that both the command and the python function have various options to generate figure or tweak the filtering thresholds. These options can be displayed using `hicstuff filter -h`

#### Matrix generation
A Hi-C sparse contact matrix can then be generated using the python submodule `hicstuff.pipeline`. The matrix can be generated in GRAAL-compatible COO format, bedgraph2 or cool format.

```python
from hicstuff import pipeline as hpi

n_frags = sum(1 for line in open(fragments_list, "r")) - 1
hpi.pairs2matrix("output_filtered.pairs", "abs_fragments_contacts_weighted.txt", 'fragments_list.txt', mat_fmt="GRAAL")
```

### File formats

* pairs files: This format is used for all intermediate files in the pipeline and is also used by `hicstuff filter`. It is a tab-separated format holding informations about Hi-C pairs. It has an [official specification](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md) defined by the 4D Nucleome data coordination and integration center.
* 2D bedgraph: This is an optional output format of `hicstuff pipeline` for the sparse matrix. It has two fragment per line, and the number of times they are found together. It has the following fields: **chr1, start1, end1, chr2, start2, end2, occurences**
    - Those files can be [loaded by cooler](https://cooler.readthedocs.io/en/latest/cli.html?highlight=load#cooler-load) using `cooler load -f bg2 <chrom.sizes>:<binsize> in.bg2.gz out.cool` where chrom.sizes is a tab delimited file with chromosome names and length on each line, and binsize is the size of bins in the matrix.
* GRAAL sparse matrix: This is a simple tab-separated file with 3 columns: **frag1, frag2, contacts**. The id columns correspond to the absolute id of the restriction fragments (0-indexed). The first row is a header containing the number of rows, number of columns and number of nonzero entries in the matrix. Example:

```
564	564	6978
0	0	3
1	2	4
1	3	3

```

* fragments_list.txt: This tab separated file provides information about restriction fragments positions, size and GC content. Note the coordinates are 0 point basepairs, unlike the pairs format, which has 1 point basepairs. Example:
   - id: 1 based restriction fragment index within chromosome.
   - chrom: Chromosome identifier. Order should be the same as in info_contigs.txt or pairs files.
   - start_pos: 0-based start of fragment, in base pairs.
   - end_pos: 0-based end of fragment, in base pairs.
   - size: Size of fragment, in base pairs.
   - gc_content: Proportion of G and C nucleotide in the fragment.
```
id	chrom	start_pos	end_pos	size	gc_content
1	seq1	0	21	21	0.5238095238095238
2	seq1	21	80	59	0.576271186440678
3	seq1	80	328	248	0.5201612903225806
```

* info_contigs.txt: This tab separated file gives information on contigs, such as number of restriction fragments and size. Example:
   - contig: Chromosome identified. Order should be the same in pairs files or fragments_list.txt.
   - length: Chromosome length, in base pairs.
   - n_frags: Number of restriction fragments in chromosome.
   - cumul_length: Cumulative length of previous chromosome, in base pairs.

```
contig	length	n_frags	cumul_length
seq1	60000	409	0
seq2	20000	155	409
```

### Contributing

All contributions are welcome, in the form of bug reports, suggestions, documentation or pull requests.
We use the [numpy standard](https://numpydoc.readthedocs.io/en/latest/format.html) for docstrings when documenting functions.

The code formatting standard we use is [black](https://github.com/psf/black), with --line-length=79 to follow PEP8 recommendations. We use `pytest` with the `pytest-doctest` and `pytest-pylint` plugins as our testing framework. Ideally, new functions should have associated unit tests, placed in the `tests` folder.
To test the code, you can run:

```bash
pytest --doctest-modules --pylint --pylint-error-types=EF --pylint-rcfile=.pylintrc hicstuff tests
```

### Citation

Please cite hicstuff using the official DOI as follows:

Cyril Matthey-Doret, Lyam Baudry, Amaury Bignaud, Axel Cournac, Remi-Montagne, Nadège Guiglielmoni, Théo Foutel Rodier and Vittore F. Scolari. 2020. hicstuff: Simple library/pipeline to generate and handle Hi-C data . Zenodo. http://doi.org/10.5281/zenodo.4066363

Bibtex entry:

```
@software{cyril_matthey_doret_2020_4066351,
  author       = {Cyril Matthey-Doret and
                  Lyam Baudry and
                  Amaury Bignaud and
                  Axel Cournac and
                  Remi-Montagne and
                  Nadège Guiglielmoni and
                  Théo Foutel-Rodier and
                  Vittore F. Scolari},
  title        = {hicstuff: Simple library/pipeline to generate and handle Hi-C data },
  month        = oct,
  year         = 2020,
  publisher    = {Zenodo},
  version      = {v2.3.1},
  doi          = {10.5281/zenodo.4066351},
  url          = {http://doi.org/10.5281/zenodo.4066363}
}
```

