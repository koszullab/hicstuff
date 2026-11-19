# Why hicstuff?

A Hi-C experiment ends as two FASTQ files. Between them and a contact map lie
a reference genome cut into restriction fragments, two aligned mates per read
pair, the 3C events that carry no information (uncut fragments, self-ligated
loops, PCR duplicates) and a choice of resolution. hicstuff runs all of it with
one command, for any organism, and keeps each step usable on its own.

It started as the merge of the [yahcp](https://github.com/baudrly/yahcp)
pipeline, the original [hicstuff](https://github.com/baudrly/hicstuff)
library, the [3C tutorial](https://github.com/axelcournac/3C_tutorial) and the
[DADE pipeline](https://github.com/scovit/dade), and is developed by the
[Koszul laboratory](https://github.com/koszullab) at Institut Pasteur.

## What it does

- **One command, from reads to matrix.** `hicstuff pipeline` aligns the reads,
  pairs the mates, assigns them to restriction fragments, filters them and
  writes a [cool](https://github.com/open2c/cooler) file, or a binned,
  balanced, multi-resolution `.mcool`, or for a metagenome a
  [metacool](https://pypi.org/project/metacooler/) file. It can also start
  from BAM or pairs files.
- **Nothing to install but hicstuff.** `pip install hicstuff` runs the whole
  pipeline in mapping mode `5spm`, the recommended one: bwa runs in-process,
  through [pybwa](https://github.com/fulcrumgenomics/pybwa). bowtie2 and
  minimap2 are there too, for the other modes, when they are installed.
- **Any genome, any protocol.** Restriction enzymes by name, several at once
  (`DpnII,HinfI` for Arima kits), fixed-size chunks for MNase or DNase
  libraries, circular chromosomes, and a minimum contig size for draft
  assemblies.
- **Restriction-fragment resolution.** The unbinned map has one bin per
  restriction fragment, the resolution
  [instaGRAAL](https://github.com/koszullab/instaGRAAL) scaffolds from, and is
  also written in instaGRAAL's own `graal` format.
- **A look at the library itself.** The 3C events each pair comes from
  (loops, uncuts, abnormal pairs), PCR duplicates, the distance law P(s) and a
  summary of the run's mapping and filtering statistics.
- **Four ways to align.** `bwa mem -5SPM`, as the 4DN, Juicer and Arima
  pipelines run it: mapping mode `5spm`, the recommended one. Then plain
  alignment, iterative alignment of truncated reads, and reads cut at their
  ligation sites.
- **Tools for the maps it writes.** `view`, `rebin`, `subsample`, `convert`
  and `scalogram` take cool, 2D bedgraph or graal matrices; `distancelaw`
  works from the pairs, and `missview` from a genome alone.
- **A Python library.** Every command is a function, and `hicstuff.hicstuff`
  holds the matrix operations they use: normalisation, trimming, despeckling,
  binning, subsampling, compartments.

## A quick look

```bash
# Reads to a balanced multi-resolution .mcool, at 1 kb and up, for an Arima library
hicstuff pipeline --mapping 5spm --enzyme DpnII,HinfI --binning 1000 --threads 8 \
    --filter --duplicates --genome genome.fa --outdir out/ R1.fq.gz R2.fq.gz

# What happened to the reads
hicstuff stats out/hicstuff_*.log

# A 2 Mb window of chromosome 1 at 10 kb
hicstuff view --region chr1:1Mb-3Mb --output chr1.png \
    out/abs_fragments_contacts_weighted.mcool::/resolutions/10000
```

The [Quickstart](quickstart.md) runs the same steps on the test data that ships
with hicstuff, and [Concepts](concepts.md) explains what fragments, bins, pairs
and 3C events are.
