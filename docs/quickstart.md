# Quickstart

This page takes two FASTQ files to a contact map, a binned multi-resolution
`.mcool` and a picture of it, with nothing installed but hicstuff. It uses the
small test data set of hicstuff's repository: 10,000 read pairs on a genome of
two sequences, 80 kb in all.

## Input files

| File | What it is |
|---|---|
| `seq.fa` | The reference genome, in FASTA format (it may be gzipped) |
| `sample.reads_for.fastq.gz` | The first mate of each read pair (R1) |
| `sample.reads_rev.fastq.gz` | The second mate, in the same order (R2) |

```bash
base=https://raw.githubusercontent.com/koszullab/hicstuff/master/test_data
curl -LO "$base/genome/seq.fa"
curl -LO "$base/sample.reads_for.fastq.gz"
curl -LO "$base/sample.reads_rev.fastq.gz"
```

The library was digested with DpnII. The enzyme matters: hicstuff cuts the
genome into the fragments the enzyme made, and assigns each read to one of
them.

## Step 1: run the pipeline

```bash
hicstuff pipeline --mapping 5spm --enzyme DpnII --genome seq.fa --outdir out \
    sample.reads_for.fastq.gz sample.reads_rev.fastq.gz
```

`--mapping 5spm` is the recommended mapping mode. It aligns the mates with
`bwa mem -5SPM`, as the 4DN, Juicer and Arima pipelines do, and records each
mate at its aligned 5' end. bwa runs in-process: no aligner needs to be
installed. [Aligning reads](guides/alignment.md#5spm) explains the mode. The
run takes a few seconds, and logs each step (excerpt):

```console
INFO :: Checking content of fastq files.
INFO :: 10000 reads found in each fastq file.
INFO :: Aligning with bwa mem -5SPM in-process, pybwa 2.3.2 (bwa 0.7.18).
INFO :: Now mapping R1 and R2 reads with bwa mem -5SPM...
INFO :: Genome digested into 564 fragments with a median length of 87.5
INFO :: 89.00% reads (single ends) mapped with Q >= 30 (17876/20000)
INFO :: 7891 pairs successfully mapped (78.91%)
INFO :: Generating matrix from pairs file out/tmp/valid_idx.pairs (7891 pairs in the file)
INFO :: Contact map generated after 0h 0m 2s
```

## Step 2: look at what it wrote

```console
$ ls out
abs_fragments_contacts_weighted.cool  hicstuff_20261117232237.log  tmp
fragments_list.txt                    info_contigs.txt             valid_idx.pairs.gz
```

| File | What it holds |
|---|---|
| `abs_fragments_contacts_weighted.cool` | The contact map, one bin per DpnII fragment |
| `fragments_list.txt` | The fragments: contig, start, end, size, GC content |
| `info_contigs.txt` | The contigs: length, number of fragments |
| `valid_idx.pairs.gz` | Every read pair kept, with the fragment of each mate |
| `hicstuff_<date>.log` | The log of the run |

The [file formats](formats.md) page describes each of them. The `.cool` file
is a [cooler](https://github.com/open2c/cooler), which cooler, cooltools,
HiGlass and the rest of the ecosystem read:

```bash
cooler info out/abs_fragments_contacts_weighted.cool
```

## Step 3: read the statistics

```bash
hicstuff stats out/hicstuff_*.log
```

```console
Sample:             hicstuff_20261117232237.log
Total read pairs:   10,000
Mapped reads:       17,876 (89.38%  of all reads)
Unmapped reads:     2,124 (10.62%  of all reads)
Recovered contacts: 7,891 (78.91%  of all read pairs)
Removed contacts:   0 (0.0% of all contacts)
  Small contigs:    0 (0.0% of all contacts)
  Filtered out:     0 (0.0% of all contacts)
    Loops:          0 (0.0% of all contacts)
    Uncuts:         0 (0.0% of all contacts)
    Weirds:         0 (0.0% of all contacts)
  PCR duplicates:   0 (0.0% of all contacts)
Final contacts:     7,891 (100.0% of all contacts)
```

Nothing was removed: the filters are off unless asked for.

## Step 4: bin, filter and balance

Most analyses want bins of a fixed size rather than restriction fragments.
`--binning` bins the map, and the pipeline then writes a multi-resolution
`.mcool`, balanced at every resolution. `--filter` removes the uninformative
3C events and `--duplicates` the PCR duplicates; `--prefix` names the files
after the sample.

```bash
hicstuff pipeline --mapping 5spm --enzyme DpnII --genome seq.fa --outdir out \
    --binning 1000 --filter --duplicates --prefix sample \
    sample.reads_for.fastq.gz sample.reads_rev.fastq.gz
```

```console
$ ls out/sample.*
out/sample.chr.tsv    out/sample.hicstuff_20261117232252.log  out/sample.valid_idx_pcrfree.pairs.gz
out/sample.frags.tsv  out/sample.mcool
```

`sample.mcool` holds the map at 1, 2, 5, 10, 20, 50 and 100 kb, each with a
`weight` column from `cooler balance`. [Running the pipeline](guides/pipeline.md)
lists the outputs of every option.

!!! note "This test library is not a real one"

    With `--filter`, nearly all its pairs are removed as uncut fragments: the
    test reads were made to exercise the code, not to look like a library.
    [Filtering 3C events](guides/filter.md) explains what the filter looks
    for, and how to check its thresholds on your own data.

## Step 5: view the map

```bash
hicstuff view --binning 5kb --normalize --output map.png \
    out/abs_fragments_contacts_weighted.cool
```

`view` bins the fragment-level map to 5 kb, balances it, and saves the
heatmap. Without `--output`, it opens a window. A resolution of the `.mcool`
is read by its URI:

```bash
hicstuff view --region seq1:10kb-40kb --output zoom.png \
    "out/sample.mcool::/resolutions/1000"
```

## The same from Python

```python
import hicstuff.io as hio
import hicstuff.hicstuff as hcs
import hicstuff.pipeline as hpi
import hicstuff.view as hcv

hpi.full_pipeline(
    genome="seq.fa",
    input1="sample.reads_for.fastq.gz",
    input2="sample.reads_rev.fastq.gz",
    mapping="5spm",
    enzyme="DpnII",
    out_dir="out",
    force=True,
)

# Any of the three formats: cool, 2D bedgraph or graal
cool = "out/abs_fragments_contacts_weighted.cool"
mat, frags, chroms = hio.flexible_hic_loader(cool)
mat = hcs.normalize_sparse(mat, norm="ICE")
hcv.plot_matrix(hcv.sparse_to_dense(mat), filename="map.png")
```

[Using the Python API](guides/python.md) goes further.

## Next steps

- [Concepts](concepts.md): fragments, bins, pairs and 3C events.
- [Running the pipeline](guides/pipeline.md): every input, output and option.
- [Aligning reads](guides/alignment.md): the four mapping modes, `5spm` first,
  and the aligners.
- [CLI reference](cli.md): every command.
