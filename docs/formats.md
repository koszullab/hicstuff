# File formats

The files hicstuff reads and writes, column by column. [Concepts](concepts.md)
explains what fragments, bins and pairs are, and
[Running the pipeline](guides/pipeline.md#outputs) which file each step
writes.

## Pairs

Every intermediate set of contacts, and the final `valid_idx*.pairs.gz` of a
run, is a pairs file, in the
[format of the 4D Nucleome consortium](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md):
a header of lines starting with `#`, then one tab-separated line per contact.

```text
## pairs format v1.0
#sorted: chr1-chr2-pos1-pos2
#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 frag1 frag2
#chromsize: seq1 60000
#chromsize: seq2 20000
NS500199:33:H3LCGBGXX:2:13208:10265:1714	seq1	22	seq1	280	+	-	1	2
NS500199:33:H3LCGBGXX:1:12309:24210:14793	seq1	26	seq1	104	+	-	1	2
```

| Column | Content |
|---|---|
| `readID` | The name of the read pair |
| `chr1`, `pos1`, `strand1` | The chromosome, position and strand of the mate that comes first in the genome, whether R1 or R2 |
| `chr2`, `pos2`, `strand2` | The same for the other mate |
| `frag1`, `frag2` | The fragment of each mate: its row in `fragments_list.txt`, counted from 0 |

Pairs have 7 columns until each mate is assigned to a fragment, and 9 after.
Positions are 1-based: the aligned 5' end of each mate in mapping mode `5spm`,
the leftmost base of its alignment otherwise (its 3' end on the reverse
strand). The `#chromsize` lines give the length of each chromosome.

The final pairs file of a run is sorted by `chr1`, `chr2`, `pos1`, `pos2`, and
compressed in the BGZF format, which pairix can index. hicstuff's commands read
uncompressed pairs files, except `distancelaw --pairs`, which reads them
gzipped or not.

## `fragments_list.txt`

The fragments of the genome, one row each, contig after contig in the order of
the FASTA file. A run with `--prefix sample` names it `sample.frags.tsv`.

```text
id	chrom	start_pos	end_pos	size	gc_content
1	seq1	0	21	21	0.5238095238095238
2	seq1	21	80	59	0.576271186440678
3	seq1	80	328	248	0.5201612903225806
```

| Column | Content |
|---|---|
| `id` | The number of the fragment on its contig, from 1 |
| `chrom` | The contig |
| `start_pos` | The start of the fragment, 0-based |
| `end_pos` | The end of the fragment, excluded: the fragment above from 0 to 21 holds bases 1 to 21 |
| `size` | `end_pos - start_pos`, in bp |
| `gc_content` | The fraction of G and C bases in the fragment |

## `info_contigs.txt`

The contigs, one row each. A run with `--prefix sample` names it
`sample.chr.tsv`.

```text
contig	length	n_frags	cumul_length
seq1	60000	409	0
seq2	20000	155	409
```

| Column | Content |
|---|---|
| `contig` | The contig, as in the FASTA file |
| `length` | Its length, in bp |
| `n_frags` | Its number of fragments |
| `cumul_length` | The number of fragments on the contigs before it: the row of its first fragment in `fragments_list.txt`, counted from 0 |

## cool and mcool

The default contact map is a [cooler](https://github.com/open2c/cooler): an
HDF5 file of its bins, its chromosomes and its non-zero pixels, which cooler,
cooltools, HiGlass and most Hi-C tools read.

- Without `--binning`, the pipeline writes a `.cool` with one bin per fragment
  (`bin-type: variable`), whose bins table has the `chrom`, `start` and `end`
  columns.
- With `--binning N`, it writes a `.mcool` holding the map at N, 2N, 5N, 10N,
  20N, 50N and 100N bp, under `/resolutions/<size>`, each with the `weight`
  column of `cooler balance`. With `--zoomify false`, it writes the map at N bp
  alone, as a balanced `.cool`.
- With `--binning bins.bed`, a BED file of bins, it writes the map of those
  bins as a balanced `.cool`, which is not zoomified.

Both store the upper triangle of the map, and raw counts. `convert`, `rebin`
and `subsample` write cool files with the `size` and `gc_content` columns as
well, when they know them, and record hicstuff's version in the metadata.

## metacool

From 50,000 contigs in the map, as a metagenome assembly may have, or with
`--matrix metacool`, the pipeline writes a `.metacool` file, made by
[metacooler](https://pypi.org/project/metacooler/)'s `create pairs`. It is a
cool file with one bin per fragment, as the default cool file without
`--binning`, or the bins of `--binning` (a bin size or a BED file of bins),
which cooler reads, with metacooler's tables of genomes and contigs on top of
it. It is neither zoomified nor balanced.

- Its contigs are not yet assigned to genomes: they are all in metacooler's
  `__unbinned__` genome, longest first, until `metacooler reassign` assigns
  them, from a binning such as metaTOR's.
- The contig of each bin is stored as an integer, where cooler tries an HDF5
  enum first, which is slow to write with many contigs (see
  [The contact map](guides/pipeline.md#the-contact-map)).
- Its metadata hold the command line of the run (`commands`) and metacooler's
  counts of the contacts within and between contigs (`qc`), which
  `metacooler metadata get -k qc` prints. They are counted on the pairs, read
  by read, and split the contacts within a contig at 1001 bp: those whose
  reads are more than 1 kb apart are long-range (`intra_contig_long`, and per
  contig the `qc_intra_contig_long` column of the contigs table), the contacts
  metaTOR's `mge` reads from the map. hicstuff sets no threshold: it is
  metacooler's default for counts on the pairs, 1001 bp from metacooler 0.4.0.
  `metacooler qc` does not read these counts: it counts the contacts again,
  from the map's pixels, between the midpoints of the fragments, and
  `metacooler qc --store` would replace the counts on the pairs with those.

## Graal matrices

The format of [instaGRAAL](https://github.com/koszullab/instaGRAAL): a sparse
matrix of fragments, written as `abs_fragments_contacts_weighted.txt` (or
`<prefix>.mat.tsv`), next to `fragments_list.txt` and `info_contigs.txt`.

```text
564	564	1030
1	2	2
2	2	17
2	3	18
```

The first line holds the number of rows, of columns and of non-zero entries.
Each other line holds two fragments, as rows of `fragments_list.txt` counted
from 0, and the number of contacts between them. The matrix carries no
coordinates: commands that read one need its `fragments_list.txt` (`--frags`),
and for some its `info_contigs.txt` (`--chroms`).

## 2D bedgraph

One line per pair of bins with contacts, tab-separated, without a header:

```text
seq1	21	80	seq1	80	328	2
seq1	80	328	seq1	80	328	17
seq1	80	328	seq1	328	473	18
```

The columns are `chrom1`, `start1`, `end1`, `chrom2`, `start2`, `end2` and the
number of contacts. cooler reads it, given its bins. Those of the pipeline's
2D bedgraph are the fragments, which `fragments_list.txt` lists:

```bash
tail -n +2 fragments_list.txt | cut -f 2-4 > fragments.bed
cooler load -f bg2 fragments.bed abs_fragments_contacts_weighted.bg2 map.cool
```

## Distance law

The table of `pipeline --distance-law` and `distancelaw --outputfile-tabl`:

```text
## hicstuff v3.3.0
## distance_law
## columns: start_bp	p(s)	chrom
1	1e-09	seq1
2	1e-09	seq1
```

One line per distance bin, of each chromosome, arm or sample: the distance at
which the bin starts, in bp, the value of the distance law, and the name of
the chromosome, arm or sample. The pipeline writes the values before they are
normalised; `distancelaw` writes them normalised, and averaged with
`--average`.

## Centromeres

The file of `--centromeres`: one line, with the position of the centromere of
each chromosome, in bp, separated by spaces, in the order of the chromosomes
of the genome.

```text
150000 238000 114000
```

## Log

`hicstuff_<date>.log`, in the output directory of a run: a header, then one
line per message, with its time and level.

```text
## hicstuff: v3.3.0 log file
## date: 2026-11-17 23:22:37
## cmd: hicstuff pipeline --mapping 5spm --enzyme DpnII --genome genome.fa --outdir out/ R1.fq.gz R2.fq.gz
## enzyme: DpnII
## input1: R1.fq.gz
## input2: R2.fq.gz
## ref: genome.fa
## aligner: bwa mem -5SPM in-process, pybwa 2.3.2 (bwa 0.7.18)
---
2026-11-17,23:22:37 :: INFO :: Checking content of fastq files.
2026-11-17,23:22:37 :: INFO :: 10000 reads found in each fastq file.
```

The `cmd` line is the command that ran the pipeline, quoted to paste in a
shell. A run started from Python, by `full_pipeline`, has none unless its
`cmd` argument gives one. With several pairs of FASTQ files, `input1` and
`input2` list the R1 and R2 files, comma-separated. The `aligner` line is written when bwa runs
in-process, or in mapping mode `5spm`. `hicstuff stats` reads the statistics
of the run from the log.
