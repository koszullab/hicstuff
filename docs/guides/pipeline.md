# Running the pipeline

`hicstuff pipeline` takes a Hi-C library from its reads, or from aligned reads
or pairs, to a contact map. This page lists its inputs, the files it writes,
and the options of each step. The [CLI reference](../cli.md#pipeline) lists
every option in one table, and [Concepts](../concepts.md#the-pipeline) draws
the steps.

```bash
hicstuff pipeline --mapping 5spm --enzyme DpnII --genome genome.fa --outdir out/ \
    R1.fq.gz R2.fq.gz
```

`--mapping 5spm` is the recommended mapping mode: see [Alignment](#alignment).

## Inputs

`--start-stage` (`-S`) says what the positional inputs are:

| `--start-stage` | Inputs | What the run does |
|---|---|---|
| `fastq` (default) | The R1 and R2 reads, in two FASTQ files, gzipped or not | Every step |
| `bam` | Two BAM files, one per mate, each aligned as single-end reads | From pairing the mates on |
| `pairs` | One [pairs file](../formats.md#pairs) of 7 columns, uncompressed | From assigning fragments on |
| `pairs_idx` | One pairs file of 9 columns, with `frag1` and `frag2`, uncompressed | From filtering on |

The BAM files may be in any order: files sorted by read name (with
`samtools sort -n` or `-N`) are read as they are, others are sorted by name
first, in the temporary directory.

### Several pairs of FASTQ files

The lanes or sequencing runs of a library are read as one library: give the
R1 files, then the R2 files. Each R1 file goes with the R2 file in the same
position.

```bash
# Glob patterns, which the shell expands
hicstuff pipeline -m 5spm -e DpnII -g genome.fa -o out/ lib1_*_R1.fq.gz lib1_*_R2.fq.gz

# The same patterns, quoted: hicstuff expands them
hicstuff pipeline -m 5spm -e DpnII -g genome.fa -o out/ 'lib1_*_R1.fq.gz' 'lib1_*_R2.fq.gz'

# Comma-separated lists
hicstuff pipeline -m 5spm -e DpnII -g genome.fa -o out/ \
    lib1_L001_R1.fq.gz,lib1_L002_R1.fq.gz lib1_L001_R2.fq.gz,lib1_L002_R2.fq.gz
```

- **Patterns** (`*`, `?`, `[...]`) match the files sorted by name, as the
  shell does: the R1 and R2 files of a lane, named alike, come in the same
  position. An argument that is the path of a file is that file, even if its
  name holds a comma or a pattern character.
- **More than two arguments**, as the shell gives them for unquoted patterns,
  are halves: the first half are the R1 files, the second half the R2 files.
  An odd number of them stops the command.
- **Each pair is checked** before the run starts: every file must exist and be
  given once, and the R1 and R2 files of a pair must start with the same read,
  named alike but for a final `/1` or `/2`. Files given in different orders
  stop the run. The reads of each pair of files are counted, unless
  `--skip-count`, and the log names a pair whose files have different numbers
  of reads.
- **The R1 files are concatenated** in the temporary directory, as are the R2
  files, and the run goes on from that pair of files: its outputs are those of
  the library in one pair of files. PCR duplicates (`--duplicates`) are found
  across the lanes. Gzip files are copied as they are; plain files are
  compressed when a gzip file is among them. The concatenated files take as
  much space as the inputs, and are removed with the other temporary files.

Several pairs of files are read from FASTQ files only: `--start-stage bam`
takes two BAM files, and `pairs` and `pairs_idx` one pairs file. Separate
libraries are better run one by one, since `--duplicates` would take a pair
found in two libraries for a PCR duplicate; their maps, on the same fragments
or bins, can then be summed with `cooler merge`.

!!! warning "`--start-stage pairs_idx` sorts its input in place"

    The pairs file given with `-S pairs_idx` is sorted where it is, and must
    be uncompressed: given a gzipped file, such as the `.pairs.gz` a run
    writes, the run fails after replacing it with uncompressed text. Give it a
    decompressed copy:

    ```bash
    zcat out/valid_idx.pairs.gz > valid_idx.pairs
    hicstuff pipeline -S pairs_idx -e DpnII -g genome.fa -o out2/ valid_idx.pairs
    ```

`--genome` (`-g`) is the reference, needed at every stage, since the
fragments are cut from it:

- A FASTA file, gzipped or not. bowtie2 and bwa need an index: one named after
  the FASTA (`genome.fa.1.bt2`, as `bowtie2-build genome.fa genome.fa` writes
  it, or `genome.fa.bwt`, as `bwa index genome.fa` does) is used; otherwise the
  pipeline copies the FASTA and builds one in a directory of its own in the
  temporary directory, at every run, which the cleanup removes.
- For bowtie2 or bwa, the prefix of an index without its FASTA (`genome` for
  `genome.1.bt2` or `genome.bwt`, …). The sequence is then read from the
  index: by `bowtie2-inspect` for a bowtie2 index, and by hicstuff itself,
  from the `.pac` and `.amb` files, for a bwa index.

Any other `--genome` stops the run before it writes anything.

`--enzyme` (`-e`) is the enzyme the library was digested with, several
enzymes separated by commas, or a fragment size in base pairs (5000 by
default); see [Digesting a genome](digest.md). `--circular` (`-C`) treats every
chromosome as circular, and `--size` (`-s`) leaves out the contigs shorter
than that many base pairs, with the pairs that have a mate on them.

## Outputs

The pipeline writes into `--outdir` (`-o`, the current directory by
default). `--prefix` (`-P`) names the files after a sample instead of using
instaGRAAL's file names:

| File | Without `--prefix` | With `--prefix sample` | Written |
|---|---|---|---|
| Log | `hicstuff_<date>.log` | `sample.hicstuff_<date>.log` | always |
| Fragments | `fragments_list.txt` | `sample.frags.tsv` | always |
| Contigs | `info_contigs.txt` | `sample.chr.tsv` | always |
| Contact map, cool | `abs_fragments_contacts_weighted.cool` | `sample.cool` | by default, below 50,000 contigs |
| Contact map, metacool | `abs_fragments_contacts_weighted.metacool` | `sample.metacool` | by default from 50,000 contigs, or with `--matrix metacool` |
| Contact map, binned | `abs_fragments_contacts_weighted.mcool` | `sample.mcool` | with `--binning` |
| Contact map, graal | `abs_fragments_contacts_weighted.txt` | `sample.mat.tsv` | with `--matrix graal` |
| Contact map, 2D bedgraph | `abs_fragments_contacts_weighted.bg2` | `sample.bg2` | with `--matrix bg2` |
| Pairs | `valid_idx.pairs.gz` | `sample.valid_idx.pairs.gz` | always |
| Distance law | `distance_law.txt` | `sample.distance_law.txt` | with `--distance-law` |
| Plots | `plots/*.pdf` | `plots/sample_*.pdf` | with `--plot` |

`<date>` is the time the run started, as `YYYYmmddHHMMSS`, so a log never
replaces another. The pairs file holds the pairs the map is built from: it is
`valid_idx_filtered.pairs.gz` after `--filter`, and
`valid_idx_pcrfree.pairs.gz` after `--duplicates`. It is sorted by `chr1`,
`chr2`, `pos1`, `pos2` and compressed in the BGZF format, which pairix can
index. With `--matrix none`, the run writes no contact map, and ends with this
file.

A run refuses to replace an output that exists, or a temporary file left by an
earlier run, unless given `--force` (`-F`).

### Temporary files

The intermediate files go to `--tmpdir` (`-T`), by default `tmp/` in the
output directory: the alignments (`for.bam`, `rev.bam`), the pairs before and
after each step (`valid.pairs`, `valid_idx.pairs`, …) and, in a directory of
their own, the index built for the run when the genome has none, or the
sequence read from an index without FASTA (a FASTA is read where it is, and
never copied). These files are deleted at the end of the run, unless
`--no-cleanup` (`-n`) keeps them. The `tmp/` directory stays, with what the
cleanup leaves in it: in mapping mode `cutsite`, the cut reads
(`digest_for.fq.gz`, `digest_rev.fq.gz`). To run
several samples in one output directory, give each its `--prefix`: the
temporary files take it too.

### Plots

`--plot` (`-p`) saves a PDF of each step that has one, in `plots/`:

| Plot | Shows | Written |
|---|---|---|
| `frags_hist.pdf` | The distribution of fragment lengths | from FASTQ or BAM files |
| `event_distance.pdf` | The 3C events by distance, with the filter's thresholds | with `--filter` |
| `event_distribution.pdf` | The proportion of each kind of event | with `--filter` |
| `distance_law.pdf` | The distance law and its slope | with `--distance-law` |

## Alignment

`--mapping` (`-m`) chooses how the reads are aligned. **`5spm` is the
recommended mode**: it aligns the mates with `bwa mem -5SPM`, as the 4DN,
Juicer and Arima pipelines do, in-process unless `--aligner-backend binary`
asks for the `bwa` program, and records each mate at its aligned 5' end. It is
not the default, `normal`. The other modes, `normal`, `iterative` and
`cutsite`, align with `--aligner` (`-a`): `bowtie2` (the default), `bwa` or
`minimap2`. `--threads` (`-t`) sets the threads of the aligner and of the other
steps. Each mate is aligned on its own; supplementary alignments are dropped,
and a pair is kept when both mates aligned with a mapping quality of at least
`--quality-min` (`-q`, 30 by default). [Aligning reads](alignment.md) explains
the aligners and the modes.

Before aligning, the pipeline counts the reads of both FASTQ files, and logs an
error if their numbers differ. `--skip-count` skips the count, which takes a
while on large files; the statistics of the run then count the read pairs as
half the reads aligned.

## Filtering

| Option | Removes |
|---|---|
| `--quality-min N` (`-q`) | Pairs with a mate aligned with a mapping quality below N (30 by default) |
| `--size N` (`-s`) | Contigs shorter than N bp, and the pairs with a mate on them |
| `--filter` (`-f`) | Uncuts, loops and impossible pairs, with thresholds estimated from the library |
| `--duplicates` (`-D`) | PCR duplicates: one pair is kept among those with the same coordinates |
| `--exclude chrM,2u` (`-E`) | These chromosomes, from the cool and metacool files only |

`--filter` needs a restriction enzyme: with a fragment size, no fragment is a
restriction fragment, and the run refuses it (though not `mnase` or `dnase`,
which it then filters as 100 bp fragments). [Filtering 3C events](filter.md)
explains how the thresholds are chosen, and how to check them.

## The contact map

`--matrix` (`-M`) chooses its format:

- **`cool`** (the default, below 50,000 contigs): a
  [cooler](https://github.com/open2c/cooler) with one bin per fragment, made by
  `cooler cload pairs`.
- **`metacool`** (the default from 50,000 contigs, but with `--binning`): a
  [metacooler](https://pypi.org/project/metacooler/) file with one bin per
  fragment, made by `metacooler create pairs`. It is a cool file, which cooler
  reads, with metacooler's tables of genomes and contigs. Its contigs are not
  yet assigned to genomes: `metacooler reassign` assigns them, from a binning
  such as metaTOR's. See [metacool](../formats.md#metacool). `--binning` and
  `--balancing-args`, options of the cool map, are refused with it, but for
  `--binning 0` and blank values; `--exclude` leaves contigs out of it, as out
  of the cool map.
- **`graal`**: instaGRAAL's sparse matrix, one row and column per fragment;
  with `fragments_list.txt` and `info_contigs.txt`, it is what instaGRAAL
  reads.
- **`bg2`**: a 2D bedgraph, one line per pair of fragments with contacts.
- **`none`**: no map. The run ends with its pairs file, for the tools that
  build their own contact maps from pairs, such as
  [metaTOR](https://github.com/koszullab/metaTOR) and
  [metacooler](https://pypi.org/project/metacooler/). `--binning`,
  `--balancing-args` and `--exclude`, options of the cool map, are refused with
  it, but for `--binning 0` and blank values.

cooler tries to store the contig of each bin of a cool file as an HDF5 enum, of
one member per contig. Building it takes a time that grows with the square of
the number of contigs (about 2 s at 32,000 contigs, 30 s at 117,000, four times
more for twice as many); HDF5 then refuses it past a few thousand members, and
cooler stores integers instead. It pays that time for every cooler it writes:
`cooler cload pairs` writes one for every 15 million pairs, then the map.
metacooler stores integers from the start. So with 50,000 contigs or more in
the map (those that `--size` and `--exclude` keep), as a metagenome assembly
may have, the default map is a metacool file, and the run warns of it before
it aligns reads or builds an index. With `--binning`, which the metacool map
does not take, the default map stays a cool file, as with `--matrix cool`: the
run then warns that it will be slow to write. `--matrix none` ends the run with
its pairs file, and `bg2` and `graal`, text files, do not pay that time
either.

### Binned, multi-resolution and balanced maps

`--binning N` (`-b`) counts the contacts in bins of N base pairs instead of
fragments, and applies to cool files only:

1. The pairs are binned at N bp.
2. The map is zoomified into a `.mcool` file holding N, 2N, 5N, 10N, 20N, 50N
   and 100N, and the single-resolution `.cool` is deleted (kept with
   `--no-cleanup`). `--zoomify false` (`-z false`) keeps the `.cool` alone.
3. Every resolution is balanced with `cooler balance`, which stores its
   `weight` column. It runs on `--threads` processes. `--balancing-args`
   (`-B`) passes more options, checked before the run starts:
   `-B "--max-iters 500 --ignore-diags 2"`. A resolution that fails to
   balance is logged as a warning, and the run goes on.

A resolution is read by its URI, `out/sample.mcool::/resolutions/5000`, with
hicstuff's commands as with cooler's.

```bash
# An Arima library at 1 kb and up, filtered, balanced
hicstuff pipeline --mapping 5spm --enzyme DpnII,HinfI --binning 1000 --threads 8 \
    --filter --duplicates --prefix sample --genome genome.fa --outdir out/ \
    R1.fq.gz R2.fq.gz
```

## The distance law

`--distance-law` (`-d`) writes `distance_law.txt`: the frequency of contacts
against the distance between the mates, for each chromosome, from the pairs
assigned to fragments, before `--filter` and `--duplicates`. With
`--centromeres FILE` (`-c`), it is computed for each chromosome arm, and
`--remove-centromeres N` (`-r`) leaves out N bp on each side of the
centromeres. The table is neither normalised nor averaged: `hicstuff
distancelaw` does both, and plots it. See [The distance law](distance-law.md).

## Starting again from intermediate files

```bash
# Reads aligned elsewhere with bwa mem -5SPM, one BAM file per mate
hicstuff pipeline -S bam -m 5spm -e DpnII -g genome.fa -o out/ R1.bam R2.bam

# A pairs file from another tool, uncompressed
hicstuff pipeline -S pairs -e DpnII -g genome.fa -o out/ sample.pairs

# The pairs of an earlier run, to bin them otherwise
zcat out/valid_idx.pairs.gz > valid_idx.pairs
hicstuff pipeline -S pairs_idx -e DpnII -g genome.fa -b 5000 -o out_5kb/ \
    valid_idx.pairs
```

The enzyme and genome must be those of the first run: they define the
fragments the pairs are assigned to. From BAM files, `-m 5spm` records each
mate at its aligned 5' end, as when the pipeline aligns the reads in that mode;
the other modes record the leftmost base of each alignment. Starting from pairs, the log holds no
alignment statistics, and `hicstuff stats` cannot read it.

## Reading the log

The log starts with the version, the date, the enzyme, the inputs and the
reference, then records every step: the number of reads, the proportion
aligned, the pairs removed by each filter, and how long the run took.
`hicstuff stats` summarises it:

```bash
hicstuff stats out/hicstuff_20261117221914.log
```

```console
Sample:             hicstuff_20261117221914.log
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

| Line | Counts |
|---|---|
| Total read pairs | The read pairs of the FASTQ files |
| Mapped reads, Unmapped reads | The mates aligned with a mapping quality of at least `--quality-min`, and the others: twice the read pairs in all, or more in `cutsite` mode |
| Recovered contacts | The pairs with both mates aligned |
| Small contigs | The pairs left out by `--size` |
| Filtered out | The pairs removed by `--filter`: loops, uncuts and weirds |
| PCR duplicates | The pairs removed by `--duplicates` |
| Final contacts | The pairs in the contact map |
