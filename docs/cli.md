# CLI reference

The `hicstuff` command has one subcommand per task. This page documents every
command and option; `-h/--help` prints the same in a terminal.

```console
hicstuff --help
hicstuff <command> --help
hicstuff --version
```

## Conventions

| Convention | Meaning |
|---|---|
| Contact maps | Every command that reads a map takes a cool file, a resolution of a `.mcool` (`file.mcool::/resolutions/<N>`), a 2D bedgraph or a graal matrix, recognised from the file. A graal matrix needs its `fragments_list.txt` with `--frags` for whatever uses coordinates. |
| Sizes | The bin sizes of `view`, `rebin` and `missview`, and the regions and distances of `view` and `scalogram`, take a unit: `500bp`, `10kb`, `2Mb`, in any case. Other sizes are whole numbers of bp. |
| Regions | `chr1`, or `chr1:1,000,000-2,000,000`, or `chr1:1Mb-2Mb`. |
| Output prefixes | `rebin`, `subsample` and `convert` take an output *prefix*, and add the extension of the format. |
| `-F, --force` | Replace an output that exists. In `cutsite`, `-F` is `--forward`. |
| `HICSTUFF_ALIGNER_BACKEND` | How bwa runs when `--aligner-backend` is not given: `auto`, `binary` or `python`. |
| `HICSTUFF_SORT=python` | Sort the pairs in Python even when a `sort` program is installed. |

New to the tool? Start with the [Quickstart](quickstart.md); for fragments,
bins, pairs and 3C events, see [Concepts](concepts.md).

## Run the pipeline

### `pipeline` { #pipeline }

Run the whole pipeline, from the reads (or BAM or pairs files) to a contact
map. **API:** `hicstuff.pipeline.full_pipeline`. **Guide:**
[Running the pipeline](guides/pipeline.md).

```console
hicstuff pipeline [OPTIONS] INPUT1... [INPUT2...]
```

`INPUT1` and `INPUT2` are the R1 and R2 FASTQ files, or the two BAM files of
`--start-stage bam`. `--start-stage pairs` and `pairs_idx` take one pairs file.
Several pairs of FASTQ files, such as the lanes of a library, are read as one
library: give the R1 files, then the R2 files, comma-separated or as glob
patterns. See [Several pairs of FASTQ files](guides/pipeline.md#several-pairs-of-fastq-files).

**Recommended:** `--mapping 5spm`, which aligns the mates with
`bwa mem -5SPM`, in-process, and records each mate at its aligned 5' end. See
[Mapping modes](guides/alignment.md#mapping-modes).

| Option | Default | Description |
|---|---|---|
| `-g, --genome` | *required* | The reference: a FASTA file, or for bowtie2 and bwa the prefix of its index, with or without the FASTA. |
| `-S, --start-stage` | `fastq` | What the inputs are: `fastq`, `bam`, `pairs` (7 columns, uncompressed) or `pairs_idx` (9 columns, with fragments, uncompressed). |
| `-e, --enzyme` | `5000` | Restriction enzyme, or several separated by commas (`DpnII,HinfI`); a fragment size in bp; or `mnase` or `dnase`, for 100 bp fragments. |
| `-a, --aligner` | `bowtie2` | `bowtie2`, `bwa` or `minimap2`, for the mapping modes other than `5spm`, which aligns with `bwa`. |
| `--aligner-backend` | `auto` | How bwa runs: `binary` (the bwa program), `python` (in-process, through pybwa) or `auto` (the program if it is on `PATH`, pybwa otherwise, and pybwa first in mode `5spm`). Default: `$HICSTUFF_ALIGNER_BACKEND`, else `auto`. |
| `-m, --mapping` | `normal` | Mapping mode: `5spm`, the recommended one, `normal`, `iterative` or `cutsite`. See [Mapping modes](guides/alignment.md#mapping-modes). |
| `-q, --quality-min` | `30` | Minimum mapping quality of both mates of a pair. |
| `-R, --read-len` | the first read's | Length of the longest reads, for mapping mode `iterative`. |
| `--skip-count` | off | Do not count the reads of the FASTQ files before aligning them. |
| `-C, --circular` | off | The chromosomes are circular, for the distance law. |
| `-s, --size` | `0` | Minimum contig size in bp: shorter contigs, and the pairs with a mate on them, are left out. |
| `-f, --filter` | off | Remove the uncuts, loops and impossible pairs. Needs a restriction enzyme. |
| `-D, --duplicates` | off | Remove the PCR duplicates. |
| `-M, --matfmt` | `cool` | Format of the contact map: `cool`, `bg2` or `graal`. |
| `-b, --binning` | `0` | Bin the cool map at this many bp; `0` keeps one bin per fragment. |
| `-z, --zoomify` | `true` | With `--binning`, write a multi-resolution `.mcool`; `-z false` writes one `.cool`. |
| `-B, --balancing-args` | none | More options for `cooler balance` (`-B "--max-iters 500"`), checked before the run starts. Balancing runs on `--threads` processes unless they set `--nproc`. A resolution that fails to balance is logged as a warning. |
| `-E, --exclude` | none | Chromosomes to leave out of the cool map, separated by commas (`chrM,2u`). |
| `-d, --distance-law` | off | Write the table of the distance law. |
| `-c, --centromeres` | none | File of the centromeres' positions, to compute the distance law per arm. |
| `-r, --remove-centromeres` | `0` | bp to leave out on each side of each centromere, for the distance law (`--help` says kb, but the value is read in bp). |
| `-o, --outdir` | the current directory | Output directory. |
| `-P, --prefix` | none | Name of the output files: `sample.cool`, `sample.frags.tsv`, … |
| `-T, --tmpdir` | `<outdir>/tmp` | Directory of the temporary files. |
| `-n, --no-cleanup` | off | Keep the temporary files. |
| `-p, --plot` | off | Save a plot of each step that has one, in `<outdir>/plots/`. |
| `-t, --threads` | `1` | Threads, and processes for balancing. |
| `-F, --force` | off | Replace the outputs and temporary files that exist. |

```bash
# Reads to a fragment-level cool, in the recommended mapping mode
hicstuff pipeline -m 5spm -e DpnII -g genome.fa -o out/ R1.fq.gz R2.fq.gz

# An Arima library, filtered, binned at 1 kb and up, balanced
hicstuff pipeline -m 5spm -e DpnII,HinfI -b 1000 -t 8 -f -D -P sample -g genome.fa \
    -o out/ R1.fq.gz R2.fq.gz

# Reads aligned elsewhere with bwa mem -5SPM, one BAM file per mate
hicstuff pipeline -S bam -m 5spm -e DpnII -g genome.fa -o out/ R1.bam R2.bam

# The lanes of a library, read as one library
hicstuff pipeline -m 5spm -e DpnII -g genome.fa -o out/ lib1_*_R1.fq.gz lib1_*_R2.fq.gz
```

## Processing

### `digest` { #digest }

Cut a genome into restriction fragments, and write `fragments_list.txt` and
`info_contigs.txt`. **API:** `hicstuff.digest.write_frag_info`. **Guide:**
[Digesting a genome](guides/digest.md).

```console
hicstuff digest [OPTIONS] FASTA
```

| Option | Default | Description |
|---|---|---|
| `-e, --enzyme` | *required* | Restriction enzyme, or several separated by commas; or a fragment size in bp, at least 50. |
| `-o, --outdir` | the current directory | Output directory. It must not exist, unless `--force` is given: name one. |
| `-s, --size` | `0` | Minimum contig size in bp: shorter contigs are left out. |
| `-c, --circular` | off | The genome is circular. Changes nothing in the fragments. |
| `-p, --plot` | off | Show a histogram of the fragments' lengths. |
| `-f, --figdir` | none | With `--plot`, save the histogram as `frags_hist.pdf` in this directory. |
| `-F, --force` | off | Write into an output directory that exists. |

```bash
hicstuff digest -e DpnII,HinfI -o digest/ genome.fa
```

### `filter` { #filter }

Remove the uncuts, loops and impossible pairs from a pairs file with
fragments. **API:** `hicstuff.filter.get_thresholds`,
`hicstuff.filter.filter_events`. **Guide:**
[Filtering 3C events](guides/filter.md).

```console
hicstuff filter [OPTIONS] INPUT_PAIRS OUTPUT_PAIRS
```

`INPUT_PAIRS` is an uncompressed pairs file of 9 columns, `OUTPUT_PAIRS` the
pairs kept, in the same format.

| Option | Default | Description |
|---|---|---|
| `-t, --thresholds` | estimated | The thresholds, in restriction sites, as `UNCUT-LOOP` (`4-5`). |
| `-i, --interactive` | off | Show the events by distance, and ask for the thresholds. |
| `-p, --plot` | off | Show the events by distance and by kind, or save them with `--figdir`. |
| `-f, --figdir` | none | With `--plot`, save the plots in this directory: `event_distribution.pdf`, and `event_distance.pdf` when the thresholds are estimated. With `--thresholds`, the directory must exist already. |
| `-P, --prefix` | none | The library's name, in the plots' titles. |

```bash
hicstuff filter --plot --figdir figures/ valid_idx.pairs valid_idx_filtered.pairs
hicstuff filter --thresholds 4-5 valid_idx.pairs valid_idx_filtered.pairs
```

## Alignment

### `iteralign` { #iteralign }

Align the reads of one FASTQ file iteratively: truncated to 20 bases first,
then extended by 20 bases at a time while they do not align. **API:**
`hicstuff.iteralign.iterative_align`. **Guide:**
[Aligning reads](guides/alignment.md#hicstuff-iteralign).

```console
hicstuff iteralign [OPTIONS] READS_FQ
```

| Option | Default | Description |
|---|---|---|
| `-g, --genome` | *required* | The genome's bowtie2 or bwa index, built beforehand, or its FASTA for minimap2. |
| `-o, --out-bam` | *required* | The BAM file to write. |
| `-a, --aligner` | `bowtie2` | `bowtie2`, `bwa` or `minimap2`. |
| `--aligner-backend` | `auto` | How bwa runs, as in [`pipeline`](#pipeline). |
| `-l, --min-len` | `20` | Length the reads are first truncated to. |
| `-R, --read-len` | the first read's | Length of the longest reads. |
| `-t, --threads` | `1` | Threads of the aligner. |
| `-T, --tempdir` | the current directory | Directory of the temporary files. |

```bash
hicstuff iteralign -g genome.fa -o R1.bam -t 8 R1.fq.gz
```

### `cutsite` { #cutsite }

Cut the reads of a pair of FASTQ files at the ligation sites of the enzymes,
and pair the pieces again. **API:** `hicstuff.cutsite.cut_ligation_sites`.
**Guide:** [Aligning reads](guides/alignment.md#hicstuff-cutsite).

```console
hicstuff cutsite [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `-F, --forward` | *required* | The R1 FASTQ file. |
| `-R, --reverse` | *required* | The R2 FASTQ file, its reads in the same order. |
| `-p, --prefix` | *required* | Output prefix: writes `<prefix>_R1.fq.gz` and `<prefix>_R2.fq.gz`. |
| `-e, --enzyme` | *required* | Restriction enzyme, or several separated by commas. |
| `-m, --mode` | `for_vs_rev` | The pieces paired: each of R1 with each of R2 (`for_vs_rev`), every two (`all`), or those next to each other (`pile`). |
| `-s, --seed-size` | `20` | A read is cut at a site only when the piece before it is longer than this; a shorter last piece is dropped. |
| `-t, --threads` | `1` | Threads. |

```bash
hicstuff cutsite -F R1.fq.gz -R R2.fq.gz -p cut/sample -e DpnII -t 8
```

## Matrix operations

### `rebin` { #rebin }

Write a contact map at a coarser resolution, in its format. **API:**
`hicstuff.hicstuff.bin_sparse`, `hicstuff.hicstuff.bin_bp_sparse`. **Guide:**
[Working with matrices](guides/matrices.md#hicstuff-rebin).

```console
hicstuff rebin [OPTIONS] CONTACT_MAP OUT_PREFIX
```

| Option | Default | Description |
|---|---|---|
| `-b, --binning` | `1` | Number of bins to merge (`10`), or a bin size in bp (`10kb`). |
| `-f, --frags` | none | The `fragments_list.txt` of a graal matrix, which needs it. |
| `-c, --chroms` | none | The `info_contigs.txt` of a graal matrix, which needs it. |
| `-F, --force` | off | Replace an output that exists. |

```bash
# Writes out/map_10kb.cool
hicstuff rebin -b 10kb out/abs_fragments_contacts_weighted.cool out/map_10kb
```

### `subsample` { #subsample }

Draw contacts at random from a contact map, every contact with the same
chance. **API:** `hicstuff.hicstuff.subsample_contacts`. **Guide:**
[Working with matrices](guides/matrices.md#hicstuff-subsample).

```console
hicstuff subsample [OPTIONS] CONTACT_MAP SUBSAMPLED_PREFIX
```

| Option | Default | Description |
|---|---|---|
| `-p, --prop` | `0.1` | The proportion of contacts to draw, between 0 and 1, or their number, above 1. |
| `-F, --force` | off | Replace an output that exists. |

```bash
hicstuff subsample -p 20000000 deep.cool deep_20M   # deep_20M.cool
```

### `convert` { #convert }

Write a contact map in another format. **API:**
`hicstuff.io.flexible_hic_loader`, `hicstuff.io.flexible_hic_saver`.
**Guide:** [Working with matrices](guides/matrices.md#hicstuff-convert).

```console
hicstuff convert [OPTIONS] CONTACT_MAP PREFIX
```

| Option | Default | Description |
|---|---|---|
| `-T, --to` | `cool` | The output format: `cool`, `bg2` or `graal`. |
| `-f, --frags` | none | The `fragments_list.txt` of a graal input, which needs it. |
| `-c, --chroms` | none | The `info_contigs.txt` of a graal input, for a graal output. |
| `-g, --genome` | none | The genome's FASTA, to compute the GC content of the bins of a cool or graal output. |
| `-F, --force` | off | Replace an output that exists. |

```bash
hicstuff convert -T graal out/abs_fragments_contacts_weighted.cool graal/sample
hicstuff convert -f fragments_list.txt -c info_contigs.txt \
    abs_fragments_contacts_weighted.txt sample
```

## Analysis and visualisation

### `view` { #view }

Draw a contact map, or the log2 ratio of two, as a heatmap. **API:**
`hicstuff.view.plot_matrix`. **Guide:**
[Working with matrices](guides/matrices.md#hicstuff-view).

```console
hicstuff view [OPTIONS] CONTACT_MAP [CONTACT_MAP2]
```

| Option | Default | Description |
|---|---|---|
| `-o, --output` | a window | The image to write, in the format of its extension. |
| `-b, --binning` | `1` | Number of bins to merge (`10`), or a bin size in bp (`10kb`). |
| `-f, --frags` | none | The `fragments_list.txt` of a graal matrix, for binning in bp, `--lines` and `--region`. |
| `-r, --region` | the whole map | The region to draw (`chr1:1Mb-2Mb`), or two separated by `;` to draw the contacts between them. |
| `-t, --trim` | none | Remove the bins whose number of contacts, in log scale, lies more than this many robust standard deviations (1.4826 times the median absolute deviation) below the median. |
| `-n, --normalize` | off | Balance the map (ICE) before drawing it. |
| `-N, --n-mad` | `3.0` | With `--normalize`, leave out of the balancing the bins more than this many robust standard deviations below the median, as `--trim` measures them. |
| `-d, --despeckle` | off | Replace the outlier pixels by the median of their diagonal. |
| `-T, --transform` | none | Draw `log2`, `log10`, `ln`, `sqrt` or `exp<x>` (`exp0.2`) of the contacts. |
| `-M, --max` | `99%` | Top of the colour scale: a value, or a percentile of the non-zero pixels (`99%`). |
| `-m, --min` | `0` | Bottom of the colour scale, the same way. |
| `-c, --cmap` | `Reds` | Matplotlib colour map. For a ratio, `seismic` unless a diverging one is given. |
| `-l, --lines` | off | Draw lines between chromosomes. Not with `--region`. |
| `-D, --dpi` | `300` | Resolution of the image, in dots per inch. |
| `-C, --circular` | off | The genome is circular. Changes nothing in the drawing. |

```bash
hicstuff view -b 10kb -n -r chr3:1Mb-2Mb -o chr3.png \
    out/sample.mcool::/resolutions/1000
hicstuff view -b 10kb -n -o ratio.png mutant.cool wild_type.cool
```

### `scalogram` { #scalogram }

Draw the contacts of each bin within growing distances of it. **API:**
`hicstuff.hicstuff.scalogram`. **Guide:**
[Working with matrices](guides/matrices.md#hicstuff-scalogram).

```console
hicstuff scalogram [OPTIONS] CONTACT_MAP
```

| Option | Default | Description |
|---|---|---|
| `-o, --output` | a window | The image to write. |
| `-n, --normalize` | off | Balance the map (ICE) first. |
| `-d, --despeckle` | off | Replace the outlier pixels first. |
| `-t, --threads` | `1` | Accepted, but changes nothing: despeckling runs on one thread. |
| `-i, --indices` | the whole map | The bins drawn: bin numbers (`100-200`) or a region (`chr1:1Mb-3Mb`). [Fails without `--range`](guides/matrices.md#hicstuff-scalogram) in 3.3.0. |
| `-r, --range` | every distance | The distances drawn, in bins (`0-50`) or bp (`10kb-500kb`). [Selects the wrong axis](guides/matrices.md#hicstuff-scalogram) in 3.3.0. |
| `-f, --frags` | none | The `fragments_list.txt` of a graal matrix, for regions and distances in bp. |
| `-C, --cmap` | `viridis` | Matplotlib colour map. |

```bash
hicstuff scalogram -n -o scalogram.png out/sample.mcool::/resolutions/5000
```

### `distancelaw` { #distancelaw }

Compute, normalise, average and plot the distance law. **API:**
`hicstuff.distance_law`. **Guide:** [The distance law](guides/distance-law.md).

```console
hicstuff distancelaw [OPTIONS]
```

One of `--pairs` and `--dist-tbl` is required.

| Option | Default | Description |
|---|---|---|
| `-p, --pairs` | none | A pairs file to compute the distance law from. |
| `-d, --dist-tbl` | none | Distance law tables, separated by commas, instead of `--pairs`. |
| `-f, --frags` | none | The `fragments_list.txt` of the pairs: distances are then measured between fragments, and the pairs need the `frag1` and `frag2` columns. |
| `-a, --average` | off | Average the curves of the chromosomes or arms. Needed for several tables. |
| `-b, --big-arm-only` | none | Average only the arms or chromosomes longer than this many bp. |
| `-c, --centromeres` | none | File of the centromeres' positions: one curve per arm. Needs `--pairs` and `--frags`. |
| `-r, --remove-centromeres` | `0` | bp to leave out on each side of each centromere (`--help` says kb, but the value is read in bp). |
| `-C, --circular` | off | Measure distances around circular chromosomes. |
| `-B, --base` | `1.1` | Each distance bin is this many times as wide as the one before it. |
| `-i, --inf` | `3000` | The shortest distance plotted and normalised over, in bp. |
| `-s, --sup` | the longest | The longest distance plotted and normalised over, in bp. |
| `-l, --labels` | the chromosomes, or `Sample <n>` | Names of the curves, separated by commas. |
| `-o, --outputfile-img` | a window | The image to write. |
| `-O, --outputfile-tabl` | none | Write the curves plotted as a table. |

```bash
hicstuff distancelaw -p valid_idx.pairs -f fragments_list.txt -o ps.pdf
hicstuff distancelaw -d wt/distance_law.txt,ko/distance_law.txt -a -l WT,KO -o ps.pdf
```

### `missview` { #missview }

Show which bins of a genome reads of a given length cannot reach. **Guide:**
[Working with matrices](guides/matrices.md#hicstuff-missview).

```console
hicstuff missview [OPTIONS] GENOME OUTPUT
```

`GENOME` is a FASTA file, `OUTPUT` the image to write.

| Option | Default | Description |
|---|---|---|
| `-R, --read-len` | *required* | Length of the simulated reads, in bp. |
| `-b, --binning` | `5000` | Resolution of the map, in bp (`5000` or `5kb`). |
| `-a, --aligner` | `bowtie2` | `bowtie2`, `bwa` (as mapping mode `5spm` aligns) or `minimap2`, for reads of 100 bp or more only. |
| `--aligner-backend` | `auto` | How bwa runs, as in [`pipeline`](#pipeline). |
| `-t, --threads` | `1` | Threads. |
| `-T, --tmpdir` | a new temporary directory | Directory of the temporary files, kept after the run. |
| `-F, --force` | off | Replace the temporary files that exist in `--tmpdir`. |

```bash
hicstuff missview -a bwa -R 75 -b 5kb genome.fa missing.png
```

### `stats` { #stats }

Summarise the log of a pipeline run: reads aligned, contacts kept, pairs
removed by each filter. **API:** `hicstuff.stats.get_pipeline_stats`.

```console
hicstuff stats LOG_FILE
```

```console
$ hicstuff stats out/hicstuff_20261117221914.log
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

A log of a run started from pairs holds no alignment statistics, and `stats`
stops with an error. [Running the pipeline](guides/pipeline.md#reading-the-log)
explains each line.
