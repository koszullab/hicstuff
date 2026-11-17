# Aligning reads

hicstuff aligns the two mates of each read pair separately, as single-end
reads, then pairs them again by read name. This page covers the aligners, the
four mapping modes of `hicstuff pipeline`, and the two commands that run an
alignment step on its own: `iteralign` and `cutsite`.

## The aligners

`--aligner` (`-a`) chooses one of three:

| Aligner | Command run | Genome |
|---|---|---|
| `bowtie2` (default) | `bowtie2 --very-sensitive-local -p <threads> -x <index> -U <reads>` | A FASTA, or the prefix of a bowtie2 index |
| `bwa` | `bwa mem -t <threads> -v 1 <index> <reads>`, or the same in-process | A FASTA |
| `minimap2` | `minimap2 -2 -t <threads> -ax sr <genome> <reads>` | A FASTA |

From each alignment, the supplementary records are dropped and the reads are
sorted by name; a pair is kept when both mates aligned with a mapping quality
of at least `--quality-min` (30 by default).

bowtie2 and bwa need an index of the genome. The pipeline uses one named after
the FASTA when it finds it, and otherwise builds one in its temporary
directory, at every run. On a large genome, build it once:

```bash
bowtie2-build genome.fa genome.fa   # genome.fa.1.bt2, …
bwa index genome.fa                 # genome.fa.bwt, …
```

minimap2 indexes the genome as it aligns.

## bwa in-process

With `--aligner bwa`, bwa runs either as the `bwa` program or in-process,
through [pybwa](https://github.com/fulcrumgenomics/pybwa), which pip installs
with hicstuff. `--aligner-backend` chooses, on `pipeline`, `iteralign` and
`missview`:

| `--aligner-backend` | bwa runs |
|---|---|
| `auto` (default) | As the `bwa` program if it is on `PATH`, through pybwa otherwise. In mapping mode `5spm`, through pybwa first |
| `binary` | As the `bwa` program |
| `python` | Through pybwa |

Without the option, the `HICSTUFF_ALIGNER_BACKEND` environment variable
chooses, and `auto` applies when it is not set either. The log says which ran:

```console
## aligner: bwa mem in-process, pybwa 2.3.2 (bwa 0.7.18)
```

Both give the same alignments, but for reads with several alignments of equal
score, of which they can make a different one primary: pybwa aligns the reads
of a file in batches, and bwa breaks ties with a hash of each read's position
in its batch. Beyond the first 10,000 reads of a file, the choice can differ
from the bwa program's. In our tests this moved 0.17% of the alignments, none
with a mapping quality of 30 or more. In mapping mode `5spm`, the 5'-most part
of a read is primary, and no such tie arises.

## Mapping modes

`--mapping` (`-m`) chooses how the pipeline aligns the reads:

### `normal`

Each mate is aligned in full. This is the default.

### `iterative`

A mate that spans a ligation junction is made of two loci, and may align to
neither. In iterative mode, each mate is first truncated to its first 20
bases and aligned; those that do not align with a mapping quality of at least
`--quality-min` are extended by 20 bases and aligned again, and so on up to
their full length. The length of the reads is that of the first read, unless
`--read-len` (`-R`) gives it, for files of reads of several lengths.

### `cutsite`

Each read is cut at the ligation sites of the enzyme (`GATCGATC` for DpnII:
the site, filled in and ligated to itself), and the pieces of the two mates are
paired again, each piece of R1 with each piece of R2. A read is only cut at a
site where the piece before it is longer than 20 bases, and a shorter last
piece is dropped. The new pairs are then aligned in `normal` mode.
`cutsite` needs a restriction enzyme: with a fragment size, the pipeline
warns and aligns in `iterative` mode instead.

### `5spm`

The mates are aligned with `bwa mem -5SPM`, as the 4DN, Juicer and Arima
pipelines align them: of a read cut into two parts by a ligation junction,
the 5'-most part is the primary alignment. Each mate is then recorded at its
aligned 5' end, the end of the sequenced molecule, as pairtools records it;
the other modes record the leftmost base of the alignment, which is the 3' end
of a mate on the reverse strand.

`5spm` aligns with bwa only: `--aligner` may be left out, or set to `bwa`. It
runs through pybwa unless `--aligner-backend binary` is given.

```bash
hicstuff pipeline -m 5spm -e DpnII,HinfI -t 8 -o out/ -g genome.fa R1.fq.gz R2.fq.gz
```

## `hicstuff iteralign`

Runs the iterative alignment on one FASTQ file, and writes a BAM file of the
reads it aligned:

```bash
hicstuff iteralign --genome genome.fa --out-bam R1.bam --aligner bwa --threads 8 \
    R1.fq.gz
```

Unlike the pipeline, `iteralign` builds no index: for bowtie2 and bwa,
`--genome` must name an index built beforehand, as above (`genome.fa` once
`bwa index genome.fa` has run). `--min-len` (`-l`) sets the length reads are
first truncated to (20 bases), and `--tempdir` (`-T`) where the truncated
reads are written.

## `hicstuff cutsite`

Cuts the reads of a pair of FASTQ files at the ligation sites of one or more
enzymes, and writes the new pairs to `<prefix>_R1.fq.gz` and
`<prefix>_R2.fq.gz`, for any aligner:

```bash
hicstuff cutsite --forward R1.fq.gz --reverse R2.fq.gz --prefix cut/sample \
    --enzyme DpnII,HinfI --threads 8
```

`--mode` (`-m`) chooses which pieces are paired:

| `--mode` | Pairs made from the pieces of a read pair |
|---|---|
| `for_vs_rev` (default) | Each piece of R1 with each piece of R2: the mode of the pipeline |
| `all` | Every two pieces, whichever mate they come from |
| `pile` | The pieces next to each other in the read pair |

The two files must list the reads in the same order. Each new pair keeps the
name of its read pair, with a number appended to tell them apart (`:00`,
`:01`, …). `--seed-size` (`-s`, 20) is the length a piece must exceed: a read
is only cut where the piece before the site is longer, and a shorter last piece
is dropped.
