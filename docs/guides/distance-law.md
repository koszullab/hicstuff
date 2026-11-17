# The distance law

The distance law, or P(s), is the frequency of contacts between two loci
against the genomic distance s between them. Its decay, and the slope of that
decay in log-log scale, characterise how chromosomes fold: compaction, loop
extrusion, cell-cycle stage. hicstuff computes it from the pairs, for each
chromosome or each chromosome arm, in bins whose width grows exponentially
with the distance.

## In the pipeline

```bash
hicstuff pipeline --mapping 5spm --distance-law --plot --enzyme DpnII \
    --genome genome.fa --outdir out/ R1.fq.gz R2.fq.gz
```

`--distance-law` (`-d`) writes `distance_law.txt`, from the pairs assigned to
fragments, before `--filter` and `--duplicates`. Its values are neither
normalised nor averaged; with `--plot`, `plots/distance_law.pdf` shows them
normalised. The options of the pipeline that change it:

| Option | Effect |
|---|---|
| `--centromeres FILE` (`-c`) | One curve per chromosome arm instead of per chromosome |
| `--remove-centromeres N` (`-r`) | Leaves out N bp on each side of each centromere (`--help` says kb, but the value is read in bp) |
| `--circular` (`-C`) | Measures distances around circular chromosomes: never more than half a chromosome |

The centromere file has one line, the position of each chromosome's
centromere in base pairs, separated by spaces, in the order of the
chromosomes of the genome:

```text
150000 238000 114000
```

## `hicstuff distancelaw`

`distancelaw` computes the distance law from a pairs file, or reads tables
computed before, then normalises, averages and plots them.

```bash
# From the pairs of a run, per chromosome
hicstuff distancelaw --pairs valid_idx.pairs --frags out/fragments_list.txt \
    --outputfile-img ps.pdf --outputfile-tabl ps.txt

# Several samples, each averaged over its chromosomes
hicstuff distancelaw --dist-tbl wt/distance_law.txt,mutant/distance_law.txt \
    --average --labels WT,mutant --outputfile-img ps.pdf
```

Its input is one of:

- **`--pairs FILE`** (`-p`): a pairs file, gzipped or not, such as the
  `valid_idx.pairs.gz` of a run. With `--frags` (`-f`), its
  `fragments_list.txt`, the distances are measured between the fragments of
  the mates, and the file must have the `frag1` and `frag2` columns. Without
  `--frags`, they are measured between the positions of the mates, and any
  pairs file with `#chromsize` lines in its header will do.
- **`--dist-tbl FILE[,FILE…]`** (`-d`): tables written by
  `pipeline --distance-law` or by `distancelaw --outputfile-tabl`. Several
  tables need `--average`.

| Option | Effect |
|---|---|
| `--average` (`-a`) | Averages the curves of a sample's chromosomes or arms into one |
| `--big-arm-only N` (`-b`) | Averages only the arms or chromosomes longer than N bp |
| `--centromeres FILE` (`-c`) | One curve per arm; needs `--pairs` with `--frags` |
| `--remove-centromeres N` (`-r`) | Leaves out N bp on each side of the centromeres |
| `--circular` (`-C`) | Measures distances around circular chromosomes |
| `--base B` (`-B`) | Each bin is B times as wide as the one before it (1.1) |
| `--inf N` (`-i`), `--sup N` (`-s`) | The range of distances plotted and normalised over: 3 kb to the longest distance, by default |
| `--labels A,B` (`-l`) | The names of the curves |
| `--outputfile-img FILE` (`-o`) | Saves the plot, instead of showing it |
| `--outputfile-tabl FILE` (`-O`) | Writes the curves plotted, as a table |

Each curve is normalised so that its values between `--inf` and `--sup` sum
to 1, which makes libraries of different depths comparable. The plot has two
panels: the distance law in log-log scale, and its slope.

## What is counted

Only the pairs whose two mates are on the same strand (`++` or `--`) count,
and only those within one chromosome, or one arm: no uncut or loop gives
them, so they need no filtering. With `--frags`, the distance runs between the
ends of the two fragments that were ligated: their starts for `--` pairs,
their ends for `++` pairs. Without it, it runs between the positions of the
mates. Each distance bin is divided by its area in the contact map, so that
bins of growing width compare.

## The table

```text
## hicstuff v3.3.0
## distance_law
## columns: start_bp	p(s)	chrom
1	1e-09	seq1
2	1e-09	seq1
```

One line per bin: the distance at which the bin starts, in base pairs, the
value of the distance law, and the chromosome, arm or sample. See
[File formats](../formats.md#distance-law).

## From Python

```python
import hicstuff.distance_law as hdl

xs, ps, names = hdl.get_distance_law("valid_idx.pairs", "fragments_list.txt")
ps = hdl.normalize_distance_law(xs, ps)
hdl.plot_ps_slope(xs, ps, labels=names, fig_path="ps.pdf")
```
