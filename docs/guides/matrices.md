# Working with matrices

Once a contact map is written, five commands work on it: `view` draws it,
`rebin` and `subsample` change its resolution and depth, `convert` changes its
format, and `scalogram` draws its contacts by distance. A sixth, `missview`,
looks at a genome before any experiment: which bins its reads could not reach.

All of them read a map in any of hicstuff's three formats: a cool file (or a
resolution of a `.mcool`, as `file.mcool::/resolutions/<N>`), a 2D bedgraph or
a graal matrix. A graal matrix carries no coordinates: give its
`fragments_list.txt` with `--frags` for whatever needs them.

## `hicstuff view`

```bash
hicstuff view --binning 10kb --normalize --output map.png \
    out/abs_fragments_contacts_weighted.cool
```

`view` draws the map as a heatmap, and saves it with `--output` (`-o`), at
`--dpi` dots per inch (300); without `--output`, it opens a window. Before
drawing, it applies, in this order:

| Option | Effect |
|---|---|
| `--binning N` (`-b`) | Merges every N bins into one |
| `--binning 10kb` | Bins the map at 10 kb (also `500bp`, `1Mb`); a graal matrix needs `--frags` |
| `--trim N` (`-t`) | Removes the bins whose number of contacts, in log scale, lies more than N robust standard deviations (1.4826 times the median absolute deviation) below the median |
| `--normalize` (`-n`) | Balances the map (ICE). Bins more than `--n-mad` (3) robust standard deviations below the median are left out of the balancing, and set to zero |
| `--region chr1:1Mb-2Mb` (`-r`) | Draws this region only; `chr1:1Mb-2Mb;chr2:3Mb-4Mb` draws the contacts between two regions |
| `--despeckle` (`-d`) | Replaces the outlier pixels, more than two standard deviations above the median of their diagonal, by that median |
| `--transform log10` (`-T`) | Draws a function of the contacts: `log2`, `log10`, `ln`, `sqrt` or `exp<x>` (as `exp0.2`) |

The colour scale runs from `--min` (`-m`, 0) to `--max` (`-M`, `99%`). Either
is a value, or a percentile of the non-zero pixels with `%`. `--cmap` (`-c`) is
any matplotlib colour map (`Reds`). `--lines` (`-l`) draws a line between
chromosomes; a graal matrix needs `--frags` for it, and it cannot be combined
with `--region`.

### Comparing two maps

```bash
hicstuff view --binning 10kb --normalize --output ratio.png mutant.cool wild_type.cool
```

Given a second map, `view` draws the log2 ratio of the first to the second,
after processing both the same way. They must have the same bins. The colour
scale runs from -2 to 2, and the colour map becomes `seismic` unless a
diverging one is given (`RdBu`, `coolwarm`, `PiYG`, …).

## `hicstuff rebin`

```bash
hicstuff rebin --binning 10kb out/abs_fragments_contacts_weighted.cool out/map_10kb
```

`rebin` writes the map at a coarser resolution, in the format of the input:
`out/map_10kb.cool` here. `--binning` takes a number of bins to merge (`10`) or
a size in base pairs (`10kb`). A graal matrix needs `--frags` and `--chroms`
(its `info_contigs.txt`), and gives three files: `<prefix>.mat.tsv`,
`<prefix>.frags.tsv` and `<prefix>.chr.tsv`. `--force` (`-F`) replaces an
output that exists.

## `hicstuff subsample`

```bash
hicstuff subsample --prop 0.5 out/abs_fragments_contacts_weighted.cool out/half
```

`subsample` draws contacts at random from the map: a proportion of them with
a value between 0 and 1 (`--prop 0.5`, 0.1 by default), a number of them above
1 (`--prop 1000000`). Every contact has the same chance to be drawn, so each
pixel keeps its share of the map. The output, `out/half.cool`, has the format
of the input. To compare two libraries at the same depth, subsample the
deeper one to the number of contacts of the other.

## `hicstuff convert`

```bash
hicstuff convert --to graal out/abs_fragments_contacts_weighted.cool out/graal/sample
```

`convert` writes the map in another format, `--to` (`-T`): `cool` (the
default), `bg2` or `graal`. The input format is recognised from the file. A
graal input needs its `--frags`, and its `--chroms` to write a graal output.
`--genome` (`-g`) computes the GC content of every bin from the genome's
FASTA, for the bins table of a cool output, or the fragments table of a graal
output. The contacts are truncated to whole numbers.

!!! warning "A graal output converted from a cool file"

    Its `<prefix>.chr.tsv` gives `cumul_length` in base pairs, where an
    `info_contigs.txt` counts fragments. Use the `info_contigs.txt` of the run
    with the graal matrix.

| `--to` | Writes |
|---|---|
| `cool` | `<prefix>.cool` |
| `bg2` | `<prefix>.bg2` |
| `graal` | `<prefix>.mat.tsv`, `<prefix>.frags.tsv`, `<prefix>.chr.tsv` |

## `hicstuff scalogram`

```bash
hicstuff scalogram --normalize --output scalogram.png \
    out/sample.mcool::/resolutions/5000
```

A scalogram shows, for every bin along the genome (horizontally), the contacts
it makes within a distance that grows upwards: at each height, the sum of its
contacts with the bins at most that far from it. Domain boundaries and
compaction show at the scales they act on.

!!! warning "`--indices` and `--range` select the wrong axis"

    In hicstuff 3.3.0, `--range` crops the bins drawn rather than the
    distances, and `--indices` fails unless `--range` is given too. Draw the
    whole map, or a map cut beforehand to the region of interest.

| Option | Effect |
|---|---|
| `--indices` (`-i`) | The bins drawn: a range of bin numbers (`100-200`) or a region (`chr1:1Mb-3Mb`) |
| `--range` (`-r`) | The distances drawn: in bins (`0-50`), or in base pairs (`10kb-500kb`), converted with the mean bin size |
| `--normalize` (`-n`) | Balances the map first (ICE) |
| `--despeckle` (`-d`) | Replaces outlier pixels first (`--threads` is accepted, but changes nothing) |
| `--cmap` (`-C`) | The matplotlib colour map (`viridis`) |
| `--frags` (`-f`) | The `fragments_list.txt` of a graal matrix, for regions and distances in base pairs |
| `--output` (`-o`) | The image to write; without it, a window opens |

## `hicstuff missview`

```bash
hicstuff missview --aligner bwa --read-len 75 --binning 5000 genome.fa missing.png
```

Before sequencing, `missview` shows which bins of a genome the reads would not
reach, at a read length: those whose sequence is repeated elsewhere, or
otherwise unmappable. It cuts the genome into reads of `--read-len` (`-R`)
bases starting at every position, aligns them, and builds their contact map at
`--binning` (`-b`, 5000 bp). A bin is missing when its reads aligned less often
than the reads of the whole genome did. The image draws the missing bins'
rows and columns white, and its title gives their proportion.

`--aligner` (`-a`) is `bowtie2` (default), `bwa` or `minimap2`; use minimap2
only for reads of 100 bp or more. `--threads` (`-t`), `--tmpdir` (`-T`) and
`--aligner-backend` work as in the pipeline.
