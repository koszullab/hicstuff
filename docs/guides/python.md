# Using the Python API

Every command of hicstuff is a function of its library, and the operations on
matrices that the commands share are functions too. This page shows the
pieces to put together; the [Python API reference](../api/hicstuff.md) lists
them all.

## The modules

| Module | Usual name | What it holds |
|---|---|---|
| `hicstuff.pipeline` | `hpi` | `full_pipeline`, and each of its steps |
| `hicstuff.io` | `hio` | Reading and writing matrices, pairs and genomes |
| `hicstuff.hicstuff` | `hcs` | Operations on matrices: normalising, trimming, binning, … |
| `hicstuff.digest` | `hcd` | Digesting a genome, assigning pairs to fragments |
| `hicstuff.filter` | `hcf` | Filtering 3C events |
| `hicstuff.distance_law` | `hdl` | Computing, averaging and plotting the distance law |
| `hicstuff.view` | `hcv` | Drawing matrices |
| `hicstuff.iteralign`, `hicstuff.cutsite`, `hicstuff.bwamem` | `hci`, `hcc`, `hbm` | Aligning reads |
| `hicstuff.stats` | `hcstats` | Reading the statistics of a run's log |

The functions of `hicstuff.hicstuff` are also those of the package itself:
`hicstuff.normalize_sparse` is `hicstuff.hicstuff.normalize_sparse`.

## Running the pipeline

`full_pipeline` takes the options of `hicstuff pipeline` as arguments:

```python
import hicstuff.pipeline as hpi

hpi.full_pipeline(
    genome="genome.fa",
    input1="R1.fq.gz",
    input2="R2.fq.gz",
    mapping="5spm",
    enzyme="DpnII,HinfI",
    filter_events=True,
    pcr_duplicates=True,
    binning=1000,
    threads=8,
    prefix="sample",
    out_dir="out",
)
```

| Option | Argument |
|---|---|
| `--start-stage` | `start_stage` |
| `--matfmt` | `mat_fmt` |
| `--filter` | `filter_events` |
| `--duplicates` | `pcr_duplicates` |
| `--quality-min` | `min_qual` |
| `--size` | `min_size` |
| `--remove-centromeres` | `remove_centros` |
| `--outdir`, `--tmpdir` | `out_dir`, `tmp_dir` |
| any other `--option-name` | `option_name` |

Unlike the command, `full_pipeline` does not turn `mnase` or `dnase` into a
fragment size: give `enzyme=100`. Several pairs of FASTQ files are lists of
paths, the R1 files in `input1` and the R2 files in `input2`, in the same
order: `full_pipeline` does not split commas or expand patterns, which
`split_inputs` does as the command does.

## Reading a matrix

`flexible_hic_loader` reads any of the three formats, and returns the map as a
sparse matrix of its upper triangle, with the table of its bins and the table
of its chromosomes:

```python
import hicstuff.io as hio

mat, frags, chroms = hio.flexible_hic_loader("out/sample.mcool::/resolutions/5000")

# A graal matrix needs its tables to have them
mat, frags, chroms = hio.flexible_hic_loader(
    "abs_fragments_contacts_weighted.txt",
    fragments_file="fragments_list.txt",
    chroms_file="info_contigs.txt",
)
```

`mat` is a `scipy.sparse.coo_matrix`. `frags` has the columns of
[`fragments_list.txt`](../formats.md#fragments_listtxt) (`id`, `chrom`,
`start_pos`, `end_pos`, and `size` and `gc_content` when the map has them),
and `chroms` the contigs, their lengths, numbers of bins and `cumul_length`:
the index of their first bin, or for a cool file their start in base pairs.

## Processing it

```python
import hicstuff.hicstuff as hcs

# Balance (ICE), leaving out the bins with too few contacts
balanced = hcs.normalize_sparse(mat, norm="ICE", n_mad=3)

# Merge every 5 bins, or bin at 10 kb from the bins' start positions
coarse = hcs.bin_sparse(mat, 5)
mat_10kb, starts = hcs.bin_bp_sparse(mat, frags.start_pos, 10_000)

# Draw a million contacts, or a tenth of them
sub = hcs.subsample_contacts(mat, 1_000_000)
tenth = hcs.subsample_contacts(mat, 0.1)

# Remove the outlier pixels, diagonal by diagonal
clean = hcs.despeckle_simple(balanced)

# Remove the bins with too few contacts: the bins table no longer matches
trimmed, _ = hcs.trim_sparse(mat, n_mad=3)
```

The functions of one chromosome take its block of the map: `compartments_sparse`
(the first two principal components), `distance_law_from_mat`,
`scalogram`.

```python
first = frags.chrom == "chr1"
chr1 = mat.tocsr()[first.to_numpy()][:, first.to_numpy()]
pc1, pc2 = hcs.compartments_sparse(chr1)
```

## Writing it

`flexible_hic_saver` writes a map in any of the three formats, adding the
extension:

```python
hio.flexible_hic_saver(
    balanced, "out/balanced", frags=frags, chroms=chroms, hic_fmt="cool"
)
# out/balanced.cool
```

`save_cool`, `save_sparse_matrix` and `save_bedgraph2d` write one format each.

## Drawing it

```python
import hicstuff.view as hcv

dense = hcv.sparse_to_dense(balanced, remove_diag=False)
hcv.plot_matrix(dense, filename="balanced.png", vmax=0.05, cmap="Reds")
```

Without `filename`, `plot_matrix` shows the map in a window.

## Logging

hicstuff logs to standard error, through the root logger of Python's
`logging` module, at the `INFO` level. To quiet it:

```python
from hicstuff.log import logger

logger.setLevel("WARNING")
```
