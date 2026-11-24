# Pipeline

`hicstuff.pipeline.full_pipeline` is `hicstuff pipeline`: its arguments are
the command's options (see [Using the Python API](../guides/python.md#running-the-pipeline)).
Each step is a function of its own: `align_reads`, `bam2pairs`,
`filter_pcr_dup`, `pairs2cool`, `pairs2binnedcool`, `pairs2metacool`,
`cool2mcool`, `balance`, `pairs2matrix`. `parse_binning`, `read_bins` and
`matrix_bins` read the bins of `--binning`: a bin size, or a BED file of bins.

`hicstuff.stats` reads the statistics of a run from its log, as
`hicstuff stats` does, and `hicstuff.log` holds the logger every module writes
to.

::: hicstuff.pipeline

::: hicstuff.stats

::: hicstuff.log
