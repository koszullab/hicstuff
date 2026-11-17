# Matrix operations

`hicstuff.hicstuff` holds the operations on contact maps that the commands
share: normalisation, trimming, despeckling, binning, subsampling, the
distance law of a matrix, compartments and scalograms. Most take a
`scipy.sparse` matrix (the `_sparse` functions) or a `numpy` array (the
`_dense` ones).

These functions are also those of the package: `hicstuff.normalize_sparse` is
`hicstuff.hicstuff.normalize_sparse`. [Using the Python API](../guides/python.md)
shows them at work.

::: hicstuff.hicstuff
