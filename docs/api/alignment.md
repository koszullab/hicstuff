# Alignment

Three modules align reads: `hicstuff.iteralign` truncates and extends them
(mapping mode `iterative`), `hicstuff.cutsite` cuts them at their ligation
sites (mapping mode `cutsite`), and `hicstuff.bwamem` runs bwa, as a program
or in-process through pybwa. See [Aligning reads](../guides/alignment.md).

::: hicstuff.iteralign

::: hicstuff.cutsite

::: hicstuff.bwamem
