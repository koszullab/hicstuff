# Input and output

`hicstuff.io` reads and writes contact maps in hicstuff's three formats, and
the files around them: pairs, genomes, aligner indexes.

| To | Use |
|---|---|
| Read a map in any format | `flexible_hic_loader` |
| Write a map in any format | `flexible_hic_saver` |
| Tell a map's format | `get_hic_format` |
| Read or write one format | `load_cool`, `save_cool`, `load_sparse_matrix`, `save_sparse_matrix`, `load_bedgraph2d`, `save_bedgraph2d` |
| Sort a pairs file | `sort_pairs` |
| Find an aligner's index | `check_fasta_index` |

[File formats](../formats.md) describes the files.

::: hicstuff.io
