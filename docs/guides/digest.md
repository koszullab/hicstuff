# Digesting a genome

`hicstuff digest` cuts a genome into the fragments a restriction enzyme makes,
as the pipeline does before assigning the reads to them, and writes the two
tables that describe them: `fragments_list.txt` and `info_contigs.txt`.

```bash
hicstuff digest --enzyme DpnII --outdir digest/ genome.fa
```

```console
INFO :: Genome digested into 564 fragments with a median length of 87.5
```

!!! warning "Give `--outdir`"

    `digest` refuses to write into a directory that exists unless given
    `--force`, and without `--outdir` it writes into the current directory,
    which always exists. Name a new directory, or add `--force`.

## Enzymes

`--enzyme` (`-e`) takes:

- **The name of an enzyme**, as Biopython's
  [Restriction](https://biopython.org/docs/latest/api/Bio.Restriction.html)
  module spells it: `DpnII`, `HindIII`, `MboI`, `HinfI`, … The genome is cut
  where the enzyme cuts, on the forward strand.
- **Several names, separated by commas**, for a library digested with several
  enzymes: `DpnII,HinfI` for the Arima kits. The genome is cut at the sites of
  each of them.
- **A number**, for fragments of that many base pairs, at least 50: each
  contig is cut every N bp from its start, and its last fragment is shorter.

`hicstuff pipeline` takes the same values, and two more: `mnase` and `dnase`,
which cut 100 bp fragments.

## Options

| Option | Effect |
|---|---|
| `--circular` (`-c`) | Changes nothing in the fragments; in `pipeline` and `distancelaw`, it measures distances around circular chromosomes |
| `--size N` (`-s`) | Contigs shorter than N bp are left out of both tables |
| `--plot` (`-p`) | Shows a histogram of fragment lengths, instead of logging their median |
| `--figdir DIR` (`-f`) | With `--plot`, saves the histogram as `DIR/frags_hist.pdf` instead of showing it |
| `--force` (`-F`) | Writes into an output directory that exists |

## Outputs

`fragments_list.txt` has one row per fragment, contig after contig, in the
order of the FASTA file:

```text
id	chrom	start_pos	end_pos	size	gc_content
1	seq1	0	21	21	0.5238095238095238
2	seq1	21	80	59	0.576271186440678
3	seq1	80	328	248	0.5201612903225806
```

`info_contigs.txt` has one row per contig:

```text
contig	length	n_frags	cumul_length
seq1	60000	409	0
seq2	20000	155	409
```

[File formats](../formats.md#fragments_listtxt) describes their columns. The
row of a fragment in `fragments_list.txt`, counted from 0, is its index in
pairs files and graal matrices.

## From Python

```python
import os

import hicstuff.digest as hcd

# Writes fragments_list.txt and info_contigs.txt in digest/, which must exist
os.makedirs("digest", exist_ok=True)
hcd.write_frag_info("genome.fa", "DpnII", output_dir="digest")

# The positions of the sites on one sequence, its start and end included
from Bio.Seq import Seq
hcd.get_restriction_table(Seq("AAGCCGGATCGG"), "HpaII")  # array([ 0,  4, 12])
```
