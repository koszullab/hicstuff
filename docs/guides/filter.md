# Filtering 3C events

Among the pairs of a Hi-C library, some come from fragments that were never
cut (**uncuts**) or that ligated to themselves (**loops**). Their mates map a
few restriction fragments apart, in a known orientation, and the contacts they
add between neighbouring fragments reflect the protocol, not the chromosome's
folding. [Concepts](../concepts.md#3c-events) explains the events; this page
explains how hicstuff finds and removes them.

## In the pipeline

```bash
hicstuff pipeline --mapping 5spm --filter --plot --enzyme DpnII --genome genome.fa \
    --outdir out/ R1.fq.gz R2.fq.gz
```

`--filter` (`-f`) estimates two thresholds from the library, removes the
uncuts, loops and impossible pairs, and builds the map from the pairs left,
which it writes as `valid_idx_filtered.pairs.gz`. It needs a restriction
enzyme: with a fragment size, the run refuses it. The log records the
thresholds and what was removed; on the test library of the
[Quickstart](../quickstart.md):

```console
INFO :: Filtering with thresholds: uncuts=7 loops=5
INFO :: Proportion of inter contacts: 4.83% (intra: 276, inter: 14)
INFO :: 7601 pairs discarded: Loops: 23, Uncuts: 7578, Weirds: 0
INFO :: 290 pairs kept (3.68%)
```

With `--plot`, `plots/event_distance.pdf` shows the events by distance with
the thresholds, and `plots/event_distribution.pdf` the proportion of each.

## What is removed

For every pair whose two mates are on the same chromosome, hicstuff orders the
mates by position and reads two things: the **event**, from their strands, and
the **distance** between them, as the number of restriction sites that
separate them (0 when both mates are on the same fragment).

| Event | Strands | Removed when |
|---|---|---|
| Uncut | `+-` | The distance is at most the uncut threshold |
| Loop | `-+` | The distance is at most the loop threshold |
| Weird | `++` or `--` | Both mates are on the same fragment |

All other pairs are kept, the pairs between chromosomes among them.

## How the thresholds are estimated

Far from each other, the four orientations of a pair are equally frequent:
which strand each mate is read from is chance. Near each other, uncuts add
`+-` pairs and loops `-+` pairs. The thresholds are the distances at which
`+-` and `-+` pairs stop being over-represented:

1. About a million pairs are drawn from the whole library (all of them if it
   has fewer). The draw hashes each line, so it is the same whatever the order
   of the file.
2. The pairs on one chromosome are counted by event and by distance, from 0 to
   49 restriction sites.
3. At each distance, hicstuff takes the median of the four (log) counts, and
   the median absolute deviation of all of them, as the spread to expect.
4. The uncut threshold is the shortest distance, of 2 sites or more, at which
   the `+-` count lies within one expected standard deviation of the median;
   the loop threshold is the same for `-+`.

When no such distance exists for one of them, the run stops with an error,
and `hicstuff filter --interactive` shows the counts to choose from.

## `hicstuff filter`

The same filter, on a pairs file with fragments (9 columns, uncompressed):

```bash
zcat out/valid_idx.pairs.gz > valid_idx.pairs
hicstuff filter --plot --figdir figures/ valid_idx.pairs valid_idx_filtered.pairs
```

| Option | Effect |
|---|---|
| `--thresholds 4-5` (`-t`) | Uses these thresholds, uncuts then loops, instead of estimating them |
| `--interactive` (`-i`) | Shows the counts by distance, and asks for the two thresholds |
| `--plot` (`-p`) | Shows the plots, or saves them in `--figdir` |
| `--figdir DIR` (`-f`) | With `--plot`, saves the plots in DIR: `event_distribution.pdf`, and `event_distance.pdf` when the thresholds are estimated. With `--thresholds`, DIR must exist already |
| `--prefix NAME` (`-P`) | The library's name, in the plots' titles |

### Checking the thresholds

On `event_distance.pdf`, each event is a curve of counts (log scale) against
the distance in restriction sites. The pairs removed are drawn in colour, up to
each threshold, and the pairs kept in green; the two thresholds are vertical
lines. Past them, the four curves should run together.

When they do not, choose the thresholds by eye, filter with
`hicstuff filter --thresholds`, and build the map from the filtered pairs:

```bash
hicstuff filter --thresholds 4-5 valid_idx.pairs valid_idx_filtered.pairs
hicstuff pipeline -S pairs_idx -e DpnII -g genome.fa -o out_filtered/ \
    valid_idx_filtered.pairs
```

## PCR duplicates

`hicstuff pipeline --duplicates` (`-D`) removes the PCR duplicates: among the
pairs with the same chromosomes, positions, strands and fragments for both
mates, it keeps one. It runs after `--filter`, writes
`valid_idx_pcrfree.pairs.gz`, and logs the proportion removed:

```console
INFO :: 0.0% PCR duplicates have been filtered out (0/290 pairs)
INFO :: 290 pairs remaining after removing PCR duplicates
```

## From Python

```python
import hicstuff.filter as hcf

uncut, loop = hcf.get_thresholds("valid_idx.pairs")
hcf.filter_events("valid_idx.pairs", "valid_idx_filtered.pairs", uncut, loop)
```
