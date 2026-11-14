"""Self-checks of the synthetic Hi-C generator used by the test-suite."""

import gzip
import re
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from tests import synthetic as syn


def _read_fastq(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    return {lines[i][1:]: lines[i + 1] for i in range(0, len(lines), 4)}


def test_generation_is_deterministic(synthetic_lib, tmp_path):
    syn.generate(tmp_path / "again", synthetic_lib.params)
    for name in ("genome.fa", "genome.fa.gz", "reads_R1.fq.gz", "reads_R2.fq.gz", "truth.tsv"):
        assert (tmp_path / "again" / name).read_bytes() == (synthetic_lib.root / name).read_bytes()
    other = syn.generate(tmp_path / "other", replace(synthetic_lib.params, seed=7, n_pairs=200))
    assert other.chroms != synthetic_lib.chroms
    assert len(other.truth) == 200


def test_genome(synthetic_lib):
    params = synthetic_lib.params
    assert synthetic_lib.chrom_names == ["chr1", "chr2", "chr3"]
    with open(synthetic_lib.genome) as handle:
        assert handle.read().count(">") == 3
    for name, seq in synthetic_lib.chroms.items():
        assert params.chrom_len_range[0] <= len(seq) <= params.chrom_len_range[1]
        assert set(seq) == set("ACGT")
        gc = (seq.count("G") + seq.count("C")) / len(seq)
        assert 0.38 <= gc <= 0.62
        for enzyme, (motif, _) in syn.ENZYMES.items():
            found = [m.start() for m in re.finditer(f"(?={motif})", seq)]
            assert len(found) >= 1
            assert found == list(synthetic_lib.sites[enzyme][name])
    # DpnII fragments are long enough to hold any read entirely.
    sizes = synthetic_lib.fragments("DpnII")["size"]
    assert sizes.min() >= params.dpnii_spacing[0] > params.read_len_range[1]
    both = synthetic_lib.fragments(["DpnII", "HindIII"])
    n_sites = sum(len(s) for e in syn.ENZYMES for s in synthetic_lib.sites[e].values())
    assert len(both) == n_sites + len(synthetic_lib.chroms)


def test_reads_match_genome_at_truth_coordinates(synthetic_lib):
    reads = {"1": _read_fastq(synthetic_lib.reads1), "2": _read_fastq(synthetic_lib.reads2)}
    truth = synthetic_lib.truth
    assert list(reads["1"]) == list(truth.read_id) == list(reads["2"])
    for row in truth.itertuples():
        for mate in "12":
            seq = reads[mate][row.read_id]
            low, high = synthetic_lib.params.read_len_range
            assert low <= len(seq) <= high
            chrom = getattr(row, "chr" + mate)
            if row.kind == "unmapped" and mate == "1":
                assert pd.isna(chrom)
                continue
            assert len(seq) == getattr(row, "len" + mate)
            start, alen = int(getattr(row, "start" + mate)), int(getattr(row, "alen" + mate))
            ref = synthetic_lib.chroms[chrom][start : start + alen]
            if getattr(row, "strand" + mate) == "+":
                assert seq[:alen] == ref
            else:
                assert syn.revcomp(seq) == ref
            bounds = synthetic_lib.boundaries()[chrom]
            local = int(getattr(row, "lfrag" + mate))
            # The aligned part lies inside its restriction fragment.
            assert bounds[local] <= start and start + alen <= bounds[local + 1] + syn.OVERHANG
        if row.chimeric:
            # R1 reads through the GATCGATC junction into the other fragment.
            split = row.chimeric_split
            assert row.alen1 == split < row.len1
            assert reads["1"][row.read_id][split - 4 : split + 4] == syn.JUNCTION
        elif row.kind != "unmapped":
            assert syn.JUNCTION not in reads["1"][row.read_id][: row.len1 - 4]


def test_truth_composition(synthetic_lib):
    params = synthetic_lib.params
    truth = synthetic_lib.truth
    counts = truth.kind.value_counts()
    n = params.n_pairs
    assert counts["uncut"] == round(params.uncut_rate * n)
    assert counts["loop"] == round(params.loop_rate * n)
    assert counts["duplicate"] == round(params.dup_rate * n)
    assert counts["unmapped"] == round(params.unmapped_rate * n)
    assert counts.sum() == n
    contacts = truth[truth.kind == "contact"]
    trans = (contacts.chr1 != contacts.chr2).mean()
    assert trans == pytest.approx(params.trans_fraction, abs=0.02)
    assert contacts.chimeric.mean() == pytest.approx(params.chimeric_fraction, abs=0.02)
    intra = contacts[contacts.chr1 == contacts.chr2]
    assert (intra.frag1 != intra.frag2).all()


def test_duplicates_copy_their_original(synthetic_lib):
    truth = synthetic_lib.truth.set_index("read_id")
    dups = truth[truth.kind == "duplicate"]
    originals = truth.loc[dups.dup_of]
    assert (originals.kind == "contact").all()
    cols = ["chr1", "start1", "strand1", "chr2", "start2", "strand2", "len1", "len2"]
    pd.testing.assert_frame_equal(
        dups[cols].reset_index(drop=True), originals[cols].reset_index(drop=True)
    )
    reads = _read_fastq(synthetic_lib.reads1)
    assert all(reads[d] == reads[o] for d, o in zip(dups.index, dups.dup_of))


def _ordered_type(rows):
    """Strand orientation of intra pairs once mates are ordered by position."""
    first_is_1 = rows.start1 <= rows.start2
    s1 = np.where(first_is_1, rows.strand1, rows.strand2)
    s2 = np.where(first_is_1, rows.strand2, rows.strand1)
    return pd.Series(s1, index=rows.index) + pd.Series(s2, index=rows.index)


def test_injected_events_orientation(synthetic_lib):
    truth = synthetic_lib.truth
    loops = truth[truth.kind == "loop"]
    assert (loops.chr1 == loops.chr2).all() and (loops.frag1 == loops.frag2).all()
    assert (_ordered_type(loops) == "-+").all()
    uncuts = truth[truth.kind == "uncut"]
    assert (uncuts.chr1 == uncuts.chr2).all()
    assert ((uncuts.frag1 - uncuts.frag2).abs() <= 2).all()
    assert (_ordered_type(uncuts) == "+-").all()
    # Genuine contacts come in all four orientations, equally likely.
    contacts = truth[(truth.kind == "contact") & (truth.chr1 == truth.chr2)]
    shares = _ordered_type(contacts).value_counts(normalize=True)
    assert set(shares.index) == {"++", "+-", "-+", "--"}
    assert shares.min() > 0.2


def test_truth_contacts_follow_power_law(synthetic_lib):
    """Contacts per pair of loci decay as s ** -alpha with the distance s
    between ligated fragment ends."""
    contacts = synthetic_lib.truth
    contacts = contacts[(contacts.kind == "contact") & (contacts.chr1 == contacts.chr2)]
    sizes = contacts.chr1.map(synthetic_lib.chrom_sizes)
    # Each contact at distance s stands for (chromosome length - s) pairs of loci.
    edges = np.logspace(np.log10(2_000), np.log10(9_000), 12)
    hist, _ = np.histogram(contacts.lig_dist, bins=edges, weights=1 / (sizes - contacts.lig_dist))
    centers = np.sqrt(edges[1:] * edges[:-1])
    slope = np.polyfit(np.log(centers), np.log(hist / np.diff(edges)), 1)[0]
    assert slope == pytest.approx(-synthetic_lib.params.alpha, abs=0.1)


def test_read_pairs_helper(tmp_path):
    path = tmp_path / "x.pairs"
    path.write_text("## pairs format v1.0\n#columns: readID chr1\nr1\tc1\t1\tc2\t5\t+\t-\n")
    pairs = syn.read_pairs(path)
    assert list(pairs.columns) == ["readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2"]
    assert pairs.loc[0, "pos2"] == 5
