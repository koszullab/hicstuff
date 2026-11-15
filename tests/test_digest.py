# Tests for the hicstuff digest module
# 20190402
import filecmp
import os
from os.path import join
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
import pytest
from Bio import SeqIO

from hicstuff import digest as hcd
from hicstuff import io as hio


def test_write_frag_info():
    """Test generation of fragments_list.txt and info_contigs.txt"""
    genome = NamedTemporaryFile(delete=False, mode="w")
    seq = "GGAATAGATCAAATGATCCACAGATC"
    genome.write(">seq1\n")
    genome.write(seq)
    genome.close()
    out_dir, tigs, frags = "test_data", "test_tigs", "test_frags"
    hcd.write_frag_info(
        genome.name,
        "DpnII",
        output_contigs=tigs,
        output_frags=frags,
        output_dir=out_dir,
    )
    tigs_df = pd.read_csv(join(out_dir, tigs), delimiter="\t")
    frags_df = pd.read_csv(join(out_dir, frags), delimiter="\t")

    assert tigs_df.length.tolist()[0] == len(seq)
    assert frags_df.start_pos.tolist() == [0, 6, 14, 22]

    os.unlink(genome.name)
    os.remove(join(out_dir, tigs))
    os.remove(join(out_dir, frags))


def test_attribute_fragments():
    """Test the attribution of reads to restriction fragments"""
    idx_pairs = NamedTemporaryFile(delete=False)
    idx_pairs.close()
    restriction_table = {}
    for record in SeqIO.parse("test_data/genome/seq.fa", "fasta"):
        # Get chromosome restriction table
        restriction_table[record.id] = hcd.get_restriction_table(record.seq, "DpnII")
    hcd.attribute_fragments("test_data/valid.pairs", idx_pairs.name, restriction_table)

    assert filecmp.cmp("test_data/valid_idx.pairs", idx_pairs.name)


def test_frag_gc_content_matches_gc_bins(tmp_path):
    """fragments_list.txt holds GC as a 0-1 fraction of the whole fragment,
    ambiguous bases included, as io.gc_bins (used by convert -g) does."""
    genome = tmp_path / "genome.fa"
    genome.write_text(">seq1\nGGAATAGATCAAATGATCCACAGATCNNNNGGCC\n")
    hcd.write_frag_info(str(genome), "DpnII", output_dir=str(tmp_path))
    frags = pd.read_csv(tmp_path / "fragments_list.txt", sep="\t")
    np.testing.assert_allclose(frags.gc_content, [2 / 6, 2 / 8, 4 / 8, 6 / 12])
    np.testing.assert_allclose(frags.gc_content, hio.gc_bins(str(genome), frags))


#: Contig lengths of a small genome cut in 100 bp chunks: ctgA (3 fragments)
#: is shorter than MIN_SIZE, ctgB (4 fragments) and ctgC (3) are not.
CHUNK_GENOME = {"ctgA": 250, "ctgB": 400, "ctgC": 300}
MIN_SIZE = 300
#: readID, chr1, pos1, chr2, pos2 (1-based), and the global 0-based fragment
#: indices expected once ctgA is left out (None: the pair is discarded).
CHUNK_PAIRS = [
    ("r1", "ctgB", 150, "ctgB", 350, 1, 3),
    ("r2", "ctgA", 50, "ctgB", 50, None, None),
    ("r3", "ctgA", 10, "ctgA", 210, None, None),
    ("r4", "ctgB", 400, "ctgC", 101, 3, 5),
    ("r5", "ctgC", 300, "ctgC", 1, 6, 4),
    ("r6", "ctgC", 120, "ctgA", 240, None, None),
]
#: The fragment indices of CHUNK_PAIRS with every contig, ctgA first.
ALL_CONTIGS_FRAGS = {
    "r1": (4, 6),
    "r2": (0, 3),
    "r3": (0, 2),
    "r4": (6, 8),
    "r5": (9, 7),
    "r6": (8, 2),
}
#: #chromsize lines of the pairs header: as hicstuff writes them, without the
#: short contig, in another order than the genome, or none at all.
HEADERS = {
    "header-lists-all": ["ctgA", "ctgB", "ctgC"],
    "header-omits-ctgA": ["ctgB", "ctgC"],
    "header-reordered": ["ctgC", "ctgA", "ctgB"],
    "no-chromsize": [],
}


def _chunk_inputs(tmp_path, listed):
    """FASTA of CHUNK_GENOME and a pairs file of CHUNK_PAIRS whose header
    lists the contigs in ``listed``."""
    rng = np.random.default_rng(0)
    genome = tmp_path / "genome.fa"
    genome.write_text(
        "".join(
            f">{name}\n{''.join(rng.choice(list('ACGT'), size=size))}\n"
            for name, size in CHUNK_GENOME.items()
        )
    )
    pairs = tmp_path / "chunks.pairs"
    header = [
        "## pairs format v1.0",
        "#sorted: readID",
        "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2",
    ] + [f"#chromsize: {name} {CHUNK_GENOME[name]}" for name in listed]
    body = [f"{r}\t{c1}\t{p1}\t{c2}\t{p2}\t+\t-" for r, c1, p1, c2, p2, _, _ in CHUNK_PAIRS]
    pairs.write_text("\n".join(header + body) + "\n")
    table = {
        record.id: hcd.get_restriction_table(record.seq, 100)
        for record in SeqIO.parse(str(genome), "fasta")
    }
    return genome, pairs, table


def _check_chunk_output(idx_pairs, frags_file, header_chroms, expected):
    """The output header lists header_chroms and the output pairs are those
    of ``expected`` (readID: (frag1, frag2)), whose fragment indices point, in
    the fragments list, at the fragment holding each read."""
    expected_header = [
        "## pairs format v1.0",
        "#sorted: readID",
        "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 frag1 frag2",
    ] + [f"#chromsize: {name} {CHUNK_GENOME[name]}" for name in header_chroms]
    assert hio.get_pairs_header(str(idx_pairs)) == expected_header
    out = pd.read_csv(idx_pairs, sep="\t", comment="#", header=None)
    assert out[0].tolist() == list(expected)
    assert list(zip(out[7], out[8])) == list(expected.values())
    frags = pd.read_csv(frags_file, sep="\t")
    for chrom, pos, frag in [(1, 2, 7), (3, 4, 8)]:
        hit = frags.iloc[out[frag]]
        assert (hit.chrom.values == out[chrom].values).all()
        assert (
            (hit.start_pos.values <= out[pos].values - 1)
            & (out[pos].values - 1 < hit.end_pos.values)
        ).all()


@pytest.mark.parametrize("listed", list(HEADERS.values()), ids=list(HEADERS))
def test_attribute_fragments_min_size(tmp_path, caplog, listed):
    """Contigs shorter than min_size are left out of the fragments list by
    write_frag_info: attribute_fragments leaves them out of the header and of
    the fragment indices, and discards the pairs with a mate on one of them,
    whether the pairs header lists them (pipeline started from pairs) or not
    (pairs made from BAM files, whose header has the contigs kept)."""
    genome, pairs, table = _chunk_inputs(tmp_path, listed)
    hcd.write_frag_info(str(genome), 100, min_size=MIN_SIZE, output_dir=str(tmp_path))
    idx_pairs = tmp_path / "idx.pairs"
    n_short = hcd.attribute_fragments(str(pairs), str(idx_pairs), table, min_size=MIN_SIZE)
    assert n_short == 3
    assert "3 pairs with a mate on a contig shorter than 300 bp were discarded" in caplog.text
    expected = {p[0]: (p[5], p[6]) for p in CHUNK_PAIRS if p[5] is not None}
    header_chroms = [name for name in listed if name != "ctgA"]
    _check_chunk_output(idx_pairs, tmp_path / "fragments_list.txt", header_chroms, expected)


@pytest.mark.parametrize("listed", list(HEADERS.values()), ids=list(HEADERS))
def test_attribute_fragments_follows_genome(tmp_path, caplog, listed):
    """Fragments are numbered as in the fragments list, in the order of the
    genome, whatever #chromsize lines the pairs header has: every pair is kept,
    and the header is left as it is."""
    genome, pairs, table = _chunk_inputs(tmp_path, listed)
    hcd.write_frag_info(str(genome), 100, output_dir=str(tmp_path))
    idx_pairs = tmp_path / "idx.pairs"
    assert hcd.attribute_fragments(str(pairs), str(idx_pairs), table) == 0
    assert "discarded" not in caplog.text
    _check_chunk_output(idx_pairs, tmp_path / "fragments_list.txt", listed, ALL_CONTIGS_FRAGS)
