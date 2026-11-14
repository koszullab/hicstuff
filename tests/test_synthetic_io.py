"""hicstuff.io round trips on the synthetic library."""

import bz2
import gzip
import importlib.util
import random
import shutil
import zipfile

import cooler
import numpy as np
import pandas as pd
import pytest

import hicstuff.io as hio
from tests import synthetic as syn

# ---------------------------------------------------------------------------
# Pairs files
# ---------------------------------------------------------------------------

#: How sort_pairs orders each column.
_KEY_TYPES = {"readID": str, "chr1": str, "pos1": int, "chr2": str, "pos2": int}
_KEY_TYPES.update({"frag1": int, "frag2": int})


@pytest.mark.parametrize(
    "keys",
    [
        ["readID"],
        ["chr1", "pos1", "chr2", "pos2"],
        ["chr1", "chr2", "pos1", "pos2"],
        ["frag1", "frag2"],
    ],
    ids="-".join,
)
def test_sort_pairs(fastq_run, tmp_path, keys):
    src = fastq_run().tmp("valid_idx.pairs")
    header = hio.get_pairs_header(str(src))
    with open(src) as handle:
        lines = [line for line in handle if not line.startswith("#")]
    random.Random(0).shuffle(lines)
    unsorted, out = tmp_path / "unsorted.pairs", tmp_path / "sorted.pairs"
    unsorted.write_text("\n".join(header) + "\n" + "".join(lines))

    hio.sort_pairs(str(unsorted), str(out), keys=keys, tmp_dir=str(tmp_path), threads=2)

    new_header = hio.get_pairs_header(str(out))
    assert f"#sorted: {'-'.join(keys)}" in new_header
    assert [h for h in new_header if not h.startswith("#sorted")] == [
        h for h in header if not h.startswith("#sorted")
    ]
    pairs = syn.read_pairs(out)
    assert sorted(pairs.readID) == sorted(line.split("\t")[0] for line in lines)
    order = list(zip(*(pairs[k].astype(_KEY_TYPES[k]) for k in keys)))
    assert order == sorted(order)


def test_sort_pairs_rejects_unknown_column(fastq_run, tmp_path):
    with pytest.raises(KeyError):
        hio.sort_pairs(
            str(fastq_run().tmp("valid_idx.pairs")), str(tmp_path / "x.pairs"), keys=["mapq"]
        )


# ---------------------------------------------------------------------------
# Compressed and sequence files
# ---------------------------------------------------------------------------


@pytest.fixture
def genome_copies(synthetic_lib, tmp_path):
    """The synthetic genome, plain and in each supported compression format."""
    text = synthetic_lib.genome.read_text()
    paths = {"plain": synthetic_lib.genome, "gz": synthetic_lib.genome_gz}
    paths["bz2"] = tmp_path / "genome.fa.bz2"
    with bz2.open(paths["bz2"], "wt") as handle:
        handle.write(text)
    paths["zip"] = tmp_path / "genome.zip"
    with zipfile.ZipFile(paths["zip"], "w") as archive:
        archive.writestr("genome.fa", text)
    return text, paths


@pytest.mark.parametrize("fmt", ["plain", "gz", "bz2", "zip"])
def test_read_compressed(genome_copies, fmt):
    text, paths = genome_copies
    assert hio.is_compressed(str(paths[fmt])) == (fmt != "plain")
    with hio.read_compressed(str(paths[fmt])) as handle:
        assert handle.read() == text
    assert hio.check_is_fasta(str(paths[fmt]))


def test_read_compressed_binary_and_multi_file_zip(genome_copies, tmp_path):
    text, paths = genome_copies
    with hio.read_compressed(str(paths["gz"]), "rb") as handle:
        assert handle.read() == text.encode()
    archive_path = tmp_path / "two.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("a.fa", text)
        archive.writestr("b.fa", text)
    with pytest.raises(OSError, match="single fastq"):
        hio.read_compressed(str(archive_path))


def test_check_is_fasta(tmp_path):
    (tmp_path / "empty.fa").write_text("")
    assert not hio.check_is_fasta(str(tmp_path / "empty.fa"))
    # An aligner index prefix is not a file
    assert not hio.check_is_fasta(str(tmp_path / "missing.fa"))


@pytest.mark.parametrize("threads", [1, 2])
def test_check_fastq_entries(synthetic_lib, tmp_path, threads):
    n = synthetic_lib.params.n_pairs
    assert hio.check_fastq_entries(str(synthetic_lib.reads1), threads) == n
    plain = tmp_path / "reads.fq"
    with gzip.open(synthetic_lib.reads2, "rb") as src, open(plain, "wb") as dst:
        shutil.copyfileobj(src, dst)
    assert hio.check_fastq_entries(str(plain), threads) == n


def test_check_bam_entries(fastq_run, synthetic_lib):
    run = fastq_run()
    assert hio.check_bam_entries(str(run.tmp("for.bam"))) == synthetic_lib.params.n_pairs


def test_check_fasta_index(indexed_genome, synthetic_lib):
    for mode in ("bowtie2", "bwa", "minimap2"):
        assert hio.check_fasta_index(indexed_genome, mode=mode) == indexed_genome
    # No index next to the original genome
    assert hio.check_fasta_index(str(synthetic_lib.genome), mode="bowtie2") is None
    assert hio.check_fasta_index(str(synthetic_lib.genome), mode="bwa") is None


def test_rename_genome(tmp_path, monkeypatch):
    src = tmp_path / "in.fa"
    src.write_text(">my-chr.1|x\nACGTRYKMacgtn\n>other\tdescription\nGGCC\n")
    out = tmp_path / "out.fa"
    hio.rename_genome(str(src), str(out))
    assert out.read_text() == ">my_chr1x\nACGTNNNNacgtn\n>other\nGGCC\n"
    # Default output name
    monkeypatch.chdir(tmp_path)
    hio.rename_genome("in.fa")
    assert (tmp_path / "in_renamed.fa").read_text() == out.read_text()


@pytest.mark.xfail(strict=True, reason="io.reorder_fasta writes sequence lengths, not records")
def test_io_reorder_fasta(synthetic_lib, tmp_path):
    out = tmp_path / "sorted.fa"
    hio.reorder_fasta(str(synthetic_lib.genome), str(out), threshold=15_000)
    assert out.read_text().count(">") == 2


def test_gc_bins(synthetic_lib):
    frags = synthetic_lib.fragments(["DpnII", "HindIII"])
    gc = hio.gc_bins(str(synthetic_lib.genome), frags)
    expected = [
        sum(synthetic_lib.chroms[c][s:e].count(b) for b in "GC") / (e - s)
        for c, s, e in zip(frags.chrom, frags.start_pos, frags.end_pos)
    ]
    np.testing.assert_allclose(gc, expected)
    with pytest.raises(ValueError, match="column starting with end"):
        hio.get_pos_cols(frags.drop(columns="end_pos"))


# ---------------------------------------------------------------------------
# Matrix formats
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def graal(matrices):
    files = matrices["graal"]
    mat, frags, chroms = hio.flexible_hic_loader(files["mat"], files["frags"], files["chroms"])
    return mat, frags, chroms


def test_sparse_matrix_round_trip(graal, tmp_path):
    mat = graal[0]
    path = tmp_path / "mat.tsv"
    hio.save_sparse_matrix(mat, str(path))
    again = hio.load_sparse_matrix(str(path))
    assert again.shape == mat.shape
    np.testing.assert_array_equal(again.toarray(), mat.toarray())
    with pytest.raises(ValueError, match="coo"):
        hio.save_sparse_matrix(mat.tocsr(), str(path))
    binned = hio.load_sparse_matrix(str(path), binning=2)
    assert binned.shape == (-(-mat.shape[0] // 2),) * 2
    assert binned.sum() == mat.sum()


@pytest.mark.xfail(strict=True, reason="binning='auto' divides by zero below 10,000 bins")
def test_sparse_matrix_auto_binning(matrices, synthetic_lib):
    mat = hio.load_sparse_matrix(matrices["graal"]["mat"], binning="auto")
    n = len(synthetic_lib.fragments())
    assert mat.shape == (n, n)


@pytest.mark.xfail(strict=True, reason="a text header makes load_sparse_matrix drop a contact")
def test_sparse_matrix_with_text_header(graal, tmp_path):
    mat = graal[0]
    path = tmp_path / "mat.tsv"
    entries = "".join(f"{r}\t{c}\t{v}\n" for r, c, v in zip(mat.row, mat.col, mat.data.astype(int)))
    path.write_text("id_frag_a\tid_frag_b\tn_contact\n" + entries)
    assert hio.load_sparse_matrix(str(path)).sum() == mat.sum()


def test_cool_round_trip(matrices, tmp_path):
    mat, frags, chroms = hio.load_cool(matrices["cool"]["mat"])
    out = tmp_path / "copy.cool"
    hio.save_cool(str(out), mat, frags.copy(), metadata={"origin": "synthetic"})
    again, frags2, chroms2 = hio.load_cool(str(out))
    np.testing.assert_array_equal(again.toarray(), mat.toarray())
    pd.testing.assert_frame_equal(frags2, frags)
    pd.testing.assert_frame_equal(chroms2, chroms)
    clr = cooler.Cooler(str(out))
    assert clr.info["metadata"] == {"origin": "synthetic"}
    # Columns can be added, and replaced
    hio.add_cool_column(clr, pd.Series(np.arange(len(frags))), "rank", metadata={"unit": "bin"})
    hio.add_cool_column(clr, pd.Series(np.arange(len(frags)) * 2), "rank", dtype=np.int64)
    np.testing.assert_array_equal(clr.bins()[:]["rank"], np.arange(len(frags)) * 2)
    # A table without an id column, a full symmetric matrix
    full = (mat + mat.T).tocoo()
    hio.save_cool(str(tmp_path / "full.cool"), full, frags.drop(columns="id"))
    assert cooler.Cooler(str(tmp_path / "full.cool")).pixels()[:]["count"].sum() == full.sum()


def test_bedgraph2d_round_trip(graal, matrices, tmp_path):
    mat, frags, _ = graal
    path = tmp_path / "mat.bg2"
    hio.save_bedgraph2d(mat, frags, str(path))
    again, frags2, chroms2 = hio.load_bedgraph2d(
        str(path), fragments_file=matrices["graal"]["frags"]
    )
    np.testing.assert_array_equal(again.toarray(), mat.toarray())
    assert list(chroms2.contig) == ["chr1", "chr2", "chr3"]
    # Without fragments, bins are guessed from the matrix: empty ones are lost
    guessed, guessed_frags, _ = hio.load_bedgraph2d(str(path))
    assert guessed.sum() == mat.sum()
    assert len(guessed_frags) <= len(frags)


def _fixed_bins_bg2(pairs, sizes, width, path):
    """Write the pairs binned in ``width`` bp bins as a bg2 file."""

    def bins(chroms, pos):
        start = (pos - 1) // width * width
        end = np.minimum(start + width, chroms.map(sizes))
        return chroms, start, end

    table = pd.concat(
        [
            pd.DataFrame(dict(zip(["c1", "s1", "e1"], bins(pairs.chr1, pairs.pos1)))),
            pd.DataFrame(dict(zip(["c2", "s2", "e2"], bins(pairs.chr2, pairs.pos2)))),
        ],
        axis=1,
    )
    counts = table.value_counts().reset_index().sort_values(["c1", "s1", "c2", "s2"])
    counts.to_csv(path, sep="\t", header=False, index=False)


def test_fixed_bins_bedgraph2d(final_pairs, synthetic_lib, tmp_path):
    """bg2 files with fixed-size bins are reloaded with all their bins."""
    width = 1_000
    pairs = syn.read_pairs(final_pairs)
    sizes = synthetic_lib.chrom_sizes
    # Drop the last, shorter bin of each chromosome (see the xfail below)
    full = {c: size // width * width for c, size in sizes.items()}
    pairs = pairs[(pairs.pos1 <= pairs.chr1.map(full)) & (pairs.pos2 <= pairs.chr2.map(full))]
    path = tmp_path / "fixed.bg2"
    _fixed_bins_bg2(pairs, sizes, width, path)
    mat, frags, chroms = hio.flexible_hic_loader(str(path))
    assert mat.sum() == len(pairs)
    assert set(frags["end_pos"] - frags["start_pos"]) == {width}
    assert list(chroms.contig) == ["chr1", "chr2", "chr3"]


@pytest.mark.xfail(
    strict=True, reason="fixed-bin bg2 fails to load unless chromosome sizes are bin multiples"
)
def test_fixed_bins_bedgraph2d_with_partial_last_bin(final_pairs, synthetic_lib, tmp_path):
    path = tmp_path / "fixed.bg2"
    pairs = syn.read_pairs(final_pairs)
    _fixed_bins_bg2(pairs, synthetic_lib.chrom_sizes, 1_000, path)
    mat, _, _ = hio.flexible_hic_loader(str(path))
    assert mat.sum() == len(pairs)


def test_flexible_hic_saver(graal, tmp_path, caplog):
    mat, frags, chroms = graal
    hio.flexible_hic_saver(mat, str(tmp_path / "bare"), hic_fmt="graal")
    assert "Could not create fragments_list.txt" in caplog.text
    assert "Could not create info_contigs.txt" in caplog.text
    assert hio.load_sparse_matrix(str(tmp_path / "bare.mat.tsv")).sum() == mat.sum()
    with pytest.raises(ValueError, match="Unknown output format"):
        hio.flexible_hic_saver(mat, str(tmp_path / "x"), frags=frags, hic_fmt="hic")
    # The graal loader warns when tables are missing
    loaded, no_frags, no_chroms = hio.flexible_hic_loader(str(tmp_path / "bare.mat.tsv"))
    assert no_frags is None and no_chroms is None
    assert "fragments_file was not provided" in caplog.text


def test_get_hic_format(matrices, tmp_path):
    assert hio.get_hic_format(matrices["graal"]["mat"]) == "graal"
    assert hio.get_hic_format(matrices["bg2"]["mat"]) == "bg2"
    assert hio.get_hic_format(matrices["cool"]["mat"]) == "cool"
    assert hio.get_hic_format("x.mcool::/resolutions/1000") == "cool"
    odd = tmp_path / "odd.tsv"
    odd.write_text("a\tb\tc\td\n")
    with pytest.raises(ValueError, match="Unkown file format"):
        hio.get_hic_format(str(odd))


# ---------------------------------------------------------------------------
# DADE and redis helpers
# ---------------------------------------------------------------------------


@pytest.mark.xfail(strict=True, reason="to_dade_matrix tests the truth value of an array")
def test_to_dade_matrix(tmp_path):
    M = np.array([[5, 8, 3], [8, 3, 5], [3, 5, 5]])
    annotations = np.array([0, 10, 20, 30])
    A = hio.to_dade_matrix(M, annotations=annotations, filename=str(tmp_path / "m.dade"))
    np.testing.assert_array_equal(A[1:, 1:], M)


def test_dade_to_graal(tmp_path):
    header = ['"RST"', "0~chr1~0~1~101", "1~chr1~1~101~251", "2~chr2~0~1~51"]
    rows = [
        ["0~chr1~0~1~101", "5"],
        ["1~chr1~1~101~251", "2", "4"],
        ["2~chr2~0~1~51", "0", "1", "3"],
    ]
    dade = tmp_path / "m.dade"
    dade.write_text("\n".join("\t".join(r) for r in [header] + rows) + "\n")
    outputs = {k: str(tmp_path / f"{k}.txt") for k in ("mat", "contigs", "frags")}
    hio.dade_to_graal(
        str(dade),
        output_matrix=outputs["mat"],
        output_contigs=outputs["contigs"],
        output_frags=outputs["frags"],
    )
    contigs = pd.read_csv(outputs["contigs"], sep="\t")
    assert contigs.to_dict("list") == {
        "contig": ["chr1", "chr2"],
        "length": [250, 50],
        "n_frags": [2, 1],
        "cumul_length": [0, 2],
    }
    frags = pd.read_csv(outputs["frags"], sep="\t")
    assert list(frags.start_pos) == [0, 100, 0]
    assert list(frags["size"]) == [100, 150, 50]


@pytest.mark.skipif(importlib.util.find_spec("redis") is not None, reason="redis is installed")
def test_redis_helpers_require_redis(tmp_path):
    with pytest.raises(SystemExit):
        hio.load_into_redis(str(tmp_path / "m.txt"))
    with pytest.raises(SystemExit):
        hio.load_from_redis("key")
