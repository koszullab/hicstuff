# Tests for the hicstuff filter module.
# 20190409

import hashlib as hashlib
import os as os
import time
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
import pytest

import hicstuff.distance_law as hcdl

fragments_file = "test_data/fragments_list.txt"
fragments = pd.read_csv(fragments_file, sep="\t", header=0, usecols=[0, 1, 2, 3])
centro_file = "test_data/centromeres.txt"
pairs_reads_file = "test_data/valid_idx_filtered.pairs"
distance_law_file = "test_data/distance_law.txt"
test_xs, test_ps, labels = hcdl.import_distance_law(distance_law_file)


def hash_file(filename):
    """Computes the MD5 hash of a file's content"""
    md5 = hashlib.md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            md5.update(chunk)
    return md5.hexdigest()


def test_distancelaw_cli() -> None:
    """Test the distance law CLI."""
    os.makedirs("test_out", exist_ok=True)
    pairs_reads_file = "test_data/valid_idx_filtered.pairs"
    os.system(
        f"hicstuff distancelaw --pairs {pairs_reads_file} -o test_out/ps.pdf -O test_out/ps.tsv"
    )
    assert os.path.exists("test_out/ps.pdf")
    pairs_reads_file = "test_data/valid_idx_filtered.pairs.gz"
    os.system(
        f"hicstuff distancelaw --average --pairs {pairs_reads_file} -o test_out/ps2.pdf -O test_out/ps2.tsv"
    )
    assert os.path.exists("test_out/ps2.pdf")


def test_export_distance_law():
    """Test exporting the distance law table files"""
    # Test if enable to create an out_dir by default.
    hcdl.export_distance_law(test_xs, test_ps, names=["seq1", "seq2"])
    # assert that file exists
    assert os.path.exists("distance_law.txt")
    # Test error message if not the same numbers of names and chromosomes.
    test = False
    try:
        hcdl.export_distance_law(test_xs, test_ps, names=["seq1"])
    except SystemExit:
        test = True
    assert test


def test_import_distance_law():
    """Test importing distance law table files"""
    xs = hcdl._logbins_xs(fragments, [60000, 20000])
    assert np.all(np.isclose(test_xs[0], xs[0], rtol=0.0001))
    assert np.all(np.isclose(test_xs[1], xs[1], rtol=0.0001))
    assert len(test_ps) == 2 and len(labels) == 2 and len(test_xs) == len(test_ps)
    assert np.isclose(sum(test_ps[0]), 3.0341050807866947e-05, rtol=1e-08)
    assert np.isclose(np.sum(test_ps[1]), 0.00010561980134394403, rtol=1e-08)
    assert list(np.unique(labels[1])) == ["seq2"]


def test_get_chr_segment_bins_index():
    """Test getting the index values of the starting positions of the
    arm/chromosome."""
    # Test with centromeres positions.
    chr_segment_bins = hcdl._get_chr_segment_bins_index(fragments, centro_file)
    assert chr_segment_bins == [0, 129, 129, 409, 409, 474, 474, 564]
    # Test without centromeres positions.
    chr_segment_bins = hcdl._get_chr_segment_bins_index(fragments)
    assert chr_segment_bins == [0, 409, 409, 564]
    # Test with centromeres positions and remove the centromeres.
    chr_segment_bins = hcdl._get_chr_segment_bins_index(fragments, centro_file, 1000)
    assert chr_segment_bins == [0, 121, 134, 409, 409, 463, 480, 564]
    # Test warning message if not the same numbers of chromsome and centromeres.
    hcdl._get_chr_segment_bins_index(fragments.iloc[0:409, :], centro_file, 1000)


def test_get_chr_segment_length():
    """Test getting the length of the arms/chromosomes."""
    chr_length = hcdl._get_chr_segment_length(fragments, [0, 129, 129, 409, 409, 474, 474, 564])
    assert chr_length == [19823, 40177, 9914, 10086]


def test_logbins_xs():
    """Test of the function making the logbins."""
    # Test with default values.
    xs = hcdl._logbins_xs(fragments, [60000, 20000])
    assert len(xs) == 2
    assert np.all(xs[0] == np.unique(np.logspace(0, 115, num=116, base=1.1, dtype=int)))
    # Test changing base.
    xs = hcdl._logbins_xs(fragments, [60000, 20000], base=1.5)
    assert np.all(xs[0] == np.unique(np.logspace(0, 27, num=28, base=1.5, dtype=int)))
    # Test with the circular option.
    xs = hcdl._logbins_xs(fragments, [60000, 20000], circular=True)
    assert np.all(xs[0] == np.unique(np.logspace(0, 108, num=109, base=1.1, dtype=int)))


def test_get_names():
    """Test getting names from a fragment file function."""
    # Test with the centromers option
    names = hcdl._get_names(fragments, [0, 200, 200, 409, 409, 522, 522, 564])
    assert names == ["seq1_left", "seq1_right", "seq2_left", "seq2_right"]
    # Test without the centromers option
    names = hcdl._get_names(fragments, [0, 409, 409, 564])
    assert names == ["seq1", "seq2"]


def test_get_distance_law():
    """Test the general distance_law function."""
    # Create a temporary file.
    distance_law = NamedTemporaryFile("w", delete=False)
    # Test with default parameters.
    hcdl.get_distance_law(pairs_reads_file, fragments_file, out_file=distance_law.name)
    assert os.path.exists(distance_law.name)
    # Test the circular option.
    hcdl.get_distance_law(
        pairs_reads_file, fragments_file, out_file=distance_law.name, circular=True
    )
    assert os.path.exists(distance_law.name)
    # Test the centromere option.
    hcdl.get_distance_law(
        pairs_reads_file,
        fragments_file,
        centro_file=centro_file,
        out_file=distance_law.name,
    )
    assert os.path.exists(distance_law.name)
    # Test error conditions if centromere adn circular.
    test = False
    try:
        hcdl.get_distance_law(
            pairs_reads_file,
            fragments_file,
            centro_file=centro_file,
            out_file=distance_law.name,
            circular=True,
        )
    except SystemExit:
        test = True
    assert test
    # Unlink the temporary file
    os.unlink(distance_law.name)


def test_get_distance_law_matches_reference_table(tmp_path):
    """The fragments path reproduces test_data/distance_law.txt exactly."""
    out = tmp_path / "distance_law.txt"
    hcdl.get_distance_law(pairs_reads_file, fragments_file, out_file=str(out))

    def body(path):
        with open(path) as handle:
            return [line for line in handle if not line.startswith("#")]

    assert body(out) == body(distance_law_file)


def _write_contigs(directory, n_contigs, n_pairs, seed=0):
    """Write a fragments list of ``n_contigs`` 200 bp contigs, of two 100 bp
    fragments each, and an indexed pairs file of ``n_pairs`` pairs joining the
    two fragments of a random contig on one strand, which both paths count.

    Returns the pairs and fragments paths and the contig of each pair.
    """
    directory.mkdir()
    rng = np.random.default_rng(seed)
    names = [f"ctg{i:05d}" for i in range(n_contigs)]
    frags = directory / "fragments_list.txt"
    with open(frags, "w") as handle:
        handle.write("id\tchrom\tstart_pos\tend_pos\tsize\tgc_content\n")
        for name in names:
            handle.write(f"1\t{name}\t0\t100\t100\t0.5\n2\t{name}\t100\t200\t100\t0.5\n")
    contig = rng.integers(0, n_contigs, n_pairs)
    strand = np.where(rng.random(n_pairs) < 0.5, "+", "-")
    pairs = directory / "valid_idx.pairs"
    with open(pairs, "w") as handle:
        handle.write("## pairs format v1.0\n")
        handle.write("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 frag1 frag2\n")
        handle.writelines(f"#chromsize: {name} 200\n" for name in names)
        handle.writelines(
            f"r{i}\t{names[c]}\t10\t{names[c]}\t150\t{s}\t{s}\t{2 * c}\t{2 * c + 1}\n"
            for i, (c, s) in enumerate(zip(contig, strand))
        )
    return str(pairs), str(frags), contig


@pytest.mark.parametrize(
    "with_frags, per_pair_function",
    [(True, "_get_pairs_distance"), (False, "_get_pairs_distance_pos")],
    ids=["fragments", "chromsizes"],
)
def test_distance_law_time_per_pair_does_not_grow_with_contigs(
    tmp_path, monkeypatch, with_frags, per_pair_function
):
    """Each pair used to find its contig with np.searchsorted on a list, which
    copies the list into an array on every call: 1 M pairs over 20,000 contigs
    did not finish in 10 minutes. The time per pair must not depend on the
    number of contigs: it was about 20 times higher over 5,000 contigs than
    over 10.

    A pair's time is the interval between two consecutive calls of the
    per-pair function: one turn of the loop over the pairs, without the
    per-contig setup before it. Their median ignores the turns that another
    process or the garbage collector interrupted.
    """
    original = getattr(hcdl, per_pair_function)
    calls = []

    def stamped(*args, **kwargs):
        calls.append(time.perf_counter())
        return original(*args, **kwargs)

    monkeypatch.setattr(hcdl, per_pair_function, stamped)
    n_pairs = 2_000

    def per_pair(n_contigs):
        del calls[:]
        pairs, frags, contig = _write_contigs(tmp_path / str(n_contigs), n_contigs, n_pairs)
        _, ps, names = hcdl.get_distance_law(pairs, frags if with_frags else None)
        assert len(names) == n_contigs and len(calls) == n_pairs
        return np.median(np.diff(calls)), ps, contig

    few, _, _ = per_pair(10)
    many, ps, contig = per_pair(5_000)
    assert many < 5 * few, (
        f"{many * 1e6:.1f} us per pair over 5,000 contigs, {few * 1e6:.1f} over 10"
    )
    # Every pair counts in its own contig, in the same bin (one distance).
    counted = np.array(ps).sum(axis=1)
    np.testing.assert_allclose(
        counted / counted.sum() * n_pairs, np.bincount(contig, minlength=5_000)
    )


def test_normalize_distance_law():
    """Test function making the average of distance law."""
    # Test normal conditions.
    inf, sup = 1000, 20000
    normed_ps = hcdl.normalize_distance_law(test_xs, test_ps, inf)
    assert len(normed_ps) == 2
    for x, p in zip(test_xs, normed_ps):
        inf_idx = np.searchsorted(x, inf)
        assert np.isclose(sum(p[inf_idx:]), 1.0, rtol=0.05)

    normed_ps = hcdl.normalize_distance_law(test_xs, test_ps, inf, sup)
    for x, p in zip(test_xs, normed_ps):
        inf_idx, sup_idx = np.searchsorted(x, [inf, sup])
        assert np.isclose(sum(p[inf_idx:sup_idx]), 1.0, rtol=0.05)
    # Test sanity check of the size of ps and xs
    test = False
    try:
        hcdl.normalize_distance_law(test_xs, test_ps[0], 1000)
    except SystemExit:
        test = True
    assert test
    # Test warning message if all values in the interval of normalization are zeros.
    hcdl.normalize_distance_law(test_xs, test_ps, 19000)


def test_average_distance_law():
    """Test function making the average of distance law."""
    # Test basic conditions.
    average_xs, average_ps = hcdl.average_distance_law(test_xs, test_ps, 15000, False)
    assert np.all(average_xs == test_xs[0])
    assert np.isclose(sum(average_ps), 6.8119e-05, rtol=10e-7)
    assert np.isclose(np.std(average_ps), 3.5803e-06, rtol=10e-7)
    # Test with big_arm_only.
    average_xs, average_ps = hcdl.average_distance_law(test_xs, test_ps, 25000, True)
    assert np.all(average_xs == test_xs[0])
    assert np.all(average_ps == test_ps[0])
    # Test error if big_arm_only value bigger than the bigger arm.
    test = False
    try:
        hcdl.average_distance_law(test_xs, test_ps, 150000, True)
    except SystemExit:
        test = True
    assert test


def test_slope_distance_law():
    """Test function calculating the slope of the distance law."""
    slope = hcdl.slope_distance_law(test_xs, test_ps)
    assert len(slope) == 2
    assert np.isclose(sum(slope[0]), 18.9329, rtol=0.0001) and np.isclose(
        sum(slope[1]), -2.7459, rtol=0.0001
    )
    assert np.isclose(np.std(slope[0]), 3.9226, rtol=0.0001) and np.isclose(
        np.std(slope[1]), 5.0451, rtol=0.0001
    )
