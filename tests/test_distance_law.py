# Tests for the hicstuff filter module.
# 20190409

import hicstuff.distance_law as hcdl
from tempfile import NamedTemporaryFile
import pandas as pd
import numpy as np
import os as os
import hashlib as hashlib

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


def test_export_distance_law():
    """Test exporting the distance law table files"""
    # Test if enable to create an out_dir by default.
    hcdl.export_distance_law(test_xs, test_ps, names=["seq1", "seq2"])
    assert hash_file("distance_law.txt") == hash_file(distance_law_file)
    os.remove("distance_law.txt")
    # Test error message if not the same numbers of names and chromosomes.
    test = False
    try:
        hcdl.export_distance_law(test_xs, test_ps, names=["seq1"])
    except SystemExit:
        test = True
    assert test


def test_import_distance_law():
    """Test importing distance law table files"""
    xs = hcdl.logbins_xs(fragments, [60000, 20000])
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
    chr_segment_bins = hcdl.get_chr_segment_bins_index(fragments, centro_file)
    assert chr_segment_bins == [0, 129, 129, 409, 409, 474, 474, 564]
    # Test without centromeres positions.
    chr_segment_bins = hcdl.get_chr_segment_bins_index(fragments)
    assert chr_segment_bins == [0, 409, 409, 564]
    # Test with centromeres positions and remove the centromeres.
    chr_segment_bins = hcdl.get_chr_segment_bins_index(fragments, centro_file, 1000)
    assert chr_segment_bins == [0, 121, 134, 409, 409, 463, 480, 564]
    # Test warning message if not the same numbers of chromsome and centromeres.
    hcdl.get_chr_segment_bins_index(fragments.iloc[0:409, :], centro_file, 1000)


def test_get_chr_segment_length():
    """Test getting the length of the arms/chromosomes."""
    chr_length = hcdl.get_chr_segment_length(
        fragments, [0, 129, 129, 409, 409, 474, 474, 564]
    )
    assert chr_length == [19823, 40177, 9914, 10086]


def test_logbins_xs():
    """Test of the function making the logbins."""
    # Test with default values.
    xs = hcdl.logbins_xs(fragments, [60000, 20000])
    assert len(xs) == 2
    assert np.all(
        xs[0] == np.unique(np.logspace(0, 115, num=116, base=1.1, dtype=np.int))
    )
    # Test changing base.
    xs = hcdl.logbins_xs(fragments, [60000, 20000], base=1.5)
    assert np.all(
        xs[0] == np.unique(np.logspace(0, 27, num=28, base=1.5, dtype=np.int))
    )
    # Test with the circular option.
    xs = hcdl.logbins_xs(fragments, [60000, 20000], circular=True)
    assert np.all(
        xs[0] == np.unique(np.logspace(0, 108, num=109, base=1.1, dtype=np.int))
    )


def test_get_names():
    """Test getting names from a fragment file function."""
    # Test with the centromers option
    names = hcdl.get_names(fragments, [0, 200, 200, 409, 409, 522, 522, 564])
    assert names == ["seq1_left", "seq1_rigth", "seq2_left", "seq2_rigth"]
    # Test without the centromers option
    names = hcdl.get_names(fragments, [0, 409, 409, 564])
    assert names == ["seq1", "seq2"]


def test_get_distance_law():
    """Test the general distance_law function."""
    # Create a temporary file.
    distance_law = NamedTemporaryFile("w", delete=False)
    # Test with default parameters.
    hcdl.get_distance_law(pairs_reads_file, fragments_file, out_file=distance_law.name)
    assert hash_file(distance_law.name) == hash_file("test_data/distance_law.txt")
    # Test the circular option.
    hcdl.get_distance_law(
        pairs_reads_file, fragments_file, out_file=distance_law.name, circular=True
    )
    assert hash_file(distance_law.name) == "9925666f7e2013c86ecad9e611cc5382"
    # Test the centromere option.
    hcdl.get_distance_law(
        pairs_reads_file,
        fragments_file,
        centro_file=centro_file,
        out_file=distance_law.name,
    )
    assert hash_file(distance_law.name) == "89ab6fb872601a47c6f2ad5be2e3cdcf"
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
