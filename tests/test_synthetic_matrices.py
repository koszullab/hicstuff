"""Matrix operations (hicstuff.hicstuff, hicstuff.view) on synthetic contact maps."""

import numpy as np
import pytest
from Bio import SeqIO
from matplotlib import pyplot as plt
from scipy.sparse import coo_matrix, triu

import hicstuff.hicstuff as hcs
import hicstuff.io as hio
import hicstuff.view as hcv


@pytest.fixture(scope="module")
def contact_map(matrices):
    """Fragment-level upper-triangle map, fragments and chromosome table."""
    files = matrices["graal"]
    return hio.flexible_hic_loader(files["mat"], files["frags"], files["chroms"])


@pytest.fixture(scope="module")
def chrom_map(contact_map):
    """Upper-triangle map of the longest chromosome alone."""
    mat, frags, _ = contact_map
    idx = np.flatnonzero(frags.chrom == "chr2")
    sub = mat.tocsr()[idx[0] : idx[-1] + 1, idx[0] : idx[-1] + 1]
    return triu(sub, format="coo").astype(float)


def _full(mat):
    dense = mat.toarray()
    return dense + dense.T - np.diag(np.diag(dense))


def test_binning_preserves_contacts(contact_map):
    mat, frags, _ = contact_map
    for factor in (2, 3, 7):
        binned = hcs.bin_sparse(mat, subsampling_factor=factor)
        assert binned.shape == (-(-mat.shape[0] // factor),) * 2
        assert binned.sum() == mat.sum()
        dense = hcs.bin_dense(mat.toarray(), subsampling_factor=factor)
        np.testing.assert_array_equal(dense, binned.toarray())
    positions = frags.start_pos.values
    binned, pos = hcs.bin_bp_sparse(mat, positions, bin_len=2_000)
    assert binned.sum() == mat.sum()
    assert binned.shape[0] == len(pos)
    dense, pos_dense = hcs.bin_bp_dense(mat.toarray(), positions, bin_len=2_000)
    assert dense.sum() == mat.sum()
    np.testing.assert_array_equal(pos_dense, pos)


def test_normalize_sparse(contact_map):
    mat = contact_map[0].astype(float)
    good = hcs.get_good_bins(mat, n_mad=3.0)
    ice = hcs.normalize_sparse(mat, norm="ICE", n_mad=3.0)
    sums = hcs.sum_mat_bins(ice)
    # ICE equalises the coverage of every detectable bin
    np.testing.assert_allclose(sums[good], 1.0, rtol=0.05)
    scn = hcs.normalize_sparse(mat, norm="SCN")
    assert scn.shape == mat.shape and np.isfinite(scn.data).all()
    doubled = hcs.normalize_sparse(mat, norm=lambda m: m * 2)
    assert doubled.sum() == 2 * mat.sum()
    with pytest.raises(Exception, match="Unknown norm"):
        hcs.normalize_sparse(mat, norm="KR")
    assert hcv.normalize(mat, norm="SCN").shape == mat.shape


def test_normalize_dense(chrom_map, capsys):
    dense = _full(chrom_map)
    scn = hcs.normalize_dense(dense, norm="SCN")
    nonzero = dense.sum(axis=1) > 0
    np.testing.assert_allclose(scn.sum(axis=1)[nonzero], 1.0, rtol=0.05)
    np.testing.assert_allclose(scn, scn.T)
    frag = hcs.normalize_dense(dense, norm="frag", order=2, iterations=5)
    assert frag.shape == dense.shape
    np.testing.assert_allclose(hcs.normalize_dense(dense, norm=lambda m: m), dense)
    # mirnylib is not installed: falls back on SCN
    np.testing.assert_allclose(hcs.normalize_dense(dense, norm="mirnylib"), scn)
    assert "mirnylib" in capsys.readouterr().out
    with pytest.raises(Exception, match="Unknown norm"):
        hcs.normalize_dense(dense, norm="KR")


def test_mad_and_trimming(contact_map):
    mat, frags, _ = contact_map
    data = mat.tocoo().data
    assert hcs.mad(mat) == pytest.approx(np.median(np.abs(data - np.median(data))))
    assert hcs.mad(mat.toarray(), axis=-1) == hcs.mad(mat.toarray(), axis=1)
    good = hcs.get_good_bins(mat, n_mad=2.0)
    symmetric = hcs.get_good_bins(mat, n_mad=2.0, symmetric=True)
    assert symmetric.sum() <= good.sum()
    chrom_starts = np.flatnonzero(np.diff(frags.start_pos) < 0) + 1
    trimmed, starts = hcs.trim_sparse(mat, n_mad=2.0, chrom_start=chrom_starts)
    assert trimmed.shape == (good.sum(), good.sum())
    assert len(starts) == 2
    np.testing.assert_array_equal(hcs.trim_dense(_full(mat), n_mad=2.0), _full(mat)[good][:, good])


def test_despeckle(contact_map):
    mat = contact_map[0]
    spiked = mat.tolil()
    spiked[10, 12] = 10_000
    cleaned = hcs.despeckle_simple(spiked.tocoo(), th2=2)
    assert cleaned[10, 12] < 100
    with pytest.raises(AttributeError):
        hcs.despeckle_simple(mat.toarray())


def test_subsample_contacts(contact_map):
    mat = contact_map[0]
    np.random.seed(0)
    half = hcs.subsample_contacts(mat, 0.5)
    assert half.sum() == int(0.5 * mat.sum())
    assert hcs.subsample_contacts(mat, 100).sum() == 100
    with pytest.raises(ValueError):
        hcs.subsample_contacts(mat, -5)


def test_scalogram(chrom_map):
    dense = _full(chrom_map)
    for circ in (False, True):
        scalo = hcs.scalogram(dense, circ=circ, max_range=10)
        assert scalo.shape == dense.shape
        # The first column sums each bin with itself only
        np.testing.assert_array_equal(scalo[:, 0], np.diag(dense))
    assert hcs.scalogram(dense.tolist()).shape == dense.shape
    with pytest.raises(ValueError, match="not square"):
        hcs.scalogram(dense[:, :-1])


def test_compartments(chrom_map):
    dense = _full(chrom_map)
    for normalize in (True, False):
        pc1, pc2 = hcs.compartments(dense, normalize=normalize)
        assert pc1.shape == pc2.shape == (dense.shape[0],)
    with pytest.raises(ValueError, match="not square"):
        hcs.compartments(dense[:, :-1])
    for normalize in (True, False):
        pc1, pc2 = hcs.compartments_sparse(chrom_map, normalize=normalize)
        assert pc1.shape == (dense.shape[0],)
    coeffs = hcs.corrcoef_sparse(coo_matrix(dense).tocsr(), B=coo_matrix(dense).tocsr())
    np.testing.assert_allclose(np.diag(coeffs), 1.0)


def test_distance_law_from_matrix(chrom_map, capsys):
    dense = _full(chrom_map)
    idx, avg = hcs.distance_law_from_mat(dense, log_bins=False)
    assert len(idx) == len(avg) == dense.shape[0]
    assert avg[0] == pytest.approx(np.diag(dense).mean())
    bins, log_avg = hcs.distance_law_from_mat(dense, indices=np.arange(20), base=1.5)
    assert len(bins) == len(log_avg) < dense.shape[0]
    # Too few bins for the log base
    idx, avg = hcs.distance_law_from_mat(np.ones((2, 2)), base=2)
    assert list(idx) == [0, 1]
    assert "Not enough bins" in capsys.readouterr().out


def test_plot_matrix(chrom_map, tmp_path):
    dense = _full(chrom_map)
    out = tmp_path / "map.png"
    hcv.plot_matrix(dense, filename=str(out), title="chr2", dpi=30, chrom_starts=[10])
    assert out.stat().st_size > 0
    fig, ax = plt.subplots()
    hcv.plot_matrix(dense, ax=ax, title="chr2", vmax=5, chrom_starts=[10])
    assert ax.get_title() == "chr2"
    assert len(ax.images) == 1 and len(ax.lines) == 2
    hcv.plot_matrix(dense)  # interactive display is a no-op with Agg
    np.testing.assert_array_equal(hcv.sparse_to_dense(chrom_map, remove_diag=False), dense)


def test_scaffold_tools(synthetic_lib, tmp_path):
    sizes = sorted(synthetic_lib.chrom_sizes.values(), reverse=True)
    assert hcv.scaffold_distribution(str(synthetic_lib.genome), threshold=0) == sizes
    assert hcv.scaffold_distribution(str(synthetic_lib.genome), plot=False, threshold=sizes[1]) == [
        sizes[0]
    ]
    out = tmp_path / "sorted.fa"
    hcv.reorder_fasta(str(synthetic_lib.genome), str(out), threshold=sizes[2])
    records = list(SeqIO.parse(str(out), "fasta"))
    assert [len(r) for r in records] == sizes[:2]
