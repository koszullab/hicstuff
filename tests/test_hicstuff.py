#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""hicstuff testing

Basic tests for functions in the hicstuff library.
"""

import random
import numpy as np
import pytest
import hicstuff.hicstuff as hcs
from scipy.sparse import coo_matrix, triu
from inspect import signature, getmembers, isfunction

SIZE_PARAMETERS = ("matrix_size", [5, 10, 20, 50, 100])


def _gen_matrices(size, full_dense=False):
    """
    Make random dense and sparse matrices of given size.
    Parameters
    ----------
    size : int
        The desired number of bins in the matrix.
    full_dense : bool
        If True, the full dense matrix is returned. Otherwise,
        only the upper triangle is returned.
    Returns
    -------
    m_d : numpy.array
        Random dense matrix of size X size.
    m_s : scipy.sparse.coo_matrix
        Random sparse matrix of size X size.
    """
    # pylint: disable=no-member
    m_d = np.array(np.triu(np.random.random((size, size))))
    m_s = coo_matrix(m_d)
    if full_dense:
        m_d += m_d.T - np.diag(np.diag(m_d))
    return m_d, m_s


@pytest.mark.filterwarnings("ignore::PendingDeprecationWarning")
@pytest.mark.parametrize(*SIZE_PARAMETERS)
def test_bin(matrix_size):
    """Test subsample binning on sparse and dense matrices.
    Note: Normal behaviour is to create small bin at the end of chromosome with
    remaining fragments.
    """
    M_d, M_s = _gen_matrices(matrix_size)
    n = M_d.shape[0]
    subsample = max(1, np.random.randint(n // 2))
    remain = 0 if n % subsample == 0 else 1
    exp_n = (n // subsample) + remain
    B_d = hcs.bin_dense(M_d, subsample)
    B_s = hcs.bin_sparse(M_s, subsample)
    # Expected dimensions ?
    assert B_d.shape[0] == exp_n
    assert B_s.shape[0] == exp_n
    # Number of contacts remains the same ?
    assert np.isclose(B_d.sum(), M_d.sum(), rtol=0.0001)
    assert np.isclose(B_s.sum(), M_s.sum(), rtol=0.0001)


@pytest.mark.filterwarnings("ignore::PendingDeprecationWarning")
@pytest.mark.parametrize(*SIZE_PARAMETERS)
def test_bin_bp(matrix_size):
    """
    Test basepair binning on sparse and dense matrices.
    Note: Normal behaviour is to create small bin at the end of chromosome with
    remaining fragments.
    """
    M_d, M_s = _gen_matrices(matrix_size)
    n = M_d.shape[0]
    binsize = 3
    # Split into 2 Chromosomes and simulate genomic positions
    size_ratio = np.random.random_sample()
    size1, size2 = int(n * size_ratio), int(n * (1 - size_ratio))
    if size1 + size2 < n:
        size2 += 1
    pos = np.concatenate([np.array(range(size1)), np.array(range(size2))])
    B_d, _ = hcs.bin_bp_dense(M_d, pos, bin_len=binsize)
    B_s, _ = hcs.bin_bp_sparse(M_s, pos, bin_len=binsize)
    exp_n = len(np.unique(pos[:size1] // binsize)) + len(
        np.unique(pos[size1:] // binsize)
    )
    # Expected dimensions ?
    assert B_d.shape[0] == exp_n
    assert B_s.shape[0] == exp_n
    # Number of contacts remains the same ?
    assert np.isclose(B_d.sum(), M_d.sum(), rtol=0.0001)
    assert np.isclose(B_s.sum(), M_s.sum(), rtol=0.0001)


@pytest.mark.filterwarnings("ignore::PendingDeprecationWarning")
@pytest.mark.parametrize(*SIZE_PARAMETERS)
def test_norm(matrix_size):
    """Test matrix normalization

    Check whether a SCN-normalized matrix has all vectors
    summing to one. Tests both the sparse and dense algorithms.
    """
    M_d, M_s = _gen_matrices(matrix_size, full_dense=True)
    N_d = hcs.normalize_dense(M_d, "SCN", iterations=50)
    N_s = hcs.normalize_sparse(M_s, "ICE", iterations=50, n_mad=1000)
    assert np.isclose(N_d.sum(axis=1), np.ones(matrix_size), rtol=0.0001).all()
    assert np.isclose(hcs.sum_mat_bins(N_s), np.ones(matrix_size), rtol=0.0001).all()
    assert np.isclose(triu(coo_matrix(N_d)).data, N_s.data, rtol=0.000001).all()


@pytest.mark.filterwarnings("ignore::PendingDeprecationWarning")
@pytest.mark.parametrize(*SIZE_PARAMETERS)
def test_trim(matrix_size):
    """
    Generate a random matrix and introduce outlier bins. Check if the correct
    number of bins are trimmed.
    """
    M_d, _ = _gen_matrices(matrix_size)
    
    # Compute thresholds
    sums = M_d.sum(axis=1)
    mad_sums = hcs.mad(M_d, axis=1)
    min_val = np.median(sums) + 3 * mad_sums
    max_val = np.median(sums) + 3 * mad_sums
    # Choose random bins
    trim_bins = np.random.randint(0, M_d.shape[0], 2)
    # Add potential pre-existing outlier bins
    trim_bins = np.append(
        trim_bins, np.where((sums <= min_val) | (sums >= max_val))[0]
    )
    trim_bins = np.unique(trim_bins)
    # Set bins to outlier values
    M_d[:, trim_bins] = M_d[trim_bins, :] = random.choice([min_val, max_val])
    # Compute trimming thresholds again
    sums = M_d.sum(axis=1)
    mad_sums = hcs.mad(M_d, axis=1)
    min_val = np.median(sums) + 3 * mad_sums
    max_val = np.median(sums) + 3 * mad_sums
    # Define bins that need to be trimmed
    trim_bins = np.where((sums <= min_val) | (sums >= max_val))[0]
    trim_shape = M_d.shape[0] - len(trim_bins)
    # Compare expected shape with results
    M_s = coo_matrix(M_d)
    T_d = hcs.trim_dense(M_d, s_min=min_val, s_max=max_val)
    assert T_d.shape[0] == trim_shape
    T_s = hcs.trim_sparse(M_s, s_min=min_val, s_max=max_val)
    assert T_s.shape[0] == trim_shape


@pytest.mark.skip(
    reason="Cannot work unless functions are annotated or a "
    "list of functions taking dense matrices as input is"
    "provided"
)
@pytest.mark.filterwarnings("ignore::PendingDeprecationWarning")
@pytest.mark.parametrize(*SIZE_PARAMETERS)
def test_basic_one_argument_functions(matrix_size):
    """Check all functions

    Generate an NxN matrix at random and feed it to all functions
    with only one argument. This is meant to catch very fundamental
    errors and also facilitate runtime type guessing with MonkeyType.

    Parameters
    ----------
    matrix_size : int
        The size of the random matrix to use for the tests.
    """
    functions_list = getmembers(hcs, isfunction)
    M_d, _ = _gen_matrices(matrix_size)
    for _, func in functions_list:
        params = signature(func).parameters
        if func.__defaults__ is None:
            nb_defaults = 0
        else:
            nb_defaults = len(func.__defaults__)
        annot = params[list(params.keys())[0]].annotation

        if len(params) == 1 or len(params) - nb_defaults == 1:
            try:
                assert func(M_d).any()
            except (ValueError, TypeError, AttributeError):
                pass


@pytest.mark.filterwarnings("ignore::PendingDeprecationWarning")
@pytest.mark.parametrize(*SIZE_PARAMETERS)
def test_corrcoef_sparse(matrix_size):
    """
    Checks if the corrcoeff sparse function yields same results
    as numpy's corrcoeff.
    """
    M_d, M_s = _gen_matrices(matrix_size)
    C_d = np.corrcoef(M_d)
    C_s = hcs.corrcoef_sparse(M_s)
    assert np.isclose(C_s, C_d, rtol=0.0001).all()


@pytest.mark.filterwarnings("ignore::PendingDeprecationWarning")
@pytest.mark.parametrize(*SIZE_PARAMETERS)
def test_compartments_sparse(matrix_size):
    """
    Checks if the eigenvectors obtained by the sparse method match what is
    returned by the dense method.
    """
    # Note: Using full dense matrix for hcs.compartments, but upper triangle
    # spars matrix for hcs.compartments_sparse because the transpose is added in
    # the function.
    M_d, M_s = _gen_matrices(matrix_size, full_dense=True)
    pc1_d, pc2_d = hcs.compartments(M_d, normalize=False)
    pc1_s, pc2_s = hcs.compartments_sparse(M_s, normalize=False)
    assert np.isclose(np.abs(pc1_d), np.abs(pc1_s), rtol=0.01).all()
    assert np.isclose(np.abs(pc2_d), np.abs(pc2_s), rtol=0.01).all()


@pytest.mark.parametrize(*SIZE_PARAMETERS)
def test_distance_law(matrix_size):
    M_d, M_s = _gen_matrices(matrix_size)

    I_d, D_d = hcs.distance_law_from_mat(M_d)
    I_s, D_s = hcs.distance_law_from_mat(M_s)

    assert I_d.shape == I_s.shape
    assert D_d.shape == D_s.shape
