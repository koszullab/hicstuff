#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Common Hi-C functions

A bunch of handy functions for processing Hi-C data
(mainly in the form of matrices):

    -Normalizations
    -Interpolations
    -Filters
    -Removing artifacts
    -Quick sum-pooling (aka 'binning') in sparse and dense form
    -Simple models with parameter estimation
    -Computing best-matching 3D structures
    -Various metrics in use among Hi-C people for eyecandy purposes
     (directional index, domainograms, etc.)

These functions are meant to be simple and relatively quick
as-is implementations of procedures described in Hi-C papers.
"""

import numpy as np
import string
import collections
import itertools
import warnings
import functools
from scipy.sparse import coo_matrix, csr_matrix, lil_matrix
from scipy.sparse.linalg import eigsh
from scipy.linalg import eig, toeplitz
import scipy.sparse as sparse
from scipy.sparse import issparse
import copy
import random
import pandas as pd
import sys
from hicstuff.log import logger


def distance_law_from_mat(matrix, indices=None, log_bins=True, base=1.1):
    """Compute distance law as a function of the genomic coordinate aka P(s).
    Bin length increases exponentially with distance if log_bins is True. Works
    on dense and sparse matrices. Less precise than the one from the pairs.

    Parameters
    ----------
    matrix : numpy.array or scipy.sparse.coo_matrix
        Hi-C contact map of the chromosome on which the distance law is
        calculated.
    indices : None or numpy array
        List of indices on which to compute the distance law. For example
        compartments or expressed genes.
    log_bins : bool
        Whether the distance law should be computed on exponentially larger
        bins.

    Returns
    -------
    numpy array of floats :
        The start index of each bin.
    numpy array of floats :
        The distance law computed per bin on the diagonal

    Examples
    --------
    >>> import numpy as np
    >>> mat = np.array([[3, 2, 1], [2, 3, 2], [1, 2, 3]])
    >>> idx, avg = distance_law_from_mat(mat, log_bins=False)
    >>> idx
    array([0, 1, 2])
    >>> avg
    array([3., 2., 1.])
    """

    n = min(matrix.shape)
    included_bins = np.zeros(n, dtype=bool)
    if indices is None:
        included_bins[:] = True
    else:
        included_bins[indices] = True
    D = np.array(
        [
            np.average(matrix.diagonal(j)[included_bins[: n - j]])
            for j in range(n)
        ]
    )
    if not log_bins:
        return np.array(range(len(D))), D
    else:
        n_bins = int(np.log(n) / np.log(base) + 1)
        logbin = np.unique(
            np.logspace(0, n_bins - 1, num=n_bins, base=base, dtype=np.int)
        )
        logbin = np.insert(logbin, 0, 0)
        logbin[-1] = min(n, logbin[-1])
        if n < logbin.shape[0]:
            print("Not enough bins. Increase logarithm base.")
            return np.array(range(len(D))), D
        logD = np.array(
            [
                np.average(D[logbin[i - 1] : logbin[i]])
                for i in range(1, len(logbin))
            ]
        )
        return logbin[:-1], logD


def despeckle_simple(B, th2=2, threads=1):
    """Single-chromosome despeckling

    Simple speckle removing function on a single chromomsome. It also works
    for multiple chromosomes but trends may be disrupted.

    Parameters
    ----------
    B : scipy.sparse.csr
        The input matrix to despeckle, in sparse (csr) format.
    th2 : float
        The number of standard deviations above the mean beyond which
        despeckling should be performed
    threads : int
        The number of CPU processes on which the function can run in parallel.

    Returns
    -------
    array_like
        The despeckled matrix, in the same format it was given.
    """
    try:
        if B.getformat() != "csr":
            B = csr_matrix(B)
    except AttributeError:
        print("Error: You must provide a sparse matrix in csr format.")
        raise

    A = copy.copy(B)
    n1 = A.shape[0]
    medians = np.zeros(n1)
    stds = np.zeros(n1)
    # Faster structure for editing values
    A = lil_matrix(A)
    for u in range(n1):
        diag = B.diagonal(u)
        medians[u] = np.median(diag)
        stds[u] = np.std(diag)
    for nw in range(n1):
        diag = A.diagonal(nw)
        diag[diag > medians[nw] + th2 * stds[nw]] = medians[nw]
        A.setdiag(diag, nw)

    return csr_matrix(A)


def despeckle_local(M, stds=2, width=2):
    """Replace outstanding values (above stds standard deviations)
    in a matrix by the average of a surrounding window of desired width.
    """

    N = np.array(M, dtype=np.float64)
    n, m = M.shape
    for i, j in itertools.product(
        range(width, n - width), range(width, m - width)
    ):
        square = M[i - width : i + width, j - width : j + width]
        avg = np.average(square)
        std = np.std(square)
        if M[i, j] >= avg + stds * std:
            N[i, j] = avg
    return (N + N.T) / 2


def bin_dense(M, subsampling_factor=3):
    """
    Wraps the bin_sparse function to apply it on dense matrices. Bins are merged
    by groups of N to produce a lower resolution matrix.

    Parameters
    ----------
    M : numpy.array
        2D array containing the Hi-C contact map
    subsampling_factor : int
        The number of bins to include in each group (subsample).

    Returns
    -------
    out_M : numpy.array
        The subsamples matrix, with a resolution lower than the input by a defined factor.
    """
    S = coo_matrix(M)
    out_S = bin_sparse(S, subsampling_factor=subsampling_factor)
    out_M = out_S.todense()

    return out_M


def bin_sparse(M, subsampling_factor=3):
    """
    Bins a sparse matrix by combining bins into groups of user defined size. Binsize
    is independent of genomic coordinates. Remaining rows and cols are put into a
    smaller bin at the end.

    Parameters
    ----------
    M : scipy.sparse.coo_matrix
        The input Hi-C matrix in a sparse format.
    subsampling_factor : int
        The number of bins to include in each group (subsample).

    Returns
    -------
    scipy.sparse.coo_matrix
        The subsamples matrix, with a resolution lower than the input by a defined factor.
        
    """

    N = M.tocoo()
    n, m = N.shape
    row, col, data = N.row, N.col, N.data

    # Divide row and column indices - duplicate coordinates are added in
    # sparse matrix construction
    remain_m = 0 if m % subsampling_factor == 0 else 1
    remain_n = 0 if n % subsampling_factor == 0 else 1
    binned_row = row // subsampling_factor
    binned_col = col // subsampling_factor
    binned_n = (n // subsampling_factor) + remain_n
    binned_m = (m // subsampling_factor) + remain_m

    # Sum data over duplicate entries
    binned = pd.DataFrame({"row": binned_row, "col": binned_col, "dat": data})
    binned = binned.groupby(["row", "col"], sort=False).sum().reset_index()
    return coo_matrix(
        (binned.dat, (binned.row, binned.col)), shape=(binned_n, binned_m)
    )


def bin_matrix(M, subsampling_factor=3):
    """Bin either sparse or dense matrices.
    """

    try:
        from scipy.sparse import issparse

        if issparse(M):
            return bin_sparse(M, subsampling_factor=subsampling_factor)
        else:
            raise ImportError
    except ImportError:
        return bin_dense(M, subsampling_factor=subsampling_factor)


def bin_annotation(annotation=None, subsampling_factor=3):
    """Perform binning on genome annotations such as contig information or bin
    positions.
    """

    if annotation is None:
        annotation = np.array([])
    n = len(annotation)
    binned_positions = [
        annotation[i] for i in range(n) if i % subsampling_factor == 0
    ]
    if len(binned_positions) == 0:
        binned_positions.append(0)
    return np.array(binned_positions)


def bin_measurement(measurement=None, subsampling_factor=3):
    """Perform binning on genome-wide measurements by summing each component
    in a window of variable size (subsampling_factor).
    """

    subs = int(subsampling_factor)
    if measurement is None:
        measurement = np.array([])
    n = len(measurement)
    binned_measurement = [
        measurement[i - subs + 1 : i].sum()
        for i in range(n)
        if i % subs == 0 and i > 0
    ]
    return np.array(binned_measurement)


def build_pyramid(M, subsampling_factor=3):
    """Iterate over a given number of times on matrix M
    so as to compute smaller and smaller matrices with bin_dense.
    """

    subs = int(subsampling_factor)
    if subs < 1:
        raise ValueError(
            "Subsampling factor needs to be an integer greater than 1."
        )
    N = [M]
    while min(N[-1].shape) > 1:
        N.append(bin_matrix(N[-1], subsampling_factor=subs))
    return N


def bin_bp_dense(M, positions, bin_len=10000):
    """Perform binning with a fixed genomic length in
    base pairs on a dense matrix. Fragments will be binned such
    that their total length is closest to the specified input.
    If a contig list is specified, binning will be performed
    such that fragments never overlap two contigs. Fragments longer
    than bin size will not be split, which can result in larger bins.
    The last smaller bin of the chromosome will be merged with the
    previous one.

    Parameters
    ----------
    M : 2D numpy array of ints or floats
        The Hi-C matrix to bin in dense format
    positions : numpy array of int
        List of 0-based basepair start positions of fragments bins
    bin_len : int
        Bin length in basepairs

    Returns
    -------
    2D numpy array of ints of floats :
        Binned matrix
    numpy array of ints :
        List of binned fragments positions in basepair
    """
    # Just converting to sparse and passing to sparse function
    S = coo_matrix(M)
    out_S, out_pos = bin_bp_sparse(S, positions, bin_len=bin_len)
    out_M = out_S.todense()

    return out_M, out_pos


def bin_bp_sparse(M, positions, bin_len=10000):
    """
    Perform binning with a fixed genomic length in
    base pairs on a sparse matrix. Fragments will be binned such
    that their total length is closest to the specified input.
    Binning will be performed such that fragments never overlap two
    contigs. Fragments longer than bin size will not be split, which
    can result in larger bins. The last smaller bin of the chromosome
    will be merged with the previous one.

    Parameters
    ----------
    M : sparse numpy matrix
        Hi-C contact matrix in sparse format.
    positions : numpy array of ints
        Start positions of fragments in the matrix, in base pairs.
    bin_len : int
        Desired length of bins, in base pairs

    Returns
    -------
    sparse scipy sparse coo_matrix:
        The binned sparse matrix in COO format.
    list of ints:
        The new bin start positions.
    """

    r = M.tocoo()
    # Get fragments where new chromosome starts (positions reset)
    chromstart = np.where(positions == 0)[0]
    chromend = np.append(chromstart[1:], len(positions))
    chromlen = chromend - chromstart
    # Assign a chromosome to each fragment
    chroms = np.repeat(range(len(chromlen)), chromlen)
    # Get binned positions
    positions = positions // bin_len
    frags = np.transpose(np.array([chroms, positions], dtype=np.int64))
    # Keep track of index fragments
    frag_idx = range(frags.shape[0])
    # Unique bin coordinates to create
    unique_bins = np.unique(frags, axis=0)
    # Check if some bins are missing (happens if a single
    # fragment should contain multiple bins)
    bins_jumps = (unique_bins[1:, 1] - unique_bins[:-1, 1]) - 1
    # Compute number of missing bins to add (no restriction site in bin)
    missing_bins = np.where(bins_jumps > 0)[0]
    n_missing_bins = np.sum(bins_jumps[bins_jumps > 0])
    # Compute correct number of bins to create
    n_bins = unique_bins.shape[0] + n_missing_bins
    # Initialise output fragment list (post binning)
    out_pos = np.zeros((n_bins, 1))
    row = copy.copy(r.row)
    col = copy.copy(r.col)
    # unique_bin_No: Number of bins w/ unique restriction fragments
    # actual_bin_No: Number of bins in total (including missing ones
    # sharing the same fragment)
    unique_bin_No, actual_bin_No = 0, 0
    # Match empty missing bins added with the original bin sharing
    # the same fragment
    added_bins = {}
    bin_per_frag = {}
    # Use (chr, bin) as grouping key (coord) and indices of fragments
    # belonging to current bin (bin_frags)
    for coords, bin_frags in itertools.groupby(
        frag_idx, lambda x: tuple(frags[x, :])
    ):

        bin_frags = list(bin_frags)
        first_frag, last_frag = bin_frags[0], bin_frags[-1] + 1
        # Pool row/col number by bin
        row[
            np.where((r.row >= first_frag) & (r.row < last_frag))
        ] = actual_bin_No
        col[
            np.where((r.col >= first_frag) & (r.col < last_frag))
        ] = actual_bin_No
        # Get bin position in basepair
        out_pos[actual_bin_No] = coords[1] * bin_len
        # Multiple bins to create in same fragment (rare)
        if unique_bin_No in missing_bins:
            # Number of basepairs to shift inside fragment for new bins
            curr_shift = 0
            # Subsequent bins belong to same frag as this one
            orig_bin = copy.copy(actual_bin_No)
            # Shifting bin index (introducing empty bins) for each bin in same frag
            for _ in range(bins_jumps[unique_bin_No]):
                curr_shift += bin_len
                actual_bin_No += 1
                # Remember bin coords and #bin /frag to fill contacts later
                added_bins[actual_bin_No] = orig_bin
                bin_per_frag[orig_bin] = bin_per_frag.get(unique_bin_No, 0) + 1
                out_pos[actual_bin_No] = coords[1] * bin_len + curr_shift
        unique_bin_No += 1
        actual_bin_No += 1
    row[np.where(r.row >= last_frag)] = actual_bin_No - 1
    col[np.where(r.col >= last_frag)] = actual_bin_No - 1
    # Sum data of duplicate row/col pairs
    # (i.e. combine contacts of all fragments in same bin)
    binned = coo_matrix(
        (r.data, (row, col)), shape=(actual_bin_No, actual_bin_No)
    )
    binned.sum_duplicates()
    binned.eliminate_zeros()

    return (binned, out_pos)


def bin_exact_bp_dense(M, positions, bin_len=10000):
    """Perform the kb-binning procedure with total bin lengths being exactly
    set to that of the specified input. Fragments overlapping two potential
    bins will be split and related contact counts will be divided according

    Parameters
    ----------
    to overlap proportions in each bin.
    M : 2D numpy array of ints or floats
        The Hi-C matrix to bin in dense format
    positions : numpy array of int
        List of basepair start positions of fragments bins
    bin_len : int
        Bin length in basepairs

    Returns
    -------
    2D numpy array of ints of floats :
        Binned matrix
    list :
        List of binned fragments
    """
    units = positions / bin_len
    n = len(positions)
    idx = [
        i for i in range(n - 1) if np.ceil(units[i]) < np.ceil(units[i + 1])
    ]
    m = len(idx) - 1
    N = np.zeros((m, m))
    remainders = [0] + [np.abs(units[i] - units[i + 1]) for i in range(m)]
    for i in range(m):
        N[i] = np.array(
            [
                (
                    M[idx[j] : idx[j + 1], idx[i] : idx[i + 1]].sum()
                    - remainders[j] * M[i][j]
                    + remainders[j + 1] * M[i + 1][j]
                )
                for j in range(m)
            ]
        )
    return N


def mad(M, axis=None):
    """
    Computes median absolute deviation of matrix bins sums.

    Parameters
    ----------
    M : scipy sparse coo_matrix
        Sparse matrix in COO format.

    axis: int
        Compute MAD on rows if 0, on columns if 1 or on all pixels if None.
        If axis is None, MAD is computed only on nonzero pixels.
    Returns
    -------
    float:
        MAD estimator of matrix bin sums
    """
    # Compute median on nonzero data values
    # otherwise, median is 0 if sufficiently sparse
    if axis is None:
        if issparse(M):
            r = M.tocoo()
            dist = r.data
        else:
            dist = M

    else:
        if axis < 0:
            axis += 2
        dist = np.array(M.sum(axis=axis, dtype=np.float)).flatten()

    return np.median(np.absolute(dist - np.median(dist)))


def get_good_bins(M, n_mad=2.0, s_min=None, s_max=None, symmetric=False):
    """
    Filters out bins with outstanding sums using median and MAD
    of the log transformed distribution of bin sums. Only filters
    weak outlier bins unless `symmetric` is set to True.

    Parameters
    ----------
    M : scipy sparse coo_matrix
        Input sparse matrix representing the Hi-C contact map.

    n_mad : float
        Minimum number of median absolut deviations around median in the
        bin sums distribution at which bins will be filtered out.

    s_min : float
        Optional fixed threshold value for bin sum below which bins should
        be filtered out.

    s_max: float
        Optional fixed threshold value for bin sum above which bins should
        be filtered out.
    symmetric : bool
        If set to true, filters out outliers on both sides of the distribution.
        Otherwise, only filters out bins on the left side (weak bins).
    Returns
    -------
    numpy array of bool :
        A 1D numpy array whose length is the number of bins in the matrix and
        values indicate if bins values are within the acceptable range (1)
        or considered outliers (0).
    """
    r = M.tocoo()
    with np.errstate(divide="ignore", invalid="ignore"):
        bins = sum_mat_bins(r)
        bins[bins == 0] = 1
        norm = np.log10(bins)
        median = np.median(norm)
        sigma = 1.4826 * mad(norm)

    if s_min is None:
        s_min = median - n_mad * sigma
    if s_max is None:
        s_max = median + n_mad * sigma

    if symmetric:
        filter_bins = (norm > s_min) * (norm < s_max)
    else:
        filter_bins = norm > s_min

    return filter_bins


def trim_dense(M, n_mad=3, s_min=None, s_max=None):
    """By default, return a matrix stripped of component
    vectors whose sparsity (i.e. total contact count on a
    single column or row) deviates more than specified number
    of standard deviations from the mean. Boolean variables
    s_min and s_max act as absolute fixed values which override
    such behaviour when specified.

    Parameters
    ----------
    M : 2D numpy array of floats
        Dense Hi-C contact matrix
    n_mad : int
        Minimum number of standard deviation by which a the sum of
        contacts in a component vector must deviate from the mean
        to be trimmed.
    s_min : float
        Fixed minimum value below which the component vectors will
        be trimmed.
    s_max : float
        Fixed maximum value above which the component vectors will
        be trimmed.

    Returns
    -------
    numpy 2D array of floats :
        The input matrix, stripped of outlier component vectors.
    """

    S = coo_matrix(M)
    S_out = trim_sparse(S, n_mad=n_mad, s_min=s_min, s_max=s_max)
    M_out = S_out.todense()
    return M_out


def trim_sparse(M, n_mad=3, s_min=None, s_max=None):
    """Apply the trimming procedure to a sparse matrix.

    Parameters
    ----------
    M : scipy.sparse.coo_matrix
        Sparse Hi-C contact map
    n_mad : int
        Minimum number of median absolute deviations by which a the sum of
        contacts in a component vector must deviate from the median
        to be trimmed.
    s_min : float
        Fixed minimum value below which the component vectors will
        be trimmed.
    s_max : float
        Fixed maximum value above which the component vectors will
        be trimmed.

    Returns
    -------
     scipy coo_matrix of floats :
        The input sparse matrix, stripped of outlier component vectors.
    """
    r = M.tocoo()
    f = get_good_bins(M, n_mad, s_min, s_max)
    miss_bins = np.cumsum(1 - f)
    # Mapping pre- and post- trimming indices of bins, post = -1 means delete
    # Note: There is probably a more efficient way than a dictionary for that
    miss_map = {old: old - offset for old, offset in enumerate(miss_bins)}
    # Indices of cells that will be kept
    indices = np.where(f[r.row] & f[r.col])
    # Remove sparse rows and shift indices accordingly
    rows = [miss_map[i] for i in r.row[indices]]
    cols = [miss_map[j] for j in r.col[indices]]
    data = r.data[indices]
    size = max(max(rows, default=-1), max(cols, default=-1)) + 1
    N = coo_matrix((data, (rows, cols)), shape=(size, size))
    return N


def normalize_dense(M, norm="SCN", order=1, iterations=40):
    """Apply one of the many normalization types to input dense
    matrix. Will also apply any callable norms such as a user-made
    or a lambda function.
    NOTE: Legacy function for dense maps

    Parameters
    ----------
    M : 2D numpy array of floats
    norm : str
        Normalization procedure to use. Can be one of "SCN",
        "mirnylib", "frag" or "global". Can also be a user-
        defined function.
    order : int
        Defines the type of vector norm to use. See numpy.linalg.norm
        for details.
    iterations : int
        Iterations parameter when using an iterative normalization
        procedure.

    Returns
    -------
    2D numpy array of floats :
        Normalized dense matrix.
    """

    s = np.array(M, np.float64)
    floatorder = np.float64(order)

    if norm == "SCN":
        for _ in range(0, iterations):

            sumrows = s.sum(axis=1)
            maskrows = (sumrows != 0)[:, None] * (sumrows != 0)[None, :]
            sums_row = sumrows[:, None] * np.ones(sumrows.shape)[None, :]
            s[maskrows] = 1.0 * s[maskrows] / sums_row[maskrows]

            sumcols = s.sum(axis=0)
            maskcols = (sumcols != 0)[:, None] * (sumcols != 0)[None, :]
            sums_col = sumcols[None, :] * np.ones(sumcols.shape)[:, None]
            s[maskcols] = 1.0 * s[maskcols] / sums_col[maskcols]

    elif norm == "mirnylib":
        try:
            from mirnylib import numutils as ntls

            s = ntls.iterativeCorrection(s, iterations)[0]
        except ImportError as e:
            print(str(e))
            print("I can't find mirnylib.")
            print(
                "Please install it from "
                "https://bitbucket.org/mirnylab/mirnylib"
            )
            print("I will use default norm as fallback.")
            return normalize_dense(M, order=order, iterations=iterations)

    elif norm == "frag":
        for _ in range(1, iterations):
            s_norm_x = np.linalg.norm(s, ord=floatorder, axis=0)
            s_norm_y = np.linalg.norm(s, ord=floatorder, axis=1)
            s_norm = np.tensordot(s_norm_x, s_norm_y, axes=0)
            s[s_norm != 0] = 1.0 * s[s_norm != 0] / s_norm[s_norm != 0]

    elif callable(norm):
        s = norm(M)

    else:
        raise Exception(
            'Unknown norm, please specify one of ("mirnylib", "SCN", "frag")'
        )

    return (s + s.T) / 2


def normalize_sparse(M, norm="SCN", iterations=40, n_mad=3.0):
    """Applies a normalization type to a sparse matrix.

    Parameters
    ----------
    M : scipy.sparse.csr_matrix of floats
    norm : str or callable
        Normalization procedure to use. Can be one of "SCN" or
        "ICE". Can also be a user-defined function.
    iterations : int
        Iterations parameter when using an iterative normalization
        procedure.
    n_mad : float
        Maximum number of median absolute deviations of bin sums to allow for 
        including bins in the normalization procedure. Bins more than `n_mad`
        mads below the median are excluded. Bins excluded from normalisation
        are set to 0.

    Returns
    -------
    scipy.sparse.csr_matrix of floats :
        Normalized sparse matrix.
    """
    # Making full symmetric matrix if not symmetric already (e.g. upper triangle)
    r = M.astype(np.float64)
    good_bins = get_good_bins(M, n_mad=n_mad)
    # Set values in non detectable bins to 0
    # For faster masking of bins, mask bins using dot product with an identity
    # matrix where bad bins have been masked on the diagonal
    # E.g. if removing the second bin (row and column):
    # 1 0 0     9 6 5     1 0 0     9 0 5
    # 0 0 0  X  6 8 7  X  0 0 0  =  0 0 0
    # 0 0 1     6 7 8     0 0 1     6 0 8
    mask_mat = sparse.eye(r.shape[0])
    mask_mat.data[0][~good_bins] = 0
    r = mask_mat.dot(r).dot(mask_mat)
    r = coo_matrix(r)
    r.eliminate_zeros()
    if norm == "ICE":
        # Row and col indices of each nonzero value in matrix
        row_indices, col_indices = r.nonzero()
        for _ in range(iterations):
            # Symmetric matrix: rows and cols have identical sums
            bin_sums = sum_mat_bins(r)
            # Normalize bin sums by the median sum of detectable bins for stability
            bin_sums /= np.median(bin_sums[good_bins])
            # Divide each nonzero value by the product of the sums of
            # their respective rows and columns.
            r.data /= bin_sums[row_indices] * bin_sums[col_indices]
        bin_sums = sum_mat_bins(r)
        # Scale to 1
        r.data = r.data * (1 / np.median(bin_sums[good_bins]))
    elif norm == "SCN":
        # Similar to ICE, but division is done sequentially by row and then column
        # sums instead of using product.
        row_indices, col_indices = r.nonzero()
        for i in range(iterations):
            bin_sums = sum_mat_bins(r)
            r.data /= bin_sums[row_indices]
            bin_sums = sum_mat_bins(r)
            r.data /= bin_sums[col_indices]

    elif callable(norm):
        r = norm(M)

    else:
        raise Exception('Unknown norm, please specify one of ("ICE", "SCN")')

    return r


def sum_mat_bins(mat):
    """
    Compute the sum of matrices bins (i.e. rows or columns) using
    only the upper triangle, assuming symmetrical matrices.

    Parameters
    ----------
    mat : scipy.sparse.csr_matrix
        Contact map in sparse format, either in upper triangle or
        full matrix.
    
    Returns
    -------
    numpy.array :
        1D array of bin sums.
    """
    # Equivalaent to row or col sum on a full matrix
    # Note: mat.sum returns a 'matrix' object. A1 extracts the 1D flat array
    # from the matrix
    return mat.sum(axis=0).A1 + mat.sum(axis=1).A1 - mat.diagonal(0)


def GC_partial(portion: str):
    """Manually compute GC content percentage in a DNA string, taking
    ambiguous values into account (according to standard IUPAC notation).

    Parameters
    ----------
    portion : str
        DNA sequence on which GC content is computed.

    Returns
    -------
    float :
        The percentage of GC in the input string.
    """

    sequence_count = collections.Counter(portion)
    gc = (
        sum([sequence_count[i] for i in "gGcCsS"])
        + sum([sequence_count[i] for i in "DdHh"]) / 3.0
        + 2 * sum([sequence_count[i] for i in "VvBb"]) / 3.0
        + sum([sequence_count[i] for i in "NnYyRrKkMm"]) / 2.0
    ) / len(portion)
    return 0 or 100 * gc


def GC_wide(genome: str, window=1000):
    """Compute GC across a window of given length.

    Parameters
    ----------
    genome : str
        The genome on which GC content will be computed.
    window : int
        The window size in which GC content is measured.

    Note
    ----
    Biopython is required.
    """

    from Bio import SeqIO

    with open(genome) as handle:
        sequence = "".join(
            [str(record.seq) for record in SeqIO.parse(handle, "fasta")]
        )

    n = len(sequence)
    for i in range(0, n, window):
        portion = sequence[i : min(i + window, n)]
        yield GC_partial(portion)


def split_genome(genome, chunk_size=10000):
    """Split genome into chunks of fixed size (save the last one).
    """

    chunks = []
    from Bio import SeqIO

    with open(genome) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = record.seq
            n = len(sequence)
            chunks += [
                str(sequence[i : min(i + chunk_size, n)])
                for i in range(0, n, chunk_size)
            ]
    return np.array(chunks)


def directional(M, window=None, circ=False, extrapolate=True, log=True):
    """From a symmetrical matrix M of size n, return a vector d whose each
    component d[i] is a T-test of two samples represented by vectors of size
    window on either side of the i-th pixel on the diagonal. Edge elements may
    be extrapolated based on the vector size reduction, except in the case of
    circular genomes. If they aren't, d will be of size n - 2*(window-1)
    instead of n.
    """

    # Sanity checks
    if not type(M) is np.ndarray:
        M = np.array(M)

    if M.shape[0] != M.shape[1]:
        raise ValueError("Matrix is not square.")

    try:
        n = min(M.shape)
    except AttributeError:
        n = M.size

    # Default window argument
    if window is None:
        window = max(n // 100, 5)

    if window >= n:
        raise ValueError("Please choose a smaller window size.")

    try:
        from scipy.stats import ttest_rel
    except ImportError as e:
        print("Scipy not")
        print(str(e))
        raise
    if log:
        N = np.zeros((n, n))
        N[M > 0] = np.log(M[M > 0])
    else:
        N = M

    if circ:
        d = [
            ttest_rel(
                np.array(list(N[i, i - window :]) + list(N[i, :i])),
                N[i, i : i + window],
            )[0]
            for i in range(window)
        ]
    elif extrapolate:
        d = [ttest_rel(N[i, 0:i], N[i, i : 2 * i])[0] for i in range(window)]
    else:
        d = []

    d += [
        ttest_rel(N[i, i - window : i], N[i, i : i + window])[0]
        for i in range(window, n - window)
    ]

    if circ:
        d += [
            ttest_rel(
                N[i, i - window : i],
                np.array((list(N[i, : i - n + window]) + list(N[i, i:]))),
            )[0]
            for i in range(n - window, n)
        ]
    elif extrapolate:
        d += [
            ttest_rel(
                N[i, i - window : i],
                (np.array(list(N[i, i:]) + list(N[i, : window - (n - i)]))),
            )[0]
            for i in range(n - window, n)
        ]

    return d


def domainogram(M, window=None, circ=False, extrapolate=True):
    """From a symmetrical matrix M of size n, return a vector d whose each
    component d[i] is the total sum of a square of 2*window+1 size centered on
    the i-th main diagonal element. Edge elements may be extrapolated based on
    the square size reduction (i.e. for window = 4, the first component will be
    equal to the first diagonal pixel multiplied by 81, the second one will be
    equal to the first 2x2 square on the diagonal multiplied by 81/4, etc.),
    except in the case of circular genomes. If they aren't, d will be of size
    n - 2*(window-1) instead of n.
    """

    # Sanity checks
    if not type(M) is np.ndarray:
        M = np.array(M)

    if M.shape[0] != M.shape[1]:
        raise ValueError("Matrix is not square.")

    try:
        n = min(M.shape)
    except AttributeError:
        n = M.size

    # Default window argument
    if window is None:
        window = max(n // 100, 5)

    if window >= n:
        raise ValueError("Please choose a smaller window size.")

    if circ:
        d = [
            (
                np.sum(M[-i + window :, -i + window :])
                + np.sum(M[: i - window + 1, : i - window + 1])
                for i in range(window)
            )
        ]
    elif extrapolate:
        d = [
            (
                np.sum(M[0 : 2 * i + 1, 0 : 2 * i + 1])
                * ((2 * window + 1) ** 2.0)
                / ((2 * i + 1) ** 2.0)
            )
            for i in range(window)
        ]
    else:
        d = []

    d += [
        np.sum(M[i - window : i + window + 1, i - window : i + window + 1])
        for i in range(window, n - window)
    ]

    if circ:
        d += [
            M[i:, i:].sum() + M[: n - i, n - i].sum()
            for i in range(n - window, n)
        ]
    elif extrapolate:
        d += [
            M[i - window :, i - window :].sum()
            * ((2 * window + 1) ** 2.0)
            / ((2 * (n - i) + 1) ** 2.0)
            for i in range(n - window, n)
        ]

    return np.array(d)


def from_structure(structure):
    """Return contact data from a 3D structure (in pdb format).
    """

    try:
        from Bio import PDB

        if isinstance(structure, str):
            p = PDB.PDBParser()
            structure = p.get_structure("S", structure)
        if isinstance(structure, PDB.Structure.Structure):
            for _ in structure.get_chains():
                atoms = [
                    np.array(atom.get_coord())
                    for atom in structure.get_atoms()
                ]
    except ImportError:
        print("Biopython not found.")
        raise

    atoms = np.array(structure)
    try:
        import scipy

        D = scipy.spatial.distance.pdist(atoms, "euclidean")
        D = scipy.spatial.distance.squareform(D)
    except ImportError:
        print("Scipy not found.")
        raise
    m = np.max(1 / D[D != 0])
    M = np.zeros(D.shape)
    M[D != 0] = 1 / D[D != 0]
    M[D == 0] = m
    return M


def largest_connected_component(matrix):
    """Compute the adjacency matrix of the largest connected component of the
    graph whose input matrix is adjacent.
    """

    try:
        import scipy.sparse

        n, components = scipy.sparse.csgraph.connected_components(
            matrix, directed=False
        )
        print("I found " + str(n) + " connected components.")
        component_dist = collections.Counter(components)
        print("Distribution of components: " + str(component_dist))
        most_common, _ = component_dist.most_common(1)[0]
        ilcc = components == most_common
        return matrix[:, ilcc][ilcc]

    except ImportError as e:
        print("I couldn't find scipy which is needed for graph routines.")
        print(str(e))
        print("Returning input matrix as fallback.")
        return matrix


def to_structure(matrix, alpha=1):
    """Compute best matching 3D genome structure from underlying input matrix
    using ShRec3D-derived method from Lesne et al., 2014.

    Link: https://www.ncbi.nlm.nih.gov/pubmed/25240436

    The method performs two steps: first compute distance matrix by treating
    contact data as an adjacency graph (of weights equal to a power law
    function of the data), then embed the resulting distance matrix into
    3D space.

    The alpha parameter influences the weighting of contacts: if alpha < 1
    long-range interactions are prioritized; if alpha >> 1 short-range
    interactions have more weight wahen computing the distance matrix.
    """

    connected = largest_connected_component(matrix)
    distances = to_distance(connected, alpha)
    n, m = connected.shape
    bary = np.sum(np.triu(distances, 1)) / (n ** 2)  # barycenters
    d = np.array(np.sum(distances ** 2, 0) / n - bary)  # distances to origin
    gram = np.array(
        [
            (d[i] + d[j] - distances[i][j] ** 2) / 2
            for i, j in itertools.product(range(n), range(m))
        ]
    ).reshape(n, m)
    normalized = gram / np.linalg.norm(gram, "fro")

    try:
        symmetric = np.array(
            (normalized + normalized.T) / 2, dtype=np.longfloat
        )  # just in case
    except AttributeError:
        symmetric = np.array((normalized + normalized.T) / 2)

    from scipy import linalg

    eigen_values, eigen_vectors = linalg.eigh(symmetric)
    if not (eigen_values >= 0).all():
        warnings.warn("Negative eigen values were found.")
    idx = eigen_values.argsort()[-3:][::-1]
    values = eigen_values[idx]
    vectors = eigen_vectors[:, idx]
    coordinates = vectors * np.sqrt(values)
    return coordinates


def get_missing_bins(original, trimmed):
    """Retrieve indices of a trimmed matrix with respect to the original matrix.
    Fairly fast but is only correct if diagonal values are different, which is
    always the case in practice.
    """

    original_diag = np.diag(original)
    trimmed_diag = np.diag(trimmed)
    index = []
    m = min(original.shape)
    for j in range(min(trimmed.shape)):
        k = 0
        while original_diag[j + k] != trimmed_diag[j] and k < 2 * m:
            k += 1
        index.append(k + j)
    return np.array(index)


def to_pdb(
    structure,
    filename,
    contigs=None,
    annotations=None,
    indices=None,
    special_bins=None,
):
    """From a structure (or matrix) generate the corresponding pdb file
    representing each chain as a contig/chromosome and filling the occupancy
    field with a custom annotation. If the matrix has been trimmed somewhat,
    remaining indices may be specified.
    """

    n = len(structure)
    letters = (
        string.ascii_uppercase
        + string.ascii_lowercase
        + string.digits
        + string.punctuation
    ) * int(n / 94 + 1)
    if contigs is None:
        contigs = np.ones(n + 1)
    if annotations is None:
        annotations = np.zeros(n + 1)
    if indices is None:
        indices = range(n + 1)
    if special_bins is None:
        special_bins = np.zeros(n + 1, dtype=int)

    structure_shapes_match = structure.shape[0] == structure.shape[1]
    if isinstance(structure, np.ndarray) and structure_shapes_match:
        structure = to_structure(structure)

    X, Y, Z = (structure[:, i] for i in range(3))
    Xmax, Ymax, Zmax = (np.max(np.abs(Xi)) for Xi in (X, Y, Z))
    X *= 100.0 / Xmax
    Y *= 100.0 / Ymax
    Z *= 100.0 / Zmax

    X = np.around(X, 3)
    Y = np.around(Y, 3)
    Z = np.around(Z, 3)

    reference = ["OW", "OW", "CE", "TE", "tR"]
    with open(filename, "w") as f:
        for i in range(1, n):

            line = "ATOM"  # 1-4 "ATOM"
            line += "  "  # 5-6 unused
            line += str(i).rjust(5)  # 7-11 atom serial number
            line += " "  # 12 unused
            line += reference[special_bins[i]].rjust(4)  # 13-16 atom name
            line += " "  # 17 alternate location indicator
            line += "SOL"  # 18-20 residue name
            line += " "  # 21 unused
            line += letters[
                int(contigs[indices[i]] - 1)
            ]  # 22 chain identifier
            line += str(i).rjust(4)  # 23-26 residue sequence number
            line += " "  # 27 code for insertion of residues
            line += "   "  # 28-30 unused
            line += str(X[i]).rjust(8)  # 31-38 X orthogonal Å coordinate
            line += str(Y[i]).rjust(8)  # 39-46 Y orthogonal Å coordinate
            line += str(Z[i]).rjust(8)  # 47-54 Z orthogonal Å coordinate
            line += "1.00".rjust(6)  # 55-60 Occupancy
            # 61-66 Temperature factor
            line += str(annotations[i - 1]).rjust(6)
            line += "      "  # 67-72 unused
            line += "    "  # 73-76 segment identifier
            line += "O".rjust(2)  # 77-78 element symbol
            line += "\n"
            f.write(line)


def matrix_to_pdb(
    matrix,
    filename,
    contigs=None,
    annotations=None,
    indices=None,
    special_bins=None,
    alpha=1,
):
    """Convert a matrix to a PDB file, shortcutting the intermediary generated
    structure.
    """
    to_pdb(
        to_structure(matrix, alpha=alpha),
        filename=filename,
        contigs=contigs,
        annotations=annotations,
        indices=indices,
        special_bins=special_bins,
    )


def to_distance(matrix, alpha=1):
    """Compute distance matrix from contact data by applying a negative power
    law (alpha) to its nonzero pixels, then interpolating on the zeroes using a
    shortest-path algorithm.
    """
    matrix = np.array(matrix)
    try:
        import scipy.sparse
    except ImportError as e:
        print("Scipy not found.")
        print(str(e))
        raise

    if callable(alpha):
        distance_function = alpha
    else:
        try:
            a = np.float64(alpha)

            def distance_function(x):
                return 1 / (x ** (1 / a))

        except TypeError:
            print("Alpha parameter must be callable or an array-like")
            raise

    if hasattr(matrix, "getformat"):
        distances = scipy.sparse.coo_matrix(matrix)
        distances.data = distance_function(distances.data)
    else:
        distances = np.zeros(matrix.shape)
        distances[matrix != 0] = distance_function(1 / matrix[matrix != 0])

    return scipy.sparse.csgraph.floyd_warshall(distances, directed=False)


def distance_to_contact(D, alpha=1):
    """Compute contact matrix from input distance matrix. Distance values of
    zeroes are given the largest contact count otherwise inferred non-zero
    distance values.
    """

    if callable(alpha):
        distance_function = alpha
    else:
        try:
            a = np.float64(alpha)

            def distance_function(x):
                return 1 / (x ** (1 / a))

        except TypeError:
            print("Alpha parameter must be callable or an array-like")
            raise
        except ZeroDivisionError:
            raise ValueError("Alpha parameter must be non-zero")

    m = np.max(distance_function(D[D != 0]))
    M = np.zeros(D.shape)
    M[D != 0] = distance_function(D[D != 0])
    M[D == 0] = m
    return M


def subsample_contacts(M, n_contacts):
    """Bootstrap sampling of contacts in a sparse Hi-C map.

    Parameters
    ----------
    M : scipy.sparse.coo_matrix
        The input Hi-C contact map in sparse format.
    n_contacts : float
        The number of contacts to be sampled if larger than one.
        The proportion of contacts to be sampled if between 0 and 1.
    Returns
    -------
    scipy.sparse.coo_matrix
        A new matrix with a fraction of the original contacts.
    """
    try:
        if n_contacts <= 1 and n_contacts > 0:
            n_contacts *= M.data.sum()
        elif n_contacts < 0:
            logger.error("n_contacts must be strictly positive")
    except ValueError as e:
        logger.error("n_contacts must be a float")
        raise e
    S = M.data.copy()
    # Match cell idx to cumulative number of contacts
    cum_counts = np.cumsum(S)
    # Total number of contacts to sample
    tot_contacts = int(cum_counts[-1])

    # Sample desired number of contacts from the range(0, n_contacts) array
    sampled_contacts = np.random.choice(
        int(tot_contacts), size=int(n_contacts), replace=False
    )

    # Get indices of sampled contacts in the cum_counts array
    idx = np.searchsorted(cum_counts, sampled_contacts, side="right")

    # Bin those indices to the same dimensions as matrix data to get counts
    sampled_counts = np.bincount(idx, minlength=S.shape[0])

    # Get nonzero values to build new sparse matrix
    nnz_mask = sampled_counts > 0
    sampled_counts = sampled_counts[nnz_mask].astype(np.float64)
    sampled_rows = M.row[nnz_mask]
    sampled_cols = M.col[nnz_mask]

    return coo_matrix(
        (sampled_counts, (sampled_rows, sampled_cols)),
        shape=(M.shape[0], M.shape[1]),
    )


def shortest_path_interpolation(matrix, alpha=1, strict=True):
    """Perform interpolation on a matrix's data by using ShRec's shortest-path
    procedure backwards and forwards. This replaces zeroes with corresponding
    shortest-path based counts and may have the additional effect of 'blurring'
    the matrix somewhat. If strict is set to True, only zeroes are replaced
    this way.

    Also known as Boost-Hi-C (https://www.ncbi.nlm.nih.gov/pubmed/30615061)
    """
    matrix = np.array(matrix, np.float64)
    contacts = distance_to_contact(
        to_distance(matrix, alpha=alpha), alpha=alpha
    )
    if not strict:
        return contacts
    else:
        M = np.copy(matrix)
        M[matrix == 0] = contacts[matrix == 0]
        return M


def pdb_to_structure(filename):
    """Import a structure object from a PDB file.
    """

    try:
        from Bio import PDB
    except ImportError:
        print("I can't import Biopython which is needed to handle PDB files.")
        raise
    p = PDB.PDBParser()
    structure = p.get_structure("S", filename)
    for _ in structure.get_chains():
        atoms = [np.array(atom.get_coord()) for atom in structure.get_atoms()]
    return atoms


def noise(matrix):
    """Just a quick function to make a matrix noisy using a standard Poisson
    distribution (contacts are treated as rare events).
    """

    D = shortest_path_interpolation(matrix, strict=True)
    return np.random.poisson(lam=D)


def flatten_positions_to_contigs(positions):
    """Flattens and converts a positions array to a contigs array, if
    applicable.
    """

    if isinstance(positions, np.ndarray):
        flattened_positions = positions.flatten()
    else:
        try:
            flattened_positions = np.array(
                [pos for contig in positions for pos in contig]
            )
        except TypeError:
            flattened_positions = np.array(positions)

    if (np.diff(positions) == 0).any() and not (0 in set(positions)):
        warnings.warn("I detected identical consecutive nonzero values.")
        return positions

    n = len(flattened_positions)
    contigs = np.ones(n)
    counter = 0
    for i in range(1, n):
        if positions[i] == 0:
            counter += 1
            contigs[i] += counter
        else:
            contigs[i] = contigs[i - 1]
    return contigs


def simple_distance_diagonal_law(matrix, circular=False):
    if not circular:
        n = len(matrix)
        return np.array([np.average(np.diagonal(matrix, j)) for j in range(n)])
    else:
        n = len(matrix)
        return [
            (
                np.average(np.diagonal(matrix, j))
                + np.average(np.diagonal(matrix, n - j))
            )
            / 2.0
            for j in range(n)
        ]


def distance_diagonal_law(matrix, positions=None, circular=False):
    """Compute a distance law trend using the contact averages of equal
    distances. Specific positions can be supplied if needed.
    """

    n = min(matrix.shape)
    if positions is None:
        return simple_distance_diagonal_law(matrix, circular=circular)
    else:
        contigs = positions_to_contigs(positions)

    def is_intra(i, j):
        return contigs[i] == contigs[j]

    max_intra_distance = max((len(contigs == u) for u in set(contigs)))

    intra_contacts = []
    inter_contacts = [
        np.average(np.diagonal(matrix, j))
        for j in range(max_intra_distance, n)
    ]
    for j in range(max_intra_distance):
        D = np.diagonal(matrix, j)
        for i in range(len(D)):
            diagonal_intra = []
            if is_intra(i, j):
                diagonal_intra.append(D[i])
        #            else:
        #                diagonal_inter.append(D[i])
        #        inter_contacts.append(np.average(np.array(diagonal_inter)))
        intra_contacts.append(np.average(np.array(diagonal_intra)))

    intra_contacts.extend(inter_contacts)

    return [positions, np.array(intra_contacts)]


def rippe_parameters(matrix, positions, lengths=None, init=None, circ=False):
    """Estimate parameters from the model described in Rippe et al., 2001.
    """

    n, _ = matrix.shape

    if lengths is None:
        lengths = np.abs(np.diff(positions))

    measurements, bins = [], []
    for i in range(n):
        for j in range(1, i):
            mean_length = (lengths[i] + lengths[j]) / 2.0
            if positions[i] < positions[j]:
                d = (
                    (positions[j] - positions[i] - lengths[i]) + mean_length
                ) / 1000.0
            else:
                d = (
                    (positions[i] - positions[j] - lengths[j]) + mean_length
                ) / 1000.0

            bins.append(np.abs(d))
            measurements.append(matrix[i, j])
    parameters = estimate_param_rippe(measurements, bins, init=init, circ=circ)
    print(parameters)
    return parameters[0]


def estimate_param_rippe(measurements, bins, init=None, circ=False):
    """Perform least square optimization needed for the rippe_parameters function.
    """

    # Init values
    DEFAULT_INIT_RIPPE_PARAMETERS = [1.0, 9.6, -1.5]
    d = 3.0

    def log_residuals(p, y, x):
        kuhn, lm, slope, A = p
        rippe = (
            np.log(A)
            + np.log(0.53)
            - 3 * np.log(kuhn)
            + slope * (np.log(lm * x) - np.log(kuhn))
            + (d - 2) / ((np.power((lm * x / kuhn), 2) + d))
        )
        err = y - rippe

        return err

    def peval(x, param):

        if circ:
            l_cont = x.max()
            n = param[1] * x / param[0]
            n0 = param[1] * x[0] / param[0]
            n_l = param[1] * l_cont / param[0]
            s = n * (n_l - n) / n_l
            s0 = n0 * (n_l - n0) / n_l
            norm_lin = param[3] * (
                0.53
                * (param[0] ** -3.0)
                * np.power(n0, (param[2]))
                * np.exp((d - 2) / ((np.power(n0, 2) + d)))
            )

            norm_circ = param[3] * (
                0.53
                * (param[0] ** -3.0)
                * np.power(s0, (param[2]))
                * np.exp((d - 2) / ((np.power(s0, 2) + d)))
            )

            rippe = (
                param[3]
                * (
                    0.53
                    * (param[0] ** -3.0)
                    * np.power(s, (param[2]))
                    * np.exp((d - 2) / ((np.power(s, 2) + d)))
                )
                * norm_lin
                / norm_circ
            )

        else:

            rippe = param[3] * (
                0.53
                * (param[0] ** -3.0)
                * np.power((param[1] * x / param[0]), (param[2]))
                * np.exp(
                    (d - 2) / ((np.power((param[1] * x / param[0]), 2) + d))
                )
            )

        return rippe

    if init is None:
        init = DEFAULT_INIT_RIPPE_PARAMETERS

    A = np.sum(measurements)

    p0 = (p for p in init), A
    from scipy.optimize import leastsq

    plsq = leastsq(log_residuals, p0, args=(np.log(measurements), bins))

    y_estim = peval(bins, plsq[0])
    kuhn_x, lm_x, slope_x, A_x = plsq[0]
    plsq_out = [kuhn_x, lm_x, slope_x, d, A_x]

    np_plsq = np.array(plsq_out)

    if np.any(np.isnan(np_plsq)) or slope_x >= 0:
        warnings.warn("Problem in parameters estimation")
        plsq_out = p0

    return plsq_out, y_estim


def null_model(
    matrix,
    positions=None,
    lengths=None,
    model="uniform",
    noisy=False,
    circ=False,
    sparsity=False,
):
    """Attempt to compute a 'null model' of the matrix given a model
    to base itself on.
    """

    n, m = matrix.shape
    positions_supplied = True
    if positions is None:
        positions = range(n)
        positions_supplied = False
    if lengths is None:
        lengths = np.diff(positions)

    N = np.copy(matrix)

    contigs = np.array(positions_to_contigs(positions))

    def is_inter(i, j):
        return contigs[i] != contigs[j]

    diagonal = np.diag(matrix)

    if model == "uniform":
        if positions_supplied:
            trans_contacts = np.array(
                [
                    matrix[i, j]
                    for i, j in itertools.product(range(n), range(m))
                    if is_inter(i, j)
                ]
            )
            mean_trans_contacts = np.average(trans_contacts)
        else:
            mean_trans_contacts = np.average(matrix) - diagonal / len(diagonal)

        N = np.random.poisson(lam=mean_trans_contacts, size=(n, m))
        np.fill_diagonal(N, diagonal)

    elif model == "distance":
        distances = distance_diagonal_law(matrix, positions)
        N = np.array(
            [
                [distances[min(abs(i - j), n)] for i in range(n)]
                for j in range(n)
            ]
        )

    elif model == "rippe":

        trans_contacts = np.array(
            [
                matrix[i, j]
                for i, j in itertools.product(range(n), range(m))
                if is_inter(i, j)
            ]
        )
        mean_trans_contacts = np.average(trans_contacts)
        kuhn, lm, slope, d, A = rippe_parameters(matrix, positions, circ=circ)

        def jc(s, frag):
            dist = s - circ * (s ** 2) / lengths[frag]
            computed_contacts = (
                0.53
                * A
                * (kuhn ** (-3.0))
                * (dist ** slope)
                * np.exp((d - 2) / (dist + d))
            )
            return np.maximum(computed_contacts, mean_trans_contacts)

        for i in range(n):
            for j in range(n):
                if not is_inter(i, j) and i != j:
                    posi, posj = positions[i], positions[j]
                    N[i, j] = jc(np.abs(posi - posj) * lm / kuhn, frag=j)
                else:
                    N[i, j] = mean_trans_contacts

    if sparsity:
        contact_sum = matrix.sum(axis=0)
        n = len(contact_sum)
        try:
            from Bio.Statistics import lowess

            trend = lowess.lowess(
                np.array(range(n), dtype=np.float64), contact_sum, f=0.03
            )
        except ImportError:
            expected_size = int(np.amax(contact_sum) / np.average(contact_sum))
            w = min(max(expected_size, 20), 100)
            trend = np.array(
                [np.average(contact_sum[i : min(i + w, n)]) for i in range(n)]
            )

        cov_score = np.sqrt((trend - np.average(trend)) / np.std(trend))

        N = ((N * cov_score).T) * cov_score

    if noisy:
        if callable(noisy):
            noise_function = noisy
        return noise_function(N)
    else:
        return N


def model_norm(
    matrix, positions=None, lengths=None, model="uniform", circ=False
):

    N = null_model(
        matrix,
        positions,
        lengths,
        model,
        noisy=False,
        circ=circ,
        sparsity=True,
    )
    return matrix / shortest_path_interpolation(N, strict=True)


def trim_structure(struct, filtering="cube", n=2):
    """Remove outlier 'atoms' (aka bins) from a structure.
    """

    X, Y, Z = (struct[:, i] for i in range(3))

    if filtering == "sphere":
        R = (np.std(X) ** 2 + np.std(Y) ** 2 + np.std(Z) ** 2) * (n ** 2)
        f = (X - np.mean(X)) ** 2 + (Y - np.mean(Y)) ** 2 + (
            Z - np.mean(Z)
        ) ** 2 < R

    if filtering == "cube":
        R = min(np.std(X), np.std(Y), np.std(Z)) * n
        f = np.ones(len(X))
        for C in (X, Y, Z):
            f *= np.abs(C - np.mean(C)) < R

    if filtering == "percentile":
        f = np.ones(len(X))
        for C in (X, Y, Z):
            f *= np.abs(C - np.mean(C)) < np.percentile(
                np.abs(C - np.mean(C)), n
            )

    return np.array([X[f], Y[f], Z[f]])


def scalogram(M, circ=False, max_range=False):
    """Computes so-called 'scalograms' used to easily
    visualize contacts at different distance scales.
    Edge cases have been painstakingly taken
    care of.

    Parameters
    ----------
    M1 : array_like
        The input contact map
    circ : bool
        Whether the contact map's reference genome is
        circular. Default is False.
    max_range : bool or int
        The maximum scale to be computed on the matrix.
        Default is False, which means the maximum possible
        range (len(M) // 2) will be taken.

    Returns
    -------
    N : array_like
        The output scalogram. Values that can't be computed
        due to edge issues, or being beyond max_range will
        be zero. In a non-circular matrix, this will result
        with a 'cone-shaped' contact map.
    """

    # Sanity checks
    if not type(M) is np.ndarray:
        M = np.array(M)

    if M.shape[0] != M.shape[1]:
        raise ValueError("Matrix is not square.")

    try:
        n = min(M.shape)
    except AttributeError:
        n = len(M)
    N = np.zeros(M.shape)
    if not max_range:
        max_range = M.shape[0] // 2
    for i in range(n):
        for j in range(max_range):
            if i + j < n and i >= j:
                N[i, j] = M[i, i - j : i + j + 1].sum()
            elif not circ and i + j < n and i < j:
                N[i, j] = M[i, i : i + j + 1].sum() * 2
            elif not circ and i + j >= n:
                N[i, j] = M[i, i - j : i + 1].sum() * 2
            elif circ and i + j < n and i < j:
                N[i, j] = M[i, i - j :].sum() + M[i, : i + j + 1].sum()
            elif circ and i >= j and i + j >= n:
                N[i, j] = M[i, i - j :].sum() + M[i, : i + j - n + 1].sum()
            elif circ and i < j and i + j >= n:
                N[i, j] = (
                    M[i, i - j :].sum()
                    + M[i, :].sum()
                    + M[i, : i + j - n + 1].sum()
                )
    return N


def asd(M1, M2):
    """Compute a Fourier transform based distance
    between two matrices.

    Inspired from Galiez et al., 2015
    (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4535829/)

    Parameters
    ----------
    M1 : array_like
        The first (normalized) input matrix.
    M2 : array_like
        The second (normalized) input matrix

    Returns
    -------
    asd : numpy.float64
        The matrix distance
    """

    from scipy.fftpack import fft2

    spectra1 = np.abs(fft2(M1))
    spectra2 = np.abs(fft2(M2))

    return np.linalg.norm(spectra2 - spectra1)


def compartments(M, normalize=True):
    """A/B compartment analysis

    Perform a PCA-based A/B compartment analysis on a normalized, single
    chromosome contact map. The results are two vectors whose values (negative
    or positive) should presumably correlate with the presence of 'active'
    vs. 'inert' chromatin.

    Parameters
    ----------
    M : array_like
        The input, normalized contact map. Must be a single chromosome.
    normalize : bool
        Whether to normalize the matrix beforehand.

    Returns
    -------
    PC1 : numpy.ndarray
        A vector representing the first component.
    PC2 : numpy.ndarray
        A vector representing the second component.
    """

    n = M.shape[0]
    if type(M) is not np.ndarray:
        M = np.array(M)

    if M.shape[0] != M.shape[1]:
        raise ValueError("Matrix is not square.")

    if normalize:
        N = normalize_dense(M)
    else:
        N = np.copy(M)
    # Computation of genomic distance law matrice:
    dist_mat = np.zeros((n, n))
    _, dist_vals = distance_law_from_mat(N, log_bins=False)
    for i in range(n):
        for j in range(n):
            dist_mat[i, j] = dist_vals[abs(j - i)]

    N /= dist_mat
    # Computation of the correlation matrice:
    N = np.corrcoef(N)
    N[np.isnan(N)] = 0.0

    # Computation of eigen vectors:
    (eig_val, eig_vec) = eig(N)
    PC1 = eig_vec[:, 0]
    PC2 = eig_vec[:, 1]
    return PC1, PC2


def corrcoef_sparse(A, B=None):
    """
    Computes correlation coefficient on sparse matrices 

    Parameters
    ----------
    A : scipy.sparse.csr_matrix
        The matrix on which to compute the correlation.
    B: scipy.sparse.csr_matrix
        An optional second matrix. If provided, the correlation between A and B
        is computed.

    Returns
    -------
    scipy.sparse.csr_matrix
        The correlation matrix.
    """
    M = A.copy()
    if B is not None:
        M = sparse.vstack((A, B), format="csr")

    A = A.astype(np.float64)
    n = A.shape[1]
    # Compute the covariance matrix
    rowsum = A.sum(axis=1)
    centering = rowsum.dot(rowsum.T) / n
    C = (A.dot(A.T) - coo_matrix(centering)) / (n - 1)
    d = C.diagonal()
    coeffs = C / np.sqrt(np.outer(d, d))

    return coeffs


def compartments_sparse(M, normalize=True):
    """A/B compartment analysis

    Performs a detrending of the power law followed by a PCA-based A/B
    compartment analysis on a sparse, normalized, single chromosome contact map.
    The results are two vectors whose values (negative or positive) should
    presumably correlate with the presence of 'active' vs. 'inert' chromatin.

    Parameters
    ----------
    M : array_like
        The input, normalized contact map. Must be a single chromosome. Values
        are assumed to be only the upper triangle of a symmetrix matrix.
    normalize : bool
        Whether to normalize the matrix beforehand.
    mask : array of bool
        An optional boolean mask indicating which bins should be used

    Returns
    -------
    pr_comp : numpy.ndarray
        An array containing the N first principal component
    """
    if normalize:
        N = normalize_sparse(M, norm="SCN")
    else:
        N = copy.copy(M)
    N = N.tocoo()
    # Detrend by the distance law
    dist_bins, dist_vals = distance_law_from_mat(N, log_bins=False)
    N.data /= dist_vals[abs(N.row - N.col)]
    N = N.tocsr()
    # Make matrix symmetric (in case of upper triangle)
    if (abs(N - N.T) > 1e-10).nnz != 0:
        N = N + N.T
        N.setdiag(N.diagonal() / 2)
        N.eliminate_zeros()
    # Compute covariance matrix on full matrix
    N = N.tocsr()
    N = corrcoef_sparse(N)
    N[np.isnan(N)] = 0.0
    # Extract eigen vectors and eigen values
    [eigen_vals, pr_comp] = eig(N)

    return pr_comp[:, 0], pr_comp[:, 1]


def remove_intra(M, contigs, mask):
    """Remove intrachromosomal contacts

    Given a contact map and a list attributing each position
    to a given chromosome, set all contacts within each
    chromosome or contig to zero. Useful to perform
    calculations on interchromosomal contacts only.

    Parameters
    ----------
    M : array_like
        The initial contact map
    contigs : list or array_like
        A 1D array whose value at index i reflect the contig
        label of the row i in the matrix M. The length of
        the array must be equal to the (identical) shape
        value of the matrix.

    Returns
    -------
    N : numpy.ndarray
        The output contact map with no intrachromosomal contacts
    """

    N = np.copy(M)
    n = len(N)

    assert n == len(contigs)

    # Naive implmentation for now
    for (i, j) in itertools.product(range(n), range(n)):
        if contigs[i] == contigs[j]:
            N[i, j] = 0

    return N


def remove_inter(M, contigs):
    """Remove interchromosomal contacts

    Given a contact map and a list attributing each position
    to a given chromosome, set all contacts between each
    chromosome or contig to zero. Useful to perform
    calculations on intrachromosomal contacts only.

    Parameters
    ----------
    M : array_like
        The initial contact map
    contigs : list or array_like
        A 1D array whose value at index i reflect the contig
        label of the row i in the matrix M. The length of
        the array must be equal to the (identical) shape
        value of the matrix.

    Returns
    -------
    N : numpy.ndarray
        The output contact map with no interchromosomal contacts
    """

    N = np.copy(M)
    n = len(N)

    assert n == len(contigs)

    # Naive implmentation for now
    for (i, j) in itertools.product(range(n), range(n)):
        if contigs[i] != contigs[j]:
            N[i, j] = 0

    return N


def positions_to_contigs(positions):
    """Label contigs according to relative positions

    Given a list of positions, return an ordered list
    of labels reflecting where the positions array started
    over (and presumably a new contig began).

    Parameters
    ----------
    positions : list or array_like
        A piece-wise ordered list of integers representing
        positions

    Returns
    -------
    contig_labels : numpy.ndarray
        The list of contig labels

    """

    contig_labels = np.zeros_like(positions)

    contig_index = 0
    for i, p in enumerate(positions):
        if p == 0:
            contig_index += 1
        contig_labels[i] = contig_index

    return contig_labels


def contigs_to_positions(contigs, binning=10000):
    """Build positions from contig labels

    From a list of contig labels and a binning parameter,
    build a list of positions that's essentially a
    concatenation of linspaces with step equal to the
    binning.

    Parameters
    ----------
    contigs : list or array_like
        The list of contig labels, must be sorted.
    binning : int, optional
        The step for the list of positions. Default is 10000.

    Returns
    -------
    positions : numpy.ndarray
        The piece-wise sorted list of positions
    """

    positions = np.zeros_like(contigs)

    index = 0
    for _, chunk in itertools.groupby(contigs):
        l = len(chunk)
        positions[index : index + l] = np.arange(list(chunk)) * binning
        index += l

    return positions


def split_matrix(M, contigs):
    """Split multiple chromosome matrix

    Split a labeled matrix with multiple chromosomes
    into unlabeled single-chromosome matrices. Inter chromosomal
    contacts are discarded.

    Parameters
    ----------
    M : array_like
        The multiple chromosome matrix to be split
    contigs : list or array_like
        The list of contig labels
    """

    index = 0
    for _, chunk in itertools.groupby(contigs):
        l = len(chunk)
        yield M[index : index + l, index : index + l]
        index += l


def distance_law(
    size,
    prefactor=10000,
    gamma1=-0.5,
    gamma2=-1.5,
    inter=0.01,
    transition=None,
):
    """Generate a theoretical matrix from the usual P(s) model with two
    exponents for short-scale and large-scale modes and a sigmoid to represent
    the transition inbetween.

    Parameters
    ----------
    size : int 
        Size of the matrix (which will be of shape (size, size))
    prefactor : float
        Prefactor that's analogous to the coverage, by default 10000
    gamma1 : float, optional
        Exponent for the short-scale mode, by default -0.5
    gamma2 : float, optional
        Exponent for the large-scale mode, by default -1.5
    inter : float, optional
        Value for inter-chromosomal contacts that also represents the minimum
        value for intra-chromosomal contacts, by default 0.01
    transition : int, optional
        Coordinate of the transition between scale modes, by default 1/10 of
        the size

    Returns
    -------
    numpy.ndarray
        A symmetrical Toeplitz matrix whose each diagonal represents a value of
        the P(s) model
    """

    if transition is None:
        transition = size // 10

    def P(s, A=10000, gamma=-0.5, inter=0.01):
        return np.fmax(inter, A * (s ** gamma))

    P1 = functools.partial(P, gamma=gamma1)
    P2 = functools.partial(P, gamma=gamma2)

    def sigmoid(s, a=transition):  # ad hoc function to smoothen the transition
        return 1.0 / (1 + np.exp(-(s - transition) / a))

    s = np.arange(
        1, size + 1
    )  # Don't start from 0 to avoid weirdness and NaNs

    return toeplitz((1 - sigmoid(s)) * P1(s) + sigmoid(s) * P2(s))
