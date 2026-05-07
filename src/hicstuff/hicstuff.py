#!/usr/bin/env python

"""Common Hi-C functions

A bunch of handy functions for processing Hi-C data
(mainly in the form of matrices):

* Normalizations
* Interpolations
* Filters
* Removing artifacts
* Quick sum-pooling (aka 'binning') in sparse and dense form
* Simple models with parameter estimation
* Computing best-matching 3D structures
* Various metrics in use among Hi-C people for eyecandy purposes (directional index, domainograms, etc.)

These functions are meant to be simple and relatively quick
as-is implementations of procedures described in Hi-C papers.
"""

import copy
import itertools

import numpy as np
import pandas as pd
import scipy.sparse as sparse
from scipy.linalg import eig
from scipy.sparse import coo_matrix, csr_matrix, issparse, lil_matrix

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
    D = np.array([np.nanmean(matrix.diagonal(j)[included_bins[: n - j]]) for j in range(n)])
    if not log_bins:
        return np.array(range(len(D))), D
    else:
        n_bins = int(np.log(n) / np.log(base) + 1)
        logbin = np.unique(np.logspace(0, n_bins - 1, num=n_bins, base=base, dtype=int))
        logbin = np.insert(logbin, 0, 0)
        logbin[-1] = min(n, logbin[-1])
        if n < logbin.shape[0]:
            print("Not enough bins. Increase logarithm base.")
            return np.array(range(len(D))), D
        logD = np.array([np.nanmean(D[logbin[i - 1] : logbin[i]]) for i in range(1, len(logbin))])
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
    return coo_matrix((binned.dat, (binned.row, binned.col)), shape=(binned_n, binned_m))


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
    for coords, bin_frags in itertools.groupby(frag_idx, lambda x: tuple(frags[x, :])):
        bin_frags = list(bin_frags)
        first_frag, last_frag = bin_frags[0], bin_frags[-1] + 1
        # Pool row/col number by bin
        row[np.where((r.row >= first_frag) & (r.row < last_frag))] = actual_bin_No
        col[np.where((r.col >= first_frag) & (r.col < last_frag))] = actual_bin_No
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
    binned = coo_matrix((r.data, (row, col)), shape=(actual_bin_No, actual_bin_No))
    binned.sum_duplicates()
    binned.eliminate_zeros()

    return (binned, out_pos)


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
        dist = np.array(M.sum(axis=axis, dtype=float)).flatten()

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
    S_out, _ = trim_sparse(S, n_mad=n_mad, s_min=s_min, s_max=s_max)
    M_out = S_out.todense()
    return M_out


def trim_sparse(M, n_mad=3, s_min=None, s_max=None, chrom_start=None):
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
    lines : bool
        Either to return the offset of the chromosomes for lines plotting.

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
    chrom_start_offset = None
    if chrom_start is not None:
        chrom_start_offset = [miss_map[start] for start in chrom_start]
    # Indices of cells that will be kept
    indices = np.where(f[r.row] & f[r.col])
    # Remove sparse rows and shift indices accordingly
    rows = [miss_map[i] for i in r.row[indices]]
    cols = [miss_map[j] for j in r.col[indices]]
    data = r.data[indices]
    size = max(max(rows, default=-1), max(cols, default=-1)) + 1
    N = coo_matrix((data, (rows, cols)), shape=(size, size))
    return N, chrom_start_offset


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
            print("Please install it from https://bitbucket.org/mirnylab/mirnylib")
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
        raise Exception('Unknown norm, please specify one of ("mirnylib", "SCN", "frag")')

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
        for _i in range(iterations):
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
    sampled_contacts = np.random.choice(int(tot_contacts), size=int(n_contacts), replace=False)

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
    if type(M) is not np.ndarray:
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
                N[i, j] = M[i, i - j :].sum() + M[i, :].sum() + M[i, : i + j - n + 1].sum()
    return N


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
    eig_val, eig_vec = eig(N)
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
    A.copy()
    if B is not None:
        sparse.vstack((A, B), format="csr")

    A = A.astype(np.float64)
    n = A.shape[1]
    # Compute the covariance matrix
    rowsum = A.sum(axis=1)
    centering = rowsum.dot(rowsum.T) / n
    C = (A.dot(A.T) - coo_matrix(centering)) / (n - 1)
    d = np.asarray(C.diagonal()).flatten()
    coeffs = np.asarray(C.todense()) / np.sqrt(np.outer(d, d))

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
