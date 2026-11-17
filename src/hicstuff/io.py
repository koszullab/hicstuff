#!/usr/bin/env python3

import bz2
import collections
import contextlib
import functools
import gzip
import heapq
import io
import itertools
import os
import re
import subprocess as sp
import sys
import tempfile
import zipfile
from decimal import Decimal
from os.path import exists, join
from random import getrandbits
from shutil import which

import cooler
import numpy as np
import pandas as pd
import pysam as ps
import scipy.stats as ss
from Bio import SeqIO, SeqUtils
from scipy.sparse import coo_matrix, tril, triu

import hicstuff.hicstuff as hcs
from hicstuff import __version__
from hicstuff.log import logger

DEFAULT_MAX_MATRIX_SHAPE = 10000

DEFAULT_FRAGMENTS_LIST_FILE_NAME = "fragments_list.txt"
DEFAULT_INFO_CONTIGS_FILE_NAME = "info_contigs.txt"
DEFAULT_SPARSE_MATRIX_FILE_NAME = "abs_fragments_contacts_weighted.txt"


def _cols_to_sparse(sparse_array, shape=None, dtype=np.float64):
    """
    Make a coordinate based sparse matrix from columns.
    Convert (3, n) shaped arrays to a sparse matrix. The first
    one acts as the row coordinates, the second one as the
    column coordinates, the third one as the data points
    for each pair. If duplicates are found, the data points are
    added.

    Parameters
    ----------
    sparse_array : array_like
        An array with exactly three columns representing the sparse
        matrix data in coordinate format.
    shape : tuple of int
        The total number of rows and columns in the matrix. Will be estimated
        from nonzero values if omitted.
    dtype : type, optional
        The type of data being loaded. Default is numpy.float64

    Example
    -------

        >>> import numpy as np
        >>> row, col = np.array([1, 2, 3]), np.array([3, 2, 1])
        >>> data = np.array([4, 5, 6])
        >>> M = np.array([row, col, data]).T
        >>> S = _cols_to_sparse(M)
        >>> print(S.todense())
        [[0. 0. 0. 0.]
         [0. 0. 0. 4.]
         [0. 0. 5. 0.]
         [0. 6. 0. 0.]]
    """

    row = sparse_array[:, 0]
    col = sparse_array[:, 1]
    data = sparse_array[:, 2]
    if shape is None:
        n = int(max(np.amax(row), np.amax(col))) + 1
        shape = (n, n)

    S = coo_matrix((data, (row, col)), shape=shape, dtype=dtype)
    return S


def load_sparse_matrix(mat_path, binning=1, dtype=np.float64):
    """Load sparse matrix

    Load a text file matrix into a sparse matrix object. The expected format is
    a 3 column file where columns are row_number, col_number, value. The first
    line consists of 3 values representing the total number of rows, columns
    and nonzero values.

    Parameters
    ----------
    mat_path : file, str or pathlib.Path
        The input matrix file in instagraal format.
    binning : int or "auto"
        The binning to perform. If "auto", binning will
        be automatically inferred so that the matrix size
        will not go beyond (10000, 10000) in shape. That
        can be changed by modifying the DEFAULT_MAX_MATRIX_SHAPE
        value. Default is 1, i.e. no binning is performed
    dtype : type, optional
        The type of data being loaded. Default is numpy.float64

    Returns
    -------
    sparse_mat : scipy.sparse.coo_matrix
        The output (sparse) matrix in COOrdinate format.

    Examples
    --------
    >>> S = load_sparse_matrix('test_data/mat_5kb.tsv', binning=1)
    >>> S.data[:10]
    array([84.,  2.,  3.,  2.,  1.,  1., 50.,  1., 66.,  1.])
    >>> S.shape
    (16, 16)
    """
    try:
        # Numeric header: the number of rows, columns and nonzero entries
        cols_arr = np.loadtxt(mat_path, delimiter="\t", dtype=dtype, ndmin=2)
        shape = (int(cols_arr[0, 0]), int(cols_arr[0, 1]))
        cols_arr = cols_arr[1:, :]
    except ValueError:
        # Text header (column names): infer the shape from the entries
        cols_arr = np.loadtxt(mat_path, delimiter="\t", dtype=dtype, skiprows=1, ndmin=2)
        shape = None

    sparse_mat = _cols_to_sparse(cols_arr, shape=shape, dtype=dtype)

    if binning == "auto":
        # Smallest factor that keeps the matrix within the maximum shape
        num_bins = max(sparse_mat.shape)
        subsampling_factor = max(1, -(-num_bins // DEFAULT_MAX_MATRIX_SHAPE))
    else:
        subsampling_factor = binning
    sparse_mat = hcs.bin_sparse(sparse_mat, subsampling_factor=subsampling_factor)
    return sparse_mat


def save_sparse_matrix(s_mat, path):
    """Save a sparse matrix

    Saves a sparse matrix object into tsv format.

    Parameters
    ----------
    s_mat : scipy.sparse.coo_matrix
        The sparse matrix to save on disk
    path : str
        File path where the matrix will be stored
    """
    if s_mat.format != "coo":
        raise ValueError("Sparse matrix must be in coo format")
    dtype = s_mat.dtype
    fmt = "%i" if isinstance(dtype, int) else "%.10e"
    sparse_arr = np.vstack([s_mat.row, s_mat.col, s_mat.data]).T

    np.savetxt(
        path,
        sparse_arr,
        header=f"{s_mat.shape[0]}\t{s_mat.shape[1]}\t{s_mat.nnz}",
        comments="",
        fmt=fmt,
        delimiter="\t",
    )


def load_pos_col(path, colnum, header=1, dtype=np.int64):
    """
    Loads a single column of a TSV file with header into a numpy array.

    Parameters
    ----------
    path : str
        The path of the TSV file to load.
    colnum : int
        The 0-based index of the column to load.
    header : int
        Number of line to skip. By default the header is a single line.

    Returns
    -------
    numpy.array :
        A 1D numpy array with the

    Examples
    --------
    >>> load_pos_col('test_data/mat_5kb.tsv', 0)[:10]
    array([0, 0, 0, 0, 0, 0, 1, 1, 2, 2])
    """
    pos_arr = np.genfromtxt(
        path,
        delimiter="\t",
        usecols=(colnum,),
        skip_header=header,
        dtype=dtype,
    )
    return pos_arr


def generate_temp_dir(path):
    """Temporary directory generation

    Generates a temporary file with a random name at the input path.

    Parameters
    ----------
    path : str
        The path at which the temporary directory will be created.

    Returns
    -------
    str
        The path of the newly created temporary directory.
    """
    exist = True
    while exist:
        # Keep trying random directory names if they already exist
        directory = str(hex(getrandbits(32)))[2:]
        full_path = join(path, directory)
        if not exists(full_path):
            exist = False
    try:
        os.makedirs(full_path)
    except PermissionError:
        raise PermissionError(
            f"The temporary directory cannot be created in {path}. "
            "Make sure you have write permission."
        ) from None
    return full_path


def _check_cooler(fun):
    """Decorates function `fun` to check if cooler is available.."""

    @functools.wraps(fun)
    def wrapped(*args, **kwargs):
        try:
            import cooler

            fun.__globals__["cooler"] = cooler
        except ImportError:
            logger.error(
                f"The cooler package is required to use {fun.__name__}, please install it first"
            )
            raise ImportError("The cooler package is required.") from None
        return fun(*args, **kwargs)

    return wrapped


@_check_cooler
def add_cool_column(clr, column, column_name, table_name="bins", metadata=None, dtype=None):
    """
    Adds a new column to a loaded Cooler store. If the column exists,
    it is replaced. This will affect the .cool file.

    Parameters
    ----------
    clr : Cooler object
        A Cooler store.
    column : pandas Series
        The column to add to the cooler.
    column_name : str
        The name of the column to add.
    table_name : str
        The name of the table to which the column should be added.
        Defaults to the "bins" table.
    metadata : dict
        A dictionary of metadata to associate with the new column.
    """
    if metadata is None:
        metadata = {}
    with clr.open("r+") as c:
        if column_name in c[table_name]:
            del c[table_name][column_name]
        h5opts = dict(compression="gzip", compression_opts=6)
        c[table_name].create_dataset(column_name, data=column, dtype=dtype, **h5opts)
        c[table_name][column_name].attrs.update(metadata)


@_check_cooler
def load_cool(cool):
    """
    Reads a cool file into memory and parses it into graal style tables.

    Parameters
    ----------
    cool : str
        Path to the input .cool file.

    Returns
    -------
    mat : scipy coo_matrix
        Hi-C contact map in COO format.
    frags : pandas DataFrame
        Table off bins matching the matrix. Corresponds to the content of the fragments_list.txt file.
    chroms : pandas DataFrame
        Table of chromosome informations.
    """
    c = cooler.Cooler(cool)  # pylint: disable=undefined-variable
    frags = c.bins()[:]
    chroms = c.chroms()[:]
    mat = c.pixels()[:]
    frags.rename(columns={"start": "start_pos", "end": "end_pos"}, inplace=True)
    frags["id"] = frags.groupby("chrom", sort=False).cumcount() + 1
    # Try loading hicstuff-specific columns
    try:
        frags = frags[["id", "chrom", "start_pos", "end_pos", "size", "gc_content"]]
    # If absent, only load standard columns
    except KeyError:
        frags = frags[["id", "chrom", "start_pos", "end_pos"]]

    chroms["cumul_length"] = chroms.length.shift(1).fillna(0).cumsum().astype(int)
    n_frags = c.bins()[:].groupby("chrom", sort=False).count().start
    chroms["n_frags"] = chroms.merge(n_frags, right_index=True, left_on="name", how="left").start
    chroms.rename(columns={"name": "contig"}, inplace=True)
    # One row per bin, including empty bins at the end of the genome
    mat = coo_matrix((mat["count"], (mat.bin1_id, mat.bin2_id)), shape=c.shape)

    return mat, frags, chroms


@_check_cooler
def save_cool(cool_out, mat, frags, metadata=None):
    """
    Writes a .cool file from graal style tables.

    Parameters
    ----------
    cool_out : str
        Path to the output cool file.
    mat : scipy coo_matrix
        The Hi-C contact matrix in sparse COO format.
    frags : pandas DataFrame
        The graal style 'fragments_list' table.
    metadata : dict
        Potential metadata to associate with the cool file.
    """
    if metadata is None:
        metadata = {}
    up_tri = False
    # Check if symmetric matrix is symmetric
    # (i.e. only upper triangle or full mat)
    if (abs(mat - mat.T) > 1e-10).nnz != 0:
        up_tri = True
    # Drop useless column
    try:
        bins = frags.drop("id", axis=1)
    except KeyError:
        bins = frags
    # Get column names right
    bins.rename(
        columns={"seq": "chrom", "start_pos": "start", "end_pos": "end"},
        inplace=True,
    )
    mat_dict = {"bin1_id": mat.row, "bin2_id": mat.col, "count": mat.data}
    pixels = pd.DataFrame(mat_dict)
    cooler.create_cooler(  # pylint: disable=undefined-variable
        cool_out,
        bins,
        pixels,
        metadata=metadata,
        symmetric_upper=up_tri,
        triucheck=False,
    )


def read_compressed(filename, mode="r"):
    """Read compressed file

    Opens the file in read mode with appropriate decompression algorithm.

    Parameters
    ----------
    filename : str
        The path to the input file
    mode : str
        "r" to read text, "rb" to read bytes.

    Returns
    -------
    file-like object
        The handle to access the input file's content

    """

    # Standard header bytes for diff compression formats
    comp_bytes = {
        b"\x1f\x8b\x08": "gz",
        b"\x42\x5a\x68": "bz2",
        b"\x50\x4b\x03\x04": "zip",
    }

    max_len = max(len(x) for x in comp_bytes)

    def file_type(filename):
        """Guess file type

        Compare header bytes with those in the file and return type.
        """
        with open(filename, "rb") as f:
            file_start = f.read(max_len)
        for magic, filetype in comp_bytes.items():
            if file_start.startswith(magic):
                return filetype
        return "uncompressed"

    # Open file with appropriate function
    mode_map = {"r": "rt", "rb": "rb"}
    comp = file_type(filename)
    if comp == "gz":
        return gzip.open(filename, mode_map[mode])
    elif comp == "bz2":
        return bz2.open(filename, mode_map[mode])
    elif comp == "zip":
        zip_arch = zipfile.ZipFile(filename, "r")
        if len(zip_arch.namelist()) > 1:
            raise OSError("Only a single fastq file must be in the zip archive.")
        else:
            zip_content = zip_arch.open(zip_arch.namelist()[0], "r")
            if mode == "rb":
                return zip_content
            # ZipFile opens as bytes by default, using io to read as text
            return io.TextIOWrapper(zip_content, encoding="utf-8")
    else:
        return open(filename, mode)


def is_compressed(filename):
    """Check compression status

    Check if the input file is compressed from the first bytes.

    Parameters
    ----------
    filename : str
        The path to the input file

    Returns
    -------
    bool
        True if the file is compressed, False otherwise.
    """

    # Standard header bytes for diff compression formats
    comp_bytes = {
        b"\x1f\x8b\x08": "gz",
        b"\x42\x5a\x68": "bz2",
        b"\x50\x4b\x03\x04": "zip",
    }
    max_len = max(len(x) for x in comp_bytes)
    with open(filename, "rb") as f:
        file_start = f.read(max_len)
    for magic, _ in comp_bytes.items():
        if file_start.startswith(magic):
            return True
    return False


def from_dade_matrix(filename, header=False):
    """Load a DADE matrix

    Load a numpy array from a DADE matrix file, optionally
    returning bin information from the header. Header data
    processing is delegated downstream.

    Parameters
    ----------
    filename : str, file or pathlib.Path
        The name of the file containing the DADE matrix.
    header : bool
        Whether to return as well information contained
        in the header. Default is False.

    Example
    -------
        >>> import numpy as np
        >>> import tempfile
        >>> lines = [['RST', 'chr1~0', 'chr1~10', 'chr2~0', 'chr2~30'],
        ...          ['chr1~0', '5'],
        ...          ['chr1~10', '8', '3'],
        ...          ['chr2~0', '3', '5', '5'],
        ...          ['chr2~30', '5', '10', '11', '2']
        ...          ]
        >>> formatted = ["\\t".join(l) + "\\n" for l in lines ]
        >>> dade = tempfile.NamedTemporaryFile(mode='w')
        >>> for fm in formatted:
        ...     dade.write(fm)
        34
        9
        12
        13
        18
        >>> dade.flush()
        >>> M, h = from_dade_matrix(dade.name, header=True)
        >>> dade.close()
        >>> print(M)
        [[ 5.  8.  3.  5.]
         [ 8.  3.  5. 10.]
         [ 3.  5.  5. 11.]
         [ 5. 10. 11.  2.]]

        >>> print(h)
        ['chr1~0', 'chr1~10', 'chr2~0', 'chr2~30']

    See https://github.com/scovit/DADE for more details about Dade.
    """

    A = pd.read_csv(filename, sep="\t", header=None)
    A.fillna("0", inplace=True)
    M, headers = np.array(A.iloc[1:, 1:], dtype=np.float64), A.iloc[0, :]
    matrix = M + M.T - np.diag(np.diag(M))
    if header:
        return matrix, headers.tolist()[1:]
    else:
        return matrix


def to_dade_matrix(M, annotations=None, filename=None):
    """Returns a Dade matrix from input numpy matrix. Any annotations are added
    as header. If filename is provided and valid, said matrix is also saved
    as text.

    Parameters
    ----------
    M : numpy.ndarray
        The square matrix to convert.
    annotations : array_like of numbers, optional
        Labels of the header row and column, one more than the number of rows
        (the first one goes in the corner). The header is left at 0 if omitted.
    filename : str, optional
        Where to save the matrix, as integers.

    Returns
    -------
    numpy.ndarray :
        The matrix with its header row and column.
    """

    n, m = M.shape
    A = np.zeros((n + 1, m + 1))
    A[1:, 1:] = M
    if annotations is not None and len(annotations) > 0:
        A[0, :] = annotations
        A[:, 0] = annotations
    if filename:
        try:
            np.savetxt(filename, A, fmt="%i")
            logger.info(f"I saved input matrix in dade format as {str(filename)}")
        except ValueError as e:
            logger.warning("I couldn't save input matrix.")
            logger.warning(str(e))

    return A


def load_into_redis(filename):
    """Load a file into redis

    Load a matrix file and sotre it in memory with redis.
    Useful to pass around huge datasets from scripts to
    scripts and load them only once.

    Inspired from https://gist.github.com/alexland/ce02d6ae5c8b63413843

    Parameters
    ----------
    filename : str, file or pathlib.Path
        The file of the matrix to load.

    Returns
    -------
    key : str
        The key of the dataset needed to retrieve it from redis.
    """
    try:
        import time

        from redis import StrictRedis as redis
    except ImportError:
        print(
            "Error! Redis does not appear to be installed in your system.",
            file=sys.stderr,
        )
        exit(1)

    M = np.genfromtxt(filename, dtype=None)
    array_dtype = str(M.dtype)
    m, n = M.shape
    M = M.ravel().tostring()
    database = redis(host="localhost", port=6379, db=0)
    key = f"{int(time.time())}|{array_dtype}#{m}#{n}"

    database.set(key, M)

    return key


def load_from_redis(key):
    """Retrieve a dataset from redis

    Retrieve a cached dataset that was stored in redis
    with the input key.

    Parameters
    ----------
    key : str
        The key of the dataset that was stored in redis.

    Returns
    -------
    M : numpy.ndarray
        The retrieved dataset in array format.
    """

    try:
        from redis import StrictRedis as redis
    except ImportError:
        print(
            "Error! Redis does not appear to be installed in your system.",
            file=sys.stderr,
        )
        exit(1)

    database = redis(host="localhost", port=6379, db=0)

    try:
        M = database.get(key)
    except KeyError:
        print(
            "Error! No dataset was found with the supplied key.",
            file=sys.stderr,
        )
        exit(1)

    array_dtype, n, m = key.split("|")[1].split("#")

    M = np.fromstring(M, dtype=array_dtype).reshape(int(n), int(m))
    return M


def dade_to_graal(
    filename,
    output_matrix=DEFAULT_SPARSE_MATRIX_FILE_NAME,
    output_contigs=DEFAULT_INFO_CONTIGS_FILE_NAME,
    output_frags=DEFAULT_FRAGMENTS_LIST_FILE_NAME,
    output_dir=None,
):
    """Convert a matrix from DADE format (https://github.com/scovit/dade)
    to a graal-compatible format. Since DADE matrices contain both fragment
    and contact information all files are generated at the same time.
    """

    with open(output_matrix, "w") as sparse_file:
        sparse_file.write("id_frag_a\tid_frag_b\tn_contact\n")
        with open(filename) as file_handle:
            first_line = file_handle.readline()
            for row_index, line in enumerate(file_handle):
                dense_row = np.array(line.split("\t")[1:], dtype=np.int32)
                for col_index in np.nonzero(dense_row)[0]:
                    line_to_write = f"{row_index}\t{col_index}\t{dense_row[col_index]}\n"
                    sparse_file.write(line_to_write)

    header = first_line.split("\t")
    bin_type = header[0]
    if bin_type == '"RST"':
        logger.info("I detected fragment-wise binning")
    elif bin_type == '"BIN"':
        logger.info("I detected fixed size binning")
    else:
        logger.warning(f"Sorry, I don't understand this matrix's binning: I read {str(bin_type)}")

    header_data = [
        header_elt.replace("'", "").replace('"', "").replace("\n", "").split("~")
        for header_elt in header[1:]
    ]

    (
        global_frag_ids,
        contig_names,
        local_frag_ids,
        frag_starts,
        frag_ends,
    ) = np.array(list(zip(*header_data)))

    frag_starts = frag_starts.astype(np.int32) - 1
    frag_ends = frag_ends.astype(np.int32) - 1
    frag_lengths = frag_ends - frag_starts

    total_length = len(global_frag_ids)

    with open(output_contigs, "w") as info_contigs:
        info_contigs.write("contig\tlength\tn_frags\tcumul_length\n")

        cumul_length = 0

        for contig in collections.OrderedDict.fromkeys(contig_names):
            length_tig = np.sum(frag_lengths[contig_names == contig])
            n_frags = collections.Counter(contig_names)[contig]
            line_to_write = f"{contig}\t{length_tig}\t{n_frags}\t{cumul_length}\n"
            info_contigs.write(line_to_write)
            cumul_length += n_frags

    with open(output_frags, "w") as fragments_list:
        fragments_list.write("id\tchrom\tstart_pos\tend_pos\tsize\tgc_content\n")
        bogus_gc = 0.5

        for i in range(total_length):
            line_to_write = f"{int(local_frag_ids[i]) + 1}\t{contig_names[i]}\t{frag_starts[i]}\t{frag_ends[i]}\t{frag_lengths[i]}\t{bogus_gc}\n"
            fragments_list.write(line_to_write)


def load_bedgraph2d(filename, bin_size=None, fragments_file=None):
    """
    Loads matrix and fragment information from a 2D bedgraph file. Note this
    function assumes chromosomes are ordered in alphabetical. order

    Parameters
    ----------
    filename : str
        Path to the bedgraph2D file.
    bin_size : int
        The size of bins in the case of fixed bin size.
    fragments_file : str
        Path to a fragments file to explicitely provide fragments positions.
        If the matrix does not have fixed bin size, this prevents errors.

    Returns
    -------
    mat : scipy.sparse.coo_matrix
        The Hi-C contact map as the upper triangle of a symetric matrix, in
        sparse format.
    frags : pandas.DataFrame
        The list of fragments/bin present in the matrix with their genomic
        positions.
    """
    bed2d = pd.read_csv(filename, sep="\t", header=None)
    chrom_sizes = {}
    if bin_size is not None:
        # If bin size if provided, retrieve chromosome lengths, this will be
        # used when regenerating bin coordinates
        chroms_left = bed2d[[3, 5]]
        chroms_left.columns = [0, 2]
        chroms = pd.concat([bed2d[[0, 2]], chroms_left]).groupby([0], sort=False).max()
        for chrom, size in zip(chroms.index, np.array(chroms)):
            chrom_sizes[chrom] = size[0]
    elif fragments_file is None:
        logger.warning(
            "Please be aware that not all information can be restored from a "
            "bg2 file without fixed bin size; fragments without any contact "
            "will be lost"
        )
    # Get all possible fragment chrom-positions into an array
    frag_pos = np.vstack([np.array(bed2d[[0, 1, 2]]), np.array(bed2d[[3, 4, 5]])])
    # Sort by position (least important, col 1)
    frag_pos = frag_pos[frag_pos[:, 1].argsort(kind="mergesort")]
    # Then by chrom (most important, col 0)
    frag_pos = frag_pos[frag_pos[:, 0].argsort(kind="mergesort")]
    # Get unique names for fragments (chrom+pos)
    ordered_frag_pos = pd.DataFrame(frag_pos).drop_duplicates().reset_index(drop=True)
    frag_pos_a = bed2d[[0, 1]].apply(lambda x: tuple(x), axis=1)
    frag_pos_b = bed2d[[3, 4]].apply(lambda x: tuple(x), axis=1)
    # If fragments file is provided, use fragments positions to indices mapping
    if fragments_file is not None:
        frags = pd.read_csv(fragments_file, delimiter="\t")
        frag_map = frags.apply(lambda x: (str(x.chrom), x.start_pos), axis=1)
        frag_map = {f_name: f_idx for f_idx, f_name in enumerate(frag_map)}
    # If fixed fragment size available, use it to reconstruct original
    # fragments ID (even if they are absent from the bedgraph file).
    elif bin_size is not None:
        frag_map = {}
        chrom_frags = []
        for chrom, size in chrom_sizes.items():
            prev_frags = len(frag_map)
            for bin_id, bin_pos in enumerate(range(0, size, bin_size)):
                frag_map[(chrom, bin_pos)] = bin_id + prev_frags
            # The last bin is shorter unless size is a multiple of bin_size
            starts = np.arange(0, size, bin_size)
            chrom_frags.append(
                pd.DataFrame(
                    {
                        "id": range(1, len(starts) + 1),
                        "chrom": np.repeat(chrom, len(starts)),
                        "start_pos": starts,
                        "end_pos": np.minimum(starts + bin_size, size),
                    }
                )
            )
        frags = pd.concat(chrom_frags, axis=0).reset_index(drop=True)
        frags.insert(loc=3, column="size", value=frags.end_pos - frags.start_pos)
    # If None available, guess fragments indices from bedgraph (potentially wrong)
    else:
        frag_map = {(v[0], v[1]): i for i, v in ordered_frag_pos.iloc[:, [0, 1]].iterrows()}
        frags = ordered_frag_pos.copy()
        frags[3] = frags.iloc[:, 2] - frags.iloc[:, 1]
        frags.insert(loc=0, column="id", value=0)
        frags.id = frags.groupby([0], sort=False).cumcount() + 1
        frags.columns = ["id", "chrom", "start_pos", "end_pos", "size"]
    # Match bin indices to their names
    frag_id_a = np.array(list(map(lambda x: frag_map[x], frag_pos_a)))
    frag_id_b = np.array(list(map(lambda x: frag_map[x], frag_pos_b)))
    contacts = np.array(bed2d.iloc[:, 6].tolist())
    # Use index to build matrix
    n_frags = len(frag_map.keys())
    mat = coo_matrix((contacts, (frag_id_a, frag_id_b)), shape=(n_frags, n_frags))

    # Get size of each chromosome in basepairs
    chromsizes = frags.groupby("chrom", sort=False).apply(lambda x: np.int64(max(x.end_pos)))
    chrom_bins = frags.groupby("chrom", sort=False).size()
    # Shift chromsizes by one to get starting bin, first one is zero
    # Make chromsize cumulative to get start bin of each chrom
    # Get chroms into a 1D array of bin starts
    chrom_start = chrom_bins.shift(1, fill_value=0).cumsum()
    chroms = pd.DataFrame(
        {
            "contig": chromsizes.index,
            "length": chromsizes.values,
            "n_frags": chrom_bins,
            "cumul_length": chrom_start,
        }
    )
    return mat, frags, chroms


def flexible_hic_loader(mat, fragments_file=None, chroms_file=None, quiet=False):
    """
    Wrapper function to load COO, bg2 or cool input and return the same output.
    COO formats requires fragments_file and chroms_file options. bg2 format can
    infer bin_size if fixed. When providing a bg2 matrix with uneven fragments
    length, one should provide fragments_file as well or empty bins will be
    truncated from the output.

    Parameters
    ----------
    mat : str
        Path to the matrix in graal, bedgraph2 or cool format.
    fragments_file : str or None
        Path to the file with fragments information (fragments_list.txt).
        Only required if the matrix is in graal format.
    chroms_file : str or None
        Path to the file with chromosome information (info_contigs.txt). Only required
        if the matrix is in graal format.
    quiet : bool
        If True, will silence warnings for empty outputs.

    Returns
    -------
    mat : scipy.sparse.coo_matrix
        Sparse upper triangle Hi-C matrix.
    frags : pandas.DataFrame or None
        Table of fragment informations. None if information was not provided.
    chroms : pandas.DataFrame or None
        Table of chromosomes/contig information. None if information was not provided.
    """
    hic_format = get_hic_format(mat)
    # Load cool based on file extension
    if hic_format == "cool":
        mat, frags, chroms = load_cool(mat)
    # Use the first line to determine COO / bg2 format
    elif hic_format == "bg2":
        # Use the frags file to define bins if available
        if fragments_file is not None:
            mat, frags, chroms = load_bedgraph2d(mat, fragments_file=fragments_file)
        else:
            # Guess if bin size is fixed based on MAD
            bg2 = pd.read_csv(mat, sep="\t")
            sizes = np.array(bg2.iloc[:, 2] - bg2.iloc[:, 1])
            size_mad = ss.median_abs_deviation(sizes, scale="normal")
            # Use only the bg2
            if size_mad > 0:
                mat, frags, chroms = load_bedgraph2d(mat)
                logger.warning(
                    "Input is a bedgraph2d file with uneven bin size, "
                    "but no fragments_file was provided. Empty bins will "
                    "be missing from the output. To avoid this, provide a "
                    "fragments file."
                )
            # Use fixed bin size
            else:
                mat, frags, chroms = load_bedgraph2d(mat, bin_size=int(np.median(sizes)))

    elif hic_format == "graal":
        mat = load_sparse_matrix(mat)
        try:
            frags = pd.read_csv(fragments_file, sep="\t")
        except ValueError:
            if not quiet:
                logger.warning(
                    "fragments_file was not provided when "
                    "loading a matrix in COO/graal format. frags will be None."
                )
            frags = None
        try:
            chroms = pd.read_csv(chroms_file, sep="\t")
        except ValueError:
            if not quiet:
                logger.warning(
                    "chroms_file was not provided when "
                    "loading a matrix in COO/graal format. chroms will be None."
                )

            chroms = None

    else:
        raise ValueError(f"Unknown input format: {hic_format}")

    # Ensure the matrix is upper triangle symmetric
    if mat.shape[0] == mat.shape[1]:
        if (abs(mat - mat.T) > 1e-10).nnz > 0:
            mat = mat + tril(mat, k=-1).T
        mat = triu(mat, format="coo")

    return mat, frags, chroms


def get_hic_format(mat):
    """Returns the format of the input Hi-C matrix

    Parameters
    ----------
    mat : str
        Path to the input Hi-C matrix
    Returns
    -------
    str :
        Hi-C format string. One of graal, bg2, cool
    """
    if mat.endswith(".cool") or mat.count(".mcool::/") == 1 or mat.count(".cool::/") == 1:
        hic_format = "cool"
    else:
        # Use the first line to determine COO / bg2 format
        ncols = len(open(mat).readline().split("\t"))
        if ncols == 7:
            hic_format = "bg2"
        elif ncols == 3:
            hic_format = "graal"
        else:
            raise ValueError("Unkown file format")
    return hic_format


def flexible_hic_saver(
    mat,
    out_prefix,
    frags=None,
    chroms=None,
    hic_fmt="graal",
    quiet=False,
):
    """
    Saves objects to the desired Hi-C file format.

    Parameters
    ----------
    mat : scipy.sparse.coo_matrix
        Hi-C contact map.
    out_prefix : str
        Output path without extension (the extension is added based on hic_fmt).
    frags : pandas.DataFrame or None
        Table of fragments informations.
    chroms : pandas.DataFrame or None
        Table of chromosomes / contigs informations.
    hic_fmt : str
        Output format. Can be one of graal for graal-compatible COO format, bg2 for
        2D bedgraph format, or cool for cooler compatible format.
    """
    if hic_fmt == "graal":
        save_sparse_matrix(mat, out_prefix + ".mat.tsv")
        try:
            frags.to_csv(out_prefix + ".frags.tsv", sep="\t", index=False)
        except AttributeError:
            if not quiet:
                logger.warning("Could not create fragments_list.txt from input files")
        try:
            chroms.to_csv(out_prefix + ".chr.tsv", sep="\t", index=False)
        except AttributeError:
            if not quiet:
                logger.warning("Could not create info_contigs.txt from input files")
    elif hic_fmt == "cool":
        frag_sizes = frags.end_pos - frags.start_pos
        size_mad = np.median(frag_sizes - np.median(frag_sizes))
        bin_type = "variable" if size_mad else "fixed"
        try:
            save_cool(
                out_prefix + ".cool",
                mat,
                frags,
                metadata={"hicstuff": __version__, "bin-type": bin_type},
            )
        except NameError:
            NameError("frags is required to save a cool file")
    elif hic_fmt == "bg2":
        try:
            save_bedgraph2d(mat, frags, out_prefix + ".bg2")
        except NameError:
            NameError("frags is required to save a bg2 file")
    else:
        raise ValueError(f"Unknown output format: {hic_fmt}")


def save_bedgraph2d(mat, frags, out_path):
    """
    Given a sparse matrix and a corresponding list of fragments, save a file
    in 2D bedgraph format.

    Parameters
    ----------
    mat : scipy.sparse.coo_matrix
        The sparse contact map.
    frags : pandas.DataFrame
        A structure containing the annotations for each matrix bin. Should
        correspond to the content of the fragments_list.txt file.

    """

    mat_df = pd.DataFrame({"row": mat.row, "col": mat.col, "data": mat.data.astype(int)})
    # Merge fragments with matrix based on row indices to annotate rows
    merge_mat = mat_df.merge(frags, left_on="row", right_index=True, how="left", suffixes=("", ""))
    # Rename annotations to assign to frag1
    merge_mat.rename(
        columns={"chrom": "chr1", "start_pos": "start1", "end_pos": "end1"},
        inplace=True,
    )
    # Do the same operation for cols (frag2)
    merge_mat = merge_mat.merge(
        frags,
        left_on="col",
        right_index=True,
        how="left",
        suffixes=("_0", "_2"),
    )
    merge_mat.rename(
        columns={"chrom": "chr2", "start_pos": "start2", "end_pos": "end2"},
        inplace=True,
    )
    # Select only relevant columns in correct order
    bg2 = merge_mat.loc[:, ["chr1", "start1", "end1", "chr2", "start2", "end2", "data"]]
    bg2.to_csv(out_path, header=None, index=False, sep="\t")


def get_pos_cols(df):
    """Get column names representing chromosome, start and end column
    from a dataframe. Allows flexible names.

    Parameters
    ----------
    df : pd.DataFrame

    Returns
    -------
    tuple of str:
        Tuple containing the names of the chromosome, start and end
        columns in the input dataframe.

    Examples
    --------
    >>> import pandas as pd
    >>> d = [1, 2, 3]
    >>> df = pd.DataFrame(
    ...     {'Chromosome': d, 'Start': d, 'End': d, 'species': d}
    ... )
    >>> get_pos_cols(df)
    ('Chromosome', 'Start', 'End')
    >>> df = pd.DataFrame(
    ...     {'id': d, 'chr': d, 'start_bp': d, 'end_bp': d}
    ... )
    >>> get_pos_cols(df)
    ('chr', 'start_bp', 'end_bp')
    """

    # Get case insensitive column names
    cnames = df.columns
    inames = cnames.str.lower()

    def _col_getter(cols, pat):
        """Helper to get column index from the start of its name"""
        mask = cols.str.startswith(pat)
        if not np.any(mask):
            raise ValueError(f"Could not find a column starting with {pat}")
        idx = np.flatnonzero(mask)[0]
        return idx

    chrom_col = cnames[_col_getter(inames, "chr")]
    start_col = cnames[_col_getter(inames, "start")]
    end_col = cnames[_col_getter(inames, "end")]
    return chrom_col, start_col, end_col


def gc_bins(genome_path, frags):
    """Generate GC content annotation for bins using input genome.

    Parameters
    ----------
    genome_path : str
        Path the the genome file in FASTA format.
    frags : pandas.DataFrame
        Table containing bin segmentation of the genome.
        Required columns: chrom, start, end.

    Returns
    -------
    numpy.ndarray of floats:
        GC content per bin, in the range [0.0, 1.0].
    """
    # Grab columns by name (catch start, Start, star_pos, etc)
    chrom_col, start_col, end_col = get_pos_cols(frags)
    # Fill the gc array by chromosome
    gc_bins = np.zeros(frags.shape[0], dtype=float)
    for rec in SeqIO.parse(genome_path, "fasta"):
        mask = frags[chrom_col] == rec.id
        # Define coordinates of each bin
        starts = frags.loc[mask, start_col]
        ends = frags.loc[mask, end_col]
        # Slice chromosome sequence based on bins
        seqs = [str(rec.seq)[s:e] for s, e in zip(starts, ends)]
        # Fill GC values for bins in the chromosome
        idx = np.flatnonzero(mask)
        gc_bins[idx] = np.array(list(map(lambda s: SeqUtils.gc_fraction(s, "ignore"), seqs)))

    return gc_bins


#: Oldest GNU sort given --parallel: 8.23 fixed a deadlock of multi-threaded
#: sort introduced in 8.6 (coreutils NEWS).
_GNU_SORT_PARALLEL = (8, 23)


def _sort_implementation(version_text):
    """Identify a UNIX sort from what ``sort --version`` printed.

    Parameters
    ----------
    version_text : str
        Standard output of ``sort --version``. The first line names the
        implementation: "sort (GNU coreutils) 9.4", "sort (uutils coreutils)
        0.10.0" (the default sort of Ubuntu 26.04), "2.3-FreeBSD" or
        "2.3-Apple (199)". busybox sort has no --version and prints nothing.

    Returns
    -------
    name : str
        "GNU", "uutils", "BSD" or "unknown".
    version : tuple of int or None
        (major, minor) of GNU and uutils sort, else None.
    """
    first_line = version_text.strip().split("\n", 1)[0]
    # The version follows the package name in parentheses
    version = re.match(r"\s*(\d+)\.(\d+)", first_line.rpartition(")")[2])
    version = tuple(map(int, version.groups())) if version else None
    if "(GNU coreutils)" in first_line and version:
        return "GNU", version
    if "(uutils coreutils)" in first_line:
        return "uutils", version
    if re.match(r"\d+\.\d+-(\w*BSD|Apple)\b", first_line):
        return "BSD", None
    return "unknown", None


def _sort_parallel_ok(version_text):
    """Whether sort_pairs can give sort --parallel, warning when it cannot.

    Parameters
    ----------
    version_text : str
        Standard output of ``sort --version``.

    Returns
    -------
    bool
        True for GNU sort 8.23 or later, and for uutils sort.
    """
    name, version = _sort_implementation(version_text)
    if name == "uutils" or (name == "GNU" and version >= _GNU_SORT_PARALLEL):
        return True
    if name == "GNU":
        logger.warning(
            "GNU sort version is %s but %s or later is required for parallel sort. "
            "Sorting without --parallel.",
            ".".join(map(str, version)),
            ".".join(map(str, _GNU_SORT_PARALLEL)),
        )
    elif name == "BSD":
        # BSD sort uses a single thread when it reads stdin, as here
        logger.warning("Using BSD sort instead of GNU sort, sorting on a single thread.")
    return False


#: Environment variable that makes sort_pairs sort in Python, whatever sort
#: program is installed, when set to "python".
SORT_ENV = "HICSTUFF_SORT"


def _sort_version():
    """What ``sort --version`` prints, if sort_pairs runs the sort program.

    sort_pairs runs GNU, uutils and BSD sort. It sorts in Python when
    HICSTUFF_SORT is "python", when there is no sort on PATH, and when the
    sort on PATH is another one: busybox sort, which has none of the long
    options sort_pairs passes and sorts in memory, or Windows' sort.exe.

    Returns
    -------
    str or None
        The standard output of ``sort --version``, or None to sort in Python.
    """
    if os.environ.get(SORT_ENV, "").strip().lower() == "python":
        logger.info(f"Sorting pairs in Python ({SORT_ENV}=python).")
        return None
    if which("sort") is None:
        logger.info("No sort program on PATH: sorting pairs in Python.")
        return None
    # stderr is dropped: a sort without --version (busybox) prints its usage.
    try:
        version_text = sp.run(["sort", "--version"], stdout=sp.PIPE, stderr=sp.DEVNULL).stdout
    except OSError:
        version_text = b""
    version_text = version_text.decode(errors="replace")
    if _sort_implementation(version_text)[0] == "unknown":
        logger.info(
            "Unknown sort (sort --version printed %r): sorting pairs in Python.",
            version_text.strip().split("\n", 1)[0],
        )
        return None
    return version_text


#: Bytes that dictionary order (sort -d) ignores in the C locale: all but
#: ASCII letters, digits and blanks.
_NOT_DICTIONARY = bytes(
    c for c in range(256) if not (chr(c).isascii() and (chr(c).isalnum() or chr(c) in " \t"))
)
#: A field as sort without -t sees it: the blanks before it, then non-blanks.
_SORT_FIELD = re.compile(rb"[ \t]*[^ \t]+")
#: The number that sort -n reads at the start of a key in the C locale: after
#: blanks, a minus sign, digits and a fraction, all optional. No plus sign,
#: exponent or thousands separator.
_SORT_NUMBER = re.compile(rb"[ \t]*(-?)([0-9]*)(?:\.([0-9]*))?")
#: Memory that a line takes while Python sorts it, on top of its length: the
#: line object, its key and their references. Measured on pairs lines: 180
#: to 440 bytes, depending on the keys.
_SORT_LINE_OVERHEAD = 450
#: Most temporary files merged at once, under the limit of 256 open files
#: that macOS sets by default.
_SORT_MERGE_FILES = 100
#: Bytes in each unit of sort -S, kibibytes by default.
_SORT_UNITS = {"": 1 << 10, "b": 1, "k": 1 << 10, "m": 1 << 20, "g": 1 << 30, "t": 1 << 40}


def _sort_buffer_bytes(buffer):
    """The bytes in a sort -S size: a number and a unit (b, K, M, G or T,
    kibibytes by default), or a percentage of the physical memory."""
    match = re.fullmatch(r"\s*(\d+)([%bkmgt]?)", str(buffer), flags=re.IGNORECASE)
    if match is None:
        raise ValueError(f"Invalid sort buffer size {buffer!r}: give a number and a unit (2G).")
    number, unit = int(match.group(1)), match.group(2).lower()
    if unit == "%":
        return number * os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") // 100
    return number * _SORT_UNITS[unit]


def _sort_number(text):
    """The value that sort -n reads at the start of a key in the C locale, 0
    when there is no number."""
    sign, integer, fraction = _SORT_NUMBER.match(text).groups()
    integer = int(integer or b"0")
    fraction = (fraction or b"").rstrip(b"0")
    value = Decimal(f"{integer}.{fraction.decode()}") if fraction else integer
    return -value if sign else value


def _sort_fields(content, n_fields):
    """The first n_fields fields of a line as sort without -t sees them, each
    with the blanks before it. A field after the last is the blanks that end
    the line, and those after it are empty."""
    ends = [match.end() for match in itertools.islice(_SORT_FIELD.finditer(content), n_fields)]
    fields, start = [], 0
    for i in range(n_fields):
        end = ends[i] if i < len(ends) else len(content)
        fields.append(content[start:end])
        start = end
    return fields


def _sort_key(sort_keys):
    """A key function that orders lines as ``LC_ALL=C sort`` with these
    -kN,N options orders them: key by key, then byte by byte.

    Without -t, sort splits a line at blanks, spaces and tabs, and a field
    begins with the blanks before it: field 2 of "r1\tchr1\t5" is "\tchr1".
    A dictionary key (d) keeps letters, digits and blanks only, and a
    numeric key (n) is compared by value. The lines given to the key end
    with their newline, which no key includes.
    """
    specs = []
    for option in sort_keys:
        field, kind = re.fullmatch(r"-k(\d+),\1([dn]?)", option).groups()
        specs.append((int(field) - 1, kind))
    n_fields = max(field for field, _ in specs) + 1

    def any_line(content):
        fields = _sort_fields(content, n_fields)
        values = []
        for field, kind in specs:
            text = fields[field]
            if kind == "n":
                digits = text.lstrip(b" \t")
                values.append(int(digits) if digits.isdigit() else _sort_number(text))
            elif kind == "d":
                values.append(text.translate(None, _NOT_DICTIONARY))
            else:
                values.append(text)
        # sort's last resort: the whole lines, byte by byte
        values.append(content)
        return tuple(values)

    # Most lines hold non-blank fields between single tabs: their key is
    # built from the line split at tabs, in one expression, which is several
    # times faster than the loop above
    terms = []
    for field, kind in specs:
        text = f"p[{field}]" if field == 0 else f'(b"\\t" + p[{field}])'
        if kind == "n":
            terms.append(f"(int(p[{field}]) if p[{field}].isdigit() else _sort_number(p[{field}]))")
        elif kind == "d":
            terms.append(f"{text}.translate(None, _NOT_DICTIONARY)")
        else:
            terms.append(text)
    source = (
        "def key(line):\n"
        "    content = line[:-1]\n"
        f'    p = content.split(b"\\t", {n_fields})\n'
        f"    if (len(p) < {n_fields} or b' ' in content or b'\\t\\t' in content\n"
        "            or content[:1] == b'\\t' or content[-1:] == b'\\t'):\n"
        "        return any_line(content)\n"
        f"    return ({', '.join(terms)}, content)\n"
    )
    namespace = {
        "any_line": any_line,
        "_sort_number": _sort_number,
        "_NOT_DICTIONARY": _NOT_DICTIONARY,
    }
    exec(source, namespace)
    return namespace["key"]


def _write_sorted(paths, key, output):
    """Write the lines of sorted files to an open binary output, merged."""
    with contextlib.ExitStack() as stack:
        runs = [stack.enter_context(open(path, "rb", buffering=1 << 18)) for path in paths]
        output.writelines(heapq.merge(*runs, key=key))


def _sort_pairs_python(in_file, out_file, sort_keys, tmp_dir=None, buffer="2G"):
    """Append the pairs of in_file to out_file in the order that
    ``LC_ALL=C sort`` with sort_keys gives them, sorting in Python.

    The pairs are sorted in chunks that fit in buffer, written to temporary
    files in tmp_dir unless they all fit in one, and merged. The temporary
    files are removed, also on error.
    """
    key = _sort_key(sort_keys)
    budget = _sort_buffer_bytes(buffer)
    chunk, size, runs, created = [], 0, [], []

    def new_run():
        """A new temporary file, opened for writing, and its path."""
        fd, path = tempfile.mkstemp(prefix="hicstuff_sort_", suffix=".pairs", dir=tmp_dir)
        created.append(path)
        return os.fdopen(fd, "wb"), path

    try:
        with read_compressed(in_file, "rb") as pairs:
            for line in pairs:
                # Skip the header, and any other comment line
                if line.startswith(b"#"):
                    continue
                if not line.endswith(b"\n"):
                    line += b"\n"
                chunk.append(line)
                size += len(line) + _SORT_LINE_OVERHEAD
                if size >= budget:
                    chunk.sort(key=key)
                    run, path = new_run()
                    with run:
                        run.writelines(chunk)
                    runs.append(path)
                    chunk, size = [], 0
        chunk.sort(key=key)
        if runs and chunk:
            run, path = new_run()
            with run:
                run.writelines(chunk)
            runs.append(path)
            chunk = []
        # Merge just enough files first that the rest can be merged at once
        while len(runs) > _SORT_MERGE_FILES:
            first = min(_SORT_MERGE_FILES, len(runs) - _SORT_MERGE_FILES + 1)
            group, runs = runs[:first], runs[first:]
            run, path = new_run()
            with run:
                _write_sorted(group, key, run)
            runs.append(path)
            for done in group:
                os.remove(done)
        with open(out_file, "ab") as output:
            if runs:
                _write_sorted(runs, key, output)
            else:
                output.writelines(chunk)
    finally:
        for path in created:
            with contextlib.suppress(FileNotFoundError):
                os.remove(path)


def sort_pairs(in_file, out_file, keys, tmp_dir=None, threads=1, buffer="2G"):
    """
    Sort a pairs file in batches using UNIX sort, or in Python.

    The pairs are sorted by GNU, uutils or BSD sort when one is on PATH, and
    in Python otherwise, or when the environment variable HICSTUFF_SORT is
    "python": both write the same file. The Python sort is single-threaded.

    Parameters
    ----------
    in_file : str
        Path to the unsorted input file. It can be compressed (gzip, bgzip,
        bzip2 or zip).
    out_file : str
        Path to the sorted output file. It is not compressed.
    keys : list of str
        list of columns to use as sort keys. Each column can be one of readID,
        chr1, pos1, chr2, pos2, strand1, strand2, frag1, frag2. Key priorities
        are according to the order in the list. Chromosomes sort byte by byte
        (chr1, chr10, chr2, as pairtools sort), positions and fragments
        numerically, read names and strands in dictionary order. Read names
        that dictionary order finds equal but that differ (r.1 and r1) are
        ordered byte by byte, so that each one forms one block.
    tmp_dir : str
        Path to the directory where temporary files will be created. Defaults
        to sort's own: $TMPDIR, else /tmp (/var/tmp for BSD sort), and to
        Python's for the Python sort.
    threads : int
        Number of parallel sorting threads, given with --parallel to GNU sort
        8.23 or later and to uutils sort. Other sorts choose their own.
    buffer : str
        Buffer size used for sorting. Consists of a number and a unit.

    Raises
    ------
    ValueError
        If a key is not one of the columns above, or, for the Python sort, the
        buffer size is not one sort accepts. Nothing is written then.
    RuntimeError
        If sort fails.
    """

    # sort runs in the C locale, where a plain key compares byte by byte and
    # finds two values equal only if they are identical: each chromosome
    # then forms one block. Version order (-V) would not: from coreutils 9.1
    # it ties contig_01 with contig_1, whose pairs then interleave and cannot
    # be indexed by pairix. Byte order is also the order of pairtools sort,
    # which pairtools merge expects. Dictionary order (-d) ignores
    # punctuation, so it ties r.1 with r1, and + with -: a plain comparison
    # of the same column follows it and decides between such values.
    key_map = {
        "readID": ["-k1,1d", "-k1,1"],
        "chr1": ["-k2,2"],
        "pos1": ["-k3,3n"],
        "chr2": ["-k4,4"],
        "pos2": ["-k5,5n"],
        "strand1": ["-k6,6d", "-k6,6"],
        "strand2": ["-k7,7d", "-k7,7"],
        "frag1": ["-k8,8n"],
        "frag2": ["-k9,9n"],
    }

    # Transform column names to the corresponding sort keys, before anything
    # is written
    unknown = [k for k in keys if k not in key_map]
    if unknown:
        raise ValueError(
            f"Unknown sort key(s) {', '.join(unknown)}: choose among {', '.join(key_map)}."
        )
    sort_keys = [opt for k in keys for opt in key_map[k]]
    sort_version = _sort_version()
    if sort_version is None:
        # Before anything is written
        _sort_buffer_bytes(buffer)
    # Rewrite header with new sorting order
    header = get_pairs_header(in_file)
    with open(out_file, "w") as output:
        for line in header:
            if line.startswith("#sorted"):
                output.write("#sorted: {}\n".format("-".join(keys)))
            else:
                output.write(line + "\n")

    if sort_version is None:
        _sort_pairs_python(in_file, out_file, sort_keys, tmp_dir=tmp_dir, buffer=buffer)
        return

    # Sort pairs and append to file. The pairs are streamed to sort from
    # Python, through a decompressing reader if the input is compressed.
    # Pass the number of threads only to a sort known to accept --parallel.
    parallel_ok = _sort_parallel_ok(sort_version)
    sort_cmd = ["sort", f"-S {buffer}"] + sort_keys
    if tmp_dir is not None:
        sort_cmd.append(f"--temporary-directory={tmp_dir}")
    if parallel_ok:
        sort_cmd.append(f"--parallel={threads}")
    with read_compressed(in_file, "rb") as pairs, open(out_file, "ab") as output:
        sort_proc = sp.Popen(
            sort_cmd,
            stdin=sp.PIPE,
            stdout=output,
            bufsize=1 << 20,
            env={**os.environ, "LC_ALL": "C"},
        )
        try:
            # Skip the header, and any other comment line
            sort_proc.stdin.writelines(line for line in pairs if not line.startswith(b"#"))
        except BrokenPipeError:
            pass  # sort exited early: its exit code is checked below
        except BaseException:
            # The input could not be read: do not let sort finish on partial input
            sort_proc.kill()
            raise
        finally:
            try:
                sort_proc.stdin.close()
            except BrokenPipeError:
                pass
            sort_proc.wait()
    if sort_proc.returncode != 0:
        raise RuntimeError(
            f"Sorting pairs from {in_file} into {out_file} failed "
            f"(sort exit code {sort_proc.returncode})."
        )


def get_pairs_header(pairs):
    r"""Retrieves the header of a .pairs file and stores lines into a list.

    Parameters
    ----------
    pairs : str
        Path to the pairs file. It can be compressed (gzip, bgzip, bzip2 or
        zip).

    Returns
    -------
    header : list of str
        A list of header lines found, in the same order they appear in pairs.

    Examples
    --------
        >>> import os
        >>> from tempfile import NamedTemporaryFile
        >>> p = NamedTemporaryFile('w', delete=False)
        >>> p.writelines(["## pairs format v1.0\n", "#sorted: chr1-chr2\n", "abcd\n"])
        >>> p.close()
        >>> h = get_pairs_header(p.name)
        >>> for line in h:
        ...     print([line])
        ['## pairs format v1.0']
        ['#sorted: chr1-chr2']
        >>> os.unlink(p.name)
    """
    with read_compressed(pairs) as pairs:
        # Store header lines into a list
        header = []
        line = pairs.readline()
        while line.startswith("#"):
            header.append(line.rstrip())
            line = pairs.readline()

    return header


def reorder_fasta(genome, output=None, threshold=100000):
    """Reorder and trim a fasta file

    Sort a fasta file by record lengths, optionally trimming the smallest ones.


    Parameters
    ----------
    genome : str, file or pathlib.Path
        The genome scaffold file (or file handle)
    output : str, file or pathlib.Path
        The output file to write to
    threshold : int, optional
        The size below which scaffolds are discarded, by default 100000
    """

    if output is None:
        output = f"{os.path.splitext(genome)[0]}_renamed.fa"

    handle = SeqIO.parse(genome, "fasta")
    handle_to_write = sorted((u for u in handle if len(u) > threshold), key=len, reverse=True)
    SeqIO.write(handle_to_write, output, "fasta")


def rename_genome(genome, output=None, ambiguous=True):
    """Rename genome and slugify headers

    Rename genomes according to a simple naming scheme; this
    is mainly done to avoid special character weirdness.

    Parameters
    ----------
    genome : file, str or pathlib.Path
        The input genome to be renamed and slugify.
    output : file, str or pathlib.Path
        The output genome to be written into. Default is <base>_renamed.fa,
        where <base> is genome_in without its extension.
    ambiguous : bool
        Whether to retain ambiguous non-N bases, otherwise rename them to Ns.
        Default is True.
    """

    if output is None:
        output = f"{os.path.splitext(genome)[0]}_renamed.fa"

    with open(output, "w") as output_handle:
        for record in SeqIO.parse(genome, "fasta"):
            # Replace hyphens, tabs and whitespace with underscores
            new_record_id = record.id.replace(" ", "_")
            new_record_id = new_record_id.replace("-", "_")
            new_record_id = new_record_id.replace("\t", "_")

            # Remove anything that's weird, i.e. not alphanumeric
            # or an underscore
            new_record_id = re.sub("[^_A-Za-z0-9]+", "", new_record_id)
            header = f">{new_record_id}\n"

            new_seq = re.sub("[^ATGCatgcNn]", "N", str(record.seq))

            output_handle.write(header)
            output_handle.write(f"{new_seq}\n")


#: Files of an aligner's index, after its prefix, for each layout that the
#: aligner finds from the prefix alone: bowtie2-build writes .bt2 files, or
#: .bt2l ones for a large index; bwa index writes <prefix>.*, or <prefix>.64.*
#: with -6.
_INDEX_SUFFIXES = {
    "bowtie2": [
        [f"{name}.{ext}" for name in (".1", ".2", ".3", ".4", ".rev.1", ".rev.2")]
        for ext in ("bt2", "bt2l")
    ],
    "bwa": [
        [f"{infix}.{ext}" for ext in ("amb", "ann", "bwt", "pac", "sa")] for infix in ("", ".64")
    ],
}


def check_fasta_index(ref, mode="bowtie2"):
    """
    Checks for the existence of a bowtie2 or bwa index based on the reference
    file name.

    The index must be named after the reference path as given, which is how
    bowtie2 and bwa find it: <ref>.1.bt2 to <ref>.rev.2.bt2 (or .bt2l) for
    bowtie2, <ref>.amb, .ann, .bwt, .pac and .sa for bwa. The reference path
    is kept as is: relative (./ or ../), absolute, or with dots in its name.

    Parameters
    ----------
    ref : str or pathlib.Path
        Path to the reference genome, or to the index prefix.
    mode : str
        The alignment software used to build the index. bowtie2 or bwa. If any
        other value is given, the function returns the reference path.

    Returns
    -------
    index : str
        The bowtie2 or bwa index prefix, i.e. ref. None if no index was found
    """
    ref = os.fspath(ref)
    if mode not in _INDEX_SUFFIXES:
        return ref
    # Index files are named by appending to the prefix: it is never parsed,
    # so leading ./ or ../ and dots in names are kept.
    for suffixes in _INDEX_SUFFIXES[mode]:
        if all(os.path.isfile(ref + suffix) for suffix in suffixes):
            return ref
    return None


#: Program that builds the index of each aligner that needs one.
_INDEX_BUILDERS = {"bowtie2": "bowtie2-build", "bwa": "bwa"}


def alignment_tools(aligner, index):
    """
    Programs that aligning reads with an aligner runs.

    Parameters
    ----------
    aligner : str
        The aligner: bowtie2, bwa or minimap2.
    index : str or None
        The index prefix found by check_fasta_index, or None if the index of
        the reference must be built first.

    Returns
    -------
    list of str
        The aligner, and the program that builds its index if there is none
        (bowtie2-build for bowtie2; bwa builds its own).
    """
    tools = [aligner]
    builder = _INDEX_BUILDERS.get(aligner)
    if index is None and builder is not None and builder != aligner:
        tools.append(builder)
    return tools


def require_tools(tools, purpose, alternative=None):
    """
    Check that programs are on PATH before running them.

    Parameters
    ----------
    tools : list of str
        Names of the programs.
    purpose : str
        What they are needed for, completing "needed ...", for the error
        message.
    alternative : str or None
        What can be done instead of installing them, completing "install
        them, or ...", for the error message.

    Raises
    ------
    ImportError
        If any of the programs is not on PATH, naming the missing ones.
    """
    missing = [tool for tool in tools if which(tool) is None]
    if missing:
        msg = (
            f"{' and '.join(missing)} {'is' if len(missing) == 1 else 'are'} not "
            f"installed or not on PATH, and needed {purpose}"
        )
        if alternative is None:
            msg += "."
        else:
            msg += f": install {'it' if len(missing) == 1 else 'them'}, or {alternative}."
        logger.error(msg)
        raise ImportError(msg)


def check_is_fasta(in_file):
    """
    Checks whether input file is in fasta format.

    Parameters
    ----------
    in_file : str
        Path to the input file.

    Returns
    -------
    bool :
        True if the input is in fasta format, False otherwise
    """
    try:
        with read_compressed(in_file) as handle:
            fasta = any(SeqIO.parse(handle, "fasta"))
    # Biopython >= 1.85 raises ValueError on text that is not FASTA instead of
    # yielding no record (UnicodeDecodeError, for binary input, is a ValueError)
    except (FileNotFoundError, ValueError):
        fasta = False

    return fasta


def check_fastq_entries(in_file, threads=1) -> int:
    """
    Check how many reads are in the input fastq.

    Parameters
    ----------
    in_file : str
        Path to the input file.
    threads : int, optional
        Number of decompression threads to use when the file is gzipped and
        ``pigz`` is installed.  Values > 1 enable parallel decompression.
        Default is 1, which falls back to ``zcat``.

    Returns
    -------
    int :
        How many reads listed in the input fastq
    """

    with open(in_file, "rb") as f:
        is_gzip = f.read(2) == b"\x1f\x8b"
    if is_gzip:
        if threads > 1 and which("pigz") is not None:
            decomp = sp.Popen(
                ["pigz", "-dc", "-p", str(threads), in_file],
                stdout=sp.PIPE,
                stderr=sp.DEVNULL,
            )
        else:
            decomp = sp.Popen(["zcat", in_file], stdout=sp.PIPE, stderr=sp.DEVNULL)
        wc = sp.Popen(["wc", "-l"], stdin=decomp.stdout, stdout=sp.PIPE)
        decomp.stdout.close()
        n_lines = int(wc.communicate()[0].strip())
    else:
        n_lines = int(sp.check_output(["wc", "-l", in_file]).split()[0])
    n_reads = int(n_lines) / 4
    return int(n_reads)


def check_bam_entries(in_file):
    """
    Check how many reads are in the input bam

    Parameters
    ----------
    in_file : str
        Path to the input file.

    Returns
    -------
    int :
        How many reads listed in the input bam
    """

    return int(ps.view("-c", in_file))


#: Mate suffix at the end of a read name.
_MATE_SUFFIX = re.compile(r"/[12]$")


def strip_mate_suffix(name):
    """
    Read name without its mate suffix, ``/1`` or ``/2``, if it has one.

    Aligners do not all keep the suffix in the names they write: bwa mem
    drops it, while bowtie2 and minimap2 keep it when aligning single reads.
    Names compared without it match whichever aligner wrote them.

    Parameters
    ----------
    name : str
        Read name, as in a FASTQ or BAM file.

    Returns
    -------
    str :
        The name without a trailing ``/1`` or ``/2``.

    Examples
    --------
    >>> strip_mate_suffix("SRR1.7/1")
    'SRR1.7'
    >>> strip_mate_suffix("read/12")
    'read/12'
    """
    return _MATE_SUFFIX.sub("", name)


#: Python code that runs ``samtools <command> <args>`` with the samtools that
#: pysam bundles, given the command and its arguments on the command line:
#: like the samtools program, it reads from its standard input and writes to
#: its standard output when given "-" as file name, and exits with samtools'
#: error message if the command fails. The working directory is left out of
#: the module search path, so that no file there can stand in for a module.
#: Ctrl-C stops it at once, as it stops the samtools program: Python would
#: only raise KeyboardInterrupt once samtools returns.
_SAMTOOLS_CHILD = """\
import sys
if sys.path and sys.path[0] == "":
    del sys.path[0]
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)
import pysam
try:
    getattr(pysam, sys.argv[1])(*sys.argv[2:], catch_stdout=False)
except pysam.SamtoolsError as err:
    sys.exit(err.value)
"""


def samtools_pipe(*commands, stdin=None):
    """
    Run samtools commands with the samtools that pysam bundles, each reading
    the output of the previous one, as ``samtools view ... | samtools sort ...``
    would.

    Each command runs in its own child process, so that it can stream from and
    to a pipe: pysam runs samtools commands within the calling process, which
    has a single standard input and output.

    Parameters
    ----------
    *commands : list
        Each command: the samtools command (view, sort...) followed by its
        arguments. All but the last write to their standard output, which the
        next one reads with "-" (or no input file) as input.
    stdin : file object or None
        Standard input of the first command, such as the standard output of
        another process. It is closed in this process once handed over, so
        that the writing process stops if the first command does.

    Raises
    ------
    RuntimeError
        If any command fails, naming the failed commands.
    """
    processes = []
    try:
        for i, command in enumerate(commands):
            piped = i < len(commands) - 1
            processes.append(
                sp.Popen(
                    [sys.executable, "-c", _SAMTOOLS_CHILD, *map(str, command)],
                    stdin=stdin,
                    stdout=sp.PIPE if piped else None,
                )
            )
            if stdin is not None:
                stdin.close()
            stdin = processes[-1].stdout
        returncodes = [process.wait() for process in processes]
    except BaseException:
        # Interrupted, by Ctrl-C or otherwise: stop the commands, as
        # subprocess.run does, rather than leave them running.
        for process in processes:
            process.kill()
            process.wait()
        raise
    failed = [
        f"samtools {command[0]} (exit code {returncode})"
        for command, returncode in zip(commands, returncodes)
        if returncode != 0
    ]
    if failed:
        raise RuntimeError(f"{', '.join(failed)} failed: see samtools' error message above.")
