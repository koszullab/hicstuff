#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
import zipfile
import bz2
import io
import os
import functools
import sys
import numpy as np
import pandas as pd
import collections
import subprocess as sp
from scipy.sparse import tril, triu
import pathlib
import re
from os.path import join, exists
from random import getrandbits
from scipy.sparse import coo_matrix
from Bio import SeqIO, SeqUtils
import hicstuff.hicstuff as hcs
from hicstuff.log import logger
from hicstuff.version import __version__
import scipy.stats as ss

DEFAULT_MAX_MATRIX_SHAPE = 10000

DEFAULT_FRAGMENTS_LIST_FILE_NAME = "fragments_list.txt"
DEFAULT_INFO_CONTIGS_FILE_NAME = "info_contigs.txt"
DEFAULT_SPARSE_MATRIX_FILE_NAME = "abs_fragments_contacts_weighted.txt"


def _cols_to_sparse(sparse_array, shape=None, dtype=np.float64):
    """
    Make a coordinate based sparse matrix from columns.
    Convert (3, n) shaped arrays to a sparse matrix. The fist
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
        cols_arr = np.loadtxt(mat_path, delimiter="\t", dtype=dtype)
        shape = (int(cols_arr[0, 0]), int(cols_arr[0, 1]))
    except ValueError:
        cols_arr = np.loadtxt(
            mat_path, delimiter="\t", dtype=dtype, skiprows=1
        )
        shape = None

    # Get values into an array without the header. Use the header to give size.
    sparse_mat = _cols_to_sparse(cols_arr[1:, :], shape=shape, dtype=dtype)

    if binning == "auto":
        num_bins = max(sparse_mat.shape) + 1
        subsampling_factor = num_bins // DEFAULT_MAX_MATRIX_SHAPE
    else:
        subsampling_factor = binning
    sparse_mat = hcs.bin_sparse(
        sparse_mat, subsampling_factor=subsampling_factor
    )
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
        ValueError("Sparse matrix must be in coo format")
    dtype = s_mat.dtype
    fmt = "%i" if dtype == int else "%.10e"
    sparse_arr = np.vstack([s_mat.row, s_mat.col, s_mat.data]).T

    np.savetxt(
        path,
        sparse_arr,
        header="{nrows}\t{ncols}\t{nonzero}".format(
            nrows=s_mat.shape[0], ncols=s_mat.shape[1], nonzero=s_mat.nnz
        ),
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
            "The temporary directory cannot be created in {}. "
            "Make sure you have write permission.".format(path)
        )
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
                "The cooler package is required to use {0}, please install it first".format(
                    fun.__name__
                )
            )
            raise ImportError("The cooler package is required.")
        return fun(*args, **kwargs)

    return wrapped


@_check_cooler
def add_cool_column(
    clr, column, column_name, table_name="bins", metadata={}, dtype=None
):
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
    with clr.open("r+") as c:
        if column_name in c[table_name]:
            del c[table_name][column_name]
        h5opts = dict(compression="gzip", compression_opts=6)
        c[table_name].create_dataset(
            column_name, data=column, dtype=dtype, **h5opts
        )
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
    frags.rename(
        columns={"start": "start_pos", "end": "end_pos"}, inplace=True
    )
    frags["id"] = frags.groupby("chrom", sort=False).cumcount() + 1
    # Try loading hicstuff-specific columns
    try:
        frags = frags[
            ["id", "chrom", "start_pos", "end_pos", "size", "gc_content"]
        ]
    # If absent, only load standard columns
    except KeyError:
        frags = frags[["id", "chrom", "start_pos", "end_pos"]]

    chroms["cumul_length"] = (
        chroms.length.shift(1).fillna(0).cumsum().astype(int)
    )
    n_frags = c.bins()[:].groupby("chrom", sort=False).count().start
    chroms["n_frags"] = chroms.merge(
        n_frags, right_index=True, left_on="name", how="left"
    ).start
    chroms.rename(columns={"name": "contig"}, inplace=True)
    n = int(max(np.amax(mat.bin1_id), np.amax(mat.bin2_id))) + 1
    shape = (n, n)
    mat = coo_matrix((mat["count"], (mat.bin1_id, mat.bin2_id)), shape=shape)

    return mat, frags, chroms


@_check_cooler
def save_cool(cool_out, mat, frags, metadata={}):
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


def read_compressed(filename):
    """Read compressed file

    Opens the file in read mode with appropriate decompression algorithm.

    Parameters
    ----------
    filename : str
        The path to the input file

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
    comp = file_type(filename)
    if comp == "gz":
        return gzip.open(filename, "rt")
    elif comp == "bz2":
        return bz2.open(filename, "rt")
    elif comp == "zip":
        zip_arch = zipfile.ZipFile(filename, "r")
        if len(zip_arch.namelist()) > 1:
            raise IOError(
                "Only a single fastq file must be in the zip archive."
            )
        else:
            # ZipFile opens as bytes by default, using io to read as text
            zip_content = zip_arch.open(zip_arch.namelist()[0], "r")
            return io.TextIOWrapper(zip_content, encoding="utf-8")
    else:
        return open(filename, "r")


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


def to_dade_matrix(M, annotations="", filename=None):
    """Returns a Dade matrix from input numpy matrix. Any annotations are added
    as header. If filename is provided and valid, said matrix is also saved
    as text.
    """

    n, m = M.shape
    A = np.zeros((n + 1, m + 1))
    A[1:, 1:] = M
    if not annotations:
        annotations = np.array(["" for _ in n], dtype=str)
    A[0, :] = annotations
    A[:, 0] = annotations.T
    if filename:
        try:
            np.savetxt(filename, A, fmt="%i")
            logger.info(
                "I saved input matrix in dade format as {0}".format(
                    str(filename)
                )
            )
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
        from redis import StrictRedis as redis
        import time
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
    key = "{0}|{1}#{2}#{3}".format(int(time.time()), array_dtype, m, n)

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
    output_frags=DEFAULT_SPARSE_MATRIX_FILE_NAME,
    output_dir=None,
):
    """Convert a matrix from DADE format (https://github.com/scovit/dade)
    to a graal-compatible format. Since DADE matrices contain both fragment
    and contact information all files are generated at the same time.
    """

    with open(output_matrix, "w") as sparse_file:
        sparse_file.write("id_frag_a\tid_frag_b\tn_contact")
        with open(filename) as file_handle:
            first_line = file_handle.readline()
            for row_index, line in enumerate(file_handle):
                dense_row = np.array(line.split("\t")[1:], dtype=np.int32)
                for col_index in np.nonzero(dense_row)[0]:
                    line_to_write = "{}\t{}\t{}\n".format(
                        row_index, col_index, dense_row[col_index]
                    )
                    sparse_file.write(line_to_write)

    header = first_line.split("\t")
    bin_type = header[0]
    if bin_type == '"RST"':
        logger.info("I detected fragment-wise binning")
    elif bin_type == '"BIN"':
        logger.info("I detected fixed size binning")
    else:
        logger.warning(
            (
                "Sorry, I don't understand this matrix's "
                "binning: I read {}".format(str(bin_type))
            )
        )

    header_data = [
        header_elt.replace("'", "")
        .replace('"', "")
        .replace("\n", "")
        .split("~")
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
            line_to_write = "%s\t%s\t%s\t%s\n" % (
                contig,
                length_tig,
                n_frags,
                cumul_length,
            )
            info_contigs.write(line_to_write)
            cumul_length += n_frags

    with open(output_frags, "w") as fragments_list:

        fragments_list.write(
            "id\tchrom\tstart_pos\tend_pos" "\tsize\tgc_content\n"
        )
        bogus_gc = 0.5

        for i in range(total_length):
            line_to_write = "%s\t%s\t%s\t%s\t%s\t%s\n" % (
                int(local_frag_ids[i]) + 1,
                contig_names[i],
                frag_starts[i],
                frag_ends[i],
                frag_lengths[i],
                bogus_gc,
            )
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
        chroms = (
            pd.concat([bed2d[[0, 2]], chroms_left])
            .groupby([0], sort=False)
            .max()
        )
        for chrom, size in zip(chroms.index, np.array(chroms)):
            chrom_sizes[chrom] = size[0]
    elif fragments_file is None:
        logger.warning(
            "Please be aware that not all information can be restored from a "
            "bg2 file without fixed bin size; fragments without any contact "
            "will be lost"
        )
    # Get all possible fragment chrom-positions into an array
    frag_pos = np.vstack(
        [np.array(bed2d[[0, 1, 2]]), np.array(bed2d[[3, 4, 5]])]
    )
    # Sort by position (least important, col 1)
    frag_pos = frag_pos[frag_pos[:, 1].argsort(kind="mergesort")]
    # Then by chrom (most important, col 0)
    frag_pos = frag_pos[frag_pos[:, 0].argsort(kind="mergesort")]
    # Get unique names for fragments (chrom+pos)
    ordered_frag_pos = (
        pd.DataFrame(frag_pos).drop_duplicates().reset_index(drop=True)
    )
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
            n_bins = size // bin_size
            chrom_frags.append(
                pd.DataFrame(
                    {
                        "id": range(1, n_bins + 1),
                        "chrom": np.repeat(chrom, n_bins),
                        "start_pos": range(0, size, bin_size),
                        "end_pos": range(bin_size, size + bin_size, bin_size),
                    }
                )
            )
        frags = pd.concat(chrom_frags, axis=0).reset_index(drop=True)
        frags.insert(
            loc=3, column="size", value=frags.end_pos - frags.start_pos
        )
    # If None available, guess fragments indices from bedgraph (potentially wrong)
    else:
        frag_map = {
            (v[0], v[1]): i
            for i, v in ordered_frag_pos.iloc[:, [0, 1]].iterrows()
        }
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
    mat = coo_matrix(
        (contacts, (frag_id_a, frag_id_b)), shape=(n_frags, n_frags)
    )

    # Get size of each chromosome in basepairs
    chromsizes = frags.groupby("chrom", sort=False).apply(
        lambda x: np.int64(max(x.end_pos))
    )
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


def flexible_hic_loader(
    mat, fragments_file=None, chroms_file=None, quiet=False
):
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
    if hic_format == "bg2":
        # Use the frags file to define bins if available
        if fragments_file is not None:
            mat, frags, chroms = load_bedgraph2d(
                mat, fragments_file=fragments_file
            )
        else:
            # Guess if bin size is fixed based on MAD
            bg2 = pd.read_csv(mat, sep="\t")
            sizes = np.array(bg2.iloc[:, 2] - bg2.iloc[:, 1])
            size_mad = ss.median_absolute_deviation(sizes)
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
                mat, frags, chroms = load_bedgraph2d(
                    mat, bin_size=int(np.median(sizes))
                )

    elif hic_format == "graal":
        try:
            mat = load_sparse_matrix(mat)
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
    if (
        mat.endswith(".cool")
        or mat.count(".mcool::/") == 1
        or mat.count(".cool::/") == 1
    ):
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
    mat, out_prefix, frags=None, chroms=None, hic_fmt="graal", quiet=False,
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
            frags.to_csv(out_prefix + ".frag.tsv", sep="\t", index=False)
        except AttributeError:
            if not quiet:
                logger.warning(
                    "Could not create fragments_list.txt from input files"
                )
        try:
            chroms.to_csv(out_prefix + ".chr.tsv", sep="\t", index=False)
        except AttributeError:
            if not quiet:
                logger.warning(
                    "Could not create info_contigs.txt from input files"
                )
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
        raise ValueError("Unknown output format: {0}".format(hic_fmt))


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

    mat_df = pd.DataFrame(
        {"row": mat.row, "col": mat.col, "data": mat.data.astype(int)}
    )
    # Merge fragments with matrix based on row indices to annotate rows
    merge_mat = mat_df.merge(
        frags, left_on="row", right_index=True, how="left", suffixes=("", "")
    )
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
    bg2 = merge_mat.loc[
        :, ["chr1", "start1", "end1", "chr2", "start2", "end2", "data"]
    ]
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
        gc_bins[idx] = np.array(list(map(SeqUtils.GC, seqs))) / 100.0

    return gc_bins


def sort_pairs(in_file, out_file, keys, tmp_dir=None, threads=1, buffer="2G"):
    """
    Sort a pairs file in batches using UNIX sort.

    Parameters
    ----------
    in_file : str
        Path to the unsorted input file
    out_file : str
        Path to the sorted output file.
    keys : list of str
        list of columns to use as sort keys. Each column can be one of readID,
        chr1, pos1, chr2, pos2, frag1, frag2. Key priorities are according to
        the order in the list.
    tmp_dir : str
        Path to the directory where temporary files will be created. Defaults
        to current directory.
    threads : int
        Number of parallel sorting threads.
    buffer : str
        Buffer size used for sorting. Consists of a number and a unit.
    """
    # TODO: Write a pure python implementation to drop GNU coreutils depencency,
    # could be inspired from: https://stackoverflow.com/q/14465154/8440675

    # Check if UNIX sort version supports parallelism
    parallel_ok = True
    sort_ver = sp.Popen(["sort", "--version"], stdout=sp.PIPE)
    sort_ver = (
        sort_ver.communicate()[0]
        .decode()
        .split("\n")[0]
        .split(" ")[-1]
        .split(".")
    )
    # If so, specify threads, otherwise don't mention it in the command line
    try:
        sort_ver = list(map(int, sort_ver))
        if sort_ver[0] < 8 or (sort_ver[0] == 8 and sort_ver[1] < 23):
            logger.warning(
                "GNU sort version is {0} but >8.23 is required for parallel "
                "sort. Sorting on a single thread.".format(
                    ".".join(map(str, sort_ver))
                )
            )
            parallel_ok = False
    # BSD sort has a different format and will throw error upon parsing. It does
    # not support parallel processes anyway.
    except ValueError:
        logger.warning(
            "Using BSD sort instead of GNU sort, sorting on a single thread."
        )
        parallel_ok = False

    key_map = {
        "readID": "-k1,1d",
        "chr1": "-k2,2V",
        "pos1": "-k3,3n",
        "chr2": "-k4,4V",
        "pos2": "-k5,5n",
        "strand1": "-k6,6d",
        "strand2": "-k7,7d",
        "frag1": "-k8,8n",
        "frag2": "-k9,9n",
    }

    # transform column names to corresponding sort keys
    try:
        sort_keys = map(lambda k: key_map[k], keys)
    except KeyError:
        print("Unkown column name.")
        raise
    # Rewrite header with new sorting order
    header = get_pairs_header(in_file)
    with open(out_file, "w") as output:
        for line in header:
            if line.startswith("#sorted"):
                output.write("#sorted: {0}\n".format("-".join(keys)))
            else:
                output.write(line + "\n")

    # Sort pairs and append to file.
    with open(out_file, "a") as output:
        grep_proc = sp.Popen(["grep", "-v", "^#", in_file], stdout=sp.PIPE)
        sort_cmd = ["sort", "-S %s" % buffer] + list(sort_keys)
        if tmp_dir is not None:
            sort_cmd.append("--temporary-directory={0}".format(tmp_dir))
        if parallel_ok:
            sort_cmd.append("--parallel={0}".format(threads))
        sort_proc = sp.Popen(sort_cmd, stdin=grep_proc.stdout, stdout=output)
        sort_proc.communicate()


def get_pairs_header(pairs):
    r"""Retrieves the header of a .pairs file and stores lines into a list.

    Parameters
    ----------
    pairs : str or file object
        Path to the pairs file.

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
    # Open file if needed
    with open(pairs, "r") as pairs:
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
        output = "{}_renamed.fa".format(genome.split(".")[0])

    handle = SeqIO.parse(genome, "fasta")
    handle_to_write = sorted(
        (len(u) for u in handle if len(u) > threshold), reverse=True
    )
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
        output = "{}_renamed.fa".format(genome.split(".")[0])

    with open(output, "w") as output_handle:
        for record in SeqIO.parse(genome, "fasta"):

            # Replace hyphens, tabs and whitespace with underscores
            new_record_id = record.id.replace(" ", "_")
            new_record_id = new_record_id.replace("-", "_")
            new_record_id = new_record_id.replace("\t", "_")

            # Remove anything that's weird, i.e. not alphanumeric
            # or an underscore
            new_record_id = re.sub("[^_A-Za-z0-9]+", "", new_record_id)
            header = ">{}\n".format(new_record_id)

            new_seq = re.sub("[^ATGCatgcNn]", "N", str(record.seq))

            output_handle.write(header)
            output_handle.write("{}\n".format(new_seq))


def check_fasta_index(ref, mode="bowtie2"):
    """
    Checks for the existence of a bowtie2 or bwa index based on the reference
    file name.

    Parameters
    ----------
    ref : str
        Path to the reference genome.
    mode : str
        The alignment software used to build the index. bowtie2 or bwa. If any
        other value is given, the function returns the reference path.

    Returns
    -------
    index : str
        The bowtie2 or bwa index basename. None if no index was found
    """
    ref = pathlib.Path(ref)
    if mode == "bowtie2":
        # Bowtie2 should have 6 index files
        bt2_idx_files = list(ref.parent.glob("{}*bt2*".format(ref.name)))
        index = None if len(bt2_idx_files) < 6 else bt2_idx_files
    elif mode == "bwa":
        refdir = str(ref.parent)
        refdir_files = os.listdir(refdir)
        bwa_idx_files = [
            join(refdir, f)
            for f in refdir_files
            if re.search(r".*\.(sa|pac|bwt|ann|amb)$", f)
        ]
        index = None if len(bwa_idx_files) < 5 else bwa_idx_files
    else:
        index = [ref]
    if index is not None:
        # Convert the PosixPath objects to strings and get the longest common prefix to obtain
        # index basename (without the dot)
        index = os.path.commonprefix(list(map(str, index))).strip(".")
    return index


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
        with open(in_file, "r") as handle:
            fasta = any(SeqIO.parse(handle, "fasta"))
    except FileNotFoundError:
        fasta = False

    return fasta
