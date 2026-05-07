"""
ARCHIVED EXPERIMENTAL FUNCTIONS

This module contains experimental code that is not maintained or supported.
The functions in this module were research explorations that did not make it
into the stable API. Use at your own risk as functionality may change or be
removed without notice.

The code is preserved here for reference purposes only.
"""

import collections
import itertools
import string
import warnings

import numpy as np

from hicstuff import hicstuff as hcs


def despeckle_local(M, stds=2, width=2):
    """Replace outstanding values (above stds standard deviations)
    in a matrix by the average of a surrounding window of desired width.
    """

    N = np.array(M, dtype=np.float64)
    n, m = M.shape
    for i, j in itertools.product(range(width, n - width), range(width, m - width)):
        square = M[i - width : i + width, j - width : j + width]
        avg = np.average(square)
        std = np.std(square)
        if M[i, j] >= avg + stds * std:
            N[i, j] = avg
    return (N + N.T) / 2


def bin_matrix(M, subsampling_factor=3):
    """Bin either sparse or dense matrices."""

    try:
        from scipy.sparse import issparse

        if issparse(M):
            return hcs.bin_sparse(M, subsampling_factor=subsampling_factor)
        else:
            raise ImportError
    except ImportError:
        return hcs.bin_dense(M, subsampling_factor=subsampling_factor)


def bin_annotation(annotation=None, subsampling_factor=3):
    """Perform binning on genome annotations such as contig information or bin
    positions.
    """

    if annotation is None:
        annotation = np.array([])
    n = len(annotation)
    binned_positions = [annotation[i] for i in range(n) if i % subsampling_factor == 0]
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
        measurement[i - subs + 1 : i].sum() for i in range(n) if i % subs == 0 and i > 0
    ]
    return np.array(binned_measurement)


def build_pyramid(M, subsampling_factor=3):
    """Iterate over a given number of times on matrix M
    so as to compute smaller and smaller matrices with bin_dense.
    """

    subs = int(subsampling_factor)
    if subs < 1:
        raise ValueError("Subsampling factor needs to be an integer greater than 1.")
    N = [M]
    while min(N[-1].shape) > 1:
        N.append(bin_matrix(N[-1], subsampling_factor=subs))
    return N


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
    idx = [i for i in range(n - 1) if np.ceil(units[i]) < np.ceil(units[i + 1])]
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
        sequence = "".join([str(record.seq) for record in SeqIO.parse(handle, "fasta")])

    n = len(sequence)
    for i in range(0, n, window):
        portion = sequence[i : min(i + window, n)]
        yield GC_partial(portion)


def split_genome(genome, chunk_size=10000):
    """Split genome into chunks of fixed size (save the last one)."""

    chunks = []
    from Bio import SeqIO

    with open(genome) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = record.seq
            n = len(sequence)
            chunks += [str(sequence[i : min(i + chunk_size, n)]) for i in range(0, n, chunk_size)]
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
    if type(M) is not np.ndarray:
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
        ttest_rel(N[i, i - window : i], N[i, i : i + window])[0] for i in range(window, n - window)
    ]

    if circ:
        d += [
            ttest_rel(
                N[i, i - window : i],
                np.array(list(N[i, : i - n + window]) + list(N[i, i:])),
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
    if type(M) is not np.ndarray:
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
        d += [M[i:, i:].sum() + M[: n - i, n - i].sum() for i in range(n - window, n)]
    elif extrapolate:
        d += [
            M[i - window :, i - window :].sum()
            * ((2 * window + 1) ** 2.0)
            / ((2 * (n - i) + 1) ** 2.0)
            for i in range(n - window, n)
        ]

    return np.array(d)


def from_structure(structure):
    """Return contact data from a 3D structure (in pdb format)."""

    try:
        from Bio import PDB

        if isinstance(structure, str):
            p = PDB.PDBParser()
            structure = p.get_structure("S", structure)
        if isinstance(structure, PDB.Structure.Structure):
            for _ in structure.get_chains():
                atoms = [np.array(atom.get_coord()) for atom in structure.get_atoms()]
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

        n, components = scipy.sparse.csgraph.connected_components(matrix, directed=False)
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
    interactions have more weight when computing the distance matrix.
    """

    connected = largest_connected_component(matrix)
    distances = to_distance(connected, alpha)
    n, m = connected.shape
    bary = np.sum(np.triu(distances, 1)) / (n**2)  # barycenters
    d = np.array(np.sum(distances**2, 0) / n - bary)  # distances to origin
    gram = np.array(
        [(d[i] + d[j] - distances[i][j] ** 2) / 2 for i, j in itertools.product(range(n), range(m))]
    ).reshape(n, m)
    normalized = gram / np.linalg.norm(gram, "fro")

    try:
        symmetric = np.array((normalized + normalized.T) / 2, dtype=np.longfloat)  # just in case
    except AttributeError:
        symmetric = np.array((normalized + normalized.T) / 2)

    from scipy import linalg

    eigen_values, eigen_vectors = linalg.eigh(symmetric)
    if not (eigen_values >= 0).all():
        warnings.warn("Negative eigen values were found.", stacklevel=2)
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
        string.ascii_uppercase + string.ascii_lowercase + string.digits + string.punctuation
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
    print(structure)
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
            line += letters[int(contigs[indices[i]] - 1)]  # 22 chain identifier
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
            raise ValueError("Alpha parameter must be non-zero") from None

    m = np.max(distance_function(D[D != 0]))
    M = np.zeros(D.shape)
    M[D != 0] = distance_function(D[D != 0])
    M[D == 0] = m
    return M


def shortest_path_interpolation(matrix, alpha=1, strict=True):
    """Perform interpolation on a matrix's data by using ShRec's shortest-path
    procedure backwards and forwards. This replaces zeroes with corresponding
    shortest-path based counts and may have the additional effect of 'blurring'
    the matrix somewhat. If strict is set to True, only zeroes are replaced
    this way.

    Also known as Boost-Hi-C (https://www.ncbi.nlm.nih.gov/pubmed/30615061)
    """
    matrix = np.array(matrix, np.float64)
    contacts = distance_to_contact(to_distance(matrix, alpha=alpha), alpha=alpha)
    if not strict:
        return contacts
    else:
        M = np.copy(matrix)
        M[matrix == 0] = contacts[matrix == 0]
        return M


def pdb_to_structure(filename):
    """Import a structure object from a PDB file."""

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
            flattened_positions = np.array([pos for contig in positions for pos in contig])
        except TypeError:
            flattened_positions = np.array(positions)

    if (np.diff(positions) == 0).any() and 0 not in set(positions):
        warnings.warn("I detected identical consecutive nonzero values.", stacklevel=2)
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
            (np.average(np.diagonal(matrix, j)) + np.average(np.diagonal(matrix, n - j))) / 2.0
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

    max_intra_distance = max(len(contigs == u) for u in set(contigs))

    intra_contacts = []
    inter_contacts = [np.average(np.diagonal(matrix, j)) for j in range(max_intra_distance, n)]
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
    """Estimate parameters from the model described in Rippe et al., 2001."""

    n, _ = matrix.shape

    if lengths is None:
        lengths = np.abs(np.diff(positions))

    measurements, bins = [], []
    for i in range(n):
        for j in range(1, i):
            mean_length = (lengths[i] + lengths[j]) / 2.0
            if positions[i] < positions[j]:
                d = ((positions[j] - positions[i] - lengths[i]) + mean_length) / 1000.0
            else:
                d = ((positions[i] - positions[j] - lengths[j]) + mean_length) / 1000.0

            bins.append(np.abs(d))
            measurements.append(matrix[i, j])
    parameters = estimate_param_rippe(measurements, bins, init=init, circ=circ)
    print(parameters)
    return parameters[0]


def estimate_param_rippe(measurements, bins, init=None, circ=False):
    """Perform least square optimization needed for the rippe_parameters function."""

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
            + (d - 2) / (np.power((lm * x / kuhn), 2) + d)
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
                * np.exp((d - 2) / (np.power(n0, 2) + d))
            )

            norm_circ = param[3] * (
                0.53
                * (param[0] ** -3.0)
                * np.power(s0, (param[2]))
                * np.exp((d - 2) / (np.power(s0, 2) + d))
            )

            rippe = (
                param[3]
                * (
                    0.53
                    * (param[0] ** -3.0)
                    * np.power(s, (param[2]))
                    * np.exp((d - 2) / (np.power(s, 2) + d))
                )
                * norm_lin
                / norm_circ
            )

        else:
            rippe = param[3] * (
                0.53
                * (param[0] ** -3.0)
                * np.power((param[1] * x / param[0]), (param[2]))
                * np.exp((d - 2) / (np.power((param[1] * x / param[0]), 2) + d))
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
        warnings.warn("Problem in parameters estimation", stacklevel=2)
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
                [matrix[i, j] for i, j in itertools.product(range(n), range(m)) if is_inter(i, j)]
            )
            mean_trans_contacts = np.average(trans_contacts)
        else:
            mean_trans_contacts = np.average(matrix) - diagonal / len(diagonal)

        N = np.random.poisson(lam=mean_trans_contacts, size=(n, m))
        np.fill_diagonal(N, diagonal)

    elif model == "distance":
        distances = distance_diagonal_law(matrix, positions)
        N = np.array([[distances[min(abs(i - j), n)] for i in range(n)] for j in range(n)])

    elif model == "rippe":
        trans_contacts = np.array(
            [matrix[i, j] for i, j in itertools.product(range(n), range(m)) if is_inter(i, j)]
        )
        mean_trans_contacts = np.average(trans_contacts)
        kuhn, lm, slope, d, A = rippe_parameters(matrix, positions, circ=circ)

        def jc(s, frag):
            dist = s - circ * (s**2) / lengths[frag]
            computed_contacts = (
                0.53 * A * (kuhn ** (-3.0)) * (dist**slope) * np.exp((d - 2) / (dist + d))
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

            trend = lowess.lowess(np.array(range(n), dtype=np.float64), contact_sum, f=0.03)
        except ImportError:
            expected_size = int(np.amax(contact_sum) / np.average(contact_sum))
            w = min(max(expected_size, 20), 100)
            trend = np.array([np.average(contact_sum[i : min(i + w, n)]) for i in range(n)])

        cov_score = np.sqrt((trend - np.average(trend)) / np.std(trend))

        N = ((N * cov_score).T) * cov_score

    if noisy:
        if callable(noisy):
            noise_function = noisy
        else:
            noise_function = noise
        return noise_function(N)
    else:
        return N


def model_norm(matrix, positions=None, lengths=None, model="uniform", circ=False):

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
    """Remove outlier 'atoms' (aka bins) from a structure."""

    X, Y, Z = (struct[:, i] for i in range(3))

    if filtering == "sphere":
        R = (np.std(X) ** 2 + np.std(Y) ** 2 + np.std(Z) ** 2) * (n**2)
        f = (X - np.mean(X)) ** 2 + (Y - np.mean(Y)) ** 2 + (Z - np.mean(Z)) ** 2 < R

    if filtering == "cube":
        R = min(np.std(X), np.std(Y), np.std(Z)) * n
        f = np.ones(len(X))
        for C in (X, Y, Z):
            f *= np.abs(C - np.mean(C)) < R

    if filtering == "percentile":
        f = np.ones(len(X))
        for C in (X, Y, Z):
            f *= np.abs(C - np.mean(C)) < np.percentile(np.abs(C - np.mean(C)), n)

    return np.array([X[f], Y[f], Z[f]])


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
    for i, j in itertools.product(range(n), range(n)):
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
    for i, j in itertools.product(range(n), range(n)):
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
        items = list(chunk)
        el = len(items)
        positions[index : index + el] = np.arange(el) * binning
        index += el

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
        el = len(chunk)
        yield M[index : index + el, index : index + el]
        index += el
