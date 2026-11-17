"""bwa mem through the bwa program, or in-process through pybwa.

hicstuff aligns reads with bwa (``--aligner bwa``) through one of two
backends: the ``bwa`` program, or pybwa, a dependency of hicstuff that runs
bwa's ``index`` and ``mem`` within Python
(https://github.com/fulcrumgenomics/pybwa). With pybwa, aligning reads needs
no program at all.

pybwa writes the primary records of the bwa program, read by read, once read
names are trimmed as bwa trims them. It aligns the reads of a file
``BATCH_SIZE`` at a time, however, and bwa breaks ties between
alignments of equal score with a hash of each read's position among the
reads aligned together, which restarts with every batch. Of a read cut into
two halves of equal score, the half that is primary can then differ from the
bwa program's choice: 0.02 to 0.25% of the pairs of real libraries. With -5
(mapping mode 5spm), the 5'-most half is primary and no such tie arises.

This module contains the following functions:

- resolve_backend
- load_index
- build_index
- align
- iterative_align
- align_pairs
"""

import contextlib
import importlib.metadata
import os
import subprocess as sp
import tempfile
from shutil import which

import pysam as ps

import hicstuff.io as hio
from hicstuff.log import logger

#: Environment variable that chooses the bwa backend when none is given.
BACKEND_ENV = "HICSTUFF_ALIGNER_BACKEND"
#: How bwa can run: chosen automatically, the bwa program, or pybwa.
BACKENDS = ("auto", "binary", "python")
#: Reads aligned per call to pybwa: their alignments take about 220 MB, as
#: those of the bwa program. Larger batches take more memory and break fewer
#: ties as the bwa program does (see the module's documentation).
BATCH_SIZE = 10_000
#: The bwa that pybwa 2.3 runs.
PYBWA_BWA_VERSION = "0.7.18"
#: Length of the reference above which ``bwa index`` builds the BWT with the
#: bwtsw algorithm, and below which with the is algorithm.
BWTSW_MIN_LENGTH = 50_000_000
#: Files of a bwa index, named after the reference.
INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")
#: Options of mapping mode 5spm: the 5'-most alignment of a read is its
#: primary (-5), mates are neither rescued (-S) nor paired (-P), and the other
#: alignments of a read cut in several parts are secondary (-M).
FIVE_SPM = "-5SPM"


def pybwa_version():
    """The version of pybwa, or None if it cannot be imported."""
    try:
        import pybwa  # noqa: F401
    except ImportError:
        return None
    return importlib.metadata.version("pybwa")


def describe(backend, five_spm=False):
    """How bwa mem runs with a backend, as the logs name it."""
    command = f"bwa mem {FIVE_SPM}" if five_spm else "bwa mem"
    if backend == "python":
        return f"{command} in-process, pybwa {pybwa_version()} (bwa {PYBWA_BWA_VERSION})"
    return f"{command}, bwa program"


def resolve_backend(requested=None, prefer_python=False):
    """The backend that aligns reads with bwa.

    Parameters
    ----------
    requested : str or None
        "auto", "binary" (the bwa program) or "python" (pybwa). None reads the
        environment variable HICSTUFF_ALIGNER_BACKEND, and is auto if it is
        not set.
    prefer_python : bool
        Whether auto prefers pybwa to the bwa program, as mapping mode 5spm
        does. Otherwise auto prefers the bwa program, so that a run aligns its
        reads as it did before pybwa.

    Returns
    -------
    str
        "binary" or "python".

    Raises
    ------
    ValueError
        If requested is not one of the backends.
    ImportError
        If the backend requested cannot run, or, for auto, neither can,
        naming what to install.
    """
    if requested is None:
        requested = os.environ.get(BACKEND_ENV) or "auto"
    requested = requested.strip().lower()
    if requested not in BACKENDS:
        raise ValueError(
            f"Unknown aligner backend {requested!r}: choose among {', '.join(BACKENDS)}."
        )
    found = {"binary": which("bwa") is not None, "python": pybwa_version() is not None}
    if requested == "auto":
        order = ("python", "binary") if prefer_python else ("binary", "python")
        for backend in order:
            if found[backend]:
                return backend
        msg = (
            "Neither the bwa program nor pybwa is installed, and one of them is needed to "
            "align reads with bwa: install pybwa (pip install pybwa), "
            "or the bwa program."
        )
    elif found[requested]:
        return requested
    elif requested == "binary":
        msg = (
            "bwa is not installed or not on PATH, and needed to align reads with "
            "--aligner-backend binary: install it, or use --aligner-backend python, which "
            "runs bwa in-process with pybwa."
        )
    else:
        msg = (
            "pybwa is not installed, and needed to align reads with --aligner-backend python: "
            "install it (pip install pybwa), or use --aligner-backend "
            "binary, which runs the bwa program."
        )
    logger.error(msg)
    raise ImportError(msg)


def _quiet():
    """Make bwa print errors only, as ``bwa -v 1``, also when it loads an
    index: this is a global of the bwa library."""
    from pybwa import BwaVerbosity, set_bwa_verbosity

    set_bwa_verbosity(BwaVerbosity.ERROR)


def build_index(fasta, prefix):
    """Build the bwa index of a fasta file at prefix with pybwa.

    The BWT is built by the algorithm ``bwa index`` chooses: bwtsw for
    references longer than 50 Mb, is otherwise. The index files are the same
    whatever the algorithm, which only sets the time and memory taken.
    """
    from pybwa import BwaIndex, BwaIndexBuildMethod

    with ps.FastxFile(fasta) as records:
        length = sum(len(record.sequence) for record in records)
    # pybwa 2.3 passes its methods to bwa off by one, which hicstuff pins:
    # RB2 runs bwtsw and BWTSW runs is.
    method = BwaIndexBuildMethod.RB2 if length > BWTSW_MIN_LENGTH else BwaIndexBuildMethod.BWTSW
    _quiet()
    BwaIndex.index(fasta, method=method, prefix=prefix)


def _write_sequence_dict(prefix):
    """Write the sequence dictionary that pybwa reads with the index at
    prefix, from the index itself: the @SQ lines of the bwa program's SAM
    header, with the names and lengths of the reference's sequences."""
    with open(f"{prefix}.ann") as ann:
        n_sequences = int(ann.readline().split()[1])
        lines = []
        for _ in range(n_sequences):
            # "gi name annotation", then "offset length ambiguous_bases"
            name = ann.readline().split()[1]
            length = int(ann.readline().split()[1])
            lines.append(f"@SQ\tSN:{name}\tLN:{length}\n")
    with open(f"{prefix}.dict", "w") as sequence_dict:
        sequence_dict.writelines(lines)


def load_index(genome, tmp_dir):
    """Load the bwa index of a reference with pybwa, building it first if
    there is none.

    pybwa reads an index with a sequence dictionary named after the index
    prefix with its last suffix replaced (genome.dict for genome.fa). Rather
    than write one next to the reference, where another file may have that
    name, the index files are linked into a new directory of tmp_dir, next
    to a dictionary written from the index. An index that must be built is
    built there.

    Parameters
    ----------
    genome : str
        Path to the reference fasta, or to the prefix of its bwa index.
    tmp_dir : str
        Directory in which the index directory is created.

    Returns
    -------
    pybwa.BwaIndex
        The index, loaded in memory.
    """
    from pybwa import BwaIndex

    prefix = os.path.join(hio.generate_temp_dir(tmp_dir), "index")
    found = hio.check_fasta_index(genome, mode="bwa")
    if found is None:
        build_index(genome, prefix)
    else:
        for suffix in INDEX_SUFFIXES:
            os.symlink(os.path.abspath(found + suffix), prefix + suffix)
    _write_sequence_dict(prefix)
    _quiet()
    return BwaIndex(prefix)


@contextlib.contextmanager
def loaded_index(genome, tmp_dir):
    """An index loaded by :func:`load_index` for the block, whose files are
    removed when the block ends, or an index already loaded.

    Parameters
    ----------
    genome : str or pybwa.BwaIndex
        Path to the reference or the prefix of its index, or an index.
    tmp_dir : str
        Directory in which a temporary directory holds the index files.
    """
    if not isinstance(genome, (str, os.PathLike)):
        yield genome
        return
    with tempfile.TemporaryDirectory(dir=tmp_dir) as directory:
        yield load_index(genome, directory)


def _bwa_name(name):
    """A read name as bwa writes it: without a final "/" and digit, which bwa
    takes for a mate number."""
    if len(name) > 2 and name[-2] == "/" and name[-1] in "0123456789":
        return name[:-2]
    return name


def _queries(records):
    """The reads of fastq records as bwa reads them, named as it names them."""
    return [
        ps.FastxRecord(name=_bwa_name(r.name), sequence=r.sequence, quality=r.quality)
        for r in records
    ]


def _batches(path, size):
    """The reads of a fastq file, plain or compressed, size at a time."""
    with ps.FastxFile(path) as reads:
        batch = []
        for read in reads:
            batch.append(read)
            if len(batch) == size:
                yield _queries(batch)
                batch = []
        if batch:
            yield _queries(batch)


def _options(threads, five_spm=False):
    """The options of ``bwa mem -t threads``, and with five_spm those of -5SPM."""
    from pybwa import BwaMemOptions

    options = BwaMemOptions(threads=threads)
    if five_spm:
        options.query_coord_as_primary = True
        options.skip_mate_rescue = True
        options.skip_pairing = True
        options.short_split_as_secondary = True
        # bwa mem -5 keeps the mapping quality of supplementary alignments,
        # which pybwa 2.3 does not set with query_coord_as_primary
        options.keep_mapq_for_supplementary = True
    return options.finalize()


def _header(index, command):
    """The SAM header of the alignments of an index: its sequences, and a
    @PG line for bwa run in-process with this command line."""
    header = index.header.to_dict()
    header["PG"] = [
        {
            "ID": "bwa",
            "PN": "bwa",
            "VN": f"{PYBWA_BWA_VERSION} (pybwa {pybwa_version()})",
            "CL": command,
        }
    ]
    return ps.AlignmentHeader.from_dict(header)


def _sort_by_name(unsorted, out_bam, threads):
    """Sort a BAM file by read name, as ``samtools sort -n`` does, and remove it."""
    try:
        ps.sort("-n", "-@", str(threads), "-O", "BAM", "-o", out_bam, unsorted)
    finally:
        os.remove(unsorted)


def align(reads, index, out_bam, threads=1, five_spm=False, batch_size=BATCH_SIZE):
    """Align the reads of a fastq file with bwa mem in-process, and write
    their alignments to a BAM file sorted by read name.

    The records written are those that hicstuff keeps of the bwa program's:
    all but the supplementary alignments, and with five_spm (``bwa mem
    -5SPM``) only the primary ones, since -M makes split alignments
    secondary.

    Parameters
    ----------
    reads : str
        Path to the fastq file, plain or compressed.
    index : pybwa.BwaIndex
        The index, from :func:`load_index`.
    out_bam : str
        Path to the BAM file to write.
    threads : int
        Threads of bwa mem and of BAM compression and sorting.
    five_spm : bool
        Align with -5SPM.
    batch_size : int
        Reads aligned together.
    """
    from pybwa import BwaMem

    aligner = BwaMem(index=index)
    options = _options(threads, five_spm)
    command = f"bwa mem {FIVE_SPM + ' ' if five_spm else ''}-t {threads} -v 1 <index> {reads}"
    unsorted = out_bam + ".unsorted"
    with ps.AlignmentFile(unsorted, "wb", header=_header(index, command), threads=threads) as out:
        for batch in _batches(reads, batch_size):
            for records in aligner.align(batch, options):
                for record in records:
                    if not (record.is_supplementary or (five_spm and record.is_secondary)):
                        out.write(record)
    _sort_by_name(unsorted, out_bam, threads)


def _primary(records):
    """The primary record of a read: bwa writes one per read."""
    return next(r for r in records if not (r.is_secondary or r.is_supplementary))


def _truncated(query, length):
    """A read truncated to its first length bases."""
    quality = None if query.quality is None else query.quality[:length]
    return ps.FastxRecord(name=query.name, sequence=query.sequence[:length], quality=quality)


def iterative_align(
    reads, index, out_bam, threads=1, min_len=20, min_qual=30, read_len=None, batch_size=BATCH_SIZE
):
    """Align the reads of a fastq file iteratively with bwa mem in-process,
    as :func:`hicstuff.iteralign.iterative_align` does with the bwa program,
    and write one record per read to a BAM file sorted by read name.

    Reads are truncated to min_len bases and aligned; those whose primary
    alignment is unique enough (flag 0 or 16, mapping quality of min_qual or
    more) are kept, and the others are extended by 20 bases and aligned
    again, until the reads are aligned at full length. The rounds run in
    memory, batch_size reads at a time, without temporary files.

    Parameters
    ----------
    reads : str
        Path to the fastq file, plain or compressed.
    index : pybwa.BwaIndex
        The index, from :func:`load_index`.
    out_bam : str
        Path to the BAM file to write.
    threads : int
        Threads of bwa mem and of BAM compression and sorting.
    min_len : int
        Length of the reads in the first round.
    min_qual : int
        Mapping quality from which an alignment is kept.
    read_len : int
        Length of the reads, that of the first read if None.
    batch_size : int
        Reads aligned together.
    """
    from pybwa import BwaMem

    if read_len is None:
        with ps.FastxFile(reads) as fastq:
            first = next(iter(fastq), None)
        read_len = 0 if first is None else len(first.sequence)
    if read_len > min_len:
        first_len = min_len
    else:
        logger.warning("min_len is longer than the reads. Iterative mapping will have no effect.")
        first_len = read_len
    # The truncated rounds, then the reads at full length, as iteralign
    lengths = list(range(first_len, read_len + 1, 20))
    lengths.append(lengths[-1] + 20)

    aligner = BwaMem(index=index)
    options = _options(threads)
    command = f"bwa mem -t {threads} -v 1 <index> {reads} (iterative)"
    left = [0] * len(lengths)
    total = 0
    unsorted = out_bam + ".unsorted"
    with ps.AlignmentFile(unsorted, "wb", header=_header(index, command), threads=threads) as out:
        for batch in _batches(reads, batch_size):
            total += len(batch)
            pending = batch
            for i, length in enumerate(lengths):
                last = i == len(lengths) - 1
                aligned = aligner.align([_truncated(q, length) for q in pending], options)
                still = []
                for query, records in zip(pending, aligned):
                    primary = _primary(records)
                    # Kept as iteralign keeps it, and written in the last round anyway
                    kept = primary.flag in (0, 16) and primary.mapping_quality >= min_qual
                    if kept or last:
                        out.write(primary)
                    if not kept:
                        still.append(query)
                left[i] += len(still)
                pending = still
    # The messages of iteralign, once all batches are aligned
    logger.info(f"{total} reads to parse")
    for i, length in enumerate(lengths):
        if i < len(lengths) - 1:
            again = "" if i == 0 else " again"
            logger.info(f"Truncating unaligned reads to {length}bp and mapping{again}.")
        else:
            logger.info(f"Trying to map unaligned reads at full length ({read_len}bp).")
        logger.info(f"{left[i]} reads left to map.")
    logger.info(f"{total - left[-1]} reads aligned / {total} total reads.")
    _sort_by_name(unsorted, out_bam, threads)


def align_pairs(reads1, reads2, genome, bam1, bam2, threads=1, backend="python", tmp_dir=None):
    """Align read pairs with ``bwa mem -5SPM``, and write each mate's primary
    alignments to a BAM file sorted by read name.

    With the bwa program, the pairs are aligned in one call, whose output is
    split by mate. With pybwa, each mate is aligned on its own, which gives
    the same primary alignments: -S and -P leave mates to their own
    alignments.

    Parameters
    ----------
    reads1, reads2 : str
        Paths to the fastq files of the first and second mates, in the same
        order, with the same names.
    genome : str or pybwa.BwaIndex
        The bwa index prefix, or with pybwa an index from :func:`load_index`.
    bam1, bam2 : str
        Paths to the BAM files to write, for each mate.
    threads : int
        Threads of bwa mem and of BAM compression and sorting.
    backend : str
        "binary" or "python".
    tmp_dir : str
        Directory of the index that pybwa loads, when genome is a path.

    Raises
    ------
    RuntimeError
        If the bwa program fails, for instance on mates named differently.
        Nothing is left of its output.
    """
    if backend == "python":
        with loaded_index(genome, tmp_dir) as index:
            align(reads1, index, bam1, threads=threads, five_spm=True)
            align(reads2, index, bam2, threads=threads, five_spm=True)
        return

    command = ["bwa", "mem", "-v", "1", FIVE_SPM, "-t", str(threads), genome, reads1, reads2]
    unsorted = [bam1 + ".unsorted", bam2 + ".unsorted"]
    try:
        with tempfile.TemporaryFile() as stderr:
            with sp.Popen(command, stdout=sp.PIPE, stderr=stderr) as proc:
                try:
                    with ps.AlignmentFile(proc.stdout, "r") as sam, contextlib.ExitStack() as stack:
                        outs = [
                            stack.enter_context(
                                ps.AlignmentFile(path, "wb", template=sam, threads=threads)
                            )
                            for path in unsorted
                        ]
                        for record in sam:
                            if not (record.is_secondary or record.is_supplementary):
                                outs[record.is_read2].write(record)
                except (OSError, ValueError) as error:
                    # bwa stopped before it wrote a SAM header: its error tells why
                    proc.kill()
                    reading_error = error
                else:
                    reading_error = None
            if proc.returncode != 0 or reading_error is not None:
                stderr.seek(0)
                lines = stderr.read().decode(errors="replace").strip().splitlines()
                reason = lines[-1] if lines else str(reading_error)
                raise RuntimeError(
                    f"bwa mem {FIVE_SPM} failed (exit code {proc.returncode}): {reason}"
                )
        for path, out_bam in zip(unsorted, (bam1, bam2)):
            _sort_by_name(path, out_bam, threads)
    finally:
        for path in unsorted:
            with contextlib.suppress(FileNotFoundError):
                os.remove(path)
