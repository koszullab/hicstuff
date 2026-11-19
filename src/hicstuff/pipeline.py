"""
Handle generation of graal-compatible contact maps from fastq files.
cmdoret, 20190322
"""

import csv
import glob
import gzip
import logging
import multiprocessing
import os
import pathlib
import re
import shlex
import shutil as st
import subprocess as sp
import sys
import time
from datetime import datetime
from os.path import join
from shutil import which

import click
import cooler
import matplotlib
import numpy as np
import pandas as pd
import pysam as ps
from Bio import SeqIO
from dateutil.relativedelta import relativedelta

import hicstuff.bwamem as hbm
import hicstuff.cutsite as hcc
import hicstuff.digest as hcd
import hicstuff.distance_law as hcdl
import hicstuff.filter as hcf
import hicstuff.io as hio
import hicstuff.iteralign as hci
import hicstuff.log as hcl
import hicstuff.stats as hcs
from hicstuff import __version__
from hicstuff.log import logger


def align_reads(
    reads,
    genome,
    out_bam,
    tmp_dir=None,
    threads=1,
    aligner="bowtie2",
    iterative=False,
    min_qual=30,
    read_len=None,
    backend=None,
):
    """
    Select and call correct alignment method and generate logs accordingly.
    Alignments are filtered so that there at most one alignment per read.
    Alignments are sorted by read name in the output file.

    Parameters
    ----------
    reads : str
        Path to the fastq file with Hi-C reads.
    genome : str or pybwa.BwaIndex
        Path to the genome bowtie2/bwa index prefix if using bowtie2 or bwa, or to the
        fasta if using minimap2. With bwa in-process, an index loaded by
        hicstuff.bwamem.load_index also works, so that it is loaded once for
        several files; a path is loaded, or indexed if it has no index.
    out_bam : str
        Path to the output BAM file containing mapped Hi-C reads.
    tmp_dir : str
        Path where temporary files are stored.
    threads : int
        Number of threads to run alignments in parallel.
    aligner : bool
        Use minimap2 or bwa instead of bowtie2.
    iterative : bool
        Wether to use the iterative mapping procedure (truncating reads and
        extending iteratively)
    min_qual : int
        Minimum mapping quality required to keep an alignment during iterative
        mapping.
    read_len : int
        Maximum read length to expect in the fastq file. Optionally used in iterative
        alignment mode. Estimated from the first read by default. Useful if input fastq
        is a composite of different read lengths.
    backend : str
        How bwa runs: "binary" (the bwa program), "python" (in-process with
        pybwa) or "auto", for which see hicstuff.bwamem.resolve_backend. None
        reads the environment variable HICSTUFF_ALIGNER_BACKEND.
    """
    if tmp_dir is None:
        tmp_dir = os.getcwd()
    if aligner == "bwa" and hbm.resolve_backend(backend) == "python":
        with hbm.loaded_index(genome, tmp_dir) as index:
            if iterative:
                hbm.iterative_align(
                    reads, index, out_bam, threads=threads, min_qual=min_qual, read_len=read_len
                )
            else:
                hbm.align(reads, index, out_bam, threads=threads)
        return
    tmp_bam = out_bam + ".tmp"

    if iterative:
        iter_tmp_dir = hio.generate_temp_dir(tmp_dir)
        hci.iterative_align(
            reads,
            tmp_dir=iter_tmp_dir,
            ref=genome,
            n_cpu=threads,
            bam_out=tmp_bam,
            min_qual=min_qual,
            aligner=aligner,
            read_len=read_len,
        )
        st.rmtree(iter_tmp_dir)
        # Remove supplementary alignments
        ps.view(
            *("-F", "2048", "-h", "-@", str(threads), "-O", "BAM", tmp_bam, "-o", out_bam),
            catch_stdout=False,
        )
    else:
        hio.require_tools(
            [aligner if aligner in ("minimap2", "bwa") else "bowtie2"], "to align reads"
        )
        if aligner == "minimap2":
            map_cmd = "minimap2 -2 -t {threads} -ax sr {fasta} {fastq} > {sam}"
        elif aligner == "bwa":
            map_cmd = "bwa mem -t {threads} -v 1 {index} {fastq} > {sam}"
        else:
            map_cmd = "bowtie2 --very-sensitive-local -p {threads} -x {index} -U {fastq} > {sam}"
        map_args = {
            "threads": threads,
            "sam": tmp_bam,
            "fastq": reads,
            "fasta": genome,
            "index": genome,
        }
        sp.call(map_cmd.format(**map_args), shell=True)
        # Remove supplementary alignments and sort reads by name
        hio.samtools_pipe(
            ["view", "-F", "2048", "-h", "-@", threads, tmp_bam],
            ["sort", "-n", "-@", threads, "-o", out_bam, "-"],
        )
    os.remove(tmp_bam)


#: Digit runs in a read name, which ``samtools sort -n`` compares by value.
_DIGIT_RUNS = re.compile(r"([0-9]+)")


def _natural_key(name):
    """Sort key ordering read names as ``samtools sort -n`` does.

    samtools compares names character by character, except that runs of
    digits found at the same place in both names are compared by value
    (``strnum_cmp`` in samtools' ``bam_sort.c``): ``r9 < r10``, although
    ``"r9" > "r10"`` as strings. Leading zeros do not count: ``r010`` and
    ``r10`` are equal under this key, as under samtools 1.17 and later.
    samtools 1.16 and earlier put ``r010`` before ``r10``; only names whose
    numbers are zero-padded to different widths are ordered differently.

    Parameters
    ----------
    name : str
        Read name.

    Returns
    -------
    tuple
        Text runs at even positions and the values of digit runs at odd ones.
        A text run followed by digits ends with ``"0"``: samtools compares a
        digit with a non-digit character as characters, and every digit
        compares alike against a non-digit.

    Examples
    --------
    >>> sorted(["r10", "r9", "r1"], key=_natural_key)
    ['r1', 'r9', 'r10']
    >>> _natural_key("r010") == _natural_key("r10")
    True
    """
    parts = _DIGIT_RUNS.split(name)
    key = []
    for i in range(1, len(parts), 2):
        key.append(parts[i - 1] + "0")
        key.append(int(parts[i]))
    key.append(parts[-1])
    return tuple(key)


def _name_order_key(bam):
    """Sort key of the read name order of a BAM file sorted by name:
    :func:`_natural_key`, or ``str`` for lexicographical order.

    ``samtools sort -n`` sorts names in natural order, ``samtools sort -N``
    (samtools 1.19 and later) in lexicographical order, which it records as
    ``SS:queryname:lexicographical`` in the ``@HD`` line.
    """
    first_line = str(bam.header).split("\n", 1)[0]
    if first_line.startswith("@HD") and "SS:queryname:lexicographical" in first_line.split("\t"):
        return str
    return _natural_key


class _NameOrderError(Exception):
    """A BAM file is not sorted in the read name order it is walked in."""


class _NameSortedBam:
    """Walk an open BAM file, found at ``path``, read by read.

    Keeps the current read, its name without ``/1`` or ``/2`` suffix, whether
    it passes ``min_qual``, and, once :meth:`order_key` has computed it, the
    ``name_key`` of that name followed by ``key_suffix``. ``counts`` tallies
    the reads read ("total") and those passing ``min_qual`` ("mapped").
    """

    def __init__(self, bam, path, name_key, min_qual, counts):
        self.bam = bam
        self.path = path
        self.name_key = name_key
        self.key_suffix = ""
        self.min_qual = min_qual
        self.counts = counts
        self.last_key = None
        self.advance()

    def advance(self):
        """Move to the next read; ``read`` and ``name`` are None at the end."""
        self.read = next(self.bam, None)
        self.key = None
        if self.read is None:
            self.name, self.passed = None, False
            return
        self.name = hio.strip_mate_suffix(self.read.query_name)
        self.passed = self.read.mapping_quality >= self.min_qual
        self.counts["total"] += 1
        self.counts["mapped"] += self.passed

    def order_key(self):
        """Sort key of the current read name, which must not come before the
        last one computed in this file: :class:`_NameOrderError` otherwise."""
        if self.key is None:
            self.key = self.name_key(self.name + self.key_suffix)
            if self.last_key is not None and self.key < self.last_key:
                raise _NameOrderError(self.path)
            self.last_key = self.key
        return self.key


def _write_pairs(bam1, bam2, out_pairs, info_contigs, name_key, min_qual, five_prime=False):
    """Write the pairs of two BAM files sorted by read name in the order of
    ``name_key`` (see :func:`bam2pairs`), walking both files together.

    Mates are matched by name without ``/1`` or ``/2`` suffix. When the names
    differ, the read whose name comes first in the order of ``name_key`` has
    no mate left in the other file and is skipped. Both files compare a name
    alike: without suffix, or with ``/1`` if both files have suffixes, so that
    names are compared as they were sorted (``r1.5/1`` comes before ``r1/1``).
    Each name compared is checked to come no earlier than the previous one
    compared in the same file: had a read been skipped while its mate was
    still ahead, the mate would break that order when reached, so
    :class:`_NameOrderError` is raised rather than the pair lost.

    Returns
    -------
    dict
        Numbers of reads ("total"), of reads with MAPQ >= min_qual ("mapped")
        and of reads without a mate in the other file ("unmatched").
    """
    # Generate header lines
    format_version = "## pairs format v1.0\n"
    sorting = "#sorted: readID\n"
    cols = "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n"
    # Chromosome order will be identical in info_contigs and pair files
    # UP031 Use format specifiers instead of percent format
    chroms = pd.read_csv(info_contigs, sep="\t").apply(
        lambda x: f"#chromsize: {x.contig} {x.length}\n", axis=1
    )
    counts = {"total": 0, "mapped": 0, "unmatched": 0}
    with ps.AlignmentFile(bam1, "rb") as file1, ps.AlignmentFile(bam2, "rb") as file2:
        forward = _NameSortedBam(file1, bam1, name_key, min_qual, counts)
        reverse = _NameSortedBam(file2, bam2, name_key, min_qual, counts)
        # Compare names as both files sorted them: with a /1 suffix if both
        # files have suffixes, since r1.5/1 comes before r1/1 but r1 before r1.5.
        if all(
            bam.read is not None and bam.read.query_name != bam.name for bam in (forward, reverse)
        ):
            forward.key_suffix = reverse.key_suffix = "/1"
        with open(out_pairs, "w") as pairs:
            pairs.writelines([format_version, sorting, cols] + chroms.tolist())
            pairs_writer = csv.writer(pairs, delimiter="\t")

            def position(end):
                # 1-based: the aligned 5' end of a reverse-strand mate is the
                # last base of its alignment, its leftmost base the first one
                if five_prime and end.is_reverse:
                    return end.reference_end
                return end.reference_start + 1

            def write_pair(name, end1, end2):
                pos1, pos2 = position(end1), position(end2)
                # Flipping to get upper triangle
                if (
                    end1.reference_id == end2.reference_id and pos1 > pos2
                ) or end1.reference_id > end2.reference_id:
                    end1, end2, pos1, pos2 = end2, end1, pos2, pos1
                pairs_writer.writerow(
                    [
                        name,
                        end1.reference_name,
                        pos1,
                        end2.reference_name,
                        pos2,
                        "-" if end1.is_reverse else "+",
                        "-" if end2.is_reverse else "+",
                    ]
                )

            # Iterate on both BAM simultaneously, skipping reads whose mate is
            # missing from the other file, until both files are exhausted.
            while forward.read is not None or reverse.read is not None:
                if forward.name == reverse.name:
                    # Keep only pairs where both reads have good quality
                    if forward.passed and reverse.passed:
                        write_pair(forward.name, forward.read, reverse.read)
                    forward.advance()
                    reverse.advance()
                    continue
                key1 = None if forward.read is None else forward.order_key()
                key2 = None if reverse.read is None else reverse.order_key()
                if key1 is not None and key1 == key2:
                    # Different names equal in that order (numbers zero-padded
                    # differently, r010 and r10) come in no set order: match
                    # every read with that key by name.
                    tied = []
                    for bam in (forward, reverse):
                        reads = {}
                        while bam.read is not None and bam.order_key() == key1:
                            reads.setdefault(bam.name, []).append((bam.read, bam.passed))
                            bam.advance()
                        tied.append(reads)
                    for name, ends1 in tied[0].items():
                        ends2 = tied[1].pop(name, [])
                        counts["unmatched"] += abs(len(ends1) - len(ends2))
                        for (end1, passed1), (end2, passed2) in zip(ends1, ends2):
                            if passed1 and passed2:
                                write_pair(name, end1, end2)
                    counts["unmatched"] += sum(len(ends) for ends in tied[1].values())
                    continue
                # Skip the read that comes first: its mate, if any, would have
                # come before the other read.
                counts["unmatched"] += 1
                if key2 is None or (key1 is not None and key1 < key2):
                    forward.advance()
                else:
                    reverse.advance()
    return counts


def bam2pairs(
    bam1, bam2, out_pairs, info_contigs, min_qual=30, tmp_dir=None, threads=1, five_prime=False
):
    """
    Make a .pairs file from two Hi-C bam files sorted by read names.
    The Hi-C mates are matched by read identifier, without /1 or /2 suffix,
    walking both files in the order they are sorted in: ``samtools sort -n``
    (natural) order, or ``samtools sort -N`` (lexicographical) order if both
    headers say so. Reads whose mate is missing from the other file are
    skipped. Files that turn out not to be in that order are walked in the
    other one, as other tools (Picard, sambamba sort -n) sort names in
    lexicographical order without saying so. Files in neither order are
    sorted again with ``samtools sort -n`` in a temporary directory, so that
    no pair is lost.
    Pairs where at least one reads maps with MAPQ below min_qual threshold are
    discarded. Pairs are sorted by readID and stored in upper triangle (first
    pair higher).

    Parameters
    ----------
    bam1 : str
        Path to the name-sorted BAM file with aligned Hi-C forward reads.
    bam2 : str
        Path to the name-sorted BAM file with aligned Hi-C reverse reads.
    out_pairs : str
        Path to the output space-separated .pairs file with columns
        readID, chr1 pos1 chr2 pos2 strand1 strand2
    info_contigs : str
        Path to the info contigs file, to get info on chromosome sizes and order.
    min_qual : int
        Minimum mapping quality required to keep a Hi-C pair.
    tmp_dir : str
        Directory in which to sort BAM files that are in neither read name
        order. Defaults to the directory of out_pairs.
    threads : int
        Number of threads to sort them with.
    five_prime : bool
        Whether each mate's position is its aligned 5' end, the first base of
        the read that aligns, at the end of the sequenced molecule, as
        pairtools records it. Otherwise it is the leftmost base of the
        alignment: the 5' end of a forward-strand mate, but the 3' end of a
        reverse-strand one.
    """
    with ps.AlignmentFile(bam1, "rb") as forward, ps.AlignmentFile(bam2, "rb") as reverse:
        orders = {_name_order_key(forward), _name_order_key(reverse)}
    # Mates are matched by walking both files in the order they are sorted in
    declared = orders.pop() if len(orders) == 1 else _natural_key
    order_names = {_natural_key: "samtools sort -n", str: "samtools sort -N"}
    for name_key in (declared, str if declared is _natural_key else _natural_key):
        try:
            counts = _write_pairs(
                bam1, bam2, out_pairs, info_contigs, name_key, min_qual, five_prime
            )
            break
        except _NameOrderError as err:
            # Pairs may have been missed: start again
            os.remove(out_pairs)
            unsorted = str(err)
            logger.info(
                "%s is not sorted by read name as %s sorts it.", unsorted, order_names[name_key]
            )
    else:
        sort_dir = hio.generate_temp_dir(tmp_dir or os.path.dirname(os.path.abspath(out_pairs)))
        logger.warning(
            "%s is in the read name order of neither samtools sort -n nor samtools sort -N: "
            "sorting both BAM files with samtools sort -n in %s.",
            unsorted,
            sort_dir,
        )
        sorted_bams = [join(sort_dir, "for.bam"), join(sort_dir, "rev.bam")]
        try:
            for bam, sorted_bam in zip((bam1, bam2), sorted_bams):
                ps.sort("-n", "-@", str(threads), "-o", sorted_bam, bam)
            counts = _write_pairs(
                *sorted_bams, out_pairs, info_contigs, _natural_key, min_qual, five_prime
            )
        except _NameOrderError as still_unsorted:
            os.remove(out_pairs)
            bam = bam1 if str(still_unsorted) == sorted_bams[0] else bam2
            raise ValueError(
                f"Cannot match the mates of {bam1} and {bam2}: the read names of "
                f"{bam}, sorted with samtools sort -n (samtools "
                f"{ps.__samtools_version__}), are still not in the order expected."
            ) from None
        finally:
            st.rmtree(sort_dir)
    if counts["unmatched"] > 0:
        logger.warning("%d reads were only present in one BAM file.", counts["unmatched"])
    logger.info(
        f"{100 * round(counts['mapped'] / counts['total'], 2):.2f}% reads (single ends) mapped with Q >= {min_qual} ({counts['mapped']}/{counts['total']})"
    )


def generate_log_header(log_path, input1, input2, genome, enzyme, aligner=None, cmd=None):
    hcl.set_file_handler(log_path, formatter=logging.Formatter(""))
    logger.info("## hicstuff: v%s log file", __version__)
    logger.info("## date: %s", time.strftime("%Y-%m-%d %H:%M:%S"))
    if cmd is not None:
        logger.info("## cmd: %s", cmd)
    logger.info("## enzyme: %s", str(enzyme))
    logger.info("## input1: %s ", input1)
    logger.info("## input2: %s", input2)
    logger.info("## ref: %s", genome)
    if aligner is not None:
        logger.info("## aligner: %s", aligner)
    logger.info("---")
    hcl.set_file_handler(log_path, formatter=hcl.logfile_formatter)


#: Keys of the coordinate sort that comes before PCR duplicate removal. The
#: strands are keys too: filter_pcr_dup compares each pair with the previous
#: one only, so a pair at the same positions on other strands must not sort
#: between two duplicates.
COORD_SORT_KEYS = ["chr1", "pos1", "chr2", "pos2", "strand1", "strand2"]

# The formats of the contact matrix, and none for no matrix.
MATRIX_FORMATS = ("graal", "bg2", "cool", "none")

# From this many contigs in the cool matrix, the pipeline warns that it is slow
# to write: cooler tries to store the contig of each bin as an HDF5 enum of one
# member per contig, which takes a time that grows with the square of their
# number, for every cooler it writes (about 2 s at 32,000 contigs, 30 s at
# 117,000), then is refused past a few thousand members, and stores integers.
MANY_CONTIGS = 50_000


def _excluded(exclude):
    """The names of the contigs that exclude leaves out of the cool matrix:
    None, or comma-separated names, as pairs2cool reads them; or a list of
    names, which the formats other than cool, which do not read it, take."""
    if not exclude:
        return set()
    return set(exclude.split(",")) if isinstance(exclude, str) else set(exclude)


def filter_pcr_dup(pairs_idx_file, filtered_file):
    """
    Filter out PCR duplicates from a pairs file sorted by COORD_SORT_KEYS,
    using overrrepresented exact coordinates. If multiple fragments have two
    reads with the exact same coordinates, only one of those fragments is kept.
    Parameters
    ----------
    pairs_idx_file : str
        Path to an indexed pairs file containing the Hi-C reads.
    filtered_file : str
        Path to the output pairs file after removing duplicates.
    """
    # Keep count of how many reads are filtered
    filter_count = 0
    reads_count = 0
    # Store header lines
    header = hio.get_pairs_header(pairs_idx_file)
    with open(pairs_idx_file) as pairs, open(filtered_file, "w") as filtered:
        # Copy header lines to filtered file
        for head_line in header:
            filtered.write(head_line + "\n")
            next(pairs)

        # Use csv methods to easily access columns
        paircols = [
            "readID",
            "chr1",
            "pos1",
            "chr2",
            "pos2",
            "strand1",
            "strand2",
            "frag1",
            "frag2",
        ]
        # Columns used for comparison of coordinates
        coord_cols = [col for col in paircols if col != "readID"]
        pair_reader = csv.DictReader(pairs, delimiter="\t", fieldnames=paircols)
        filt_writer = csv.DictWriter(filtered, delimiter="\t", fieldnames=paircols)

        # Initialise a variable to store coordinates of reads in previous pair
        prev = {k: 0 for k in paircols}
        for pair in pair_reader:
            reads_count += 1
            # If coordinates are the same as before, skip pair
            if all(pair[pair_var] == prev[pair_var] for pair_var in coord_cols):
                filter_count += 1
                continue
            # Else write pair and store new coordinates as previous
            else:
                filt_writer.writerow(pair)
                prev = pair
        logger.info(
            f"{100 * round(filter_count / reads_count, 3):.1f}% PCR duplicates have been filtered out ({filter_count}/{reads_count} pairs)"
        )
        logger.info(f"{reads_count - filter_count} pairs remaining after removing PCR duplicates")


def sort_and_compress_pairs(pairs_file, out_file, tmp_dir=None, threads=1):
    """
    Sort a pairs file by chr1, chr2, pos1, pos2 (the standard order for
    pairs files) and compress it in the BGZF format, so that pairix can index
    it: with bgzip, on several threads, if it is installed, otherwise with
    pysam, on one thread.

    Parameters
    ----------
    pairs_file : str
        Path to the input, uncompressed pairs file.
    out_file : str
        Path to the output sorted and compressed pairs file.
    tmp_dir : str
        Directory where temporary files will be created.
    threads : int
        Number of threads used for sorting and bgzip compression.
    """
    sorted_tmp = out_file + ".tmp"
    hio.sort_pairs(
        pairs_file,
        sorted_tmp,
        keys=["chr1", "chr2", "pos1", "pos2"],
        tmp_dir=tmp_dir,
        threads=threads,
    )
    try:
        if which("bgzip") is not None:
            with open(out_file, "wb") as out:
                sp.run(["bgzip", "-c", "-@", str(threads), sorted_tmp], stdout=out, check=True)
        else:
            ps.tabix_compress(sorted_tmp, out_file, force=True)
    finally:
        os.remove(sorted_tmp)


def _cooler_command(*names):
    """
    Return a `cooler` command by name, e.g. ``_cooler_command("cload", "pairs")``.

    The command is looked up through cooler's command-line entry point
    (``cooler.cli:cli``, the `cooler` program), so that hicstuff relies on
    cooler's command names and options, not on the modules defining them.
    """
    from cooler.cli import cli

    command = cli
    for name in names:
        command = command.commands[name]
    return command


def _run_cooler(command, args, what):
    """
    Run a `cooler` command in this process, as `cooler` would run it from the
    command line, and raise if it fails.

    Running cooler's own command, rather than its Python functions, keeps the
    command-line semantics (defaults, checks, and the free-form
    `balancing_args`) without needing a `cooler` program on PATH. cooler's
    progress records (one per chunk or balancing iteration) are kept out of
    hicstuff's log; its warnings are not.

    Parameters
    ----------
    command : click.Command
        The cooler command, e.g. ``_cooler_command("balance")``.
    args : list of str
        Its command-line arguments.
    what : str
        What is being run, for the error message.
    """
    cooler_logger = logging.getLogger("cooler")
    level = cooler_logger.level
    cooler_logger.setLevel(logging.WARNING)
    status = None
    try:
        status = command.main(args=list(args), prog_name="cooler", standalone_mode=False)
    except SystemExit as err:
        # Some commands exit rather than raise, e.g. on non-convergence
        status = err.code
    except click.exceptions.Abort as err:
        # click turns an interruption, but also an EOFError (e.g. from a
        # truncated gzip input), into Abort
        if isinstance(err.__cause__, KeyboardInterrupt):
            raise err.__cause__ from None
        # No local name for the cause: it would keep cooler's frames, and
        # with cooler < 0.11 its temporary files, alive until garbage collection
        raise RuntimeError(f"{what} failed: {err.__cause__ or err!r}") from (err.__cause__ or err)
    except Exception as err:
        raise RuntimeError(f"{what} failed: {err}") from err
    finally:
        cooler_logger.setLevel(level)
    # Without standalone mode, click returns the status of a `ctx.exit(status)`
    # rather than exiting
    if isinstance(status, (int, str)) and status not in (0, ""):
        raise RuntimeError(f"{what} failed (exit status {status}).")


def _run_cooler_program(names, args, what):
    """
    Run a `cooler` command in a new Python interpreter (`python -m cooler`),
    and raise if it fails, as :func:`_run_cooler` does in this process.

    Its standard output is hicstuff's; its log is kept and the last line of
    it given in the error.

    Parameters
    ----------
    names : list of str
        The names of the cooler command, e.g. ``["balance"]``.
    args : list of str
        Its command-line arguments.
    what : str
        What is being run, for the error message.
    """
    proc = sp.run(
        [sys.executable, "-m", "cooler", *names, *args],
        stderr=sp.PIPE,
        text=True,
    )
    if proc.returncode != 0:
        lines = [line for line in proc.stderr.splitlines() if line.strip()]
        reason = f": {lines[-1].strip()}" if lines else ""
        raise RuntimeError(f"{what} failed (exit status {proc.returncode}){reason}")


def _run_cooler_cload_pairs(bins, pairs_file, cool_file):
    """
    Run `cooler cload pairs` on a hicstuff pairs file and fail loudly if
    cooler does not succeed.

    Parameters
    ----------
    bins : str
        Bins argument passed to cooler: either a BED-like bins file or a
        ``chromsizes:binsize`` string.
    pairs_file : str
        Path to the pairs file containing input contact data.
    cool_file : str
        Path to the output cool file name to generate.
    """
    # Field numbers are one-based and follow the hicstuff pairs columns:
    # readID chr1 pos1 chr2 pos2 strand1 strand2 [frag1 frag2]
    args = ["--chrom1", "2", "--pos1", "3", "--chrom2", "4", "--pos2", "5"]
    _run_cooler(
        _cooler_command("cload", "pairs"),
        args + [bins, pairs_file, cool_file],
        f"`cooler cload pairs` generating {cool_file} from {pairs_file}",
    )


def pairs2cool(pairs_file, cool_file, bins_file, exclude):
    """
    Make a cooler file from the pairs file. See: https://github.com/mirnylab/cooler/ for more informations.

    Parameters
    ----------

    pairs_file : str
        Path to the pairs file containing input contact data.
    cool_file : str
        Path to the output cool file name to generate.
    bins_file : str
        Path to the file containing genomic segmentation information. (fragments_list.txt).
    """

    # Exclude some chromosomes from bins
    bins_tmp = bins_file + ".cooler"
    bins = pd.read_csv(bins_file, sep="\t", usecols=[1, 2, 3], skiprows=1, header=None)
    if exclude is not None:
        bins = bins[~bins[1].isin(exclude.split(","))]
    bins.to_csv(bins_tmp, sep="\t", header=False, index=False)

    # Make cool
    try:
        _run_cooler_cload_pairs(bins_tmp, pairs_file, cool_file)
    finally:
        os.remove(bins_tmp)


def pairs2binnedcool(pairs_file, cool_file, binning, info_contigs, exclude):
    """
    Make a *binned* cooler file from the pairs file. See: https://github.com/mirnylab/cooler/ for more informations.

    Parameters
    ----------

    pairs_file : str
        Path to the pairs file containing input contact data.
    cool_file : str
        Path to the output cool file name to generate.
    binning : int
        If mat_fmt is set to "cool", the cool file will be further binned to
        this resolution.
    info_contigs : pathlib.Path or str
        The path to the contigs info file
    """

    # Exclude some chromosomes from bins
    chroms_tmp = info_contigs + ".chroms"
    chroms = pd.read_csv(info_contigs, sep="\t", usecols=[0, 1], skiprows=1, header=None)
    if exclude is not None:
        chroms = chroms[~chroms[0].isin(exclude.split(","))]
    chroms.to_csv(chroms_tmp, index=False, sep="\t", header=False)

    # Run `cooler cload pairs`
    try:
        _run_cooler_cload_pairs(chroms_tmp + ":" + str(binning), pairs_file, cool_file)
    finally:
        os.remove(chroms_tmp)


def cool2mcool(cool_file, output_file):
    """
    Zoomify a *binned* cooler file. See: https://github.com/mirnylab/cooler/ for more informations.

    Parameters
    ----------

    cool_file : str
        Path to an existing binned .cool file.
    output_file : str
        Path to the new multi-resolution .mcool file.
    """

    # Get resolutions
    clr = cooler.Cooler(cool_file)
    res = clr.binsize
    multires = [res, res * 2, res * 5, res * 10, res * 20, res * 50, res * 100]

    # Run cooler zoomify
    cooler.zoomify_cooler(clr.filename, output_file, multires, chunksize=10000000)


def _balancing_args(balancing_args, threads=1):
    """
    Split `balancing_args` into `cooler balance` arguments, the way a shell
    would, and check them against cooler's own option definitions, so that
    a typo fails before the pipeline runs rather than after.

    Unless `balancing_args` sets it, the number of processes (--nproc) is
    `threads` rather than cooler's default of 8.

    Parameters
    ----------
    balancing_args : str or None
        Extra arguments passed to `cooler balance`.
    threads : int
        Number of processes used to balance, unless set in `balancing_args`.

    Returns
    -------
    list of str
    """
    invalid = f"Invalid balancing arguments {balancing_args!r}"
    try:
        args = shlex.split(balancing_args) if balancing_args else []
    except ValueError as err:  # e.g. an unclosed quote
        raise ValueError(f"{invalid}: {err}") from None
    command = _cooler_command("balance")
    try:
        ctx = command.make_context("cooler balance", args + ["COOL_PATH"])
    except click.ClickException as err:
        raise ValueError(f"{invalid}: {err.format_message()}") from None
    except click.exceptions.Exit:
        # --help: click printed cooler's help and asks to exit
        raise ValueError(f"{invalid}: they ask for help, printed above") from None
    if ctx.get_parameter_source("nproc") is click.core.ParameterSource.DEFAULT:
        args = ["--nproc", str(threads)] + args
    return args


def balance(cool_file, balancing_args, threads=1):
    """
    Balance a *binned* cooler file, every resolution of a .mcool file, with
    `cooler balance`. See: https://github.com/mirnylab/cooler/ for more informations.

    Balancing never stops the pipeline. When `cooler balance` fails at a
    resolution, a warning carrying cooler's error is logged and that
    resolution stays as it was; when the weights it stores did not
    converge, a warning says so. "Matrix balanced" is logged only for the
    resolutions that got weights.

    Parameters
    ----------

    cool_file : str
        Path to an existing binned .(m)cool file.
    balancing_args : str
        Extra arguments passed to `cooler balance` (default: None)
    threads : int
        Number of processes used to balance, unless set in `balancing_args`.

    Returns
    -------
    list of str
        The URIs of the resolutions that hold balancing weights.
    """
    args = _balancing_args(balancing_args, threads)
    command = _cooler_command("balance")
    params = command.make_context("cooler balance", args + ["COOL_PATH"]).params
    # The weight column cooler writes: --name, "weight" by default
    name = params["name"]
    # cooler's process pool starts its workers with Python's default start
    # method. Unless it is fork, every worker imports the script that called
    # hicstuff again (spawn on macOS, forkserver on Linux from Python 3.14),
    # which then fails unless that script guards its code with
    # `if __name__ == "__main__":`. Balance in a new `python -m cooler`
    # process then, whose main module is guarded.
    new_process = params["nproc"] > 1 and multiprocessing.get_start_method() != "fork"
    if cooler.fileops.is_multires_file(cool_file):
        uris = [cool_file + "::" + path for path in cooler.fileops.list_coolers(cool_file)]
    else:
        uris = [cool_file]
    balanced = []
    for uri in uris:
        try:
            what = f"`cooler balance` of {uri}"
            if new_process:
                _run_cooler_program(["balance"], args + [uri], what)
            else:
                _run_cooler(command, args + [uri], what)
        except RuntimeError as err:
            logger.warning("Continuing without balancing: %s", err)
            continue
        with cooler.Cooler(uri).open("r") as grp:
            weights = grp["bins"].get(name)
            # One flag per chromosome with --cis-only or --trans-only
            converged = weights is None or bool(np.all(weights.attrs.get("converged", True)))
        if weights is None:
            # e.g. --stdout, or --convergence-policy discard without convergence
            logger.warning("`cooler balance` stored no %r column in %s.", name, uri)
            continue
        if not converged:
            logger.warning(
                "`cooler balance` did not converge for %s: its %r column is "
                "stored with converged=False.",
                uri,
                name,
            )
        balanced.append(uri)
    if balanced:
        logger.info(
            "Matrix balanced with cooler balance%s.",
            "" if len(balanced) == len(uris) else f" ({len(balanced)} of {len(uris)} resolutions)",
        )
    return balanced


def pairs2matrix(pairs_file, mat_file, fragments_file, mat_fmt="graal", threads=1, tmp_dir=None):
    """Generate the matrix by counting the number of occurences of each
    combination of restriction fragments in a pairs file.

    Parameters
    ----------
    pairs_file : str
        Path to a Hi-C pairs file, with frag1 and frag2 columns added.
    mat_file : str
        Path where the matrix will be written.
    fragments_file : str
        Path to the fragments_list.txt file. Used to know total
        matrix size in case some observations are not observed at the end.
    mat_fmt : str
        The format to use when writing the matrix. Can be graal or bg2 format.
    threads : int
        Number of threads to use in parallel.
    tmp_dir : str
        Temporary directory for sorting files. If None given, will use the system default.
    """
    # Number of fragments is N lines in frag list - 1 for the header
    n_frags = sum(1 for line in open(fragments_file)) - 1
    frags = pd.read_csv(fragments_file, delimiter="\t")

    def write_mat_entry(frag1, frag2, contacts):
        """Write a single sparse matrix entry in either graal or bg2 format"""
        if mat_fmt == "graal":
            mat.write("\t".join(map(str, [frag1, frag2, contacts])) + "\n")
        elif mat_fmt == "bg2":
            frag1, frag2 = int(frag1), int(frag2)
            mat.write(
                "\t".join(
                    map(
                        str,
                        [
                            frags.chrom[frag1],
                            frags.start_pos[frag1],
                            frags.end_pos[frag1],
                            frags.chrom[frag2],
                            frags.start_pos[frag2],
                            frags.end_pos[frag2],
                            contacts,
                        ],
                    )
                )
                + "\n"
            )

    pre_mat_file = mat_file + ".pre.pairs"
    hio.sort_pairs(
        pairs_file,
        pre_mat_file,
        keys=["frag1", "frag2"],
        threads=threads,
        tmp_dir=tmp_dir,
    )
    header_size = len(hio.get_pairs_header(pre_mat_file))
    with open(pre_mat_file) as pairs, open(mat_file, "w") as mat:
        # Skip header lines
        for _ in range(header_size):
            next(pairs)
        prev_pair = None  # Pairs identified by [frag1, frag2]
        n_occ = 0  # Number of occurences of each frag combination
        n_nonzero = 0  # Total number of nonzero matrix entries
        n_pairs = 0  # Total number of pairs entered in the matrix
        pairs_reader = csv.reader(pairs, delimiter="\t")
        # First line contains nrows, ncols and number of nonzero entries.
        # Number of nonzero entries is unknown for now
        if mat_fmt == "graal":
            mat.write("\t".join(map(str, [n_frags, n_frags, "-"])) + "\n")
        for pair in pairs_reader:
            # Fragment ids are field 8 and 9
            curr_pair = [pair[7], pair[8]]
            # Increment number of occurences for fragment pair
            if prev_pair == curr_pair:
                n_occ += 1
            # Write previous pair and start a new one
            else:
                if n_occ > 0:
                    write_mat_entry(prev_pair[0], prev_pair[1], n_occ)
                    n_nonzero += 1
                    n_pairs += n_occ
                prev_pair = curr_pair
                n_occ = 1
        # Write the last value
        if n_occ > 0:
            write_mat_entry(prev_pair[0], prev_pair[1], n_occ)
            n_nonzero += 1
            n_pairs += n_occ

    # Edit header line to fill number of nonzero entries inplace in graal header
    if mat_fmt == "graal":
        with open(mat_file) as mat, open(pre_mat_file, "w") as tmp_mat:
            header = mat.readline()
            header = header.replace("-", str(n_nonzero))
            tmp_mat.write(header)
            st.copyfileobj(mat, tmp_mat)
        # Replace the matrix file with the one with corrected header
        os.rename(pre_mat_file, mat_file)
    else:
        os.remove(pre_mat_file)

    logger.info(
        "%d pairs used to build a contact map of %d bins with %d nonzero entries.",
        n_pairs,
        n_frags,
        n_nonzero,
    )


def check_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    return which(name) is not None


#: Characters that make an input of ``hicstuff pipeline`` a glob pattern
_GLOB_CHARS = re.compile(r"[*?[]")


def expand_inputs(arg):
    """The files that an input argument of ``hicstuff pipeline`` names.

    An argument that is a path that exists is that file, whatever characters
    its name holds. Otherwise, it is a comma-separated list of paths, each of
    which may be a glob pattern (with ``*``, ``?`` or ``[...]``), replaced by
    the files it matches, sorted by name.

    Parameters
    ----------
    arg : str
        A path, comma-separated paths, or glob patterns.

    Returns
    -------
    list of str
        The paths of the files, in the order given.

    Raises
    ------
    ValueError
        If a pattern matches no file, or a list holds an empty item.

    Examples
    --------
    >>> expand_inputs("lane1_R1.fq.gz,lane2_R1.fq.gz")
    ['lane1_R1.fq.gz', 'lane2_R1.fq.gz']
    """
    if os.path.exists(arg):
        return [arg]
    files = []
    for item in arg.split(","):
        if not item:
            raise ValueError(f"Empty path in the input list {arg!r}.")
        if _GLOB_CHARS.search(item) is None:
            files.append(item)
            continue
        matches = sorted(glob.glob(item))
        if not matches:
            raise ValueError(f"No file matches the input pattern {item!r}.")
        files.extend(matches)
    return files


def split_inputs(args, start_stage="fastq"):
    """The R1 and R2 files given by the input arguments of ``hicstuff pipeline``.

    Each argument is expanded by :func:`expand_inputs`. From fastq or BAM
    files, two arguments are the R1 files and the R2 files. More arguments,
    as the shell gives them for unquoted glob patterns such as
    ``lib1_*_R1.fq.gz lib1_*_R2.fq.gz``, must be even in number: the first
    half are the R1 files, the second half the R2 files. From pairs, every
    argument is a pairs file.

    Parameters
    ----------
    args : list of str
        The input arguments, as given on the command line.
    start_stage : str
        The start stage of the pipeline: "fastq", "bam", "pairs" or
        "pairs_idx".

    Returns
    -------
    tuple of (list of str, list of str)
        The R1 files and the R2 files, in matching order; the files and an
        empty list from pairs, or from a single argument.

    Raises
    ------
    ValueError
        If the R1 and R2 files cannot be told apart, or an argument names no
        file (see :func:`expand_inputs`).

    Examples
    --------
    >>> split_inputs(["a_R1.fq.gz,b_R1.fq.gz", "a_R2.fq.gz,b_R2.fq.gz"])
    (['a_R1.fq.gz', 'b_R1.fq.gz'], ['a_R2.fq.gz', 'b_R2.fq.gz'])
    >>> split_inputs(["a_R1.fq.gz", "b_R1.fq.gz", "a_R2.fq.gz", "b_R2.fq.gz"])
    (['a_R1.fq.gz', 'b_R1.fq.gz'], ['a_R2.fq.gz', 'b_R2.fq.gz'])
    """
    args = list(args)
    if start_stage in ("pairs", "pairs_idx"):
        return [f for arg in args for f in expand_inputs(arg)], []
    if len(args) > 2 and len(args) % 2:
        raise ValueError(
            f"{len(args)} input files: give the R1 files, then the R2 files, as many of each."
        )
    half = max(len(args) // 2, 1)
    forward = [f for arg in args[:half] for f in expand_inputs(arg)]
    reverse = [f for arg in args[half:] for f in expand_inputs(arg)]
    return forward, reverse


def _input_files(value):
    """The list of files an input of full_pipeline gives: None, a path, or a
    list of paths."""
    if value is None:
        return []
    if isinstance(value, (str, os.PathLike)):
        return [os.fspath(value)]
    return [os.fspath(path) for path in value]


def _is_gzip(path):
    """Whether a file starts as a gzip file."""
    with open(path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def _first_read_name(path):
    """The name of the first read of a fastq file, without the /1 or /2 that
    some files end it with, or None if the file holds no read."""
    with hio.read_compressed(path) as handle:
        header = handle.readline()
    if not header.startswith("@"):
        return None
    fields = header[1:].split(maxsplit=1)
    if not fields:
        return None
    return re.sub(r"/[12]$", "", fields[0])


def check_fastq_pairs(forward, reverse):
    """Check that several pairs of fastq files can be read as one library.

    Every file must exist and be given once, and the R1 and R2 files of each
    pair must start with the same read: Illumina names both mates of a read
    pair alike, with an optional /1 or /2.

    Parameters
    ----------
    forward : list of str
        The R1 files.
    reverse : list of str
        The R2 files, in the order of their R1 files.

    Raises
    ------
    ValueError
        If a file is missing or given twice, or the files of a pair do not
        start with the same read.
    """
    files = forward + reverse
    missing = [path for path in files if not os.path.exists(path)]
    if missing:
        raise ValueError(f"Input file not found: {', '.join(missing)}")
    seen = {}
    for path in files:
        real = os.path.realpath(path)
        if real in seen:
            raise ValueError(f"Input file given twice: {seen[real]} and {path}.")
        seen[real] = path
    for fq1, fq2 in zip(forward, reverse):
        name1, name2 = _first_read_name(fq1), _first_read_name(fq2)
        if name1 != name2:
            raise ValueError(
                f"{fq1} and {fq2} are not the R1 and R2 files of the same reads: they start "
                f"with reads {name1} and {name2}. Give the R1 files, then the R2 files, in "
                "the same order."
            )


def concatenate_fastq(files, out_file):
    """Concatenate fastq files into one, gzip-compressed if any of them is.

    The gzip files are copied as they are, which makes a valid gzip file of
    several members, read as one by Python, htslib and the aligners. With a
    gzip file among them, the plain files are compressed; otherwise the
    output is a plain file. A plain file that does not end with a newline
    gets one.

    Parameters
    ----------
    files : list of str
        The fastq files, plain or gzip-compressed, in the order to write them.
    out_file : str
        Path of the concatenated file.
    """
    gzipped = [_is_gzip(path) for path in files]
    with open(out_file, "wb") as out:
        for path, is_gzip in zip(files, gzipped):
            with open(path, "rb") as src:
                if is_gzip:
                    st.copyfileobj(src, out, 1 << 20)
                elif any(gzipped):
                    # A temporary file, read again soon: compressed fast
                    with gzip.GzipFile(fileobj=out, mode="wb", compresslevel=1) as dst:
                        _copy_lines(src, dst)
                else:
                    _copy_lines(src, out)


def _copy_lines(src, dst):
    """Copy a text file, adding the newline it may lack at its end."""
    last = b""
    for block in iter(lambda: src.read(1 << 20), b""):
        dst.write(block)
        last = block[-1:]
    if last not in (b"", b"\n"):
        dst.write(b"\n")


def full_pipeline(
    genome,
    input1,
    input2=None,
    aligner=None,
    centromeres=None,
    exclude=None,
    circular=False,
    distance_law=False,
    enzyme=5000,
    filter_events=False,
    force=False,
    mapping="normal",
    mat_fmt="cool",
    binning=0,
    zoomify=True,
    balancing_args=None,
    min_qual=30,
    min_size=0,
    no_cleanup=False,
    out_dir=None,
    pcr_duplicates=False,
    plot=False,
    prefix=None,
    read_len=None,
    remove_centros=None,
    start_stage="fastq",
    threads=1,
    tmp_dir=None,
    skip_count=False,
    aligner_backend=None,
    cmd=None,
):
    """
    Run the whole hicstuff pipeline. Starting from fastq files and a genome to
    obtain a contact matrix.

    Parameters
    ----------
    genome : str
        Path to the bowtie2/bwa index prefix if using bowtie2/bwa or to the genome
        in fasta format if using minimap2.
    input1 : str or list of str
        Path to the Hi-C reads in fastq format (forward), the aligned Hi-C reads
        in BAM format, or the pairs file, depending on the value of start_stage.
        From fastq, a list of files: the R1 files of several pairs of fastq
        files (lanes or runs of a library), read as one library.
    input2 : str or list of str
        Path to the Hi-C reads in fastq format (reverse), the aligned Hi-C reads
        in BAM format, or None, depending on the value of start_stage. From
        fastq, a list of files: the R2 files, in the order of their R1 files.
        Each pair must start with the same read; the R1 files are concatenated
        in the temporary directory, as are the R2 files, and aligned as one
        pair of files.
    binning : int
        If mat_fmt is set to "cool", the cool file will be further binned to
        this resolution.
    zoomify : bool
        Whether to zoomify binned cool matrix (only used if mat_fmt == "cool" and binning is set)
    balancing_args : str
        Extra arguments passed to `cooler balance` (default: None) (only used if
        mat_fmt == "cool" and binning is set)
    enzyme : int or strtest_data/genome/seq.fa
    circular : bool
        Use if the genome is circular.
    out_dir : str or None
        Path where output files should be written. Current directory by default.
    tmp_dir : str or None
        Path where temporary files will be written. Creates a "tmp" folder in
        out_dir by default.
    plot : bool
        Whether plots should be generated at different steps of the pipeline.
        Plots are saved in a "plots" directory inside out_dir.
    min_qual : int
        Minimum mapping quality required to keep a pair of Hi-C reads.
    min_size : int
        Minimum contig size required to keep it. Shorter contigs are left out
        of the fragments list, and pairs with a mate on one of them are
        discarded.
    threads : int
        Number of threads to use for parallel operations.
    no_cleanup : bool
        Whether temporary files should be deleted at the end of the pipeline.
    mapping : str
        normal|iterative|cutsite|5spm. Use normal, iterative, cutsite or 5spm
        mapping. "normal": Normal alignement. "iterative": Truncates and extends
        reads until unambiguous alignment. "cutsite": Digest reads at religation
        sites and build new pairs from the fragments created. "5spm": Align the
        mates with ``bwa mem -5SPM``, which makes the 5'-most alignment of a read
        its primary one, as the 4DN, Juicer and Arima pipelines do; it aligns
        with bwa, in-process with pybwa unless aligner_backend says otherwise.
        The pairs of mode 5spm record each mate at its aligned 5' end, those of
        the other modes at the leftmost base of its alignment.
    filter_events : bool
        Filter spurious or uninformative 3C events. Requires a restriction enzyme.
    force : bool
        If True, overwrite existing files with the same name as output.
    prefix : str or None
        Choose a common name for output files instead of default graal names.
    start_stage : str
        Step at which the pipeline should start. Can be "fastq", "bam", "pairs"
        or "pairs_idx". With starting from bam allows to skip alignment and start
        from named-sorted bam files. With
        "pairs", a single pairs file is given as input, and with "pairs_idx", the
        pairs in the input must already be attributed to fragments and fragment
        attribution is skipped.
    mat_fmt : str
        Select the output matrix format. Can be either "bg2" for the
        bedgraph2 format, "cool" for Mirnylab's cool format, graal for a
        plain text COO format compatible with Koszullab's instagraal software,
        or "none" for no matrix: the run ends with its pairs file, which tools
        that build their own contact maps read (metaTOR, metacooler). Options
        of the cool matrix, binning other than 0, and balancing_args and
        exclude that are not blank, are refused with it.
    aligner : str
        Read alignment software to use. Can be either "minimap2", "bwa" or
        "bowtie2". Defaults to bowtie2, and to bwa in mapping mode 5spm, which
        only aligns with bwa.
    pcr_duplicates : bool
        If True, PCR duplicates will be filtered based on genomic positions.
        Pairs where both reads have exactly the same coordinates are considered
        duplicates and only one of those will be conserved.
    distance_law : bool
        If True, generates a distance law file with the values of the probabilities
        to have a contact between two distances for each chromosomes or arms if the
        file with the positions has been given. The values are not normalized, or
        averaged.
    centromeres : None or str
        If not None, path of file with Positions of the centromeres separated by a
        space and in the same order than the chromosomes.
    exclude : None or str
        If not None, the name of the chromosomes to remove (e.g. "2u,chrM")
    read_len : int
        Maximum read length to expect in the fastq file. Optionally used in iterative
        alignment mode. Estimated from the first read by default. Useful if input fastq
        is a composite of different read lengths.
    remove_centros : None or int
        If the distance law is computed, this is the number of kb that will be removed
        around the centromere position given by in the centromere file.
    skip_count : bool
        If True, skip counting reads in the input FASTQ files. Useful for very
        large files where counting is slow.
    aligner_backend : str
        How bwa runs: "binary" (the bwa program), "python" (in-process with
        pybwa) or "auto". None reads the environment variable
        HICSTUFF_ALIGNER_BACKEND, and is auto if it is not set. auto runs the
        bwa program if it is on PATH and pybwa otherwise, and in mapping mode
        5spm pybwa first.
    cmd : str
        The command line that ran the pipeline, written to the log header after
        the date. None leaves it out.
    """
    if mapping == "5spm":
        if aligner not in (None, "bwa"):
            msg = (
                f"Mapping mode 5spm aligns with bwa mem -5SPM, not with {aligner}: "
                "leave out the aligner, or set it to bwa."
            )
            logger.error(msg)
            raise ValueError(msg)
        aligner = "bwa"
    elif aligner is None:
        aligner = "bowtie2"
    if aligner_backend is not None and aligner != "bwa" and aligner_backend.lower() == "python":
        msg = f"pybwa runs bwa only: --aligner-backend python cannot align with {aligner}."
        logger.error(msg)
        raise ValueError(msg)
    if aligner not in ("bowtie2", "minimap2", "bwa"):
        logger.error("Incompatible aligner software, choose bowtie2, minimap2 or bwa.")
        raise ValueError("aligner should be either bowtie2, minimap2 or bwa.")
    if mat_fmt not in MATRIX_FORMATS:
        msg = f"Unknown matrix format {mat_fmt!r}: choose graal, bg2, cool or none."
        logger.error(msg)
        raise ValueError(msg)
    if mat_fmt == "cool":
        try:
            import cooler
        except ImportError:
            logger.error(
                "The cooler package is require to return matrix in cool format, please install it first."
            )
            raise ImportError("The cooler package is required.") from None
        if binning > 0:
            # Fail now on invalid balancing arguments, not after the alignment
            _balancing_args(balancing_args)
    elif mat_fmt == "none":
        cool_only = [
            name
            for name, given in [
                ("binning", binning != 0),
                ("balancing_args", bool(balancing_args and balancing_args.strip())),
                ("exclude", any(name.strip() for name in _excluded(exclude))),
            ]
            if given
        ]
        if cool_only:
            options = "is an option" if len(cool_only) == 1 else "are options"
            msg = (
                f"The matrix format none writes no matrix: {', '.join(cool_only)} "
                f"{options} of the cool matrix."
            )
            logger.error(msg)
            raise ValueError(msg)

    # Pipeline can start from 3 input types
    start_time = datetime.now()
    stages = {"fastq": 0, "bam": 1, "pairs": 2, "pairs_idx": 3}
    stage_name, start_stage = start_stage, stages[start_stage]

    # Check if the genome is an index or fasta file
    idx = hio.check_fasta_index(genome, mode=aligner)
    is_fasta = hio.check_is_fasta(genome)
    # The genome is a fasta file, or the prefix of a bowtie2 or bwa index from
    # which its sequence is read, for the digestion and the fragments.
    # minimap2 reads a fasta file only.
    if not (is_fasta or (idx is not None and aligner != "minimap2")):
        logger.error(
            "The genome must be a fasta file, or the prefix of a %s index: %s",
            "bowtie2 or bwa" if aligner == "minimap2" else aligner,
            genome,
        )
        sys.exit(1)
    # Check that the third parties the run uses can be run: the aligner only
    # aligns reads from fastq files. bwa runs as a program or in-process.
    backend = None
    if start_stage == 0:
        if aligner == "bwa":
            backend = hbm.resolve_backend(aligner_backend, prefer_python=mapping == "5spm")
        if backend != "python":
            hio.require_tools(
                hio.alignment_tools(aligner, idx),
                f"to align the reads of start stage {stage_name}",
                alternative=(
                    "run with -a bwa, which aligns in-process with pybwa"
                    if aligner != "bwa" and hbm.pybwa_version() is not None
                    else None
                ),
            )
    if aligner == "bowtie2" and not is_fasta:
        hio.require_tools(
            ["bowtie2-inspect"], f"to extract the genome sequence from the bowtie2 index {genome}"
        )

    # Check if the number of input files is correct
    for_files, rev_files = _input_files(input1), _input_files(input2)
    if start_stage <= 1:
        if not for_files or not rev_files:
            logger.error("You must provide 2 input files when --start-stage is fastq or bam.")
            sys.exit(1)
        if len(for_files) != len(rev_files):
            logger.error(
                "%d R1 files and %d R2 files: each R1 file needs its R2 file.",
                len(for_files),
                len(rev_files),
            )
            sys.exit(1)
        if start_stage == 1 and len(for_files) > 1:
            logger.error(
                "Several pairs of input files are read from fastq only: give one BAM file per mate."
            )
            sys.exit(1)
    else:
        if rev_files or len(for_files) != 1:
            logger.error(
                "You must provide a single input file when --start-stage is pairs or pairs_idx."
            )
            sys.exit(1)
    # Several pairs of fastq files are concatenated into one pair of files
    several_pairs = len(for_files) > 1
    if several_pairs:
        try:
            check_fastq_pairs(for_files, rev_files)
        except ValueError as error:
            logger.error(str(error))
            sys.exit(1)
    input1 = ",".join(for_files)
    input2 = ",".join(rev_files) if rev_files else None
    # sanitize enzyme
    enzyme = str(enzyme)
    # Remember whether fragments_file has been generated during this run
    fragments_updated = False

    if out_dir is None:
        out_dir = os.getcwd()

    if tmp_dir is None:
        tmp_dir = join(out_dir, "tmp")

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    # Define figures output paths
    if plot:
        fig_dir = join(out_dir, "plots")
        os.makedirs(fig_dir, exist_ok=True)
        if prefix:
            frag_plot = join(fig_dir, prefix + "_frags_hist.pdf")
            dist_plot = join(fig_dir, prefix + "_event_distance.pdf")
            pie_plot = join(fig_dir, prefix + "_event_distribution.pdf")
            distance_law_plot = join(fig_dir, prefix + "_distance_law.pdf")
        else:
            frag_plot = join(fig_dir, "frags_hist.pdf")
            dist_plot = join(fig_dir, "event_distance.pdf")
            pie_plot = join(fig_dir, "event_distribution.pdf")
            distance_law_plot = join(fig_dir, "distance_law.pdf")
        matplotlib.use("Agg")
    else:
        fig_dir = None
        dist_plot = pie_plot = frag_plot = None

    # Use current time for logging and to identify files
    now = time.strftime("%Y%m%d%H%M%S")

    def _tmp_file(fname):
        if prefix:
            fname = prefix + "." + fname
        full_path = join(tmp_dir, fname)
        if not force and os.path.exists(full_path):
            raise OSError(f"Temporary file {full_path} already exists. Use --force to overwrite")
        return full_path

    def _out_file(fname):
        if prefix:
            fname = prefix + "." + fname
        full_path = join(out_dir, fname)
        if not force and os.path.exists(full_path):
            raise OSError(f"Output file {full_path} already exists. Use --force to overwrite")

        return full_path

    # Define temporary file names
    log_file = _out_file("hicstuff_" + now + ".log")
    bam1 = _tmp_file("for.bam")
    bam2 = _tmp_file("rev.bam")
    pairs = _tmp_file("valid.pairs")
    pairs_idx = _tmp_file("valid_idx.pairs")
    pairs_filtered = _tmp_file("valid_idx_filtered.pairs")
    pairs_pcr = _tmp_file("valid_idx_pcrfree.pairs")
    reads_cat = []
    if several_pairs:
        reads_cat = [
            _tmp_file(f"{mate}.fq.gz" if any(map(_is_gzip, files)) else f"{mate}.fq")
            for mate, files in (("for", for_files), ("rev", rev_files))
        ]

    # Enable file logging
    hcl.set_file_handler(log_file)
    # The logs of bwa in-process, or in mapping mode 5spm, name the aligner;
    # other runs log what they did before these existed.
    generate_log_header(
        log_file,
        input1,
        input2,
        genome,
        enzyme,
        aligner=(
            hbm.describe(backend, five_spm=mapping == "5spm")
            if backend == "python" or mapping == "5spm"
            else None
        ),
        cmd=cmd,
    )

    # If defautl `mat_fmt` or set `mat_fmt=cool`, notify the user
    if mat_fmt == "cool":
        try:
            import cooler

            logger.info(
                "The default output format is now `.cool`. The Hi-C "
                "matrix will be generated with cooler v%s "
                "(Abdennur & Mirny, Bioinformatics 2020).",
                cooler.__version__,
            )
        except:
            raise

    # The genome has two roles: fasta is the file its sequences are read from
    # (digestion, fragments), index the prefix of the aligner's index. What
    # this run writes from the genome goes to a directory of its own, which
    # cleanup removes: the sequence read from an index without fasta, the
    # index built from a fasta file without one, or the files pybwa loads.
    # Nothing is written next to the genome, and a fasta file is not copied.
    genome_dir = None
    if (start_stage == 0 and (backend == "python" or idx is None)) or not is_fasta:
        genome_dir = hio.generate_temp_dir(tmp_dir)
    if is_fasta:
        fasta = genome
    elif aligner == "bwa":
        # A bwa index without fasta: its .pac file holds the sequence
        fasta = join(genome_dir, "genome.fa.gz")
        with gzip.open(fasta, "wb", compresslevel=6) as dst:
            hbm.write_fasta(genome, dst)
    else:
        # A bowtie2 index without fasta: bowtie2-inspect reads the sequence
        fasta = join(genome_dir, "genome.fa.gz")
        bt2fa = sp.Popen(
            ["bowtie2-inspect", genome],
            stdout=sp.PIPE,
            stderr=sp.PIPE,
        )
        # At gzip's default level, as `gzip -c` compressed it
        with gzip.open(fasta, "wb", compresslevel=6) as dst:
            st.copyfileobj(bt2fa.stdout, dst, 1 << 20)
        _, bt2err = bt2fa.communicate()
        # bowtie2-inspect still has return code 0 when crashing, need to
        # actively look for error in stderr
        if re.search(r"[Ee]rror", bt2err.decode()):
            logger.error(bt2err)
            logger.error(
                "bowtie2-inspect has failed, make sure you provided "
                "the path to the bowtie2 index without the extension."
            )
            sys.exit(1)

    # Check for spaces in fasta headers and issue error if found, and count
    # the contigs of the cool matrix (those --size and --exclude keep), before
    # an index is built: its warning comes first.
    excluded = _excluded(exclude)
    n_contigs = 0
    for record in SeqIO.parse(hio.read_compressed(fasta), "fasta"):
        if len(record.seq) >= int(min_size) and record.id not in excluded:
            n_contigs += 1
        if " " in record.id:
            logger.error("Sequence identifiers contain spaces. Please clean the input genome.")
    if mat_fmt == "cool" and n_contigs >= MANY_CONTIGS:
        logger.warning(
            f"{n_contigs:,} contigs in the cool matrix: it will be slow to write, as cooler "
            "takes a time that grows with the square of the number of contigs for every cooler "
            "it writes (about 30 s at 117,000 contigs, four times more for twice as many). "
            "--matfmt none ends the run with its pairs file, which metaTOR and metacooler read; "
            "--matfmt bg2 and graal write a text matrix."
        )

    # A fasta file without bowtie2 or bwa index is indexed in genome_dir when
    # reads are aligned, read where it is, plain or compressed. pybwa indexes
    # it when it loads it, below.
    index = idx
    if idx is None and start_stage == 0:
        if aligner == "bwa" and os.path.isfile(genome + ".bwt.2bit.64"):
            logger.warning(
                "Found a bwa-mem2 index at %s, which bwa and pybwa cannot read: "
                "a bwa index is built instead.",
                genome,
            )
        if backend != "python":
            index = join(genome_dir, "genome")
            logger.info(
                "%s index not found at %s, building a temporary one in %s.",
                aligner,
                genome,
                genome_dir,
            )
            if aligner == "bowtie2":
                sp.run(["bowtie2-build", "-q", fasta, index], stderr=sp.PIPE)
            else:
                sp.run(["bwa", "index", "-p", index, fasta], stderr=sp.PIPE)

    # Define output file names (tsv files)
    if prefix:
        fragments_list = _out_file("frags.tsv")
        info_contigs = _out_file("chr.tsv")
        if mat_fmt != "none":
            mat = _out_file("mat.tsv")
        # If matrix has a different format, give it the right extension
        if mat_fmt not in ("graal", "none"):
            mat = _out_file(mat_fmt)
    else:
        # Default graal file names
        fragments_list = _out_file("fragments_list.txt")
        info_contigs = _out_file("info_contigs.txt")
        if mat_fmt != "none":
            mat = _out_file("abs_fragments_contacts_weighted.txt")
        if mat_fmt not in ("graal", "none"):
            mat = _out_file("abs_fragments_contacts_weighted." + mat_fmt)
    # Define what input files are given
    if start_stage == 0:
        reads1, reads2 = for_files[0], rev_files[0]
    elif start_stage == 1:
        bam1, bam2 = input1, input2
    elif start_stage == 2:
        pairs = input1
    elif start_stage == 3:
        pairs_idx = input1

    # Perform genome alignment
    nreads_input1 = 0
    if start_stage == 0:
        if skip_count:
            logger.info("Skipping read count check (--skip-count).")
        else:
            # Check number of reads in both fastqs, of each pair of files
            logger.info("Checking content of fastq files.")
            nreads_input2 = 0
            counts_differ = False
            for fq1, fq2 in zip(for_files, rev_files):
                nreads1 = hio.check_fastq_entries(fq1, threads)
                nreads2 = hio.check_fastq_entries(fq2, threads)
                nreads_input1 += nreads1
                nreads_input2 += nreads2
                counts_differ = counts_differ or nreads1 != nreads2
                if several_pairs and nreads1 != nreads2:
                    logger.error(
                        "Fastq files %s and %s do not have the same number of reads: %d and %d.",
                        fq1,
                        fq2,
                        nreads1,
                        nreads2,
                    )
                elif several_pairs:
                    logger.info("%d reads in each of %s and %s.", nreads1, fq1, fq2)
            if counts_differ:
                logger.error("Fastq files do not have the same number of reads.")
            else:
                logger.info(f"{int(nreads_input1)} reads found in each fastq file.")

        if several_pairs:
            logger.info(
                "Concatenating %d pairs of fastq files into %s and %s.",
                len(for_files),
                *reads_cat,
            )
            for files, out in zip((for_files, rev_files), reads_cat):
                concatenate_fastq(files, out)
            reads1, reads2 = reads_cat

        # Define mapping choice (default normal):
        if mapping == "normal":
            iterative = False
        elif mapping == "iterative":
            iterative = True
        elif mapping == "cutsite":
            # If no enzyme given use iterative alignment.
            try:
                int(enzyme)
                logger.warning(
                    "No enzyme has been given. Can't map using cutsite, iterative mapping will be used instead."
                )
                iterative = True
            # If cutsite enabled and enzyme given, cut the reads before making a
            # normal alignment.
            except ValueError:
                iterative = False
                digest_for = _tmp_file("digest_for.fq.gz")
                digest_rev = _tmp_file("digest_rev.fq.gz")
                hcc.cut_ligation_sites(
                    fq_for=reads1,
                    fq_rev=reads2,
                    digest_for=digest_for,
                    digest_rev=digest_rev,
                    enzyme=enzyme,
                    mode="for_vs_rev",
                    seed_size=20,
                    n_cpu=threads,
                )
                reads1, reads2 = digest_for, digest_rev
        elif mapping != "5spm":
            logger.error("mapping must be either normal, iterative, cutsite or 5spm.")
            raise ValueError

        # bwa in-process loads the index once, for both mates, and builds it
        # from the fasta if there is none
        if backend == "python":
            logger.info("Aligning with %s.", hbm.describe(backend, five_spm=mapping == "5spm"))
            index = hbm.load_index(fasta if index is None else index, genome_dir)
        if mapping == "5spm":
            logger.info("Now mapping R1 and R2 reads with bwa mem -5SPM...")
            hbm.align_pairs(
                reads1, reads2, index, bam1, bam2, threads=threads, backend=backend, tmp_dir=tmp_dir
            )
        else:
            logger.info("Now separately mapping R1 and R2 reads...")
            for reads, bam in ((reads1, bam1), (reads2, bam2)):
                align_reads(
                    reads,
                    index,
                    bam,
                    tmp_dir=tmp_dir,
                    threads=threads,
                    aligner=aligner,
                    iterative=iterative,
                    min_qual=min_qual,
                    read_len=read_len,
                    backend=backend,
                )

    # Detect if multiple enzymes are given
    if re.search(",", enzyme):
        enzyme = enzyme.split(",")

    # Starting from bam files
    if start_stage <= 1:
        # Check number of reads in both fastqs
        if bam1 == input1:
            logger.info("Checking content of bam files.")
            nreads_input1 = hio.check_bam_entries(bam1)
            nreads_input2 = hio.check_bam_entries(bam2)
            if nreads_input1 != nreads_input2:
                logger.error("Bam files do not have the same number of reads.")
            else:
                logger.info(f"{nreads_input1} reads found in each bam file.")

        fragments_updated = True
        # Generate info_contigs and fragments_list output files
        hcd.write_frag_info(
            fasta,
            enzyme,
            min_size=min_size,
            circular=circular,
            output_contigs=info_contigs,
            output_frags=fragments_list,
        )

        # Log fragment size distribution
        hcd.frag_len(frags_file_name=fragments_list, plot=plot, fig_path=frag_plot)

        # Make pairs file (readID, chr1, chr2, pos1, pos2, strand1, strand2).
        # Mode 5spm records the mates' aligned 5' ends.
        bam2pairs(
            bam1,
            bam2,
            pairs,
            info_contigs,
            min_qual=min_qual,
            tmp_dir=tmp_dir,
            threads=threads,
            five_prime=mapping == "5spm",
        )

    # Pairs with a mate on a contig shorter than min_size, which has no
    # fragment in the fragments list
    small_contig_pairs = 0
    # Starting from pairs file
    if start_stage <= 2:
        restrict_table = {}
        for record in SeqIO.parse(hio.read_compressed(fasta), "fasta"):
            # Get chromosome restriction table
            restrict_table[record.id] = hcd.get_restriction_table(
                record.seq, enzyme, circular=circular
            )

        # Add fragment index to pairs (readID, chr1, pos1, chr2,
        # pos2, strand1, strand2, frag1, frag2), without the pairs on contigs
        # shorter than min_size
        small_contig_pairs = hcd.attribute_fragments(
            pairs, pairs_idx, restrict_table, min_size=min_size
        )

    # Sort pairs file by coordinates (and strands, for PCR duplicate removal)
    # for next steps
    hio.sort_pairs(
        pairs_idx,
        pairs_idx + ".sorted",
        keys=COORD_SORT_KEYS,
        threads=threads,
        tmp_dir=tmp_dir,
    )
    os.rename(pairs_idx + ".sorted", pairs_idx)

    # Count total pairs, including those discarded on small contigs
    tot_pairs = small_contig_pairs
    with open(pairs_idx) as file:
        for line in file:
            if line.startswith("#"):
                continue
            else:
                tot_pairs += 1
    if nreads_input1 != 0:
        logger.info(
            f"{tot_pairs} pairs successfully mapped ({round(100 * tot_pairs / (nreads_input1), 2)}%)"
        )
    else:
        logger.info(f"{tot_pairs} pairs successfully mapped")

    # Filter pairs if requested
    if filter_events:
        uncut_thr, loop_thr = hcf.get_thresholds(
            pairs_idx, plot_events=plot, fig_path=dist_plot, prefix=prefix
        )
        hcf.filter_events(
            pairs_idx,
            pairs_filtered,
            uncut_thr,
            loop_thr,
            plot_events=plot,
            fig_path=pie_plot,
            prefix=prefix,
        )
        use_pairs = pairs_filtered
    else:
        use_pairs = pairs_idx

    # Generate fragments file if it has not been already
    if not fragments_updated:
        hcd.write_frag_info(
            fasta,
            enzyme,
            min_size=min_size,
            circular=circular,
            output_contigs=info_contigs,
            output_frags=fragments_list,
        )

    # Generate distance law table if enabled
    if distance_law:
        out_distance_law = _out_file("distance_law.txt")
        if remove_centros is None:
            remove_centros = 0
        remove_centros = int(remove_centros)
        x_s, p_s, _ = hcdl.get_distance_law(
            pairs_idx,
            fragments_list,
            centro_file=centromeres,
            base=1.1,
            out_file=out_distance_law,
            circular=circular,
            rm_centro=remove_centros,
        )
        # Generate distance law figure is plots are enabled
        if plot:
            # Retrieve chrom labels from distance law file
            _, _, chr_labels = hcdl.import_distance_law(out_distance_law)
            chr_labels = [lab[0] for lab in chr_labels]
            chr_labels_idx = np.unique(chr_labels, return_index=True)[1]
            chr_labels = [chr_labels[index] for index in sorted(chr_labels_idx)]
            p_s = hcdl.normalize_distance_law(x_s, p_s)
            hcdl.plot_ps_slope(x_s, p_s, labels=chr_labels, fig_path=distance_law_plot)

    # Filter out PCR duplicates if requested
    if pcr_duplicates:
        filter_pcr_dup(use_pairs, pairs_pcr)
        use_pairs = pairs_pcr

    # Build matrix from pairs.
    if mat_fmt == "cool":
        # Log which pairs file is being used and how many pairs are listed
        pairs_count = 0
        with open(use_pairs) as file:
            for line in file:
                if line.startswith("#"):
                    continue
                else:
                    pairs_count += 1
        logger.info(
            "Generating matrix from pairs file %s (%d pairs in the file) ",
            use_pairs,
            pairs_count,
        )

        # Name matrix file in .cool
        mat = os.path.splitext(mat)[0] + ".cool"

        # If binning is **not** set, parse the pairs into a **un-binned** cool
        if binning == 0:
            ## THIS NEEDS TO BE FIXED AT SOME POINT
            pairs2cool(use_pairs, mat, fragments_list, exclude)

        # If binning is set, proceed to bin the pairs instead
        else:
            cool_file = os.path.splitext(mat)[0] + ".cool"
            pairs2binnedcool(use_pairs, cool_file, binning, info_contigs, exclude)
            mat = cool_file

            # If zoomify == True, zoomify binned cool
            if zoomify:
                mcool_file = os.path.splitext(mat)[0] + ".mcool"
                cool2mcool(mat, mcool_file)
                mat = mcool_file

            # Balance binned matrix
            balance(mat, balancing_args, threads=threads)
    elif mat_fmt == "none":
        logger.info("No matrix (--matfmt none): the run ends with its pairs file.")
    else:
        pairs2matrix(
            use_pairs,
            mat,
            fragments_list,
            mat_fmt=mat_fmt,
            threads=threads,
            tmp_dir=tmp_dir,
        )

    # Get stats on the pipeline
    try:
        logger.info("Fetching mapping and pairing stats")
        stats = hcs.get_pipeline_stats(log_file)
        logger.info(stats)
    except IndexError:
        logger.warning("IndexError. Stats not compiled.")
        pass

    # Move final pairs file to main dir.
    p = pathlib.Path(use_pairs).absolute()
    out_dir_abs = pathlib.Path(out_dir).absolute()
    pairsf = out_dir_abs / p.name
    # The final pairs are the input file itself when starting from indexed
    # pairs without filtering: it must be neither moved nor deleted.
    is_input = p == pathlib.Path(input1).absolute()
    # The file may already be in out_dir, even when named through another
    # path (tmp_dir == out_dir, a symbolic link to out_dir).
    in_out_dir = pairsf.exists() and os.path.samefile(p, pairsf)
    # Only rename (move) a temporary file; otherwise copy, unless the file
    # is already in out_dir.
    if not in_out_dir:
        if p.parent == out_dir_abs / "tmp" and not is_input:
            p.rename(pairsf)
        else:
            st.copy2(str(p), str(pairsf))

    # Sort and compress final pairs file
    sort_and_compress_pairs(str(pairsf), str(pairsf) + ".gz", tmp_dir=tmp_dir, threads=threads)
    if not (is_input and in_out_dir):
        os.remove(pairsf)

    # Clean temporary files
    if not no_cleanup:
        tempfiles = [
            pairs,
            pairs_idx,
            pairs_filtered,
            pairs_pcr,
            bam1,
            bam2,
            # The concatenated fastq files of several pairs of files
            *reads_cat,
        ]
        # The genome sequence and the index this run wrote, if any
        if genome_dir is not None:
            st.rmtree(genome_dir)
        # Do not delete files that were given as input
        try:
            tempfiles.remove(input1)
            tempfiles.remove(input2)
        except ValueError:
            pass

        # Delete single-resolution matrix if `--zoomify` is set
        if binning > 0 and zoomify:
            try:
                tempfiles.append(cool_file)
            except ValueError:
                pass

        # Remove the rest of tempfiles
        for file in tempfiles:
            try:
                os.remove(file)
            except FileNotFoundError:
                pass

    end_time = datetime.now()
    duration = relativedelta(end_time, start_time)
    made = "Pairs" if mat_fmt == "none" else "Contact map"
    logger.info(f"{made} generated after {duration.hours}h {duration.minutes}m {duration.seconds}s")
