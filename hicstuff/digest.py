#!/usr/bin/env python3
# coding: utf-8

"""Genome digestion

Functions used to write auxiliary instagraal compatible
sparse matrices.
"""

from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, Analysis
import os, sys, csv
import re
import collections
import copy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from hicstuff.log import logger
import hicstuff.io as hio

DEFAULT_FRAGMENTS_LIST_FILE_NAME = "fragments_list.txt"
DEFAULT_INFO_CONTIGS_FILE_NAME = "info_contigs.txt"
DEFAULT_SPARSE_MATRIX_FILE_NAME = "abs_fragments_contacts_weighted.txt"
DEFAULT_KB_BINNING = 1
DEFAULT_THRESHOLD_SIZE = 0
# Most used enzyme for eukaryotes
DEFAULT_ENZYME = "DpnII"
# If using evenly-sized chunks instead of restriction
# enzymes, they shouldn't be too short
DEFAULT_MIN_CHUNK_SIZE = 50


def write_frag_info(
    fasta,
    enzyme,
    min_size=DEFAULT_THRESHOLD_SIZE,
    circular=False,
    output_contigs=DEFAULT_INFO_CONTIGS_FILE_NAME,
    output_frags=DEFAULT_FRAGMENTS_LIST_FILE_NAME,
    output_dir=None,
):
    """Digest and write fragment information

    Write the fragments_list.txt and info_contigs.txt that are necessary for
    instagraal to run.

    Parameters
    ----------
    fasta : pathlib.Path or str
        The path to the reference genome
    enzyme : str, int or list of str
        If a string, must be the name of an enzyme (e.g. DpnII) and the genome
        will be cut at the enzyme's restriction sites. If a number, the genome
        will be cut uniformly into chunks with length equal to that number. A
        list of enzymes can also be specified if using multiple enzymes.
    min_size : float, optional
        Size below which shorter contigs are discarded. Default is 0, i.e. all
        contigs are retained.
    circular : bool, optional
        Whether the genome is circular. Default is False.
    output_contigs : str, optional
        The name of the file with contig info. Default is info_contigs.txt
    output_frags : str, optional
        The name of the file with fragment info. Default is fragments_list.txt
    output_dir : [type], optional
        The path to the output directory, which will be created if not already
        existing. Default is the current directory.
    """

    records = SeqIO.parse(fasta, "fasta")

    try:
        info_contigs_path = os.path.join(output_dir, output_contigs)
        frag_list_path = os.path.join(output_dir, output_frags)
    except TypeError:
        info_contigs_path = output_contigs
        frag_list_path = output_frags

    with open(info_contigs_path, "w") as info_contigs:

        info_contigs.write("contig\tlength\tn_frags\tcumul_length\n")

        with open(frag_list_path, "w") as fragments_list:

            fragments_list.write(
                "id\tchrom\tstart_pos" "\tend_pos\tsize\tgc_content\n"
            )

            total_frags = 0

            for record in records:
                contig_seq = record.seq
                contig_name = record.id
                contig_length = len(contig_seq)
                if contig_length < int(min_size):
                    continue

                sites = get_restriction_table(
                    contig_seq, enzyme, circular=circular
                )
                fragments = (
                    contig_seq[sites[i] : sites[i + 1]]
                    for i in range(len(sites) - 1)
                )
                n_frags = 0

                current_id = 1
                start_pos = 0
                for frag in fragments:
                    frag_length = len(frag)
                    if frag_length > 0:
                        end_pos = start_pos + frag_length
                        gc_content = SeqUtils.GC(frag) / 100.0

                        current_fragment_line = "%s\t%s\t%s\t%s\t%s\t%s\n" % (
                            current_id,
                            contig_name,
                            start_pos,
                            end_pos,
                            frag_length,
                            gc_content,
                        )

                        fragments_list.write(current_fragment_line)

                        try:
                            assert (current_id == 1 and start_pos == 0) or (
                                current_id > 1 and start_pos > 0
                            )
                        except AssertionError:
                            logger.error((current_id, start_pos))
                            raise
                        start_pos = end_pos
                        current_id += 1
                        n_frags += 1

                current_contig_line = "%s\t%s\t%s\t%s\n" % (
                    contig_name,
                    contig_length,
                    n_frags,
                    total_frags,
                )
                total_frags += n_frags
                info_contigs.write(current_contig_line)


def attribute_fragments(pairs_file, idx_pairs_file, restriction_table):
    """
    Writes the indexed pairs file, which has two more columns than the input
    pairs file corresponding to the restriction fragment index of each read.
    Note that pairs files have 1bp point positions whereas restriction table
    has 0bp point poisitions.

    Parameters
    ----------
    pairs_file: str
        Path the the input pairs file. Consists of 7 tab-separated
        columns: readID, chr1, pos1, chr2, pos2, strand1, strand2
    idx_pairs_file: str
        Path to the output indexed pairs file. Consists of 9 white space
        separated columns: readID, chr1, pos1, chr2, pos2, strand1, strand2,
        frag1, frag2. frag1 and frag2 are 0-based restriction fragments based
        on whole genome.
    restriction_table: dict
        Dictionary with chromosome identifiers (str) as keys and list of
        positions (int) of restriction sites as values.
    """

    # NOTE: Bottlenecks here are 1. binary search in find_frag and 2. writerow
    # 1. could be reduced by searching groups of N frags in parallel and 2. by
    # writing N frags simultaneously using a single call of writerows.

    # Parse and update header section
    pairs_header = hio.get_pairs_header(pairs_file)
    header_size = len(pairs_header)
    chrom_order = []
    with open(idx_pairs_file, "w") as idx_pairs:
        for line in pairs_header:
            # Add new column names to header
            if line.startswith("#columns"):
                line = line.rstrip() + " frag1 frag2"
            if line.startswith("#chromsize"):
                chrom_order.append(line.split()[1])
            idx_pairs.write(line + "\n")

    # Get number of fragments per chrom to allow genome-based indices
    shift_frags = {}
    prev_frags = 0
    for rank, chrom in enumerate(chrom_order):
        if rank > 0:
            # Note the "-1" because there are nfrags + 1 sites in rest table
            prev_frags += len(restriction_table[chrom_order[rank - 1]]) - 1
        # Idx of each chrom's frags will be shifted by n frags in previous chroms
        shift_frags[chrom] = prev_frags

    missing_contigs = set()
    # Attribute pairs to fragments and append them to output file (after header)
    with open(pairs_file, "r") as pairs, open(
        idx_pairs_file, "a"
    ) as idx_pairs:
        # Skip header lines
        for _ in range(header_size):
            next(pairs)

        # Define input and output fields
        pairs_cols = [
            "readID",
            "chr1",
            "pos1",
            "chr2",
            "pos2",
            "strand1",
            "strand2",
        ]
        idx_cols = pairs_cols + ["frag1", "frag2"]

        # Use csv reader / writer to automatically parse columns into a dict
        pairs_reader = csv.DictReader(
            pairs, fieldnames=pairs_cols, delimiter="\t"
        )
        pairs_writer = csv.DictWriter(
            idx_pairs, fieldnames=idx_cols, delimiter="\t"
        )

        for pair in pairs_reader:
            # Get the 0-based indices of corresponding restriction fragments
            # Deducing 1 from pair position to get it into 0bp point
            pair["frag1"] = find_frag(
                int(pair["pos1"]) - 1, restriction_table[pair["chr1"]]
            )
            pair["frag2"] = find_frag(
                int(pair["pos2"]) - 1, restriction_table[pair["chr2"]]
            )
            # Shift fragment indices to make them genome-based instead of
            # chromosome-based
            try:
                pair["frag1"] += shift_frags[pair["chr1"]]
            except KeyError:
                missing_contigs.add(pair["chr1"])
            try:
                pair["frag2"] += shift_frags[pair["chr2"]]
            except KeyError:
                missing_contigs.add(pair["chr2"])

            # Write indexed pairs in the new file
            pairs_writer.writerow(pair)

        if missing_contigs:
            logger.warning(
                "Pairs on the following contigs were discarded as "
                "those contigs are not listed in the paris file header. "
                "This is normal if you filtered out small contigs: %s"
                % " ".join(list(missing_contigs))
            )


def get_restriction_table(seq, enzyme, circular=False):
    """
    Get the restriction table for a single genomic sequence.

    Parameters
    ----------
    seq : Seq object
        A biopython Seq object representing a chromosomes or contig.
    enzyme : int, str or list of str
        The name of the restriction enzyme used, or a list of restriction
        enzyme names. Can also be an integer, to digest by fixed chunk size.
    circular : bool
        Wether the genome is circular.

    Returns
    -------
    numpy.array:
        List of restriction fragment boundary positions for the input sequence.

    >>> from Bio.Seq import Seq
    >>> get_restriction_table(Seq("AAGCCGGATCGG"),"HpaII")
    array([ 0,  4, 12])
    >>> get_restriction_table(Seq("AA"),["HpaII", "MluCI"])
    array([0, 2])
    >>> get_restriction_table(Seq("AA"),"aeiou1")
    Traceback (most recent call last):
        ...
    ValueError: aeiou1 is not a valid restriction enzyme.
    >>> get_restriction_table("AA","HpaII")
    Traceback (most recent call last):
        ...
    TypeError: Expected Seq or MutableSeq instance, got <class 'str'> instead

    """
    chrom_len = len(seq)
    wrong_enzyme = "{} is not a valid restriction enzyme.".format(enzyme)
    # Restriction batch containing the restriction enzyme
    try:
        enz = [enzyme] if isinstance(enzyme, str) else enzyme
        cutter = RestrictionBatch(enz)
    except (TypeError, ValueError):
        try:
            cutter = max(int(enzyme), DEFAULT_MIN_CHUNK_SIZE)
        except ValueError:
            raise ValueError(wrong_enzyme)

    # Conversion from string type to restriction type
    if isinstance(cutter, int):
        sites = [i for i in range(0, chrom_len, cutter)]
        if sites[-1] < chrom_len:
            sites.append(chrom_len)
    else:
        # Find sites of all restriction enzymes given
        ana = Analysis(cutter, seq, linear=not circular)
        sites = ana.full()
        # Gets all sites into a single flat list with 0-based index
        sites = [site - 1 for enz in sites.values() for site in enz]
        # Sort by position and allow first add start and end of seq
        sites.sort()
        sites.insert(0, 0)
        sites.append(chrom_len)

    return np.array(sites)


def find_frag(pos, r_sites):
    """
    Use binary search to find the index of a chromosome restriction fragment
    corresponding to an input genomic position.
    Parameters
    ----------
    pos : int
        Genomic position, in base pairs.
    r_sites : list
        List of genomic positions corresponding to restriction sites.
    Returns
    -------
    int
        The 0-based index of the restriction fragment to which the position belongs.

    >>> find_frag(15, [0, 20, 30])
    0
    >>> find_frag(15, [10, 20, 30])
    Traceback (most recent call last):
        ...
    ValueError: The first position in the restriction table is not 0.
    >>> find_frag(31, [0, 20, 30])
    Traceback (most recent call last):
        ...
    ValueError: Read position is larger than last entry in restriction table.

    """
    if r_sites[0] != 0:
        raise ValueError(
            "The first position in the restriction table is not 0."
        )
    if pos > r_sites[-1]:
        raise ValueError(
            "Read position is larger than last entry in restriction table."
        )
    # binary search for the index of the read
    index = max(np.searchsorted(r_sites, pos, side="right") - 1, 0)
    # Last site = end of the chrom, index of last fragment is last site - 1
    index = min(len(r_sites) - 2, index)

    return index


def frag_len(
    frags_file_name=DEFAULT_FRAGMENTS_LIST_FILE_NAME,
    output_dir=None,
    plot=False,
    fig_path=None,
):
    """
    logs summary statistics of fragment length distribution based on an
    input fragment file. Can optionally show a histogram instead
    of text summary.
    Parameters
    ----------
    frags_file_name : str
        Path to the output list of fragments.
    output_dir : str
        Directory where the list should be saved.
    plot : bool
        Wether a histogram of fragment length should be shown.
    fig_path : str
        If a path is given, the figure will be saved instead of shown.
    """

    try:
        frag_list_path = os.path.join(output_dir, frags_file_name)
    except TypeError:
        frag_list_path = frags_file_name
    frags = pd.read_csv(frag_list_path, sep="\t")
    nfrags = frags.shape[0]
    med_len = frags["size"].median()
    nbins = 40
    if plot:
        fig, ax = plt.subplots()
        _, _, _ = ax.hist(frags["size"], bins=nbins)

        ax.set_xlabel("Fragment length [bp]")
        ax.set_ylabel("Log10 number of fragments")
        ax.set_title("Distribution of restriction fragment length")
        ax.set_yscale("log", basey=10)
        ax.annotate(
            "Total fragments: {}".format(nfrags),
            xy=(0.95, 0.95),
            xycoords="axes fraction",
            fontsize=12,
            horizontalalignment="right",
            verticalalignment="top",
        )
        ax.annotate(
            "Median length: {}".format(med_len),
            xy=(0.95, 0.90),
            xycoords="axes fraction",
            fontsize=12,
            horizontalalignment="right",
            verticalalignment="top",
        )
        # Tweak spacing to prevent clipping of ylabel
        fig.tight_layout()
        if fig_path:
            plt.savefig(fig_path)
        else:
            plt.show()
        plt.clf()
    else:
        logger.info(
            "Genome digested into {0} fragments with a median "
            "length of {1}".format(nfrags, med_len)
        )


def gen_enzyme_religation_regex(enzyme):
    """Return a regex which corresponds to all possible religation sites given a
    set of enzyme.
    Parameters:
    -----------
    enzyme : str
        String that contains the names of the enzyme separated by a comma.
    Returns:
    --------
    re.Pattern :
        Regex that corresponds to all possible ligation sites given a set of
        enzyme.
    Examples:
    ---------
    >>> gen_enzyme_religation_regex('HpaII')
    re.compile('CCGCGG')
    >>> gen_enzyme_religation_regex('HpaII,MluCI')
    re.compile('AATTAATT|AATTCGG|CCGAATT|CCGCGG')
    """

    # Split the str on the comma to separate the different enzymes.
    enzyme = enzyme.split(",")

    # Check on Biopython dictionnary the enzyme.
    rb = RestrictionBatch(enzyme)

    # Initiation:
    give_list = []
    accept_list = []
    ligation_list = []

    # Iterates on the enzymes.
    for enz in rb:

        # Extract restriction sites and look for cut sites.
        site = enz.elucidate()
        fw_cut = site.find("^")
        rev_cut = site.find("_")

        # Process "give" site. Remove N on the left (useless).
        give_site = site[:rev_cut].replace("^", "")
        while give_site[0] == "N":
            give_site = give_site[1:]
        give_list.append(give_site)

        # Process "accept" site. Remove N on the rigth (useless).
        accept_site = site[fw_cut + 1 :].replace("_", "")
        while accept_site[-1] == "N":
            accept_site = accept_site[:-1]
        accept_list.append(accept_site)

    # Iterates on the two list to build all the possible HiC ligation sites.
    for give_site in give_list:
        for accept_site in accept_list:
            # Replace "N" by "." for regex searching of the sites
            ligation_list.append((give_site + accept_site).replace("N", "."))
            ligation_list.append(
                str(Seq(give_site + accept_site).reverse_complement()).replace(
                    "N", "."
                )
            )

    # Build the regex for any ligation sites.
    pattern = "|".join(sorted(list(set(ligation_list))))
    return re.compile(pattern)
