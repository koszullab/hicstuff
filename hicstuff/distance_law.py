#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import sys
import matplotlib.pyplot as plt
import warnings
from scipy import ndimage
from matplotlib import cm
import pandas as pd
import os as os
import csv as csv
import hicstuff.io as hio
from hicstuff.log import logger


def export_distance_law(xs, ps, names, out_file=None):
    """ Export the x(s) and p(s) from two list of numpy.ndarrays to a table
    in txt file with three columns separated by a tabulation. The first column
    contains the x(s), the second the p(s) and the third the name of the arm or
    chromosome. The file is created in the directory given by outdir or the
    current file if no file is given.

    Parameters
    ----------
    xs : list of numpy.ndarray
        The list of the start position of logbins of each p(s) in base pairs.
    ps : list of numpy.ndarray
        The list of p(s).
    names : list of string
        List containing the names of the chromosomes/arms/conditions of the p(s)
        values given.
    out_file : str or None
        Path where output file should be written. ./distance_law.txt by default.

    Return
    ------
    txt file:
        File with three columns separated by a tabulation. The first column
        contains the x(s), the second the p(s) and the third the name of the arm
        or chromosome. The file is creates in the output file given or the
        default one if none given.
    """
    # ./distance_law.txt as out_file if no out_file is given.
    if out_file is None:
        out_file = os.getcwd() + "/distance_law.txt"
    # Sanity check: as many chromosomes/arms as ps
    if len(xs) != len(names):
        logger.error("Number of chromosomes/arms and number of p(s) list differ.")
        sys.exit(1)
    # Create the file and write it
    f = open(out_file, "w")
    for i in range(len(xs)):
        for j in range(len(xs[i])):
            line = (
                str(format(xs[i][j], "g"))
                + "\t"
                + str(format(ps[i][j], "g"))
                + "\t"
                + names[i]
                + "\n"
            )
            f.write(line)
    f.close()


def import_distance_law(distance_law_file):
    """ Import the table created by export_distance_law and return the list of
    x(s) and p(s) in the order of the chromosomes.

    Parameters
    ----------
    distance_law_file : string
        Path to the file containing three columns : the x(s), the p(s), and the
        chromosome/arm name.

    Returns
    -------
    list of numpy.ndarray :
        The start coordinate of each bin one array per chromosome or arm.
    list of numpy.ndarray :
        The distance law probabilities corresponding of the bins of the
        previous list.
    list of numpy.ndarray :
        The names of the arms/chromosomes corresponding to the previous
        list.
    """
    file = pd.read_csv(
        distance_law_file,
        sep="\t",
        header=None,
        dtype={"a": np.int32, "b": np.float32, "c": str},
    )
    names_idx = np.unique(file.iloc[:, 2], return_index=True)[1]
    names = [file.iloc[:, 2][index] for index in sorted(names_idx)]
    xs = [None] * len(names)
    ps = [None] * len(names)
    labels = [None] * len(names)
    for i in range(len(names)):
        subfile = file[file.iloc[:, 2] == names[i]]
        xs[i] = np.array(subfile.iloc[:, 0])
        ps[i] = np.array(subfile.iloc[:, 1])
        labels[i] = np.array(subfile.iloc[:, 2])
    return xs, ps, labels


def get_chr_segment_bins_index(fragments, centro_file=None, rm_centro=0):
    """Get the index positions of the start and end bins of different 
    chromosomes, or arms if the centromers position have been given from the
    fragments file made by hicstuff.
    
    Parameters
    ----------
    fragments : pandas.DataFrame
        Table containing in the first coulum the ID of the fragment, in the 
        second the names of the chromosome in the third and fourth the start 
        position and the end position of the fragment. The file have no header.
        (File like the 'fragments_list.txt' from hicstuff)
    centro_file : None or str
        None or path to a file with the genomic positions of the centromers 
        sorted as the chromosomes separated by a space. The file have only one 
        line.
    rm_centro : int
        If a value is given, will remove the contacts close the centromeres.
        It will remove as many kb as the argument given. Default is zero.
        
    Returns
    -------
    list of floats :
        The start and end indices of chromosomes/arms to compute the distance
        law on each chromosome/arm separately.
    """
    # Get bins where chromosomes start
    chr_start_bins = np.where(fragments == 0)[0]
    # Create a list of same length for the end of the bins
    chr_end_bins = np.zeros(len(chr_start_bins))
    # Get bins where chromsomes end
    for i in range(len(chr_start_bins) - 1):
        chr_end_bins[i] = chr_start_bins[i + 1]
    chr_end_bins[-1] = len(fragments.iloc[:, 0])
    # Combine start and end of bins in a single array. Values are the id of the
    # bins
    chr_segment_bins = np.sort(np.concatenate((chr_start_bins, chr_end_bins)))
    if centro_file is not None:
        # Read the file of the centromers
        with open(centro_file, "r", newline="") as centro:
            centro = csv.reader(centro, delimiter=" ")
            centro_pos = next(centro)
        # Sanity check: as many chroms as centromeres
        if len(chr_start_bins) != len(centro_pos):
            logger.warning(
                "Number of chromosomes and centromeres differ, centromeres position are not taking into account."
            )
            centro_file = None
    if centro_file is not None:
        # Get bins of centromeres
        centro_bins = np.zeros(2 * len(centro_pos))
        for i in range(len(chr_start_bins)):
            if (i + 1) < len(chr_start_bins):
                subfrags = fragments[chr_start_bins[i] : chr_start_bins[i + 1]]
            else:
                subfrags = fragments[chr_start_bins[i] :]
            # index of last fragment starting before centro in same chrom
            centro_bins[2 * i] = chr_start_bins[i] + max(
                np.where(
                    subfrags["start_pos"][:] // (int(centro_pos[i]) - rm_centro) == 0
                )[0]
            )
            centro_bins[2 * i + 1] = chr_start_bins[i] + max(
                np.where(
                    subfrags["start_pos"][:] // (int(centro_pos[i]) + rm_centro) == 0
                )[0]
            )
        # Combine centro and chrom bins into a single array. Values are the id
        # of the bins started and ending the arms.
        chr_segment_bins = np.sort(
            np.concatenate((chr_start_bins, chr_end_bins, centro_bins))
        )
    return list(chr_segment_bins)


def get_chr_segment_length(fragments, chr_segment_bins):
    """Compute a list of the length of the different objects (arm or 
    chromosome) given by chr_segment_bins.
    
    Parameters
    ----------
    fragments : pandas.DataFrame
        Table containing in the first coulum the ID of the fragment, in the 
        second the names of the chromosome in the third and fourth the start 
        position and the end position of the fragment. The file have no header.
        (File like the 'fragments_list.txt' from hicstuff)
    chr_segment_bins : list of floats
        The start and end indices of chromosomes/arms to compute the distance
        law on each chromosome/arm separately.
        
    Returns
    -------
    list of numpy.ndarray:
        The length in base pairs of each chromosome or arm.
    """
    chr_segment_length = [None] * int(len(chr_segment_bins) / 2)
    # Iterate in chr_segment_bins in order to obtain the size of each chromosome/arm
    for i in range(len(chr_segment_length)):
        # Obtain the size of the chromosome/arm, the if loop is to avoid the
        # case of arms where the end position of the last fragments doesn't
        # mean the size of arm. If it's the right arm we have to start to count the
        # size from the beginning of the arm.
        if fragments["start_pos"].iloc[int(chr_segment_bins[2 * i])] == 0:
            n = fragments["end_pos"].iloc[int(chr_segment_bins[2 * i + 1]) - 1]
        else:
            n = (
                fragments["end_pos"].iloc[int(chr_segment_bins[2 * i + 1]) - 1]
                - fragments["start_pos"].iloc[int(chr_segment_bins[2 * i])]
            )
        chr_segment_length[i] = n
    return chr_segment_length


def logbins_xs(fragments, chr_segment_length, base=1.1, circular=False):
    """Compute the logbins of each chromosome/arm in order to have theme to
    compute distance law. At the end you will have bins of increasing with a 
    logspace with the base of the value given in base.
    
    Parameters
    ----------
    fragments : pandas.DataFrame
        Table containing in the first coulum the ID of the fragment, in the 
        second the names of the chromosome in the third and fourth the start 
        position and the end position of the fragment. The file have no header.
        (File like the 'fragments_list.txt' from hicstuff)
    chr_segment_length: list of floats
        List of the size in base pairs of the different arms or chromosomes.
    base : float
        Base use to construct the logspace of the bins, 1.1 by default.
    circular : bool
        If True, calculate the distance as the chromosome is circular. Default 
        value is False.
        
    Returns
    -------
    list of numpy.ndarray :
        The start coordinate of each bin one array per chromosome or arm.
    """
    # Create the xs array and a list of the length of the chromosomes/arms
    xs = [None] * len(chr_segment_length)
    # Iterate in chr_segment_bins in order to make the logspace
    for i in range(len(chr_segment_length)):
        n = chr_segment_length[i]
        # if the chromosome is circular the mawimum distance between two reads
        # are divided by two
        if circular:
            n /= 2
        n_bins = int(np.log(n) / np.log(base))
        # For each chromosome/arm compute a logspace to have the logbin
        # equivalent to the size of the arms and increasing size of bins
        xs[i] = np.unique(
            np.logspace(0, n_bins, num=n_bins + 1, base=base, dtype=np.int)
        )
    return xs


def circular_distance_law(distance, chr_segment_length, chr_bin):
    """Recalculate the distance to return the distance in a circular chromosome
    and not the distance between the two genomic positions.

    Parameters
    ----------
    chr_segment_bins : list of floats
        The start and end indices of chromosomes/arms to compute the distance
        law on each chromosome/arm separately.
    chr_segment_length: list of floats
        List of the size in base pairs of the different arms or chromosomes.
    distance : int
        Distance between two fragments with a contact.

    Returns
    -------
    int :
        The real distance in the chromosome circular and not the distance 
        between two genomic positions

    Examples
    --------
    >>> circular_distance_law(7500, [2800, 9000], 1)
    1500
    >>> circular_distance_law(1300, [2800, 9000], 0)
    1300
    >>> circular_distance_law(1400, [2800, 9000], 0)
    1400
    """
    chr_len = chr_segment_length[chr_bin]
    if distance > chr_len / 2:
        distance = chr_len - distance
    return distance


def get_pairs_distance(
    line, fragments, chr_segment_bins, chr_segment_length, xs, ps, circular=False
):
    """From a line of a pair reads file, filter -/+ or +/- reads, keep only the 
    reads in the same chromosome/arm and compute the distance of the the two
    fragments. It modify the input ps in order to count or not the line given. 
    It will add one in the logbin corresponding to the distance.

    Parameters
    ----------
    line : OrderedDict
        Line of a pair reads file with the these keys readID, chr1, pos1, chr2,
        pos2, strand1, strand2, frag1, frag2. The values are in a dictionnary.
    fragments : pandas.DataFrame
        Table containing in the first coulum the ID of the fragment, in the 
        second the names of the chromosome in the third and fourth the start 
        position and the end position of the fragment. The file have no header.
        (File like the 'fragments_list.txt' from hicstuff)
    chr_segment_bins : list of floats
        The start and end indices of chromosomes/arms to compute the distance
        law on each chromosome/arm separately.
    chr_segment_length: list of floats
        List of the size in base pairs of the different arms or chromosomes.
    xs : list of lists
        The start coordinate of each bin one array per chromosome or arm.
    ps : list of lists
        The sum of contact already count. xs and ps should have the same 
        dimensions.
    circular : bool
        If True, calculate the distance as the chromosome is circular. Default 
        value is False.
    """
    # Check this is a pairs_idx file and not simple pairs
    if line['frag1'] is None:
        logger.error(
            'Input pairs file must have frag1 and frag2 columns. In hicstuff '
            'pipeline, this is the "valid_idx.pairs" file.'
        )
    # We only keep the event +/+ or -/-. This is done to avoid to have any
    # event of uncut which are not possible in these events. We can remove the
    # good events of +/- or -/+ because we don't need a lot of reads to compute
    # the distance law and if we eliminate these reads we do not create others
    # biases as they should have the same distribution.
    if line["strand1"] == line["strand2"]:
        # Find in which chromosome/arm are the fragment 1 and 2.
        chr_bin1 = (
            np.searchsorted(chr_segment_bins, int(line["frag1"]), side="right") - 1
        )
        chr_bin2 = (
            np.searchsorted(chr_segment_bins, int(line["frag2"]), side="right") - 1
        )
        # We only keep the reads with the two fragments in the same chromosome
        # or arm.
        if chr_bin1 == chr_bin2:
            # Remove the contacts in the centromeres if centro_remove
            if chr_bin1 % 2 == 0:
                chr_bin1 = int(chr_bin1 / 2)
                # For the reads -/-, the fragments should be religated with both
                # their start position (position in the left on the genomic
                # sequence, 5'). For the reads +/+ it's the contrary. We compute
                # the distance as the distance between the two extremities which
                # are religated.
                if line["strand1"] == "-":
                    distance = abs(
                        np.array(fragments["start_pos"][int(line["frag1"])])
                        - np.array(fragments["start_pos"][int(line["frag2"])])
                    )
                if line["strand1"] == "+":
                    distance = abs(
                        np.array(fragments["end_pos"][int(line["frag1"])])
                        - np.array(fragments["end_pos"][int(line["frag2"])])
                    )
                if circular:
                    distance = circular_distance_law(
                        distance, chr_segment_length, chr_bin1
                    )
                xs_temp = xs[chr_bin1][:]
                # Find the logbins in which the distance is and add one to the sum
                # of contact.
                ps_indice = np.searchsorted(xs_temp, distance, side="right") - 1
                ps[chr_bin1][ps_indice] += 1


def get_names(fragments, chr_segment_bins):
    """Make a list of the names of the arms or the chromosomes.

    Parameters
    ----------
    fragments : pandas.DataFrame
        Table containing in the first coulum the ID of the fragment, in the 
        second the names of the chromosome in the third and fourth the start 
        position and the end position of the fragment. The file have no header.
        (File like the 'fragments_list.txt' from hicstuff)
    chr_segment_bins : list of floats
        The start position of chromosomes/arms to compute the distance law on 
        each chromosome/arm separately.

    Returns
    -------
    list of floats : 
        List of the labels given to the curves. It will be the name of the arms
        or chromosomes.
    """
    # Get the name of the chromosomes.
    chr_names_idx = np.unique(fragments.iloc[:, 1], return_index=True)[1]
    chr_names = [fragments.iloc[:, 1][index] for index in sorted(chr_names_idx)]
    # Case where they are separate in chromosomes
    if len(chr_segment_bins) / 2 != len(chr_names):
        names = []
        for chr in chr_names:
            names.append(str(chr) + "_left")
            names.append(str(chr) + "_rigth")
        chr_names = names
    return chr_names


def get_distance_law(
    pairs_reads_file,
    fragments_file,
    centro_file=None,
    base=1.1,
    out_file=None,
    circular=False,
    rm_centro=0,
):
    """Compute distance law as a function of the genomic coordinate aka P(s).
    Bin length increases exponentially with distance. Works on pairs file 
    format from 4D Nucleome Omics Data Standards Working Group. If the genome 
    is composed of several chromosomes and you want to compute the arms 
    separately, provide a file with the positions of centromers. Create a file 
    with three coulumns separated by a tabulation. The first column contains 
    the xs, the second the ps and the third the name of the arm or chromosome. 
    The file is create in the directory given in outdir or in the current 
    directory if no directory given.

    Parameters
    ----------
    pairs_reads_file : string
        Path of a pairs file format from 4D Nucleome Omics Data Standards 
        Working Group with the 8th and 9th coulumns are the ID of the fragments
        of the reads 1 and 2.
    fragments_file : path
        Path of a table containing in the first column the ID of the fragment,
        in the second the names of the chromosome in the third and fourth 
        the start position and the end position of the fragment. The file have 
        no header. (File like the 'fragments_list.txt' from hicstuff)
    centro_file : None or str
        None or path to a file with the genomic positions of the centromers 
        sorted as the chromosomes separated by a space. The file have only one 
        line.
    base : float
        Base use to construct the logspace of the bins - 1.1 by default.
    out_file : None or str
        Path of the output file. If no path given, the output is returned.
    circular : bool
        If True, calculate the distance as the chromosome is circular. Default 
        value is False. Cannot be True if centro_file is not None
    rm_centro : int
        If a value is given, will remove the contacts close the centromeres.
        It will remove as many kb as the argument given. Default is None.

    Returns
    -------
    xs : list of numpy.ndarray
        Basepair coordinates of log bins used to compute distance law.
    ps : list of numpy.ndarray
        Contacts value, in arbitrary units, at increasingly long genomic ranges
        given by xs.
    names : list of strings
        Names of chromosomes that are plotted
    """
    # Sanity check : centro_fileition should be None if chromosomes are
    # circulars (no centromeres is circular chromosomes).
    if circular and centro_file != None:
        logger.error("Chromosomes cannot have a centromere and be circular")
        sys.exit(1)
    # Import third columns of fragments file
    fragments = pd.read_csv(fragments_file, sep="\t", header=0, usecols=[0, 1, 2, 3])
    # Calculate the indice of the bins to separate into chromosomes/arms
    chr_segment_bins = get_chr_segment_bins_index(fragments, centro_file, rm_centro)
    # Calculate the length of each chromosoms/arms
    chr_segment_length = get_chr_segment_length(fragments, chr_segment_bins)
    xs = logbins_xs(fragments, chr_segment_length, base, circular)
    # Create the list of p(s) with one array for each chromosome/arm and each
    # array contain as many values as in the logbin
    ps = [None] * len(chr_segment_length)
    for i in range(len(xs)):
        ps[i] = [0] * len(xs[i])
    # Read the pair reads file
    with open(pairs_reads_file, "r", newline="") as reads:
        # Remove the line of the header
        header_length = len(hio.get_pairs_header(pairs_reads_file))
        for i in range(header_length):
            next(reads)
        # Reads all the others lines and put the values in a dictionnary with
        # the keys : 'readID', 'chr1', 'pos1', 'chr2', 'pos2', 'strand1',
        # 'strand2', 'frag1', 'frag2'
        reader = csv.DictReader(
            reads,
            fieldnames=[
                "readID",
                "chr1",
                "pos1",
                "chr2",
                "pos2",
                "strand1",
                "strand2",
                "frag1",
                "frag2",
            ],
            delimiter="\t",
        )
        for line in reader:
            # Iterate in each line of the file after the header
            get_pairs_distance(
                line, fragments, chr_segment_bins, chr_segment_length, xs, ps, circular
            )
    # Divide the number of contacts by the area of the logbin
    for i in range(len(xs)):
        n = chr_segment_length[i]
        for j in range(len(xs[i]) - 1):
            # Use the area of a trapezium to know the area of the logbin with n
            # the size of the matrix.
            ps[i][j] /= ((2 * n - xs[i][j + 1] - xs[i][j]) / 2) * (
                (1 / np.sqrt(2)) * (xs[i][j + 1] - xs[i][j])
            )
            # print(
            #    ((2 * n - xs[i][j + 1] - xs[i][j]) / 2)
            #    * ((1 / np.sqrt(2)) * (xs[i][j + 1] - xs[i][j]))
            # )
        # Case of the last logbin which is an isosceles rectangle triangle
        # print(ps[i][-5:-1], ((n - xs[i][-1]) ** 2) / 2)
        ps[i][-1] /= ((n - xs[i][-1]) ** 2) / 2
    names = get_names(fragments, chr_segment_bins)
    if out_file:
        export_distance_law(xs, ps, names, out_file)
    return xs, ps, names


def normalize_distance_law(xs, ps, inf=3000, sup=None):
    """Normalize the distance in order to have the sum of the ps values between
    'inf' (default value is 3kb) until the end of the array equal to one and
    limit the effect of coverage between two conditions/chromosomes/arms when
    you compare them together. If we have a list of ps, it will normalize until
    the length of the shorter object or the value of sup, whichever is smaller.

    Parameters
    ----------
    xs : list of numpy.ndarray
        list of logbins corresponding to the ps.
    ps : list of numpy.ndarray
        Average ps or list of ps of the chromosomes/arms. xs and ps have to 
        have the same shape.
    inf : integer
        Inferior value of the interval on which, the normalization is applied.
    sup : integer
        Superior value of the interval on which, the normalization is applied.

    Returns
    -------
    list of numpy.ndarray :
        List of ps each normalized separately.
    """
    # Sanity check: xs and ps have the same dimension
    if np.shape(xs) != np.shape(ps):
        logger.error("xs and ps should have the same dimension.")
        sys.exit(1)
    # Define the length of shortest chromosomes as a lower bound for the sup boundary
    min_xs = len(min(xs, key=len))
    normed_ps = [None] * len(ps)
    if sup is None:
        sup = np.inf
    for chrom_id, chrom_ps in enumerate(ps):
        # Iterate on the different ps to normalize each of theme separately
        chrom_sum = 0
        # Change the last value to have something continuous because the last
        # one is much bigger (computed on matrix corner = triangle instead of trapezoid).
        chrom_ps[-1] = chrom_ps[-2]
        for bin_id, bin_value in enumerate(chrom_ps):
            # Compute normalization factor based on values between inf and sup
            # Sup will be whatever is smaller between user-provided sup and length of
            # the shortest chromosome
            if (xs[chrom_id][bin_id] > inf) and (xs[chrom_id][bin_id] < sup) and (bin_id < min_xs):
                chrom_sum += bin_value
        if chrom_sum == 0:
            chrom_sum += 1
            logger.warning("No values of p(s) in one segment")
        # Make the normalisation
        normed_ps[chrom_id] = np.array(ps[chrom_id]) / chrom_sum
    return normed_ps


def average_distance_law(xs, ps, sup, big_arm_only=False):
    """Compute the average distance law between the file the different distance
    law of the chromosomes/arms.

    Parameters
    ----------
    xs : list of numpy.ndarray
        The list of logbins.
    ps : list of lists of floats
        The list of numpy.ndarray.
    sup : int
        Value given to set the minimum size of the chromosomes/arms to make the
        average.
    big_arm_only : bool
        By default False. If True, will only take into account the arms/chromosomes 
        longer than the value of sup. Sup mandatory if set.

    Returns
    -------
    numpy.ndarray :
        List of the xs with the max length.
    numpy.ndarray :
        List of the average_ps.
    """
    # Find longest chromosome / arm and make two arrays of this length for the
    # average distance law and remove the last value.
    xs = max(xs, key=len)
    max_length = len(xs)
    ps_values = np.zeros(max_length)
    ps_occur = np.zeros(max_length)
    for chrom_ps in ps:
        # Iterate on ps in order to calculate the number of occurences (all the
        # chromossomes/arms are not as long as the longest one) and the sum of
        # the values of distance law.
        # Change the last value to have something continuous because the last
        # one is much bigger.
        chrom_ps[-1] = chrom_ps[-2]
        # Sanity check : sup strictly inferior to maw length arms.
        if big_arm_only:
            if sup >= xs[-1]:
                logger.error(
                    "sup have to be inferior to the max length of arms/chromsomes if big arm only set"
                )
                sys.exit(1)
            if sup <= xs[len(chrom_ps) - 1]:
                ps_occur[: len(chrom_ps)] += 1
                ps_values[: len(chrom_ps)] += chrom_ps
        else:
            ps_occur[: len(chrom_ps)] += 1
            ps_values[: len(chrom_ps)] += chrom_ps
    # Make the mean
    averaged_ps = ps_values / ps_occur
    return xs, averaged_ps


def slope_distance_law(xs, ps):
    """Compute the slope of the loglog curve of the ps as the 
    [log(ps(n+1)) - log(ps(n))] / [log(n+1) - log(n)].
    Compute only list of ps, not list of array.

    Parameters
    ----------
    xs : list of numpy.ndarray
        The list of logbins.
    ps : list of numpy.ndarray
        The list of ps.

    Returns
    -------
    list of numpy.ndarray :
        The slope of the distance law. It will be shorter of one value than the
        ps given initially.
    """
    slope = [None] * len(ps)
    for i in range(len(ps)):
        ps[i][ps[i] == 0] = 10 ** (-9)
        # Compute the slope
        slope_temp = np.log(np.array(ps[i][1:]) / np.array(ps[i][:-1])) / np.log(
            np.array(xs[i][1:]) / np.array(xs[i][:-1])
        )
        # The 1.8 is the intensity of the normalisation, it could be adapted.
        slope_temp[slope_temp == np.nan] = 10 ** (-15)
        slope[i] = ndimage.filters.gaussian_filter1d(slope_temp, 1.8)
    return slope


def get_ylim(xs, curve, inf, sup):
    """Compute the max and the min of the list of list between the inferior
    (inf) and superior (sup) bounds.

    Parameters
    ----------
    xs : list of numpy.ndarray
        The list of the logbins starting position in basepair.
    curve : list of numpy.ndarray
        A list of numpy.ndarray from which you want to extract minimum and
        maximum values in a given interval.
    inf : int
        Inferior limit of the interval in basepair.
    sup : int
        Superior limit of the interval in basepair.

    Returns
    -------
    min_tot : float
        Minimum value of the list of list in this interval.
    max_tot : float
        Maximum value of the list of list in this interval.

    Examples
    --------
    >>> get_ylim([np.array([1, 4, 15]), np.array([1, 4, 15, 40])],
    ... [np.array([5.5, 3.2, 17.10]), np.array([24, 32, 1.111, 18.5])],
    ... 2,
    ... 15
    ... )
    (1.111, 32.0)
    """
    # Create a list in order to add all the interesting values
    flatten_list = []
    for i, logbins in enumerate(xs):
        # Iterate on xs.
        # Search for the minimum index corresponding to the smallest bin
        # superior or equal to inf (in pair base).
        min_value = min(logbins[logbins >= inf], default=0)
        min_index = np.where(logbins == min_value)[0]
        # Skip chromosome if total size smaller than inf
        if len(min_index) == 0:
            continue
        # Search for the maximum index corresponding to the biggest bin
        # inferior or equal to sup (in pair base).
        max_value = max(logbins[logbins <= sup])
        max_index = np.where(logbins == max_value)[0]
        # Add the values in the interval in the flattened list.
        if int(max_index) != len(logbins) - 1:
            max_index += 1
        for j in range(int(min_index), int(max_index)):
            flatten_list.append(curve[i][j])
    # Caluclate the min and the max of this list.
    min_tot = min(flatten_list)
    max_tot = max(flatten_list)
    return min_tot, max_tot


def plot_ps_slope(xs, ps, labels, fig_path=None, inf=3000, sup=None):
    """Compute two plots, one with the different distance law of each 
    arm/chromosome in loglog scale and one with the slope (derivative) of these
    curves. Generate a svg file with savefig.

    Parameters
    ----------
    xs : list of numpy.ndarray
        The list of the logbins of each ps.
    ps : list of numpy.ndarray
        The list of ps.
    labels_file : list of strings
        Names of the different curves in the order in which they are given.
    fig_path : str
        Path where the figure will be created. If None (default), the plot is
        shown in an interactive window.
    inf : int 
        Value of the mimimum x of the window of the plot. Have to be strictly
        positive. By default 3000.
    sup : int 
        Value of the maximum x of the window of the plot. By default None.
        
    """
    # Give the max value for sup if no value have been attributed
    if sup is None:
        sup = max(max(xs, key=len))

    # Compute slopes from the curves
    slope = slope_distance_law(xs, ps)
    # Make the plot of distance law
    # Give a range of color
    cols = iter(cm.rainbow(np.linspace(0, 1, len(ps))))
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, figsize=(18, 10))
    plt.subplots_adjust(left=0.05, right=0.85, top=0.93, bottom=0.07)
    ax1.set_xlabel("Distance (pb)", fontsize="x-large")
    ax1.set_ylabel("P(s)", fontsize="x-large")
    ax1.set_title("Distance law", fontsize="xx-large")
    ylim = get_ylim(xs, ps, inf, sup)
    ax1.set_ylim(0.9 * ylim[0], 1.1 * ylim[1])
    for i in range(len(ps)):
        # Iterate on the different distance law array and take them by order of
        # size in order to have the color scale equivalent to the size scale
        col = next(cols)
        ax1.loglog(xs[i], ps[i], label=labels[i])
    # Make the same plot with the slope
    cols = iter(cm.rainbow(np.linspace(0, 1, len(slope))))
    ax2.set_xlabel("Distance (pb)", fontsize="x-large")
    ax2.set_ylabel("Slope", fontsize="x-large")
    ax2.set_title("Slope of the distance law", fontsize="xx-large")
    ax2.set_xlim([inf, sup])
    ylim = get_ylim(xs, slope, inf, sup)
    ax2.set_ylim(1.1 * ylim[0], 0.9 * ylim[1])
    xs2 = [None] * len(xs)
    for i in range(len(slope)):
        xs2[i] = xs[i][:-1]
        col = next(cols)
        ax2.semilogx(xs2[i], slope[i], label=labels[i], subsx=[2, 3, 4, 5, 6, 7, 8, 9])
    ax2.legend(loc="upper left", bbox_to_anchor=(1.02, 1.00), ncol=1, fontsize="large")
    # Save the figure in svg
    if fig_path is not None:
        plt.savefig(fig_path)
    else:
        plt.show()
