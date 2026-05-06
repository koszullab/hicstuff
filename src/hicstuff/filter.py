#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Hi-C event filtering

Analyse the contents of a 3C library and filters spurious events such as loops
and uncuts to improve the overall signal. Filtering consists in excluding +-
and -+ Hi-C pairs if their reads are closer than a threshold in minimum number
of restriction fragments. This threshold represents the distance at which the
abundance of these events deviate significantly from the rest of the library.
It is estimated automatically by default using the median absolute deviation of
pairs at longer distances, but can also be entered manually.

The program takes a 2D BED file as input with the following fields:
chromA startA endA indiceA strandA chromB startB endB indiceB strandB
Each line of the file represents a Hi-C pair with reads A and B. The indices
are 0-based and represent the restriction fragment to which reads are
attributed.

@author: cmdoret (reimplementation of Axel KournaK's code)
"""

import sys
import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt
from hicstuff.log import logger


def process_read_pair(line):
    r"""Process and order read pairs in a .pairs record.

    Takes a read pair (line) from a .pairs file as input, reorders the pair
    so that read 1 in intrachromosomal pairs always has the smallest genomic
    coordinate.

    Parameters
    ----------
    line : str
        Record from a pairs file with columns: readID, chr1, pos1, chr2,
        pos2, strand1, strand2, frag1, frag2

    Returns
    -------
    dict
        Dictionary with reordered pair where each column from the line as an
        entry. The number of restriction sites separating reads in the pair and
        the type of event are also stored in the dictionary, under keys
        "nsites" and "type".

    Examples
    --------
        >>> d = process_read_pair("readX\ta\t1\tb\t20\t-\t-\t1\t3")
        >>> for u in sorted(d.items()):
        ...     print(u)
        ('chr1', 'a')
        ('chr2', 'b')
        ('frag1', 1)
        ('frag2', 3)
        ('nsites', 2)
        ('pos1', 1)
        ('pos2', 20)
        ('readID', 'readX')
        ('strand1', '-')
        ('strand2', '-')
        ('type', 'inter')
        >>> d = process_read_pair('readY\ta\t20\ta\t10\t-\t+\t2\t1')
        >>> [d[x] for x in sorted(d.keys())]
        ['a', 'a', 1, 2, 1, 10, 20, 'readY', '+', '-', '+-']
    """
    # Split line by whitespace
    p = line.split("\t")
    if len(p) != 9:
        raise ValueError(
            "Your input file does not have 9 columns. Make sure "
            "the file has the readID and twice the following 4 fields "
            "(once for each read in the pair): chr pos strand frag."
        )
    # Saving each column as a dictionary entry with meaningful key
    cols = [
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
    p = {cols[i]: p[i] for i in range(len(cols))}
    # Transforming numeric columns to int
    for col in ["pos1", "pos2", "frag1", "frag2"]:
        p[col] = int(p[col])

    # invert records for intrachromosomal pairs where rec2 comes before
    # rec1 in genomic coordinates
    if p["chr1"] == p["chr2"] and p["pos2"] < p["pos1"]:
        p["strand1"], p["strand2"] = p["strand2"], p["strand1"]
        p["pos1"], p["pos2"] = p["pos2"], p["pos1"]
        p["frag1"], p["frag2"] = p["frag2"], p["frag1"]

    # Number of restriction sites separating reads in the pair
    p["nsites"] = abs(p["frag2"] - p["frag1"])
    # Get event type
    if p["chr1"] == p["chr2"]:
        p["type"] = "".join([p["strand1"], p["strand2"]])
    else:
        p["type"] = "inter"

    return p


def get_thresholds(
    in_dat, interactive=False, plot_events=False, fig_path=None, prefix=None
):
    """Guess distance threshold for event filtering

    Analyse the events in the first million of Hi-C pairs in the library, plot
    the occurrences of each event type according to number of restriction
    fragments, and ask user interactively for the minimum threshold for uncuts
    and loops.

    Parameters
    ----------
    in_dat: str
        Path to the .pairs file containing Hi-C pairs.
    interactive: bool
        If True, plots are diplayed and thresholds are required interactively.
    plot_events : bool
        Whether to show the plot
    fig_path : str
        Path where the figure will be saved. If None, the figure will be
        diplayed interactively.
    prefix : str
        If the library has a name, it will be shown on plots.
    
    Returns
    -------
    dictionary
        dictionary with keys "uncuts" and "loops" where the values are the
        corresponding thresholds entered by the user.
    """
    thr_uncut = None
    thr_loop = None
    max_sites = 50
    # Map of event -> legend name of event for intrachromosomal pairs.
    legend = {
        "++": "++ (weird)",
        "--": "-- (weird)",
        "+-": "+- (uncuts)",
        "-+": "-+ (loops)",
    }
    colors = {"++": "#222222", "+-": "r", "--": "#666666", "-+": "tab:orange"}
    n_events = {event: np.zeros(max_sites) for event in legend}
    i = 0
    # open the file for reading (just the first 1 000 000 lines)
    with open(in_dat, "r") as pairs:
        for line in pairs:
            # Skip header lines
            if line.startswith("#"):
                continue
            i += 1
            # Only use the first million pairs to estimate thresholds
            if i == 1000000:
                break
            # Process Hi-C pair into a dictionary
            p = process_read_pair(line)
            # Type of event and number of restriction site between reads
            etype = p["type"]
            nsites = p["nsites"]
            # Count number of events for intrachrom pairs
            if etype != "inter" and nsites < max_sites:
                n_events[etype][nsites] += 1

    def plot_event(n_events, legend, name):
        """Plot the frequency of a given event types over distance."""
        plt.xlim([-0.5, 15])
        plt.plot(
            range(n_events[name].shape[0]),
            n_events[name],
            "o-",
            label=legend[name],
            linewidth=2.0,
            c=colors[name],
        )

    if interactive:
        # PLot:
        try:
            plt.figure(0)
            for event in legend:
                plot_event(n_events, legend, event)
            plt.grid()
            plt.xlabel("Number of restriction fragment(s)")
            plt.ylabel("Number of events")
            plt.yscale("log")
            plt.legend()
            plt.show(block=False)

        except Exception:
            logger.error(
                "Unable to show plots, skipping figure generation. Perhaps "
                "there is no Xserver running ? (might be due to windows "
                "environment). Try running without the interactive option."
            )

        # Asks the user for appropriate thresholds
        print(
            "Please enter the number of restriction fragments separating "
            "reads in a Hi-C pair below or at which loops and "
            "uncuts events will be excluded\n",
            file=sys.stderr,
        )
        thr_uncut = int(input("Enter threshold for the uncuts events (+-):"))
        thr_loop = int(input("Enter threshold for the loops events (-+):"))
        try:
            plt.clf()
        except Exception:
            pass
    else:
        # Estimate thresholds from data
        for event in n_events:
            fixed = n_events[event]
            fixed[fixed == 0] = 1
            n_events[event] = fixed

        all_events = np.log(np.array(list(n_events.values())))
        # Compute median occurences at each restriction sites
        event_med = np.median(all_events, axis=0)
        # Compute MAD, to have a robust estimator of the expected deviation
        # from median at long distances
        mad = np.median(abs(all_events - event_med))
        exp_stdev = mad / 0.67449
        # Iterate over sites, from furthest to frag+2
        for site in range(max_sites)[:1:-1]:
            # For uncuts and loops, keep the last (closest) site where the
            # deviation from other events <= expected_stdev
            if (
                abs(np.log(n_events["+-"][site]) - event_med[site])
                <= exp_stdev
            ):
                thr_uncut = site
            if (
                abs(np.log(n_events["-+"][site]) - event_med[site])
                <= exp_stdev
            ):
                thr_loop = site
        if thr_uncut is None or thr_loop is None:
            raise ValueError(
                "The threshold for loops or uncut could not be estimated. "
                "Please try running with -i to investigate the problem."
            )
        logger.info(
            "Filtering with thresholds: uncuts={0} loops={1}".format(
                thr_uncut, thr_loop
            )
        )
        if plot_events:
            try:
                plt.figure(1)
                plt.xlim([-0.5, 15])
                # Draw colored lines for events to discard
                plt.plot(
                    range(0, thr_uncut + 1),
                    n_events["+-"][: thr_uncut + 1],
                    "o-",
                    c=colors["+-"],
                    label=legend["+-"],
                )
                plt.plot(
                    range(0, thr_loop + 1),
                    n_events["-+"][: thr_loop + 1],
                    "o-",
                    c=colors["-+"],
                    label=legend["-+"],
                )
                plt.plot(
                    range(0, 2),
                    n_events["--"][:2],
                    "o-",
                    c=colors["--"],
                    label=legend["--"],
                )
                plt.plot(
                    range(0, 2),
                    n_events["++"][:2],
                    "o-",
                    c=colors["++"],
                    label=legend["++"],
                )
                # Draw black lines for events to keep
                plt.plot(
                    range(thr_uncut, n_events["+-"].shape[0]),
                    n_events["+-"][thr_uncut:],
                    "o-",
                    range(thr_loop, n_events["-+"].shape[0]),
                    n_events["-+"][thr_loop:],
                    "o-",
                    range(1, n_events["--"].shape[0]),
                    n_events["--"][1:],
                    "o-",
                    range(1, n_events["++"].shape[0]),
                    n_events["++"][1:],
                    "o-",
                    label="kept",
                    linewidth=2.0,
                    c="g",
                )
                plt.grid()
                plt.xlabel("Number of restriction site(s)")
                plt.ylabel("Number of events")
                plt.yscale("log")
                # Remove duplicate "kept" entries in legend
                handles, labels = plt.gca().get_legend_handles_labels()
                by_label = OrderedDict(zip(labels, handles))
                plt.legend(by_label.values(), by_label.keys())
                # Show uncut and loop threshold as vertical lines
                plt.axvline(x=thr_loop, color=colors["-+"])
                plt.axvline(x=thr_uncut, color=colors["+-"])

                if prefix:
                    plt.title(
                        "Library events by distance in {}".format(prefix)
                    )
                plt.tight_layout()
                if fig_path:
                    plt.savefig(fig_path)
                else:
                    plt.show(block=False)
                # plt.clf()

            except Exception:
                logger.error(
                    "Unable to show plots, skipping figure generation. Is "
                    "an X server running? (might be due to windows "
                    "environment). Try running without the plot option."
                )
    return thr_uncut, thr_loop


def filter_events(
    in_dat,
    out_filtered,
    thr_uncut,
    thr_loop,
    plot_events=False,
    fig_path=None,
    prefix=None,
):
    """Filter events (loops, uncuts and weirds)

    Filter out spurious intrachromosomal Hi-C pairs from input file. +- pairs
    with reads closer or at the uncut threshold and -+ pairs with reads closer
    or at the loop thresholds are excluded from the ouput file. -- and ++ pairs
    with both mates on the same fragments are also discarded. All others are
    written.

    Parameters
    ----------
    in_dat : file object
        File handle in read mode to the 2D BED file containing Hi-C pairs.
    out_filtered : file object
        File handle in write mode the output filtered 2D BED file.
    thr_uncut : int
        Minimum number of restriction sites between reads to keep an
        intrachromosomal +- pair.
    thr_loop : int
        Minimum number of restriction sites between reads to keep an
        intrachromosomal -+ pair.
    plot_events : bool
        If True, a plot showing the proportion of each type of event will be
        shown after filtering.
    fig_path : str
        Path where the figure will be saved. If None, figure is displayed
        interactively.
    prefix : str
        If the library has a name, it will be shown on plots.
    """
    n_uncuts = 0
    n_loops = 0
    n_weirds = 0
    lrange_intra = 0
    lrange_inter = 0

    # open the files for reading and writing
    with open(in_dat, "r") as pairs, open(out_filtered, "w") as filtered:
        for line in pairs:  # iterate over each line
            # Copy header lines to output
            if line.startswith("#"):
                filtered.write(line)
                continue

            p = process_read_pair(line)
            line_to_write = (
                "\t".join(
                    map(
                        str,
                        (
                            p["readID"],
                            p["chr1"],
                            p["pos1"],
                            p["chr2"],
                            p["pos2"],
                            p["strand1"],
                            p["strand2"],
                            p["frag1"],
                            p["frag2"],
                        ),
                    )
                )
                + "\n"
            )
            if p["chr1"] == p["chr2"]:
                # Do not report ++ and -- pairs on the same fragment (impossible)
                if p["frag1"] == p["frag2"] and p["strand1"] == p["strand2"]:
                    n_weirds += 1
                elif p["nsites"] <= thr_loop and p["type"] == "-+":
                    n_loops += 1
                elif p["nsites"] <= thr_uncut and p["type"] == "+-":
                    n_uncuts += 1
                else:
                    lrange_intra += 1
                    filtered.write(line_to_write)

            if p["chr1"] != p["chr2"]:
                lrange_inter += 1
                filtered.write(line_to_write)

    if lrange_inter > 0:
        ratio_inter = round(
            100 * lrange_inter / float(lrange_intra + lrange_inter), 2
        )
    else:
        ratio_inter = 0

    # Log quick summary of operation results
    kept = lrange_intra + lrange_inter
    discarded = n_loops + n_uncuts + n_weirds
    total = kept + discarded
    logger.info(
        "Proportion of inter contacts: {0}% (intra: {1}, "
        "inter: {2})".format(ratio_inter, lrange_intra, lrange_inter)
    )
    logger.info(
        "{0} pairs discarded: Loops: {1}, Uncuts: {2}, Weirds: {3}".format(
            discarded, n_loops, n_uncuts, n_weirds
        )
    )
    logger.info(
        "{0} pairs kept ({1}%)".format(
            kept, round(100 * kept / (kept + discarded), 2)
        )
    )

    # Visualize summary if requested by user
    if plot_events:
        try:
        # Plot: make a square figure and axes to plot a pieChart:
            plt.figure(2, figsize=(6, 6))
            # The slices will be ordered and plotted counter-clockwise.
            fracs = [n_uncuts, n_loops, n_weirds, lrange_intra, lrange_inter]
            # Format labels to include event names and proportion
            labels = list(
                map(
                    lambda x: (x[0] + ": %.2f%%") % (100 * x[1] / total),
                    [
                        ("Uncuts", n_uncuts),
                        ("Loops", n_loops),
                        ("Weirds", n_weirds),
                        ("3D intra", lrange_intra),
                        ("3D inter", lrange_inter),
                    ],
                )
            )
            colors = ["salmon", "lightskyblue", "yellow", "palegreen", "plum"]
            patches, _ = plt.pie(fracs, colors=colors, startangle=90)
            plt.legend(
                patches,
                labels,
                loc='upper left',
                bbox_to_anchor=(-0.1, 1.),
            )
            if prefix:
                plt.title(
                    "Distribution of library events in {}".format(prefix),
                    bbox={"facecolor": "1.0", "pad": 5},
                )
            plt.text(
                0.3,
                1.15,
                "Threshold Uncuts = " + str(thr_uncut),
                fontdict=None,
            )
            plt.text(
                0.3,
                1.05,
                "Threshold Loops = " + str(thr_loop),
                fontdict=None,
            )

            plt.text(
                -1.5,
                -1.2,
                "Total number of reads = " + str(total),
                fontdict=None,
            )
            plt.text(
                -1.5,
                -1.3,
                "Ratio inter/(intra+inter) = " + str(ratio_inter) + "%",
                fontdict=None,
            )
            percentage = round(
                100
                * float(lrange_inter + lrange_intra)
                / (n_loops + n_uncuts + n_weirds + lrange_inter + lrange_intra)
            )
            plt.text(
                -1.5,
                -1.4,
                "Selected reads = {0}%".format(percentage),
                fontdict=None,
            )
            if fig_path:
                plt.savefig(fig_path)
            else:
                plt.show()
            plt.clf()
        except Exception:
            logger.error(
                "Unable to show plots. Perhaps there is no Xserver running ?"
                "(might be due to windows environment) skipping figure "
                "generation."
            )
