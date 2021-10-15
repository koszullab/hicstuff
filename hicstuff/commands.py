#!/usr/bin/env python3
# coding: utf-8

"""Abstract command classes for hicstuff

This module contains all classes related to hicstuff
commands:

    -iteralign (iterative mapping)
    -digest (genome chunking)
    -cutsite (preprocess fastq by cutting reads into digestion products)
    -filter (Hi-C 'event' sorting: loops, uncuts, weird
     and 'true contacts')
    -view (map visualization)
    -pipeline (whole contact map generation)
    -distancelaw (Analysis tool and plot for the distance law)

Running 'pipeline' implies running 'digest', but not
iteralign or filter unless specified, because they can
take up a lot of time for dimnishing returns.

Note
----
Structure based on Rémy Greinhofer (rgreinho) tutorial on subcommands in
docopt : https://github.com/rgreinho/docopt-subcommands-example
cmdoret, 20181412

Raises
------
NotImplementedError
    Will be raised if AbstractCommand is called for
    some reason instead of one of its children.
ValueError
    Will be raised if an incorrect chunking method (e.g.
    not an enzyme or number or invalid range view is
    specified.
"""
import re
import sys, os, shutil
import tempfile
from os.path import join, dirname
from matplotlib import pyplot as plt
from matplotlib import cm
from docopt import docopt
import pandas as pd
import numpy as np
import pysam as ps
import glob
import copy
from Bio import SeqIO
import hicstuff.view as hcv
import hicstuff.hicstuff as hcs
import hicstuff.cutsite as hcc
import hicstuff.digest as hcd
import hicstuff.iteralign as hci
import hicstuff.filter as hcf
from hicstuff.version import __version__
import hicstuff.io as hio
from hicstuff.log import logger
import hicstuff.pipeline as hpi
import hicstuff.distance_law as hcdl

DIVERGENT_CMAPS = [
    "PiYG",
    "PRGn",
    "BrBG",
    "PuOr",
    "RdGy",
    "RdBu",
    "RdYlBu",
    "RdYlGn",
    "Spectral",
    "coolwarm",
    "bwr",
    "seismic",
]


class AbstractCommand:
    """Abstract base command class

    Base class for the commands from which
    other hicstuff commands derive.
    """

    def __init__(self, command_args, global_args):
        """Initialize the commands"""
        self.args = docopt(self.__doc__, argv=command_args)
        self.global_args = global_args
        # Map Hi-C format to file extension
        self.fmt2ext = {"cool": ".cool", "bg2": ".bg2", "graal": ".mat.tsv"}

    def execute(self):
        """Execute the commands"""
        raise NotImplementedError

    def check_output_path(self, path, force=False):
        """Throws error if the output file exists. Create required file tree otherwise."""
        # Get complete output filename and prevent overwriting unless force is enabled
        if not force and os.path.exists(path):
            raise IOError(
                "Output file already exists. Use --force to overwrite"
            )
        if dirname(path):
            os.makedirs(dirname(path), exist_ok=True)


class Iteralign(AbstractCommand):
    """Iterative mapping command

    Truncate reads from a fastq file to 20 basepairs and iteratively extend and
    re-align the unmapped reads to optimize the proportion of uniquely aligned
    reads in a 3C library.

    usage:
        iteralign [--aligner=bowtie2] [--threads=1] [--min-len=20] [--read-len=INT]
                  [--tempdir=DIR] --out-bam=FILE --genome=FILE <reads.fq>

    arguments:
        reads.fq                Fastq file containing the reads to be aligned

    options:
        -g, --genome=FILE        The genome on which to map the reads. Must be
                                 the path to the bowtie2/bwa index if using bowtie2/bwa
                                 or to the genome in fasta format if using minimap2.
        -t, --threads=INT        Number of parallel threads allocated for the
                                 alignment [default: 1].
        -T, --tempdir=DIR        Temporary directory. Defaults to current
                                 directory.
        -a, --aligner=bowtie2    Choose alignment software between bowtie2,
                                 minimap2 or bwa. minimap2 should only be used for
                                 reads > 100 bp. [default: bowtie2]
        -l, --min-len=INT        Length to which the reads should be
                                 truncated [default: 20].
        -o, --out-bam=FILE       Path where the alignment will be written in
                                 BAM format.
        -R, --read-len=INT       Read length in input FASTQ file. If not provided,
                                 this is estimated from the first read in the file.
    """

    def execute(self):
        read_len = self.args["--read-len"]

        if read_len is not None:
            read_len = int(read_len)

        if not self.args["--tempdir"]:
            self.args["--tempdir"] = "."

        temp_directory = hio.generate_temp_dir(self.args["--tempdir"])

        hci.iterative_align(
            self.args["<reads.fq>"],
            temp_directory,
            self.args["--genome"],
            self.args["--threads"],
            self.args["--out-bam"],
            aligner=self.args["--aligner"],
            min_len=int(self.args["--min-len"]),
            read_len=read_len,
        )
        # Deletes the temporary folder
        shutil.rmtree(temp_directory)


class Digest(AbstractCommand):
    """Genome chunking command

    Digests a fasta file into fragments based on a restriction enzyme or a
    fixed chunk size. Generates two output files into the target directory
    named "info_contigs.txt" and "fragments_list.txt"

    usage:
        digest [--plot] [--figdir=FILE] [--force] [--circular] [--size=0]
               [--outdir=DIR] --enzyme=ENZ <fasta>

    arguments:
        fasta                     Fasta file to be digested

    options:
        -c, --circular                  Specify if the genome is circular.
        -e, --enzyme=ENZ[,ENZ2,...]     A restriction enzyme or an integer
                                        representing fixed chunk sizes (in bp).
                                        Multiple comma-separated enzymes can
                                        be given.
        -F, --force                     Write even if the output file already exists.
        -s, --size=INT                  Minimum size threshold to keep
                                        fragments. [default: 0]
        -o, --outdir=DIR                Directory where the fragments and
                                        contigs files will be written.
                                        Defaults to current directory.
        -p, --plot                      Show a histogram of fragment length
                                        distribution after digestion.
        -f, --figdir=FILE               Path to directory of the output figure.
                                        By default, the figure is only shown
                                        but not saved.

    output:
        fragments_list.txt: information about restriction fragments (or chunks)
        info_contigs.txt: information about contigs or chromosomes

    """

    def execute(self):
        # If circular is not specified, change it from None to False
        if not self.args["--circular"]:
            self.args["--circular"] = False
        if not self.args["--outdir"]:
            self.args["--outdir"] = os.getcwd()
        # Create output directory if it does not exist
        if os.path.exists(self.args["--outdir"]):
            if not self.args["--force"]:
                raise IOError(
                    "Output directory already exists. Use --force to overwrite"
                )
        else:
            os.makedirs(self.args["--outdir"], exist_ok=True)
        if self.args["--figdir"]:
            figpath = join(self.args["--figdir"], "frags_hist.pdf")
        else:
            figpath = None
        # Split into a list if multiple enzymes given
        enzyme = self.args["--enzyme"]
        if re.search(r",", enzyme):
            enzyme = enzyme.split(",")

        hcd.write_frag_info(
            self.args["<fasta>"],
            enzyme,
            self.args["--size"],
            output_dir=self.args["--outdir"],
            circular=self.args["--circular"],
        )

        hcd.frag_len(
            output_dir=self.args["--outdir"],
            plot=self.args["--plot"],
            fig_path=figpath,
        )


class Cutsite(AbstractCommand):
    """Cutsite command

    Generates new gzipped fastq files from original fastq. The function will cut
    the reads at their religation sites and creates new pairs of reads with the
    different fragments obtained after cutting at the digestion sites.

    There are three choices to how combine the fragments. 1. "for_vs_rev": All
    the combinations are made between one forward fragment and one reverse
    fragment. 2. "all": All 2-combinations are made. 3. "pile": Only
    combinations between adjacent fragments in the initial reads are made.

    usage:
        cutsite --forward=FILE --reverse=FILE  --prefix=STR --enzyme=STR
        [--mode=for_vs_rev] [--seed-size=20] [--threads=1]

    options:
        -1, --forward=FILE      Fastq file containing the forward reads to 
                                digest.
        -2, --reverse=FILE      Fastq file containing the reverse reads to 
                                digest.
        -p, --prefix=STR        Prefix of the path where to write the digested 
                                gzipped fastq files. Filenames will be added the
                                suffix "_{1,2}.fq.gz".
        -e, --enzyme=STR        The list of restriction enzyme used to digest
                                the genome separated by a comma. Example:
                                HpaII,MluCI.
        -m, --mode=STR          Digestion mode. There are three possibilities:
                                "for_vs_rev", "all" and "pile". The first one
                                "for_vs_rev" makes all possible contact between
                                fragments from forward read versus the fragments
                                of the reverse reads. The second one "all"
                                consist two make all pairs of fragments 
                                possible. The third one "pile" will make the
                                contacts only with the adjacent fragments.
                                [Default: for_vs_rev]
        -s, --seed-size=INT     Minimum size of a read. (i.e. seed size used 
                                in mapping as reads smaller won't be mapped.)
                                [Default: 20]
        -t, --threads=INT       Number of parallel threads allocated for the
                                alignement. [Default: 1]
    """

    def execute(self):

        # Check for mandatory options
        for option in ["--prefix", "--forward", "--reverse"]:
            if self.args[option] is None:
                raise ValueError(f"{option} is mandatory.")

        prefix = self.args["--prefix"]
        # Create output directory if it does not exist
        if dirname(prefix):
            os.makedirs(dirname(prefix), exist_ok=True)
        output_for = prefix + "_1.fq.gz"
        output_rev = prefix + "_2.fq.gz"
        # Digestion of the reads.
        logger.info("Digestion of the reads:")
        logger.info("Enzyme used: {0}".format(self.args["--enzyme"]))
        logger.info(
            "Mode used to cut the reads: {0}".format(self.args["--mode"])
        )
        hcc.cut_ligation_sites(
            self.args["--forward"],
            self.args["--reverse"],
            output_for,
            output_rev,
            enzyme=self.args["--enzyme"],
            mode=self.args["--mode"],
            seed_size=int(self.args["--seed-size"]),
            n_cpu=int(self.args["--threads"]),
        )


class Filter(AbstractCommand):
    """Mapping event filtering command

    Filters spurious 3C events such as loops and uncuts from the library based
    on a minimum distance threshold automatically estimated from the library by
    default. Can also plot 3C library statistics.

    usage:
        filter [--interactive | --thresholds INT-INT] [--plot]
               [--figdir FILE] [--prefix STR] <input> <output>

    arguments:
        input       2D BED file containing coordinates of Hi-C interacting
                    pairs, the index of their restriction fragment and their
                    strands.
        output      Path to the filtered file, in the same format as the input.

    options:
        -f, --figdir=DIR                  Path to the output figure directory.
                                          By default, the figure is only shown
                                          but not saved.
        -i, --interactive                 Interactively shows plots and asks
                                          for thresholds.
        -p, --plot                        Shows plots of library composition
                                          and 3C events abundance.
        -P, --prefix STR                  If the library has a name, it will
                                          be shown on the figures.
        -t, --thresholds=INT-INT          Manually defines integer values for
                                          the thresholds in the order
                                          [uncut, loop]. Reads above those values
                                          are kept.
    """

    def execute(self):
        figpath = None
        if self.args["--thresholds"]:
            # Thresholds supplied by user beforehand
            uncut_thr, loop_thr = self.args["--thresholds"].split("-")
            try:
                uncut_thr = int(uncut_thr)
                loop_thr = int(loop_thr)
            except ValueError:
                logger.error(
                    "You must provide integer numbers for the thresholds."
                )
        else:
            # Threshold defined at runtime
            if self.args["--figdir"]:
                figpath = join(self.args["--figdir"], "event_distance.pdf")
                if not os.path.exists(self.args["--figdir"]):
                    os.makedirs(self.args["--figdir"])
            uncut_thr, loop_thr = hcf.get_thresholds(
                self.args["<input>"],
                interactive=self.args["--interactive"],
                plot_events=self.args["--plot"],
                fig_path=figpath,
                prefix=self.args["--prefix"],
            )
        # Filter library and write to output file
        figpath = None
        if self.args["--figdir"]:
            figpath = join(self.args["--figdir"], "event_distribution.pdf")

        hcf.filter_events(
            self.args["<input>"],
            self.args["<output>"],
            uncut_thr,
            loop_thr,
            plot_events=self.args["--plot"],
            fig_path=figpath,
            prefix=self.args["--prefix"],
        )


class View(AbstractCommand):
    """Contact map visualization command

    Visualize a Hi-C matrix file as a heatmap of contact frequencies. Allows to
    tune visualisation by binning and normalizing the matrix, and to save the
    output image to disk. If no output is specified, the output is displayed.

    usage:
        view [--binning=1] [--despeckle] [--frags FILE] [--trim INT] [--n-mad 3.0] [--lines]
             [--normalize] [--min=0] [--max=99%] [--output=IMG] [--cmap=Reds] [--dpi=300]
             [--transform=STR] [--circular] [--region=STR] <contact_map> [<contact_map2>]

    arguments:
        contact_map             Sparse contact matrix in bg2, cool or graal format
        contact_map2            Sparse contact matrix in bg2, cool or graal format,
                                if given, the log ratio of contact_map/contact_map2
                                will be shown.


    options:
        -b, --binning=INT[bp|kb|Mb|Gb]   Rebin the matrix. If no unit is given, bins will
                                         be merged by groups of INT. If a unit is given,
                                         bins of that size will be generated. [default: 1]
        -c, --cmap=STR                   The name of a matplotlib colormap to
                                         use for the matrix. [default: Reds]
        -C, --circular                   Use if the genome is circular.
        -d, --despeckle                  Remove sharp increases in long range
                                         contact by averaging surrounding
                                         values.
        -D, --dpi=INT                    Map resolution in DPI (dots per inch). [default: 300]
        -f, --frags=FILE                 Required for bp binning and chromosome lines.
                                         Tab-separated file with headers, containing
                                         fragments start position in the 3rd
                                         column, as generated by hicstuff
                                         pipeline.
        -T, --transform=STR              Apply a mathematical transformation to pixel values
                                         to improve visibility of long range signals. Possible
                                         values are: log2, log10, ln, sqrt, exp0.2.
        -l, --lines                      Add dotted lines marking separation between chromosomes
                                         or contigs. Requires --frags.
        -M, --max=INT                    Saturation threshold. Maximum pixel
                                         value is set to this number. Can be
                                         followed by % to use a percentile of
                                         nonzero pixels in the contact
                                         map. [default: 99%]
        -m, --min=INT                    Minimum of the colorscale, works
                                         identically to --max. [default: 0]
        -N, --n-mad=INT                  Number of median absolute deviations (MAD) from the median
                                         of log bin sums allowed to keep bins in the normalization
                                         procedure [default: 3.0].
        -n, --normalize                  Should ICE normalization be performed
                                         before rendering the matrix ?
        -o, --output=FILE                Name of the image file where the view is stored.
        -r, --region=STR[;STR]           Only view a region of the contact map.
                                         Regions are specified as UCSC strings.
                                         (e.g.:chr1:1000-12000). If only one
                                         region is given, it is viewed on the
                                         diagonal. If two regions are given,
                                         The contacts between both are shown.
        -t, --trim=INT                   Trims outlier rows/columns from the
                                         matrix if the sum of their contacts
                                         deviates from the mean by more than
                                         INT standard deviations.
    """

    def data_transform(self, dense_map, operation="log10"):
        """
        Apply a mathematical operation on a dense Hi-C map. Valid
        operations are: log2, log10, ln, sqrt, exp0.2
        """
        ops = {
            "log10": np.log10,
            "log2": np.log2,
            "ln": np.log,
            "sqrt": np.sqrt,
        }
        if operation in ops:
            return ops[operation](dense_map)
        elif re.match(r"exp", operation):
            splitop = operation.split("exp")
            exp_val = float(splitop[1])
            return dense_map ** exp_val
        elif hasattr(np, operation) and callable(np.__dict__[operation]):
            logger.warning("Using built-in numpy callable: %s", operation)
            return np.__dict__[operation](dense_map)
        else:
            raise TypeError("Supplied transform function is not supported.")

    def process_matrix(self, sparse_map):
        """
        Performs any combination of binning, normalisation, log transformation,
        trimming and subsetting based on the attributes of the instance class.
        """
        # BINNING
        if self.binning > 1:
            if self.bp_unit:
                self.pos = self.frags.iloc[:, 2]
                binned_map, binned_pos = hcs.bin_bp_sparse(
                    M=sparse_map, positions=self.pos, bin_len=self.binning
                )
                # Get bin numbers of chromosome starts
                binned_start = np.append(
                    np.where(binned_pos == 0)[0], len(binned_pos)
                )
                # Get bin length of each chromosome
                num_binned = binned_start[1:] - binned_start[:-1]
                # Get unique chromosome names without losing original order
                # (numpy.unique sorts output)
                chr_names_idx = np.unique(
                    self.frags.iloc[:, 1], return_index=True
                )[1]
                chr_names = [
                    self.frags.iloc[index, 1]
                    for index in sorted(chr_names_idx)
                ]
                binned_chrom = np.repeat(chr_names, num_binned)
                binned_frags = pd.DataFrame(
                    {"chrom": binned_chrom, "start_pos": binned_pos[:, 0]}
                )
                binned_frags["end_pos"] = binned_frags.groupby("chrom")[
                    "start_pos"
                ].shift(-1)
                chrom_ends = self.frags.groupby("chrom").end_pos.max()
                # Fill ends of chromosome bins with actual chromosome length
                for cn in chrom_ends.index:
                    binned_frags.end_pos[
                        np.isnan(binned_frags.end_pos)
                        & (binned_frags.chrom == cn)
                    ] = chrom_ends[cn]

            else:
                # Note this is a basic binning procedure, chromosomes are
                # not taken into account -> last few fragments of a chrom
                # are merged with the first few of the next
                binned_map = hcs.bin_sparse(
                    M=sparse_map, subsampling_factor=self.binning
                ) 
                if self.frags:
                    binned_frags = self.frags.iloc[:: self.binning, :]
                    binned_frags = binned_frags.reset_index(drop=True)
                    # Since matrix binning ignores chromosomes, we
                    # have to do the same procedure with fragments
                    # we just correct the coordinates to start at 0
                    def shift_min(x):
                        try:
                            x[x == min(x)] = 0
                        except ValueError:
                            pass
                        return x

                    binned_frags.start_pos = binned_frags.groupby(
                        "chrom", sort=False
                    ).start_pos.apply(shift_min)
                else:
                    binned_frags = self.frags

        else:
            binned_map = sparse_map
            binned_frags = self.frags

        # TRIMMING
        if self.args["--trim"]:
            try:
                trim_std = float(self.args["--trim"])
            except ValueError:
                logger.error(
                    "You must specify a number of standard deviations for "
                    "trimming"
                )
                raise
            binned_map = hcs.trim_sparse(binned_map, n_mad=trim_std)

        # NORMALIZATION
        if self.args["--normalize"]:
            binned_map = hcs.normalize_sparse(
                binned_map, norm="ICE", n_mad=float(self.args["--n-mad"])
            )

        # ZOOM REGION
        if self.args["--region"]:
            if self.args["--lines"]:
                raise NotImplementedError(
                    "Chromosome lines are currently incompatible with a region zoom"
                )
            if self.frags is None:
                logger.error(
                    "A fragment file must be provided to subset "
                    "genomic regions. See hicstuff view --help"
                )
                sys.exit(1)
            # Load chromosomes and positions from fragments list
            reg_pos = binned_frags[['chrom', 'start_pos']]
            region = self.args["--region"]
            if ";" in region:
                # 2 input regions: zoom anywhere in matrix
                self.symmetric = False
                reg1, reg2 = region.split(";")
                reg1 = parse_ucsc(reg1, reg_pos)
                reg2 = parse_ucsc(reg2, reg_pos)
            else:
                # Only 1 input region: zoom on diagonal
                region = parse_ucsc(region, reg_pos)
                reg1 = reg2 = region
            binned_map = binned_map.tocsr()
            binned_map = binned_map[reg1[0] : reg1[1], reg2[0] : reg2[1]]
            binned_map = binned_map.tocoo()

        return binned_map, binned_frags

    def execute(self):

        input_map = self.args["<contact_map>"]
        hic_fmt = hio.get_hic_format(input_map)
        cmap = self.args["--cmap"]
        # Switch to a divergent colormap for plotting ratios
        if (
            self.args["<contact_map2>"] is not None
            and cmap not in DIVERGENT_CMAPS
        ):
            # In case user specified a custom cmap incompatible with ratios
            if cmap != "Reds":
                logger.warning(
                    "You chose a non-divergent colormap. Valid divergent "
                    "cmaps are:\n\t{}".format(" ".join(DIVERGENT_CMAPS))
                )
            logger.info(
                "Defaulting to seismic colormap for ratios. You can pick "
                "another divergent colormap if you wish."
            )
            cmap = "seismic"
        self.bp_unit = False
        bin_str = self.args["--binning"].upper()
        self.symmetric = True
        transform = self.args["--transform"]
        try:
            # Subsample binning
            self.binning = int(bin_str)
        except ValueError:
            if re.match(r"^[0-9]+[KMG]?B[P]?$", bin_str):
                if hic_fmt == "graal" and not self.args["--frags"]:
                    logger.error(
                        "A fragment file must be provided to perform "
                        "basepair binning. See hicstuff view --help"
                    )
                    sys.exit(1)
                # Load positions from fragments list
                self.binning = parse_bin_str(bin_str)
                self.bp_unit = True
            else:
                logger.error(
                    "Please provide an integer or basepair value for binning."
                )
                raise
        sparse_map, self.frags, _ = hio.flexible_hic_loader(
            input_map, fragments_file=self.args["--frags"], quiet=True
        )
        output_file = self.args["--output"]
        processed_map, frags = self.process_matrix(sparse_map)
        # If 2 matrices given compute log ratio
        if self.args["<contact_map2>"]:
            sparse_map2, _, _ = hio.flexible_hic_loader(
                self.args["<contact_map2>"],
                fragments_file=self.args["--frags"],
                quiet=True,
            )
            processed_map2, _ = self.process_matrix(sparse_map2)
            if sparse_map2.shape != sparse_map.shape:
                logger.error(
                    "You cannot compute the ratio of matrices with "
                    "different dimensions"
                )
            # Get log of values for both maps
            processed_map.data = np.log2(processed_map.data)
            processed_map2.data = np.log2(processed_map2.data)
            # Note: Taking diff of logs instead of log of ratio because sparse
            # mat division yields dense matrix in current implementation.
            # Changing base to 2 afterwards.
            processed_map = processed_map.tocsr() - processed_map2.tocsr()
            processed_map = processed_map.tocoo()
            processed_map.data[np.isnan(processed_map.data)] = 0.0
            # Log transformation done already
            transform = False

        if self.args["--despeckle"]:
            processed_map = hcs.despeckle_simple(processed_map)
        try:
            if self.symmetric:
                dense_map = hcv.sparse_to_dense(
                    processed_map, remove_diag=False
                )
            else:
                dense_map = processed_map.toarray()

            def set_v(v, mat):
                if "%" in v:
                    try:
                        valid_pixels = (mat > 0) & (np.isfinite(mat))
                        val = np.percentile(
                            mat[valid_pixels], float(v.strip("%"))
                        )
                    # No nonzero / finite value
                    except IndexError:
                        val = 0
                else:
                    val = float(v)
                return val

            dense_map = dense_map.astype(float)
            self.vmax = set_v(self.args["--max"], dense_map)
            self.vmin = set_v(self.args["--min"], dense_map)
            if self.args["<contact_map2>"]:
                self.vmin, self.vmax = -2, 2
            # Log transform the map and the colorscale limits if needed
            if transform:
                dense_map = self.data_transform(dense_map, transform)
                # self.vmin = np.percentile(dense_map[np.isfinite(dense_map)], 1)
                # self.vmax = self.data_transform(self.vmax, transform)
                self.vmax = set_v(self.args["--max"], dense_map)
                self.vmin = set_v(self.args["--min"], dense_map)
            else:
                # Set 0 values in matrix to NA
                dense_map[dense_map == 0] = np.inf
            # Get chromosome coordinates if required
            if self.args["--lines"]:
                chrom_starts = np.where(np.diff(frags.start_pos) < 0)[0] + 1
            else:
                chrom_starts = None
            # Display NA values in white
            current_cmap = cm.get_cmap().copy()
            current_cmap.set_bad(color=current_cmap(0))

            hcv.plot_matrix(
                dense_map,
                filename=output_file,
                vmin=self.vmin,
                vmax=self.vmax,
                dpi=int(self.args["--dpi"]),
                cmap=cmap,
                chrom_starts=chrom_starts,
            )
        except MemoryError:
            logger.error("contact map is too large to load, try binning more")


class Pipeline(AbstractCommand):
    """Whole (end-to-end) contact map generation command

    Entire Pipeline to process fastq files into a Hi-C matrix. Uses all the
    individual components of hicstuff.

    usage:
        pipeline [--aligner=bowtie2] [--centromeres=FILE] [--circular] [--distance-law]
                 [--duplicates] [--enzyme=5000] [--filter] [--force] [--mapping=normal]
                 [--matfmt=graal] [--no-cleanup] [--outdir=DIR] [--plot] [--prefix=PREFIX]
                 [--quality-min=30] [--read-len=INT] [--remove-centromeres=0] [--size=0]
                 [--start-stage=fastq] [--threads=1] [--tmpdir=DIR] --genome=FILE <input1> [<input2>]

    arguments:
        input1:             Forward fastq file, if start_stage is "fastq", sam
                            file for aligned forward reads if start_stage is
                            "bam", or a .pairs file if start_stage is "pairs".
        input2:             Reverse fastq file, if start_stage is "fastq", sam
                            file for aligned reverse reads if start_stage is
                            "bam", or nothing if start_stage is "pairs".


    options:
        -a, --aligner=STR             Alignment software to use. Can be either
                                      bowtie2, minimap2 or bwa. minimap2 should
                                      only be used for reads > 100 bp.
                                      [default: bowtie2]
        -c, --centromeres=FILE        Positions of the centromeres separated by
                                      a space and in the same order than the
                                      chromosomes. Discordant with the circular
                                      option.
        -C, --circular                Enable if the genome is circular.
                                      Discordant with the centromeres option.
        -d, --distance-law            If enabled, generates a distance law file
                                      with the values of the probabilities to
                                      have a contact between two distances for
                                      each chromosomes or arms if the file with
                                      the positions has been given. The values
                                      are not normalized, or averaged.
        -D, --duplicates              Filter out PCR duplicates based on read
                                      positions.
        -e, --enzyme={STR|INT}        Restriction enzyme or "mnase" if a string,
                                      or chunk size (i.e. resolution) if a number.
                                      Can also be multiple comma-separated
                                      enzymes. [default: 5000]
        -f, --filter                  Filter out spurious 3C events (loops and
                                      uncuts) using hicstuff filter. Requires
                                      "-e" to be a restriction enzyme or mnase,
                                      not a chunk size. For more informations, see
                                      Cournac et al. BMC Genomics, 2012.
        -F, --force                   Write even if the output file already exists.
        -g, --genome=FILE             Reference genome to map against. Path to
                                      the bowtie2/bwa index if using bowtie2/bwa,
                                      or to a FASTA file if using minimap2.
        -m, --mapping=STR             normal|iterative|cutsite. Parameter of 
                                      mapping. "normal": Directly map reads 
                                      without any process. "iterative": Map 
                                      reads iteratively using iteralign, by 
                                      truncating reads to 20bp and then 
                                      repeatedly extending to align them.
                                      "cutsite": Cut reads at the religation 
                                      sites of the given enzyme using cutsite, 
                                      create new pairs of reads and then align 
                                      them ; enzyme is required [default: normal].
        -M, --matfmt=STR              The format of the output sparse matrix.
                                      Can be "bg2" for 2D Bedgraph format,
                                      "cool" for Mirnylab's cooler software, or
                                      "graal" for graal-compatible plain text
                                      COO format. [default: graal]
        -n, --no-cleanup              If enabled, intermediary BED files will
                                      be kept after generating the contact map.
                                      Disabled by defaut.
        -o, --outdir=DIR              Output directory. Defaults to the current
                                      directory.
        -p, --plot                    Generates plots in the output directory
                                      at different steps of the pipeline.
        -P, --prefix=STR              Overrides default filenames and prefixes all
                                      output files with a custom name.
        -q, --quality-min=INT         Minimum mapping quality for selecting
                                      contacts. [default: 30].
        -r, --remove-centromeres=INT  Integer. Number of kb that will be remove around
                                      the centromere position given by in the centromere
                                      file. [default: 0]
        -R, --read-len=INT            Maximum read length in the fastq file. Optionally
                                      used in iterative alignment mode. Estimated from
                                      the first read by default. Useful if input fastq
                                      is a composite of different read lengths.
        -s, --size=INT                Minimum size threshold to consider
                                      contigs. Keep all contigs by default.
                                      [default: 0]
        -S, --start-stage=STR         Define the starting point of the pipeline
                                      to skip some steps. Default is "fastq" to
                                      run from the start. Can also be "bam" to
                                      skip the alignment, "pairs" to start from a
                                      single pairs file or "pairs_idx" to skip
                                      fragment attribution and only build the
                                      matrix. [default: fastq]
        -t, --threads=INT             Number of threads to allocate.
                                      [default: 1].
        -T, --tmpdir=DIR              Directory for storing intermediary BED
                                      files and temporary sort files. Defaults
                                      to the output directory.

    output:
        abs_fragments_contacts_weighted.txt: the sparse contact map
        fragments_list.txt: information about restriction fragments (or chunks)
        info_contigs.txt: information about contigs or chromosomes
        hicstuff.log: details and statistics about the run.
    """

    def execute(self):

        if self.args["--filter"] and self.args["--enzyme"].isdigit():
            raise ValueError(
                "You cannot filter without specifying a restriction enzyme."
            )
        elif self.args["--enzyme"] in ("mnase", "dnase"):
            logger.info("## Enzyme provided is 'mnase', setting bin-size to 100bp")
            self.args["--enzyme"] = 100

        if not self.args["--outdir"]:
            self.args["--outdir"] = os.getcwd()

        if self.args["--matfmt"] not in ("graal", "bg2", "cool"):
            logger.error("matfmt must be either bg2, cool or graal.")
            raise ValueError

        read_len = self.args["--read-len"]
        if read_len is not None:
            read_len = int(read_len)

        hpi.full_pipeline(
            genome=self.args["--genome"],
            input1=self.args["<input1>"],
            input2=self.args["<input2>"],
            aligner=self.args["--aligner"],
            centromeres=self.args["--centromeres"],
            circular=self.args["--circular"],
            distance_law=self.args["--distance-law"],
            enzyme=self.args["--enzyme"],
            filter_events=self.args["--filter"],
            force=self.args["--force"],
            mapping=self.args["--mapping"],
            mat_fmt=self.args["--matfmt"],
            min_qual=int(self.args["--quality-min"]),
            min_size=int(self.args["--size"]),
            no_cleanup=self.args["--no-cleanup"],
            out_dir=self.args["--outdir"],
            pcr_duplicates=self.args["--duplicates"],
            plot=self.args["--plot"],
            prefix=self.args["--prefix"],
            read_len=read_len,
            remove_centros=self.args["--remove-centromeres"],
            start_stage=self.args["--start-stage"],
            threads=int(self.args["--threads"]),
            tmp_dir=self.args["--tmpdir"],
        )


class Scalogram(AbstractCommand):
    """
    Generate a scalogram.

    usage:
        scalogram [--cmap=viridis] [--centromeres=FILE] [--frags=FILE] [--range=INT-INT]
                  [--threads=1] [--output=FILE] [--normalize]
                  [--indices=INT-INT] [--despeckle] <contact_map>

    argument:
        <contact_map> The sparse Hi-C contact matrix.

    options:
        -C, --cmap=STR                     The matplotlib colormap to use for
                                           the plot. [default: viridis]
        -d, --despeckle                    Remove speckles (artifactual spots)
                                           from the matrix.
        -f, --frags=FILE                   Fragments_list.txt file providing mapping
                                           between genomic coordinates and bin IDs.
        -i, --indices=INT-INT              The range of bin numbers of the matrix to
                                           use for the plot. Can also be given in
                                           UCSC style genomic coordinates (requires -f).
                                           E.g. chr1:1Mb-10Mb.
        -o, --output=FILE                  Output file where the plot should be
                                           saved. Plot is only displayed by
                                           default.
        -n, --normalize                    Normalize the matrix first.
        -r, --range=INT-INT                The range of contact distance to look
                                           at. No limit by default. Values in
                                           basepairs by default but a unit can
                                           be specified (kb, Mb, ...).
        -t, --threads=INT                  Parallel processes to run in for
                                           despeckling. [default: 1]
    """

    def execute(self):
        mat, frags, _ = hio.flexible_hic_loader(
            self.args["<contact_map>"], fragments_file=self.args["--frags"]
        )
        if frags is not None:
            # If fragments_list.txt is provided, load chrom start and end columns
            frags = pd.read_csv(
                self.args["--frags"], delimiter="\t", usecols=(1, 2, 3)
            )
        if self.args["--range"]:
            shortest, longest = self.args["--range"].split("-")
            # If range given in number of bins
            try:
                shortest, longest = int(shortest), int(longest)
            # If range given in genomic scale
            except ValueError:
                shortest, longest = (
                    parse_bin_str(shortest),
                    parse_bin_str(longest),
                )
                # Use average bin size to convert genomic scale to number of bins
                avg_res = (frags.end_pos - frags.start_pos).mean()
                shortest, longest = (
                    int(shortest // avg_res),
                    int(longest // avg_res),
                )

        if self.args["--indices"]:
            start, end = self.args["--indices"].split("-")
            # If given in bin numbers
            try:
                start = int(start)
                end = int(end)
            # If given in genomic coordinates
            except ValueError:
                start, end = parse_ucsc(
                    self.args["--indices"],
                    frags.loc[:, ["chrom", "start_pos"]],
                )

        output_file = self.args["--output"]
        # good_bins = np.array(range(S.shape[0]))
        S = mat.tocsr()
        if not self.args["--range"]:
            shortest = 0
            longest = S.shape[0]

        if self.args["--normalize"]:
            # good_bins = np.where(hcs.get_good_bins(S, n_std=3) == 1)[0]
            S = hcs.normalize_sparse(S, norm="ICE")
            S = S.tocsr()
        if self.args["--despeckle"]:
            S = hcs.despeckle_simple(S, threads=int(self.args["--threads"]))

        # Cropping matrix before transforming to dense to reduce memory overhead
        # Note we leave a margin equal to longest range so that all windows can be computed
        if self.args["--indices"]:
            crop_inf, crop_sup = (
                max(0, start - longest),
                min(S.shape[0], end + longest),
            )
            crop_later = longest
            S = S[crop_inf:crop_sup, crop_inf:crop_sup]
        else:
            crop_later = 0

        D = hcv.sparse_to_dense(S)
        D = np.fliplr(np.rot90(hcs.scalogram(D), k=-1))
        # Crop the margin left previously to get actual indices on dimenstion 0
        # and focus scale to --range on dimension 1
        plt.contourf(
            D[crop_later : D.shape[1] - crop_later, shortest:longest],
            cmap=self.args["--cmap"],
        )
        if output_file:
            plt.savefig(output_file)
        else:
            plt.show()


class Rebin(AbstractCommand):
    """
    Rebins a Hi-C matrix and modifies its fragment and chrom files accordingly.
    Output files are in the same format as the input files (cool, graal or bg2).
    usage:
        rebin [--binning=1] [--frags=FILE] [--force] [--chroms=FILE] <contact_map> <out_prefix>

    arguments:
        contact_map             Sparse contact matrix in graal, cool or bg2 format.
        out_prefix              Prefix path (without extension) for the output files.

    options:
        -b, --binning=INT[bp|kb|Mb|Gb]   Subsampling factor or fix value in
                                         basepairs to use for binning
                                         [default: 1].
        -f, --frags=FILE                 Tab-separated file with headers,
                                         containing fragments start position in
                                         the 3rd column. This is the file
                                         "fragments_list.txt" generated by
                                         hicstuff pipeline. Required for graal
                                         matrices and recommended for bg2.
        -F, --force                      Write even if the output file already exists.
        -c, --chroms=FILE                Tab-separated with headers, containing
                                         chromosome names, size, number of
                                         restriction fragments. This is the file
                                         "info_contigs.txt" generated by hicstuff
                                         pipeline.
    """

    def execute(self):
        prefix = self.args["<out_prefix>"]
        bin_str = self.args["--binning"].upper()
        hic_fmt = hio.get_hic_format(self.args["<contact_map>"])
        # Get complete output filename and prevent overwriting unless --force is enabled
        out_name = prefix + self.fmt2ext[hic_fmt]
        self.check_output_path(out_name, force=self.args["--force"])
        # Load positions from fragments list and chromosomes from chrom file
        map_path = self.args["<contact_map>"]
        hic_map, frags, chromlist = hio.flexible_hic_loader(
            map_path,
            fragments_file=self.args["--frags"],
            chroms_file=self.args["--chroms"],
        )
        if hic_fmt == "graal" and (frags is None or chromlist is None):
            raise ValueError(
                "You must provide a chroms file and a fragments file "
                "when rebinning a matrix in graal format. (hint: the "
                "files info_contigs.txt and fragments_list.txt)"
            )
        # Create output directory if it does not exist
        if dirname(prefix):
            os.makedirs(dirname(prefix), exist_ok=True)
        bp_unit = False
        try:
            # Subsample binning
            binning = int(bin_str)
        except ValueError:
            # Basepair binning: determine bin size
            if re.match(r"^[0-9]+[KMG]?B[P]?$", bin_str):
                binning = parse_bin_str(bin_str)
                bp_unit = True
            else:
                logger.error(
                    "Please provide an integer or basepair value for binning."
                )
                raise
        chromnames = np.unique(frags.chrom)

        if bp_unit:
            # Basepair binning: Perform binning
            hic_map, _ = hcs.bin_bp_sparse(hic_map, frags.start_pos, binning)
            for chrom in chromnames:
                chrom_mask = frags.chrom == chrom
                # For all chromosomes, get new bin start positions
                bin_id = frags.loc[chrom_mask, "start_pos"] // binning
                frags.loc[chrom_mask, "id"] = bin_id + 1
                frags.loc[chrom_mask, "start_pos"] = binning * bin_id
                bin_ends = binning * bin_id + binning
                # Do not allow bin ends to be larger than chrom size
                try:
                    chromsize = chromlist.length[
                        chromlist.contig == chrom
                    ].values[0]
                except AttributeError:
                    chromsize = chromlist["length_kb"][
                        chromlist.contig == chrom
                    ].values[0]
                bin_ends[bin_ends > chromsize] = chromsize
                frags.loc[frags.chrom == chrom, "end_pos"] = bin_ends

            # Account for special cases where restriction fragments are larger than
            # bin size, resulting in missing bins (i.e. jumps in bin ids)
            id_diff = (
                np.array(frags.loc[:, "id"])[1:]
                - np.array(frags.loc[:, "id"])[:-1]
            )
            # Normal jump is 1, new chromosome (reset id) is < 0, abnormal is > 1
            # Get panda indices of abnormal jumps
            jump_frag_idx = np.where(id_diff > 1)[0]
            add_bins = id_diff - 1
            # Need to insert [jump] bins after indices with abnormal [jump]
            miss_bins = [None] * np.sum(add_bins[jump_frag_idx])
            miss_bin_id = 0
            for idx in jump_frag_idx:
                jump_size = add_bins[idx]
                for j in range(1, jump_size + 1):
                    # New generated bins will be given attributes based on the previous bin
                    # e.g. if 2 missing bins between bins 2 and 5:
                    # id[3] = id[2] + 1 * 1 and id[4] = id[2] + 1 * 2
                    miss_bins[miss_bin_id] = {
                        "id": frags.loc[idx, "id"] + 1 * j,
                        "chrom": frags.loc[idx, "chrom"],
                        "start_pos": frags.loc[idx, "start_pos"] + binning * j,
                        "end_pos": frags.loc[idx, "end_pos"] + binning * j,
                        "size": binning,
                        "gc_content": np.NaN,
                    }
                    miss_bin_id += 1
                    # Shift bins row idx to allow

            # Give existing bins spaced row idx to allow inserting missing bins
            idx_shift = copy.copy(id_diff)
            idx_shift[idx_shift < 1] = 1
            existing_bins_idx = np.cumsum(idx_shift)
            # Prepend first bin (lost when computing diff)
            existing_bins_idx = np.insert(existing_bins_idx, 0, 0)
            # Add missing bins to original table, and sort by idx
            # missing bins are "holes" in the continuous range of existing bins
            missing_bins_idx = sorted(
                set(range(existing_bins_idx[0], existing_bins_idx[-1]))
                - set(existing_bins_idx)
            )
            miss_bins_df = pd.DataFrame(
                miss_bins, columns=frags.columns, index=missing_bins_idx
            )
            frags["tmp_idx"] = existing_bins_idx
            miss_bins_df["tmp_idx"] = missing_bins_idx
            frags = pd.concat([frags, miss_bins_df], axis=0, sort=False)
            frags.sort_values("tmp_idx", axis=0, inplace=True)
            frags.drop("tmp_idx", axis=1, inplace=True)

        else:
            # Subsample binning
            hic_map = hcs.bin_sparse(hic_map, binning)
            # Use index for binning, but keep 1-indexed.
            # Exception when binning is 1 (no binning) where no need to shift
            shift_id = 0 if binning == 1 else 1
            frags.id = (frags.id // binning) + shift_id

        # Save original columns order
        col_ordered = list(frags.columns)
        # Get new start and end position for each bin
        frags = frags.groupby(["chrom", "id"], sort=False)
        positions = frags.agg({"start_pos": "min", "end_pos": "max"})
        positions.reset_index(inplace=True)
        # Compute mean for all added features in each index bin
        # Normally only other feature is GC content
        try:
            features = frags.agg("mean")
            features.reset_index(inplace=True)
            # set new bins positions
            frags = features
            frags["start_pos"] = 0
            frags["end_pos"] = 0
            frags.loc[:, positions.columns] = positions
        except pd.core.base.DataError:
            frags = positions
        frags["size"] = frags.end_pos - frags.start_pos
        cumul_bins = 0
        for chrom in chromnames:
            # Update cumulative length column in chromlist
            chrom_frags = frags.chrom == chrom
            n_bins = frags.start_pos[chrom_frags].shape[0]
            chromlist.loc[chromlist.contig == chrom, "n_frags"] = n_bins
            chromlist.loc[
                chromlist.contig == chrom, "cumul_length"
            ] = cumul_bins
            cumul_bins += n_bins
            # Adjust each chromosome's last bin end to match chromsize
            last_frag_end = frags.loc[chrom_frags, 'end_pos'].max()
            chromlen = chromlist.loc[chromlist.contig == chrom, 'length'].values[0]
            frags.loc[chrom_frags & (frags.end_pos == last_frag_end), 'end_pos'] = chromlen


        # Keep original column order
        frags = frags.reindex(columns=col_ordered)
        # Write 3 binned output files
        hio.flexible_hic_saver(
            hic_map,
            self.args["<out_prefix>"],
            frags=frags,
            chroms=chromlist,
            hic_fmt=hic_fmt,
        )


class Subsample(AbstractCommand):
    """
    Subsample contacts from a Hi-C matrix. Probability of sampling is proportional
    to the intensity of the bin.
    usage:
        subsample  [--prop=0.1] [--force] <contact_map> <subsampled_prefix>

    arguments:
        contact_map             Sparse contact matrix in graal, bg2 or cool format.
        subsampled_prefix       Path without extension to the output map in the same
                                format as the input containing only a fraction of the
                                contacts.

    options:
        -p, --prop=FLOAT        Proportion of contacts to sample from the input matrix
                                if between 0 and 1. Raw number of contacts to keep if
                                superior to 1. [default: 0.1]
        -F, --force             Write even if the output file already exists.
    """

    def execute(self):
        hic_fmt = hio.get_hic_format(self.args["<contact_map>"])
        prefix = self.args["<subsampled_prefix>"]
        # Get complete output filename and prevent overwriting unless --force is enabled
        out_name = prefix + self.fmt2ext[hic_fmt]
        self.check_output_path(out_name, force=self.args["--force"])
        mat, frags, _ = hio.flexible_hic_loader(
            self.args["<contact_map>"], quiet=True
        )
        subsampled = hcs.subsample_contacts(mat, float(self.args["--prop"]))
        subsampled = subsampled.tocoo()
        hio.flexible_hic_saver(
            subsampled,
            prefix,
            frags=frags,
            hic_fmt=hic_fmt,
            quiet=True,
        )


class Convert(AbstractCommand):
    """
    Convert between different Hi-C dataformats. Currently supports tsv (graal),
    bedgraph2D (bg2) and cooler (cool). Input format is automatically inferred.

    usage:
        convert [--frags=FILE] [--chroms=FILE] [--force] [--genome=FILE]
                [--to=cool] <contact_map> <prefix>

    arguments:
        contact_map               The file containing the contact frequencies.
        prefix                    The prefix path for output files. An extension
                                  will be added to the files depending on the
                                  output format.

    options:
        -f, --frags=FILE          Tab-separated file with headers,
                                  containing columns id, chrom, start_pos,
                                  end_pos size. This is the file
                                  "fragments_list.txt" generated by
                                  hicstuff pipeline. Required for graal
                                  matrices and recommended for bg2.
        -F, --force               Write even if the output file already exists.
        -g, --genome=FILE         Optional genome file used to add a GC content
                                  column to the fragments table. This is
                                  required to generate instagraal-compatible
                                  files.
        -c, --chroms=FILE         Tab-separated with headers, containing
                                  columns contig, length, n_frags, cumul_length.
                                  This is the file "info_contigs.txt" generated
                                  by hicstuff pipeline.
        -T, --to=STR              The format to which files should be
                                  converted. [default: cool]
    """

    def execute(self):
        out_fmt = self.args["--to"]
        mat_path = self.args["<contact_map>"]
        frags_path = self.args["--frags"]
        genome_path = self.args["--genome"]
        chroms_path = self.args["--chroms"]
        prefix = self.args["<prefix>"]

        # Get complete output filename and prevent overwriting unless --force is enabled
        out_name = prefix + self.fmt2ext[out_fmt]
        self.check_output_path(out_name, force=self.args["--force"])

        # Load
        mat, frags, chroms = hio.flexible_hic_loader(
            mat_path,
            fragments_file=frags_path,
            chroms_file=chroms_path,
            quiet=True,
        )

        # Modify fragments for instagraal compatibility
        # Add fragments size column
        chrom_col, start_col, end_col = hio.get_pos_cols(frags)
        size = frags[end_col] - frags[start_col]
        if "size" not in frags.columns:
            frags = frags.join(pd.DataFrame({"size": size}))
        # If genome was provided, add gc_content column
        if genome_path:
            gc = hio.gc_bins(genome_path, frags)
            frags = frags.join(pd.DataFrame({"gc_content": gc}))

        # Write
        mat = mat.astype(int)
        hio.flexible_hic_saver(
            mat=mat,
            out_prefix=prefix,
            frags=frags,
            chroms=chroms,
            hic_fmt=out_fmt,
        )


class Distancelaw(AbstractCommand):
    """Distance law tools.
    Take the distance law file from hicstuff and can average it, normalize it compute the
    slope of the curve and plot it.

    usage:
        distancelaw [--average] [--big-arm-only=INT] [--base=1.1] [--centromeres=FILE]
                    [--circular] [--frags=FILE] [--inf=INT] [--outputfile-img=FILE]
                    [--outputfile-tabl=FILE] [--labels=DIR] [--sup=INT]
                    [--remove-centromeres=0] (--pairs=FILE | --dist-tbl=FILE1[,FILE2,...])

    options:
        -a, --average                       If set, calculate the average of the distance
                                            law of the different chromosomes/arms in each
                                            condition. If two file given average is
                                            mandatory.
        -b, --big-arm-only=INT              Integer. It given will take only the arms bigger
                                            than the value given as argument.
        -B, --base=FLOAT                    Float corresponding of the base of the log to
                                            make the logspace which will slice the genomes
                                            in logbins. These slices will be in basepairs
                                            unit. [default is 1.1]
        -c, --centromeres=FILE              Positions of the centromeres separated by
                                            a space and in the same order as the
                                            chromosomes. This allows to plot chromosomal arms
                                            separately. Note this will only work with --pairs
                                            input, as the distance law needs to be recomputed.
                                            Incompatible with the circular option.
        -C, --circular                      Enable if the genome is circular. Discordant
                                            with the centromeres option.
        -d, --dist-tbl=FILE1[,FILE2,...]    Directory to the file or files containing the
                                            compute distance law. File should have the same
                                            format than the ones made by hicstuff pipeline.
        -f, --frags=FILE                    Tab-separated file with headers, containing
                                            columns id, chrom, start_pos, end_pos size.
                                            This is the file "fragments_list.txt" generated by
                                            hicstuff pipeline. Required if pairs are given.
        -i, --inf=INT                       Inferior born to plot the distance law. By
                                            default the value is 3000 bp (3 kb). Have to
                                            be strictly positive.
        -l, --labels=STR1,STR2...           List of string of the labels for the plot
                                            separated by a coma. If no labels given, give
                                            the names "Sample 1", "Sample 2"...
        -o, --outputfile-img=FILE           Output file. Format must be compatible with
                                            plt.savefig. Default : None.
        -O, --outputfile-tabl=TABLE         Output file. Default : None.
        -p, --pairs=FILE                    Pairs file. Format from 4D Nucleome Omics Data
                                            Standards Working Group with the 8th and 9th
                                            coulumns are the ID of the fragments of the
                                            reads 1 and 2. Only add if no distance_law table
                                            given. It will compute the table from these pairs
                                            and the fragments from the fragments file.
        -r, --remove-centromeres=INT        Integer. Number of kb that will be remove around
                                            the centromere position given by in the centromere
                                            file. [default: 0]
        -s, --sup=INT                       Superior born to plot the distance law. By
                                            default the value is the maximum length of all
                                            the dataset given. Also if big arm only set, it
                                            will be the minimum size of the arms/chromosomes
                                            taken to make the average.
    """

    def execute(self):
        # Give no file as output_file_img if no given.
        if self.args["--outputfile-img"]:
            output_file_img = self.args["--outputfile-img"]
        else:
            output_file_img = None
        # Compute the table of distance law if pairs given
        if self.args["--pairs"]:
            # Sanity check : frags mandatory if pairs given.
            if not self.args["--frags"] or self.args["--dist-tbl"]:
                logger.error(
                    "You have to give fragments and/or not give table of the distance law if pairs file given."
                )
                sys.exit(1)
            pairs = self.args["--pairs"]
            fragments = self.args["--frags"]
            # Give no file as output_file_tabl if no given.
            if self.args["--outputfile-tabl"]:
                output_file_tabl = self.args["--outputfile-tabl"]
            else:
                output_file_tabl = None
            # Check if centromeres file given and if the centromeres have to be removed
            if self.args["--centromeres"]:
                centromeres = self.args["--centromeres"]
                rm_centro = int(self.args["--remove-centromeres"])
            else:
                centromeres = None
                rm_centro = 0
            # Check if circular condition given
            if self.args["--circular"]:
                circular = self.args["--circular"]
            else:
                circular = None
            # Check logarithm base
            if self.args["--base"]:
                base = float(self.args["--base"])
            else:
                base = 1.1
            xs, ps, names = [None], [None], [None]
            xs[0], ps[0], names[0] = hcdl.get_distance_law(
                pairs_reads_file=pairs,
                fragments_file=fragments,
                centro_file=centromeres,
                base=base,
                out_file=output_file_tabl,
                circular=circular,
                rm_centro=rm_centro,
            )
            length_files = 1
        else:
            # Put in a list the path or the different paths given.
            distance_law_file = self.args["--dist-tbl"]
            distance_law_files = distance_law_file.split(",")
            length_files = len(distance_law_files)
            # Make new lists for the modified distance law.
            xs = [None] * length_files
            ps = [None] * length_files
            names = [None] * length_files
            # Iterate on the different file given by the user.
            for i in range(length_files):
                xs[i], ps[i], names[i] = hcdl.import_distance_law(
                    distance_law_files[i]
                )
            names = [name[0] for name in names]
        # Put the inf and sup according to the arguments given.
        if self.args["--inf"]:
            inf = int(self.args["--inf"])
        else:
            inf = 3000
        if self.args["--sup"]:
            sup = int(self.args["--sup"])
        else:
            sup = max(max(xs[0], key=len))
        # Add the option big army only.
        if self.args["--big-arm-only"]:
            big_arm_only = True
            arm_sup = int(self.args["--big-arm-only"])
        else:
            big_arm_only = False
            arm_sup = sup
        # Sanity check : Average mandatory if more than one file.
        if not self.args["--average"] and length_files > 1:
            logger.error("You have to average if more than one file.")
            sys.exit(1)
        # Iterate on the different file given by the user.
        for i in range(length_files):
            # Make the average if enabled
            if self.args["--average"]:
                xs[i], ps[i] = hcdl.average_distance_law(
                    xs[i], ps[i], arm_sup, big_arm_only
                )
                # If not average, we should to remove one level of list to have the good dimension.
        if not self.args["--average"]:
            names = names[0]
            xs = xs[0]
            ps = ps[0]
        # Normalize and make the derivative
        ps = hcdl.normalize_distance_law(xs, ps, inf, arm_sup)

        # Gave new names for the different samples.
        if self.args["--labels"]:
            labels = self.args["--labels"].split(",")
        else:
            if length_files == 1 and not self.args["--average"]:
                labels = []
                if len(names) > 1:
                    for chr_names in names:
                        labels.append(chr_names)
                else:
                    labels = names
            else:
                labels = []
                for i in range(length_files):
                    labels.append("Sample " + str(i))
                    
        # Make the plot if enabled, if not average plot the different arms or
        # chromosomes with the initial names else plot the different conditions
        # with the names labels.
        hcdl.plot_ps_slope(xs, ps, labels, output_file_img, inf, sup)

        # Export the new table if required.
        if self.args["--outputfile-tabl"]:
            hcdl.export_distance_law(
                xs,
                ps,
                labels, 
                self.args["--outputfile-tabl"]
            )


class Missview(AbstractCommand):
    """
    Previews bins that will be missing in a Hi-C map with a given read length by
    finding repetitive regions in the genome.

    usage:
        missview [--aligner=bowtie2] [--force] [--binning=5000]
                 [--threads=1] [--tmpdir=STR] --read-len=INT <genome> <output>

    arguments:
        genome               Genome file in fasta format.
        output               Path to the output image.

    options:
        -a, --aligner=STR    The read alignment software to use. Can be either
                             bowtie2, minimap2 or bwa. minimap2 should only be
                             used for reads > 100 bp. [default: bowtie2]
        -b, --binning=INT    Resolution to use to preview the Hi-C map. [default: 5000]
        -F, --force          Write even if the output file already exists.
        -R, --read-len=INT   Write even if the output file already exists.
        -t, --threads=INT    Number of CPUs to use in parallel. [default: 1]
        -T, --tmpdir=STR     Directory where temporary files will be generated.
    """

    def execute(self):
        aligner = self.args["--aligner"]
        force = self.args["--force"]
        genome = self.args["<genome>"]
        out = self.args["<output>"]
        resolution = parse_bin_str(self.args["--binning"])
        read_len = int(self.args["--read-len"])
        threads = int(self.args["--threads"])
        tmp_dir = self.args["--tmpdir"]
        if tmp_dir is None:
            tmp_dir = tempfile.TemporaryDirectory().name
        # Simulate reads and save into a fastq file
        phred = "F" * read_len
        tmp_fq = join(tmp_dir, "simulated_reads.fq")
        tmp_bam = join(tmp_dir, "simulated_reads.bam")
        self.check_output_path(tmp_fq, force=force)
        logger.info(
            "Simulating reads by splitting the genome into %i bp chunks",
            read_len,
        )
        with open(tmp_fq, "w") as fq_handle:
            for rec in SeqIO.parse(genome, "fasta"):
                for i in range(len(rec.seq) - read_len):
                    fq_handle.write("@NS_SIM_%s_%i\n" % (rec.id, i))
                    fq_handle.write(str(rec.seq[i : i + read_len]))
                    fq_handle.write("\n+\n" + phred + "\n")
        # Map reads to genome
        hpi.align_reads(
            tmp_fq,
            genome,
            tmp_bam,
            tmp_dir=tmp_dir,
            threads=threads,
            aligner=aligner,
        )
        # Sort alignments by name
        ps.sort(
            "-@",
            str(threads),
            "-n",
            "-O",
            "BAM",
            "-o",
            tmp_bam + ".sorted",
            tmp_bam,
        )
        shutil.move(tmp_bam + ".sorted", tmp_bam)
        # Run the standard pipeline with, using twice the forward reads.
        # This will generate a diagonal-only matrix
        hpi.full_pipeline(
            genome,
            tmp_bam,
            tmp_bam,
            start_stage="bam",
            aligner=aligner,
            enzyme=resolution,
            force=force,
            threads=threads,
            out_dir=tmp_dir,
            tmp_dir=tmp_dir,
        )
        # Plot and save matrix
        mat_path = join(tmp_dir, "abs_fragments_contacts_weighted.txt")
        mat = hio.load_sparse_matrix(mat_path)
        # Check which bins are not at the median (i.e. drop in mapping rate)
        log_content = open(glob.glob(join(tmp_dir, "*.log"))[0]).read()
        # Get (int rounded) percentage of reads mapped and convert to proportion
        prop_mapped = (
            int(re.search(r".*INFO :: ([0-9]*)% reads.*", log_content)[1])
            / 100
        )
        logger.info(
            "Bins with less than %s mapped reads will be considered undetectable",
            str(100 * prop_mapped) + "%",
        )
        unmappable = mat.diagonal(0) < prop_mapped * resolution
        mappable_mat = np.ones(mat.shape)
        mappable_mat[unmappable, :] = 0
        mappable_mat[:, unmappable] = 0
        hcv.plot_matrix(
            mappable_mat,
            filename=out,
            title=" %.3f%% missing bins for %s \nwith %i bp reads at resolution %i."
            % (
                round(100 * sum(unmappable) / len(unmappable), 3),
                os.path.basename(genome),
                read_len,
                resolution,
            ),
            dpi=600,
            vmax=2,
            cmap="Greys",
        )
        logger.info("Output image saved at %s.", out)


def parse_bin_str(bin_str):
    """Bin string parsing

    Take a basepair binning string as input and converts it into
    corresponding basepair values.

    Parameters
    ----------
    bin_str : str
        A basepair region (e.g. 150KB). Unit can be BP, KB, MB, GB.

    Example
    -------

        >>> parse_bin_str("150KB")
        150000
        >>> parse_bin_str("0.1mb")
        100000

    Returns
    -------
    binning : int
        The number of basepair corresponding to the binning string.
    """
    try:
        binning = int(bin_str)
    except ValueError:
        bin_str = bin_str.upper()
        binsuffix = {"B": 1, "K": 1000, "M": 1e6, "G": 1e9}
        unit_pos = re.search(r"[KMG]?B[P]?$", bin_str).start()
        bp_unit = bin_str[unit_pos:]
        # Extract unit and multiply accordingly for fixed bp binning
        binning = int(float(bin_str[:unit_pos]) * binsuffix[bp_unit[0]])

    return binning


def parse_ucsc(ucsc_str, bins):
    """
    Take a UCSC region in UCSC notation and a list of bin chromosomes and
    positions (in basepair) and converts it to range of bins.

    Parameters
    ----------
    ucsc_str : str
        The region string in UCSC notation (e.g. chr1:1000-2000)
    bins : pandas.DataFrame
        Dataframe of two columns containing the chromosome and start
        position of each bin. Each row must be one bin.

    Returns
    -------
    coord : tuple
        A tuple containing the bin range containing in the requested region.
    """
    if ":" in ucsc_str:
        chrom, bp = ucsc_str.split(":")
        bp = bp.replace(",", "").upper()
        start, end = bp.split("-")
        start, end = parse_bin_str(start), parse_bin_str(end)
        # Make absolute bin index (independent of chrom)
        bins["id"] = bins.index
        chrombins = bins.loc[bins.iloc[:, 0] == chrom, :]
        start = max([start, 1])
        start = max(chrombins.id[chrombins.iloc[:, 1] < start])
        end = max(chrombins.id[chrombins.iloc[:, 1] < end])
    else:
        chrom = ucsc_str
        # Make absolute bin index (independent of chrom)
        bins[2] = bins.index
        chrombins = bins.loc[bins.iloc[:, 0] == chrom, :]
        try:
            start = min(chrombins[2])
            end = max(chrombins[2])
        except ValueError:
            logger.error("Invalid chromosome")
            raise
    coord = (int(start), int(end))
    return coord
