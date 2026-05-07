"""hicstuff command-line interface."""

import copy
import glob
import os
import re
import shutil
import tempfile
from os.path import dirname, join

import numpy as np
import pandas as pd
import pysam as ps
import rich_click as click
from Bio import SeqIO
from matplotlib import cm
from matplotlib import pyplot as plt

import hicstuff.cutsite as hcc
import hicstuff.digest as hcd
import hicstuff.distance_law as hcdl
import hicstuff.filter as hcf
import hicstuff.hicstuff as hcs
import hicstuff.io as hio
import hicstuff.iteralign as hci
import hicstuff.pipeline as hpi
import hicstuff.stats as hcstats
import hicstuff.view as hcv
from hicstuff import __version__
from hicstuff.log import logger

# ---------------------------------------------------------------------------
# rich-click configuration
# ---------------------------------------------------------------------------
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.MAX_WIDTH = 100
click.rich_click.ERRORS_SUGGESTION = (
    "Run [bold cyan]hicstuff --help[/bold cyan] for usage information."
)
click.rich_click.COMMAND_GROUPS = {
    "hicstuff": [
        {
            "name": "Alignment",
            "commands": ["iteralign", "cutsite"],
        },
        {
            "name": "Processing",
            "commands": ["digest", "filter", "pipeline"],
        },
        {
            "name": "Matrix operations",
            "commands": ["rebin", "subsample", "convert"],
        },
        {
            "name": "Analysis & Visualization",
            "commands": ["view", "scalogram", "distancelaw", "missview", "stats"],
        },
    ]
}

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

# Map Hi-C format to file extension
_FMT2EXT = {"cool": ".cool", "bg2": ".bg2", "graal": ".mat.tsv"}


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _parse_bin_str(bin_str: str) -> int:
    """Convert a basepair binning string to an integer number of basepairs.

    Examples
    --------
    >>> _parse_bin_str("150KB")
    150000
    >>> _parse_bin_str("0.1mb")
    100000
    """
    try:
        return int(bin_str)
    except ValueError:
        bin_str = bin_str.upper().strip("P").strip("B")
    try:
        return int(bin_str)
    except ValueError:
        binsuffix = {"K": 1000, "M": 1e6, "G": 1e9}
        unit_pos = re.search(r"[KMG]?$", bin_str).start()
        bp_unit = bin_str[unit_pos:]
        return int(float(bin_str[:unit_pos]) * binsuffix[bp_unit[0]])


def _parse_ucsc(ucsc_str: str, bins: pd.DataFrame) -> tuple:
    """Convert a UCSC region string to a bin-index range.

    Parameters
    ----------
    ucsc_str : str
        Region in UCSC notation (e.g. ``chr1:1000-2000``).
    bins : pandas.DataFrame
        Two-column DataFrame: chromosome name and start position per bin.

    Returns
    -------
    tuple
        ``(start_bin, end_bin)`` integer bin indices.
    """
    if ":" in ucsc_str:
        chrom, bp = ucsc_str.split(":")
        bp = bp.replace(",", "").upper()
        start, end = bp.split("-")
        start, end = _parse_bin_str(start), _parse_bin_str(end)
        bins = bins.copy()
        bins["id"] = bins.index
        chrombins = bins.loc[bins.iloc[:, 0] == chrom, :]
        start = max(start, 1)
        start = max(chrombins.id[chrombins.iloc[:, 1] < start])
        end = max(chrombins.id[chrombins.iloc[:, 1] < end])
    else:
        chrom = ucsc_str
        bins = bins.copy()
        bins["id"] = bins.index
        chrombins = bins.loc[bins.iloc[:, 0] == chrom, :]
        try:
            start = min(chrombins["id"])
            end = max(chrombins["id"])
        except ValueError:
            logger.error("Invalid chromosome: %s", chrom)
            raise
    return (int(start), int(end))


def _check_output_path(path: str, force: bool = False) -> None:
    """Raise OSError if the output path already exists and force is False."""
    if not force and os.path.exists(path):
        raise OSError("Output file already exists. Use --force to overwrite.")
    if dirname(path):
        os.makedirs(dirname(path), exist_ok=True)


def _data_transform(dense_map, operation: str):
    """Apply a mathematical transformation to a dense Hi-C map.

    Supported operations: ``log2``, ``log10``, ``ln``, ``sqrt``,
    ``exp<float>`` (e.g. ``exp0.2``).
    """
    ops = {
        "log10": np.log10,
        "log2": np.log2,
        "ln": np.log,
        "sqrt": np.sqrt,
    }
    if operation in ops:
        return ops[operation](dense_map)
    if re.match(r"^exp", operation):
        exp_val = float(operation.split("exp")[1])
        return dense_map**exp_val
    if hasattr(np, operation) and callable(getattr(np, operation)):
        logger.warning("Using built-in numpy callable: %s", operation)
        return getattr(np, operation)(dense_map)
    raise TypeError(f"Supplied transform '{operation}' is not supported.")


# ---------------------------------------------------------------------------
# CLI group
# ---------------------------------------------------------------------------


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(version=__version__, prog_name="hicstuff")
def cli():
    """Simple Hi-C pipeline for generating and manipulating contact matrices."""


# ---------------------------------------------------------------------------
# iteralign
# ---------------------------------------------------------------------------


@cli.command("iteralign")
@click.argument("reads_fq")
@click.option(
    "-g", "--genome", required=True, metavar="FILE", help="Reference genome or aligner index path."
)
@click.option("-o", "--out-bam", required=True, metavar="FILE", help="Output BAM file path.")
@click.option(
    "-t",
    "--threads",
    default="1",
    show_default=True,
    metavar="INT",
    help="Number of parallel threads.",
)
@click.option(
    "-T",
    "--tempdir",
    default=None,
    metavar="DIR",
    help="Temporary directory (default: current directory).",
)
@click.option(
    "-a",
    "--aligner",
    default="bowtie2",
    show_default=True,
    metavar="STR",
    help="Aligner: bowtie2, minimap2 or bwa.",
)
@click.option(
    "-l",
    "--min-len",
    default=20,
    show_default=True,
    type=int,
    metavar="INT",
    help="Minimum truncated read length.",
)
@click.option(
    "-R",
    "--read-len",
    default=None,
    type=int,
    metavar="INT",
    help="Maximum read length (estimated from first read if omitted).",
)
def iteralign(reads_fq, genome, out_bam, threads, tempdir, aligner, min_len, read_len):
    """Iteratively align reads to a reference genome.

    Truncates reads to 20 bp then iteratively extends and re-aligns unmapped
    reads to maximise the proportion of uniquely aligned reads in a 3C library.
    """
    if tempdir is None:
        tempdir = "."
    temp_directory = hio.generate_temp_dir(tempdir)
    try:
        hci.iterative_align(
            reads_fq,
            temp_directory,
            genome,
            threads,
            out_bam,
            aligner=aligner,
            min_len=min_len,
            read_len=read_len,
        )
    finally:
        shutil.rmtree(temp_directory)


# ---------------------------------------------------------------------------
# digest
# ---------------------------------------------------------------------------


@cli.command("digest")
@click.argument("fasta")
@click.option(
    "-e",
    "--enzyme",
    required=True,
    metavar="ENZ[,ENZ2,...]",
    help="Restriction enzyme name(s) or fixed chunk size in bp.",
)
@click.option(
    "-o",
    "--outdir",
    default=None,
    metavar="DIR",
    help="Output directory (default: current directory).",
)
@click.option(
    "-s",
    "--size",
    default=0,
    show_default=True,
    type=int,
    metavar="INT",
    help="Minimum fragment size to keep.",
)
@click.option("-c", "--circular", is_flag=True, help="Genome is circular.")
@click.option("-p", "--plot", is_flag=True, help="Show fragment length distribution histogram.")
@click.option(
    "-f", "--figdir", default=None, metavar="DIR", help="Directory to save the distribution figure."
)
@click.option("-F", "--force", is_flag=True, help="Overwrite existing output directory.")
def digest(fasta, enzyme, outdir, size, circular, plot, figdir, force):
    """Digest a genome FASTA into restriction fragments.

    Writes ``fragments_list.txt`` and ``info_contigs.txt`` to the output directory.
    """
    if outdir is None:
        outdir = os.getcwd()
    if os.path.exists(outdir):
        if not force:
            raise OSError("Output directory already exists. Use --force to overwrite.")
    else:
        os.makedirs(outdir, exist_ok=True)
    figpath = join(figdir, "frags_hist.pdf") if figdir else None
    enzyme_parsed = enzyme.split(",") if "," in enzyme else enzyme
    hcd.write_frag_info(fasta, enzyme_parsed, size, output_dir=outdir, circular=circular)
    hcd.frag_len(output_dir=outdir, plot=plot, fig_path=figpath)


# ---------------------------------------------------------------------------
# cutsite
# ---------------------------------------------------------------------------


@cli.command("cutsite")
@click.option("-F", "--forward", required=True, metavar="FILE", help="Forward reads FASTQ file.")
@click.option("-R", "--reverse", required=True, metavar="FILE", help="Reverse reads FASTQ file.")
@click.option(
    "-p",
    "--prefix",
    required=True,
    metavar="STR",
    help="Output prefix (suffixed with _R1.fq.gz / _R2.fq.gz).",
)
@click.option(
    "-e", "--enzyme", required=True, metavar="STR", help="Comma-separated restriction enzyme(s)."
)
@click.option(
    "-m",
    "--mode",
    default="for_vs_rev",
    show_default=True,
    metavar="STR",
    type=click.Choice(["for_vs_rev", "all", "pile"]),
    help="Fragment pairing mode.",
)
@click.option(
    "-s",
    "--seed-size",
    default=20,
    show_default=True,
    type=int,
    metavar="INT",
    help="Minimum read size after cutting.",
)
@click.option(
    "-t",
    "--threads",
    default=1,
    show_default=True,
    type=int,
    metavar="INT",
    help="Number of parallel threads.",
)
def cutsite(forward, reverse, prefix, enzyme, mode, seed_size, threads):
    """Preprocess FASTQ files by cutting reads at religation sites.

    Generates gzipped FASTQ files with reads cut at ligation junctions,
    creating new fragment pairs for mapping.
    """
    if dirname(prefix):
        os.makedirs(dirname(prefix), exist_ok=True)
    logger.info("Digestion of reads with enzyme: %s, mode: %s", enzyme, mode)
    hcc.cut_ligation_sites(
        forward,
        reverse,
        prefix + "_R1.fq.gz",
        prefix + "_R2.fq.gz",
        enzyme=enzyme,
        mode=mode,
        seed_size=seed_size,
        n_cpu=threads,
    )


# ---------------------------------------------------------------------------
# filter
# ---------------------------------------------------------------------------


@cli.command("filter")
@click.argument("input_pairs")
@click.argument("output_pairs")
@click.option("-f", "--figdir", default=None, metavar="DIR", help="Directory for output figures.")
@click.option(
    "-i",
    "--interactive",
    is_flag=True,
    help="Ask for thresholds interactively after showing plots.",
)
@click.option(
    "-p", "--plot", is_flag=True, help="Show library composition and 3C event abundance plots."
)
@click.option(
    "-P", "--prefix", default=None, metavar="STR", help="Library name displayed on figures."
)
@click.option(
    "-t",
    "--thresholds",
    default=None,
    metavar="INT-INT",
    help="Manual thresholds as UNCUT-LOOP (e.g. 4-5).",
)
def filter(input_pairs, output_pairs, figdir, interactive, plot, prefix, thresholds):
    """Filter spurious Hi-C events (loops and uncuts) from a pairs file."""
    if thresholds:
        try:
            uncut_thr, loop_thr = (int(x) for x in thresholds.split("-"))
        except ValueError as err:
            raise click.BadParameter(
                "Thresholds must be two integers separated by '-', e.g. 4-5.",
                param_hint="--thresholds",
            ) from err
    else:
        figpath = None
        if figdir:
            os.makedirs(figdir, exist_ok=True)
            figpath = join(figdir, "event_distance.pdf")
        uncut_thr, loop_thr = hcf.get_thresholds(
            input_pairs,
            interactive=interactive,
            plot_events=plot,
            fig_path=figpath,
            prefix=prefix,
        )
    figpath = join(figdir, "event_distribution.pdf") if figdir else None
    hcf.filter_events(
        input_pairs,
        output_pairs,
        uncut_thr,
        loop_thr,
        plot_events=plot,
        fig_path=figpath,
        prefix=prefix,
    )


# ---------------------------------------------------------------------------
# view
# ---------------------------------------------------------------------------


@cli.command("view")
@click.argument("contact_map")
@click.argument("contact_map2", required=False, default=None)
@click.option(
    "-b",
    "--binning",
    default="1",
    show_default=True,
    metavar="INT[bp|kb|Mb]",
    help="Merge bins by factor or generate fixed-size bins.",
)
@click.option(
    "-c",
    "--cmap",
    default="Reds",
    show_default=True,
    metavar="STR",
    help="Matplotlib colormap name.",
)
@click.option("-C", "--circular", is_flag=True, help="Genome is circular.")
@click.option("-d", "--despeckle", is_flag=True, help="Remove speckle artefacts.")
@click.option(
    "-D", "--dpi", default=300, show_default=True, type=int, metavar="INT", help="Output image DPI."
)
@click.option(
    "-f",
    "--frags",
    default=None,
    metavar="FILE",
    help="fragments_list.txt (required for bp binning and --lines).",
)
@click.option(
    "-T",
    "--transform",
    default=None,
    metavar="STR",
    help="Pixel transform: log2, log10, ln, sqrt, exp<val>.",
)
@click.option(
    "-l", "--lines", is_flag=True, help="Draw chromosome separator lines (requires --frags)."
)
@click.option(
    "-M",
    "--max",
    "vmax",
    default="99%",
    show_default=True,
    metavar="NUM[%]",
    help="Colorscale maximum (percentile with %).",
)
@click.option(
    "-m",
    "--min",
    "vmin",
    default="0",
    show_default=True,
    metavar="NUM[%]",
    help="Colorscale minimum.",
)
@click.option(
    "-N",
    "--n-mad",
    default=3.0,
    show_default=True,
    type=float,
    metavar="FLOAT",
    help="MAD threshold for ICE normalization bin filtering.",
)
@click.option("-n", "--normalize", is_flag=True, help="Perform ICE normalization before rendering.")
@click.option(
    "-o",
    "--output",
    default=None,
    metavar="FILE",
    help="Output image path (display interactively if omitted).",
)
@click.option(
    "-r",
    "--region",
    default=None,
    metavar="STR[;STR]",
    help="UCSC region to zoom into (e.g. chr1:1000-12000).",
)
@click.option(
    "-t",
    "--trim",
    default=None,
    type=float,
    metavar="FLOAT",
    help="Trim bins deviating by more than this many MADs.",
)
def view(
    contact_map,
    contact_map2,
    binning,
    cmap,
    circular,
    despeckle,
    dpi,
    frags,
    transform,
    lines,
    vmax,
    vmin,
    n_mad,
    normalize,
    output,
    region,
    trim,
):
    """Visualize a Hi-C matrix as a heatmap."""
    hic_fmt = hio.get_hic_format(contact_map)

    # Switch to divergent colormap for ratio plots
    if contact_map2 is not None and cmap not in DIVERGENT_CMAPS:
        if cmap != "Reds":
            logger.warning(
                "Non-divergent colormap selected for ratio plot. Divergent options: %s",
                " ".join(DIVERGENT_CMAPS),
            )
        logger.info(
            "Defaulting to seismic colormap for ratio. Specify another divergent colormap if desired."
        )
        cmap = "seismic"

    bin_str = binning.upper()
    symmetric = True
    try:
        binning_val = int(bin_str)
        bp_unit = False
    except ValueError as err:
        if re.match(r"^[0-9]+[KMG]?B[P]?$", bin_str):
            if hic_fmt == "graal" and frags is None:
                raise click.UsageError(
                    "A fragment file (--frags) is required for basepair binning of graal matrices."
                ) from err
            binning_val = _parse_bin_str(bin_str)
            bp_unit = True
        else:
            raise click.BadParameter(
                f"Invalid binning value '{binning}'. Use an integer or a basepair string (e.g. 1kb).",
                param_hint="--binning",
            ) from None

    sparse_map, frags_df, _ = hio.flexible_hic_loader(contact_map, fragments_file=frags, quiet=True)

    def _process_matrix(sparse_map_in):
        nonlocal symmetric
        # Binning
        if binning_val > 1:
            if bp_unit:
                pos = frags_df.iloc[:, 2]
                binned_map, binned_pos = hcs.bin_bp_sparse(
                    M=sparse_map_in, positions=pos, bin_len=binning_val
                )
                binned_start = np.append(np.where(binned_pos == 0)[0], len(binned_pos))
                num_binned = binned_start[1:] - binned_start[:-1]
                chr_names_idx = np.unique(frags_df.iloc[:, 1], return_index=True)[1]
                chr_names = [frags_df.iloc[index, 1] for index in sorted(chr_names_idx)]
                binned_chrom = np.repeat(chr_names, num_binned)
                binned_frags = pd.DataFrame({"chrom": binned_chrom, "start_pos": binned_pos[:, 0]})
                binned_frags["end_pos"] = binned_frags.groupby("chrom")["start_pos"].shift(-1)
                chrom_ends = frags_df.groupby("chrom").end_pos.max()
                for cn in chrom_ends.index:
                    binned_frags.loc[
                        np.isnan(binned_frags.end_pos) & (binned_frags.chrom == cn), "end_pos"
                    ] = chrom_ends[cn]
            else:
                binned_map = hcs.bin_sparse(M=sparse_map_in, subsampling_factor=binning_val)
                if frags_df is not None:
                    binned_frags = frags_df.iloc[::binning_val, :].reset_index(drop=True)

                    def _shift_min(x):
                        try:
                            x[x == min(x)] = 0
                        except ValueError:
                            pass
                        return x

                    binned_frags.start_pos = binned_frags.groupby(
                        "chrom", sort=False
                    ).start_pos.apply(_shift_min)
                else:
                    binned_frags = frags_df
        else:
            binned_map = sparse_map_in
            binned_frags = frags_df

        # Chromosome separator lines
        chrom_starts = None
        if lines and binned_frags is not None:
            chrom_starts = np.where(np.diff(binned_frags.start_pos) < 0)[0] + 1

        # Trimming
        if trim is not None:
            binned_map, chrom_starts = hcs.trim_sparse(
                binned_map, n_mad=trim, chrom_start=chrom_starts
            )

        # Normalization
        if normalize:
            binned_map = hcs.normalize_sparse(binned_map, norm="ICE", n_mad=n_mad)

        # Region zoom
        if region:
            if lines:
                raise NotImplementedError("Chromosome lines are incompatible with a region zoom.")
            if frags_df is None:
                raise click.UsageError(
                    "A fragment file (--frags) is required to zoom into a genomic region."
                )
            reg_pos = binned_frags[["chrom", "start_pos"]]
            if ";" in region:
                symmetric = False
                reg1_str, reg2_str = region.split(";")
                reg1 = _parse_ucsc(reg1_str, reg_pos)
                reg2 = _parse_ucsc(reg2_str, reg_pos)
            else:
                reg = _parse_ucsc(region, reg_pos)
                reg1 = reg2 = reg
            binned_map = binned_map.tocsr()[reg1[0] : reg1[1], reg2[0] : reg2[1]].tocoo()

        return binned_map, chrom_starts

    processed_map, chrom_starts = _process_matrix(sparse_map)

    # If a second matrix was provided, compute the log2 ratio
    if contact_map2 is not None:
        sparse_map2, _, _ = hio.flexible_hic_loader(contact_map2, fragments_file=frags, quiet=True)
        processed_map2, _ = _process_matrix(sparse_map2)
        if sparse_map2.shape != sparse_map.shape:
            raise click.UsageError("Cannot compute ratio of matrices with different dimensions.")
        processed_map.data = np.log2(processed_map.data)
        processed_map2.data = np.log2(processed_map2.data)
        processed_map = (processed_map.tocsr() - processed_map2.tocsr()).tocoo()
        processed_map.data[np.isnan(processed_map.data)] = 0.0
        transform = None

    if despeckle:
        processed_map = hcs.despeckle_simple(processed_map)

    try:
        if symmetric:
            dense_map = hcv.sparse_to_dense(processed_map, remove_diag=False)
        else:
            dense_map = processed_map.toarray()

        def _set_v(v, mat):
            if "%" in str(v):
                try:
                    valid_pixels = (mat > 0) & (np.isfinite(mat))
                    return np.percentile(mat[valid_pixels], float(str(v).strip("%")))
                except IndexError:
                    return 0
            return float(v)

        dense_map = dense_map.astype(float)
        vmax_val = _set_v(vmax, dense_map)
        vmin_val = _set_v(vmin, dense_map)
        if contact_map2 is not None:
            vmin_val, vmax_val = -2, 2
        if transform:
            dense_map = _data_transform(dense_map, transform)
            vmax_val = _set_v(vmax, dense_map)
            vmin_val = _set_v(vmin, dense_map)
        else:
            dense_map[dense_map == 0] = np.inf
        current_cmap = cm.get_cmap().copy()
        current_cmap.set_bad(color=current_cmap(0))

        hcv.plot_matrix(
            dense_map,
            filename=output,
            vmin=vmin_val,
            vmax=vmax_val,
            dpi=dpi,
            cmap=cmap,
            chrom_starts=chrom_starts,
        )
    except MemoryError:
        logger.error("Contact map is too large to load. Try binning more.")


# ---------------------------------------------------------------------------
# pipeline
# ---------------------------------------------------------------------------


@cli.command("pipeline")
@click.argument("input1")
@click.argument("input2", required=False, default=None)
@click.option(
    "-g", "--genome", required=True, metavar="FILE", help="Reference genome or aligner index."
)
@click.option(
    "-a",
    "--aligner",
    default="bowtie2",
    show_default=True,
    metavar="STR",
    help="Aligner: bowtie2, minimap2 or bwa.",
)
@click.option(
    "-B",
    "--balancing-args",
    default=None,
    metavar="STR",
    help="Extra arguments passed to `cooler balance`.",
)
@click.option(
    "-b",
    "--binning",
    default=0,
    show_default=True,
    type=int,
    metavar="INT",
    help="Bin the cool matrix to this resolution (bp). 0 means no binning.",
)
@click.option(
    "-c", "--centromeres", default=None, metavar="FILE", help="Centromere positions file."
)
@click.option("-C", "--circular", is_flag=True, help="Genome is circular.")
@click.option("-d", "--distance-law", is_flag=True, help="Generate a distance law output file.")
@click.option("-D", "--duplicates", is_flag=True, help="Filter PCR duplicates.")
@click.option(
    "-e",
    "--enzyme",
    default="5000",
    show_default=True,
    metavar="{STR|INT}",
    help="Restriction enzyme name, 'mnase'/'dnase', or chunk size in bp.",
)
@click.option(
    "-E",
    "--exclude",
    default=None,
    metavar="STR",
    help="Comma-separated chromosomes to exclude (e.g. chrM,2u).",
)
@click.option(
    "-f",
    "--filter",
    "filter_events",
    is_flag=True,
    help="Filter spurious 3C events (loops and uncuts).",
)
@click.option("-F", "--force", is_flag=True, help="Overwrite existing output files.")
@click.option(
    "-m",
    "--mapping",
    default="normal",
    show_default=True,
    metavar="STR",
    type=click.Choice(["normal", "iterative", "cutsite"]),
    help="Mapping mode.",
)
@click.option(
    "-M",
    "--matfmt",
    default="cool",
    show_default=True,
    metavar="STR",
    type=click.Choice(["graal", "bg2", "cool"]),
    help="Output matrix format.",
)
@click.option("-n", "--no-cleanup", is_flag=True, help="Keep intermediary files.")
@click.option(
    "-o",
    "--outdir",
    default=None,
    metavar="DIR",
    help="Output directory (default: current directory).",
)
@click.option("-p", "--plot", is_flag=True, help="Generate plots at pipeline steps.")
@click.option("-P", "--prefix", default=None, metavar="STR", help="Prefix for all output files.")
@click.option(
    "-q",
    "--quality-min",
    default=30,
    show_default=True,
    type=int,
    metavar="INT",
    help="Minimum mapping quality.",
)
@click.option(
    "-r",
    "--remove-centromeres",
    default=0,
    show_default=True,
    type=int,
    metavar="INT",
    help="kb to remove around centromere positions.",
)
@click.option(
    "-R",
    "--read-len",
    default=None,
    type=int,
    metavar="INT",
    help="Maximum read length (estimated from first read if omitted).",
)
@click.option(
    "-s",
    "--size",
    default=0,
    show_default=True,
    type=int,
    metavar="INT",
    help="Minimum contig size threshold.",
)
@click.option(
    "-S",
    "--start-stage",
    default="fastq",
    show_default=True,
    metavar="STR",
    type=click.Choice(["fastq", "bam", "pairs", "pairs_idx"]),
    help="Pipeline start stage.",
)
@click.option(
    "-t",
    "--threads",
    default=1,
    show_default=True,
    type=int,
    metavar="INT",
    help="Number of parallel threads.",
)
@click.option("-T", "--tmpdir", default=None, metavar="DIR", help="Temporary directory.")
@click.option(
    "-z",
    "--zoomify",
    default=True,
    show_default=True,
    type=bool,
    metavar="BOOL",
    help="Generate multi-resolution .mcool from the binned cool matrix.",
)
@click.option(
    "--skip-count",
    is_flag=True,
    default=False,
    help="Skip the read-count check on input FASTQ files.",
)
def pipeline(
    input1,
    input2,
    genome,
    aligner,
    balancing_args,
    binning,
    centromeres,
    circular,
    distance_law,
    duplicates,
    enzyme,
    exclude,
    filter_events,
    force,
    mapping,
    matfmt,
    no_cleanup,
    outdir,
    plot,
    prefix,
    quality_min,
    remove_centromeres,
    read_len,
    size,
    start_stage,
    threads,
    tmpdir,
    zoomify,
    skip_count,
):
    """Run the full Hi-C pipeline from FASTQ to contact matrix.

    \b
    Example — generate a multi-resolution .mcool for Arima Hi-C:

        hicstuff pipeline --enzyme "DpnII,HinfI" --binning 1000 --threads 8 \\
            --genome ref.fa R1.fq.gz R2.fq.gz
    """
    if filter_events and enzyme.isdigit():
        raise click.UsageError("Cannot filter events without specifying a restriction enzyme.")
    if enzyme in ("mnase", "dnase"):
        logger.info("Enzyme is '%s', setting chunk size to 100 bp.", enzyme)
        enzyme = "100"
    if outdir is None:
        outdir = os.getcwd()
    hpi.full_pipeline(
        genome=genome,
        input1=input1,
        input2=input2,
        aligner=aligner,
        centromeres=centromeres,
        circular=circular,
        exclude=exclude,
        distance_law=distance_law,
        enzyme=enzyme,
        filter_events=filter_events,
        force=force,
        mapping=mapping,
        mat_fmt=matfmt,
        binning=binning,
        zoomify=zoomify,
        balancing_args=balancing_args,
        min_qual=quality_min,
        min_size=size,
        no_cleanup=no_cleanup,
        out_dir=outdir,
        pcr_duplicates=duplicates,
        plot=plot,
        prefix=prefix,
        read_len=read_len,
        remove_centros=remove_centromeres,
        start_stage=start_stage,
        threads=threads,
        tmp_dir=tmpdir,
        skip_count=skip_count,
    )


# ---------------------------------------------------------------------------
# scalogram
# ---------------------------------------------------------------------------


@cli.command("scalogram")
@click.argument("contact_map")
@click.option(
    "-C", "--cmap", default="viridis", show_default=True, metavar="STR", help="Matplotlib colormap."
)
@click.option("-d", "--despeckle", is_flag=True, help="Remove speckle artefacts before plotting.")
@click.option(
    "-f",
    "--frags",
    default=None,
    metavar="FILE",
    help="fragments_list.txt for coordinate conversion.",
)
@click.option(
    "-i",
    "--indices",
    default=None,
    metavar="INT-INT",
    help="Bin range or UCSC coordinates to display.",
)
@click.option(
    "-o",
    "--output",
    default=None,
    metavar="FILE",
    help="Output image path (display interactively if omitted).",
)
@click.option("-n", "--normalize", is_flag=True, help="ICE-normalize the matrix before plotting.")
@click.option(
    "-r",
    "--range",
    "range_str",
    default=None,
    metavar="INT-INT",
    help="Contact distance range to display.",
)
@click.option(
    "-t",
    "--threads",
    default=1,
    show_default=True,
    type=int,
    metavar="INT",
    help="Parallel threads for despeckling.",
)
def scalogram(contact_map, cmap, despeckle, frags, indices, output, normalize, range_str, threads):
    """Generate a scalogram from a Hi-C contact matrix."""
    mat, frags_df, _ = hio.flexible_hic_loader(contact_map, fragments_file=frags)
    if frags_df is not None and frags is not None:
        frags_df = pd.read_csv(frags, delimiter="\t", usecols=(1, 2, 3))

    shortest, longest = 0, None
    if range_str:
        shortest_str, longest_str = range_str.split("-")
        try:
            shortest, longest = int(shortest_str), int(longest_str)
        except ValueError:
            shortest = _parse_bin_str(shortest_str)
            longest = _parse_bin_str(longest_str)
            avg_res = (frags_df.end_pos - frags_df.start_pos).mean()
            shortest = int(shortest // avg_res)
            longest = int(longest // avg_res)

    start, end = None, None
    if indices:
        try:
            start, end = (int(x) for x in indices.split("-"))
        except ValueError:
            start, end = _parse_ucsc(indices, frags_df.loc[:, ["chrom", "start_pos"]])

    S = mat.tocsr()
    if longest is None:
        longest = S.shape[0]
    if normalize:
        S = hcs.normalize_sparse(S, norm="ICE").tocsr()
    if despeckle:
        S = hcs.despeckle_simple(S, threads=threads)

    if start is not None and end is not None:
        crop_inf = max(0, start - longest)
        crop_sup = min(S.shape[0], end + longest)
        crop_later = longest
        S = S[crop_inf:crop_sup, crop_inf:crop_sup]
    else:
        crop_later = 0

    D = hcv.sparse_to_dense(S)
    D = np.fliplr(np.rot90(hcs.scalogram(D), k=-1))
    plt.contourf(D[crop_later : D.shape[1] - crop_later, shortest:longest], cmap=cmap)
    if output:
        plt.savefig(output)
    else:
        plt.show()


# ---------------------------------------------------------------------------
# rebin
# ---------------------------------------------------------------------------


@cli.command("rebin")
@click.argument("contact_map")
@click.argument("out_prefix")
@click.option(
    "-b",
    "--binning",
    default="1",
    show_default=True,
    metavar="INT[bp|kb|Mb]",
    help="Binning factor or basepair bin size.",
)
@click.option("-f", "--frags", default=None, metavar="FILE", help="fragments_list.txt file.")
@click.option("-c", "--chroms", default=None, metavar="FILE", help="info_contigs.txt file.")
@click.option("-F", "--force", is_flag=True, help="Overwrite existing output.")
def rebin(contact_map, out_prefix, binning, frags, chroms, force):
    """Rebin a Hi-C matrix to a coarser resolution.

    Output files use the same format as the input (cool, bg2 or graal).
    """
    hic_fmt = hio.get_hic_format(contact_map)
    out_name = out_prefix + _FMT2EXT[hic_fmt]
    _check_output_path(out_name, force=force)

    hic_map, frags_df, chromlist = hio.flexible_hic_loader(
        contact_map,
        fragments_file=frags,
        chroms_file=chroms,
    )
    if hic_fmt == "graal" and (frags_df is None or chromlist is None):
        raise click.UsageError(
            "Graal format requires --frags (fragments_list.txt) and --chroms (info_contigs.txt)."
        )
    if dirname(out_prefix):
        os.makedirs(dirname(out_prefix), exist_ok=True)

    bin_str = binning.upper()
    try:
        binning_val = int(bin_str)
        bp_unit = False
    except ValueError:
        if re.match(r"^[0-9]+[KMG]?B[P]?$", bin_str):
            binning_val = _parse_bin_str(bin_str)
            bp_unit = True
        else:
            raise click.BadParameter(
                f"Invalid binning '{binning}'.", param_hint="--binning"
            ) from None

    chromnames = np.unique(frags_df.chrom)

    if bp_unit:
        hic_map, _ = hcs.bin_bp_sparse(hic_map, frags_df.start_pos, binning_val)
        for chrom in chromnames:
            chrom_mask = frags_df.chrom == chrom
            bin_id = frags_df.loc[chrom_mask, "start_pos"] // binning_val
            frags_df.loc[chrom_mask, "id"] = bin_id + 1
            frags_df.loc[chrom_mask, "start_pos"] = binning_val * bin_id
            bin_ends = binning_val * bin_id + binning_val
            try:
                chromsize = chromlist.length[chromlist.contig == chrom].values[0]
            except AttributeError:
                chromsize = chromlist["length_kb"][chromlist.contig == chrom].values[0]
            bin_ends[bin_ends > chromsize] = chromsize
            frags_df.loc[frags_df.chrom == chrom, "end_pos"] = bin_ends

        id_diff = np.array(frags_df.loc[:, "id"])[1:] - np.array(frags_df.loc[:, "id"])[:-1]
        jump_frag_idx = np.where(id_diff > 1)[0]
        add_bins = id_diff - 1
        miss_bins = [None] * np.sum(add_bins[jump_frag_idx])
        miss_bin_id = 0
        for idx in jump_frag_idx:
            jump_size = add_bins[idx]
            for j in range(1, jump_size + 1):
                miss_bins[miss_bin_id] = {
                    "id": frags_df.loc[idx, "id"] + j,
                    "chrom": frags_df.loc[idx, "chrom"],
                    "start_pos": frags_df.loc[idx, "start_pos"] + binning_val * j,
                    "end_pos": frags_df.loc[idx, "end_pos"] + binning_val * j,
                    "size": binning_val,
                    "gc_content": np.nan,
                }
                miss_bin_id += 1

        idx_shift = copy.copy(id_diff)
        idx_shift[idx_shift < 1] = 1
        existing_bins_idx = np.cumsum(idx_shift)
        existing_bins_idx = np.insert(existing_bins_idx, 0, 0)
        missing_bins_idx = sorted(
            set(range(existing_bins_idx[0], existing_bins_idx[-1])) - set(existing_bins_idx)
        )
        miss_bins_df = pd.DataFrame(miss_bins, columns=frags_df.columns, index=missing_bins_idx)
        for col in frags_df.columns:
            if frags_df[col].dtype.name != "category":
                try:
                    miss_bins_df[col] = pd.to_numeric(miss_bins_df[col])
                except (ValueError, TypeError):
                    pass
        frags_df["tmp_idx"] = existing_bins_idx
        miss_bins_df["tmp_idx"] = missing_bins_idx
        frags_df = pd.concat([frags_df, miss_bins_df], axis=0, sort=False)
        frags_df.sort_values("tmp_idx", axis=0, inplace=True)
        frags_df.drop("tmp_idx", axis=1, inplace=True)
    else:
        hic_map = hcs.bin_sparse(hic_map, binning_val)
        shift_id = 0 if binning_val == 1 else 1
        frags_df.id = (frags_df.id // binning_val) + shift_id

    col_ordered = list(frags_df.columns)
    frags_df = frags_df.groupby(["chrom", "id"], sort=False, observed=True)
    positions = frags_df.agg({"start_pos": "min", "end_pos": "max"})
    positions.reset_index(inplace=True)
    try:
        features = frags_df.agg("mean")
        features.reset_index(inplace=True)
        frags_df = features
        frags_df["start_pos"] = 0
        frags_df["end_pos"] = 0
        frags_df.loc[:, positions.columns] = positions
    except pd.errors.DataError:
        frags_df = positions
    frags_df["size"] = frags_df.end_pos - frags_df.start_pos
    cumul_bins = 0
    for chrom in chromnames:
        chrom_frags = frags_df.chrom == chrom
        n_bins = frags_df.start_pos[chrom_frags].shape[0]
        chromlist.loc[chromlist.contig == chrom, "n_frags"] = n_bins
        chromlist.loc[chromlist.contig == chrom, "cumul_length"] = cumul_bins
        cumul_bins += n_bins
        last_frag_end = frags_df.loc[chrom_frags, "end_pos"].max()
        chromlen = chromlist.loc[chromlist.contig == chrom, "length"].values[0]
        frags_df.loc[chrom_frags & (frags_df.end_pos == last_frag_end), "end_pos"] = chromlen

    frags_df = frags_df.reindex(columns=col_ordered)
    hio.flexible_hic_saver(hic_map, out_prefix, frags=frags_df, chroms=chromlist, hic_fmt=hic_fmt)


# ---------------------------------------------------------------------------
# subsample
# ---------------------------------------------------------------------------


@cli.command("subsample")
@click.argument("contact_map")
@click.argument("subsampled_prefix")
@click.option(
    "-p",
    "--prop",
    default=0.1,
    show_default=True,
    type=float,
    metavar="FLOAT",
    help="Proportion (0–1) or raw contact count (>1) to keep.",
)
@click.option("-F", "--force", is_flag=True, help="Overwrite existing output.")
def subsample(contact_map, subsampled_prefix, prop, force):
    """Subsample contacts from a Hi-C matrix.

    Sampling probability is proportional to the intensity of each bin.
    """
    hic_fmt = hio.get_hic_format(contact_map)
    out_name = subsampled_prefix + _FMT2EXT[hic_fmt]
    _check_output_path(out_name, force=force)
    mat, frags_df, _ = hio.flexible_hic_loader(contact_map, quiet=True)
    subsampled = hcs.subsample_contacts(mat, prop).tocoo()
    hio.flexible_hic_saver(
        subsampled, subsampled_prefix, frags=frags_df, hic_fmt=hic_fmt, quiet=True
    )


# ---------------------------------------------------------------------------
# convert
# ---------------------------------------------------------------------------


@cli.command("convert")
@click.argument("contact_map")
@click.argument("prefix")
@click.option("-f", "--frags", default=None, metavar="FILE", help="fragments_list.txt file.")
@click.option("-c", "--chroms", default=None, metavar="FILE", help="info_contigs.txt file.")
@click.option("-F", "--force", is_flag=True, help="Overwrite existing output.")
@click.option(
    "-g",
    "--genome",
    default=None,
    metavar="FILE",
    help="Genome FASTA to compute GC content column.",
)
@click.option(
    "-T",
    "--to",
    "out_fmt",
    default="cool",
    show_default=True,
    metavar="STR",
    type=click.Choice(["graal", "bg2", "cool"]),
    help="Output format.",
)
def convert(contact_map, prefix, frags, chroms, force, genome, out_fmt):
    """Convert between Hi-C matrix formats (graal, bg2, cool).

    Input format is automatically inferred.
    """
    out_name = prefix + _FMT2EXT[out_fmt]
    _check_output_path(out_name, force=force)
    mat, frags_df, chroms_df = hio.flexible_hic_loader(
        contact_map,
        fragments_file=frags,
        chroms_file=chroms,
        quiet=True,
    )
    chrom_col, start_col, end_col = hio.get_pos_cols(frags_df)
    size = frags_df[end_col] - frags_df[start_col]
    if "size" not in frags_df.columns:
        frags_df = frags_df.join(pd.DataFrame({"size": size}))
    if genome:
        gc = hio.gc_bins(genome, frags_df)
        frags_df = frags_df.join(pd.DataFrame({"gc_content": gc}))
    hio.flexible_hic_saver(
        mat=mat.astype(int),
        out_prefix=prefix,
        frags=frags_df,
        chroms=chroms_df,
        hic_fmt=out_fmt,
    )


# ---------------------------------------------------------------------------
# distancelaw
# ---------------------------------------------------------------------------


@cli.command("distancelaw")
@click.option(
    "-p", "--pairs", "pairs_file", default=None, metavar="FILE", help="Input indexed pairs file."
)
@click.option(
    "-d",
    "--dist-tbl",
    default=None,
    metavar="FILE1[,FILE2,...]",
    help="Pre-computed distance law file(s), comma-separated.",
)
@click.option("-f", "--frags", default=None, metavar="FILE", help="fragments_list.txt file.")
@click.option("-a", "--average", is_flag=True, help="Average distance law across chromosomes/arms.")
@click.option(
    "-b",
    "--big-arm-only",
    default=None,
    type=int,
    metavar="INT",
    help="Only use arms larger than this value (bp).",
)
@click.option(
    "-B",
    "--base",
    default=1.1,
    show_default=True,
    type=float,
    metavar="FLOAT",
    help="Log base for genomic bin spacing.",
)
@click.option(
    "-c", "--centromeres", default=None, metavar="FILE", help="Centromere positions file."
)
@click.option("-C", "--circular", is_flag=True, help="Genome is circular.")
@click.option(
    "-i",
    "--inf",
    default=3000,
    show_default=True,
    type=int,
    metavar="INT",
    help="Minimum distance to plot (bp).",
)
@click.option(
    "-l",
    "--labels",
    default=None,
    metavar="STR1,STR2,...",
    help="Comma-separated sample labels for the plot.",
)
@click.option("-o", "--outputfile-img", default=None, metavar="FILE", help="Output image path.")
@click.option("-O", "--outputfile-tabl", default=None, metavar="FILE", help="Output table path.")
@click.option(
    "-r",
    "--remove-centromeres",
    default=0,
    show_default=True,
    type=int,
    metavar="INT",
    help="kb to remove around centromere positions.",
)
@click.option(
    "-s", "--sup", default=None, type=int, metavar="INT", help="Maximum distance to plot (bp)."
)
def distancelaw(
    pairs_file,
    dist_tbl,
    frags,
    average,
    big_arm_only,
    base,
    centromeres,
    circular,
    inf,
    labels,
    outputfile_img,
    outputfile_tabl,
    remove_centromeres,
    sup,
):
    """Analyse and plot the Hi-C distance law (P(s) curve)."""
    if pairs_file and dist_tbl:
        raise click.UsageError("Cannot use --pairs and --dist-tbl simultaneously.")
    if not pairs_file and not dist_tbl:
        raise click.UsageError("Provide either --pairs or --dist-tbl.")

    if pairs_file:
        xs = [None]
        ps_list = [None]
        names = [None]
        xs[0], ps_list[0], names[0] = hcdl.get_distance_law(
            pairs_reads_file=pairs_file,
            fragments_file=frags if frags else None,
            centro_file=centromeres,
            base=base,
            out_file=outputfile_tabl,
            circular=circular if circular else None,
            rm_centro=remove_centromeres,
        )
        length_files = 1
    else:
        distance_law_files = dist_tbl.split(",")
        length_files = len(distance_law_files)
        xs = [None] * length_files
        ps_list = [None] * length_files
        names = [None] * length_files
        for i in range(length_files):
            xs[i], ps_list[i], names[i] = hcdl.import_distance_law(distance_law_files[i])
        names = [name[0] for name in names]

    if sup is None:
        sup = max(max(xs[0], key=len))
    arm_sup = big_arm_only if big_arm_only is not None else sup

    if not average and length_files > 1:
        raise click.UsageError("--average is required when providing more than one file.")

    for i in range(length_files):
        if average:
            xs[i], ps_list[i] = hcdl.average_distance_law(
                xs[i], ps_list[i], arm_sup, big_arm_only is not None
            )

    if not average:
        names = names[0]
        xs = xs[0]
        ps_list = ps_list[0]

    ps_list = hcdl.normalize_distance_law(xs, ps_list, inf, arm_sup)

    if labels:
        plot_labels = labels.split(",")
    elif length_files == 1 and not average:
        plot_labels = list(names) if len(names) > 1 else names
    else:
        plot_labels = [f"Sample {i}" for i in range(length_files)]

    hcdl.plot_ps_slope(xs, ps_list, plot_labels, outputfile_img, inf, sup)
    if outputfile_tabl:
        hcdl.export_distance_law(xs, ps_list, plot_labels, outputfile_tabl)


# ---------------------------------------------------------------------------
# missview
# ---------------------------------------------------------------------------


@cli.command("missview")
@click.argument("genome")
@click.argument("output")
@click.option(
    "-R", "--read-len", required=True, type=int, metavar="INT", help="Simulated read length (bp)."
)
@click.option(
    "-a",
    "--aligner",
    default="bowtie2",
    show_default=True,
    metavar="STR",
    help="Aligner: bowtie2, minimap2 or bwa.",
)
@click.option(
    "-b",
    "--binning",
    default="5000",
    show_default=True,
    metavar="INT",
    help="Resolution for the preview map.",
)
@click.option("-F", "--force", is_flag=True, help="Overwrite existing output.")
@click.option(
    "-t",
    "--threads",
    default=1,
    show_default=True,
    type=int,
    metavar="INT",
    help="Number of parallel threads.",
)
@click.option("-T", "--tmpdir", default=None, metavar="DIR", help="Temporary directory.")
def missview(genome, output, read_len, aligner, binning, force, threads, tmpdir):
    """Preview unmappable Hi-C bins for a given read length.

    Simulates short reads across the genome and identifies genomic bins that
    will be systematically missing due to repetitive sequences.
    """
    resolution = _parse_bin_str(str(binning))
    if tmpdir is None:
        tmpdir = tempfile.TemporaryDirectory().name
    phred = "F" * read_len
    tmp_fq = join(tmpdir, "simulated_reads.fq")
    tmp_bam = join(tmpdir, "simulated_reads.bam")
    _check_output_path(tmp_fq, force=force)

    logger.info("Simulating reads by splitting the genome into %i bp chunks", read_len)
    with open(tmp_fq, "w") as fq_handle:
        for rec in SeqIO.parse(genome, "fasta"):
            for i in range(len(rec.seq) - read_len):
                fq_handle.write(
                    f"@NS_SIM_{rec.id}_{i}\n{str(rec.seq[i : i + read_len])}\n+\n{phred}\n"
                )

    hpi.align_reads(tmp_fq, genome, tmp_bam, tmp_dir=tmpdir, threads=threads, aligner=aligner)
    ps.sort("-@", str(threads), "-n", "-O", "BAM", "-o", tmp_bam + ".sorted", tmp_bam)
    shutil.move(tmp_bam + ".sorted", tmp_bam)

    hpi.full_pipeline(
        genome,
        tmp_bam,
        tmp_bam,
        start_stage="bam",
        aligner=aligner,
        enzyme=resolution,
        force=force,
        threads=threads,
        out_dir=tmpdir,
        tmp_dir=tmpdir,
    )

    mat_path = join(tmpdir, "abs_fragments_contacts_weighted.txt")
    mat = hio.load_sparse_matrix(mat_path)
    log_files = glob.glob(join(tmpdir, "*.log"))
    log_content = open(log_files[0]).read()
    prop_mapped = int(re.search(r".*INFO :: ([0-9]*)% reads.*", log_content)[1]) / 100
    logger.info(
        "Bins with less than %s%% mapped reads will be considered undetectable.",
        100 * prop_mapped,
    )
    unmappable = mat.diagonal(0) < prop_mapped * resolution
    mappable_mat = np.ones(mat.shape)
    mappable_mat[unmappable, :] = 0
    mappable_mat[:, unmappable] = 0
    hcv.plot_matrix(
        mappable_mat,
        filename=output,
        title=(
            f"{100 * sum(unmappable) / len(unmappable):.3f}% missing bins"
            f" for {os.path.basename(genome)}"
            f" with {read_len} bp reads at resolution {resolution}."
        ),
        dpi=600,
        vmax=2,
        cmap="Greys",
    )
    logger.info("Output image saved at %s.", output)


# ---------------------------------------------------------------------------
# stats
# ---------------------------------------------------------------------------


@cli.command("stats")
@click.argument("log_file")
def stats(log_file):
    """Extract mapping statistics from a hicstuff pipeline log file."""
    pipeline_stats = hcstats.get_pipeline_stats(log_file)
    hcstats.print_pipeline_stats(pipeline_stats)
