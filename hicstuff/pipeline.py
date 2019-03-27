"""
Handle generation of GRAAL-compatible contact maps from fastq files.
cmdoret, 20190322
"""
import os, sys, time, csv
from datetime import datetime
from dateutil.relativedelta import relativedelta
import shutil as st
import itertools
from os.path import join
import subprocess as sp
from Bio import SeqIO
import pandas as pd
import hicstuff.digest as hcd
import hicstuff.iteralign as hci
import hicstuff.filter as hcf
import hicstuff.io as hio
import pysam as ps
from hicstuff.log import logger


def align_reads(
    reads,
    genome,
    out_sam,
    tmp_dir=None,
    threads=1,
    minimap2=False,
    iterative=False,
    min_qual=30,
):
    """
    Select and call correct alignment method and generate logs accordingly.
    Alignments are filtered so that there at most one alignment per read.

    Parameters
    ----------
    reads : str
        Path to the fastq file with Hi-C reads.
    genome : str
        Path to the genome in fasta format
    out_sam : str
        Path to the output SAM file containing mapped Hi-C reads.
    tmp_dir : str
        Path where temporary files are stored.
    threads : int
        Number of threads to run alignments in parallel.
    aligner : bool
        Use minimap2 instead of bowtie2.
    iterative : bool
        Wether to use the iterative mapping procedure (truncating reads and
        extending iteratively)
    min_qual : int
        Minimum mapping quality required to keep an alignment during iterative
        mapping.
    """
    if tmp_dir is None:
        tmp_dir = os.getcwd()
    index = None
    tmp_sam = out_sam + ".tmp"

    if iterative:
        hci.temp_directory = hci.generate_temp_dir(tmp_dir)
        hci.iterative_align(
            reads,
            tmp_dir=tmp_dir,
            ref=genome,
            n_cpu=threads,
            sam_out=tmp_sam,
            min_qual=min_qual,
        )
    else:
        if minimap2:
            map_cmd = "minimap2 -2 -t {threads} -ax sr {fasta} {fastq} > {sam}"
        else:
            index = hci.check_bt2_index(genome)
            map_cmd = "bowtie2 --very-sensitive-local -p {threads} -x {index} -U {fastq} > {sam}"
        map_args = {
            "threads": threads,
            "sam": tmp_sam,
            "fastq": reads,
            "fasta": genome,
            "index": index,
        }
        sp.call(map_cmd.format(**map_args), shell=True)

    # Remove supplementary alignments
    # TODO: replace sp.call with ps.view command. It currently has a bug
    # preventing output redirection.
    sp.call(
        "samtools view -F 2048 -h -@ {threads} -o {out} {tmp}".format(
            tmp=tmp_sam, threads=threads, out=out_sam
        ),
        shell=True,
    )
    os.remove(tmp_sam)


def sam2pairs(sam1, sam2, out_pairs, info_contigs, min_qual=30):
    """
    Make a .pairs file from two Hi-C sam files sorted by read names.
    The Hi-C mates are matched by read identifier. Pairs where at least one
    reads maps with MAPQ below  min_qual threshold are discarded. Pairs are
    sorted by readID and stored in upper triangle (first pair higher).

    Parameters
    ----------
    sam1 : str
        Path to the name-sorted SAM file with aligned Hi-C forward reads.
    sam2 : str
        Path to the name-sorted SAM file with aligned Hi-C reverse reads.
    out_pairs : str
        Path to the output space-separated .pairs file with columns 
        readID, chr1 pos1 chr2 pos2 strand1 strand2
    info_contigs : str
        Path to the info contigs file, to get info on chromosome sizes and order.
    min_qual : int
        Minimum mapping quality required to keep a Hi=C pair.
    """
    forward = ps.AlignmentFile(sam1)
    reverse = ps.AlignmentFile(sam2)

    # Generate header lines
    format_version = "## pairs format v1.0\n"
    sorting = "#sorted: readID\n"
    cols = "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n"
    # Chromosome order will be identical in info_contigs and pair files
    chroms = pd.read_csv(info_contigs, sep="\t").apply(
        lambda x: "#chromsize: %s %d\n" % (x.contig, x.length), axis=1
    )

    with open(out_pairs, "w") as pairs:
        pairs.writelines([format_version, sorting, cols, *chroms])
        pairs_writer = csv.writer(pairs, delimiter=" ")
        # Iterate on both SAM simultaneously
        for end1, end2 in itertools.zip_longest(forward, reverse):
            # Keep only pairs where both reads have good quality
            if (
                end1.mapping_quality >= min_qual
                and end2.mapping_quality >= min_qual
            ):
                if end1.query_name == end2.query_name:
                    if (
                        end1.reference_start > end2.reference_start
                        or end1.reference_id > end2.reference_id
                    ):
                        end1, end2 = end2, end1
                    pairs_writer.writerow(
                        [
                            end1.query_name,
                            end1.reference_name,
                            end1.reference_start + 1,
                            end2.reference_name,
                            end2.reference_start + 1,
                            "+" if end1.is_reverse else "-",
                            "+" if end2.is_reverse else "-",
                        ]
                    )
                else:
                    print(
                        "Error: Reads do not match between SAM files. "
                        "Verify both files are name-sorted and do not have "
                        "supplementary alignments."
                    )
                    sys.exit(1)
    pairs.close()


def pairs2matrix(pairs_file, mat_file, n_frags, mat_format="GRAAL", threads=1):
    """Generate the matrix by counting the number of occurences of each
    combination of restriction fragments in a 2D BED file.

    Parameters
    ----------
    pairs_file : str
        Path to a Hi-C pairs file, with frag1 and frag2 columns added.
    mat_file : str
        Path where the matrix will be written.
    n_frags : int
        Total number of restriction fragments in genome. Used to know total
        matrix size in case some observations are not observed at the end.
    mat_format : str
        The format to use when writing the matrix. Can be GRAAL or cooler format.
    threads : int
        Number of threads to use in parallel.
    """
    pre_mat_file = mat_file + ".pre.pairs"
    hio.sort_pairs(
        pairs_file, pre_mat_file, keys=["frag1", "frag2"], threads=threads
    )
    header_size = len(hio.get_pairs_header(pre_mat_file))
    time.sleep(
        1
    )  # Crashes if no sleep. File pointer must not be closed. Why ??
    with open(pre_mat_file, "r") as pairs, open(mat_file, "w") as mat:
        # Skip header lines
        for _ in range(header_size):
            next(pairs)
        prev_pair = ["0", "0"]  # Pairs identified by [frag1, frag2]
        n_occ = 0  # Number of occurences of each frag combination
        n_nonzero = 0  # Total number of nonzero matrix entries
        pairs_reader = csv.reader(pairs, delimiter=" ")
        # First line contains nrows, ncols and number of nonzero entries.
        # Number of nonzero entries is unknown for now
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
                    mat.write(
                        "\t".join(
                            map(str, [prev_pair[0], prev_pair[1], n_occ])
                        )
                        + "\n"
                    )
                prev_pair = curr_pair
                n_occ = 1
                n_nonzero += 1
        # Write the last value
        mat.write(
            "\t".join(map(str, [curr_pair[0], curr_pair[1], n_occ])) + "\n"
        )
        n_nonzero += 1
    # Edit header line to fill number of nonzero entries inplace
    with open(mat_file) as mat, open(pre_mat_file, "w") as tmp_mat:
        header = mat.readline()
        header = header.replace("-", str(n_nonzero))
        tmp_mat.write(header)
        st.copyfileobj(mat, tmp_mat)

    # Replace the matrix file with the one with corrected header
    os.rename(pre_mat_file, mat_file)


def full_pipeline(
    genome,
    input1,
    input2=None,
    enzyme=5000,
    circular=False,
    out_dir=None,
    tmp_dir=None,
    plot=False,
    min_qual=30,
    min_size=0,
    threads=1,
    no_cleanup=False,
    iterative=False,
    filter_events=False,
    prefix=None,
    start_stage="fastq",
    bedgraph=False,
    minimap2=False,
):
    """
    Run the whole hicstuff pipeline. Starting from fastq files and a genome to
    obtain a contact matrix.
    
    Parameters
    ----------
    genome : str
        Path to the genome in fasta format.
    input1 : str
        Path to the Hi-C reads in fastq format (forward), the aligned Hi-C reads
        in SAM format, or the pairs file, depending on the value of start_stage.
    input2 : str
        Path to the Hi-C reads in fastq format (forward), the aligned Hi-C reads
        in SAM format, or None, depending on the value of start_stage.
    enzyme : int or str
        Name of the enzyme used for the digestion (e.g "DpnII"). If an integer
        is used instead, the fragment attribution will be done directly using a
        fixed chunk size.
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
        Minimum fragment size required to keep a restriction fragment.
    threads : int
        Number of threads to use for parallel operations.
    no_cleanup : bool
        Whether temporary files should be deleted at the end of the pipeline.
    iterative : bool
        Use iterative mapping. Truncates and extends reads until unambiguous
        alignment.
    filter_events : bool
        Filter spurious or uninformative 3C events. Requires a restriction enzyme.
    prefix : str or None
        Choose a common name for output files instead of default GRAAL names.
    start_stage : str
        Step at which the pipeline should start. Can be "fastq", "sam" or "pairs".
    bedgraph : bool
        Use the cooler-compatible bedgraph2 format instead of GRAAL format when
        writing the matrix
    minimap2 : bool
        Use minimap2 instead of bowtie2 for read alignment.
    """
    # Pipeline can start from 3 input types
    start_time = datetime.now()
    stages = {"fastq": 0, "sam": 1, "pairs": 2}
    start_stage = stages[start_stage]

    if out_dir is None:
        out_dir = os.getcwd()

    if tmp_dir is None:
        tmp_dir = join(out_dir, "tmp")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    # Define figures output paths
    if plot:
        fig_dir = join(out_dir, "plots")
        frag_plot = join(fig_dir, "frags_hist.pdf")
        dist_plot = join(fig_dir, "event_distance.pdf")
        pie_plot = join(fig_dir, "event_distribution.pdf")
    else:
        fig_dir = None
        dist_plot = pie_plot = frag_plot = None

    # Use current time for logging and to identify files
    now = time.strftime("%Y%m%d%H%M%S")

    def _tmp_file(fname):
        if prefix:
            fname = prefix + "." + fname
        return join(tmp_dir, fname)

    def _out_file(fname):
        if prefix:
            fname = prefix + "." + fname
        return join(out_dir, fname)

    # Define temporary file names
    log_file = _out_file("hicstuff_" + now + ".log")
    sam1 = _tmp_file("for.sam")
    sam2 = _tmp_file("rev.sam")
    pairs = _tmp_file("valid.pairs")
    pairs_idx = _tmp_file("valid_idx.pairs")
    pairs_filtered = _tmp_file("valid_idx_filtered.pairs")

    # Define output file names
    if prefix:
        fragments_list = _out_file("mat.tsv")
        info_contigs = _out_file("chr.tsv")
        mat = _out_file("mat.tsv")
    else:
        # Default GRAAL file names
        fragments_list = _out_file("fragments_list.txt")
        info_contigs = _out_file("info_contigs.txt")
        mat = _out_file("abs_fragments_contacts_weighted.txt")

    # Define what input files are given
    if start_stage == 0:
        reads1, reads2 = input1, input2
    elif start_stage == 1:
        sam1, sam2 = input1, input2
    elif start_stage == 2:
        pairs_idx = input1

    # Perform genome alignment
    if start_stage == 0:
        align_reads(
            reads1,
            genome,
            sam1,
            tmp_dir=tmp_dir,
            threads=threads,
            minimap2=minimap2,
            iterative=iterative,
            min_qual=min_qual,
        )
        align_reads(
            reads2,
            genome,
            sam2,
            tmp_dir=tmp_dir,
            threads=threads,
            minimap2=minimap2,
            iterative=iterative,
            min_qual=min_qual,
        )
        # Sort alignments by read name
        ps.sort(
            "-@", str(threads), "-n", "-O", "SAM", "-o", sam1 + ".sorted", sam1
        )
        st.move(sam1 + ".sorted", sam1)
        ps.sort(
            "-@", str(threads), "-n", "-O", "SAM", "-o", sam2 + ".sorted", sam2
        )
        st.move(sam2 + ".sorted", sam2)

    if start_stage < 2:

        # Generate info_contigs and fragments_list output files
        hcd.write_frag_info(
            genome,
            enzyme,
            min_size=min_size,
            circular=circular,
            output_contigs=info_contigs,
            output_frags=fragments_list,
        )

        # Log fragment size distribution
        hcd.frag_len(
            frags_file_name=fragments_list, plot=plot, fig_path=frag_plot
        )

        # Make pairs file (readID, chr1, chr2, pos1, pos2, strand1, strand2)
        sam2pairs(sam1, sam2, pairs, info_contigs, min_qual=min_qual)

        restrict_table = {}
        for record in SeqIO.parse(genome, "fasta"):
            # Get chromosome restriction table
            restrict_table[record.id] = hcd.get_restriction_table(
                record.seq, enzyme, circular=circular
            )

        # Add fragment index to pairs (readID, chr1, pos1, chr2,
        # pos2, strand1, strand2, frag1, frag2)
        hcd.attribute_fragments(pairs, pairs_idx, restrict_table)

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

    # Build matrix from pairs.
    mat_format = "cooler" if bedgraph else "GRAAL"
    # Number of fragments is N lines in frag list - 1 for the header
    n_frags = sum(1 for line in open(fragments_list, "r")) - 1
    pairs2matrix(
        use_pairs, mat, n_frags, mat_format=mat_format, threads=threads
    )

    # Clean temporary files
    if not no_cleanup:
        tempfiles = [pairs, pairs_idx, pairs_filtered, sam1, sam2]
        # Do not delete files that were given as input
        try:
            tempfiles.remove(input1)
            tempfiles.remove(input2)
        except ValueError:
            pass
        for file in tempfiles:
            try:
                os.remove(file)
            except FileNotFoundError:
                pass

    end_time = datetime.now()
    duration = relativedelta(end_time, start_time)
    logger.info(
        "Contact map generated after {h}h {m}m {s}s".format(
            h=duration.hours, m=duration.minutes, s=duration.seconds
        )
    )