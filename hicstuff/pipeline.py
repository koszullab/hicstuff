"""
Handle generation of graal-compatible contact maps from fastq files.
cmdoret, 20190322
"""
import os, time, csv, sys, re
from datetime import datetime
from dateutil.relativedelta import relativedelta
import shutil as st
import itertools
from shutil import which
import logging
from os.path import join
import subprocess as sp
from Bio import SeqIO
import pandas as pd
import numpy as np
import pysam as ps
import hicstuff.digest as hcd
import hicstuff.iteralign as hci
import hicstuff.filter as hcf
import hicstuff.io as hio
import hicstuff.distance_law as hcdl
import hicstuff.cutsite as hcc
import matplotlib
import pathlib
from hicstuff.version import __version__
import hicstuff.log as hcl
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
):
    """
    Select and call correct alignment method and generate logs accordingly.
    Alignments are filtered so that there at most one alignment per read.
    Alignments are sorted by read name in the output file.

    Parameters
    ----------
    reads : str
        Path to the fastq file with Hi-C reads.
    genome : str
        Path to the genome bowtie2/bwa index prefix if using bowtie2 or bwa, or to the 
        fasta if using minimap2.
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
    """
    if tmp_dir is None:
        tmp_dir = os.getcwd()
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
        sp.call(
            "samtools view -F 2048 -h -@ {threads} -O BAM {tmp} -o {out}".format(
                threads=threads, tmp=tmp_bam, out=out_bam
            ),
            shell=True,
        )
    else:
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
        sp.call("samtools view -F 2048 -h -@ {threads} {tmp} | samtools sort -n -@ {threads} -o {out} -".format(
                tmp=tmp_bam, threads=threads, out=out_bam
            ),
            shell=True,
        )
    os.remove(tmp_bam)


def bam2pairs(bam1, bam2, out_pairs, info_contigs, min_qual=30):
    """
    Make a .pairs file from two Hi-C bam files sorted by read names.
    The Hi-C mates are matched by read identifier. Pairs where at least one
    reads maps with MAPQ below  min_qual threshold are discarded. Pairs are
    sorted by readID and stored in upper triangle (first pair higher).

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
    """
    forward = ps.AlignmentFile(bam1, "rb")
    reverse = ps.AlignmentFile(bam2, "rb")

    # Generate header lines
    format_version = "## pairs format v1.0\n"
    sorting = "#sorted: readID\n"
    cols = "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n"
    # Chromosome order will be identical in info_contigs and pair files
    chroms = pd.read_csv(info_contigs, sep="\t").apply(
        lambda x: "#chromsize: %s %d\n" % (x.contig, x.length), axis=1
    )
    with open(out_pairs, "w") as pairs:
        pairs.writelines([format_version, sorting, cols] + chroms.tolist())
        pairs_writer = csv.writer(pairs, delimiter="\t")
        n_reads = {"total": 0, "mapped": 0}
        # Remember if some read IDs were missing from either file
        unmatched_reads = 0
        # Remember if all reads in one bam file have been read
        exhausted = [False, False]
        # Iterate on both BAM simultaneously
        end_regex = re.compile(r'/[12]$')
        for end1, end2 in itertools.zip_longest(forward, reverse):
            # Remove end-specific suffix if any
            end1.query_name = re.sub(end_regex, '', end1.query_name)
            end2.query_name = re.sub(end_regex, '', end2.query_name)
            # Both file still have reads
            # Check if reads pass filter
            try:
                end1_passed = end1.mapping_quality >= min_qual
            # Happens if end1 bam file has been exhausted
            except AttributeError:
                exhausted[0] = True
                end1_passed = False
            try:
                end2_passed = end2.mapping_quality >= min_qual
            # Happens if end2 bam file has been exhausted
            except AttributeError:
                exhausted[1] = True
                end2_passed = False
            # Skip read if mate is not present until they match or reads
            # have been exhausted
            while sum(exhausted) == 0 and end1.query_name != end2.query_name:
                # Get next read and check filters again
                # Count single-read iteration
                unmatched_reads += 1
                n_reads["total"] += 1
                if end1.query_name < end2.query_name:
                    try:
                        end1 = next(forward)
                        end1_passed = end1.mapping_quality >= min_qual
                    # If EOF is reached in BAM 1
                    except (StopIteration, AttributeError):
                        exhausted[0] = True
                        end1_passed = False
                    n_reads["mapped"] += end1_passed
                elif end1.query_name > end2.query_name:
                    try:
                        end2 = next(reverse)
                        end2_passed = end2.mapping_quality >= min_qual
                    # If EOF is reached in BAM 2
                    except (StopIteration, AttributeError):
                        exhausted[1] = True
                        end2_passed = False
                    n_reads["mapped"] += end2_passed

            # 2 reads processed per iteration, unless one file is exhausted
            n_reads["total"] += 2 - sum(exhausted)
            n_reads["mapped"] += sum([end1_passed, end2_passed])
            # Keep only pairs where both reads have good quality
            if end1_passed and end2_passed:

                # Flipping to get upper triangle
                if (
                    end1.reference_id == end2.reference_id
                    and end1.reference_start > end2.reference_start
                ) or end1.reference_id > end2.reference_id:
                    end1, end2 = end2, end1
                pairs_writer.writerow(
                    [
                        end1.query_name,
                        end1.reference_name,
                        end1.reference_start + 1,
                        end2.reference_name,
                        end2.reference_start + 1,
                        "-" if end1.is_reverse else "+",
                        "-" if end2.is_reverse else "+",
                    ]
                )
    pairs.close()
    if unmatched_reads > 0:
        logger.warning(
            "%d reads were only present in one BAM file. Make sure you sorted reads by name before running the pipeline.",
            unmatched_reads,
        )
    logger.info(
        "{perc_map}% reads (single ends) mapped with Q >= {qual} ({mapped}/{total})".format(
            total=n_reads["total"],
            mapped=n_reads["mapped"],
            perc_map=round(100 * n_reads["mapped"] / n_reads["total"]),
            qual=min_qual,
        )
    )


def generate_log_header(log_path, input1, input2, genome, enzyme):
    hcl.set_file_handler(log_path, formatter=logging.Formatter(""))
    logger.info("## hicstuff: v%s log file", __version__)
    logger.info("## date: %s", time.strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("## enzyme: %s", str(enzyme))
    logger.info("## input1: %s ", input1)
    logger.info("## input2: %s", input2)
    logger.info("## ref: %s", genome)
    logger.info("---")
    hcl.set_file_handler(log_path, formatter=hcl.logfile_formatter)


def filter_pcr_dup(pairs_idx_file, filtered_file):
    """
    Filter out PCR duplicates from a coordinate-sorted pairs file using
    overrrepresented exact coordinates. If multiple fragments have two reads
    with the exact same coordinates, only one of those fragments is kept.
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
    with open(pairs_idx_file, "r") as pairs, open(filtered_file, "w") as filtered:
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
            "%d%% PCR duplicates have been filtered out (%d / %d pairs) "
            % (100 * round(filter_count / reads_count, 3), filter_count, reads_count)
        )


def pairs2cool(pairs_file, cool_file, bins_file):
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

    # Make bins file compatible with cooler cload
    bins_tmp = bins_file + ".cooler"
    bins = pd.read_csv(bins_file, sep="\t", usecols=[1, 2, 3], skiprows=1, header=None)
    bins.to_csv(bins_tmp, sep="\t", header=False, index=False)

    cooler_cmd = "cooler cload pairs -c1 2 -p1 3 -p2 4 -c2 5 {bins} {pairs} {cool}"
    cool_args = {"bins": bins_tmp, "pairs": pairs_file, "cool": cool_file}
    sp.call(cooler_cmd.format(**cool_args), shell=True)
    os.remove(bins_tmp)


def pairs2matrix(
    pairs_file, mat_file, fragments_file, mat_fmt="graal", threads=1, tmp_dir=None
):
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
    n_frags = sum(1 for line in open(fragments_file, "r")) - 1
    frags = pd.read_csv(fragments_file, delimiter="\t")

    def write_mat_entry(frag1, frag2, contacts):
        """Write a single sparse matrix entry in either graal or bg2 format"""
        if mat_fmt == "graal":
            mat.write("\t".join(map(str, [frag1, frag2, n_occ])) + "\n")
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
    with open(pre_mat_file, "r") as pairs, open(mat_file, "w") as mat:
        # Skip header lines
        for _ in range(header_size):
            next(pairs)
        prev_pair = ["0", "0"]  # Pairs identified by [frag1, frag2]
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
                prev_pair = curr_pair
                n_pairs += n_occ
                n_occ = 1
                n_nonzero += 1
        # Write the last value
        write_mat_entry(curr_pair[0], curr_pair[1], n_occ)
        n_nonzero += 1
        n_pairs += 1

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


def full_pipeline(
    genome,
    input1,
    input2=None,
    aligner="bowtie2",
    centromeres=None,
    circular=False,
    distance_law=False,
    enzyme=5000,
    filter_events=False,
    force=False,
    mapping="normal",
    mat_fmt="graal",
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
):
    """
    Run the whole hicstuff pipeline. Starting from fastq files and a genome to
    obtain a contact matrix.

    Parameters
    ----------
    genome : str
        Path to the bowtie2/bwa index prefix if using bowtie2/bwa or to the genome 
        in fasta format if using minimap2.
    input1 : str
        Path to the Hi-C reads in fastq format (forward), the aligned Hi-C reads
        in BAM format, or the pairs file, depending on the value of start_stage.
    input2 : str
        Path to the Hi-C reads in fastq format (forward), the aligned Hi-C reads
        in BAM format, or None, depending on the value of start_stage.
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
        Minimum contig size required to keep it.
    threads : int
        Number of threads to use for parallel operations.
    no_cleanup : bool
        Whether temporary files should be deleted at the end of the pipeline.
    mapping : str
        normal|iterative|cutsite. Use normal, iterative or cutsite mapping. 
        "normal": Normal alignement. "iterative": Truncates and extends reads 
        until unambiguous alignment. "cutsite": Digest reads at religation sites 
        and build new pairs from the fragments created.
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
        bedgraph2 format, "cool" for Mirnylab's cool format, or graal for a
        plain text COO format compatible with Koszullab's instagraal software.
    aligner : str
        Read alignment software to use. Can be either "minimap2", "bwa" or "bowtie2".
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
    read_len : int
        Maximum read length to expect in the fastq file. Optionally used in iterative
        alignment mode. Estimated from the first read by default. Useful if input fastq
        is a composite of different read lengths.
    remove_centros : None or int
        If the distance law is computed, this is the number of kb that will be removed
        around the centromere position given by in the centromere file.
    """
    # Check if third parties can be run
    if aligner in ("bowtie2", "minimap2", "bwa"):
        if check_tool(aligner) is None:
            logger.error("%s is not installed or not on PATH", aligner)
            sys.exit(1)
    else:
        logger.error("Incompatible aligner software, choose bowtie2, minimap2 or bwa")
        sys.exit(1)
    if check_tool("samtools") is None:
        logger.error("Samtools is not installed or not on PATH")
        sys.exit(1)

    # Pipeline can start from 3 input types
    start_time = datetime.now()
    stages = {"fastq": 0, "bam": 1, "pairs": 2, "pairs_idx": 3}
    start_stage = stages[start_stage]

    # Check if the number of input files is correct
    if start_stage <= 1:
        if input2 is None:
            logger.error(
                "You must provide 2 input files when --start-stage is fastq " "or bam."
            )
            sys.exit(1)
    else:
        if input2 is not None:
            logger.error(
                "You must provide a single input file when --start-stage is "
                "pairs or pairs_idx."
            )
            sys.exit(1)
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
            raise IOError(
                "Temporary file {} already exists. Use --force to overwrite".format(
                    full_path
                )
            )
        return full_path

    def _out_file(fname):
        if prefix:
            fname = prefix + "." + fname
        full_path = join(out_dir, fname)
        if not force and os.path.exists(full_path):
            raise IOError(
                "Output file {} already exists. Use --force to overwrite".format(
                    full_path
                )
            )

        return full_path

    # Define temporary file names
    log_file = _out_file("hicstuff_" + now + ".log")
    tmp_genome = _tmp_file("genome.fasta")
    bam1 = _tmp_file("for.bam")
    bam2 = _tmp_file("rev.bam")
    pairs = _tmp_file("valid.pairs")
    pairs_idx = _tmp_file("valid_idx.pairs")
    pairs_filtered = _tmp_file("valid_idx_filtered.pairs")
    pairs_pcr = _tmp_file("valid_idx_pcrfree.pairs")

    # Enable file logging
    hcl.set_file_handler(log_file)
    generate_log_header(log_file, input1, input2, genome, enzyme)

    # If the user chose bowtie2 and supplied an index, extract fasta from it
    # For later steps of the pipeline (digestion / frag attribution)
    # Check if the genome is an index or fasta file
    idx = hio.check_fasta_index(genome, mode=aligner)
    is_fasta = hio.check_is_fasta(genome)
    
    # Different aligners accept different files. Make sure the input format is good.
    # Note bowtie2 can extract fasta from the index, but bwa cannot
    sane_input = {
            'bowtie2': is_fasta or idx,
            'minimap2': is_fasta, 
            'bwa': is_fasta
    }

    if not sane_input[aligner]:
        logger.error("You must provide either a fasta or bowtie2 index prefix as genome")

    # Just use the input genome is it is indexed
    if is_fasta and idx:
        fasta = genome
    # Otherwise copy it in tmpdir for indexing, unless the input is a bt2 index, in which
    # case fasta will be extracted later from it.
    else:
        if is_fasta:
            st.copy(genome, tmp_genome)
            genome = tmp_genome
        fasta = tmp_genome
        

    # Bowtie2-specific feature: extract fasta from the index
    if aligner == 'bowtie2' and not is_fasta:
        # Index is present, extract fasta file from it
        bt2fa = sp.Popen(
            ["bowtie2-inspect", genome],
            stdout=open(tmp_genome, "w"),
            stderr=sp.PIPE,
        )
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

    # Build index with bowtie2 / bwa if required
    if idx is None and aligner in ['bowtie2', 'bwa']:
        if aligner == 'bowtie2':
            index_cmd = ["bowtie2-build", '-q', fasta, fasta]
        elif aligner == 'bwa':
            index_cmd = ['bwa', 'index', fasta]
        # We only need the index if the user provided fastq input
        if start_stage == 0:
            # If no index present assume input is fasta, copy it in tmp and
            # index it (to avoid conflict between instances)
            logger.info(
                "%s index not found at %s, generating "
                "a local temporary index.", aligner, genome
            )
            sp.run(index_cmd, stderr=sp.PIPE)

    # Check for spaces in fasta headers and issue error if found
    for record in SeqIO.parse(fasta, "fasta"):
        if " " in record.id:
            logger.error(
                "Sequence identifiers contain spaces. Please clean the input genome."
            )
    # Define output file names (tsv files)
    if prefix:
        fragments_list = _out_file("frags.tsv")
        info_contigs = _out_file("chr.tsv")
        mat = _out_file("mat.tsv")
        # If matrix has a different format, give it the right extension
        if mat_fmt != "graal":
            mat = _out_file(mat_fmt)
    else:
        # Default graal file names
        fragments_list = _out_file("fragments_list.txt")
        info_contigs = _out_file("info_contigs.txt")
        mat = _out_file("abs_fragments_contacts_weighted.txt")
        if mat_fmt != "graal":
            mat = _out_file("abs_fragments_contacts_weighted." + mat_fmt)
    # Define what input files are given
    if start_stage == 0:
        reads1, reads2 = input1, input2
    elif start_stage == 1:
        bam1, bam2 = input1, input2
    elif start_stage == 2:
        pairs = input1
    elif start_stage == 3:
        pairs_idx = input1
 
    # Perform genome alignment
    if start_stage == 0:
        
        # Define mapping choice (default normal):
        if mapping == "normal":
            iterative = False
        elif mapping == "iterative":
            iterative = True   
        elif mapping == "cutsite":
            # If no enzyme given use iterative alignment.
            try:
                int(enzyme)
                logger.warning("No enzyme has been given. Can't map using cutsite, iterative mapping will be used instead.")
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
        else:
            logger.error("mapping must be either normal, iterative or cutsite.")
            raise ValueError
        
        align_reads(
            reads1,
            genome,
            bam1,
            tmp_dir=tmp_dir,
            threads=threads,
            aligner=aligner,
            iterative=iterative,
            min_qual=min_qual,
            read_len=read_len,
        )
        align_reads(
            reads2,
            genome,
            bam2,
            tmp_dir=tmp_dir,
            threads=threads,
            aligner=aligner,
            iterative=iterative,
            min_qual=min_qual,
            read_len=read_len,
        )

    # Detect if multiple enzymes are given
    if re.search(",", enzyme):
        enzyme = enzyme.split(",")
        
    # Starting from bam files
    if start_stage <= 1:

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

        # Make pairs file (readID, chr1, chr2, pos1, pos2, strand1, strand2)
        bam2pairs(bam1, bam2, pairs, info_contigs, min_qual=min_qual)

    # Starting from pairs file
    if start_stage <= 2:
        restrict_table = {}
        for record in SeqIO.parse(fasta, "fasta"):
            # Get chromosome restriction table
            restrict_table[record.id] = hcd.get_restriction_table(
                record.seq, enzyme, circular=circular
            )

        # Add fragment index to pairs (readID, chr1, pos1, chr2,
        # pos2, strand1, strand2, frag1, frag2)
        hcd.attribute_fragments(pairs, pairs_idx, restrict_table)

    # Sort pairs file by coordinates for next steps
    hio.sort_pairs(
        pairs_idx,
        pairs_idx + ".sorted",
        keys=["chr1", "pos1", "chr2", "pos2"],
        threads=threads,
        tmp_dir=tmp_dir,
    )
    os.rename(pairs_idx + ".sorted", pairs_idx)

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
        # Name matrix file in .cool
        cool_file = os.path.splitext(mat)[0] + ".cool"
        pairs2cool(use_pairs, cool_file, fragments_list)
    else:
        pairs2matrix(
            use_pairs,
            mat,
            fragments_list,
            mat_fmt=mat_fmt,
            threads=threads,
            tmp_dir=tmp_dir,
        )

    # Clean temporary files
    if not no_cleanup:
        tempfiles = [
            pairs,
            pairs_idx,
            pairs_filtered,
            bam1,
            bam2,
            pairs_pcr,
            tmp_genome,
        ]
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
