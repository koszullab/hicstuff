#!/usr/bin/env python3
# coding: utf-8
"""Iterative alignment

Aligns iteratively reads from a 3C fastq file: reads
are trimmed with a range-sweeping number of basepairs
and each read set generated this way is mapped onto
the reference genome. This may result in a small
increase of properly mapped reads.

@author: Remi Montagne & cmdoret
"""

import os
import sys
import glob
import re
import subprocess as sp
import pysam as ps
import shutil as st
import hicstuff.io as hio
import contextlib
from hicstuff.log import logger
from os.path import join




def iterative_align(
    fq_in,
    tmp_dir,
    ref,
    n_cpu,
    bam_out,
    aligner="bowtie2",
    min_len=20,
    min_qual=30,
    read_len=None,
):
    """Iterative alignment

    Aligns reads iteratively reads of fq_in with bowtie2, minimap2 or bwa. Reads are
    truncated to the 20 first nucleotides and unmapped reads are extended by 20
    nucleotides and realigned on each iteration.

    Parameters
    ----------

    fq_in : str
        Path to input fastq file to align iteratively.
    tmp_dir : str
        Path where temporary files should be written.
    ref : str
        Path to the reference genome if Minimap2 is used for alignment.
        Path to the index genome if Bowtie2/bwa is used for alignment. 
    n_cpu : int
        The number of CPUs to use for the iterative alignment.
    bam_out : str
        Path where the final alignment should be written in BAM format.
    aligner : str
        Choose between minimap2, bwa or bowtie2 for the alignment.
    min_len : int
        The initial length of the fragments to align.
    min_qual : int
        Minimum mapping quality required to keep Hi-C pairs.
    read_len : int
        Read length in the fasta file. If set to None, the length of the first read
        is used. Set this value to the longest read length in the file if you have
        different read lengths.
        
    Examples
    --------

    iterative_align(fq_in='example_for.fastq', ref='example_bt2_index', bam_out='example_for.bam', aligner="bowtie2")
    iterative_align(fq_in='example_for.fastq', ref='example_genome.fa', bam_out='example_for.bam', aligner="minimap2")
    """
    # set with the name of the unaligned reads :
    remaining_reads = set()
    total_reads = 0
    # Store path of SAM containing aligned reads at each iteration.
    iter_out = []

    # If there is already a file with the same name as the output file,
    # remove it. Otherwise, ignore.
    with contextlib.suppress(FileNotFoundError):
        try:
            os.remove(bam_out)
        except IsADirectoryError:
            logger.error("You need to give the BAM output file, not a folder.")
            raise

    # Bowtie only accepts uncompressed fastq: uncompress it into a temp file
    if aligner == "bowtie2" and hio.is_compressed(fq_in):
        uncomp_path = join(tmp_dir, os.path.basename(fq_in) + ".tmp")
        with hio.read_compressed(fq_in) as inf:
            with open(uncomp_path, "w") as uncomp:
                st.copyfileobj(inf, uncomp)
    else:
        uncomp_path = fq_in

    # throw error if index does not exist
    index = hio.check_fasta_index(ref, mode=aligner)
    if index is None:
        logger.error(
            "Reference index is missing, please build the {} ".format(aligner), 
            "index first."
        )
        sys.exit(1)
    # Counting reads
    with hio.read_compressed(uncomp_path) as inf:
        for _ in inf:
            total_reads += 1
    total_reads /= 4

    # Use first read to guess read length if not provided.
    if read_len is None:
        with hio.read_compressed(uncomp_path) as inf:
            # Skip first line (read header)
            size = inf.readline()
            # Stripping newline from sequence line.
            read_len = len(inf.readline().rstrip())

    # initial length of the fragments to align
    # In case reads are shorter than provided min_len
    if read_len > min_len:
        n = min_len
    else:
        logger.warning(
            "min_len is longer than the reads. Iterative mapping will have no effect."
        )
        n = read_len
    logger.info("{0} reads to parse".format(int(total_reads)))

    first_round = True
    # iterative alignment per se
    while n <= read_len:
        logger.info(
            "Truncating unaligned reads to {size}bp and mapping{again}.".format(
                size=int(n), again="" if first_round else " again"
            )
        )
        iter_out += [join(tmp_dir, "trunc_{0}.bam".format(str(n)))]
        # Generate a temporary input fastq file with the n first nucleotids
        # of the reads.
        truncated_reads = truncate_reads(
            tmp_dir, uncomp_path, remaining_reads, n, first_round
        )

        # Align the truncated reads on reference genome
        temp_alignment = join(tmp_dir, "temp_alignment.bam")
        map_args = {
            "fa": ref,
            "cpus": n_cpu,
            "fq": truncated_reads,
            "idx": index,
            "bam": temp_alignment,
        }
        if re.match(r"^(minimap[2]?|mm[2]?)$", aligner, flags=re.IGNORECASE):
            cmd = "minimap2 -x sr -a -t {cpus} {fa} {fq}".format(**map_args)
        elif re.match(r"^(bwa)$", aligner, flags=re.IGNORECASE):
            cmd = "bwa mem -t {cpus} -v 1 {idx} {fq}".format(**map_args)
        elif re.match(r"^(bowtie[2]?|bt[2]?)$", aligner, flags=re.IGNORECASE):
            cmd = (
                "bowtie2 -x {idx} -p {cpus}" " --quiet --very-sensitive {fq}"
            ).format(**map_args)
        else:
            raise ValueError("Unknown aligner. Select bowtie2, minimap2 or bwa.")

        map_process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
        sort_process = sp.Popen(
            "samtools sort -n -@ {cpus} -O BAM -o {bam}".format(**map_args),
            shell=True,
            stdin=map_process.stdout,
        )
        out, err = sort_process.communicate()

        # filter the reads: the reads whose truncated end was aligned are written
        # to the output file.
        # The reads whose truncated end was not aligned are kept for the next round.
        remaining_reads = filter_bamfile(temp_alignment, iter_out[-1], min_qual)

        n += 20
        first_round = False

    # one last round without trimming
    logger.info(
        "Trying to map unaligned reads at full length ({0}bp).".format(int(read_len))
    )

    truncated_reads = truncate_reads(
        tmp_dir,
        infile=uncomp_path,
        unaligned_set=remaining_reads,
        trunc_len=n,
        first_round=first_round,
    )
    if aligner == "minimap2" or aligner == "Minimap2":
        cmd = "minimap2 -x sr -a -t {cpus} {fa} {fq}".format(
            fa=ref, cpus=n_cpu, fq=truncated_reads
        )
    elif aligner == "bwa" or aligner == "Bwa" or aligner == "BWA":
        cmd = "bwa mem -v 1 -t {cpus} {idx} {fq}".format(
            idx=index, cpus=n_cpu, fq=truncated_reads
        )
    else:
        cmd = ("bowtie2 -x {idx} -p {cpus} --quiet " "--very-sensitive {fq}").format(
            idx=index, cpus=n_cpu, fq=truncated_reads
        )
    map_process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    # Keep reads sorted by name
    sort_process = sp.Popen(
        "samtools sort -n -@ {cpus} -O BAM -o {bam}".format(
            cpus=n_cpu, bam=temp_alignment
        ),
        shell=True,
        stdin=map_process.stdout,
    )
    out, err = sort_process.communicate()
    iter_out += [join(tmp_dir, "trunc_{0}.bam".format(str(n)))]
    remaining_reads = filter_bamfile(temp_alignment, iter_out[-1], min_qual)

    # Report unaligned reads as well
    iter_out += [join(tmp_dir, "unaligned.bam")]
    temp_bam = ps.AlignmentFile(temp_alignment, "rb", check_sq=False)
    unmapped = ps.AlignmentFile(iter_out[-1], "wb", template=temp_bam)
    for r in temp_bam:
        # Do not write supplementary alignments (keeping 1 alignment/read)
        if r.query_name in remaining_reads and not r.is_supplementary:
            unmapped.write(r)
    unmapped.close()
    temp_bam.close()

    # Merge all aligned reads and unmapped reads into a single bam
    ps.merge("-n", "-O", "BAM", "-@", str(n_cpu), bam_out, *iter_out)
    logger.info(
        "{0} reads aligned / {1} total reads.".format(
            int(total_reads - len(remaining_reads)), int(total_reads)
        )
    )

    return 0


def truncate_reads(tmp_dir, infile, unaligned_set, trunc_len, first_round):
    """Trim read ends

    Writes the n first nucleotids of each sequence in infile to an auxiliary.
    file in the temporary folder.

    Parameters
    ----------

    tmp_dir : str
        Path to the temporary folder.
    infile : str
        Path to the fastq file to truncate.
    unaligned_set : set
        Contains the names of all reads that did not map unambiguously in
        previous rounds.
    trunc_len : int
        The number of basepairs to keep in each truncated sequence.
    first_round : bool
        If this is the first round, truncate all reads without checking mapping.
    
    Returns
    -------

    str :
        Path to the output fastq file containing truncated reads.
    """

    outfile = "{0}/truncated.fastq".format(tmp_dir)
    with ps.FastxFile(infile, "r") as inf, open(outfile, "w") as outf:
        for entry in inf:
            # If the read did not align in previous round or this is the first round
            if (entry.name in unaligned_set) or first_round:
                entry.sequence = entry.sequence[:trunc_len]
                entry.quality = entry.quality[:trunc_len]
                outf.write(str(entry) + "\n")
    return outfile


def filter_bamfile(temp_alignment, filtered_out, min_qual=30):
    """Filter alignment BAM files

    Reads all the reads in the input BAM alignment file.
    Write reads to the output file if they are aligned with a good
    quality, otherwise add their name in a set to stage them for the next round
    of alignment.

    Parameters
    ----------
    temp_alignment : str
        Path to the input temporary alignment.
    outfile : str
        Path to the output filtered temporary alignment.
    min_qual : int
        Minimum mapping quality required to keep a Hi-C pair.
    
    Returns
    -------
    set:
        Contains the names reads that did not align.
    """
    # Check the quality and status of each aligned fragment.
    # Write the ones with good quality in the final output file.
    # Keep those that do not map unambiguously for the next round.

    unaligned = set()
    temp_bam = ps.AlignmentFile(temp_alignment, "rb", check_sq=False)
    outf = ps.AlignmentFile(filtered_out, "wb", template=temp_bam)
    for r in temp_bam:
        if r.flag in [0, 16] and r.mapping_quality >= min_qual:
            outf.write(r)
        else:
            unaligned.add(r.query_name)

    logger.info("{0} reads left to map.".format(len(unaligned)))
    temp_bam.close()
    outf.close()

    return unaligned
