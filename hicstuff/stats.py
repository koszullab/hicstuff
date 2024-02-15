#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re 
from os.path import join
import json

def get_pipeline_stats(prefix, out_dir, log_file):
    """Get stats after pipeline execution.

    Parameters
    ----------
    prefix : str
        The prefix used to create output files by the pipeline.
    out_dir : str
        The prefix used to create output files by the pipeline.
    log_file : str
        Path to hicstuff log file.

    Returns
    -------
    txt file:
        A single text file containing stats about hicstuff pipeline execution: 
        - Number of sequenced pairs from fastqs
        - Number (% of total reads) of unmapped reads 
        - Number (% of total reads) of mapped reads 
        - Number (% of total reads) of pairs 
        - Number (% of total pairs) of filtered pairs 
            - Number (% of total pairs) of loops
            - Number (% of total pairs) of uncuts
            - Number (% of total pairs) of abnormal (-- and ++)
        - Number (% of total pairs) of deduplicated pairs [Number (% of total pairs) of PCR duplicates]
        - From all pairs used in contact matrix: 
            - Number (% of matrix) of cis pairs
            - Number (% of matrix) of trans pairs
            - Trans ratio
    """

    with open(log_file) as file:
        log_lines = [line.rstrip() for line in file]

    # 1. Number of sequenced pairs from fastqs
    fastq_pairs = [s for s in log_lines if re.search("reads found in eajch fastq file.", s)][0]
    fastq_pairs = re.sub(".*INFO :: ", "", fastq_pairs)
    fastq_pairs = re.sub(" reads found in each fastq file.*", "", fastq_pairs)
    fastq_pairs = int(float(fastq_pairs))
    
    # 2. Number (% of total) of (un)mapped reads
    tot_mapped = [s for s in log_lines if re.search("mapped with Q ", s)][0]
    tot_mapped = re.sub(".*Q >= 30 \(", "", tot_mapped)
    tot_mapped = re.sub("/.*", "", tot_mapped)
    tot_mapped = int(float(tot_mapped))
    tot_unmapped = fastq_pairs*2 - tot_mapped
    pct_mapped = round(tot_mapped/(fastq_pairs*2)*100, 2)
    pct_unmapped = round(tot_unmapped/(fastq_pairs*2)*100, 2)

    # 3. Number (% of total) of pairs
    tot_pairs = [s for s in log_lines if re.search("pairs successfully mapped", s)][0]
    tot_pairs = re.sub(".*INFO :: ", "", tot_pairs)
    tot_pairs = re.sub(" pairs successfully.*", "", tot_pairs)
    tot_pairs = int(float(tot_pairs))
    pct_pairs = round(tot_pairs/fastq_pairs*100, 2)
    
    # 4. Number (% of total) of filtered pairs
    filtered_pairs = [s for s in log_lines if re.search("pairs discarded:", s)]
    if (len(filtered_pairs) > 0): 
        filtered_pairs = filtered_pairs[0]
        loops_pairs = re.sub(".*Loops: ", "", filtered_pairs)
        loops_pairs = re.sub(", Uncuts:.*", "", loops_pairs)
        loops_pairs = int(float(loops_pairs))
        uncuts_pairs = re.sub(".*Uncuts: ", "", filtered_pairs)
        uncuts_pairs = re.sub(", Weirds:.*", "", uncuts_pairs)
        uncuts_pairs = int(float(uncuts_pairs))
        abnormal_pairs = re.sub(".*Weirds: ", "", filtered_pairs)
        abnormal_pairs = int(float(abnormal_pairs))
        filtered_pairs = re.sub(".*INFO :: ", "", filtered_pairs)
        filtered_pairs = re.sub(" pairs discarded.*", "", filtered_pairs)
        filtered_pairs = int(float(filtered_pairs))
    else: 
        loops_pairs=0
        uncuts_pairs=0
        abnormal_pairs=0
        filtered_pairs=0
    pct_filtered = round(filtered_pairs/tot_pairs*100, 2)
    pct_loops_pairs = round(loops_pairs/tot_pairs*100, 2)
    pct_uncuts_pairs = round(uncuts_pairs/tot_pairs*100, 2)
    pct_abnormal_pairs = round(abnormal_pairs/tot_pairs*100, 2)

    # 5. Number (% of total) of PCR dups pairs
    PCR_pairs = [s for s in log_lines if re.search("PCR duplicates have been filtered", s)]
    if (len(PCR_pairs) > 0): 
        PCR_pairs = PCR_pairs[0]
        PCR_pairs = re.sub(".*have been filtered out \(", "", PCR_pairs)
        PCR_pairs = re.sub(" / .*", "", PCR_pairs)
        PCR_pairs = int(float(PCR_pairs))
    else: 
        PCR_pairs = 0
    pct_PCR = round(PCR_pairs/tot_pairs*100, 2)

    # 6. Number (%) of final pairs
    kept_pairs=tot_pairs-filtered_pairs-PCR_pairs
    pct_kept = round(kept_pairs/tot_pairs*100, 2)

    # Open the log file and read its contents
    stats_file_path = join(out_dir, prefix + ".stats.txt")
    with open(stats_file_path, 'w') as stats_file:
        stats_file.write("Sample:             {prefix}\n".format(prefix = prefix))
        stats_file.write("Total read pairs:   {reads}\n".format(reads = "{:,}".format(fastq_pairs)))
        stats_file.write("Mapped reads:       {tot_mapped} ({pct_mapped}%  of all reads)\n".format(
                tot_mapped = "{:,}".format(tot_mapped), 
                pct_mapped = pct_mapped
            )
        )
        stats_file.write("Unmapped reads:     {tot_unmapped} ({pct_unmapped}%  of all reads)\n".format(
            tot_unmapped = "{:,}".format(tot_unmapped), pct_unmapped = pct_unmapped
        ))
        stats_file.write("Recovered contacts: {tot_pairs} ({pct_pairs}%  of all read pairs)\n".format(
            tot_pairs = "{:,}".format(tot_pairs), pct_pairs = pct_pairs
        ))
        stats_file.write("Final contacts:     {kept_pairs} ({pct_kept}% of all contacts)\n".format(
            kept_pairs = "{:,}".format(kept_pairs), pct_kept = pct_kept
        ))
        stats_file.write("  Filtered out:     {filtered_pairs} ({pct_filtered}% of all contacts)\n".format(
            filtered_pairs = "{:,}".format(filtered_pairs), pct_filtered = pct_filtered
        ))
        stats_file.write("    Loops:          {loops_pairs} ({pct_loops_pairs}% of all contacts)\n".format(
            loops_pairs = "{:,}".format(loops_pairs), pct_loops_pairs = pct_loops_pairs
        ))
        stats_file.write("    Uncuts:         {uncuts_pairs} ({pct_uncuts_pairs}% of all contacts)\n".format(
            uncuts_pairs = "{:,}".format(uncuts_pairs), pct_uncuts_pairs = pct_uncuts_pairs
        ))
        stats_file.write("    Weirds:         {abnormal_pairs} ({pct_abnormal_pairs}% of all contacts)\n".format(
            abnormal_pairs = "{:,}".format(abnormal_pairs), pct_abnormal_pairs = pct_abnormal_pairs
        ))
        stats_file.write("  PCR duplicates:   {PCR_pairs} ({pct_PCR}% of all contacts)\n".format(
            PCR_pairs = "{:,}".format(PCR_pairs), pct_PCR = pct_PCR
        ))

    stats = {
        'Total read pairs': fastq_pairs,
        'Mapped reads': tot_mapped,
        'Unmapped reads': tot_unmapped,
        'Recovered contacts': tot_pairs,
        'Final contacts': kept_pairs,
        'Filtered out': filtered_pairs,
        'Loops': loops_pairs,
        'Uncuts': uncuts_pairs,
        'Weirds': abnormal_pairs,
        'PCR duplicates': PCR_pairs
    }
    stats_json_path = join(out_dir, prefix + ".stats.json")
    with open(stats_json_path, 'w') as json_file:
        json.dump(stats, json_file, indent=4)

