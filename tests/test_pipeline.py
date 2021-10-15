# Test functions for the pipeline submodule

from tempfile import NamedTemporaryFile
import os, shutil
import pandas as pd
import pytest
import filecmp
import numpy as np
import hicstuff.pipeline as hpi


MAPPING_PARAMETERS = ("mapping", ['normal', 'iterative'])
ALIGNER_PARAMETERS = ("aligner", ['bowtie2', 'bwa', 'minimap2'])


def test_sam2pairs():
    ...


def test_pairs2mat():
    ...


def test_filter_pcr_dup():
    """Test if PCR duplicates are removed correctly"""
    dup_rm = NamedTemporaryFile(mode="w", delete=False)
    dup_rm.close()
    test_pairs = NamedTemporaryFile(mode="w", delete=False)
    test_rm = NamedTemporaryFile(mode="w", delete=False)
    lnum = 0
    # Copy the test valid_idx file, but generate PCR dups of the pair at line 50
    with open("test_data/valid_idx.pairs", "r") as pairs:
        for line in pairs:
            test_pairs.write(line)
            if lnum == 50:
                # Making 30 duplicates of this pair
                for i in range(30):
                    test_pairs.write(line)
            lnum += 1
    test_pairs.close()

    # Remove duplicates from the original pairs file and from the artificially
    # amplified one
    hpi.filter_pcr_dup("test_data/valid_idx.pairs", dup_rm.name)
    hpi.filter_pcr_dup(test_pairs.name, test_rm.name)

    # Check if duplicates have been removed correctly (both files are identical
    # after PCR filter)
    assert filecmp.cmp(dup_rm.name, test_rm.name)
    os.unlink(test_pairs.name)
    os.unlink(dup_rm.name)
    os.unlink(test_rm.name)


def test_full_pipeline_frags():
    mapping_to_enzyme = {
        'normal': "DpnII",
        'cutsite': "DpnII,HpaII"
    }
    for mapping, enzyme in mapping_to_enzyme.items():
        hpi.full_pipeline(
            input1="test_data/sample.reads_for.fastq.gz",
            input2="test_data/sample.reads_rev.fastq.gz",
            genome="test_data/genome/seq",
            enzyme=enzyme,
            mapping=mapping,
            out_dir="test_out",
            plot=True,
            pcr_duplicates=True,
            filter_events=True,
            no_cleanup=True,
            force=True,
        )


@pytest.mark.parametrize(*MAPPING_PARAMETERS)
@pytest.mark.parametrize(*ALIGNER_PARAMETERS)
def test_full_pipeline_bin(mapping, aligner):
    """Crash Test for the whole pipeline"""
    start_input = {
        'fastq': [
            "test_data/sample.reads_for.fastq.gz",
            "test_data/sample.reads_rev.fastq.gz",
        ],
        'bam': ['test_out/tmp/for.bam', 'test_out/tmp/rev.bam'],
        'pairs': ['test_out/tmp/valid.pairs', None],
        'pairs_idx': ['test_out/tmp/valid_idx.pairs', None]
    }
    for stage, [in1, in2] in start_input.items():
        # Indexed or non-indexed genome
        # for genome in ['test_data/genome/seq', 'test_data/genome/seq.fa']:
        hpi.full_pipeline(
            input1=in1,
            input2=in2,
            genome="test_data/genome/seq.fa",
            enzyme=5000,
            out_dir="test_out2",
            aligner=aligner,
            mapping=mapping,
            prefix="test",
            distance_law=True,
            start_stage=stage,
            mat_fmt="cooler",
            force=True,
        )
    shutil.rmtree("test_out2/")
    if mapping == 'iterative' and aligner == 'minimap2':
        shutil.rmtree("test_out/")
