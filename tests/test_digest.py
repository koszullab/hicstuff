# Tests for the hicstuff digest module
# 20190402
import filecmp
import os
from os.path import join
from tempfile import NamedTemporaryFile

import pandas as pd
from Bio import SeqIO

from hicstuff import digest as hcd


def test_write_frag_info():
    """Test generation of fragments_list.txt and info_contigs.txt"""
    genome = NamedTemporaryFile(delete=False, mode="w")
    seq = "GGAATAGATCAAATGATCCACAGATC"
    genome.write(">seq1\n")
    genome.write(seq)
    genome.close()
    out_dir, tigs, frags = "test_data", "test_tigs", "test_frags"
    hcd.write_frag_info(
        genome.name,
        "DpnII",
        output_contigs=tigs,
        output_frags=frags,
        output_dir=out_dir,
    )
    tigs_df = pd.read_csv(join(out_dir, tigs), delimiter="\t")
    frags_df = pd.read_csv(join(out_dir, frags), delimiter="\t")

    assert tigs_df.length.tolist()[0] == len(seq)
    assert frags_df.start_pos.tolist() == [0, 6, 14, 22]

    os.unlink(genome.name)
    os.remove(join(out_dir, tigs))
    os.remove(join(out_dir, frags))


def test_attribute_fragments():
    """Test the attribution of reads to restriction fragments"""
    idx_pairs = NamedTemporaryFile(delete=False)
    idx_pairs.close()
    restriction_table = {}
    for record in SeqIO.parse("test_data/genome/seq.fa", "fasta"):
        # Get chromosome restriction table
        restriction_table[record.id] = hcd.get_restriction_table(record.seq, "DpnII")
    hcd.attribute_fragments("test_data/valid.pairs", idx_pairs.name, restriction_table)

    assert filecmp.cmp("test_data/valid_idx.pairs", idx_pairs.name)
