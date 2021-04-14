# Tests for CLI tools of hicstuff
# Commands are simply run to test for crashes.
# TODO: add tests to check for output contents
import hicstuff.commands as hcmd
import os
import pytest
from pathlib import Path
import shutil as su

# Use global variables for input files
GRAAL = "test_data/abs_fragments_contacts_weighted.txt"
BG2 = "test_data/mat.bg2"
COOL = "test_data/mat.cool"
FRAG = "test_data/fragments_list.txt"
CHROM = "test_data/info_contigs.txt"
OUT = "test_cli"
os.makedirs(OUT, exist_ok=True)
MATS = ("mat", [GRAAL, BG2, COOL])


@pytest.mark.parametrize(*MATS)
def test_view(mat):
    args = (
        "-b 500bp -c Reds -d -f {0} -T log2 -n -t 2 -m 0.98 -r seq1:100-50000 "
        + "-o {1}/test.png {2}"
    ).format(FRAG, OUT, mat)
    proc = hcmd.View(args.split(" "), {})
    proc.execute()


def test_pipeline():
    args = (
        "-e DpnII -t 1 -f -D -d -m iterative -n -P test -o {0} -g test_data/genome/seq "
        + "test_data/sample.reads_for.fastq.gz test_data/sample.reads_rev.fastq.gz"
    ).format(OUT)
    proc = hcmd.Pipeline(args.split(" ") + ["-F"], {})
    proc.execute()
    with pytest.raises(IOError):
        proc = hcmd.Pipeline(args.split(" "), {})
        proc.execute()


@pytest.mark.parametrize(*MATS)
def test_rebin(mat):
    args = "-b 1kb -f {0} -c {1} {2} {3}".format(
        FRAG, CHROM, mat, str(Path(OUT) / "rebinned")
    )
    proc = hcmd.Rebin(args.split(" ") + ["-F"], {})
    proc.execute()
    with pytest.raises(IOError):
        proc = hcmd.Rebin(args.split(" "), {})
        proc.execute()


def test_convert():
    args = "-f {0} -c {1} {2} {3}".format(
        FRAG, CHROM, GRAAL, str(Path(OUT) / "converted")
    )
    proc = hcmd.Convert(args.split(" ") + ["-F"], {})
    proc.execute()
    with pytest.raises(IOError):
        proc = hcmd.Convert(args.split(" "), {})
        proc.execute()


def test_distancelaw():
    args = (
        "-a -o test.png -d test_data/distance_law.txt"
        + " -c test_data/centromeres.txt -b 10000 -r 1000 -B 1.1"
    )
    proc = hcmd.Distancelaw(args.split(" "), {})
    proc.execute()


def test_distance_law_2():
    args = (
        "-p test_data/valid_idx_filtered.pairs -f {0} -C"
        + " -O {1}/test_distance_law.txt -i 500 -s 45000"
    ).format(FRAG, OUT)
    proc = hcmd.Distancelaw(args.split(" "), {})
    proc.execute()


def test_iteralign():
    args = (
        "-g test_data/genome/seq -t 1 -T tmp -l 30"
        + " -o {0}/test.sam test_data/sample.reads_for.fastq.gz"
    ).format(OUT)
    proc = hcmd.Iteralign(args.split(" "), {})
    proc.execute()


def test_digest():
    args = "-e DpnII -p -f {0} -o {0} test_data/genome/seq.fa".format(OUT)
    su.rmtree(OUT)
    proc = hcmd.Digest(args.split(" "), {})
    proc.execute()
    # Should fail, directory already exists
    with pytest.raises(IOError):
        proc = hcmd.Digest(args.split(" "), {})
        proc.execute()
    # Should succeed with --force flag
    args = "-F " + args
    proc = hcmd.Digest(args.split(" "), {})
    proc.execute()


def test_filter():
    args = "-f {0} -p test_data/valid_idx.pairs {0}/valid_idx_filtered.pairs".format(
        OUT
    )
    proc = hcmd.Filter(args.split(" "), {})
    proc.execute()


def test_scalogram():
    args = "-C viridis -n -t 1 -o {0}/scalo.png {1}".format(OUT, GRAAL)
    proc = hcmd.Scalogram(args.split(" "), {})
    proc.execute()


@pytest.mark.parametrize(*MATS)
def test_subsample(mat):
    args = "-p 0.5 {0} {1}".format(mat, str(Path(OUT) / "subsampled"))
    proc = hcmd.Subsample(args.split(" ") + ["-F"], {})
    proc.execute()
    with pytest.raises(IOError):
        proc = hcmd.Subsample(args.split(" "), {})
        proc.execute()


@pytest.mark.parametrize("mode", ["for_vs_rev", "all", "pile"])
def test_cutsite(mode):
    arg_vals = {
        "FASTQ_FOR": "test_data/sample.reads_for.fastq.gz",
        "FASTQ_REV": "test_data/sample.reads_rev.fastq.gz",
        "OUT": "test_data/digested",
    }
    args = (
        "-1 {FASTQ_FOR} -2 {FASTQ_REV} -e DpnII,HinfI -m {0} -p {OUT} -t 8"
    ).format(mode, **arg_vals)
    proc = hcmd.Cutsite(args.split(" "), {})
    proc.execute()
