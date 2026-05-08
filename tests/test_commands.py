# Tests for CLI tools of hicstuff
# Commands are run via click.testing.CliRunner to test for crashes.
import os
import shutil as su
from pathlib import Path

import pytest
from click.testing import CliRunner

from hicstuff.cli import cli

# Use global variables for input files
GRAAL = "test_data/abs_fragments_contacts_weighted.txt"
BG2 = "test_data/mat.bg2"
COOL = "test_data/mat.cool"
FRAG = "test_data/fragments_list.txt"
CHROM = "test_data/info_contigs.txt"
OUT = "test_cli"
os.makedirs(OUT, exist_ok=True)
MATS = ("mat", [GRAAL, BG2, COOL])


@pytest.fixture
def runner():
    return CliRunner()


@pytest.mark.parametrize(*MATS)
def test_view(runner, mat):
    result = runner.invoke(
        cli,
        [
            "view",
            "-b",
            "500bp",
            "-c",
            "Reds",
            "-d",
            "-f",
            FRAG,
            "-T",
            "log2",
            "-n",
            "-t",
            "2",
            "-m",
            "0.98",
            "-r",
            "seq1:100-50000",
            "-o",
            f"{OUT}/test.png",
            mat,
        ],
    )
    assert result.exit_code == 0, result.output


def test_pipeline(runner):
    result = runner.invoke(
        cli,
        [
            "pipeline",
            "-e",
            "DpnII",
            "-t",
            "1",
            "-f",
            "-D",
            "-d",
            "-m",
            "iterative",
            "-n",
            "-P",
            "test",
            "-o",
            OUT,
            "-g",
            "test_data/genome/seq",
            "-F",
            "test_data/sample.reads_for.fastq.gz",
            "test_data/sample.reads_rev.fastq.gz",
        ],
    )
    assert result.exit_code == 0, result.output
    # Should fail without --force when output already exists
    result = runner.invoke(
        cli,
        [
            "pipeline",
            "-e",
            "DpnII",
            "-t",
            "1",
            "-f",
            "-D",
            "-d",
            "-m",
            "iterative",
            "-n",
            "-P",
            "test",
            "-o",
            OUT,
            "-g",
            "test_data/genome/seq",
            "test_data/sample.reads_for.fastq.gz",
            "test_data/sample.reads_rev.fastq.gz",
        ],
    )
    assert result.exit_code != 0


@pytest.mark.parametrize(*MATS)
def test_rebin(runner, mat):
    out_prefix = str(Path(OUT) / "rebinned")
    result = runner.invoke(
        cli,
        ["rebin", "-b", "1kb", "-f", FRAG, "-c", CHROM, "-F", mat, out_prefix],
    )
    assert result.exit_code == 0, result.output
    # Should fail without --force when output already exists
    result = runner.invoke(
        cli,
        ["rebin", "-b", "1kb", "-f", FRAG, "-c", CHROM, mat, out_prefix],
    )
    assert result.exit_code != 0


def test_convert(runner):
    out_prefix = str(Path(OUT) / "converted")
    result = runner.invoke(
        cli,
        ["convert", "-f", FRAG, "-c", CHROM, "-F", GRAAL, out_prefix],
    )
    assert result.exit_code == 0, result.output
    # Should fail without --force when output already exists
    result = runner.invoke(
        cli,
        ["convert", "-f", FRAG, "-c", CHROM, GRAAL, out_prefix],
    )
    assert result.exit_code != 0


def test_distancelaw(runner):
    result = runner.invoke(
        cli,
        [
            "distancelaw",
            "-a",
            "-o",
            "test.png",
            "-d",
            "test_data/distance_law.txt",
            "-c",
            "test_data/centromeres.txt",
            "-b",
            "10000",
            "-r",
            "1000",
            "-B",
            "1.1",
        ],
    )
    assert result.exit_code == 0, result.output


def test_distance_law_2(runner):
    result = runner.invoke(
        cli,
        [
            "distancelaw",
            "-p",
            "test_data/valid_idx_filtered.pairs",
            "-f",
            FRAG,
            "-C",
            "-O",
            f"{OUT}/test_distance_law.txt",
            "-i",
            "500",
            "-s",
            "45000",
        ],
    )
    assert result.exit_code == 0, result.output


def test_iteralign(runner):
    result = runner.invoke(
        cli,
        [
            "iteralign",
            "-g",
            "test_data/genome/seq",
            "-t",
            "1",
            "-T",
            "tmp",
            "-l",
            "30",
            "-o",
            f"{OUT}/test.bam",
            "test_data/sample.reads_for.fastq.gz",
        ],
    )
    assert result.exit_code == 0, result.output


def test_digest(runner):
    su.rmtree(OUT, ignore_errors=True)
    result = runner.invoke(
        cli,
        ["digest", "-e", "DpnII", "-p", "-f", OUT, "-o", OUT, "test_data/genome/seq.fa"],
    )
    assert result.exit_code == 0, result.output
    # Should fail when output directory already exists without --force
    result = runner.invoke(
        cli,
        ["digest", "-e", "DpnII", "-p", "-f", OUT, "-o", OUT, "test_data/genome/seq.fa"],
    )
    assert result.exit_code != 0
    # Should succeed with --force
    result = runner.invoke(
        cli,
        ["digest", "-e", "DpnII", "-p", "-f", OUT, "-o", OUT, "-F", "test_data/genome/seq.fa"],
    )
    assert result.exit_code == 0, result.output


def test_filter(runner):
    result = runner.invoke(
        cli,
        [
            "filter",
            "-f",
            OUT,
            "-p",
            "test_data/valid_idx.pairs",
            f"{OUT}/valid_idx_filtered.pairs",
        ],
    )
    assert result.exit_code == 0, result.output


def test_scalogram(runner):
    result = runner.invoke(
        cli,
        ["scalogram", "-C", "viridis", "-n", "-t", "1", "-o", f"{OUT}/scalo.png", GRAAL],
    )
    assert result.exit_code == 0, result.output


@pytest.mark.parametrize(*MATS)
def test_subsample(runner, mat):
    out_prefix = str(Path(OUT) / "subsampled")
    result = runner.invoke(
        cli,
        ["subsample", "-p", "0.5", "-F", mat, out_prefix],
    )
    assert result.exit_code == 0, result.output
    # Should fail without --force when output already exists
    result = runner.invoke(
        cli,
        ["subsample", "-p", "0.5", mat, out_prefix],
    )
    assert result.exit_code != 0


@pytest.mark.parametrize("mode", ["for_vs_rev", "all", "pile"])
def test_cutsite(runner, mode):
    result = runner.invoke(
        cli,
        [
            "cutsite",
            "-F",
            "test_data/sample.reads_for.fastq.gz",
            "-R",
            "test_data/sample.reads_rev.fastq.gz",
            "-e",
            "DpnII,HinfI",
            "-m",
            mode,
            "-p",
            "test_data/digested",
            "-t",
            "8",
        ],
    )
    assert result.exit_code == 0, result.output


def test_stats(runner):
    log_file = f"{OUT}/hicstuff.log"
    if not os.path.exists(log_file):
        pytest.skip(f"Log file {log_file} not found; run test_pipeline first.")
    result = runner.invoke(cli, ["stats", log_file])
    assert result.exit_code == 0, result.output
