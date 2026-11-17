"""Cutting reads at ligation sites: reading the reads, the writer thread and
its errors."""

import gzip
import importlib
import subprocess as sp
import sys
import threading

import pysam
import pytest

import hicstuff.cutsite as hcc
import hicstuff.digest as hcd

FOR = "test_data/sample.reads_for.fastq.gz"
REV = "test_data/sample.reads_rev.fastq.gz"

UNGUARDED_SCRIPT = """\
import multiprocessing
import sys

multiprocessing.set_start_method("spawn", force=True)
import hicstuff.cutsite as hcc

print("script body runs", flush=True)
hcc.cut_ligation_sites(*sys.argv[1:5], "DpnII", "for_vs_rev", 20, 4)
"""


def test_cutsite_needs_no_main_guard(tmp_path):
    """A script without `if __name__ == "__main__":` can cut reads where
    processes are spawned (macOS, and forkserver from Python 3.14)."""
    script = tmp_path / "unguarded.py"
    script.write_text(UNGUARDED_SCRIPT)
    out_for, out_rev = tmp_path / "for.fq.gz", tmp_path / "rev.fq.gz"
    proc = sp.run(
        [sys.executable, str(script), *map(str, [FOR, REV, out_for, out_rev])],
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.count("script body runs") == 1
    with gzip.open(out_for, "rt") as f_for, gzip.open(out_rev, "rt") as f_rev:
        assert sum(1 for _ in f_for) == sum(1 for _ in f_rev) > 40000


def test_cutsite_raises_writer_error(tmp_path):
    """A writer failure is raised, rather than leaving the run blocked."""
    missing = tmp_path / "missing" / "dir"
    threads = threading.active_count()
    with pytest.raises(FileNotFoundError):
        hcc.cut_ligation_sites(
            FOR,
            REV,
            str(missing / "for.fq.gz"),
            str(missing / "rev.fq.gz"),
            "DpnII",
            "for_vs_rev",
            20,
            2,
        )
    assert threading.active_count() == threads


def test_cutsite_stops_writer_on_unsorted_reads(tmp_path):
    """Reads that are not in the same order stop the run, and the writer
    thread with it."""
    shuffled = tmp_path / "shuffled.fq.gz"
    with gzip.open(REV, "rt") as f:
        lines = f.read().splitlines(keepends=True)
    records = [lines[i : i + 4] for i in range(0, len(lines), 4)]
    with gzip.open(shuffled, "wt") as f:
        f.writelines(line for record in records[::-1] for line in record)
    threads = threading.active_count()
    with pytest.raises(SystemExit):
        hcc.cut_ligation_sites(
            FOR,
            str(shuffled),
            str(tmp_path / "f.fq.gz"),
            str(tmp_path / "r.fq.gz"),
            "DpnII",
            "for_vs_rev",
            20,
            2,
        )
    assert threading.active_count() == threads


#: Reads with an Illumina comment, a tab, a mate suffix, no bases, and
#: lowercase bases, with the (name, sequence, quality) expected of each.
ODD_READS = (
    "@r1 1:N:0:ATCACG\nACGTAC\n+\nIIIIII\n"
    "@r2\tBC:Z:AAA\nGGCC\n+r2\nJJJJ\n"
    "@r3/1\nAC\n+\n@I\n"
    "@r4\n\n+\n\n"
    "@r5\nacgtn\n+\nIIIII\n"
)
ODD_RECORDS = [
    ("r1", "ACGTAC", "IIIIII"),
    ("r2", "GGCC", "JJJJ"),
    ("r3/1", "AC", "@I"),
    ("r4", "", ""),
    ("r5", "acgtn", "IIIII"),
]


@pytest.mark.parametrize("compression", ["plain", "gzip", "bgzf"])
def test_read_fastq(compression, tmp_path):
    """Reads are read by pysam, and an empty read has an empty quality."""
    path = tmp_path / "reads.fq"
    if compression == "plain":
        path.write_text(ODD_READS)
    elif compression == "gzip":
        path.write_bytes(gzip.compress(ODD_READS.encode()))
    else:
        with pysam.BGZFile(str(path), "wb") as handle:
            handle.write(ODD_READS.encode())
    assert list(hcc._read_fastq(str(path))) == ODD_RECORDS


@pytest.mark.parametrize(
    "text",
    [">r1\nACGT\n>r2\nGG\n", "@r1\nACGT\n+\nIIII\n@r2\nAC\n"],
    ids=["fasta", "truncated"],
)
def test_read_fastq_refuses_reads_without_quality(text, tmp_path):
    path = tmp_path / "reads.fq"
    path.write_text(text)
    with pytest.raises(ValueError, match="is not a fastq file: read r[12] has no quality"):
        list(hcc._read_fastq(str(path)))


def _plain_fastq(path):
    """(name, sequence, quality) of each read of a 4-line fastq file."""
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    return [(lines[i][1:].split()[0], lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4)]


@pytest.mark.parametrize("mode", ["for_vs_rev", "all", "pile"])
def test_cutsite_without_pyfastx(mode, tmp_path, monkeypatch):
    """Reads are read with pysam, not pyfastx, and cut into the pairs that the
    reads of a plain 4-line parser are cut into."""
    monkeypatch.setitem(sys.modules, "pyfastx", None)
    importlib.reload(hcc)
    out_for, out_rev = tmp_path / "for.fq.gz", tmp_path / "rev.fq.gz"
    hcc.cut_ligation_sites(FOR, REV, str(out_for), str(out_rev), "DpnII", mode, 20, 2)

    sites = hcd.gen_enzyme_religation_regex("DpnII")
    reads = list(zip(_plain_fastq(FOR), _plain_fastq(REV)))
    expected_for, expected_rev, n_pairs = "", "", 0
    for (name, seq_for, qual_for), (_, seq_rev, qual_rev) in reads:
        cut_for = hcc.cutsite_read(sites, seq_for, qual_for, 20)
        cut_rev = hcc.cutsite_read(sites, seq_rev, qual_rev, 20)
        expected_for, expected_rev, n_pairs = hcc.write_pair(
            expected_for, expected_rev, name, *cut_for, *cut_rev, mode, n_pairs
        )
    with gzip.open(out_for, "rt") as f_for, gzip.open(out_rev, "rt") as f_rev:
        assert (f_for.read(), f_rev.read()) == (expected_for, expected_rev)
    # Some reads were cut into several pairs
    assert n_pairs > len(reads)
