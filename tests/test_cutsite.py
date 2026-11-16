"""Cutting reads at ligation sites: the writer thread and its errors."""

import gzip
import subprocess as sp
import sys
import threading

import pytest

import hicstuff.cutsite as hcc

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
