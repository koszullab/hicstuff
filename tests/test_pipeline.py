# Test functions for the pipeline submodule

import filecmp
import gzip
import logging
import os
import random
import re
import shutil
import subprocess as sp
import sys
from tempfile import NamedTemporaryFile

import click
import cooler
import h5py
import numpy as np
import pandas as pd
import pysam
import pytest

import hicstuff.io as hio
import hicstuff.pipeline as hpi

MAPPING_PARAMETERS = ("mapping", ["normal", "iterative"])
ALIGNER_PARAMETERS = ("aligner", ["bowtie2", "bwa", "minimap2"])
PAIRS_FILE = "test_data/valid_idx_filtered.pairs"


def _read_pairs_coords(pairs_file):
    return pd.read_csv(
        pairs_file,
        sep="\t",
        comment="#",
        header=None,
        usecols=[1, 2, 3, 4],
        names=["chrom1", "pos1", "chrom2", "pos2"],
    )


@pytest.mark.parametrize("exclude", [None, "seq2"])
def test_pairs2binnedcool(tmp_path, exclude):
    """Binned cool pixels must match a direct count of the pairs file."""
    info_contigs = str(tmp_path / "info_contigs.txt")
    shutil.copy("test_data/info_contigs.txt", info_contigs)
    cool_file = str(tmp_path / "binned.cool")
    binning = 1000

    hpi.pairs2binnedcool(PAIRS_FILE, cool_file, binning, info_contigs, exclude)

    clr = cooler.Cooler(cool_file)
    pairs = _read_pairs_coords(PAIRS_FILE)
    if exclude is not None:
        assert exclude not in clr.chromnames
        pairs = pairs[(pairs.chrom1 != exclude) & (pairs.chrom2 != exclude)]
    # The test data must contain inter-chromosomal contacts, so that the
    # chrom2 / pos2 columns are actually exercised.
    if exclude is None:
        assert (pairs.chrom1 != pairs.chrom2).any()

    bins = clr.bins()[:]
    offsets = {chrom: bins.index[bins.chrom == chrom][0] for chrom in clr.chromnames}
    bin1 = pairs.chrom1.map(offsets) + (pairs.pos1 - 1) // binning
    bin2 = pairs.chrom2.map(offsets) + (pairs.pos2 - 1) // binning
    expected = (
        pd.DataFrame({"bin1_id": np.minimum(bin1, bin2), "bin2_id": np.maximum(bin1, bin2)})
        .value_counts()
        .sort_index()
    )
    observed = clr.pixels()[:].set_index(["bin1_id", "bin2_id"])["count"].sort_index()
    np.testing.assert_array_equal(observed.index.to_numpy(), expected.index.to_numpy())
    np.testing.assert_array_equal(observed.to_numpy(), expected.to_numpy())
    # Temporary chroms file must be cleaned up
    assert not os.path.exists(info_contigs + ".chroms")


def test_pairs2cool(tmp_path):
    """Fragment-level cool must contain every pair, cis and trans."""
    frags = str(tmp_path / "fragments_list.txt")
    shutil.copy("test_data/fragments_list.txt", frags)
    cool_file = str(tmp_path / "frags.cool")

    hpi.pairs2cool(PAIRS_FILE, cool_file, frags, None)

    pixels = cooler.Cooler(cool_file).pixels(join=True)[:]
    pairs = _read_pairs_coords(PAIRS_FILE)
    assert pixels["count"].sum() == len(pairs)
    n_trans = (pairs.chrom1 != pairs.chrom2).sum()
    assert pixels.loc[pixels.chrom1 != pixels.chrom2, "count"].sum() == n_trans
    assert not os.path.exists(frags + ".cooler")


def _data_lines(lines):
    return [line for line in lines if not line.startswith("#")]


def _check_sorted_pairs_gz(pairs_gz, expected_data_lines):
    """Check a final pairs.gz: header, content and chr1-chr2-pos1-pos2 order."""
    with gzip.open(pairs_gz, "rt") as f:
        lines = f.read().splitlines()
    assert "#sorted: chr1-chr2-pos1-pos2" in lines
    assert any(line.startswith("#columns: readID chr1 pos1 chr2 pos2") for line in lines)
    data = _data_lines(lines)
    assert sorted(data) == sorted(expected_data_lines)
    keys = [(f[1], f[3], int(f[2]), int(f[4])) for f in (line.split("\t") for line in data)]
    assert keys == sorted(keys)


@pytest.mark.parametrize("compressor", ["bgzip", "pysam"])
def test_sort_and_compress_pairs(tmp_path, hide_tools, compressor):
    if compressor == "bgzip" and shutil.which("bgzip") is None:
        pytest.skip("bgzip is not installed")
    if compressor == "pysam":
        # Without bgzip, pysam compresses the pairs
        hide_tools("bgzip")

    # Shuffle pairs so that the test actually exercises sorting
    with open(PAIRS_FILE) as f:
        lines = f.read().splitlines()
    header, data = hio.get_pairs_header(PAIRS_FILE), _data_lines(lines)
    random.Random(42).shuffle(data)
    unsorted = tmp_path / "unsorted.pairs"
    unsorted.write_text("\n".join(header + data) + "\n")
    out_gz = str(tmp_path / "sorted.pairs.gz")

    hpi.sort_and_compress_pairs(str(unsorted), out_gz, tmp_dir=str(tmp_path), threads=2)

    _check_sorted_pairs_gz(out_gz, data)
    assert not os.path.exists(out_gz + ".tmp")
    # BGZF output starts with a gzip member carrying the BGZF "BC" extra field
    with open(out_gz, "rb") as f:
        assert f.read(16)[12:14] == b"BC"


def test_sort_pairs_failure_raises(tmp_path):
    """A failing UNIX sort must raise instead of leaving a truncated file."""
    with pytest.raises(RuntimeError, match="Sorting pairs"):
        hio.sort_pairs(
            PAIRS_FILE,
            str(tmp_path / "out.pairs"),
            keys=["chr1", "pos1"],
            buffer="not_a_size",
        )


def test_full_pipeline_writes_final_pairs(tmp_path):
    """The pipeline must leave a sorted, compressed pairs file in out_dir."""
    out_dir = tmp_path / "out"
    hpi.full_pipeline(
        input1="test_data/sample.reads_for.fastq.gz",
        input2="test_data/sample.reads_rev.fastq.gz",
        genome="test_data/genome/seq",
        enzyme="DpnII",
        out_dir=str(out_dir),
        prefix="test",
        pcr_duplicates=True,
        no_cleanup=True,
        force=True,
    )
    pairs_gz = out_dir / "test.valid_idx_pcrfree.pairs.gz"
    assert pairs_gz.exists()
    with gzip.open(pairs_gz, "rt") as f:
        data = _data_lines(f.read().splitlines())
    assert len(data) > 0
    _check_sorted_pairs_gz(pairs_gz, data)
    # The uncompressed copy is removed once compression succeeded
    assert not (out_dir / "test.valid_idx_pcrfree.pairs").exists()


def test_cooler_cload_columns_match_pairs_header(tmp_path, monkeypatch):
    """Column numbers given to cooler must match the pairs `#columns:` header.

    Older cooler versions silently reordered swapped columns, so this check
    does not depend on the installed cooler version.
    """
    calls = []
    monkeypatch.setattr(hpi, "_run_cooler", lambda command, args, what: calls.append(args))
    info_contigs = str(tmp_path / "info_contigs.txt")
    shutil.copy("test_data/info_contigs.txt", info_contigs)
    hpi.pairs2binnedcool(PAIRS_FILE, str(tmp_path / "x.cool"), 1000, info_contigs, None)
    cmd = calls[0]

    with open(PAIRS_FILE) as pairs:
        header = next(line for line in pairs if line.startswith("#columns:"))
    columns = header.split()[1:]
    expected = {
        "--chrom1": columns.index("chr1") + 1,
        "--pos1": columns.index("pos1") + 1,
        "--chrom2": columns.index("chr2") + 1,
        "--pos2": columns.index("pos2") + 1,
    }
    observed = {opt: int(cmd[cmd.index(opt) + 1]) for opt in expected}
    assert observed == expected


def test_pairs2binnedcool_cooler_failure(tmp_path):
    """A failing `cooler cload pairs` must raise instead of passing silently."""
    info_contigs = str(tmp_path / "info_contigs.txt")
    shutil.copy("test_data/info_contigs.txt", info_contigs)
    bad_pairs = tmp_path / "bad.pairs"
    bad_pairs.write_text("read1\tseq1\tnot_a_position\tseq1\t100\t+\t+\n")
    cool_file = str(tmp_path / "bad.cool")

    with pytest.raises(RuntimeError, match="cooler cload pairs"):
        hpi.pairs2binnedcool(str(bad_pairs), cool_file, 1000, info_contigs, None)
    assert not os.path.exists(info_contigs + ".chroms")


def _binned_mcool(tmp_path):
    info_contigs = str(tmp_path / "info_contigs.txt")
    shutil.copy("test_data/info_contigs.txt", info_contigs)
    cool, mcool = str(tmp_path / "x.cool"), str(tmp_path / "x.mcool")
    hpi.pairs2binnedcool("test_data/valid_idx.pairs", cool, 1000, info_contigs, None)
    hpi.cool2mcool(cool, mcool)
    return cool, mcool


def test_cooler_steps_need_no_cooler_program(tmp_path, hide_tools):
    """cload and balance run in-process: no `cooler` program on PATH needed."""
    hide_tools("cooler")
    _, mcool = _binned_mcool(tmp_path)
    hpi.balance(mcool, None)
    for path in cooler.fileops.list_coolers(mcool):
        assert "weight" in cooler.Cooler(f"{mcool}::{path}").bins().columns


@pytest.mark.parametrize(
    "balancing_args", [None, "--max-iters 100  --ignore-diags 1", "--cis-only"]
)
def test_balance_matches_cooler_balance_program(tmp_path, balancing_args):
    """In-process balancing writes what `cooler balance` writes."""
    cool, _ = _binned_mcool(tmp_path)
    expected = str(tmp_path / "expected.cool")
    shutil.copy(cool, expected)
    sp.run(
        [sys.executable, "-m", "cooler", "balance", *hpi._balancing_args(balancing_args), expected],
        check=True,
        stderr=sp.DEVNULL,
    )
    hpi.balance(cool, balancing_args)
    with h5py.File(expected, "r") as f_exp, h5py.File(cool, "r") as f_obs:
        exp, obs = f_exp["bins/weight"], f_obs["bins/weight"]
        assert exp[:].tobytes() == obs[:].tobytes()
        assert dict(exp.attrs).keys() == dict(obs.attrs).keys()
        for key in exp.attrs:
            np.testing.assert_array_equal(exp.attrs[key], obs.attrs[key])
        assert (exp.compression, exp.compression_opts, exp.shuffle) == (
            obs.compression,
            obs.compression_opts,
            obs.shuffle,
        )


@pytest.mark.parametrize(
    "balancing_args, message",
    [
        ("--bogus 3", "--bogus"),
        ("--max-iters x", "max-iters"),
        ("--stdout --nproc", "nproc"),
        ("--name 'w", "closing quotation"),
        ("--help", "help"),
    ],
)
def test_invalid_balancing_args_fail_before_the_pipeline_runs(tmp_path, balancing_args, message):
    out_dir = tmp_path / "out"
    with pytest.raises(ValueError, match=message):
        hpi.full_pipeline(
            input1="test_data/sample.reads_for.fastq.gz",
            input2="test_data/sample.reads_rev.fastq.gz",
            genome="test_data/genome/seq",
            enzyme="DpnII",
            out_dir=str(out_dir),
            binning=1000,
            balancing_args=balancing_args,
        )
    assert not out_dir.exists()


def test_balancing_processes_default_to_threads():
    """No process pool in hicstuff's process unless asked for."""
    assert hpi._balancing_args(None) == ["--nproc", "1"]
    assert hpi._balancing_args("--cis-only", threads=4) == ["--nproc", "4", "--cis-only"]
    assert hpi._balancing_args("-p 2 --cis-only", threads=4) == ["-p", "2", "--cis-only"]


def _warnings(caplog):
    return [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]


def test_balance_failure_warns_and_continues(tmp_path, caplog):
    """A failing `cooler balance` logs cooler's error as a warning and leaves
    the matrix unbalanced, instead of stopping the run or passing silently."""
    caplog.set_level(logging.INFO)
    cool, _ = _binned_mcool(tmp_path)
    args = "--max-iters 1 --tol 1e-30 --min-nnz 1 --mad-max 0 --convergence-policy error"
    assert hpi.balance(cool, args) == []
    assert any(
        "Continuing without balancing" in w and "cooler balance" in w for w in _warnings(caplog)
    )
    assert "weight" not in cooler.Cooler(cool).bins().columns
    assert "Matrix balanced" not in caplog.text
    caplog.clear()
    assert hpi.balance(cool, None) == [cool]
    assert "Matrix balanced with cooler balance." in caplog.text
    assert not _warnings(caplog)
    # An existing weight column: cooler >= 0.11 raises "already exists",
    # older versions print it and exit with status 1
    caplog.clear()
    assert hpi.balance(cool, None) == []
    assert any("Continuing without balancing" in w for w in _warnings(caplog))
    assert "Matrix balanced" not in caplog.text


@pytest.mark.parametrize("scope", ["", " --cis-only"])
def test_balance_warns_when_weights_do_not_converge(tmp_path, caplog, scope):
    """Weights stored without convergence are kept, with a warning; with
    --cis-only, cooler stores one convergence flag per chromosome."""
    caplog.set_level(logging.INFO)
    cool, _ = _binned_mcool(tmp_path)
    args = "--max-iters 1 --tol 1e-30 --min-nnz 1 --mad-max 0" + scope
    assert hpi.balance(cool, args) == [cool]
    assert any("did not converge" in w for w in _warnings(caplog))
    with h5py.File(cool, "r") as f:
        assert not np.all(f["bins/weight"].attrs["converged"])


def test_balance_warns_when_no_weights_are_stored(tmp_path, caplog, capsys):
    """--stdout prints the weights: nothing is stored, and that is said."""
    caplog.set_level(logging.INFO)
    cool, _ = _binned_mcool(tmp_path)
    assert hpi.balance(cool, "--stdout --name w2") == []
    assert any("stored no 'w2' column" in w for w in _warnings(caplog))
    assert "Matrix balanced" not in caplog.text


def test_balance_counts_resolutions_that_got_weights(tmp_path, caplog):
    """When only some resolutions of a .mcool balance, the log says how many."""
    caplog.set_level(logging.INFO)
    _, mcool = _binned_mcool(tmp_path)
    paths = cooler.fileops.list_coolers(mcool)
    assert len(hpi.balance(mcool, None)) == len(paths)
    with h5py.File(mcool, "r+") as f:
        del f[paths[0] + "/bins/weight"]
    caplog.clear()
    assert hpi.balance(mcool, None) == [f"{mcool}::{paths[0]}"]
    assert f"Matrix balanced with cooler balance (1 of {len(paths)} resolutions)." in caplog.text
    assert len(_warnings(caplog)) == len(paths) - 1


def test_failed_balancing_does_not_stop_the_pipeline(tmp_path, caplog):
    """The pipeline finishes, unbalanced, when balancing fails at every
    resolution, and says so."""
    caplog.set_level(logging.INFO)
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("not a bed file\n")
    out_dir = tmp_path / "out"
    hpi.full_pipeline(
        input1="test_data/sample.reads_for.fastq.gz",
        input2="test_data/sample.reads_rev.fastq.gz",
        genome="test_data/genome/seq",
        enzyme="DpnII",
        out_dir=str(out_dir),
        binning=1000,
        balancing_args=f"--blacklist {blacklist}",
    )
    (mcool,) = out_dir.glob("*.mcool")
    paths = cooler.fileops.list_coolers(str(mcool))
    for path in paths:
        assert "weight" not in cooler.Cooler(f"{mcool}::{path}").bins().columns
    failed = [w for w in _warnings(caplog) if w.startswith("Continuing without balancing")]
    assert len(failed) == len(paths)
    assert "Matrix balanced" not in caplog.text
    assert "Contact map generated" in caplog.text
    assert list(out_dir.glob("*.pairs.gz"))


UNGUARDED_SCRIPT = """\
import multiprocessing
import sys

multiprocessing.set_start_method("spawn", force=True)
import hicstuff.pipeline as hpi

print("script body runs", flush=True)
print(hpi.balance(sys.argv[1], None, threads=2), flush=True)
"""


def test_balance_needs_no_main_guard(tmp_path):
    """A script without `if __name__ == "__main__":` can balance with several
    processes when they are spawned (macOS, and forkserver from Python 3.14)."""
    cool, _ = _binned_mcool(tmp_path)
    script = tmp_path / "unguarded.py"
    script.write_text(UNGUARDED_SCRIPT)
    proc = sp.run(
        [sys.executable, str(script), cool],
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.count("script body runs") == 1
    assert proc.stdout.splitlines()[-1] == str([cool])
    assert "weight" in cooler.Cooler(cool).bins().columns


def test_balance_in_a_new_process_matches_in_process(tmp_path, monkeypatch):
    """Without fork, several balancing processes run in `python -m cooler`,
    which writes the same weights as balancing in hicstuff's process."""
    cool, _ = _binned_mcool(tmp_path)
    expected = str(tmp_path / "expected.cool")
    shutil.copy(cool, expected)
    hpi.balance(expected, None, threads=1)
    calls = []
    run_program = hpi._run_cooler_program
    monkeypatch.setattr(hpi.multiprocessing, "get_start_method", lambda: "spawn")
    monkeypatch.setattr(hpi, "_run_cooler_program", lambda *a: calls.append(a) or run_program(*a))
    assert hpi.balance(cool, None, threads=2) == [cool]
    assert len(calls) == 1
    with h5py.File(expected, "r") as f_exp, h5py.File(cool, "r") as f_obs:
        assert f_exp["bins/weight"][:].tobytes() == f_obs["bins/weight"][:].tobytes()
    # One process, or fork: cooler runs in hicstuff's process
    hpi.balance(cool, "--force", threads=1)
    monkeypatch.setattr(hpi.multiprocessing, "get_start_method", lambda: "fork")
    hpi.balance(cool, "--force", threads=2)
    assert len(calls) == 1


def test_run_cooler_program_failure(tmp_path):
    """A failing `python -m cooler` raises with its exit status and the last
    line of its log."""
    cool, _ = _binned_mcool(tmp_path)
    with pytest.raises(RuntimeError, match=r"balancing failed \(exit status 2\): .*--bogus"):
        hpi._run_cooler_program(["balance"], ["--bogus", cool], "balancing")


def test_run_cooler_exit_status_and_interrupt():
    @click.command()
    @click.argument("code", type=int)
    def exits(code):
        sys.exit(code)

    @click.command()
    @click.argument("code", type=int)
    @click.pass_context
    def ctx_exits(ctx, code):
        ctx.exit(code)

    @click.command()
    def interrupted():
        raise KeyboardInterrupt

    @click.command()
    def truncated():
        raise EOFError("Compressed file ended before the end-of-stream marker")

    hpi._run_cooler(exits, ["0"], "exits")
    with pytest.raises(RuntimeError, match="exits failed .exit status 3"):
        hpi._run_cooler(exits, ["3"], "exits")
    # click returns, rather than raises, the status of ctx.exit()
    hpi._run_cooler(ctx_exits, ["0"], "ctx_exits")
    with pytest.raises(RuntimeError, match="ctx_exits failed .exit status 2"):
        hpi._run_cooler(ctx_exits, ["2"], "ctx_exits")
    with pytest.raises(KeyboardInterrupt):
        hpi._run_cooler(interrupted, [], "interrupted")
    # click turns an EOFError into Abort, which says nothing
    with pytest.raises(RuntimeError, match="truncated failed: EOFError"):
        hpi._run_cooler(truncated, [], "truncated")


def test_cooler_commands_are_found_by_name():
    """The commands are those of the `cooler` program, whatever their module."""
    import cooler.cli

    assert hpi._cooler_command("balance") is cooler.cli.cli.commands["balance"]
    assert hpi._cooler_command("cload", "pairs").name == "pairs"


def test_sam2pairs(): ...


def test_pairs2mat(): ...


def test_filter_pcr_dup():
    """Test if PCR duplicates are removed correctly"""
    dup_rm = NamedTemporaryFile(mode="w", delete=False)
    dup_rm.close()
    test_pairs = NamedTemporaryFile(mode="w", delete=False)
    test_rm = NamedTemporaryFile(mode="w", delete=False)
    lnum = 0
    # Copy the test valid_idx file, but generate PCR dups of the pair at line 50
    with open("test_data/valid_idx.pairs") as pairs:
        for line in pairs:
            test_pairs.write(line)
            if lnum == 50:
                # Making 30 duplicates of this pair
                for _i in range(30):
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


def test_filter_pcr_dup_after_coordinate_sort(tmp_path):
    """A pair at the same positions on other strands does not hide a duplicate."""
    header = (
        "## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 frag1 frag2\n"
    )
    lines = [
        "rA\tseq1\t10\tseq1\t50\t+\t+\t1\t2",
        "rB\tseq1\t10\tseq1\t50\t+\t-\t1\t2",
        "rC\tseq1\t10\tseq1\t50\t+\t+\t1\t2",
    ]
    unsorted, sorted_pairs = tmp_path / "unsorted.pairs", tmp_path / "sorted.pairs"
    unsorted.write_text(header + "\n".join(lines) + "\n")
    hio.sort_pairs(str(unsorted), str(sorted_pairs), keys=hpi.COORD_SORT_KEYS)
    pcr_free = tmp_path / "pcrfree.pairs"
    hpi.filter_pcr_dup(str(sorted_pairs), str(pcr_free))
    kept = [line.split("\t")[0] for line in pcr_free.read_text().splitlines() if line[0] != "#"]
    # rC duplicates rA; sorted by coordinates only, rB came between them
    assert kept == ["rA", "rB"]


def test_full_pipeline_frags():
    mapping_to_enzyme = {"normal": "DpnII", "cutsite": "DpnII,HpaII"}
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


# Testing:
# no binning (no zoomify, no balancing)
# binning: no zoomify (balancing)
# binning: no zoomify, tweaked balancing
# binning: zoomify, (balancing)
# binning: zoomify, tweaked balancing


def test_full_pipeline_frags_with_binning_balancing():
    hpi.full_pipeline(
        input1="test_data/sample.reads_for.fastq.gz",
        input2="test_data/sample.reads_rev.fastq.gz",
        genome="test_data/genome/seq",
        enzyme="DpnII,HpaII",
        mapping="normal",
        out_dir="test_out",
        plot=True,
        pcr_duplicates=True,
        filter_events=True,
        no_cleanup=True,
        force=True,
    )
    hpi.full_pipeline(
        input1="test_data/sample.reads_for.fastq.gz",
        input2="test_data/sample.reads_rev.fastq.gz",
        genome="test_data/genome/seq",
        enzyme="DpnII,HpaII",
        mapping="normal",
        out_dir="test_out",
        plot=True,
        pcr_duplicates=True,
        filter_events=True,
        no_cleanup=True,
        force=True,
        binning=1000,
        zoomify=False,
    )
    hpi.full_pipeline(
        input1="test_data/sample.reads_for.fastq.gz",
        input2="test_data/sample.reads_rev.fastq.gz",
        genome="test_data/genome/seq",
        enzyme="DpnII,HpaII",
        mapping="normal",
        out_dir="test_out",
        plot=True,
        pcr_duplicates=True,
        filter_events=True,
        no_cleanup=True,
        force=True,
        binning=1000,
        zoomify=False,
        balancing_args="--min-nnz 1 --mad-max 1",
    )
    hpi.full_pipeline(
        input1="test_data/sample.reads_for.fastq.gz",
        input2="test_data/sample.reads_rev.fastq.gz",
        genome="test_data/genome/seq",
        enzyme="DpnII,HpaII",
        mapping="normal",
        out_dir="test_out",
        plot=True,
        pcr_duplicates=True,
        filter_events=True,
        no_cleanup=True,
        force=True,
        binning=1000,
    )
    hpi.full_pipeline(
        input1="test_data/sample.reads_for.fastq.gz",
        input2="test_data/sample.reads_rev.fastq.gz",
        genome="test_data/genome/seq",
        enzyme="DpnII,HpaII",
        mapping="normal",
        out_dir="test_out",
        plot=True,
        pcr_duplicates=True,
        filter_events=True,
        no_cleanup=True,
        force=True,
        binning=1000,
        balancing_args="--min-nnz 1 --mad-max 1",
    )
    hpi.full_pipeline(
        input1="test_data/sample.reads_for.fastq.gz",
        input2="test_data/sample.reads_rev.fastq.gz",
        genome="test_data/genome/seq",
        enzyme="DpnII,HpaII",
        mapping="normal",
        out_dir="test_out",
        plot=True,
        pcr_duplicates=True,
        filter_events=True,
        no_cleanup=False,
        force=True,
        binning=1000,
        exclude="seq2",
        balancing_args="--min-nnz 1 --mad-max 1",
    )
    hpi.full_pipeline(
        input1="test_data/sample.reads_for.fastq.gz",
        input2="test_data/sample.reads_rev.fastq.gz",
        genome="test_data/genome/seq_chromnames-as-number",
        enzyme="DpnII,HpaII",
        mapping="normal",
        out_dir="test_out",
        plot=True,
        pcr_duplicates=True,
        filter_events=True,
        no_cleanup=False,
        force=True,
        binning=1000,
        distance_law=True,
        balancing_args="--min-nnz 1 --mad-max 1",
    )
    hpi.full_pipeline(
        input1="test_data/sample.reads_for.fastq.gz",
        input2="test_data/sample.reads_rev.fastq.gz",
        genome="test_data/genome/seq_chromnames-as-number",
        enzyme="DpnII,HpaII",
        mapping="normal",
        out_dir="test_out",
        plot=True,
        pcr_duplicates=True,
        filter_events=True,
        no_cleanup=False,
        force=True,
        binning=1000,
        distance_law=True,
        balancing_args="--min-nnz 1 --mad-max 1",
    )


def test_full_pipeline_frags_gzipped_genome():
    mapping_to_enzyme = {"normal": "DpnII", "cutsite": "DpnII,HpaII"}
    for mapping, enzyme in mapping_to_enzyme.items():
        hpi.full_pipeline(
            input1="test_data/sample.reads_for.fastq.gz",
            input2="test_data/sample.reads_rev.fastq.gz",
            genome="test_data/genome/seq.fa.gz",
            enzyme=enzyme,
            mapping=mapping,
            mat_fmt="cool",
            binning=1000,
            zoomify=True,
            balancing_args="--max-iters 10 --min-nnz 10",
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
    out_dir = f"test_out_{aligner}_{mapping}"
    # Copy test_data/valid_idx.pairs into a writable location so the pipeline
    # can sort it in-place without touching the read-only test_data directory.
    os.makedirs(os.path.join(out_dir, "tmp"), exist_ok=True)
    pairs_idx_copy = os.path.join(out_dir, "tmp", "valid_idx_input.pairs")
    shutil.copy("test_data/valid_idx.pairs", pairs_idx_copy)
    start_input = {
        "fastq": [
            "test_data/sample.reads_for.fastq.gz",
            "test_data/sample.reads_rev.fastq.gz",
        ],
        "bam": ["test_out/tmp/for.bam", "test_out/tmp/rev.bam"],
        "pairs": ["test_out/tmp/valid.pairs", None],
        "pairs_idx": [pairs_idx_copy, None],
    }
    for stage, [in1, in2] in start_input.items():
        # Indexed or non-indexed genome
        # for genome in ['test_data/genome/seq', 'test_data/genome/seq.fa']:
        hpi.full_pipeline(
            input1=in1,
            input2=in2,
            genome="test_data/genome/seq.fa.gz",
            enzyme=5000,
            out_dir=out_dir,
            aligner=aligner,
            mapping=mapping,
            prefix="test",
            distance_law=True,
            start_stage=stage,
            mat_fmt="cool",
            no_cleanup=True,
            force=True,
        )
    shutil.rmtree(out_dir)
    if mapping == "iterative" and aligner == "minimap2":
        shutil.rmtree("test_out/")


# ---------------------------------------------------------------------------
# bam2pairs: matching mates in name-sorted BAM files
# ---------------------------------------------------------------------------

BAM_REFS = [("chr1", 20000), ("chr2", 8000)]


def _samtools_version():
    """Version of the samtools on the PATH, as a tuple of ints; skip without one."""
    if shutil.which("samtools") is None:
        pytest.skip("samtools is not installed")
    # Only the first line: the rest is not always valid UTF-8
    out = sp.run(["samtools", "--version"], capture_output=True, check=True).stdout
    return tuple(int(x) for x in re.match(rb"samtools (\d+)\.(\d+)", out).groups())


def _write_bam(path, reads, header_hd=None):
    """Write reads, given as (name, ref_id, pos, reverse, mapq) in that order,
    to a BAM file. ref_id -1 makes an unmapped read."""
    header = {"HD": header_hd or {"VN": "1.6"}, "SQ": [{"SN": n, "LN": ln} for n, ln in BAM_REFS]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for name, ref_id, pos, reverse, mapq in reads:
            read = pysam.AlignedSegment(bam.header)
            read.query_name = name
            read.query_sequence = "A" * 20
            if ref_id < 0:
                read.flag = 4
            else:
                read.flag = 16 if reverse else 0
                read.reference_id = ref_id
                read.reference_start = pos
                read.mapping_quality = mapq
                read.cigarstring = "20M"
            bam.write(read)


def _name_sorted_bam(path, reads, flag, tmp_path):
    """Name-sorted BAM of reads, sorted by ``samtools sort <flag>``."""
    unsorted = tmp_path / f"unsorted_{path.name}"
    _write_bam(unsorted, reads)
    sp.run(["samtools", "sort", flag, "-o", str(path), str(unsorted)], check=True)
    return str(path)


def _bam_names(path):
    with pysam.AlignmentFile(str(path), "rb", check_sq=False) as bam:
        return [read.query_name for read in bam]


def test_natural_key_follows_samtools_sort_n(tmp_path):
    """Read names come out of ``samtools sort -n`` in _natural_key order.

    The names mix digit runs, with and without leading zeros, with characters
    sorting before ('!', '-', '.', '/') and after (':', 'A', '_', 'z', '~')
    digits."""
    version = _samtools_version()
    rng = random.Random(0)
    alphabet = "0000123456789" * 2 + "!-./:A_z~"
    names = {"r9", "r10", "r9.1", "r9:1", "r010", "x1", "x.1", "x:1", "a0b", "a00c"}
    while len(names) < 3000:
        names.add("".join(rng.choice(alphabet) for _ in range(rng.randint(1, 10))))
    names = sorted(names)
    rng.shuffle(names)
    sorted_names = _bam_names(
        _name_sorted_bam(
            tmp_path / "sorted.bam", [(n, -1, 0, False, 0) for n in names], "-n", tmp_path
        )
    )
    if version < (1, 17):
        # samtools 1.16 and earlier put r010 before r10, before comparing what
        # follows the numbers: leave out names with zero-padded numbers.
        padded = re.compile(r"(?<![0-9])0[0-9]")
        sorted_names = [name for name in sorted_names if not padded.search(name)]
    keys = [hpi._natural_key(name) for name in sorted_names]
    # Names equal under the key (r010, r10) may come in either order
    assert all(a <= b for a, b in zip(keys, keys[1:]))
    # The test is only meaningful if samtools' order is not the string order
    assert sorted_names != sorted(sorted_names)


def _alignments(rng, names):
    """Random alignment (contig, position, reverse strand, MAPQ) of each read,
    a tenth of them with MAPQ < 30."""
    return {
        name: (
            rng.randrange(2),
            rng.randrange(7000),
            rng.random() < 0.5,
            rng.choice([5] + [40] * 9),
        )
        for name in names
    }


def _run_bam2pairs(tmp_path, mates, **kwargs):
    """Run bam2pairs on tmp_path/for.bam and tmp_path/rev.bam, and check that
    it writes the pairs of the reads of ``mates`` (the alignments of each end,
    by read name) with MAPQ >= 30, flipped to the upper triangle."""
    info_contigs = tmp_path / "info_contigs.txt"
    info_contigs.write_text(
        "contig\tlength\tn_frags\tcumul_length\n"
        + "".join(f"{name}\t{length}\t1\t0\n" for name, length in BAM_REFS)
    )
    out = tmp_path / "out.pairs"
    hpi.bam2pairs(
        str(tmp_path / "for.bam"), str(tmp_path / "rev.bam"), str(out), str(info_contigs), **kwargs
    )
    expected = []
    for name in mates["for"].keys() & mates["rev"].keys():
        aln1, aln2 = mates["for"][name], mates["rev"][name]
        if min(aln1[3], aln2[3]) < 30:
            continue
        if aln1[:2] > aln2[:2]:
            aln1, aln2 = aln2, aln1
        expected.append(
            [name]
            + [BAM_REFS[aln1[0]][0], aln1[1] + 1, BAM_REFS[aln2[0]][0], aln2[1] + 1]
            + ["-" if aln1[2] else "+", "-" if aln2[2] else "+"]
        )
    pairs = pd.read_csv(out, sep="\t", comment="#", header=None, dtype={0: str})
    assert sorted(pairs.values.tolist()) == sorted(expected)


@pytest.mark.parametrize("suffixes", [("", ""), ("/1", "/2")], ids=["plain", "suffixed"])
@pytest.mark.parametrize("flag", ["-n", "-N"])
def test_bam2pairs_matches_mates_in_samtools_name_order(tmp_path, caplog, flag, suffixes):
    """Every pair is found when mates are missing from either BAM file, as with
    BAM files aligned with --no-unal, in the order samtools sorted them in.

    Read names r1, r2... sort differently as strings (r10 < r9) and naturally
    (r9 < r10, samtools sort -n)."""
    if _samtools_version() < (1, 19) and flag == "-N":
        pytest.skip("samtools sort -N needs samtools 1.19 or later")
    rng = random.Random(1)
    n_reads = 300
    names = [f"r{i}" for i in range(1, n_reads + 1)]
    mates = {end: _alignments(rng, names) for end in ("for", "rev")}
    # Reads absent from one file, as when unaligned reads are not written
    del_for = {name for i, name in enumerate(names, 1) if i % 7 == 3 or i % 11 == 0}
    del_rev = {name for i, name in enumerate(names, 1) if i % 5 == 2 or i in (9, 99)}
    for end, deleted, suffix in zip(("for", "rev"), (del_for, del_rev), suffixes):
        mates[end] = {name: mates[end][name] for name in names if name not in deleted}
        reads = [(name + suffix, *aln) for name, aln in mates[end].items()]
        rng.shuffle(reads)
        _name_sorted_bam(tmp_path / f"{end}.bam", reads, flag, tmp_path)

    _run_bam2pairs(tmp_path, mates)

    n_unmatched = len(del_for ^ del_rev)
    assert f"{n_unmatched} reads were only present in one BAM file" in caplog.text
    # Both files were walked once, in the order their header declares
    assert "is not sorted by read name" not in caplog.text


@pytest.mark.parametrize("suffixes", [("", ""), ("/1", "/2")], ids=["plain", "suffixed"])
def test_bam2pairs_matches_mates_whose_names_tie(tmp_path, caplog, suffixes):
    """Names whose numbers are zero-padded differently (r1, r001) are equal to
    samtools sort -n, which (from samtools 1.17) leaves them in input order:
    their mates are matched all the same. Names extended by a character sorting before '/'
    (r1, r1.5) sort before them once suffixed (r1.5/1 < r1/1), and are walked
    in that order."""
    rng = random.Random(3)
    names = sorted(
        {name for i in range(1, 121) for name in (f"r{i}", f"r{i:03d}", f"r{i}.5", f"r{i:02d}x")}
    )
    rng.shuffle(names)
    mates = {}
    for end, suffix in zip(("for", "rev"), suffixes):
        # Reads absent from one file, as when unaligned reads are not written
        kept = [name for name in names if rng.random() > 0.25]
        mates[end] = _alignments(rng, kept)
        # Both files list reads in the same input (FASTQ) order before sorting
        _write_bam(tmp_path / f"unsorted_{end}.bam", [(n + suffix, *mates[end][n]) for n in kept])
        pysam.sort("-n", "-o", str(tmp_path / f"{end}.bam"), str(tmp_path / f"unsorted_{end}.bam"))

    _run_bam2pairs(tmp_path, mates)

    n_unmatched = len(mates["for"].keys() ^ mates["rev"].keys())
    assert f"{n_unmatched} reads were only present in one BAM file" in caplog.text
    assert "sorting both BAM files" not in caplog.text


def _bams_in_order(tmp_path, order, all_mates=False):
    """Write tmp_path/for.bam and tmp_path/rev.bam, with reads r1 to r300 in
    the given order, and return the alignments of each end by read name."""
    rng = random.Random(2)
    names = [f"r{i}" for i in range(1, 301)]
    rng.shuffle(names)
    lexicographical = {"VN": "1.6", "SO": "queryname", "SS": "queryname:lexicographical"}
    mates = {}
    for end in ("for", "rev"):
        kept = [name for name in names if all_mates or rng.random() > 0.25]
        mates[end] = _alignments(rng, kept)
        reads = [(name, *mates[end][name]) for name in kept]
        path = tmp_path / f"{end}.bam"
        if order == "lexicographical":
            # As sorted by Picard: the header does not say in which order
            _write_bam(path, sorted(reads), header_hd={"VN": "1.6", "SO": "queryname"})
        elif order == "mixed" and end == "for":
            _write_bam(tmp_path / "unsorted_for.bam", reads)
            pysam.sort("-n", "-o", str(path), str(tmp_path / "unsorted_for.bam"))
        elif order == "mixed":
            _write_bam(path, sorted(reads), header_hd=lexicographical)
        elif order == "coordinate":
            reads.sort(key=lambda read: read[1:3])
            _write_bam(path, reads, header_hd={"VN": "1.6", "SO": "coordinate"})
        else:
            # In the input (FASTQ) order, as written by the aligner
            _write_bam(path, reads, header_hd={"VN": "1.6", "SO": "unsorted"})
    return mates


@pytest.mark.parametrize(
    "order, all_mates, walk",
    [
        ("lexicographical", False, "lexicographical"),
        ("mixed", False, "sorted"),
        ("coordinate", False, "sorted"),
        ("unsorted", False, "sorted"),
        ("unsorted", True, "natural"),
    ],
    ids=["lexicographical", "mixed", "coordinate", "unsorted", "unsorted-all-mates"],
)
def test_bam2pairs_finds_mates_in_bams_in_another_order(tmp_path, caplog, order, all_mates, walk):
    """Every pair is found in BAM files not in the order of samtools sort -n.

    Files in lexicographical order whose header does not say so (as sorted by
    Picard) are walked in that order. Files sorted by name in different
    orders, sorted by coordinate or not sorted are sorted again with samtools
    sort -n, when reads have no mate in the other file: files whose mates all
    follow each other can be walked as they are."""
    mates = _bams_in_order(tmp_path, order, all_mates=all_mates)
    sort_dir = tmp_path / "sort"
    sort_dir.mkdir()

    _run_bam2pairs(tmp_path, mates, tmp_dir=str(sort_dir), threads=2)

    n_unmatched = len(mates["for"].keys() ^ mates["rev"].keys())
    if all_mates:
        assert "reads were only present in one BAM file" not in caplog.text
    else:
        assert f"{n_unmatched} reads were only present in one BAM file" in caplog.text
    retried = "is not sorted by read name as samtools sort -n sorts it" in caplog.text
    resorted = "sorting both BAM files with samtools sort -n" in caplog.text
    assert (retried, resorted) == {
        "natural": (False, False),
        "lexicographical": (True, False),
        "sorted": (True, True),
    }[walk]
    assert list(sort_dir.iterdir()) == []


def _sort_in_place(*args):
    """A samtools sort that leaves the reads as they are."""
    shutil.copy(args[-1], args[-2])


def _failing_sort(*args):
    raise pysam.utils.SamtoolsError("samtools sort failed")


@pytest.mark.parametrize(
    "sort, error, message",
    [
        (_sort_in_place, ValueError, "still not in the order"),
        (_failing_sort, pysam.utils.SamtoolsError, "failed"),
    ],
    ids=["names-out-of-order", "sort-fails"],
)
def test_bam2pairs_writes_nothing_if_sorting_does_not_help(
    tmp_path, monkeypatch, sort, error, message
):
    """If BAM files sorted again with samtools sort -n are still not in its
    order, bam2pairs raises rather than write pairs that may miss some; if
    sorting them fails, the error is raised. Either way, neither pairs nor
    sorted files are left behind."""
    _bams_in_order(tmp_path, "unsorted")
    monkeypatch.setattr(hpi.ps, "sort", sort)
    with pytest.raises(error, match=message):
        hpi.bam2pairs(
            str(tmp_path / "for.bam"),
            str(tmp_path / "rev.bam"),
            str(tmp_path / "out.pairs"),
            "test_data/info_contigs.txt",
            tmp_dir=str(tmp_path),
        )
    assert sorted(path.name for path in tmp_path.iterdir()) == ["for.bam", "rev.bam"]
