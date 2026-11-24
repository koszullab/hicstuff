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
import metacooler
import numpy as np
import pandas as pd
import pysam
import pytest

import hicstuff.bwamem as hbm
import hicstuff.io as hio
import hicstuff.pipeline as hpi
import hicstuff.stats as hcs
from tests import synthetic as syn
from tests.conftest import skip_unless_can_align

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


@pytest.mark.needs("sort")
def test_sort_pairs_failure_raises(tmp_path):
    """A failing UNIX sort must raise instead of leaving a truncated file."""
    with pytest.raises(RuntimeError, match="Sorting pairs"):
        hio.sort_pairs(
            PAIRS_FILE,
            str(tmp_path / "out.pairs"),
            keys=["chr1", "pos1"],
            buffer="not_a_size",
        )


#: Files of the bowtie2 and bwa indexes at prefix "genome", and of the index
#: that pybwa builds, with the sequence dictionary it loads, at "index".
BOWTIE2_INDEX = [f"genome.{name}.bt2" for name in ("1", "2", "3", "4", "rev.1", "rev.2")]
BWA_INDEX = [f"genome.{ext}" for ext in ("amb", "ann", "bwt", "pac", "sa")]
PYBWA_INDEX = [f"index.{ext}" for ext in ("amb", "ann", "bwt", "dict", "pac", "sa")]


@pytest.mark.parametrize("genome_name", ["genome.fa", "genome.fa.gz"])
@pytest.mark.parametrize("no_cleanup", [False, True])
@pytest.mark.parametrize(
    ("aligner", "mapping", "backend", "programs", "index_files"),
    [
        ("bowtie2", "normal", None, ("bowtie2", "bowtie2-build"), BOWTIE2_INDEX),
        ("bwa", "normal", "binary", ("bwa",), BWA_INDEX),
        ("bwa", "5spm", "python", (), PYBWA_INDEX),
    ],
)
def test_genome_indexed_in_a_directory_of_the_run(
    aligner,
    mapping,
    backend,
    programs,
    index_files,
    no_cleanup,
    genome_name,
    synthetic_lib,
    tmp_path,
):
    """A fasta file without index, plain or compressed, is indexed where it
    is, once, in a directory of tmp_dir of the run, which cleanup removes. It
    is not copied, and nothing is written next to it."""
    missing = [program for program in programs if shutil.which(program) is None]
    if missing:
        pytest.skip(f"not on PATH: {', '.join(missing)}")
    if backend == "python" and hbm.pybwa_version() is None:
        pytest.skip("pybwa is not installed")
    genome = tmp_path / genome_name
    compressed = genome_name.endswith(".gz")
    shutil.copy(synthetic_lib.genome_gz if compressed else synthetic_lib.genome, genome)
    out_dir = tmp_path / "out"
    hpi.full_pipeline(
        genome=str(genome),
        input1=str(synthetic_lib.reads1),
        input2=str(synthetic_lib.reads2),
        aligner=aligner,
        aligner_backend=backend,
        mapping=mapping,
        enzyme=syn.LIBRARY_ENZYME,
        out_dir=str(out_dir),
        no_cleanup=no_cleanup,
    )
    assert sorted(os.listdir(tmp_path)) == [genome_name, "out"]
    (log,) = out_dir.glob("*.log")
    log_text = log.read_text()
    tmp_dir = out_dir / "tmp"
    if not no_cleanup:
        assert f"{aligner} index not found at {genome}, building a temporary one in" in log_text
        assert list(tmp_dir.iterdir()) == []
        return
    (genome_dir,) = [path for path in tmp_dir.iterdir() if path.is_dir()]
    # The log names where the index is built: pybwa builds it in a
    # directory of genome_dir
    built = f"{aligner} index not found at {genome}, building a temporary one in {genome_dir}"
    assert built in log_text
    # genome_dir holds one index, and no copy of the genome
    files = [path for path in genome_dir.rglob("*") if path.is_file()]
    assert sorted(path.name for path in files) == index_files
    assert len({path.parent for path in files}) == 1
    # The reads were aligned to it
    with pysam.AlignmentFile(str(tmp_dir / "for.bam")) as bam:
        assert any(not read.is_unmapped for read in bam)


@pytest.mark.parametrize("aligner", ["bowtie2", "bwa"])
@pytest.mark.parametrize(
    ("stage", "pairs"), [("pairs", "valid.pairs"), ("pairs_idx", "valid_idx.pairs")]
)
def test_genome_not_indexed_without_alignment(aligner, stage, pairs, tmp_path):
    """A run that aligns no reads neither copies nor indexes a fasta file
    without index: it reads it where it is."""
    genome = tmp_path / "genome.fa"
    shutil.copy("test_data/genome/seq.fa", genome)
    # The pairs file of stage pairs_idx is the run's own: it sorts it in place
    shutil.copy(f"test_data/{pairs}", tmp_path / pairs)
    out_dir = tmp_path / "out"
    hpi.full_pipeline(
        genome=str(genome),
        input1=str(tmp_path / pairs),
        start_stage=stage,
        aligner=aligner,
        enzyme="DpnII",
        out_dir=str(out_dir),
        no_cleanup=True,
    )
    assert sorted(os.listdir(tmp_path)) == ["genome.fa", "out", pairs]
    # No directory for the genome in tmp_dir
    assert [path for path in (out_dir / "tmp").iterdir() if path.is_dir()] == []
    (log,) = out_dir.glob("*.log")
    assert "index not found" not in log.read_text()
    assert (out_dir / "fragments_list.txt").stat().st_size > 0


@pytest.mark.parametrize(
    ("aligner", "mapping", "backend", "build"),
    [
        ("bowtie2", "normal", None, ["bowtie2-build", "-q", "{fasta}", "{prefix}"]),
        ("bwa", "normal", "binary", ["bwa", "index", "-p", "{prefix}", "{fasta}"]),
        ("bwa", "5spm", "python", None),
    ],
)
@pytest.mark.parametrize("fasta_name", ["genome.fa", "genome.fa.gz"])
def test_index_named_without_fasta_extension(
    aligner, mapping, backend, build, fasta_name, synthetic_lib, tmp_path
):
    """The index of a FASTA file named after it without its extensions
    (bowtie2-build genome.fa genome, bwa index -p genome genome.fa) is
    aligned to, and no other index is built."""
    programs = [aligner] if build is None else [aligner, build[0]]
    missing = [program for program in programs if shutil.which(program) is None]
    if backend != "python" and missing:
        pytest.skip(f"not on PATH: {', '.join(missing)}")
    if backend == "python" and hbm.pybwa_version() is None:
        pytest.skip("pybwa is not installed")
    ref_dir = tmp_path / "ref.v2"
    ref_dir.mkdir()
    fasta = ref_dir / fasta_name
    compressed = fasta_name.endswith(".gz")
    shutil.copy(synthetic_lib.genome_gz if compressed else synthetic_lib.genome, fasta)
    prefix = ref_dir / "genome"
    if build is None:
        hbm.build_index(str(synthetic_lib.genome), str(prefix))
    else:
        sp.run(
            [arg.format(fasta=synthetic_lib.genome, prefix=prefix) for arg in build],
            check=True,
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL,
        )
    out_dir = tmp_path / "out"
    hpi.full_pipeline(
        genome=str(fasta),
        input1=str(synthetic_lib.reads1),
        input2=str(synthetic_lib.reads2),
        aligner=aligner,
        aligner_backend=backend,
        mapping=mapping,
        enzyme=syn.LIBRARY_ENZYME,
        out_dir=str(out_dir),
        no_cleanup=True,
    )
    (log,) = out_dir.glob("*.log")
    log_text = log.read_text()
    assert f"Found the {aligner} index of {fasta} at {prefix}." in log_text
    assert "index not found" not in log_text
    # No index was built: pybwa only links to the index it loads
    built = [
        path
        for path in (out_dir / "tmp").rglob("*")
        if path.is_file() and not path.is_symlink() and path.suffix in (".bt2", ".bwt", ".sa")
    ]
    assert built == []
    with pysam.AlignmentFile(str(out_dir / "tmp" / "for.bam")) as bam:
        assert any(not read.is_unmapped for read in bam)


@pytest.mark.parametrize("backend", ["binary", "python"])
def test_bwa_mem2_index_is_named(backend, synthetic_lib, tmp_path):
    """A bwa-mem2 index, which bwa cannot read, is named in the log of a run
    that builds a bwa index instead."""
    if backend == "binary" and shutil.which("bwa") is None:
        pytest.skip("bwa is not installed")
    if backend == "python" and hbm.pybwa_version() is None:
        pytest.skip("pybwa is not installed")
    genome = tmp_path / "genome.fa"
    shutil.copy(synthetic_lib.genome, genome)
    # The files of a bwa-mem2 index, which only their names identify here
    for suffix in (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"):
        (tmp_path / f"genome.fa{suffix}").write_bytes(b"")
    hpi.full_pipeline(
        genome=str(genome),
        input1=str(synthetic_lib.reads1),
        input2=str(synthetic_lib.reads2),
        aligner="bwa",
        aligner_backend=backend,
        mapping="5spm",
        enzyme=syn.LIBRARY_ENZYME,
        out_dir=str(tmp_path / "out"),
    )
    (log,) = (tmp_path / "out").glob("*.log")
    log_text = log.read_text()
    assert f"Found a bwa-mem2 index at {genome}, which bwa and pybwa cannot read" in log_text
    assert f"bwa index not found at {genome}" in log_text


def test_genome_neither_fasta_nor_index(tmp_path, caplog):
    """A genome that is neither a fasta file nor an index stops the run
    before it writes anything."""
    out_dir = tmp_path / "out"
    with pytest.raises(SystemExit):
        hpi.full_pipeline(
            genome=str(tmp_path / "missing"),
            input1="test_data/valid.pairs",
            start_stage="pairs",
            aligner="bwa",
            out_dir=str(out_dir),
        )
    assert "The genome must be a fasta file, or the prefix of a bwa index" in caplog.text
    assert not out_dir.exists()


@pytest.mark.parametrize("cmd", [None, "hicstuff pipeline -g ref.fa R1.fq R2.fq"])
def test_log_header_records_the_command(tmp_path, cmd):
    """The command that ran the pipeline follows the date, when it is known."""
    log = tmp_path / "hicstuff.log"
    hpi.generate_log_header(str(log), "R1.fq", "R2.fq", "ref.fa", 5000, cmd=cmd)
    header = log.read_text().splitlines()
    assert header[1].startswith("## date: ")
    if cmd is None:
        assert not any(line.startswith("## cmd:") for line in header)
    else:
        assert header[2] == f"## cmd: {cmd}"


@pytest.mark.needs("bowtie2", "bowtie2-inspect")
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


@pytest.mark.needs("bowtie2", "bowtie2-inspect")
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


@pytest.mark.needs("bowtie2", "bowtie2-inspect")
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


@pytest.mark.needs("bowtie2", "bowtie2-inspect")
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


@pytest.mark.needs("bowtie2", "bowtie2-build")
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
@pytest.mark.needs("bowtie2", "bowtie2-inspect")
def test_full_pipeline_bin(mapping, aligner):
    """Crash Test for the whole pipeline. It starts from bam files that
    test_full_pipeline_frags writes."""
    skip_unless_can_align(aligner)
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
    """Write reads of 20 bases, given as (name, ref_id, pos, reverse, mapq)
    in that order, to a BAM file. ref_id -1 makes an unmapped read. A sixth
    element is the CIGAR string of the alignment, 20M by default."""
    header = {"HD": header_hd or {"VN": "1.6"}, "SQ": [{"SN": n, "LN": ln} for n, ln in BAM_REFS]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for name, ref_id, pos, reverse, mapq, *cigar in reads:
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
                read.cigarstring = cigar[0] if cigar else "20M"
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


def _write_info_contigs(tmp_path):
    """Write tmp_path/info_contigs.txt, of the contigs of BAM_REFS."""
    info_contigs = tmp_path / "info_contigs.txt"
    info_contigs.write_text(
        "contig\tlength\tn_frags\tcumul_length\n"
        + "".join(f"{name}\t{length}\t1\t0\n" for name, length in BAM_REFS)
    )
    return info_contigs


def _bam2pairs_rows(tmp_path, **kwargs):
    """Rows of the pairs bam2pairs writes from tmp_path/for.bam and
    tmp_path/rev.bam."""
    out = tmp_path / "out.pairs"
    hpi.bam2pairs(
        str(tmp_path / "for.bam"),
        str(tmp_path / "rev.bam"),
        str(out),
        str(_write_info_contigs(tmp_path)),
        **kwargs,
    )
    pairs = pd.read_csv(out, sep="\t", comment="#", header=None, dtype={0: str})
    return pairs.values.tolist()


def _run_bam2pairs(tmp_path, mates, **kwargs):
    """Run bam2pairs on tmp_path/for.bam and tmp_path/rev.bam, and check that
    it writes the pairs of the reads of ``mates`` (the alignments of each end,
    by read name) with MAPQ >= 30, flipped to the upper triangle. Mates are at
    their leftmost base or, with five_prime, at their 5' end: the last base of
    the 20 of a reverse-strand alignment."""

    def position(aln):
        return aln[1] + (20 if kwargs.get("five_prime") and aln[2] else 1)

    expected = []
    for name in mates["for"].keys() & mates["rev"].keys():
        aln1, aln2 = mates["for"][name], mates["rev"][name]
        if min(aln1[3], aln2[3]) < 30:
            continue
        if (aln1[0], position(aln1)) > (aln2[0], position(aln2)):
            aln1, aln2 = aln2, aln1
        expected.append(
            [name]
            + [BAM_REFS[aln1[0]][0], position(aln1), BAM_REFS[aln2[0]][0], position(aln2)]
            + ["-" if aln1[2] else "+", "-" if aln2[2] else "+"]
        )
    assert sorted(_bam2pairs_rows(tmp_path, **kwargs)) == sorted(expected)


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


def test_bam2pairs_records_aligned_5_prime_ends(tmp_path):
    """With five_prime, each mate is at its aligned 5' end, the end of the
    sequenced molecule: the first base of a forward-strand alignment and the
    last of a reverse-strand one, soft-clipped bases left out. Pairs are
    flipped to the upper triangle by these positions: the reverse mate of r1
    starts left of its forward mate, but ends right of it. Otherwise, each
    mate is at the leftmost base of its alignment."""
    _write_bam(
        tmp_path / "for.bam",
        [("r1", 0, 100, False, 40), ("r2", 1, 500, True, 40), ("r3", 1, 1000, False, 40, "5S15M")],
    )
    _write_bam(
        tmp_path / "rev.bam",
        [("r1", 0, 90, True, 40), ("r2", 0, 700, False, 40), ("r3", 1, 2000, True, 40, "15M5S")],
    )

    assert _bam2pairs_rows(tmp_path) == [
        ["r1", "chr1", 91, "chr1", 101, "-", "+"],
        ["r2", "chr1", 701, "chr2", 501, "+", "-"],
        ["r3", "chr2", 1001, "chr2", 2001, "+", "-"],
    ]
    assert _bam2pairs_rows(tmp_path, five_prime=True) == [
        ["r1", "chr1", 101, "chr1", 110, "+", "-"],
        ["r2", "chr1", 701, "chr2", 520, "+", "-"],
        ["r3", "chr2", 1001, "chr2", 2015, "+", "-"],
    ]


@pytest.mark.parametrize("order", ["lexicographical", "unsorted"])
def test_bam2pairs_records_5_prime_ends_of_bams_in_another_order(tmp_path, order):
    """Mates are at their 5' ends in BAM files walked again in another order,
    or sorted again."""
    mates = _bams_in_order(tmp_path, order)
    sort_dir = tmp_path / "sort"
    sort_dir.mkdir()

    _run_bam2pairs(tmp_path, mates, tmp_dir=str(sort_dir), five_prime=True)


def _run_from_pairs(tmp_path, name, **kwargs):
    """A pipeline run from a copy of test_data's pairs, in tmp_path/name,
    which it returns. test_data's genome has two contigs: seq1, of 60 kb, and
    seq2, of 20 kb."""
    pairs = tmp_path / f"{name}.pairs"
    shutil.copy("test_data/valid.pairs", pairs)
    out_dir = tmp_path / name
    hpi.full_pipeline(
        genome="test_data/genome/seq.fa",
        input1=str(pairs),
        start_stage="pairs",
        enzyme="DpnII",
        out_dir=str(out_dir),
        **kwargs,
    )
    return out_dir


def _outputs(out_dir):
    """The files that a run wrote in out_dir, but its log, by name: the
    contents of a compressed one uncompressed."""
    found = {}
    for path in sorted(out_dir.iterdir()):
        if path.is_file() and path.suffix != ".log":
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rb") as handle:
                found[path.name] = handle.read()
    return found


@pytest.mark.parametrize(
    "prefix, written, matrix",
    [
        (
            None,
            ["fragments_list.txt", "info_contigs.txt", "valid_idx.pairs.gz"],
            "abs_fragments_contacts_weighted.cool",
        ),
        (
            "sample",
            ["sample.chr.tsv", "sample.frags.tsv", "sample.valid_idx.pairs.gz"],
            "sample.cool",
        ),
    ],
)
def test_full_pipeline_without_matrix(prefix, written, matrix, tmp_path, caplog):
    """--matrix none writes no matrix: the run ends with its pairs file,
    which, with the fragments and contigs tables, is what a run that writes
    the cool matrix writes too, with or without a prefix."""
    caplog.set_level(logging.INFO)
    with_cool = _outputs(_run_from_pairs(tmp_path, "cool", mat_fmt="cool", prefix=prefix))
    without = _outputs(_run_from_pairs(tmp_path, "none", mat_fmt="none", prefix=prefix))
    assert sorted(without) == written
    assert sorted(set(with_cool) - set(without)) == [matrix]
    assert without == {name: with_cool[name] for name in without}
    assert "No matrix (--matrix none): the run ends with its pairs file." in caplog.text
    assert re.search(r"Pairs generated after \d+h \d+m \d+s", caplog.text)


@pytest.mark.parametrize(
    "mat_fmt, kwargs, named",
    [
        ("none", dict(binning=1000), "writes no matrix: binning is an option"),
        (
            "none",
            dict(balancing_args="--max-iters 10"),
            "writes no matrix: balancing_args is an option",
        ),
        ("none", dict(exclude="seq2"), "writes no matrix: exclude is an option"),
        (
            "none",
            dict(binning=1000, exclude="seq2"),
            "writes no matrix: binning, exclude are options",
        ),
        ("none", dict(binning=-5), "writes no matrix: binning is an option"),
        ("none", dict(exclude=["seq2"]), "writes no matrix: exclude is an option"),
        ("metacool", dict(zoomify=True), "is neither zoomified nor balanced: zoomify is an option"),
        (
            "metacool",
            dict(balancing_args="--max-iters 10"),
            "is neither zoomified nor balanced: balancing_args is an option",
        ),
        (
            "metacool",
            dict(binning=1000, zoomify=True, balancing_args="--max-iters 10", exclude="seq2"),
            "is neither zoomified nor balanced: zoomify, balancing_args are options",
        ),
    ],
)
def test_full_pipeline_refuses_cool_options(mat_fmt, kwargs, named, tmp_path, caplog):
    """Binning, zoomify and balancing arguments are options of the cool
    matrix, and contigs to exclude of the cool and metacool matrices: with
    the matrix format none, which writes no matrix, they are refused, and
    zoomify True and balancing arguments with metacool, which is binned but
    neither zoomified nor balanced, before the run writes anything, the
    reason logged as an error."""
    out_dir = tmp_path / "out"
    message = f"The matrix format {mat_fmt} {named} of the cool matrix."
    with pytest.raises(ValueError, match=f"^{re.escape(message)}$"):
        hpi.full_pipeline(
            genome="test_data/genome/seq.fa",
            input1="test_data/valid.pairs",
            start_stage="pairs",
            enzyme="DpnII",
            out_dir=str(out_dir),
            mat_fmt=mat_fmt,
            **kwargs,
        )
    assert not out_dir.exists()
    assert ("ERROR", message) in [
        (record.levelname, record.getMessage()) for record in caplog.records
    ]


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(binning=0, zoomify=False, balancing_args="", exclude=""),
        dict(binning=" 0 ", balancing_args=" ", exclude=" , "),
    ],
)
def test_full_pipeline_without_matrix_takes_options_without_effect(kwargs, tmp_path):
    """No binning (0), no zoomify, blank balancing arguments and no contig to
    exclude would shape no matrix: the matrix format none takes them, and
    writes what a run without them writes, no matrix."""
    plain = _outputs(_run_from_pairs(tmp_path, "plain", mat_fmt="none"))
    assert sorted(plain) == ["fragments_list.txt", "info_contigs.txt", "valid_idx.pairs.gz"]
    assert _outputs(_run_from_pairs(tmp_path, "given", mat_fmt="none", **kwargs)) == plain


def test_full_pipeline_takes_a_list_of_contigs_to_exclude(tmp_path):
    """exclude given as a list of names, which full_pipeline does not
    document, still runs with a format that does not read it. A guard
    against the previous commit, where bg2 never read exclude: it failed
    once the count of the contigs of the cool matrix read it as a string."""
    out_dir = _run_from_pairs(tmp_path, "run", mat_fmt="bg2", exclude=["seq2"])
    assert (out_dir / "abs_fragments_contacts_weighted.bg2").stat().st_size > 0


@pytest.mark.parametrize("mat_fmt", ["hic", "mcool", "None"])
def test_full_pipeline_refuses_an_unknown_matrix_format(mat_fmt, tmp_path, caplog):
    """A matrix format other than graal, bg2, cool, metacool and none is
    refused before the run writes anything, the reason logged as an error.
    Another name wrote an empty matrix file under its extension. None, the
    Python value, is the default format, and the name "None" is unknown."""
    out_dir = tmp_path / "out"
    message = f"Unknown matrix format {mat_fmt!r}: choose graal, bg2, cool, metacool or none."
    with pytest.raises(ValueError, match=f"^{re.escape(message)}$"):
        hpi.full_pipeline(
            genome="test_data/genome/seq.fa",
            input1="test_data/valid.pairs",
            start_stage="pairs",
            enzyme="DpnII",
            out_dir=str(out_dir),
            mat_fmt=mat_fmt,
        )
    assert not out_dir.exists()
    assert ("ERROR", message) in [
        (record.levelname, record.getMessage()) for record in caplog.records
    ]


def _slow_cool_warnings(caplog):
    return [
        record.getMessage()
        for record in caplog.records
        if record.levelname == "WARNING" and "contigs in the cool matrix" in record.getMessage()
    ]


#: The warning of a cool matrix of MANY_CONTIGS contigs or more, for
#: test_data's two contigs
SLOW_COOL = (
    "2 contigs in the cool matrix: it will be slow to write (cooler takes a time that grows "
    "with the square of the number of contigs for every cooler it writes: about 30 s at "
    "117,000 contigs, four times more for twice as many). --matrix metacool writes it with "
    "metacooler, which does not pay that time, with one bin per fragment; --matrix none ends "
    "the run with its pairs file; --matrix bg2 and graal write a text matrix."
)

#: The warning of the default matrix, metacool from MANY_CONTIGS contigs, for
#: test_data's two contigs
DEFAULT_METACOOL = (
    "2 contigs in the matrix: it will be a metacool matrix, written with metacooler, rather "
    "than a cool matrix, which cooler is slow to write with this many contigs (cooler takes a "
    "time that grows with the square of the number of contigs for every cooler it writes: "
    "about 30 s at 117,000 contigs, four times more for twice as many). --matrix cool writes "
    "the cool matrix all the same."
)


def _warnings(caplog):
    return [record.getMessage() for record in caplog.records if record.levelname == "WARNING"]


@pytest.mark.parametrize(
    "mat_fmt, many, warned",
    [
        ("cool", 2, True),
        ("cool", 3, False),
        ("bg2", 2, False),
        ("none", 2, False),
        ("metacool", 2, False),
    ],
)
def test_full_pipeline_warns_of_a_slow_cool_matrix(
    mat_fmt, many, warned, tmp_path, caplog, monkeypatch
):
    """A cool matrix of MANY_CONTIGS contigs or more makes the run warn that
    it will be slow to write, and which formats do not pay that time:
    test_data's genome has two contigs. Fewer contigs, or another format,
    make no warning."""
    monkeypatch.setattr(hpi, "MANY_CONTIGS", many)
    _run_from_pairs(tmp_path, "run", mat_fmt=mat_fmt)
    assert _slow_cool_warnings(caplog) == ([SLOW_COOL] if warned else [])
    assert DEFAULT_METACOOL not in _warnings(caplog)


@pytest.mark.parametrize(
    "kwargs, chroms",
    [
        (dict(), ["seq1", "seq2"]),
        (dict(min_size=30_000), ["seq1"]),
        (dict(exclude="seq2"), ["seq1"]),
    ],
)
def test_full_pipeline_counts_the_contigs_of_the_cool_matrix(
    kwargs, chroms, tmp_path, caplog, monkeypatch
):
    """The warning counts the contigs of the cool matrix, not those of the
    genome: --size and --exclude leave some out."""
    monkeypatch.setattr(hpi, "MANY_CONTIGS", 1)
    out_dir = _run_from_pairs(tmp_path, "run", mat_fmt="cool", **kwargs)
    (warning,) = _slow_cool_warnings(caplog)
    assert warning.startswith(f"{len(chroms)} contigs in the cool matrix: ")
    assert cooler.Cooler(str(out_dir / "abs_fragments_contacts_weighted.cool")).chromnames == chroms


@pytest.mark.parametrize(
    "mat_fmt, warning",
    [
        ("cool", "contigs in the cool matrix"),
        (None, "contigs in the matrix: it will be a metacool"),
    ],
)
def test_full_pipeline_warns_before_building_an_index(
    mat_fmt, warning, tmp_path, caplog, monkeypatch
):
    """The warnings of a slow cool matrix, and of the default metacool matrix,
    come once the genome is read, before the index of a fasta file without
    one is built, which takes long on a large assembly: the run is stopped
    there."""
    genome = tmp_path / "seq.fa"
    shutil.copy("test_data/genome/seq.fa", genome)
    monkeypatch.setattr(hpi, "MANY_CONTIGS", 2)
    run = sp.run

    class IndexBuilt(Exception):
        pass

    def stop_at_the_index(args, *rest, **kwargs):
        if args[0] == "bowtie2-build":
            raise IndexBuilt
        return run(args, *rest, **kwargs)

    monkeypatch.setattr(hpi.sp, "run", stop_at_the_index)
    # The run stops before it runs a program: none need be installed (CI's
    # pip-only job runs the tests without any).
    monkeypatch.setattr(hpi.hio, "require_tools", lambda *args, **kwargs: None)
    caplog.set_level(logging.INFO)
    with pytest.raises(IndexBuilt):
        hpi.full_pipeline(
            genome=str(genome),
            input1="test_data/sample.reads_for.fastq.gz",
            input2="test_data/sample.reads_rev.fastq.gz",
            enzyme="DpnII",
            out_dir=str(tmp_path / "out"),
            aligner="bowtie2",
            mat_fmt=mat_fmt,
        )
    messages = [record.getMessage() for record in caplog.records]
    (warned,) = [i for i, message in enumerate(messages) if warning in message]
    (indexing,) = [i for i, message in enumerate(messages) if "index not found" in message]
    assert warned < indexing


def _contacts(path):
    """The contacts of a cool or metacool matrix by pair of bins, a bin being
    its (contig, start, end), whatever the order of the contigs in the file.
    cooler names the contig of a bin of a cool matrix, and gives the number of
    that of a metacool matrix, which metacooler stores as an integer."""
    clr = cooler.Cooler(str(path))
    bins = clr.bins()[:]
    chroms = bins.chrom.astype(str)
    if bins.chrom.dtype.kind in "iu":
        chroms = np.asarray(clr.chromnames)[bins.chrom]
    keys = list(zip(chroms, bins.start, bins.end))
    contacts = {}
    for bin1, bin2, count in clr.pixels()[:][["bin1_id", "bin2_id", "count"]].itertuples(
        index=False
    ):
        pair = tuple(sorted((keys[bin1], keys[bin2])))
        contacts[pair] = contacts.get(pair, 0) + int(count)
    return contacts


@pytest.mark.parametrize(
    "prefix, matrix", [(None, "abs_fragments_contacts_weighted"), ("sample", "sample")]
)
def test_full_pipeline_metacool_matrix(prefix, matrix, tmp_path):
    """--matrix metacool writes a metacool matrix, with metacooler: the bins
    of the cool matrix, one per fragment, and its contacts. The other outputs
    are those of the cool run, byte for byte, with or without a prefix.
    metacooler checks the file, which cooler reads too; its contigs are not
    yet assigned to genomes, and the run's command is stored as the one that
    made it."""
    cool_dir = _run_from_pairs(tmp_path, "cool", mat_fmt="cool", prefix=prefix)
    cmd = "hicstuff pipeline -M metacool -S pairs"
    meta_dir = _run_from_pairs(tmp_path, "meta", mat_fmt="metacool", prefix=prefix, cmd=cmd)
    with_cool, with_metacool = _outputs(cool_dir), _outputs(meta_dir)
    assert sorted(set(with_cool) ^ set(with_metacool)) == [f"{matrix}.cool", f"{matrix}.metacool"]
    common = set(with_cool) & set(with_metacool)
    assert {name: with_metacool[name] for name in common} == {
        name: with_cool[name] for name in common
    }
    metacool = str(meta_dir / f"{matrix}.metacool")
    contacts = _contacts(metacool)
    assert contacts == _contacts(cool_dir / f"{matrix}.cool")
    assert sum(contacts.values()) == 10_000
    assert all(check.ok for check in metacooler.check(metacool))
    mc = metacooler.MetaCooler(metacool)
    assert mc.info["bin-type"] == "variable"
    assert mc.genomes.name.tolist() == ["__unbinned__"]
    assert mc.contigs.name.tolist() == ["seq1", "seq2"]
    assert mc.metadata["commands"] == [cmd]


def test_full_pipeline_metacool_matrix_excludes_contigs(tmp_path):
    """--exclude leaves its contigs out of the metacool matrix, with their
    contacts, as out of the cool matrix. A run without a command stores
    none."""
    cool_dir = _run_from_pairs(tmp_path, "cool", mat_fmt="cool", exclude="seq2")
    meta_dir = _run_from_pairs(tmp_path, "meta", mat_fmt="metacool", exclude=["seq2"])
    metacool = str(meta_dir / "abs_fragments_contacts_weighted.metacool")
    contacts = _contacts(metacool)
    assert contacts == _contacts(cool_dir / "abs_fragments_contacts_weighted.cool")
    assert {bin1[0] for bin1, _ in contacts} == {"seq1"}
    mc = metacooler.MetaCooler(metacool)
    assert mc.contigs.name.tolist() == ["seq1"]
    assert "commands" not in mc.metadata


@pytest.mark.parametrize(
    "many, kwargs, written",
    [
        (2, dict(), "metacool"),
        (3, dict(), "cool"),
        (2, dict(min_size=30_000), "cool"),
        (2, dict(exclude="seq2"), "cool"),
        (1, dict(exclude="seq2"), "metacool"),
    ],
)
def test_full_pipeline_default_matrix(many, kwargs, written, tmp_path, caplog, monkeypatch):
    """Without a matrix format, the run writes a cool matrix, or, from
    MANY_CONTIGS contigs in the matrix (those that --size and --exclude keep),
    a metacool matrix, and warns of it: test_data's genome has two contigs."""
    monkeypatch.setattr(hpi, "MANY_CONTIGS", many)
    out_dir = _run_from_pairs(tmp_path, "run", **kwargs)
    matrices = [path.name for path in out_dir.iterdir() if path.suffix in (".cool", ".metacool")]
    assert matrices == [f"abs_fragments_contacts_weighted.{written}"]
    n_contigs = 1 if kwargs else 2
    default_warnings = [
        message
        for message in _warnings(caplog)
        if message.startswith(f"{n_contigs} contigs in the")
    ]
    if written == "metacool":
        assert default_warnings == [DEFAULT_METACOOL.replace("2 contigs", f"{n_contigs} contigs")]
        assert metacooler.is_metacooler(str(out_dir / matrices[0]))
    else:
        assert default_warnings == []


#: The warning of the default matrix, metacool from MANY_CONTIGS contigs,
#: binned, for test_data's two contigs
DEFAULT_BINNED_METACOOL = DEFAULT_METACOOL.replace(
    " --matrix cool", " It is binned, but neither zoomified nor balanced. --matrix cool"
)


@pytest.mark.parametrize("zoomify", [None, False])
def test_full_pipeline_default_binned_matrix(zoomify, tmp_path, caplog, monkeypatch):
    """With binning too, the default matrix is metacool from MANY_CONTIGS
    contigs: binned, but neither zoomified nor balanced, and the run warns of
    it. zoomify left to its default, or False, is taken."""
    monkeypatch.setattr(hpi, "MANY_CONTIGS", 2)
    out_dir = _run_from_pairs(tmp_path, "run", binning=5000, zoomify=zoomify)
    assert not list(out_dir.glob("*.cool")) and not list(out_dir.glob("*.mcool"))
    metacool = str(out_dir / "abs_fragments_contacts_weighted.metacool")
    assert cooler.Cooler(metacool).binsize == 5000
    assert "weight" not in cooler.Cooler(metacool).bins().columns
    assert DEFAULT_BINNED_METACOOL in _warnings(caplog)
    assert _slow_cool_warnings(caplog) == []


def test_full_pipeline_default_binned_matrix_below_many_contigs(tmp_path, caplog, monkeypatch):
    """Below MANY_CONTIGS contigs, the default binned matrix is a cool one,
    zoomified and balanced, as before, and the run warns of no format."""
    monkeypatch.setattr(hpi, "MANY_CONTIGS", 3)
    out_dir = _run_from_pairs(tmp_path, "run", binning=5000)
    assert (out_dir / "abs_fragments_contacts_weighted.mcool").is_file()
    assert not list(out_dir.glob("*.metacool"))
    assert _slow_cool_warnings(caplog) == []
    assert DEFAULT_BINNED_METACOOL not in _warnings(caplog)


def _bed(path, rows):
    """A BED file of bins at path, of rows of (contig, start, end)."""
    path.write_text("".join(f"{chrom}\t{start}\t{end}\n" for chrom, start, end in rows))
    return str(path)


def _bins(path):
    """The bins of a cool or metacool matrix, as (contig, start, end)."""
    clr = cooler.Cooler(str(path))
    bins = clr.bins()[:]
    chroms = bins.chrom.astype(str)
    if bins.chrom.dtype.kind in "iu":
        chroms = np.asarray(clr.chromnames)[bins.chrom]
    return sorted(zip(chroms, bins.start, bins.end))


#: test_data's contigs and their lengths
LENGTHS = {"seq1": 60_000, "seq2": 20_000}

#: Bins of no single size: of 3 kb on seq2, then of 7 kb on seq1
UNEVEN_BINS = [("seq2", start, min(start + 3000, 20_000)) for start in range(0, 20_000, 3000)] + [
    ("seq1", start, min(start + 7000, 60_000)) for start in range(0, 60_000, 7000)
]


@pytest.mark.parametrize("bed", [False, True])
def test_full_pipeline_binned_metacool_matrix(bed, tmp_path):
    """-M metacool with binning, a bin size or a BED file of bins: the pairs
    are binned from their positions, as for the cool matrix. The metacool
    matrix has the bins of the cool matrix of the same binning, and its
    contacts, bin for bin, with metacooler's contact statistics counted on
    the pairs; it is neither zoomified nor balanced."""
    binning = _bed(tmp_path / "bins.bed", UNEVEN_BINS) if bed else 5000
    cool_dir = _run_from_pairs(tmp_path, "cool", mat_fmt="cool", binning=binning, zoomify=False)
    meta_dir = _run_from_pairs(tmp_path, "meta", mat_fmt="metacool", binning=binning)
    cool = cool_dir / "abs_fragments_contacts_weighted.cool"
    metacool = str(meta_dir / "abs_fragments_contacts_weighted.metacool")
    assert _bins(metacool) == _bins(cool)
    if bed:
        assert _bins(metacool) == sorted(UNEVEN_BINS)
    contacts = _contacts(metacool)
    assert contacts == _contacts(cool)
    assert sum(contacts.values()) == 10_000
    assert all(check.ok for check in metacooler.check(metacool))
    assert cooler.Cooler(metacool).binsize == (None if bed else 5000)
    assert "weight" in cooler.Cooler(str(cool)).bins().columns
    assert "weight" not in cooler.Cooler(metacool).bins().columns
    qc = metacooler.MetaCooler(metacool).metadata["qc"]
    assert (qc["source"], qc["distance_threshold"]) == ("pairs", 1001)
    assert not list(meta_dir.glob("*.mcool"))


def test_full_pipeline_metacool_matrix_of_contigs(tmp_path):
    """A BED file of one bin per contig makes the map of the contigs, as
    metaTOR reads it: its contacts are those of the fragments, summed by pair
    of contigs."""
    binning = _bed(
        tmp_path / "contigs.bed", [(name, 0, length) for name, length in LENGTHS.items()]
    )
    frag_dir = _run_from_pairs(tmp_path, "frag", mat_fmt="metacool")
    meta_dir = _run_from_pairs(tmp_path, "meta", mat_fmt="metacool", binning=binning)
    by_contig = {}
    for pair, count in _contacts(frag_dir / "abs_fragments_contacts_weighted.metacool").items():
        key = tuple(sorted((name, 0, LENGTHS[name]) for name, _, _ in pair))
        by_contig[key] = by_contig.get(key, 0) + count
    metacool = meta_dir / "abs_fragments_contacts_weighted.metacool"
    assert _contacts(metacool) == by_contig
    assert _bins(metacool) == [("seq1", 0, 60_000), ("seq2", 0, 20_000)]


@pytest.mark.parametrize(
    "kwargs", [dict(exclude="seq2"), dict(exclude=["seq2"]), dict(min_size=30_000)]
)
@pytest.mark.parametrize("mat_fmt", ["cool", "metacool"])
def test_full_pipeline_bins_on_the_contigs_of_the_matrix(mat_fmt, kwargs, tmp_path):
    """The bins of a BED file, or of a bin size, on the contigs that --exclude
    or --size leaves out are left out of the matrix, with their contacts."""
    if mat_fmt == "cool" and isinstance(kwargs.get("exclude"), list):
        pytest.skip("the cool matrix reads comma-separated contigs to exclude only")
    for binning in (_bed(tmp_path / "bins.bed", UNEVEN_BINS), 5000):
        out_dir = _run_from_pairs(
            tmp_path,
            f"{mat_fmt}{binning != 5000}",
            mat_fmt=mat_fmt,
            binning=binning,
            zoomify=False,
            **kwargs,
        )
        (matrix,) = [path for path in out_dir.iterdir() if path.suffix in (".cool", ".metacool")]
        assert {chrom for chrom, _, _ in _bins(matrix)} == {"seq1"}
        assert {pair[0][0] for pair in _contacts(matrix)} == {"seq1"}


def test_full_pipeline_bed_binned_cool_matrix(tmp_path):
    """The cool matrix of a BED file of bins is balanced, but not zoomified:
    zoomify left to its default is ignored, as a BED file gives no bin
    size."""
    out_dir = _run_from_pairs(tmp_path, "run", binning=_bed(tmp_path / "bins.bed", UNEVEN_BINS))
    assert not list(out_dir.glob("*.mcool"))
    cool = out_dir / "abs_fragments_contacts_weighted.cool"
    assert _bins(cool) == sorted(UNEVEN_BINS)
    assert "weight" in cooler.Cooler(str(cool)).bins().columns


def _refused_before_the_run(tmp_path, caplog, message, **kwargs):
    """full_pipeline refuses kwargs with message, logged as an error, before
    it writes anything."""
    out_dir = tmp_path / "out"
    with pytest.raises(ValueError, match=f"^{re.escape(message)}$"):
        hpi.full_pipeline(
            genome="test_data/genome/seq.fa",
            input1="test_data/valid.pairs",
            start_stage="pairs",
            enzyme="DpnII",
            out_dir=str(out_dir),
            **kwargs,
        )
    assert not out_dir.exists()
    assert ("ERROR", message) in [(r.levelname, r.getMessage()) for r in caplog.records]


@pytest.mark.parametrize("binning", [-5, "-5", " -5 "])
def test_full_pipeline_refuses_a_negative_binning(binning, tmp_path, caplog):
    message = "The binning must be 0, a bin size in bp or a BED file of bins, not -5."
    _refused_before_the_run(tmp_path, caplog, message, mat_fmt="cool", binning=binning)


def test_full_pipeline_refuses_a_binning_that_is_no_file(tmp_path, caplog):
    path = str(tmp_path / "nowhere.bed")
    message = (
        f"The binning {path!r} is neither a bin size in bp nor a BED file of bins: there is no "
        "such file."
    )
    _refused_before_the_run(tmp_path, caplog, message, mat_fmt="metacool", binning=path)


@pytest.mark.parametrize(
    "content, reason",
    [
        ("", "it is empty."),
        (
            "chrom\tstart\tend\nseq1\t0\t100\n",
            "its second and third columns, the starts and ends of the bins, must be integers, "
            "without a header line.",
        ),
        ("seq1\t100\t100\n", "the bin seq1:100-100 is empty."),
        ("seq1\t-1\t100\n", "the first bin of seq1 starts at -1, not 0."),
        ("seq1\t0\t100\nseq2\t100\t200\n", "the first bin of seq2 starts at 100, not 0."),
        (
            "seq1\t100\t300\nseq2\t0\t10\nseq1\t0\t200\n",
            "the bins seq1:0-200 and seq1:100-300 overlap.",
        ),
        (
            "seq1\t0\t100\nseq1\t200\t300\n",
            "the bins seq1:0-100 and seq1:200-300 leave a gap.",
        ),
    ],
)
def test_full_pipeline_refuses_a_bed_file_that_is_not_of_bins(content, reason, tmp_path, caplog):
    """A BED file of bins is read before the run writes anything: one that is
    not of bins, or whose bins do not tile their contigs from 0, without gap
    or overlap, is refused, and why."""
    bed = tmp_path / "bins.bed"
    bed.write_text(content)
    message = f"{bed} is not a BED file of bins: {reason}"
    _refused_before_the_run(tmp_path, caplog, message, mat_fmt="metacool", binning=str(bed))


def test_full_pipeline_refuses_a_bed_file_of_fewer_than_three_columns(tmp_path, caplog):
    bed = tmp_path / "bins.bed"
    bed.write_text("seq1\t0\n")
    out_dir = tmp_path / "out"
    with pytest.raises(ValueError, match=f"^{re.escape(str(bed))} is not a BED file of bins: "):
        hpi.full_pipeline(
            genome="test_data/genome/seq.fa",
            input1="test_data/valid.pairs",
            start_stage="pairs",
            enzyme="DpnII",
            out_dir=str(out_dir),
            mat_fmt="cool",
            binning=str(bed),
        )
    assert not out_dir.exists()


@pytest.mark.parametrize("mat_fmt", ["cool", None])
def test_full_pipeline_refuses_zoomify_with_a_bed_file(mat_fmt, tmp_path, caplog):
    """zoomify True is refused with the bins of a BED file, which give no bin
    size, before the run writes anything: with the cool matrix, and with the
    default one, which is cool or metacool, which is not zoomified."""
    bed = _bed(tmp_path / "bins.bed", UNEVEN_BINS)
    message = hpi._zoomify_refused_message(bed)
    _refused_before_the_run(tmp_path, caplog, message, mat_fmt=mat_fmt, binning=bed, zoomify=True)


@pytest.mark.parametrize(
    "rows, reason",
    [
        ([("seq3", 0, 100)], "are on contigs that are not in the genome: seq3."),
        (
            [(f"other{i}", 0, 100) for i in range(7)],
            "are on contigs that are not in the genome: other0, other1, other2, other3, "
            "other4, and 2 more.",
        ),
        ([("seq2", 0, 25_000)], "ends at 25,000 bp, not at the end of the contig, 20,000 bp."),
        (
            [("seq1", 0, 60_000), ("seq2", 0, 10_000)],
            "ends at 10,000 bp, not at the end of the contig, 20,000 bp.",
        ),
    ],
)
def test_full_pipeline_refuses_bins_that_are_not_of_the_genome(rows, reason, tmp_path, caplog):
    """The bins of a BED file must be on contigs of the genome, and those of a
    contig must end at its end, as they tile it: checked once the genome is
    read, before the pairs are."""
    bed = _bed(tmp_path / "bins.bed", rows)
    with pytest.raises(ValueError) as raised:
        _run_from_pairs(tmp_path, "run", mat_fmt="metacool", binning=bed)
    message = str(raised.value)
    assert message.endswith(reason)
    assert bed in message
    assert ("ERROR", message) in [(r.levelname, r.getMessage()) for r in caplog.records]
    assert _outputs(tmp_path / "run") == {}


def test_full_pipeline_refuses_bins_on_no_contig_of_the_matrix(tmp_path):
    bed = _bed(tmp_path / "bins.bed", [("seq2", 0, 20_000)])
    with pytest.raises(ValueError, match="is on a contig of the matrix"):
        _run_from_pairs(tmp_path, "run", mat_fmt="cool", binning=bed, exclude="seq2")
    assert _outputs(tmp_path / "run") == {}


@pytest.mark.parametrize(
    "kwargs, options",
    [
        (dict(balancing_args="--max-iters 10"), ["balancing_args"]),
        (dict(zoomify=True), ["zoomify"]),
        (dict(zoomify=True, balancing_args="--max-iters 10"), ["zoomify", "balancing_args"]),
    ],
)
def test_full_pipeline_refuses_cool_options_with_the_default_metacool(
    kwargs, options, tmp_path, caplog, monkeypatch
):
    """From MANY_CONTIGS contigs, the default matrix is metacool, binned or
    not: zoomify True and balancing arguments, options of the cool matrix, are
    refused with it, once the genome is read, before the pairs are, with
    DefaultMatrixError, which names them and says to give the cool format."""
    monkeypatch.setattr(hpi, "MANY_CONTIGS", 2)
    with pytest.raises(hpi.DefaultMatrixError) as raised:
        _run_from_pairs(tmp_path, "run", binning=5000, **kwargs)
    assert (raised.value.n_contigs, raised.value.options) == (2, options)
    message = hpi._default_metacool_message(2, options)
    assert str(raised.value) == message
    assert message.startswith("2 contigs in the matrix: its default format is metacool, which ")
    assert ("ERROR", message) in [(r.levelname, r.getMessage()) for r in caplog.records]
    assert _outputs(tmp_path / "run") == {}


def test_full_pipeline_warns_of_a_slow_binned_cool_matrix(tmp_path, caplog, monkeypatch):
    """The warning of a slow cool matrix says that the metacool matrix would
    have the same bins, neither zoomified nor balanced."""
    monkeypatch.setattr(hpi, "MANY_CONTIGS", 2)
    _run_from_pairs(tmp_path, "run", mat_fmt="cool", binning=5000, zoomify=False)
    binned = SLOW_COOL.replace(
        "with one bin per fragment", "with the same bins, but neither zoomified nor balanced"
    )
    assert _slow_cool_warnings(caplog) == [binned]


@pytest.mark.parametrize(
    "binning, parsed",
    [(0, 0), ("0", 0), (" 0 ", 0), (" 1000 ", 1000), (np.int64(5), 5), (5000, 5000)],
)
def test_parse_binning(binning, parsed):
    assert hpi.parse_binning(binning) == parsed
    assert type(hpi.parse_binning(binning)) is int
    assert hpi._binned(binning) is (parsed != 0)


def test_parse_binning_of_a_bed_file(tmp_path):
    bed = _bed(tmp_path / "bins.bed", UNEVEN_BINS)
    import pathlib

    assert hpi.parse_binning(pathlib.Path(bed)) == bed
    assert hpi.parse_binning(bed) == bed
    assert hpi._binned(bed)


def test_matrix_bins_cuts_the_contigs_as_cooler_does(tmp_path):
    """The bins of a bin size are cooler's: the contigs of info_contigs in its
    order, but those of exclude, each cut from 0, its last bin shorter."""
    info_contigs = tmp_path / "info_contigs.txt"
    info_contigs.write_text(
        "contig\tlength\tn_frags\tcumul_length\n"
        "b\t2500\t1\t0\nx\t900\t1\t1\na\t1000\t1\t2\n10\t1\t1\t3\n"
    )
    chromsizes = pd.Series({"b": 2500, "a": 1000, "10": 1})
    expected = cooler.util.binnify(chromsizes, 1000)
    got = hpi.matrix_bins(1000, str(info_contigs), "x")
    assert list(zip(got.chrom, got.start, got.end)) == list(
        zip(expected.chrom.astype(str), expected.start, expected.end)
    )
    bed = _bed(
        tmp_path / "bins.bed", [("a", 0, 1000), ("x", 0, 900), ("b", 1000, 2500), ("b", 0, 1000)]
    )
    got = hpi.matrix_bins(bed, str(info_contigs), ["x"])
    assert list(zip(got.chrom, got.start, got.end)) == [
        ("b", 0, 1000),
        ("b", 1000, 2500),
        ("a", 0, 1000),
    ]


def test_binned_metacool_needs_the_contigs(tmp_path):
    with pytest.raises(ValueError, match="needs info_contigs"):
        hpi.pairs2metacool(
            "test_data/valid.pairs", str(tmp_path / "x.metacool"), None, None, binning=1000
        )


def test_full_pipeline_default_metacool_matrix_from_fastq(tmp_path, caplog, monkeypatch):
    """A run from FASTQ files that writes the default metacool matrix ends
    with the statistics of the run: metacooler's log, in the run's, does not
    stop them. pybwa aligns the reads: no program is needed."""
    monkeypatch.setattr(hpi, "MANY_CONTIGS", 2)
    out_dir = tmp_path / "out"
    caplog.set_level(logging.INFO)
    hpi.full_pipeline(
        genome="test_data/genome/seq.fa",
        input1="test_data/sample.reads_for.fastq.gz",
        input2="test_data/sample.reads_rev.fastq.gz",
        enzyme="DpnII",
        aligner="bwa",
        aligner_backend="python",
        out_dir=str(out_dir),
    )
    (log,) = out_dir.glob("hicstuff_*.log")
    log_text = log.read_text()
    assert "Created metacool file" in log_text
    assert "Stats not compiled" not in log_text
    stats = hcs.get_pipeline_stats(str(log))
    contacts = _contacts(out_dir / "abs_fragments_contacts_weighted.metacool")
    assert sum(contacts.values()) == stats["Final contacts"]
