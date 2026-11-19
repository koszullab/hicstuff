"""Several pairs of fastq files, such as the lanes of a library, read as one
library by ``hicstuff pipeline``."""

import gzip
import logging
import os

import cooler
import pandas as pd
import pytest
from click.testing import CliRunner

import hicstuff.pipeline as hpi
import hicstuff.stats as hcs
from hicstuff.cli import cli
from tests.conftest import skip_unless_can_align

READS = ("test_data/sample.reads_for.fastq.gz", "test_data/sample.reads_rev.fastq.gz")
GENOME = "test_data/genome/seq.fa"
#: Reads in the first lane of the split library, of 10000
LANE1_READS = 3000


def _fastq_lines(path):
    with gzip.open(path, "rt") as handle:
        return handle.readlines()


@pytest.fixture
def lanes(tmp_path):
    """The test library split in two lanes: lib_L1_R1.fq.gz, lib_L1_R2.fq.gz,
    lib_L2_R1.fq.gz and lib_L2_R2.fq.gz."""
    lane_dir = tmp_path / "lanes"
    lane_dir.mkdir()
    for mate, path in zip(("R1", "R2"), READS):
        lines = _fastq_lines(path)
        split = 4 * LANE1_READS
        for lane, chunk in (("L1", lines[:split]), ("L2", lines[split:])):
            with gzip.open(lane_dir / f"lib_{lane}_{mate}.fq.gz", "wt") as handle:
                handle.writelines(chunk)
    return lane_dir


def _lane_files(lanes, mate):
    return [str(lanes / f"lib_{lane}_{mate}.fq.gz") for lane in ("L1", "L2")]


def _write_fastq(path, names, suffix=""):
    """A fastq file of reads with these names, each of 4 bases."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as handle:
        for name in names:
            handle.write(f"@{name}{suffix} 1:N:0:0\nACGT\n+\nFFFF\n")
    return str(path)


# ---------------------------------------------------------------------------
# expand_inputs and split_inputs: what the arguments name
# ---------------------------------------------------------------------------


def test_expand_inputs_takes_an_existing_path_as_it_is(tmp_path):
    """A file whose name holds commas or pattern characters is that file."""
    path = tmp_path / "odd,name[1]*.fq"
    path.write_text("")
    assert hpi.expand_inputs(str(path)) == [str(path)]


def test_expand_inputs_splits_commas_and_expands_patterns(tmp_path):
    for name in ("c_R1.fq", "a_R1.fq", "b_R1.fq", "a_R2.fq"):
        (tmp_path / name).write_text("")
    arg = f"{tmp_path}/*_R1.fq,{tmp_path}/z_R1.fq,{tmp_path}/a_R?.fq"
    assert hpi.expand_inputs(arg) == [
        str(tmp_path / "a_R1.fq"),
        str(tmp_path / "b_R1.fq"),
        str(tmp_path / "c_R1.fq"),
        # A path without pattern is kept, found or not
        str(tmp_path / "z_R1.fq"),
        str(tmp_path / "a_R1.fq"),
        str(tmp_path / "a_R2.fq"),
    ]


@pytest.mark.parametrize(
    "arg, message",
    [
        ("{tmp}/*_R3.fq", "No file matches the input pattern"),
        ("{tmp}/a_R1.fq,,{tmp}/b_R1.fq", "Empty path in the input list"),
    ],
)
def test_expand_inputs_refuses(tmp_path, arg, message):
    with pytest.raises(ValueError, match=message):
        hpi.expand_inputs(arg.format(tmp=tmp_path))


@pytest.mark.parametrize(
    "args, stage, expected",
    [
        (["R1.fq", "R2.fq"], "fastq", (["R1.fq"], ["R2.fq"])),
        (["a1.fq,b1.fq", "a2.fq,b2.fq"], "fastq", (["a1.fq", "b1.fq"], ["a2.fq", "b2.fq"])),
        # As the shell expands unquoted patterns: the R1 files, then the R2 files
        (["a1.fq", "b1.fq", "a2.fq", "b2.fq"], "fastq", (["a1.fq", "b1.fq"], ["a2.fq", "b2.fq"])),
        (["R1.bam", "R2.bam"], "bam", (["R1.bam"], ["R2.bam"])),
        (["R1.fq"], "fastq", (["R1.fq"], [])),
        (["valid.pairs"], "pairs", (["valid.pairs"], [])),
        (["a.pairs", "b.pairs"], "pairs_idx", (["a.pairs", "b.pairs"], [])),
    ],
)
def test_split_inputs(args, stage, expected):
    assert hpi.split_inputs(args, stage) == expected


def test_split_inputs_expands_patterns(lanes):
    assert hpi.split_inputs([f"{lanes}/lib_*_R1.fq.gz", f"{lanes}/lib_*_R2.fq.gz"]) == (
        _lane_files(lanes, "R1"),
        _lane_files(lanes, "R2"),
    )


def test_split_inputs_refuses_an_odd_number_of_files():
    with pytest.raises(ValueError, match="3 input files: give the R1 files, then the R2"):
        hpi.split_inputs(["a1.fq", "b1.fq", "a2.fq"])


# ---------------------------------------------------------------------------
# check_fastq_pairs: R1 and R2 files that go together
# ---------------------------------------------------------------------------


def test_check_fastq_pairs_accepts_the_lanes_of_a_library(lanes):
    hpi.check_fastq_pairs(_lane_files(lanes, "R1"), _lane_files(lanes, "R2"))


def test_check_fastq_pairs_reads_names_ending_in_1_and_2(tmp_path):
    """Reads named read/1 and read/2 are the mates of a pair, in plain files
    as in gzip files."""
    forward = [_write_fastq(tmp_path / f"{lane}_R1.fq", ["r1", "r2"], "/1") for lane in "ab"]
    reverse = [_write_fastq(tmp_path / f"{lane}_R2.fq.gz", ["r1", "r2"], "/2") for lane in "ab"]
    hpi.check_fastq_pairs(forward, reverse)


@pytest.mark.parametrize(
    "content, name",
    [("@r1/1 1:N:0:0\n", "r1"), ("@r1/2\n", "r1"), ("@r1/3\n", "r1/3"), ("", None), ("@\n", None)],
)
def test_first_read_name(tmp_path, content, name):
    """Without its /1 or /2; None for a file without reads, or a read without
    a name."""
    path = tmp_path / "reads.fq"
    path.write_text(content)
    assert hpi._first_read_name(str(path)) == name


def test_check_fastq_pairs_refuses_files_of_other_reads(lanes):
    forward = _lane_files(lanes, "R1")
    with pytest.raises(ValueError, match="are not the R1 and R2 files of the same reads"):
        hpi.check_fastq_pairs(forward, _lane_files(lanes, "R2")[::-1])


def test_check_fastq_pairs_refuses_a_missing_file(lanes):
    forward = _lane_files(lanes, "R1")
    reverse = [_lane_files(lanes, "R2")[0], str(lanes / "missing_R2.fq.gz")]
    with pytest.raises(ValueError, match="Input file not found: .*missing_R2.fq.gz"):
        hpi.check_fastq_pairs(forward, reverse)


def test_check_fastq_pairs_refuses_a_file_given_twice(lanes):
    """Through a symbolic link too."""
    forward = _lane_files(lanes, "R1")
    link = lanes / "link_R1.fq.gz"
    os.symlink(forward[0], link)
    with pytest.raises(ValueError, match="Input file given twice"):
        hpi.check_fastq_pairs([forward[0], str(link)], _lane_files(lanes, "R2"))


# ---------------------------------------------------------------------------
# concatenate_fastq
# ---------------------------------------------------------------------------


def test_concatenate_fastq_copies_gzip_files_as_they_are(lanes, tmp_path):
    files = _lane_files(lanes, "R1")
    out = tmp_path / "for.fq.gz"
    hpi.concatenate_fastq(files, str(out))
    assert out.read_bytes() == b"".join(open(path, "rb").read() for path in files)
    assert _fastq_lines(out) == _fastq_lines(READS[0])


def test_concatenate_fastq_of_plain_files_is_plain(tmp_path):
    """A plain file without a newline at its end gets one."""
    first = tmp_path / "a.fq"
    first.write_text("@r1\nACGT\n+\nFFFF")
    second = _write_fastq(tmp_path / "b.fq", ["r2"])
    out = tmp_path / "for.fq"
    hpi.concatenate_fastq([str(first), second], str(out))
    assert out.read_text() == "@r1\nACGT\n+\nFFFF\n@r2 1:N:0:0\nACGT\n+\nFFFF\n"


def test_concatenate_fastq_compresses_plain_files_among_gzip_files(tmp_path):
    files = [
        _write_fastq(tmp_path / "a.fq", ["r1"]),
        _write_fastq(tmp_path / "b.fq.gz", ["r2"]),
        _write_fastq(tmp_path / "c.fq", ["r3"]),
    ]
    out = tmp_path / "for.fq.gz"
    hpi.concatenate_fastq(files, str(out))
    assert "".join(_fastq_lines(out)) == "".join(
        f"@{name} 1:N:0:0\nACGT\n+\nFFFF\n" for name in ("r1", "r2", "r3")
    )


# ---------------------------------------------------------------------------
# full_pipeline and the command line
# ---------------------------------------------------------------------------


def _pairs_rows(out_dir):
    with gzip.open(out_dir / "valid_idx.pairs.gz", "rt") as handle:
        return [line for line in handle if not line.startswith("#")]


def _pixels(out_dir):
    clr = cooler.Cooler(str(out_dir / "abs_fragments_contacts_weighted.cool"))
    return clr.bins()[:], clr.pixels()[:]


def test_lanes_give_the_map_of_the_library(lanes, tmp_path):
    """The lanes of a library, given as the shell expands unquoted patterns,
    give the pairs and the map of the library in one pair of files."""
    skip_unless_can_align("bwa", GENOME)
    whole, split = tmp_path / "whole", tmp_path / "split"
    hpi.full_pipeline(
        genome=GENOME,
        input1=READS[0],
        input2=READS[1],
        aligner="bwa",
        enzyme="DpnII",
        out_dir=str(whole),
    )
    result = CliRunner().invoke(
        cli,
        [
            "pipeline",
            "-a",
            "bwa",
            "-e",
            "DpnII",
            "-g",
            GENOME,
            "-o",
            str(split),
            *_lane_files(lanes, "R1"),
            *_lane_files(lanes, "R2"),
        ],
    )
    assert result.exit_code == 0, result.output

    assert _pairs_rows(split) == _pairs_rows(whole)
    (bins, pixels), (whole_bins, whole_pixels) = _pixels(split), _pixels(whole)
    pd.testing.assert_frame_equal(bins, whole_bins)
    pd.testing.assert_frame_equal(pixels, whole_pixels)
    # The concatenated files are removed with the other temporary files
    assert not list((split / "tmp").glob("*.fq*"))

    (log,) = split.glob("hicstuff_*.log")
    log_text = log.read_text()
    r1, r2 = _lane_files(lanes, "R1"), _lane_files(lanes, "R2")
    assert f"## input1: {','.join(r1)}" in log_text
    assert f"## input2: {','.join(r2)}" in log_text
    assert f"{LANE1_READS} reads in each of {r1[0]} and {r2[0]}." in log_text
    assert f"{10000 - LANE1_READS} reads in each of {r1[1]} and {r2[1]}." in log_text
    assert "Concatenating 2 pairs of fastq files" in log_text
    assert hcs.get_pipeline_stats(str(log))["Total read pairs"] == 10000


@pytest.mark.parametrize("offset", [False, True])
def test_reads_are_counted_in_each_pair_of_files(lanes, tmp_path, caplog, monkeypatch, offset):
    """A pair of files with fewer reads in one file is named in an error, as
    the files of one pair are, and the run goes on. Pairs whose differences
    cancel out in the sum of the reads are too (offset)."""
    caplog.set_level(logging.INFO)
    r1, r2 = _lane_files(lanes, "R1"), _lane_files(lanes, "R2")
    for path in (r2[1], r1[0]) if offset else (r2[1],):
        lines = _fastq_lines(path)
        with gzip.open(path, "wt") as handle:
            handle.writelines(lines[:-4])

    def stop(files, out_file):
        raise RuntimeError("stop before the concatenation")

    monkeypatch.setattr(hpi, "concatenate_fastq", stop)
    with pytest.raises(RuntimeError, match="stop before the concatenation"):
        hpi.full_pipeline(
            genome=GENOME,
            input1=r1,
            input2=r2,
            aligner="bwa",
            enzyme="DpnII",
            out_dir=str(tmp_path / "out"),
        )
    lane2 = 10000 - LANE1_READS
    assert (
        f"Fastq files {r1[1]} and {r2[1]} do not have the same number of reads: "
        f"{lane2} and {lane2 - 1}." in caplog.text
    )
    assert "Fastq files do not have the same number of reads." in caplog.text
    assert "reads found in each fastq file" not in caplog.text


@pytest.mark.parametrize(
    "stage, input1, input2, message",
    [
        (
            "fastq",
            ["a_R1.fq", "b_R1.fq"],
            ["a_R2.fq"],
            "2 R1 files and 1 R2 files: each R1 file needs its R2 file.",
        ),
        (
            "bam",
            ["a_R1.bam", "b_R1.bam"],
            ["a_R2.bam", "b_R2.bam"],
            "Several pairs of input files are read from fastq only",
        ),
        ("pairs", ["a.pairs", "b.pairs"], None, "You must provide a single input file"),
        ("fastq", ["a_R1.fq"], None, "You must provide 2 input files"),
    ],
)
def test_full_pipeline_refuses_inputs(tmp_path, caplog, stage, input1, input2, message):
    """Before it writes anything."""
    out_dir = tmp_path / "out"
    with pytest.raises(SystemExit):
        hpi.full_pipeline(
            genome=GENOME,
            input1=input1,
            input2=input2,
            aligner="bwa",
            enzyme="DpnII",
            start_stage=stage,
            out_dir=str(out_dir),
        )
    assert message in caplog.text
    assert not out_dir.exists()


def test_full_pipeline_refuses_files_of_other_reads(lanes, tmp_path, caplog):
    out_dir = tmp_path / "out"
    with pytest.raises(SystemExit):
        hpi.full_pipeline(
            genome=GENOME,
            input1=_lane_files(lanes, "R1"),
            input2=_lane_files(lanes, "R2")[::-1],
            aligner="bwa",
            enzyme="DpnII",
            out_dir=str(out_dir),
        )
    assert "are not the R1 and R2 files of the same reads" in caplog.text
    assert not out_dir.exists()


@pytest.mark.parametrize(
    "inputs, message",
    [
        # Two R1 files and one R2 file, as the shell expands lib_*_R1.fq.gz lib_L1_R2.fq.gz
        (lambda lanes: [*_lane_files(lanes, "R1"), f"{lanes}/lib_L1_R2.fq.gz"], "3 input files"),
        (lambda lanes: [f"{lanes}/lib_*_R3.fq.gz", f"{lanes}/lib_*_R4.fq.gz"], "No file matches"),
    ],
)
def test_pipeline_command_refuses_inputs(lanes, tmp_path, inputs, message):
    """A usage error, before the run starts."""
    result = CliRunner().invoke(
        cli,
        ["pipeline", "-e", "DpnII", "-g", GENOME, "-o", str(tmp_path / "out"), *inputs(lanes)],
    )
    assert result.exit_code == 2
    assert message in result.output
    assert not (tmp_path / "out").exists()
