"""Programs hicstuff runs: a run needs only those it uses.

Programs are hidden from PATH with the ``hide_tools`` fixture, so that each
test sees which ones a run looks for, runs, or can do without.
"""

import contextlib
import gzip
import os
import random
import re
import shutil
import signal
import subprocess as sp
import sys
import time

import cooler
import pysam
import pytest
from click.testing import CliRunner

import hicstuff.bwamem as hbm
import hicstuff.io as hio
import hicstuff.iteralign as hci
import hicstuff.pipeline as hpi
from hicstuff.cli import cli
from tests import synthetic as syn

ALIGNERS = ("bowtie2", "bwa", "minimap2")
MAPPINGS = ("normal", "iterative", "cutsite")
#: Programs that align reads, or build or read an aligner's index.
ALIGNER_TOOLS = ("bowtie2", "bowtie2-build", "bowtie2-inspect", "bwa", "minimap2")


def _sorted_lines(path):
    with hio.read_compressed(str(path)) as handle:
        return sorted(line for line in handle if not line.startswith("#"))


def _missing(tools, purpose, alternative=None):
    """Regex matching the whole message of require_tools for these tools."""
    verb = "is" if len(tools) == 1 else "are"
    msg = f"{' and '.join(tools)} {verb} not installed or not on PATH, and needed {purpose}"
    if alternative is None:
        msg += "."
    else:
        msg += f": install {'it' if len(tools) == 1 else 'them'}, or {alternative}."
    return f"^{re.escape(msg)}$"


#: How a run can align reads when the aligner it was given is missing.
WITH_PYBWA = "run with -a bwa, which aligns in-process with pybwa"
#: The error of a run that aligns with bwa, with neither the bwa program nor pybwa.
NO_BWA = (
    "^Neither the bwa program nor pybwa is installed, and one of them is needed to align "
    "reads with bwa: install pybwa"
)


def _alternative():
    """The alternative that missing aligners' errors name: none without pybwa."""
    return WITH_PYBWA if hbm.pybwa_version() is not None else None


# ---------------------------------------------------------------------------
# pipeline
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("aligner", ALIGNERS)
@pytest.mark.parametrize("stage", ["bam", "pairs", "pairs_idx"])
def test_start_stage_runs_without_aligner(
    stage, aligner, fastq_run, indexed_genome, run_pipeline, hide_tools, tmp_path
):
    """Runs that start from BAM or pairs files align nothing: they need
    neither the aligner nor its index programs, and none needs samtools."""
    ref = fastq_run()
    if stage == "bam":
        inputs = (str(ref.tmp("for.bam")), str(ref.tmp("rev.bam")))
    elif stage == "pairs":
        inputs = (str(ref.tmp("valid.pairs")), None)
    else:
        # The pairs_idx stage sorts its input in place: give it a copy.
        inputs = (str(tmp_path / "valid_idx.pairs"), None)
        shutil.copy(ref.tmp("valid_idx.pairs"), inputs[0])
    hide_tools(*ALIGNER_TOOLS, "samtools")
    run = run_pipeline(
        f"{stage}_without_{aligner}",
        genome=indexed_genome,
        input1=inputs[0],
        input2=inputs[1],
        aligner=aligner,
        start_stage=stage,
        enzyme=syn.LIBRARY_ENZYME,
        mat_fmt="graal",
    )
    assert _sorted_lines(run.path("valid_idx.pairs.gz")) == _sorted_lines(
        ref.tmp("valid_idx.pairs")
    )


@pytest.mark.parametrize("mapping", MAPPINGS)
@pytest.mark.parametrize("aligner", ALIGNERS)
def test_fastq_start_needs_the_aligner(
    aligner, mapping, synthetic_lib, indexed_genome, hide_tools, tmp_path, monkeypatch
):
    """Without its aligner, a run from fastq files stops before writing
    anything, naming the aligner and the start stage that needs it, and
    running with bwa in-process if pybwa is installed. bwa needs either."""
    hide_tools(aligner)
    if aligner == "bwa":
        # bwa runs in-process when pybwa is installed
        monkeypatch.setattr(hbm, "pybwa_version", lambda: None)
        message = NO_BWA
    else:
        message = _missing([aligner], "to align the reads of start stage fastq", _alternative())
    out_dir = tmp_path / "out"
    with pytest.raises(ImportError, match=message):
        hpi.full_pipeline(
            genome=indexed_genome,
            input1=str(synthetic_lib.reads1),
            input2=str(synthetic_lib.reads2),
            aligner=aligner,
            mapping=mapping,
            enzyme=syn.LIBRARY_ENZYME,
            out_dir=str(out_dir),
        )
    assert not out_dir.exists()


def test_fastq_start_needs_bowtie2_build_without_index(synthetic_lib, hide_tools, tmp_path):
    """A genome without bowtie2 index is indexed first, with bowtie2-build."""
    genome = tmp_path / "genome.fa"
    shutil.copy(synthetic_lib.genome, genome)
    hide_tools("bowtie2-build", "bwa", "minimap2")
    out_dir = tmp_path / "out"
    with pytest.raises(
        ImportError,
        match=_missing(
            ["bowtie2-build"], "to align the reads of start stage fastq", _alternative()
        ),
    ):
        hpi.full_pipeline(
            genome=str(genome),
            input1=str(synthetic_lib.reads1),
            input2=str(synthetic_lib.reads2),
            enzyme=syn.LIBRARY_ENZYME,
            out_dir=str(out_dir),
        )
    assert not out_dir.exists()


def test_bowtie2_index_genome_needs_bowtie2_inspect(hide_tools, tmp_path):
    """The genome sequence is read from a bowtie2 index with bowtie2-inspect,
    whatever the start stage: the only aligner program such a run needs."""
    pairs = tmp_path / "valid.pairs"
    shutil.copy("test_data/valid.pairs", pairs)
    out_dir = tmp_path / "out"
    kwargs = dict(
        genome="test_data/genome/seq",
        input1=str(pairs),
        start_stage="pairs",
        enzyme="DpnII",
        mat_fmt="graal",
        out_dir=str(out_dir),
    )
    hide_tools(*ALIGNER_TOOLS, "samtools")
    purpose = "to extract the genome sequence from the bowtie2 index test_data/genome/seq"
    with pytest.raises(ImportError, match=_missing(["bowtie2-inspect"], purpose)):
        hpi.full_pipeline(**kwargs)
    assert not out_dir.exists()
    hide_tools("bowtie2", "bowtie2-build", "bwa", "minimap2", "samtools")
    hpi.full_pipeline(**kwargs)
    assert len(_sorted_lines(out_dir / "valid_idx.pairs.gz")) > 0


@pytest.mark.parametrize(
    "aligner, mapping",
    [("bowtie2", m) for m in MAPPINGS] + [("bwa", "normal"), ("minimap2", "normal")],
)
def test_fastq_run_without_samtools(aligner, mapping, fastq_run, run_pipeline, hide_tools):
    """Alignments are filtered and sorted with the samtools that pysam bundles:
    a run from fastq files runs no samtools program, and writes what a run
    with samtools on PATH writes."""
    ref = fastq_run(aligner, mapping)
    # The genome is indexed already
    hide_tools("samtools", "bowtie2-build")
    run = run_pipeline(f"{aligner}_{mapping}_without_samtools", **ref.kwargs)
    for name in ("fragments_list.txt", "info_contigs.txt"):
        assert run.path(name).read_bytes() == ref.path(name).read_bytes()
    for name in ("valid_idx_pcrfree.pairs.gz", "abs_fragments_contacts_weighted.cool"):
        assert _content(run.path(name)) == _content(ref.path(name))
    for bam in ("for.bam", "rev.bam"):
        with pysam.AlignmentFile(str(run.tmp(bam)), "rb") as handle:
            programs = [pg for pg in handle.header.to_dict()["PG"] if pg["PN"] == "samtools"]
            reads = list(handle)
        assert programs
        assert all(pysam.__samtools_version__ in pg["VN"] for pg in programs)
        assert not any(read.is_supplementary for read in reads)
        with pysam.AlignmentFile(str(ref.tmp(bam)), "rb") as handle:
            assert [read.to_string() for read in reads] == [read.to_string() for read in handle]


def _content(path):
    """Lines of a pairs file, or pixels of a cool file."""
    if str(path).endswith(".cool"):
        return cooler.Cooler(str(path)).pixels()[:].values.tolist()
    with hio.read_compressed(str(path)) as handle:
        return handle.read().splitlines()


@pytest.mark.parametrize("aligner", ALIGNERS)
def test_align_reads_needs_the_aligner(
    aligner, synthetic_lib, indexed_genome, hide_tools, no_pybwa, tmp_path
):
    hide_tools(aligner)
    message = NO_BWA if aligner == "bwa" else _missing([aligner], "to align reads")
    with pytest.raises(ImportError, match=message):
        hpi.align_reads(
            str(synthetic_lib.reads1), indexed_genome, str(tmp_path / "out.bam"), aligner=aligner
        )


# ---------------------------------------------------------------------------
# iteralign, missview
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("aligner", ALIGNERS)
def test_iterative_align_needs_the_aligner(
    aligner, synthetic_lib, indexed_genome, hide_tools, no_pybwa, tmp_path
):
    hide_tools(aligner)
    message = NO_BWA if aligner == "bwa" else _missing([aligner], "for iterative alignment")
    with pytest.raises(ImportError, match=message):
        hci.iterative_align(
            str(synthetic_lib.reads1),
            str(tmp_path),
            indexed_genome,
            1,
            str(tmp_path / "out.bam"),
            aligner=aligner,
        )


def test_iteralign_command_needs_the_aligner(synthetic_lib, hide_tools, tmp_path):
    hide_tools("minimap2")
    work = tmp_path / "work"
    work.mkdir()
    args = ["iteralign", "-a", "minimap2", "-g", str(synthetic_lib.genome), "-T", str(work)]
    result = CliRunner().invoke(
        cli, args + ["-o", str(tmp_path / "r1.bam"), str(synthetic_lib.reads1)]
    )
    assert isinstance(result.exception, ImportError)
    assert re.match(_missing(["minimap2"], "for iterative alignment"), str(result.exception))
    # The temporary directory is removed
    assert list(work.iterdir()) == []


@pytest.mark.parametrize(
    "aligner, hidden",
    [
        ("bowtie2", "bowtie2"),
        ("bowtie2", "bowtie2-build"),
        ("bwa", "bwa"),
        ("minimap2", "minimap2"),
    ],
)
def test_missview_needs_the_aligner(aligner, hidden, synthetic_lib, hide_tools, no_pybwa, tmp_path):
    """missview aligns reads it simulates, with the aligner, after building
    an index for a genome that has none: it checks for both first. bwa needs
    the bwa program or pybwa."""
    genome = tmp_path / "small.fa"
    genome.write_text(">chr3\n" + synthetic_lib.chroms["chr3"][:4000] + "\n")
    work = tmp_path / "work"
    hide_tools(hidden)
    args = ["missview", "-R", "60", "-a", aligner, "-T", str(work), str(genome)]
    result = CliRunner().invoke(cli, args + [str(tmp_path / "m.png")])
    assert isinstance(result.exception, ImportError)
    message = NO_BWA if aligner == "bwa" else _missing([hidden], "to align the simulated reads")
    assert re.match(message, str(result.exception))
    assert not (work / "simulated_reads.fq").exists()


def test_missview_without_samtools(synthetic_lib, hide_tools, tmp_path, monkeypatch):
    genome = tmp_path / "small.fa"
    genome.write_text(">chr3\n" + synthetic_lib.chroms["chr3"][:4000] + "\n")
    monkeypatch.setattr("hicstuff.view.plot_matrix", lambda *args, **kwargs: None)
    hide_tools("samtools")
    args = ["missview", "-R", "100", "-a", "minimap2", "-b", "1000", str(genome)]
    result = CliRunner().invoke(cli, args + [str(tmp_path / "m.png")])
    assert result.exit_code == 0, (result.output, result.exception)


# ---------------------------------------------------------------------------
# samtools_pipe: samtools commands of pysam, piped
# ---------------------------------------------------------------------------


def _write_sam(path, names):
    """SAM file with, for each name in that order, a primary and a
    supplementary alignment."""
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 10000}]}
    with pysam.AlignmentFile(str(path), "w", header=header) as sam:
        for i, name in enumerate(names):
            for flag in (0, 2048):
                read = pysam.AlignedSegment(sam.header)
                read.query_name = name
                read.flag = flag
                read.reference_id = 0
                read.reference_start = i
                read.mapping_quality = 40
                read.cigarstring = "20M"
                read.query_sequence = "A" * 20
                sam.write(read)


@pytest.fixture
def unsorted_sam(tmp_path):
    """SAM file of reads whose names are in no order, natural or not
    (r9 < r10 < r100), each with a supplementary alignment."""
    names = [f"r{i}" for i in range(1000)] + [f"r{i}.{j}" for i in range(5) for j in (2, 10)]
    random.Random(0).shuffle(names)
    _write_sam(tmp_path / "unsorted.sam", names)
    return tmp_path / "unsorted.sam", sorted(names, key=hpi._natural_key)


def _names(bam):
    with pysam.AlignmentFile(str(bam), "rb") as handle:
        return [(read.query_name, read.flag) for read in handle]


def test_samtools_pipe_filters_and_sorts(unsorted_sam, tmp_path, hide_tools):
    """``samtools view -F 2048 | samtools sort -n``, as the pipeline runs it
    after aligning reads, without samtools on PATH."""
    sam, names = unsorted_sam
    hide_tools("samtools")
    out = tmp_path / "sorted.bam"
    hio.samtools_pipe(
        ["view", "-F", "2048", "-h", "-@", 2, sam], ["sort", "-n", "-@", 2, "-o", out, "-"]
    )
    assert _names(out) == [(name, 0) for name in names]


def test_samtools_pipe_reads_another_process(unsorted_sam, tmp_path, hide_tools):
    """``aligner | samtools sort -n``, as iterative alignment runs it."""
    sam, names = unsorted_sam
    hide_tools("samtools")
    out = tmp_path / "sorted.bam"
    writer = sp.Popen(["cat", str(sam)], stdout=sp.PIPE)
    hio.samtools_pipe(["sort", "-n", "-o", out], stdin=writer.stdout)
    assert writer.wait() == 0
    assert writer.stdout.closed
    assert [name for name, _ in _names(out)] == [name for name in names for _ in range(2)]


def test_samtools_pipe_failure(tmp_path, capfd):
    """A failing command raises, naming every command that failed, after
    samtools' own error message."""
    missing = tmp_path / "missing.sam"
    with pytest.raises(
        RuntimeError,
        match=re.escape("samtools view (exit code 1), samtools sort (exit code 1) failed"),
    ):
        hio.samtools_pipe(["view", "-h", missing], ["sort", "-o", tmp_path / "x.bam", "-"])
    assert f'failed to open "{missing}"' in capfd.readouterr().err


def test_samtools_pipe_ignores_modules_in_working_directory(unsorted_sam, tmp_path, monkeypatch):
    """A file of the working directory named like a module cannot stand in
    for it in the child process."""
    sam, names = unsorted_sam
    work = tmp_path / "work"
    work.mkdir()
    (work / "pysam.py").write_text("raise ImportError('not pysam')\n")
    monkeypatch.chdir(work)
    hio.samtools_pipe(["view", "-F", "2048", "-h", sam], ["sort", "-n", "-o", "out.bam", "-"])
    assert _names(work / "out.bam") == [(name, 0) for name in names]


#: Runs samtools_pipe until interrupted, on a standard input that never ends:
#: the read end of a pipe, given as file descriptor.
_INTERRUPTED_PIPE = """\
import os
import sys

import hicstuff.io as hio

try:
    hio.samtools_pipe(
        ["view", "-h", "-"],
        ["sort", "-n", "-o", sys.argv[2], "-"],
        stdin=os.fdopen(int(sys.argv[1]), "rb"),
    )
except KeyboardInterrupt:
    sys.exit(42)
"""


def _running(pgid):
    """Processes of a process group that have not exited."""
    pids = []
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            with open(f"/proc/{entry}/stat") as handle:
                # After the command name: state, parent, process group...
                fields = handle.read().rsplit(")", 1)[1].split()
        except OSError:
            continue
        if int(fields[2]) == pgid and fields[0] != "Z":
            pids.append(int(entry))
    return pids


def _runs_samtools(pid):
    """Whether a process has loaded the samtools of pysam."""
    try:
        with open(f"/proc/{pid}/maps") as handle:
            return "libcsamtools" in handle.read()
    except OSError:
        return False


@pytest.mark.skipif(not os.path.isdir("/proc"), reason="reads processes from /proc")
@pytest.mark.parametrize("signalled", ["group", "hicstuff"])
def test_samtools_pipe_stops_on_ctrl_c(signalled, tmp_path):
    """Ctrl-C, which signals every process of the terminal's foreground group,
    stops the samtools commands as it stops samtools programs, and so does an
    interruption of hicstuff alone: no command is left running after it."""
    read, write = os.pipe()
    script = tmp_path / "pipe.py"
    script.write_text(_INTERRUPTED_PIPE)
    proc = sp.Popen(
        [sys.executable, str(script), str(read), str(tmp_path / "out.bam")],
        pass_fds=(read,),
        start_new_session=True,
    )
    os.close(read)
    try:
        deadline = time.monotonic() + 60
        while True:
            commands = [pid for pid in _running(proc.pid) if pid != proc.pid]
            if len(commands) == 2 and all(_runs_samtools(pid) for pid in commands):
                break
            assert time.monotonic() < deadline, "the samtools commands did not start"
            time.sleep(0.1)
        # Both commands now wait for input, in samtools
        time.sleep(0.5)
        if signalled == "group":
            os.killpg(proc.pid, signal.SIGINT)
        else:
            os.kill(proc.pid, signal.SIGINT)
        assert proc.wait(timeout=60) == 42
        assert _running(proc.pid) == []
    finally:
        # A command left running would wait for this write end to close.
        os.close(write)
        with contextlib.suppress(ProcessLookupError):
            os.killpg(proc.pid, signal.SIGKILL)
        proc.wait()


class _Interrupted(Exception):
    pass


def _interrupt(signum, frame):
    raise _Interrupted


def test_samtools_pipe_kills_commands_when_interrupted(tmp_path, monkeypatch):
    """samtools_pipe kills its commands, and waits for them, when waiting for
    them is interrupted."""
    started = []

    class Popen(sp.Popen):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            started.append(self)

    monkeypatch.setattr(hio.sp, "Popen", Popen)
    read, write = os.pipe()
    previous = signal.signal(signal.SIGALRM, _interrupt)
    signal.setitimer(signal.ITIMER_REAL, 2)
    try:
        with pytest.raises(_Interrupted):
            # The input never ends: the commands run until killed.
            hio.samtools_pipe(
                ["view", "-h", "-"],
                ["sort", "-n", "-o", tmp_path / "out.bam", "-"],
                stdin=os.fdopen(read, "rb"),
            )
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)
        os.close(write)
    assert [process.returncode for process in started] == [-signal.SIGKILL] * 2


def test_check_bam_entries_without_samtools(unsorted_sam, tmp_path, hide_tools):
    sam, names = unsorted_sam
    hide_tools("samtools")
    # Every record is counted, supplementary ones included
    assert hio.check_bam_entries(str(sam)) == 2 * len(names)


# ---------------------------------------------------------------------------
# Final pairs file: BGZF, with or without bgzip
# ---------------------------------------------------------------------------

#: The empty block that ends a BGZF file.
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def _bgzf_blocks(path):
    """Uncompressed sizes of the blocks of a BGZF file, checking that each is a
    gzip member whose extra field gives its size, as pairix and tabix need."""
    data = path.read_bytes()
    sizes = []
    pos = 0
    while pos < len(data):
        header = data[pos : pos + 18]
        # gzip magic, deflate, FEXTRA flag; a 6 byte extra field "BC" holding
        # the block size minus one
        assert header[:4] == b"\x1f\x8b\x08\x04", pos
        assert header[10:16] == b"\x06\x00BC\x02\x00", pos
        end = pos + int.from_bytes(header[16:18], "little") + 1
        assert end <= len(data)
        sizes.append(int.from_bytes(data[end - 4 : end], "little"))
        pos = end
    assert data.endswith(BGZF_EOF)
    return sizes


def test_sort_and_compress_pairs_bgzip_and_pysam_agree(tmp_path, hide_tools):
    """bgzip, when installed, and pysam otherwise, write BGZF files that hold
    the same bytes."""
    if shutil.which("bgzip") is None:
        pytest.skip("bgzip is not installed")
    pairs = "test_data/valid_idx.pairs"
    bgzip, pysam_bgzf = tmp_path / "bgzip.pairs.gz", tmp_path / "pysam.pairs.gz"
    hpi.sort_and_compress_pairs(pairs, str(bgzip), tmp_dir=str(tmp_path), threads=2)
    hide_tools("bgzip")
    hpi.sort_and_compress_pairs(pairs, str(pysam_bgzf), tmp_dir=str(tmp_path), threads=2)
    # Several blocks of data, then the empty one
    assert len(_bgzf_blocks(bgzip)) > 2
    assert len(_bgzf_blocks(pysam_bgzf)) > 2
    assert gzip.decompress(pysam_bgzf.read_bytes()) == gzip.decompress(bgzip.read_bytes())


def test_pipeline_without_external_tools(final_pairs, synthetic_lib, hide_tools, tmp_path):
    """From pairs files, the pipeline runs with none of the aligners, samtools
    or bgzip installed, and still writes BGZF pairs."""
    pairs = tmp_path / "valid.pairs"
    shutil.copy(final_pairs, pairs)
    hide_tools(*ALIGNER_TOOLS, "samtools", "bgzip")
    out_dir = tmp_path / "out"
    hpi.full_pipeline(
        genome=str(synthetic_lib.genome),
        input1=str(pairs),
        start_stage="pairs_idx",
        enzyme=syn.LIBRARY_ENZYME,
        out_dir=str(out_dir),
    )
    final = out_dir / "valid.pairs.gz"
    assert sum(_bgzf_blocks(final)) == len(gzip.decompress(final.read_bytes()))
    assert _sorted_lines(final) == _sorted_lines(final_pairs)
    assert (out_dir / "abs_fragments_contacts_weighted.cool").exists()


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "aligner, index, tools",
    [
        ("bowtie2", "ref", ["bowtie2"]),
        ("bowtie2", None, ["bowtie2", "bowtie2-build"]),
        ("bwa", None, ["bwa"]),
        ("minimap2", "ref.fa", ["minimap2"]),
    ],
)
def test_alignment_tools(aligner, index, tools):
    assert hio.alignment_tools(aligner, index) == tools


def test_require_tools_names_every_missing_tool(hide_tools):
    hide_tools("bowtie2", "bwa")
    with pytest.raises(ImportError, match=_missing(["bowtie2", "bwa"], "for this")):
        hio.require_tools(["bowtie2", "sort", "bwa"], "for this")
    hio.require_tools(["sort"], "for this")
