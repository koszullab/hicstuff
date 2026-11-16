"""Programs hicstuff runs: a run needs only those it uses.

Programs are hidden from PATH with the ``hide_tools`` fixture, so that each
test sees which ones a run looks for, runs, or can do without.
"""

import re
import shutil

import pytest
from click.testing import CliRunner

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


def _missing(tools, purpose):
    """Regex matching the whole message of require_tools for these tools."""
    verb = "is" if len(tools) == 1 else "are"
    msg = f"{' and '.join(tools)} {verb} not installed or not on PATH, and needed {purpose}."
    return f"^{re.escape(msg)}$"


# ---------------------------------------------------------------------------
# pipeline
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("aligner", ALIGNERS)
@pytest.mark.parametrize("stage", ["bam", "pairs", "pairs_idx"])
def test_start_stage_runs_without_aligner(
    stage, aligner, fastq_run, indexed_genome, run_pipeline, hide_tools, tmp_path
):
    """Runs that start from BAM or pairs files align nothing: they need
    neither the aligner nor its index programs."""
    ref = fastq_run()
    if stage == "bam":
        inputs = (str(ref.tmp("for.bam")), str(ref.tmp("rev.bam")))
    elif stage == "pairs":
        inputs = (str(ref.tmp("valid.pairs")), None)
    else:
        # The pairs_idx stage sorts its input in place: give it a copy.
        inputs = (str(tmp_path / "valid_idx.pairs"), None)
        shutil.copy(ref.tmp("valid_idx.pairs"), inputs[0])
    hide_tools(*ALIGNER_TOOLS, *([] if stage == "bam" else ["samtools"]))
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
    aligner, mapping, synthetic_lib, indexed_genome, hide_tools, tmp_path
):
    """Without its aligner, a run from fastq files stops before writing
    anything, naming the aligner and the start stage that needs it."""
    hide_tools(aligner)
    out_dir = tmp_path / "out"
    with pytest.raises(
        ImportError, match=_missing([aligner], "to align the reads of start stage fastq")
    ):
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
        ImportError, match=_missing(["bowtie2-build"], "to align the reads of start stage fastq")
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


@pytest.mark.parametrize("aligner", ALIGNERS)
def test_align_reads_needs_the_aligner(
    aligner, synthetic_lib, indexed_genome, hide_tools, tmp_path
):
    hide_tools(aligner)
    with pytest.raises(ImportError, match=_missing([aligner], "to align reads")):
        hpi.align_reads(
            str(synthetic_lib.reads1), indexed_genome, str(tmp_path / "out.bam"), aligner=aligner
        )


# ---------------------------------------------------------------------------
# iteralign, missview
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("aligner", ALIGNERS)
def test_iterative_align_needs_the_aligner(
    aligner, synthetic_lib, indexed_genome, hide_tools, tmp_path
):
    hide_tools(aligner)
    with pytest.raises(ImportError, match=_missing([aligner], "for iterative alignment")):
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
def test_missview_needs_the_aligner(aligner, hidden, synthetic_lib, hide_tools, tmp_path):
    """missview aligns reads it simulates, with the aligner, after building
    an index for a genome that has none: it checks for both first."""
    genome = tmp_path / "small.fa"
    genome.write_text(">chr3\n" + synthetic_lib.chroms["chr3"][:4000] + "\n")
    work = tmp_path / "work"
    hide_tools(hidden)
    args = ["missview", "-R", "60", "-a", aligner, "-T", str(work), str(genome)]
    result = CliRunner().invoke(cli, args + [str(tmp_path / "m.png")])
    assert isinstance(result.exception, ImportError)
    assert re.match(_missing([hidden], "to align the simulated reads"), str(result.exception))
    assert not (work / "simulated_reads.fq").exists()


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
