"""Shared fixtures: headless plotting and a synthetic Hi-C library.

The synthetic library (see ``tests/synthetic.py``) is generated once per test
session. Pipeline runs on it are cached by name, so several test modules can
assert different things about the same run without paying for it twice.
"""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

import glob  # noqa: E402
import shutil  # noqa: E402
import subprocess as sp  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

import hicstuff.bwamem as hbm  # noqa: E402
import hicstuff.io as hio  # noqa: E402
import hicstuff.log as hcl  # noqa: E402
import hicstuff.pipeline as hpi  # noqa: E402
from tests import synthetic as syn  # noqa: E402

ALIGNERS = ("bowtie2", "bwa", "minimap2")
MAPPINGS = ("normal", "iterative", "cutsite")
EXTERNAL_TOOLS = ("bowtie2", "bowtie2-build", "bwa", "minimap2")
#: The programs that each aligner runs, its index builder included.
ALIGNER_TOOLS = {
    "bowtie2": ("bowtie2", "bowtie2-build"),
    "bwa": ("bwa",),
    "minimap2": ("minimap2",),
}


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "needs(*programs): skip the test unless the programs are on PATH"
    )


@pytest.hookimpl(tryfirst=True)
def pytest_runtest_setup(item):
    """Skip a test marked ``needs`` before its fixtures are set up, unless
    the programs it names are on PATH."""
    for marker in item.iter_markers("needs"):
        missing = [tool for tool in marker.args if shutil.which(tool) is None]
        if missing:
            pytest.skip(f"not on PATH: {', '.join(missing)}")


def can_align(aligner, genome=None):
    """Whether reads can be aligned with an aligner on a genome: with the
    programs a run needs, its index builder only if the genome has no index,
    or for bwa in-process, with pybwa."""
    index = None if genome is None else hio.check_fasta_index(genome, mode=aligner)
    if all(shutil.which(tool) is not None for tool in hio.alignment_tools(aligner, index)):
        return True
    return aligner == "bwa" and hbm.pybwa_version() is not None


def skip_unless_can_align(aligner, genome=None):
    if not can_align(aligner, genome):
        pytest.skip(
            f"{aligner} cannot align reads: {', '.join(ALIGNER_TOOLS[aligner])} not on PATH"
        )


@pytest.fixture(autouse=True)
def _close_figures_and_logs():
    yield
    plt.close("all")
    # The log file of a pipeline run must not catch the logs of later tests
    hcl.remove_file_handler()


@pytest.fixture(scope="session")
def synthetic_lib(tmp_path_factory):
    """The default synthetic genome, read pairs and truth table."""
    return syn.generate(tmp_path_factory.mktemp("synthetic"))


@pytest.fixture(scope="session")
def require_tools():
    missing = [tool for tool in EXTERNAL_TOOLS if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"external tools not installed: {', '.join(missing)}")


@pytest.fixture
def hide_tools(tmp_path_factory, monkeypatch):
    """Hide programs from PATH for the rest of the test.

    ``hide_tools("samtools", "bwa")`` sets PATH to a directory of links to
    every other program found on the test's PATH, the first of each name as a
    lookup would find it, so that neither shutil.which nor a shell finds the
    hidden ones while everything else keeps running. Each call replaces the
    previous one.
    """
    path = os.environ.get("PATH", "")

    def _hide(*tools):
        shadow = tmp_path_factory.mktemp("path")
        for directory in path.split(os.pathsep):
            if not os.path.isdir(directory):
                continue
            for name in sorted(os.listdir(directory)):
                target = os.path.join(directory, name)
                link = shadow / name
                if name in tools or os.path.lexists(link):
                    continue
                if os.path.isfile(target) and os.access(target, os.X_OK):
                    link.symlink_to(target)
        monkeypatch.setenv("PATH", str(shadow))
        assert all(shutil.which(tool) is None for tool in tools)
        return shadow

    return _hide


@pytest.fixture
def no_pybwa(monkeypatch):
    """Make hicstuff find pybwa not installed, for the rest of the test."""
    import hicstuff.bwamem as hbm

    monkeypatch.setattr(hbm, "pybwa_version", lambda: None)


@pytest.fixture(scope="session")
def indexed_genome(synthetic_lib, tmp_path_factory):
    """Copy of the synthetic genome with bowtie2 and bwa indices next to it,
    named after it (genome.fa.1.bt2, genome.fa.bwt...), in its own directory.
    Each index is built if its program is on PATH, and the bwa index with
    pybwa otherwise.
    """
    genome = tmp_path_factory.mktemp("index") / "genome.fa"
    shutil.copy(synthetic_lib.genome, genome)
    if shutil.which("bowtie2-build") is not None:
        sp.run(["bowtie2-build", "-q", str(genome), str(genome)], check=True)
    if shutil.which("bwa") is not None:
        sp.run(["bwa", "index", str(genome)], check=True, stderr=sp.DEVNULL)
    elif hbm.pybwa_version() is not None:
        hbm.build_index(str(genome), str(genome))
    return str(genome)


class PipelineRun:
    """Outputs of one ``hicstuff.pipeline.full_pipeline`` call."""

    def __init__(self, out_dir, kwargs):
        self.out_dir = out_dir
        self.tmp_dir = out_dir / "tmp"
        self.kwargs = kwargs
        logs = glob.glob(str(out_dir / "*.log"))
        assert len(logs) == 1, logs
        self.log = logs[0]
        with open(self.log) as handle:
            self.log_text = handle.read()

    def path(self, name):
        prefix = self.kwargs.get("prefix")
        return self.out_dir / (f"{prefix}.{name}" if prefix else name)

    def tmp(self, name):
        prefix = self.kwargs.get("prefix")
        return self.tmp_dir / (f"{prefix}.{name}" if prefix else name)


@pytest.fixture(scope="session")
def run_pipeline(tmp_path_factory):
    """Run the pipeline once per distinct ``name`` and cache the result. A
    run from fastq files is skipped if its aligner cannot align reads."""
    cache = {}

    def _run(name, **kwargs):
        if kwargs.get("start_stage", "fastq") == "fastq":
            mapping = kwargs.get("mapping", "normal")
            aligner = kwargs.get("aligner") or ("bwa" if mapping == "5spm" else "bowtie2")
            skip_unless_can_align(aligner, kwargs.get("genome"))
        if name not in cache:
            out_dir = tmp_path_factory.mktemp(name)
            kwargs = dict(kwargs, out_dir=str(out_dir))
            try:
                hpi.full_pipeline(**kwargs)
            finally:
                hcl.remove_file_handler()
            cache[name] = PipelineRun(out_dir, kwargs)
        return cache[name]

    return _run


@pytest.fixture(scope="session")
def fastq_run(synthetic_lib, indexed_genome, run_pipeline):
    """Factory for full runs from the synthetic FASTQ files, with event
    filtering, PCR duplicate removal and the distance law enabled."""

    def _run(aligner="bowtie2", mapping="normal", **extra):
        name = f"{aligner}_{mapping}" + "".join(f"_{k}" for k in sorted(extra))
        kwargs = dict(
            genome=indexed_genome,
            input1=str(synthetic_lib.reads1),
            input2=str(synthetic_lib.reads2),
            aligner=aligner,
            mapping=mapping,
            enzyme=syn.LIBRARY_ENZYME,
            filter_events=True,
            pcr_duplicates=True,
            distance_law=True,
            no_cleanup=True,
        )
        kwargs.update(extra)
        return run_pipeline(name, **kwargs)

    return _run


@pytest.fixture(scope="session")
def final_pairs(fastq_run, tmp_path_factory):
    """Uncompressed copy of the final (filtered, deduplicated) pairs of the
    default bowtie2 run, for tests that start the pipeline from pairs."""
    import gzip

    run = fastq_run()
    dest = tmp_path_factory.mktemp("inputs") / "final.pairs"
    with gzip.open(run.path("valid_idx_pcrfree.pairs.gz"), "rt") as src:
        dest.write_text(src.read())
    return dest


@pytest.fixture(scope="session")
def matrices(final_pairs, synthetic_lib, run_pipeline, tmp_path_factory):
    """The final pairs of the default run turned into graal, bg2 and
    fragment-level cool matrices, keyed by format; each value is a dict with
    the paths of the matrix, fragments_list.txt and info_contigs.txt."""
    out = {}
    names = {
        "graal": "abs_fragments_contacts_weighted.txt",
        "bg2": "abs_fragments_contacts_weighted.bg2",
        "cool": "abs_fragments_contacts_weighted.cool",
    }
    for fmt, name in names.items():
        # The pairs_idx stage sorts its input in place: give each run a copy.
        src = tmp_path_factory.mktemp(f"in_{fmt}") / "input.pairs"
        src.write_text(final_pairs.read_text())
        run = run_pipeline(
            f"pairs_idx_{fmt}",
            genome=str(synthetic_lib.genome),
            input1=str(src),
            start_stage="pairs_idx",
            enzyme=syn.LIBRARY_ENZYME,
            mat_fmt=fmt,
        )
        out[fmt] = {
            "run": run,
            "mat": str(run.path(name)),
            "frags": str(run.path("fragments_list.txt")),
            "chroms": str(run.path("info_contigs.txt")),
        }
    return out
