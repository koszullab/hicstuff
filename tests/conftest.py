"""Shared fixtures: headless plotting and a synthetic Hi-C library.

The synthetic library (see ``tests/synthetic.py``) is generated once per test
session. Pipeline runs on it are cached by name, so several test modules can
assert different things about the same run without paying for it twice.
"""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

import glob  # noqa: E402
import logging  # noqa: E402
import shutil  # noqa: E402
import subprocess as sp  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

import hicstuff.pipeline as hpi  # noqa: E402
from tests import synthetic as syn  # noqa: E402

ALIGNERS = ("bowtie2", "bwa", "minimap2")
MAPPINGS = ("normal", "iterative", "cutsite")
EXTERNAL_TOOLS = ("bowtie2", "bowtie2-build", "bwa", "minimap2", "samtools")


def detach_file_handlers():
    """Remove the log file handlers the pipeline attaches to the root logger,
    so that a run does not keep writing into the log of a previous one."""
    root = logging.getLogger()
    for handler in root.handlers[:]:
        if isinstance(handler, logging.FileHandler):
            root.removeHandler(handler)
            handler.close()


@pytest.fixture(autouse=True)
def _close_figures_and_logs():
    yield
    plt.close("all")
    detach_file_handlers()


@pytest.fixture(scope="session")
def synthetic_lib(tmp_path_factory):
    """The default synthetic genome, read pairs and truth table."""
    return syn.generate(tmp_path_factory.mktemp("synthetic"))


@pytest.fixture(scope="session")
def require_tools():
    missing = [tool for tool in EXTERNAL_TOOLS if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"external tools not installed: {', '.join(missing)}")


@pytest.fixture(scope="session")
def indexed_genome(synthetic_lib, require_tools, tmp_path_factory):
    """Copy of the synthetic genome with bowtie2 and bwa indices next to it.

    It lives in its own directory: hicstuff detects bwa indices by listing the
    whole directory, so indices must not sit next to other FASTA files.
    """
    genome = tmp_path_factory.mktemp("index") / "genome.fa"
    shutil.copy(synthetic_lib.genome, genome)
    sp.run(["bowtie2-build", "-q", str(genome), str(genome)], check=True)
    sp.run(["bwa", "index", str(genome)], check=True, stderr=sp.DEVNULL)
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
def run_pipeline(tmp_path_factory, require_tools):
    """Run the pipeline once per distinct ``name`` and cache the result."""
    cache = {}

    def _run(name, **kwargs):
        if name not in cache:
            out_dir = tmp_path_factory.mktemp(name)
            kwargs = dict(kwargs, out_dir=str(out_dir))
            try:
                hpi.full_pipeline(**kwargs)
            finally:
                detach_file_handlers()
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
