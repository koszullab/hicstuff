"""bwa in-process with pybwa: the backend that runs, the index it loads, and
alignments that are those of the bwa program."""

import os
import random
import shutil
import subprocess as sp
import sys

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

pybwa = pytest.importorskip("pybwa")

MAPPINGS = ("normal", "iterative", "cutsite")


def _content(path):
    """Lines of a pairs file, or pixels of a cool file."""
    if str(path).endswith(".cool"):
        return cooler.Cooler(str(path)).pixels()[:].values.tolist()
    with hio.read_compressed(str(path)) as handle:
        return handle.read().splitlines()


def _records(bam):
    with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as handle:
        return [record.to_string() for record in handle]


def _needs_bwa():
    if shutil.which("bwa") is None:
        pytest.skip("the bwa program is needed to compare with")


# ---------------------------------------------------------------------------
# backend
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("program", "module", "prefer_python", "expected"),
    [
        (True, True, False, "binary"),
        (False, True, False, "python"),
        (True, False, False, "binary"),
        (True, True, True, "python"),
        (True, False, True, "binary"),
    ],
)
def test_auto_backend(program, module, prefer_python, expected, hide_tools, monkeypatch):
    """auto runs the bwa program first, so that runs align as they did, and
    pybwa first for mapping mode 5spm."""
    _needs_bwa()
    if not program:
        hide_tools("bwa")
    if not module:
        monkeypatch.setattr(hbm, "pybwa_version", lambda: None)
    assert hbm.resolve_backend("auto", prefer_python=prefer_python) == expected


def test_backend_errors_name_the_fixes(hide_tools, monkeypatch):
    hide_tools("bwa")
    with pytest.raises(
        ImportError,
        match="^bwa is not installed or not on PATH, and needed to align reads with "
        "--aligner-backend binary: install it, or use --aligner-backend python",
    ):
        hbm.resolve_backend("binary")
    monkeypatch.setattr(hbm, "pybwa_version", lambda: None)
    with pytest.raises(
        ImportError,
        match=r"^pybwa is not installed, and needed to align reads with --aligner-backend "
        r"python: install it \(pip install pybwa",
    ):
        hbm.resolve_backend("python")
    with pytest.raises(ImportError, match="^Neither the bwa program nor pybwa is installed"):
        hbm.resolve_backend("auto")
    with pytest.raises(ValueError, match="Unknown aligner backend 'gpu'"):
        hbm.resolve_backend("gpu")


def test_backend_from_the_environment(monkeypatch):
    """HICSTUFF_ALIGNER_BACKEND chooses the backend when none is given."""
    _needs_bwa()
    monkeypatch.setenv("HICSTUFF_ALIGNER_BACKEND", "Python")
    assert hbm.resolve_backend() == "python"
    assert hbm.resolve_backend(None, prefer_python=False) == "python"
    assert hbm.resolve_backend("binary") == "binary"
    monkeypatch.setenv("HICSTUFF_ALIGNER_BACKEND", "")
    assert hbm.resolve_backend() == "binary"


def test_pybwa_not_installed(monkeypatch):
    monkeypatch.setitem(sys.modules, "pybwa", None)
    assert hbm.pybwa_version() is None


# ---------------------------------------------------------------------------
# index
# ---------------------------------------------------------------------------


def test_load_index_writes_nothing_next_to_the_genome(indexed_genome, tmp_path):
    """The index of the bwa program is linked into tmp_dir with a sequence
    dictionary of names and lengths, as the bwa program's SAM header."""
    genome_dir = os.path.dirname(indexed_genome)
    before = sorted(os.listdir(genome_dir))
    index = hbm.load_index(indexed_genome, str(tmp_path))
    assert sorted(os.listdir(genome_dir)) == before
    with pysam.FastxFile(indexed_genome) as fasta:
        sequences = [{"SN": record.name, "LN": len(record.sequence)} for record in fasta]
    assert index.header.to_dict()["SQ"] == sequences
    (index_dir,) = tmp_path.iterdir()
    links = sorted(path.name for path in index_dir.iterdir() if path.is_symlink())
    assert links == [f"index{suffix}" for suffix in hbm.INDEX_SUFFIXES]
    assert (index_dir / "index.dict").read_text().startswith("@SQ\tSN:")


def test_load_index_builds_bwa_index(synthetic_lib, tmp_path):
    """A reference without index is indexed into tmp_dir, with the files that
    bwa index writes."""
    _needs_bwa()
    genome = tmp_path / "genome.fa"
    shutil.copy(synthetic_lib.genome, genome)
    work = tmp_path / "work"
    work.mkdir()
    hbm.load_index(str(genome), str(work))
    assert sorted(os.listdir(tmp_path)) == ["genome.fa", "work"]
    expected = tmp_path / "expected"
    expected.mkdir()
    sp.run(
        ["bwa", "index", "-p", str(expected / "genome"), str(genome)], check=True, stderr=sp.DEVNULL
    )
    (index_dir,) = work.iterdir()
    for suffix in hbm.INDEX_SUFFIXES:
        built = (index_dir / f"index{suffix}").read_bytes()
        assert built == (expected / f"genome{suffix}").read_bytes(), suffix


@pytest.mark.parametrize(("min_length", "method"), [(100, "RB2"), (1000, "BWTSW")])
def test_build_index_runs_the_algorithm_of_bwa_index(min_length, method, tmp_path, monkeypatch):
    """bwtsw above 50 Mb, is below: pybwa 2.3 runs them for RB2 and BWTSW."""
    methods = []

    class Index:
        @classmethod
        def index(cls, fasta, method, prefix):
            methods.append(method)

    monkeypatch.setattr(pybwa, "BwaIndex", Index)
    monkeypatch.setattr(hbm, "BWTSW_MIN_LENGTH", min_length)
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">a\n" + "ACGT" * 25 + "\n>b\n" + "ACGT" * 15 + "\n")
    hbm.build_index(str(fasta), str(tmp_path / "index"))
    assert methods == [getattr(pybwa.BwaIndexBuildMethod, method)]


def test_pybwa_rb2_runs_bwtsw(synthetic_lib, tmp_path, capfd):
    """pybwa 2.3 passes its index methods to bwa off by one, which
    build_index relies on: its RB2 runs bwtsw. A pybwa that fixes it must
    come with a change of build_index."""
    pybwa.set_bwa_verbosity(pybwa.BwaVerbosity.INFO)
    try:
        method = pybwa.BwaIndexBuildMethod.RB2
        pybwa.BwaIndex.index(synthetic_lib.genome, method=method, prefix=str(tmp_path / "index"))
    finally:
        pybwa.set_bwa_verbosity(pybwa.BwaVerbosity.ERROR)
    assert "[BWTIncCreate]" in capfd.readouterr().err


# ---------------------------------------------------------------------------
# alignments
# ---------------------------------------------------------------------------


def test_read_names_as_bwa_writes_them(indexed_genome, synthetic_lib, tmp_path):
    """Names lose a final "/" and digit, as bwa trims them, when longer than
    two characters."""
    _needs_bwa()
    names = ["r1/1", "r2/2", "r3/7", "r4/x", "/1", "r5/12", "r6"]
    with pysam.FastxFile(str(synthetic_lib.reads1)) as reads:
        records = [next(reads) for _ in names]
    fastq = tmp_path / "reads.fq"
    fastq.write_text(
        "".join(f"@{n}\n{r.sequence}\n+\n{r.quality}\n" for n, r in zip(names, records))
    )
    sam = sp.run(
        ["bwa", "mem", "-v", "1", indexed_genome, str(fastq)],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    expected = [line.split("\t")[0] for line in sam.splitlines() if not line.startswith("@")]
    out = tmp_path / "out.bam"
    hbm.align(str(fastq), hbm.load_index(indexed_genome, str(tmp_path)), str(out))
    with pysam.AlignmentFile(str(out)) as bam:
        assert sorted(record.query_name for record in bam) == sorted(expected)
    assert sorted(expected) == ["/1", "r1", "r2", "r3", "r4/x", "r5/12", "r6"]


def test_five_spm_makes_the_5_end_primary(tmp_path):
    """Of a read cut by a ligation junction, the primary alignment is the
    longest part, and with -5SPM the 5' one."""
    rng = random.Random(5)
    chroms = {name: "".join(rng.choices("ACGT", k=5000)) for name in ("chr1", "chr2")}
    genome = tmp_path / "genome.fa"
    genome.write_text("".join(f">{name}\n{seq}\n" for name, seq in chroms.items()))
    read = chroms["chr1"][1000:1050] + chroms["chr2"][2000:2070]
    fastq = tmp_path / "read.fq"
    fastq.write_text(f"@chimeric\n{read}\n+\n{'I' * len(read)}\n")
    index = hbm.load_index(str(genome), str(tmp_path))
    primaries = {}
    for five_spm in (False, True):
        out = tmp_path / f"{five_spm}.bam"
        hbm.align(str(fastq), index, str(out), five_spm=five_spm)
        with pysam.AlignmentFile(str(out)) as bam:
            records = list(bam)
        primaries[five_spm] = [r.reference_name for r in records if not r.is_secondary]
        # Supplementary records are dropped, and with -5SPM the secondary ones
        assert not any(r.is_supplementary for r in records)
        assert all(not r.is_secondary for r in records) or not five_spm
    assert primaries == {False: ["chr2"], True: ["chr1"]}


@pytest.mark.parametrize("iterative", [False, True])
def test_align_reads_loads_the_index(
    iterative, indexed_genome, synthetic_lib, tmp_path, monkeypatch
):
    """align_reads given the reference with pybwa loads its index, writes the
    alignments of the bwa program, and leaves nothing in its temporary
    directory, the working directory by default."""
    _needs_bwa()
    monkeypatch.chdir(tmp_path)
    for backend in ("binary", "python"):
        hpi.align_reads(
            str(synthetic_lib.reads1),
            indexed_genome,
            f"{backend}.bam",
            aligner="bwa",
            iterative=iterative,
            backend=backend,
        )
    assert sorted(os.listdir(tmp_path)) == ["binary.bam", "python.bam"]
    assert _records(tmp_path / "python.bam") == _records(tmp_path / "binary.bam")


def test_iteralign_and_align_pairs_remove_the_index(indexed_genome, synthetic_lib, tmp_path):
    """The index that pybwa loads for iterative_align and align_pairs is
    removed from their temporary directory."""
    work = tmp_path / "work"
    work.mkdir()
    reads = str(synthetic_lib.reads1), str(synthetic_lib.reads2)
    bam1, bam2 = str(tmp_path / "for.bam"), str(tmp_path / "rev.bam")
    hci.iterative_align(
        reads[0], str(work), indexed_genome, 1, bam1, aligner="bwa", backend="python"
    )
    hbm.align_pairs(*reads, indexed_genome, bam1, bam2, backend="python", tmp_dir=str(work))
    assert list(work.iterdir()) == []


@pytest.mark.parametrize("iterative", [False, True])
def test_align_in_batches(iterative, indexed_genome, synthetic_lib, tmp_path):
    """Reads aligned in batches have one primary record each, and those of
    MAPQ 30 or more are where they are aligned in one batch: only ties
    between equal-score alignments, which have none, can differ."""
    index = hbm.load_index(indexed_genome, str(tmp_path))
    align = hbm.iterative_align if iterative else hbm.align
    alignments = []
    for batch_size in (hbm.BATCH_SIZE, 1000):
        out = tmp_path / f"{batch_size}.bam"
        align(str(synthetic_lib.reads1), index, str(out), batch_size=batch_size)
        with pysam.AlignmentFile(str(out)) as bam:
            records = list(bam)
        assert len(records) == len({r.query_name for r in records}) == 6000
        alignments.append(
            [
                (r.query_name, r.reference_name, r.reference_start, r.cigarstring)
                for r in records
                if r.mapping_quality >= 30
            ]
        )
    assert alignments[0] == alignments[1]


def test_iterative_align_reads_shorter_than_min_len(
    indexed_genome, synthetic_lib, tmp_path, caplog
):
    """Reads no longer than min_len are aligned at their length, then at full
    length, as iteralign aligns them."""
    _needs_bwa()
    bams = {}
    for backend in ("binary", "python"):
        bams[backend] = tmp_path / f"{backend}.bam"
        work = tmp_path / backend
        work.mkdir()
        hci.iterative_align(
            str(synthetic_lib.reads1),
            str(work),
            indexed_genome,
            1,
            str(bams[backend]),
            aligner="bwa",
            min_len=500,
            backend=backend,
        )
    assert _records(bams["python"]) == _records(bams["binary"])
    assert caplog.text.count("min_len is longer than the reads") == 2


def test_iterative_align_needs_a_bwa_index(synthetic_lib, tmp_path):
    """iteralign needs the reference indexed, also in-process."""
    genome = tmp_path / "genome.fa"
    shutil.copy(synthetic_lib.genome, genome)
    with pytest.raises(SystemExit):
        hci.iterative_align(
            str(synthetic_lib.reads1),
            str(tmp_path),
            str(genome),
            1,
            str(tmp_path / "out.bam"),
            aligner="bwa",
            backend="python",
        )


@pytest.mark.parametrize("mapping", MAPPINGS)
def test_fastq_run_with_pybwa(mapping, fastq_run, run_pipeline, hide_tools):
    """Without the bwa program, -a bwa aligns in-process with pybwa, and writes
    the alignments, pairs and matrix of the bwa program."""
    _needs_bwa()
    ref = fastq_run("bwa", mapping)
    hide_tools("bwa")
    run = run_pipeline(f"bwa_{mapping}_pybwa", **ref.kwargs)
    assert f"## aligner: bwa mem in-process, pybwa {hbm.pybwa_version()}" in run.log_text
    assert "## aligner" not in ref.log_text
    for name in ("fragments_list.txt", "info_contigs.txt"):
        assert run.path(name).read_bytes() == ref.path(name).read_bytes()
    for name in ("valid_idx_pcrfree.pairs.gz", "abs_fragments_contacts_weighted.cool"):
        assert _content(run.path(name)) == _content(ref.path(name))
    for bam in ("for.bam", "rev.bam"):
        assert _records(run.tmp(bam)) == _records(ref.tmp(bam))
        with pysam.AlignmentFile(str(run.tmp(bam))) as handle:
            (bwa,) = [pg for pg in handle.header.to_dict()["PG"] if pg["ID"] == "bwa"]
        assert bwa["VN"] == f"{hbm.PYBWA_BWA_VERSION} (pybwa {hbm.pybwa_version()})"


def test_fastq_run_indexes_with_pybwa(synthetic_lib, fastq_run, run_pipeline, hide_tools, tmp_path):
    """Without the bwa program nor an index, pybwa indexes the reference."""
    ref = fastq_run("bwa", "normal")
    genome = tmp_path / "genome.fa"
    shutil.copy(synthetic_lib.genome, genome)
    hide_tools("bwa")
    run = run_pipeline("bwa_normal_pybwa_index", **dict(ref.kwargs, genome=str(genome)))
    assert "bwa index not found" in run.log_text
    assert sorted(os.listdir(tmp_path)) == ["genome.fa"]
    for name in ("valid_idx_pcrfree.pairs.gz", "abs_fragments_contacts_weighted.cool"):
        assert _content(run.path(name)) == _content(ref.path(name))


# ---------------------------------------------------------------------------
# mapping mode 5spm
# ---------------------------------------------------------------------------


def _five_spm_kwargs(synthetic_lib, indexed_genome):
    return dict(
        genome=indexed_genome,
        input1=str(synthetic_lib.reads1),
        input2=str(synthetic_lib.reads2),
        mapping="5spm",
        enzyme=syn.LIBRARY_ENZYME,
        filter_events=True,
        pcr_duplicates=True,
        no_cleanup=True,
    )


def test_5spm_backends_write_the_same_pairs(synthetic_lib, indexed_genome, run_pipeline):
    """Mapping mode 5spm aligns with bwa, in-process first: the bwa program,
    which aligns both mates in one call, writes the same pairs and matrix."""
    _needs_bwa()
    kwargs = _five_spm_kwargs(synthetic_lib, indexed_genome)
    binary = run_pipeline("5spm_binary", aligner_backend="binary", **kwargs)
    python = run_pipeline("5spm_python", **kwargs)
    assert "## aligner: bwa mem -5SPM, bwa program" in binary.log_text
    assert f"## aligner: bwa mem -5SPM in-process, pybwa {hbm.pybwa_version()}" in python.log_text
    for name in ("valid_idx_pcrfree.pairs.gz", "abs_fragments_contacts_weighted.cool"):
        assert _content(python.path(name)) == _content(binary.path(name))
    # The mates' alignments of MAPQ 30 or more are the same. Only the bwa
    # program's are paired: unmapped mates are placed next to their mate.
    for bam in ("for.bam", "rev.bam"):
        alignments = []
        for run in (binary, python):
            with pysam.AlignmentFile(str(run.tmp(bam))) as handle:
                records = list(handle)
            assert not any(r.is_secondary or r.is_supplementary for r in records)
            alignments.append(
                [
                    (r.query_name, r.reference_name, r.reference_start, r.is_reverse, r.cigarstring)
                    for r in records
                    if r.mapping_quality >= 30
                ]
            )
        assert alignments[0] == alignments[1]
        assert len(alignments[0]) > 4000


def test_5spm_only_aligns_with_bwa(synthetic_lib, indexed_genome, tmp_path):
    kwargs = _five_spm_kwargs(synthetic_lib, indexed_genome)
    out_dir = tmp_path / "out"
    with pytest.raises(ValueError, match="Mapping mode 5spm aligns with bwa mem -5SPM, not with"):
        hpi.full_pipeline(aligner="bowtie2", out_dir=str(out_dir), **kwargs)
    with pytest.raises(ValueError, match="pybwa runs bwa only"):
        hpi.full_pipeline(
            **dict(kwargs, mapping="normal"),
            aligner="minimap2",
            aligner_backend="python",
            out_dir=str(out_dir),
        )
    assert not out_dir.exists()


def test_5spm_mates_named_differently(synthetic_lib, indexed_genome, tmp_path):
    """bwa stops on mates named differently, after it has written alignments:
    the error is raised with its message, and nothing is left of them."""
    _needs_bwa()
    reads2 = tmp_path / "reads2.fq"
    with pysam.FastxFile(str(synthetic_lib.reads2)) as reads:
        records = list(reads)
    records[-1].name = "another_name"
    reads2.write_text("".join(f"{record}\n" for record in records))
    bams = [tmp_path / "for.bam", tmp_path / "rev.bam"]
    with pytest.raises(RuntimeError, match="paired reads have different names"):
        hbm.align_pairs(
            str(synthetic_lib.reads1),
            str(reads2),
            indexed_genome,
            *map(str, bams),
            backend="binary",
        )
    assert sorted(os.listdir(tmp_path)) == ["reads2.fq"]


def test_5spm_bwa_fails_before_writing(synthetic_lib, tmp_path):
    """bwa without index stops before it writes its SAM header: its error is
    raised, and nothing is written."""
    _needs_bwa()
    bams = [tmp_path / "for.bam", tmp_path / "rev.bam"]
    with pytest.raises(RuntimeError, match="fail to locate the index files"):
        hbm.align_pairs(
            str(synthetic_lib.reads1),
            str(synthetic_lib.reads2),
            str(tmp_path / "no_index"),
            *map(str, bams),
            backend="binary",
        )
    assert os.listdir(tmp_path) == []


def test_unknown_mapping_mode(synthetic_lib, indexed_genome, tmp_path, caplog):
    kwargs = dict(_five_spm_kwargs(synthetic_lib, indexed_genome), mapping="5spn", aligner="bwa")
    with pytest.raises(ValueError):
        hpi.full_pipeline(out_dir=str(tmp_path / "out"), **kwargs)
    assert "mapping must be either normal, iterative, cutsite or 5spm." in caplog.text


def test_5spm_command_line(synthetic_lib, indexed_genome, tmp_path):
    """-m 5spm is accepted in any case, and only aligns with bwa."""
    args = ["pipeline", "-m", "5SPM", "-a", "minimap2", "-g", indexed_genome]
    args += ["-o", str(tmp_path / "out"), str(synthetic_lib.reads1), str(synthetic_lib.reads2)]
    result = CliRunner().invoke(cli, args)
    assert "Mapping mode 5spm aligns with bwa mem -5SPM, not with minimap2" in str(result.exception)


# ---------------------------------------------------------------------------
# iteralign, missview
# ---------------------------------------------------------------------------


def test_iteralign_with_pybwa(synthetic_lib, indexed_genome, hide_tools, tmp_path):
    """iteralign aligns in-process the reads that the bwa program aligns."""
    _needs_bwa()
    work = tmp_path / "work"
    work.mkdir()
    bams = {}
    for backend in ("binary", "python"):
        bams[backend] = tmp_path / f"{backend}.bam"
        args = ["iteralign", "-a", "bwa", "--aligner-backend", backend, "-g", indexed_genome]
        args += ["-T", str(work), "-o", str(bams[backend]), str(synthetic_lib.reads1)]
        result = CliRunner().invoke(cli, args)
        assert result.exit_code == 0, result.output
        hide_tools("bwa")
    assert _records(bams["python"]) == _records(bams["binary"])
    assert list(work.iterdir()) == []


def test_missview_with_pybwa(synthetic_lib, hide_tools, tmp_path):
    """missview aligns its reads in-process without the bwa program, and
    writes the matrix that the bwa program gives."""
    _needs_bwa()
    genome = tmp_path / "small.fa"
    genome.write_text(">chr3\n" + synthetic_lib.chroms["chr3"][:4000] + "\n")
    matrices = []
    for backend in ("binary", "python"):
        work = tmp_path / backend
        args = ["missview", "-R", "60", "-a", "bwa", "--aligner-backend", backend]
        args += ["-T", str(work), str(genome), str(tmp_path / f"{backend}.png")]
        result = CliRunner().invoke(cli, args)
        assert result.exit_code == 0, result.output
        matrices.append((work / "abs_fragments_contacts_weighted.txt").read_text())
        hide_tools("bwa")
    assert matrices[0] == matrices[1]
