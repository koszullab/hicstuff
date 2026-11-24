"""Command-line tests on the synthetic library (via click's CliRunner)."""

import gzip
import re
import shlex
import shutil

import cooler
import metacooler
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

import hicstuff.distance_law as hcdl
import hicstuff.io as hio
import hicstuff.iteralign as hci
import hicstuff.pipeline as hpi
import hicstuff.view as hcv
from hicstuff import __version__
from hicstuff.cli import cli
from tests import synthetic as syn

COMMANDS = [
    "convert",
    "cutsite",
    "digest",
    "distancelaw",
    "filter",
    "iteralign",
    "missview",
    "pipeline",
    "rebin",
    "scalogram",
    "stats",
    "subsample",
    "view",
]


@pytest.fixture
def runner():
    return CliRunner()


def invoke(runner, *args, exit_code=0):
    result = runner.invoke(cli, [str(arg) for arg in args])
    assert result.exit_code == exit_code, (result.output, result.exception)
    return result


def text(result):
    """CLI output without rich-click's styling, box drawing and line wrapping.

    rich-click forces a terminal, hence ANSI styling, when GITHUB_ACTIONS or
    FORCE_COLOR is set, so the escape codes are stripped as well.
    """
    plain = re.sub(r"\x1b\[[0-9;]*m", "", result.output)
    return " ".join(re.sub(r"[\u2500-\u257f]", " ", plain).split())


def _pixels(path, frags=None):
    """Pixel counts of a matrix in any format, indexed by (bin1, bin2)."""
    if path.endswith(".cool"):
        pixels = cooler.Cooler(path).pixels()[:]
        return pixels.set_index(["bin1_id", "bin2_id"])["count"]
    mat, _, _ = hio.flexible_hic_loader(path, frags, quiet=True)
    return pd.Series(mat.data, index=pd.MultiIndex.from_arrays([mat.row, mat.col]))


def _dense(mat):
    """Full symmetric dense matrix from an upper-triangle sparse matrix."""
    dense = mat.toarray()
    return dense + dense.T - np.diag(np.diag(dense))


def _cool_sum(path):
    return cooler.Cooler(str(path)).pixels()[:]["count"].sum()


# ---------------------------------------------------------------------------
# Help and generic errors
# ---------------------------------------------------------------------------


def test_every_command_is_listed(runner):
    assert sorted(cli.commands) == COMMANDS
    result = invoke(runner, "--help")
    for command in COMMANDS:
        assert command in result.output
    assert __version__ in invoke(runner, "--version").output


@pytest.mark.parametrize("command", COMMANDS)
def test_command_help(runner, command):
    result = invoke(runner, command, "--help")
    assert "Usage" in result.output
    assert command in result.output


@pytest.mark.parametrize(
    "args",
    [
        ["view", "missing.bg2"],
        ["stats", "missing.log"],
        ["convert", "-T", "bg2", "missing.bg2", "out"],
        ["subsample", "missing.bg2", "out"],
        ["digest", "-e", "DpnII", "-o", "{tmp}/out", "missing.fa"],
        pytest.param(
            ["pipeline", "-g", "{genome}", "-o", "{tmp}/out", "missing_R1.fq", "missing_R2.fq"],
            marks=pytest.mark.needs("bowtie2", "bowtie2-build"),
        ),
    ],
    ids=lambda args: args[0],
)
def test_missing_input(runner, args, tmp_path, monkeypatch, synthetic_lib):
    monkeypatch.chdir(tmp_path)
    args = [a.format(tmp=tmp_path, genome=synthetic_lib.genome) for a in args]
    result = invoke(runner, *args, exit_code=1)
    assert isinstance(result.exception, (FileNotFoundError, OSError))


def test_unknown_enzyme(runner, synthetic_lib, final_pairs, tmp_path):
    result = invoke(
        runner,
        "digest",
        "-e",
        "NotAnEnzyme",
        "-o",
        tmp_path / "d",
        synthetic_lib.genome,
        exit_code=1,
    )
    assert "NotAnEnzyme is not a valid restriction enzyme" in str(result.exception)
    src = tmp_path / "input.pairs"
    shutil.copy(final_pairs, src)
    result = invoke(
        runner,
        *("pipeline", "-S", "pairs_idx", "-e", "NotAnEnzyme", "-o", tmp_path / "p"),
        *("-g", synthetic_lib.genome, src),
        exit_code=1,
    )
    assert "not a valid restriction enzyme" in str(result.exception)


# ---------------------------------------------------------------------------
# digest, cutsite, iteralign
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("enzymes", ["HindIII", "DpnII,HindIII"])
def test_digest(runner, synthetic_lib, tmp_path, enzymes):
    out, figs = tmp_path / "out", tmp_path / "figs" / "digest"
    invoke(runner, "digest", "-e", enzymes, "-o", out, synthetic_lib.genome)
    frags = pd.read_csv(out / "fragments_list.txt", sep="\t")
    expected = synthetic_lib.fragments(enzymes.split(","))
    pd.testing.assert_frame_equal(frags.drop(columns="gc_content"), expected)
    # The output directory exists: only overwritten with --force.
    result = invoke(runner, "digest", "-e", enzymes, "-o", out, synthetic_lib.genome, exit_code=1)
    assert "already exists" in str(result.exception)
    # The figure directory is created, as with filter -f
    invoke(runner, "digest", "-e", enzymes, "-o", out, "-F", "-p", "-f", figs, synthetic_lib.genome)
    assert (figs / "frags_hist.pdf").stat().st_size > 0


def test_digest_defaults_to_current_directory(runner, synthetic_lib, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    invoke(runner, "digest", "-e", "HindIII", synthetic_lib.genome, exit_code=1)
    invoke(runner, "digest", "-e", "HindIII", "-F", synthetic_lib.genome)
    assert (tmp_path / "info_contigs.txt").exists()


@pytest.mark.parametrize("mode, pairs_per_chimera", [("for_vs_rev", 2), ("all", 3), ("pile", 2)])
def test_cutsite(runner, synthetic_lib, tmp_path, mode, pairs_per_chimera):
    """Reads are cut at ligation junctions: only chimeric R1 reads contain one,
    giving two pieces that are paired with R2 according to the mode."""
    prefix = tmp_path / "sub" / "cut"
    invoke(
        runner,
        *("cutsite", "-F", synthetic_lib.reads1, "-R", synthetic_lib.reads2),
        *("-e", "DpnII", "-m", mode, "-p", prefix, "-t", "2"),
    )
    reads = {}
    for mate in ("R1", "R2"):
        with gzip.open(f"{prefix}_{mate}.fq.gz", "rt") as handle:
            reads[mate] = handle.read().splitlines()[::4]
    assert reads["R1"] == reads["R2"]
    truth = synthetic_lib.truth
    n_chimeric = int(truth.chimeric.sum())
    assert len(reads["R1"]) == len(truth) + (pairs_per_chimera - 1) * n_chimeric
    names = pd.Series(reads["R1"]).str.extract(r"^@(sim\d+):(\d+)$")
    cut = (
        set(names[0][names[1] != "01"]) if mode != "for_vs_rev" else set(names[0][names[1] != "00"])
    )
    assert cut == set(truth.read_id[truth.chimeric])


def test_cutsite_rejects_unknown_mode(runner, synthetic_lib, tmp_path):
    invoke(
        runner,
        *("cutsite", "-F", synthetic_lib.reads1, "-R", synthetic_lib.reads2),
        *("-e", "DpnII", "-m", "shuffle", "-p", tmp_path / "x"),
        exit_code=2,
    )


def test_iteralign(runner, synthetic_lib, tmp_path, monkeypatch, require_tools):
    monkeypatch.chdir(tmp_path)
    out = tmp_path / "r1.bam"
    invoke(
        runner,
        *("iteralign", "-a", "minimap2", "-g", synthetic_lib.genome, "-o", out),
        *("-l", "40", synthetic_lib.reads1),
    )
    # One record per read, mapped or not; the temporary directory is removed.
    assert hio.check_bam_entries(str(out)) == synthetic_lib.params.n_pairs
    assert [p.name for p in tmp_path.iterdir()] == ["r1.bam"]


def test_iterative_align_errors(synthetic_lib, tmp_path, caplog, require_tools):
    with pytest.raises(SystemExit):
        # No bowtie2 index next to this genome
        hci.iterative_align(
            str(synthetic_lib.reads1), str(tmp_path), str(synthetic_lib.genome), 1, "x.bam"
        )
    with pytest.raises(ValueError, match="Unknown aligner"):
        hci.iterative_align(
            str(synthetic_lib.reads1),
            str(tmp_path),
            str(synthetic_lib.genome),
            1,
            str(tmp_path / "x.bam"),
            aligner="blast",
            min_len=500,
        )
    assert "min_len is longer than the reads" in caplog.text


# ---------------------------------------------------------------------------
# pipeline, filter, stats
# ---------------------------------------------------------------------------


_WITHOUT_MATRIX = ["pipeline", "-S", "pairs", "-e", "DpnII", "-g", "test_data/genome/seq.fa"]


@pytest.mark.parametrize(
    "matrix, options, named",
    [
        ("none", ["-b", "1000"], "writes no matrix: --binning is an option"),
        ("none", ["-B", "--max-iters 10"], "writes no matrix: --balancing-args is an option"),
        ("none", ["-E", "seq2"], "writes no matrix: --exclude is an option"),
        (
            "none",
            ["-b", "1000", "-E", "seq2"],
            "writes no matrix: --binning, --exclude are options",
        ),
        ("metacool", ["-z", "true"], "is neither zoomified nor balanced: --zoomify is an option"),
        (
            "metacool",
            ["-b", "1000", "-z", "True", "-B", "--max-iters 10", "-E", "seq2"],
            "is neither zoomified nor balanced: --zoomify, --balancing-args are options",
        ),
    ],
)
def test_pipeline_refuses_cool_options(runner, matrix, options, named, tmp_path):
    """-M none writes no matrix, and -M metacool is binned but neither
    zoomified nor balanced: the options of the cool matrix, given with them,
    are a usage error, before the run starts. -M metacool takes --binning and
    --exclude."""
    out_dir = tmp_path / "out"
    args = [*_WITHOUT_MATRIX, "-M", matrix, *options, "-o", out_dir, "test_data/valid.pairs"]
    result = invoke(runner, *args, exit_code=2)
    assert f"--matrix {matrix} {named} of the cool matrix." in text(result)
    assert not out_dir.exists()


def _bins_file(tmp_path):
    """A BED file of bins of no single size: of 7 kb on seq1, 3 kb on seq2."""
    path = tmp_path / "bins.bed"
    rows = [("seq1", start, min(start + 7000, 60_000)) for start in range(0, 60_000, 7000)]
    rows += [("seq2", start, min(start + 3000, 20_000)) for start in range(0, 20_000, 3000)]
    path.write_text("".join(f"{chrom}\t{start}\t{end}\n" for chrom, start, end in rows))
    return path


@pytest.mark.parametrize("binning", ["-5", "nowhere.bed"])
def test_pipeline_refuses_a_binning(runner, binning, tmp_path):
    """--binning is 0, a bin size in bp or a BED file of bins: another value
    is a usage error, before the run starts."""
    out_dir = tmp_path / "out"
    args = [*_WITHOUT_MATRIX, "-b", binning, "-o", out_dir, "test_data/valid.pairs"]
    result = invoke(runner, *args, exit_code=2)
    assert "Invalid value for '-b' / '--binning'" in text(result)
    assert not out_dir.exists()


@pytest.mark.parametrize("options", [["-b", "5000"], ["-b", "5000", "-z", "false"], ["-b", "BED"]])
def test_pipeline_binned_metacool_matrix(runner, options, tmp_path):
    """-M metacool takes --binning, a bin size or a BED file of bins, and -z
    false: it writes one metacool matrix of those bins, neither zoomified nor
    balanced."""
    bed = _bins_file(tmp_path)
    options = [str(bed) if option == "BED" else option for option in options]
    out_dir = tmp_path / "out"
    invoke(
        runner, *_WITHOUT_MATRIX, "-M", "metacool", *options, "-o", out_dir, "test_data/valid.pairs"
    )
    matrices = sorted(path.name for path in out_dir.iterdir() if "cool" in path.suffix)
    assert matrices == ["abs_fragments_contacts_weighted.metacool"]
    clr = cooler.Cooler(str(out_dir / matrices[0]))
    assert clr.binsize == (None if options[1] == str(bed) else 5000)
    assert "weight" not in clr.bins().columns
    assert clr.pixels()[:]["count"].sum() == 10_000


def test_pipeline_bed_binned_cool_matrix(runner, tmp_path):
    """-M cool with a BED file of bins writes one balanced cool matrix: -z
    left to its default is ignored, and -z true is a usage error, before the
    run starts, with -M cool and without -M."""
    bed = _bins_file(tmp_path)
    refused = tmp_path / "refused"
    for matrix in (["-M", "cool"], []):
        args = [*_WITHOUT_MATRIX, *matrix, "-b", bed, "-z", "true", "-o", refused]
        result = invoke(runner, *args, "test_data/valid.pairs", exit_code=2)
        assert f"--zoomify is refused with the bins of {bed}" in text(result)
        assert not refused.exists()
    out_dir = tmp_path / "out"
    invoke(
        runner, *_WITHOUT_MATRIX, "-M", "cool", "-b", bed, "-o", out_dir, "test_data/valid.pairs"
    )
    matrices = sorted(path.name for path in out_dir.iterdir() if "cool" in path.suffix)
    assert matrices == ["abs_fragments_contacts_weighted.cool"]
    assert "weight" in cooler.Cooler(str(out_dir / matrices[0])).bins().columns


@pytest.mark.parametrize(
    "options, named",
    [
        (["-B", "--max-iters 10"], "--balancing-args is an option"),
        (["-z", "true"], "--zoomify is an option"),
        (["-z", "true", "-B", "--max-iters 10"], "--zoomify, --balancing-args are options"),
    ],
)
def test_pipeline_refuses_cool_options_with_the_default_metacool(
    runner, options, named, tmp_path, monkeypatch
):
    """From MANY_CONTIGS contigs, the default matrix is metacool, binned or
    not: -z true and --balancing-args, given with it, are a usage error, once
    the genome is read, that says to give --matrix cool."""
    monkeypatch.setattr(hpi, "MANY_CONTIGS", 2)
    out_dir = tmp_path / "out"
    args = [*_WITHOUT_MATRIX, "-b", "5000", *options, "-o", out_dir, "test_data/valid.pairs"]
    result = invoke(runner, *args, exit_code=2)
    assert (
        "2 contigs in the matrix: its default format is metacool, which is neither zoomified "
        f"nor balanced: {named} of the cool matrix. Give --matrix cool to write the cool matrix"
    ) in text(result)
    assert not [path for path in out_dir.iterdir() if path.is_file() and path.suffix != ".log"]


def test_pipeline_default_binned_metacool_matrix(runner, tmp_path, monkeypatch):
    """From MANY_CONTIGS contigs, --binning makes the default matrix a binned
    metacool matrix: -z left to its default, or false, is taken."""
    monkeypatch.setattr(hpi, "MANY_CONTIGS", 2)
    for name, options in [("default", []), ("no-zoomify", ["-z", "false"])]:
        out_dir = tmp_path / name
        invoke(
            runner, *_WITHOUT_MATRIX, "-b", "5000", *options, "-o", out_dir, "test_data/valid.pairs"
        )
        matrices = sorted(path.name for path in out_dir.iterdir() if "cool" in path.suffix)
        assert matrices == ["abs_fragments_contacts_weighted.metacool"]


@pytest.mark.parametrize(
    "options", [[], ["-b", "0", "-z", "False", "-B", "", "-E", ""], ["-B", " ", "-E", " , "]]
)
def test_pipeline_without_matrix(runner, options, tmp_path):
    """-M none ends the run with its pairs file. The options of the cool
    matrix are taken with values that would shape no matrix: no binning, no
    zoomify, blank balancing arguments, no contig to exclude."""
    out_dir = tmp_path / "out"
    invoke(runner, *_WITHOUT_MATRIX, "-M", "none", *options, "-o", out_dir, "test_data/valid.pairs")
    written = sorted(
        path.name for path in out_dir.iterdir() if path.is_file() and path.suffix != ".log"
    )
    assert written == ["fragments_list.txt", "info_contigs.txt", "valid_idx.pairs.gz"]


@pytest.mark.parametrize(
    "options",
    [["--matrix", "none"], ["--matfmt", "none"], ["--matfmt", "none", "-M", "none"]],
)
def test_pipeline_matrix_option(runner, options, tmp_path):
    """--matrix (-M) chooses the matrix format. --matfmt, its former name,
    still works, alone or with --matrix of the same format, but the help no
    longer shows it."""
    out_dir = tmp_path / "out"
    invoke(runner, *_WITHOUT_MATRIX, *options, "-o", out_dir, "test_data/valid.pairs")
    written = sorted(
        path.name for path in out_dir.iterdir() if path.is_file() and path.suffix != ".log"
    )
    assert written == ["fragments_list.txt", "info_contigs.txt", "valid_idx.pairs.gz"]
    help_text = text(invoke(runner, "pipeline", "--help"))
    assert "--matrix -M STR Output matrix format:" in help_text
    assert "--matfmt" not in help_text
    assert "[default: (cool, or metacool from 50,000 contigs)]" in help_text


def test_pipeline_matrix_option_with_another_former_name(runner, tmp_path):
    """--matfmt and --matrix of two formats are a usage error, before the run
    starts."""
    out_dir = tmp_path / "out"
    args = [*_WITHOUT_MATRIX, "--matfmt", "bg2", "-M", "none", "-o", out_dir]
    result = invoke(runner, *args, "test_data/valid.pairs", exit_code=2)
    assert "--matfmt is the former name of --matrix: give --matrix none alone." in text(result)
    assert not out_dir.exists()


def test_pipeline_metacool_matrix(runner, tmp_path):
    """-M metacool writes a metacool matrix, without the contigs of --exclude,
    and stores the command line of the run, as its log records it."""
    out_dir = tmp_path / "out"
    args = [*_WITHOUT_MATRIX, "-M", "metacool", "-E", "seq2", "-o", out_dir]
    invoke(runner, *args, "test_data/valid.pairs")
    mc = metacooler.MetaCooler(str(out_dir / "abs_fragments_contacts_weighted.metacool"))
    assert mc.contigs.name.tolist() == ["seq1"]
    (log,) = out_dir.glob("hicstuff_*.log")
    (cmd,) = [line for line in log.read_text().splitlines() if "## cmd: " in line]
    assert mc.metadata["commands"] == [cmd.split("## cmd: ", 1)[1]]
    # CliRunner names the program after the function: cli
    assert mc.metadata["commands"][0].startswith("cli pipeline -S pairs -e DpnII")


def test_pipeline_command(runner, synthetic_lib, final_pairs, tmp_path, monkeypatch):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    monkeypatch.chdir(tmp_path)
    args = ["pipeline", "-S", "pairs_idx", "-e", "mnase", "-M", "graal", "-n"]
    args += ["-g", synthetic_lib.genome, inputs / "input.pairs"]
    shutil.copy(final_pairs, inputs / "input.pairs")
    invoke(runner, *args)
    mat = hio.load_sparse_matrix(str(tmp_path / "abs_fragments_contacts_weighted.txt"))
    assert mat.sum() == len(syn.read_pairs(final_pairs))
    # MNase libraries are binned in 100 bp chunks
    frags = pd.read_csv(tmp_path / "fragments_list.txt", sep="\t")
    assert (frags["size"].iloc[:10] == 100).all()
    # Existing outputs are only overwritten with --force
    shutil.copy(final_pairs, inputs / "input.pairs")
    result = invoke(runner, *args, exit_code=1)
    assert isinstance(result.exception, OSError)
    invoke(runner, *args, "-F")
    # Filtering 3C events requires restriction fragments
    result = invoke(runner, *args, "-F", "-f", "-e", "5000", exit_code=2)
    assert "Cannot filter events" in text(result)


def test_pipeline_command_line(runner, monkeypatch):
    """The pipeline gets the command that ran it, quoted to paste in a shell."""
    calls = []
    monkeypatch.setattr(hpi, "full_pipeline", lambda **kwargs: calls.append(kwargs))
    args = ["pipeline", "-g", "genomes/sac cer3.fa", "R1.fq.gz", "R2.fq.gz", "--mapping", "5spm"]
    invoke(runner, *args)
    assert calls[0]["cmd"] == f"cli {shlex.join(args)}"
    assert "'genomes/sac cer3.fa'" in calls[0]["cmd"]


def test_filter_command(runner, fastq_run, tmp_path):
    pairs = fastq_run().tmp("valid_idx.pairs")
    auto, manual = tmp_path / "auto.pairs", tmp_path / "manual.pairs"
    figs = tmp_path / "figs"
    invoke(runner, "filter", "-p", "-f", figs, "-P", "synthetic", pairs, auto)
    assert (figs / "event_distance.pdf").exists()
    assert (figs / "event_distribution.pdf").exists()
    thresholds = "-".join(re.search(r"uncuts=(\d+) loops=(\d+)", fastq_run().log_text).groups())
    invoke(runner, "filter", "-t", thresholds, pairs, manual)
    assert manual.read_text() == auto.read_text()
    assert len(syn.read_pairs(manual)) < len(syn.read_pairs(pairs))
    result = invoke(runner, "filter", "-t", "three", pairs, manual, exit_code=2)
    assert "Thresholds must be two integers" in text(result)


def test_stats_command(runner, fastq_run, synthetic_lib, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    result = invoke(runner, "stats", fastq_run("bowtie2", "iterative").log)
    assert f"Total read pairs:   {synthetic_lib.params.n_pairs:,}" in result.output
    assert "PCR duplicates:" in result.output
    result = invoke(runner, "stats", fastq_run("bowtie2", "normal").log)
    assert f"Total read pairs:   {synthetic_lib.params.n_pairs:,}" in result.output


def test_stats_command_without_mapping(runner, matrices):
    result = invoke(runner, "stats", matrices["graal"]["run"].log, exit_code=1)
    assert "No read mapping statistics found" in text(result)


# ---------------------------------------------------------------------------
# view
# ---------------------------------------------------------------------------


@pytest.fixture
def plotted(monkeypatch):
    """Record the arrays and settings given to hicstuff.view.plot_matrix."""
    calls = []
    original = hcv.plot_matrix

    def _record(array, **kwargs):
        calls.append((np.array(array), kwargs))
        return original(array, **kwargs)

    monkeypatch.setattr(hcv, "plot_matrix", _record)
    return calls


@pytest.fixture(scope="module")
def cool_dense(matrices):
    mat, _, _ = hio.load_cool(matrices["cool"]["mat"])
    return _dense(mat).astype(float)


def test_view_raw_counts(runner, matrices, cool_dense, plotted, tmp_path):
    out = tmp_path / "map.png"
    invoke(runner, "view", "-D", "30", "-o", out, matrices["cool"]["mat"])
    assert out.stat().st_size > 0
    array, kwargs = plotted[0]
    # Empty pixels are not drawn (set to inf)
    np.testing.assert_array_equal(np.where(np.isinf(array), 0, array), cool_dense)
    assert kwargs["vmax"] == np.percentile(cool_dense[cool_dense > 0], 99)
    assert kwargs["vmin"] == 0


@pytest.mark.parametrize(
    "transform, func",
    [
        ("log2", np.log2),
        ("log10", np.log10),
        ("ln", np.log),
        ("sqrt", np.sqrt),
        ("exp0.5", lambda x: x**0.5),
        ("log1p", np.log1p),
    ],
)
def test_view_transforms(runner, matrices, cool_dense, plotted, tmp_path, transform, func):
    invoke(
        runner, "view", "-D", "30", "-T", transform, "-M", "95%", "-m", "1",
        "-o", tmp_path / "map.png", matrices["cool"]["mat"],
    )  # fmt: skip
    array, kwargs = plotted[0]
    with np.errstate(divide="ignore"):
        expected = func(cool_dense)
    np.testing.assert_allclose(array, expected)
    assert kwargs["vmin"] == 1


def test_view_rejects_unknown_transform(runner, matrices, tmp_path):
    result = invoke(
        runner, "view", "-T", "fourier", "-o", tmp_path / "x.png", matrices["cool"]["mat"],
        exit_code=1,
    )  # fmt: skip
    assert isinstance(result.exception, TypeError)


@pytest.mark.parametrize("fmt", ["graal", "bg2", "cool"])
@pytest.mark.parametrize(
    "options",
    [
        ["-n", "-N", "2"],
        ["-t", "3", "-d"],
        ["-b", "3"],
        ["-b", "3", "-l"],
        ["-b", "2kb", "-l"],
        ["-b", "2kb", "-r", "chr2:1000-9000"],
        ["-r", "chr1"],
        ["-r", "chr1:1000-8000;chr3:2000-6000", "-c", "Blues"],
    ],
    ids=[
        "ice",
        "trim-despeckle",
        "bin-factor",
        "bin-factor-lines",
        "bin-bp-lines",
        "bin-bp-region",
        "chrom",
        "two-regions",
    ],
)
def test_view_options(runner, matrices, plotted, synthetic_lib, tmp_path, fmt, options):
    files = matrices[fmt]
    out = tmp_path / "map.png"
    invoke(runner, "view", "-D", "30", "-f", files["frags"], *options, "-o", out, files["mat"])
    assert out.stat().st_size > 0
    ((array, kwargs),) = plotted
    frags = synthetic_lib.fragments()
    n = len(frags)
    shapes = {
        "-n": (n, n),
        "-b 3": (-(-n // 3), -(-n // 3)),
        "-r chr1": ((frags.chrom == "chr1").sum(),) * 2,
    }
    key = " ".join(options[:2])
    if key in shapes:
        assert array.shape == shapes[key]
    elif ";" in key:
        assert array.shape[0] != array.shape[1]
        assert kwargs["cmap"] == "Blues"
    else:
        assert array.shape[0] == array.shape[1] <= n
    if "-l" in options:
        assert len(kwargs["chrom_starts"]) == 2
    if key == "-b 3" and "-l" in options:
        # Each bin is named after its first fragment
        chroms = frags.chrom[::3].to_numpy()
        expected = np.flatnonzero(chroms[1:] != chroms[:-1]) + 1
        np.testing.assert_array_equal(kwargs["chrom_starts"], expected)


def test_view_region_bounds(runner, matrices, plotted, synthetic_lib, tmp_path):
    """A region (1-based, inclusive) holds every bin it overlaps, the last one included."""
    files = matrices["cool"]
    invoke(runner, "view", "-r", "chr2:1000-9000", "-f", files["frags"], "-o", tmp_path / "x.png",
           files["mat"])  # fmt: skip
    frags = synthetic_lib.fragments()
    chr2 = frags[frags.chrom == "chr2"]
    overlap = ((chr2.end_pos >= 1000) & (chr2.start_pos < 9000)).sum()
    assert plotted[0][0].shape == (overlap, overlap)


def test_cli_doctests():
    import doctest

    import hicstuff.cli

    assert doctest.testmod(hicstuff.cli).failed == 0


def test_view_ratio(runner, matrices, plotted, tmp_path, caplog):
    """The log2 ratio of a map with itself is empty and uses a divergent map."""
    out = tmp_path / "ratio.png"
    mats = matrices["cool"]["mat"], matrices["bg2"]["mat"]
    frags = matrices["bg2"]["frags"]
    invoke(runner, "view", "-D", "30", "-c", "Blues", "-f", frags, "-o", out, *mats)
    array, kwargs = plotted[0]
    assert np.isinf(array).all()
    assert (kwargs["cmap"], kwargs["vmin"], kwargs["vmax"]) == ("seismic", -2, 2)
    assert "Non-divergent colormap" in caplog.text


def test_view_binned_graal_without_fragments(runner, matrices, plotted, synthetic_lib, tmp_path):
    invoke(runner, "view", "-b", "4", "-o", tmp_path / "x.png", matrices["graal"]["mat"])
    n = len(synthetic_lib.fragments())
    assert plotted[0][0].shape == (-(-n // 4),) * 2


@pytest.mark.parametrize(
    "options, exit_code, message",
    [
        (["-b", "lots"], 2, "Invalid binning value"),
        (["-b", "1kb"], 2, "A fragment file (--frags) is required"),
        (["-r", "chr1:1-5000"], 2, "required to zoom"),
    ],
)
def test_view_usage_errors(runner, matrices, tmp_path, options, exit_code, message):
    result = invoke(
        runner, "view", *options, "-o", tmp_path / "x.png", matrices["graal"]["mat"],
        exit_code=exit_code,
    )  # fmt: skip
    assert message in text(result)


def test_view_region_with_lines_is_rejected(runner, matrices, tmp_path):
    files = matrices["cool"]
    result = invoke(
        runner, "view", "-l", "-r", "chr1", "-f", files["frags"], "-o", tmp_path / "x.png",
        files["mat"], exit_code=1,
    )  # fmt: skip
    assert isinstance(result.exception, NotImplementedError)


# ---------------------------------------------------------------------------
# scalogram, rebin, subsample, convert
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "options",
    [
        ["-i", "10-60", "-r", "2-20", "-d", "-n"],
        ["-i", "chr2:1000-9000", "-r", "1kb-5kb", "-C", "magma"],
    ],
)
def test_scalogram(runner, matrices, tmp_path, options):
    files = matrices["cool"]
    out = tmp_path / "scalo.png"
    invoke(runner, "scalogram", "-f", files["frags"], *options, "-o", out, files["mat"])
    assert out.stat().st_size > 0


@pytest.mark.parametrize(
    "fmt, binning",
    [("cool", "2kb"), ("cool", "250bp"), ("graal", "3"), ("graal", "2kb"), ("bg2", "2")],
)
def test_rebin_preserves_contacts(runner, matrices, tmp_path, fmt, binning):
    files = matrices[fmt]
    prefix = tmp_path / "sub" / "rebinned"
    invoke(
        runner, "rebin", "-b", binning, "-f", files["frags"], "-c", files["chroms"],
        files["mat"], prefix,
    )  # fmt: skip
    before, _, _ = hio.flexible_hic_loader(files["mat"], files["frags"], files["chroms"])
    ext = {"cool": ".cool", "bg2": ".bg2", "graal": ".mat.tsv"}[fmt]
    after_frags = {"graal": f"{prefix}.frags.tsv", "bg2": None, "cool": None}[fmt]
    after, frags, _ = hio.flexible_hic_loader(
        f"{prefix}{ext}", after_frags, f"{prefix}.chr.tsv", quiet=True
    )
    assert after.sum() == before.sum()
    if binning == "2kb":
        # One bin per 2 kb window up to the last restriction fragment.
        starts = pd.read_csv(files["frags"], sep="\t").groupby("chrom").start_pos.max()
        assert len(frags) == (starts // 2000 + 1).sum()
    # Existing output requires --force
    invoke(runner, "rebin", "-b", binning, files["mat"], prefix, exit_code=1)


def test_rebin_errors(runner, matrices, tmp_path):
    files = matrices["graal"]
    result = invoke(runner, "rebin", "-b", "2", files["mat"], tmp_path / "x", exit_code=2)
    assert "Graal format requires --frags" in text(result)
    result = invoke(
        runner, "rebin", "-b", "huge", "-f", files["frags"], "-c", files["chroms"], files["mat"],
        tmp_path / "x", exit_code=2,
    )  # fmt: skip
    assert "Invalid binning" in text(result)


@pytest.mark.parametrize(
    "prop, expected", [("0.5", lambda n: int(0.5 * n)), ("1000", lambda n: 1000)]
)
@pytest.mark.parametrize("fmt", ["cool", "bg2", "graal"])
def test_subsample(runner, matrices, tmp_path, fmt, prop, expected):
    files = matrices[fmt]
    prefix = tmp_path / "sub"
    invoke(runner, "subsample", "-p", prop, files["mat"], prefix)
    ext = {"cool": ".cool", "bg2": ".bg2", "graal": ".mat.tsv"}[fmt]
    before = _pixels(files["mat"], files["frags"])
    after = _pixels(f"{prefix}{ext}", files["frags"])
    assert after.sum() == expected(before.sum())
    # Subsampled contacts come from existing pixels, at most as many
    assert after.index.isin(before.index).all()
    assert (after <= before.reindex(after.index)).all()
    invoke(runner, "subsample", "-p", prop, files["mat"], prefix, exit_code=1)
    invoke(runner, "subsample", "-F", "-p", prop, files["mat"], prefix)


def test_convert(runner, matrices, synthetic_lib, tmp_path):
    files = matrices["graal"]
    reference, _, _ = hio.flexible_hic_loader(files["mat"], files["frags"], files["chroms"])
    # GC content is recomputed from the genome (see test_convert_gc_with_hicstuff_fragments)
    frags = tmp_path / "frags.tsv"
    pd.read_csv(files["frags"], sep="\t").drop(columns="gc_content").to_csv(
        frags, sep="\t", index=False
    )
    for fmt, ext in (("cool", ".cool"), ("bg2", ".bg2"), ("graal", ".mat.tsv")):
        prefix = tmp_path / fmt
        invoke(
            runner, "convert", "-T", fmt, "-g", synthetic_lib.genome, "-f", frags,
            "-c", files["chroms"], files["mat"], prefix,
        )  # fmt: skip
        mat, _, _ = hio.flexible_hic_loader(f"{prefix}{ext}", files["frags"], quiet=True)
        np.testing.assert_array_equal(mat.toarray(), reference.toarray())
    # GC content of each bin, computed from the genome
    bins = cooler.Cooler(str(tmp_path / "cool.cool")).bins()[:]
    gc = [
        sum(synthetic_lib.chroms[c][s:e].count(b) for b in "GC") / (e - s)
        for c, s, e in zip(bins.chrom, bins.start, bins.end)
    ]
    np.testing.assert_allclose(bins.gc_content, gc)
    # cool without size column, to bg2
    invoke(runner, "convert", "-T", "bg2", matrices["cool"]["mat"], tmp_path / "from_cool")
    mat, _, _ = hio.flexible_hic_loader(str(tmp_path / "from_cool.bg2"), files["frags"])
    np.testing.assert_array_equal(mat.toarray(), reference.toarray())
    invoke(
        runner, "convert", "-T", "bg2", matrices["cool"]["mat"], tmp_path / "from_cool", exit_code=1
    )


def test_convert_gc_with_hicstuff_fragments(runner, matrices, synthetic_lib, tmp_path):
    """-g recomputes the gc_content column of hicstuff's fragments_list.txt."""
    files = matrices["graal"]
    frags = pd.read_csv(files["frags"], sep="\t")
    altered = tmp_path / "frags.tsv"
    frags.assign(gc_content=-1.0).to_csv(altered, sep="\t", index=False)
    invoke(
        runner, "convert", "-T", "cool", "-g", synthetic_lib.genome, "-f", altered,
        "-c", files["chroms"], files["mat"], tmp_path / "out",
    )  # fmt: skip
    bins = cooler.Cooler(str(tmp_path / "out.cool")).bins()[:]
    np.testing.assert_allclose(bins.gc_content, frags.gc_content)


# ---------------------------------------------------------------------------
# distancelaw
# ---------------------------------------------------------------------------


def test_distancelaw_command(runner, fastq_run, tmp_path):
    run = fastq_run()
    pairs, frags = run.tmp("valid_idx.pairs"), run.path("fragments_list.txt")
    table, image = tmp_path / "ps.txt", tmp_path / "ps.png"
    invoke(
        runner, "distancelaw", "-p", pairs, "-f", frags, "-O", table, "-o", image,
        "-i", "1000", "-s", "15000",
    )  # fmt: skip
    assert image.stat().st_size > 0
    xs, ps, names = hcdl.import_distance_law(str(table))
    assert sorted({n[0] for n in names}) == ["chr1", "chr2", "chr3"]
    # Position-based distances when no fragments file is given
    image2 = tmp_path / "ps_pos.png"
    invoke(runner, "distancelaw", "-p", pairs, "-o", image2, "-O", tmp_path / "pos.txt")
    assert image2.stat().st_size > 0
    # Several tables can be compared once averaged
    image3 = tmp_path / "avg.png"
    tables = f"{table},{tmp_path / 'pos.txt'}"
    invoke(
        runner, "distancelaw", "-d", tables, "-a", "-b", "9000", "-l", "frags,positions",
        "-o", image3,
    )  # fmt: skip
    assert image3.stat().st_size > 0
    invoke(runner, "distancelaw", "-d", table, "-a", "-o", tmp_path / "avg1.png")


@pytest.mark.parametrize(
    "options, message",
    [
        ([], "Provide either --pairs or --dist-tbl"),
        (["-p", "{pairs}", "-d", "{pairs}"], "Cannot use --pairs and --dist-tbl"),
        (["-d", "{table},{table}"], "--average is required"),
    ],
    ids=["no-input", "both-inputs", "no-average"],
)
def test_distancelaw_usage_errors(runner, fastq_run, options, message):
    run = fastq_run()
    fmt = {"pairs": run.tmp("valid_idx.pairs"), "table": run.path("distance_law.txt")}
    result = invoke(runner, "distancelaw", *(o.format(**fmt) for o in options), exit_code=2)
    assert message in text(result)


# ---------------------------------------------------------------------------
# missview
# ---------------------------------------------------------------------------


def test_missview(runner, synthetic_lib, tmp_path, require_tools, plotted):
    genome = tmp_path / "small.fa"
    genome.write_text(">chr3\n" + synthetic_lib.chroms["chr3"][:4000] + "\n")
    out = tmp_path / "missview.png"
    invoke(
        runner, "missview", "-R", "100", "-a", "minimap2", "-b", "1000", "-T", tmp_path / "work",
        genome, out,
    )  # fmt: skip
    assert out.stat().st_size > 0
    # A unique sequence has no missing bin, the shorter last one included.
    array, kwargs = plotted[0]
    np.testing.assert_array_equal(array, np.ones((4, 4)))
    assert kwargs["title"].startswith("0.000% missing bins for small.fa with 100 bp reads")
    # Temporary files are kept in the requested directory, and reused with --force.
    assert (tmp_path / "work" / "abs_fragments_contacts_weighted.txt").exists()
    invoke(
        runner, "missview", "-R", "100", "-a", "minimap2", "-T", tmp_path / "work", genome, out,
        exit_code=1,
    )  # fmt: skip
    invoke(
        runner, "missview", "-R", "100", "-a", "minimap2", "-b", "2kb", "-F",
        "-T", tmp_path / "work", genome, out,
    )  # fmt: skip
    assert plotted[1][0].shape == (2, 2)


@pytest.mark.parametrize("aligner", ["bowtie2", "bwa"])
def test_missview_finds_repeats(runner, synthetic_lib, tmp_path, require_tools, plotted, aligner):
    """Bins of a sequence present twice in the genome cannot be mapped."""
    seq = synthetic_lib.chroms["chr3"]
    genome = tmp_path / "repeats.fa"
    genome.write_text(f">chr3\n{seq[:4000]}\n>copy\n{seq[1000:2000]}{seq[6000:8000]}\n")
    invoke(runner, "missview", "-R", "60", "-a", aligner, "-b", "1000", genome, tmp_path / "m.png")
    array, kwargs = plotted[0]
    missing = np.flatnonzero(array.sum(axis=0) == 0)
    # chr3:1000-2000 and copy:0-1000
    assert list(missing) == [1, 4]
    assert kwargs["title"].startswith("28.571% missing bins")
