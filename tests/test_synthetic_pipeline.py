"""End-to-end pipeline tests on the synthetic library, checked against its truth.

Every aligner (bowtie2, bwa, minimap2) is run in every mapping mode (normal,
iterative, cutsite) once per session, with event filtering, PCR duplicate
removal and the distance law enabled; the tests below then compare what the
pipeline wrote with what was simulated.
"""

import gzip
import os
import re
import shutil

import cooler
import numpy as np
import pandas as pd
import pytest

import hicstuff.distance_law as hcdl
import hicstuff.io as hio
import hicstuff.pipeline as hpi
from tests import synthetic as syn

ALIGNERS = ("bowtie2", "bwa", "minimap2")
MAPPINGS = ("normal", "iterative", "cutsite")
COMBOS = [(aligner, mapping) for aligner in ALIGNERS for mapping in MAPPINGS]
#: Fewer mates align than simulated when an aligner gives short reads a
#: MAPQ below the pipeline's threshold (30): minimap2 does so for about half
#: of the unique 50-60 bp reads. Mates of at least this length align with
#: every aligner...
LONG_READ = 70
#: ... except when minimap2 aligns them iteratively, i.e. truncated to 20, 40,
#: 60 bp... first.
MIN_LONG_READ_RECALL = {("minimap2", "iterative"): 0.8}


@pytest.fixture(scope="module", params=COMBOS, ids=["-".join(c) for c in COMBOS])
def run(request, fastq_run):
    aligner, mapping = request.param
    return fastq_run(aligner, mapping)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def load_pairs(path):
    """Pairs file as a DataFrame, with cutsite suffixes split from read IDs.

    In cutsite mode, read ``x`` cut into segments yields pairs ``x:ij``; the
    pair of the first segments (``:00``) is the simulated contact.
    """
    pairs = syn.read_pairs(path)
    parts = pairs.readID.str.extract(r"^(?P<base>[^:]+)(?::(?P<sub>\d+))?$")
    pairs["base"] = parts["base"]
    pairs["sub"] = parts["sub"].fillna("00")
    return pairs


def log_value(run, pattern):
    match = re.search(pattern, run.log_text)
    assert match is not None, f"{pattern!r} not found in {run.log}"
    return tuple(int(x) for x in match.groups())


def classify(chr1, frag1, strand1, pos1, chr2, frag2, strand2, pos2, thr_uncut, thr_loop):
    """3C event type of a pair, following hicstuff's rules: mates are ordered
    by position, ``-+`` pairs closer than the loop threshold are loops, ``+-``
    pairs closer than the uncut threshold are uncuts."""
    if chr1 != chr2:
        return "inter"
    if (frag1, pos1) > (frag2, pos2):
        frag1, strand1, frag2, strand2 = frag2, strand2, frag1, strand1
    nsites = frag2 - frag1
    if nsites == 0 and strand1 == strand2:
        return "weird"
    if strand1 + strand2 == "-+" and nsites <= thr_loop:
        return "loop"
    if strand1 + strand2 == "+-" and nsites <= thr_uncut:
        return "uncut"
    return "intra"


def truth_class(row, thr_uncut, thr_loop):
    return classify(
        row.chr1,
        int(row.frag1),
        row.strand1,
        int(row.start1),
        row.chr2,
        int(row.frag2),
        row.strand2,
        int(row.start2),
        thr_uncut,
        thr_loop,
    )


def compare_mates(observed, truth, mapping):
    """Return the read IDs of ``observed`` pairs that disagree with the truth.

    Mates are compared after sorting by (fragment, strand): chromosome,
    fragment and strand must match exactly. Forward mates must start at the
    simulated position. In iterative mode, a reverse mate is reported at the
    leftmost base of its aligned 5' prefix, so it must lie between the
    simulated leftmost base and 19 bp before its 5' end.
    """
    merged = observed.merge(truth.add_suffix("_t"), left_on="base", right_index=True)
    bad = []
    for row in merged.itertuples():
        obs = sorted(
            [
                (int(row.frag1), row.strand1, row.chr1, int(row.pos1)),
                (int(row.frag2), row.strand2, row.chr2, int(row.pos2)),
            ]
        )
        exp = sorted(
            [
                (int(row.frag1_t), row.strand1_t, row.chr1_t, int(row.start1_t), int(row.five1_t)),
                (int(row.frag2_t), row.strand2_t, row.chr2_t, int(row.start2_t), int(row.five2_t)),
            ]
        )
        for (frag, strand, chrom, pos), (t_frag, t_strand, t_chrom, start, five) in zip(obs, exp):
            if (frag, strand, chrom) != (t_frag, t_strand, t_chrom):
                ok = False
            elif strand == "+" or mapping != "iterative":
                ok = pos == start + 1
            else:
                ok = start + 1 <= pos <= five + 1 - 19
            if not ok:
                bad.append(row.base)
                break
    return bad


# ---------------------------------------------------------------------------
# Alignment and pairing
# ---------------------------------------------------------------------------


def test_valid_pairs_match_truth(run, synthetic_lib):
    aligner, mapping = run.kwargs["aligner"], run.kwargs["mapping"]
    truth = synthetic_lib.truth.set_index("read_id")
    pairs = load_pairs(run.tmp("valid_idx.pairs"))
    main = pairs[pairs["sub"] == "00"].set_index("readID")

    # Precision: no pair from the unmappable library fraction, no duplicated
    # read, and every reported mate is where it was simulated.
    assert main.base.is_unique
    assert set(main.base) <= set(synthetic_lib.mapped_truth().read_id)
    assert compare_mates(main, truth, mapping) == []

    # Recall, by category of simulated pairs.
    found = set(main.base)
    expected = synthetic_lib.mapped_truth()

    def recall(subset):
        assert len(subset) > 50
        return len(found & set(subset.read_id)) / len(subset)

    plain = expected[~expected.chimeric]
    long_reads = plain[(plain.len1 >= LONG_READ) & (plain.len2 >= LONG_READ)]
    assert recall(long_reads) >= MIN_LONG_READ_RECALL.get((aligner, mapping), 0.97)
    if aligner != "minimap2":
        assert recall(plain) >= 0.99
        assert len(main) == pytest.approx(len(expected), rel=0.07)
    # Reads spanning a ligation junction are rescued by iterative alignment
    # (which aligns their 5' end) and by cutsite (which cuts the junction).
    chimeric = expected[expected.chimeric]
    if aligner == "bwa" or (aligner == "bowtie2" and mapping != "normal"):
        assert recall(chimeric) >= 0.99
        assert len(main) == pytest.approx(len(expected), abs=len(expected) // 100)

    # The pipeline logs what it found
    assert log_value(run, r"(\d+) pairs successfully mapped") == (len(pairs),)
    n_input_pairs = synthetic_lib.params.n_pairs
    if mapping == "cutsite":
        # Each chimeric R1 is cut in two, each piece paired with R2.
        n_input_pairs += int(synthetic_lib.truth.chimeric.sum())
    (_, n_reads) = log_value(run, r"mapped with Q >= 30 \((\d+)/(\d+)\)")
    assert n_reads == 2 * n_input_pairs


def test_cutsite_extra_pairs_are_uncut_like(fastq_run, synthetic_lib):
    """The second piece of a cut chimeric R1 is paired with R2: both come from
    the same fragment and face each other."""
    pairs = load_pairs(fastq_run("bowtie2", "cutsite").tmp("valid_idx.pairs"))
    extra = pairs[pairs["sub"] != "00"]
    truth = synthetic_lib.truth.set_index("read_id")
    assert set(extra["sub"]) == {"10"}
    assert set(extra.base) == set(truth.index[truth.chimeric])
    assert (extra.frag1 == extra.frag2).all()
    assert ((extra.strand1 + extra.strand2) == "+-").all()
    # R2 is the mate of the second piece too: same fragment as in the truth.
    assert (extra.frag1.values == truth.loc[extra.base, "frag2"].values).all()


# ---------------------------------------------------------------------------
# Filtering of 3C events and PCR duplicates
# ---------------------------------------------------------------------------


def test_filter_events_match_truth(run, synthetic_lib):
    thr_uncut, thr_loop = log_value(run, r"Filtering with thresholds: uncuts=(\d+) loops=(\d+)")
    # Simulated contacts are at least 1 kb (several fragments) apart: the
    # estimated thresholds separate them from the injected events.
    assert 2 <= thr_uncut <= 6
    assert 2 <= thr_loop <= 6

    truth = synthetic_lib.truth.set_index("read_id")
    pairs = load_pairs(run.tmp("valid_idx.pairs"))
    main = pairs[pairs["sub"] == "00"]
    extra = pairs[pairs["sub"] != "00"]
    # Event type of every pair, predicted from the truth (from the observed
    # mates for the extra pairs of cut chimeric reads), keyed by read ID.
    cols = ["chr1", "frag1", "strand1", "pos1", "chr2", "frag2", "strand2", "pos2"]
    expected = pd.concat(
        [
            pd.Series(
                [
                    truth_class(row, thr_uncut, thr_loop)
                    for row in truth.loc[main.base].itertuples()
                ],
                index=main.readID.values,
                dtype=object,
            ),
            pd.Series(
                [classify(*row[1:], thr_uncut, thr_loop) for row in extra[cols].itertuples()],
                index=extra.readID.values,
                dtype=object,
            ),
        ]
    )
    counts = expected.value_counts()
    discarded, loops, uncuts, weirds = log_value(
        run, r"(\d+) pairs discarded: Loops: (\d+), Uncuts: (\d+), Weirds: (\d+)"
    )
    assert loops == counts.get("loop", 0)
    assert uncuts == counts.get("uncut", 0)
    assert weirds == counts.get("weird", 0) == 0
    assert discarded == loops + uncuts + weirds

    # The discarded pairs are exactly the predicted ones, not just as many.
    kept = load_pairs(run.tmp("valid_idx_filtered.pairs"))
    assert len(kept) == len(pairs) - discarded
    removed = set(pairs.readID) - set(kept.readID)
    assert removed == set(expected.index[expected.isin(["loop", "uncut", "weird"])])

    # Each count is the injected events of that kind, all of them discarded,
    # plus the genuine contacts (or cutsite pieces) within the threshold,
    # which hicstuff cannot tell apart from them.
    kinds = pd.Series(truth.loc[main.base, "kind"].values, index=main.readID.values)
    params = synthetic_lib.params
    for kind, rate in (("loop", params.loop_rate), ("uncut", params.uncut_rate)):
        injected = kinds.index[kinds == kind]
        assert (expected.loc[injected] == kind).all()
        assert set(injected) <= removed
        # The injected events align at their simulated rate.
        assert len(injected) / len(main) == pytest.approx(rate, rel=0.25)


def test_pcr_duplicates_removed(run, synthetic_lib):
    truth = synthetic_lib.truth.set_index("read_id")
    before = load_pairs(run.tmp("valid_idx_filtered.pairs"))
    after = load_pairs(run.path("valid_idx_pcrfree.pairs.gz"))
    # A duplicate and its original form one molecule; cutsite pieces of one
    # molecule are distinct pairs.
    molecule = pd.Series(np.where(truth.dup_of != "", truth.dup_of, truth.index), truth.index)
    before["molecule"] = molecule.loc[before.base].values + ":" + before["sub"]
    after["molecule"] = molecule.loc[after.base].values + ":" + after["sub"]

    assert after.molecule.is_unique
    assert set(after.molecule) == set(before.molecule)
    removed = len(before) - len(after)
    assert removed == len(before) - before.molecule.nunique()
    injected = truth.kind.loc[before.base].eq("duplicate").sum()
    assert 0 < removed <= injected
    assert log_value(run, r"\((\d+)/(\d+) pairs\)") == (removed, len(before))


# ---------------------------------------------------------------------------
# Matrices
# ---------------------------------------------------------------------------


def _pixel_counts(pairs, bin1, bin2):
    lo, hi = np.minimum(bin1, bin2), np.maximum(bin1, bin2)
    return pd.DataFrame({"bin1_id": lo, "bin2_id": hi}).value_counts().sort_index()


def test_cool_matches_final_pairs(run, synthetic_lib):
    pairs = load_pairs(run.path("valid_idx_pcrfree.pairs.gz"))
    clr = cooler.Cooler(str(run.path("abs_fragments_contacts_weighted.cool")))
    bins = clr.bins()[:]
    frags = synthetic_lib.fragments()
    assert clr.chromnames == synthetic_lib.chrom_names
    np.testing.assert_array_equal(bins.start, frags.start_pos)
    np.testing.assert_array_equal(bins.end, frags.end_pos)

    pixels = clr.pixels()[:]
    assert pixels["count"].sum() == len(pairs)
    observed = pixels.set_index(["bin1_id", "bin2_id"])["count"].sort_index()
    expected = _pixel_counts(pairs, pairs.frag1.values, pairs.frag2.values)
    np.testing.assert_array_equal(observed.index.to_numpy(), expected.index.to_numpy())
    np.testing.assert_array_equal(observed.to_numpy(), expected.to_numpy())

    # Inter-chromosomal contacts: as many as simulated, give or take what
    # the filters removed.
    trans = (pairs.chr1 != pairs.chr2).mean()
    assert trans == pytest.approx(synthetic_lib.params.trans_fraction, abs=0.03)


def test_final_pairs_sorted(run):
    with gzip.open(run.path("valid_idx_pcrfree.pairs.gz"), "rt") as handle:
        header = [line for line in handle if line.startswith("#")]
    assert "#sorted: chr1-chr2-pos1-pos2\n" in header
    pairs = load_pairs(run.path("valid_idx_pcrfree.pairs.gz"))
    keys = list(zip(pairs.chr1, pairs.chr2, pairs.pos1, pairs.pos2))
    assert keys == sorted(keys)


def test_digestion_matches_truth(fastq_run, synthetic_lib):
    run = fastq_run()
    frags = pd.read_csv(run.path("fragments_list.txt"), sep="\t")
    expected = synthetic_lib.fragments()
    pd.testing.assert_frame_equal(frags.drop(columns="gc_content"), expected)
    chroms = pd.read_csv(run.path("info_contigs.txt"), sep="\t")
    assert list(chroms.contig) == synthetic_lib.chrom_names
    assert list(chroms.length) == list(synthetic_lib.chrom_sizes.values())
    n_frags = expected.groupby("chrom", sort=False).size()
    assert list(chroms.n_frags) == list(n_frags)
    assert list(chroms.cumul_length) == [0] + list(np.cumsum(n_frags)[:-1])


def test_fragment_gc_content(fastq_run, synthetic_lib):
    frags = pd.read_csv(fastq_run().path("fragments_list.txt"), sep="\t")
    seqs = synthetic_lib.chroms
    gc = [
        sum(seqs[c][s:e].count(b) for b in "GC") / (e - s)
        for c, s, e in zip(frags.chrom, frags.start_pos, frags.end_pos)
    ]
    np.testing.assert_allclose(frags.gc_content, gc, rtol=1e-6)


# ---------------------------------------------------------------------------
# Distance law
# ---------------------------------------------------------------------------


def _fit_slope(xs, ps, lo=2_000, hi=9_000):
    xs, ps = np.asarray(xs, dtype=float), np.asarray(ps, dtype=float)
    keep = (xs >= lo) & (xs <= hi) & (ps > 0)
    assert keep.sum() >= 10
    return np.polyfit(np.log(xs[keep]), np.log(ps[keep]), 1)[0]


def test_distance_law_recovers_slope(fastq_run, synthetic_lib):
    """P(s) decays with the simulated exponent (between 2 and 9 kb, where
    every chromosome contributes)."""
    alpha = synthetic_lib.params.alpha
    run = fastq_run("bowtie2", "iterative")
    xs, ps, names = hcdl.import_distance_law(str(run.path("distance_law.txt")))
    assert [n[0] for n in names] == synthetic_lib.chrom_names
    avg_xs, avg_ps = hcdl.average_distance_law(xs, ps, 9_000)
    assert _fit_slope(avg_xs, avg_ps) == pytest.approx(-alpha, abs=0.2)

    # Without a fragments file, distances are taken between read positions.
    xs, ps, names = hcdl.get_distance_law(str(run.tmp("valid_idx.pairs")), None)
    assert names == synthetic_lib.chrom_names
    avg_xs, avg_ps = hcdl.average_distance_law(xs, ps, 9_000)
    assert _fit_slope(avg_xs, avg_ps) == pytest.approx(-alpha, abs=0.2)


# ---------------------------------------------------------------------------
# Other start stages, formats and options
# ---------------------------------------------------------------------------


def _sorted_lines(path):
    with hio.read_compressed(str(path)) as handle:
        return sorted(line for line in handle if not line.startswith("#"))


@pytest.mark.parametrize("stage", ["bam", "pairs"])
def test_start_stages_reproduce_fastq_run(stage, fastq_run, indexed_genome, run_pipeline):
    ref = fastq_run()
    inputs = {
        "bam": (str(ref.tmp("for.bam")), str(ref.tmp("rev.bam"))),
        "pairs": (str(ref.tmp("valid.pairs")), None),
    }[stage]
    run = run_pipeline(
        f"from_{stage}",
        genome=indexed_genome,
        input1=inputs[0],
        input2=inputs[1],
        start_stage=stage,
        enzyme=syn.LIBRARY_ENZYME,
        mat_fmt="graal",
    )
    assert _sorted_lines(run.path("valid_idx.pairs.gz")) == _sorted_lines(
        ref.tmp("valid_idx.pairs")
    )
    # Inputs are not deleted, intermediate files are.
    assert all(path is None or os.path.exists(path) for path in inputs)
    assert not run.tmp("valid_idx.pairs").exists()


def test_graal_bg2_and_cool_agree(matrices, final_pairs, synthetic_lib):
    pairs = load_pairs(final_pairs)
    n = len(synthetic_lib.fragments())
    expected = np.zeros((n, n), dtype=int)
    np.add.at(
        expected,
        (np.minimum(pairs.frag1, pairs.frag2), np.maximum(pairs.frag1, pairs.frag2)),
        1,
    )
    graal = hio.load_sparse_matrix(matrices["graal"]["mat"])
    bg2, bg2_frags, _ = hio.load_bedgraph2d(
        matrices["bg2"]["mat"], fragments_file=matrices["graal"]["frags"]
    )
    cool, cool_frags, _ = hio.load_cool(matrices["cool"]["mat"])
    for mat in (graal, bg2, cool):
        assert mat.shape == (n, n)
        np.testing.assert_array_equal(mat.toarray(), expected)
    np.testing.assert_array_equal(bg2_frags.start_pos, cool_frags.start_pos)

    # graal header: rows, columns, nonzero entries
    with open(matrices["graal"]["mat"]) as handle:
        header = handle.readline().rstrip("\n")
    assert header.split("\t")[:2] == [str(n), str(n)]


@pytest.mark.parametrize("fmt", ["graal", "bg2"])
def test_text_matrix_counts(matrices, final_pairs, synthetic_lib, fmt):
    """The graal header and the log give the number of entries and of pairs."""
    with open(matrices[fmt]["mat"]) as handle:
        lines = handle.read().splitlines()
    nnz = len(lines) - (fmt == "graal")
    if fmt == "graal":
        assert int(lines[0].split("\t")[2]) == nnz
    n_pairs, n = len(load_pairs(final_pairs)), len(synthetic_lib.fragments())
    message = f"{n_pairs} pairs used to build a contact map of {n} bins with {nnz} nonzero entries"
    assert message in matrices[fmt]["run"].log_text


@pytest.mark.parametrize("fmt", ["graal", "bg2"])
def test_text_matrix_from_empty_pairs(tmp_path, fmt):
    pairs = tmp_path / "empty.pairs"
    pairs.write_text(
        "## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 frag1 frag2\n"
    )
    shutil.copy("test_data/fragments_list.txt", tmp_path / "frags.txt")
    mat = tmp_path / "mat.txt"
    hpi.pairs2matrix(str(pairs), str(mat), str(tmp_path / "frags.txt"), mat_fmt=fmt)
    n_frags = len(pd.read_csv(tmp_path / "frags.txt", sep="\t"))
    assert mat.read_text() == (f"{n_frags}\t{n_frags}\t0\n" if fmt == "graal" else "")


def test_binned_mcool(final_pairs, synthetic_lib, run_pipeline, tmp_path):
    src = tmp_path / "input.pairs"
    src.write_text(final_pairs.read_text())
    run = run_pipeline(
        "binned_mcool",
        genome=str(synthetic_lib.genome),
        input1=str(src),
        start_stage="pairs_idx",
        enzyme=syn.LIBRARY_ENZYME,
        binning=1000,
    )
    pairs = load_pairs(final_pairs)
    mcool = str(run.path("abs_fragments_contacts_weighted.mcool"))
    # The single-resolution cool is removed after zoomify.
    assert not run.path("abs_fragments_contacts_weighted.cool").exists()
    resolutions = sorted(int(p.split("/")[-1]) for p in cooler.fileops.list_coolers(mcool))
    assert resolutions == [1000, 2000, 5000, 10000, 20000, 50000, 100000]
    sizes = synthetic_lib.chrom_sizes
    for res in resolutions:
        clr = cooler.Cooler(f"{mcool}::/resolutions/{res}")
        bins = clr.bins()[:]
        assert (bins.start % res == 0).all()
        assert len(bins) == sum(-(-size // res) for size in sizes.values())
        assert clr.pixels()[:]["count"].sum() == len(pairs)

    # Base resolution: pixels are the pairs binned at 1 kb, and balanced.
    clr = cooler.Cooler(f"{mcool}::/resolutions/1000")
    offsets = clr.offset
    bin1 = np.array([offsets(c) for c in pairs.chr1]) + (pairs.pos1.values - 1) // 1000
    bin2 = np.array([offsets(c) for c in pairs.chr2]) + (pairs.pos2.values - 1) // 1000
    expected = _pixel_counts(pairs, bin1, bin2)
    observed = clr.pixels()[:].set_index(["bin1_id", "bin2_id"])["count"].sort_index()
    np.testing.assert_array_equal(observed.index.to_numpy(), expected.index.to_numpy())
    np.testing.assert_array_equal(observed.to_numpy(), expected.to_numpy())
    assert "weight" in clr.bins().columns


def test_binned_cool_without_zoomify(final_pairs, synthetic_lib, run_pipeline, tmp_path):
    src = tmp_path / "input.pairs"
    src.write_text(final_pairs.read_text())
    run = run_pipeline(
        "binned_cool",
        genome=str(synthetic_lib.genome),
        input1=str(src),
        start_stage="pairs_idx",
        enzyme=syn.LIBRARY_ENZYME,
        binning=2000,
        zoomify=False,
        exclude="chr3",
        balancing_args="--max-iters 10 --min-nnz 1",
    )
    pairs = load_pairs(final_pairs)
    clr = cooler.Cooler(str(run.path("abs_fragments_contacts_weighted.cool")))
    assert clr.binsize == 2000
    assert clr.chromnames == ["chr1", "chr2"]
    kept = pairs[(pairs.chr1 != "chr3") & (pairs.chr2 != "chr3")]
    assert clr.pixels()[:]["count"].sum() == len(kept)
    assert "weight" in clr.bins().columns


def test_plots_and_prefix(fastq_run, run_pipeline, indexed_genome):
    """Plots are written for every step; outputs carry the prefix."""
    ref = fastq_run()
    run = run_pipeline(
        "plots",
        genome=indexed_genome,
        input1=str(ref.tmp("for.bam")),
        input2=str(ref.tmp("rev.bam")),
        start_stage="bam",
        enzyme=syn.LIBRARY_ENZYME,
        prefix="synth",
        plot=True,
        filter_events=True,
        pcr_duplicates=True,
        distance_law=True,
        mat_fmt="bg2",
    )
    plots = run.out_dir / "plots"
    for name in ("frags_hist", "event_distance", "event_distribution", "distance_law"):
        assert (plots / f"synth_{name}.pdf").stat().st_size > 0
    for name in ("frags.tsv", "chr.tsv", "bg2", "distance_law.txt"):
        assert run.path(name).exists()
    # Intermediate files are cleaned up, BAM inputs are kept.
    assert not run.tmp("valid_idx.pairs").exists()
    assert ref.tmp("for.bam").exists()


def test_cutsite_without_enzyme_falls_back_to_iterative(
    synthetic_lib, run_pipeline, indexed_genome
):
    """Fixed-size chunks have no ligation site: cutsite mode aligns
    iteratively, and fragments are the chunks."""
    run = run_pipeline(
        "cutsite_chunks",
        genome=indexed_genome,
        input1=str(synthetic_lib.reads1),
        input2=str(synthetic_lib.reads2),
        mapping="cutsite",
        enzyme=2000,
        mat_fmt="graal",
        skip_count=True,
    )
    assert "iterative mapping will be used instead" in run.log_text
    assert "Skipping read count check" in run.log_text
    frags = pd.read_csv(run.path("fragments_list.txt"), sep="\t")
    sizes = synthetic_lib.chrom_sizes
    assert len(frags) == sum(-(-size // 2000) for size in sizes.values())
    mat = hio.load_sparse_matrix(str(run.path("abs_fragments_contacts_weighted.txt")))
    pairs = load_pairs(run.path("valid_idx.pairs.gz"))
    assert mat.sum() == len(pairs) == len(synthetic_lib.mapped_truth())


def test_pipeline_rejects_bad_inputs(synthetic_lib, indexed_genome, tmp_path, monkeypatch):
    common = dict(genome=indexed_genome, out_dir=str(tmp_path / "out"), enzyme="DpnII")
    reads = dict(input1=str(synthetic_lib.reads1), input2=str(synthetic_lib.reads2))
    with pytest.raises(ValueError, match="aligner"):
        hpi.full_pipeline(aligner="blast", **reads, **common)
    with pytest.raises(SystemExit):
        hpi.full_pipeline(input1=str(synthetic_lib.reads1), **common)
    with pytest.raises(SystemExit):
        hpi.full_pipeline(start_stage="pairs", **reads, **common)
    with pytest.raises(ValueError):
        hpi.full_pipeline(mapping="psychic", force=True, **reads, **common)
    # Existing outputs and temporary files are protected unless force is set.
    (tmp_path / "out" / "fragments_list.txt").touch()
    with pytest.raises(OSError, match="Output file .* already exists"):
        hpi.full_pipeline(**reads, **common)
    (tmp_path / "out2" / "tmp").mkdir(parents=True)
    (tmp_path / "out2" / "tmp" / "for.bam").touch()
    with pytest.raises(OSError, match="Temporary file .* already exists"):
        hpi.full_pipeline(**reads, **dict(common, out_dir=str(tmp_path / "out2")))
    monkeypatch.setattr(hpi, "which", lambda tool: None if tool == "samtools" else tool)
    with pytest.raises(ImportError, match="samtools"):
        hpi.full_pipeline(**reads, **common)


@pytest.mark.parametrize("folder", [".", "tmp"])
def test_pairs_input_inside_out_dir(final_pairs, synthetic_lib, tmp_path, folder):
    """An input pairs file in out_dir, or in its tmp folder, is compressed into
    out_dir and neither moved nor deleted."""
    src = tmp_path / folder / "input.pairs"
    src.parent.mkdir(exist_ok=True)
    src.write_text(final_pairs.read_text())
    hpi.full_pipeline(
        genome=str(synthetic_lib.genome),
        input1=str(src),
        start_stage="pairs_idx",
        enzyme=syn.LIBRARY_ENZYME,
        mat_fmt="graal",
        out_dir=str(tmp_path),
    )
    assert _sorted_lines(tmp_path / "input.pairs.gz") == _sorted_lines(final_pairs)
    assert _sorted_lines(src) == _sorted_lines(final_pairs)


def test_tmp_dir_is_out_dir(final_pairs, synthetic_lib, tmp_path):
    """Temporary files may be written into out_dir itself."""
    src = tmp_path / "in" / "input.pairs"
    src.parent.mkdir()
    src.write_text(final_pairs.read_text())
    out = tmp_path / "out"
    hpi.full_pipeline(
        genome=str(synthetic_lib.genome),
        input1=str(src),
        start_stage="pairs_idx",
        enzyme=syn.LIBRARY_ENZYME,
        mat_fmt="graal",
        pcr_duplicates=True,
        out_dir=str(out),
        tmp_dir=str(out),
    )
    assert _sorted_lines(out / "valid_idx_pcrfree.pairs.gz") == _sorted_lines(final_pairs)
    assert not (out / "valid_idx_pcrfree.pairs").exists()
    assert src.exists()


@pytest.mark.parametrize("linked", ["out_dir", "tmp_dir"])
def test_final_pairs_reached_through_a_link(final_pairs, synthetic_lib, tmp_path, linked):
    """out_dir, or a tmp_dir equal to it, may be named through a symbolic link:
    a final pairs file already in out_dir is still recognised as such."""
    real, link = tmp_path / "real", tmp_path / "link"
    real.mkdir()
    link.symlink_to(real, target_is_directory=True)
    if linked == "out_dir":
        # The input pairs are in out_dir, which is named through the link
        src, final = real / "input.pairs", "input.pairs"
        dirs = dict(out_dir=str(link))
    else:
        # The final pairs are a temporary file in out_dir, as tmp_dir
        src, final = tmp_path / "input.pairs", "valid_idx_pcrfree.pairs"
        dirs = dict(out_dir=str(real), tmp_dir=str(link))
    src.write_text(final_pairs.read_text())
    hpi.full_pipeline(
        genome=str(synthetic_lib.genome),
        input1=str(src),
        start_stage="pairs_idx",
        enzyme=syn.LIBRARY_ENZYME,
        mat_fmt="graal",
        pcr_duplicates=linked == "tmp_dir",
        **dirs,
    )
    assert _sorted_lines(real / f"{final}.gz") == _sorted_lines(final_pairs)
    assert _sorted_lines(src) == _sorted_lines(final_pairs)
    assert (real / final).exists() == (linked == "out_dir")


def _bam_subset(src, dst, keep):
    import pysam

    with pysam.AlignmentFile(src, "rb") as bam, pysam.AlignmentFile(dst, "wb", template=bam) as out:
        for i, read in enumerate(bam):
            if keep(i):
                out.write(read)


def _check_read_counts(log, *bams, min_qual=30):
    """bam2pairs logs every read once, with the number passing min_qual."""
    import pysam

    mapped = total = 0
    for path in bams:
        with pysam.AlignmentFile(path, "rb") as bam:
            for read in bam:
                total += 1
                mapped += read.mapping_quality >= min_qual
    assert f"mapped with Q >= {min_qual} ({mapped}/{total})" in log


def test_bam2pairs_skips_reads_missing_from_one_file(fastq_run, tmp_path, caplog):
    """Mates are matched by name: reads absent from one BAM are skipped."""
    run = fastq_run()
    rev = str(tmp_path / "rev.bam")
    _bam_subset(str(run.tmp("rev.bam")), rev, lambda i: i % 10 != 3)
    out = tmp_path / "out.pairs"
    hpi.bam2pairs(str(run.tmp("for.bam")), rev, str(out), str(run.path("info_contigs.txt")))
    ref = syn.read_pairs(run.tmp("valid.pairs"))
    # Reads are sorted by name and named after their rank.
    n_reads = hio.check_bam_entries(str(run.tmp("rev.bam")))
    dropped = {f"sim{i:06d}" for i in range(3, n_reads, 10)}
    pairs = syn.read_pairs(out)
    assert set(pairs.readID) == set(ref.readID) - dropped
    pd.testing.assert_frame_equal(
        pairs.reset_index(drop=True),
        ref[~ref.readID.isin(dropped)].reset_index(drop=True),
    )
    assert f"{len(dropped)} reads were only present in one BAM file" in caplog.text
    _check_read_counts(caplog.text, str(run.tmp("for.bam")), rev)


@pytest.mark.parametrize("truncated", ["for", "rev"])
def test_bam2pairs_with_truncated_bam(fastq_run, tmp_path, caplog, truncated):
    """Reads after the end of the shorter BAM have no mate and are skipped."""
    run = fastq_run()
    bams = {end: str(run.tmp(f"{end}.bam")) for end in ("for", "rev")}
    n_reads = hio.check_bam_entries(bams[truncated])
    bams[truncated] = str(tmp_path / "truncated.bam")
    _bam_subset(str(run.tmp(f"{truncated}.bam")), bams[truncated], lambda i: i < n_reads - 10)
    out = tmp_path / "out.pairs"
    hpi.bam2pairs(bams["for"], bams["rev"], str(out), str(run.path("info_contigs.txt")))
    ref = syn.read_pairs(run.tmp("valid.pairs"))
    dropped = {f"sim{i:06d}" for i in range(n_reads - 10, n_reads)}
    pd.testing.assert_frame_equal(
        syn.read_pairs(out).reset_index(drop=True),
        ref[~ref.readID.isin(dropped)].reset_index(drop=True),
    )
    assert "10 reads were only present in one BAM file" in caplog.text
    _check_read_counts(caplog.text, bams["for"], bams["rev"])


def test_log_file_handler_leaves_other_handlers(tmp_path):
    """set_file_handler replaces the handler it set before, and only that one."""
    import logging

    import hicstuff.log as hcl

    before = list(hcl.logger.handlers)
    other = logging.FileHandler(tmp_path / "other.log")
    hcl.logger.addHandler(other)
    try:
        hcl.set_file_handler(str(tmp_path / "first.log"))
        hcl.logger.info("message one")
        hcl.set_file_handler(str(tmp_path / "second.log"))
        hcl.logger.info("message two")
        assert other in hcl.logger.handlers
        # So do pytest's own handlers, a FileHandler among them
        assert all(handler in hcl.logger.handlers for handler in before)
        files = [h.baseFilename for h in hcl.logger.handlers if isinstance(h, logging.FileHandler)]
        assert sorted(os.path.basename(f) for f in files if str(tmp_path) in f) == [
            "other.log",
            "second.log",
        ]
    finally:
        hcl.logger.removeHandler(other)
        other.close()
    assert "message one" in (tmp_path / "first.log").read_text()
    assert "message two" not in (tmp_path / "first.log").read_text()
    assert "message two" in (tmp_path / "second.log").read_text()
    assert "message one" in (tmp_path / "other.log").read_text()
    assert "message two" in (tmp_path / "other.log").read_text()
