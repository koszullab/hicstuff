"""hicstuff.stats on the logs of pipeline runs on the synthetic library."""

import json
import re

import pytest

import hicstuff.stats as hcs
from tests import synthetic as syn


def _n_pairs(path):
    return len(syn.read_pairs(path))


def test_pipeline_stats_match_outputs(fastq_run, synthetic_lib):
    run = fastq_run("bowtie2", "iterative")
    stats = hcs.get_pipeline_stats(run.log)

    n = synthetic_lib.params.n_pairs
    unmappable = int((synthetic_lib.truth.kind == "unmapped").sum())
    assert stats["Total read pairs"] == n
    # Every genomic read aligns: only the random mates are unmapped.
    assert stats["Mapped reads"] == 2 * n - unmappable
    assert stats["Unmapped reads"] == unmappable

    valid = _n_pairs(run.tmp("valid_idx.pairs"))
    filtered = _n_pairs(run.tmp("valid_idx_filtered.pairs"))
    final = _n_pairs(run.path("valid_idx_pcrfree.pairs.gz"))
    assert stats["Recovered contacts"] == valid == len(synthetic_lib.mapped_truth())
    assert stats["Filtered out"] == valid - filtered
    assert stats["Filtered out"] == stats["Loops"] + stats["Uncuts"] + stats["Weirds"]
    assert stats["Loops"] >= (synthetic_lib.truth.kind == "loop").sum()
    assert stats["Uncuts"] >= (synthetic_lib.truth.kind == "uncut").sum()
    assert stats["PCR duplicates"] == filtered - final
    assert stats["Final contacts"] == final
    assert stats["Removed contacts"] == valid - final


def test_stats_without_filtering(fastq_run, run_pipeline, synthetic_lib, indexed_genome):
    """Without event filtering nor PCR duplicate removal, nothing is removed."""
    run = run_pipeline(
        "iterative_unfiltered",
        genome=indexed_genome,
        input1=str(synthetic_lib.reads1),
        input2=str(synthetic_lib.reads2),
        mapping="iterative",
        enzyme=syn.LIBRARY_ENZYME,
        mat_fmt="graal",
    )
    stats = hcs.get_pipeline_stats(run.log)
    for key in ("Filtered out", "Loops", "Uncuts", "Weirds", "PCR duplicates", "Removed contacts"):
        assert stats[key] == 0
    assert stats["Final contacts"] == stats["Recovered contacts"]


@pytest.mark.parametrize(
    "mapping, extra",
    [
        ("normal", {}),
        ("iterative", {}),
        ("cutsite", {}),
        ("normal", {"skip_count": True}),
        ("cutsite", {"skip_count": True}),
    ],
    ids=["normal", "iterative", "cutsite", "normal-skip-count", "cutsite-skip-count"],
)
def test_stats_of_every_mapping_mode(fastq_run, synthetic_lib, tmp_path, mapping, extra):
    run = fastq_run("bowtie2", mapping, **extra)
    stats = hcs.get_pipeline_stats(run.log)
    n = synthetic_lib.params.n_pairs
    assert stats["Total read pairs"] == n
    # cutsite aligns the pieces of reads cut at ligation sites
    pieces = re.search(r"Number of pairs after digestion: (\d+)", run.log_text)
    n_aligned = 2 * (int(pieces[1]) if mapping == "cutsite" else n)
    assert stats["Mapped reads"] + stats["Unmapped reads"] == n_aligned
    assert stats["Unmapped reads"] >= (synthetic_lib.truth.kind == "unmapped").sum()
    assert stats["Recovered contacts"] == _n_pairs(run.tmp("valid_idx.pairs"))
    assert stats["Final contacts"] == _n_pairs(run.path("valid_idx_pcrfree.pairs.gz"))
    # The pipeline logs the stats it compiled
    assert "Stats not compiled" not in run.log_text
    # Mapped reads are a share of the reads aligned
    out = tmp_path / "stats.txt"
    hcs.write_pipeline_stats(stats, str(out))
    pct = float(re.search(r"Mapped reads: .*\(([\d.]+)%", out.read_text())[1])
    assert pct == round(100 * stats["Mapped reads"] / n_aligned, 2)


def test_stats_need_read_mapping(matrices):
    """A run started from pairs has no mapping statistics."""
    with pytest.raises(IndexError, match="No read mapping statistics"):
        hcs.get_pipeline_stats(matrices["graal"]["run"].log)


def _stats(**values):
    keys = [
        "Total read pairs",
        "Mapped reads",
        "Unmapped reads",
        "Recovered contacts",
        "Final contacts",
        "Removed contacts",
        "Filtered out",
        "Loops",
        "Uncuts",
        "Weirds",
        "PCR duplicates",
    ]
    stats = {"Sample": "lib"}
    stats.update({key: values.get(key, 0) for key in keys})
    return stats


def test_write_pipeline_stats(tmp_path, capsys, monkeypatch):
    stats = _stats(
        **{
            "Total read pairs": 2000,
            "Mapped reads": 3000,
            "Unmapped reads": 1000,
            "Recovered contacts": 1000,
            "Final contacts": 750,
            "Removed contacts": 250,
            "Filtered out": 200,
            "Loops": 50,
            "Uncuts": 150,
            "PCR duplicates": 50,
        }
    )
    txt, js = tmp_path / "stats.txt", tmp_path / "stats.json"
    hcs.write_pipeline_stats(stats, str(txt))
    hcs.write_pipeline_stats(stats, str(js))
    assert json.loads(js.read_text()) == stats
    lines = txt.read_text().splitlines()
    assert lines[0].split() == ["Sample:", "lib"]
    assert "Total read pairs:   2,000" in lines
    assert "Mapped reads:       3,000 (75.0%  of all reads)" in lines
    assert "Recovered contacts: 1,000 (50.0%  of all read pairs)" in lines
    assert "    Uncuts:         150 (15.0% of all contacts)" in lines
    assert "Final contacts:     750 (75.0% of all contacts)" in lines

    monkeypatch.chdir(tmp_path)
    hcs.print_pipeline_stats(stats)
    assert capsys.readouterr().out.strip() == txt.read_text().strip()
    # The temporary file used for printing is removed
    assert sorted(p.name for p in tmp_path.iterdir()) == ["stats.json", "stats.txt"]


def test_write_empty_stats(tmp_path):
    out = tmp_path / "empty.txt"
    hcs.write_pipeline_stats(_stats(), str(out))
    text = out.read_text()
    assert "Total read pairs:   0" in text
    assert text.count("(0.0%") == 10
