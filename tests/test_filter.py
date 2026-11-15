# Tests for the hicstuff filter module.
# 20190402

import hashlib
import os
from tempfile import NamedTemporaryFile

import numpy as np

from hicstuff import filter as hcf

PAIRS = "test_data/valid_idx.pairs"


def hash_file(filename):
    """Computes the MD5 hash of a file's content"""
    md5 = hashlib.md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            md5.update(chunk)
    return md5.hexdigest()


def test_get_threshold():
    """Test estimation of filtering threshold and figure generation."""
    fig_file = "test_event_dist.png"
    thr_uncut, thr_loop = hcf.get_thresholds(
        "test_data/valid_idx.pairs", plot_events=True, fig_path=fig_file
    )
    # Check if the estimated threshold are identical to expected values.
    assert thr_uncut == 6
    assert thr_loop == 5

    # Monkey-patch input to simulate stdin when testing interactive mode
    input_values = [6, 5]
    hcf.input = lambda x: input_values.pop(0)
    thr_uncut, thr_loop = hcf.get_thresholds(
        "test_data/valid_idx.pairs",
        interactive=True,
        plot_events=True,
        fig_path=fig_file,
    )
    # Remove figure if it was created
    if os.path.isfile(fig_file):
        os.remove(fig_file)


def test_filter_pairs():
    """Test pairs file filtering function."""
    filt_pairs = NamedTemporaryFile("w", delete=False)
    filt_pairs.close()
    fig_file = "test_piechart.png"

    hcf.filter_events(
        "test_data/valid_idx.pairs",
        filt_pairs.name,
        6,
        5,
        plot_events=True,
        fig_path=fig_file,
    )

    # Test if the filtered pairs file mathes expectations
    assert hash_file("test_data/valid_idx_filtered.pairs") == hash_file(filt_pairs.name)

    # Remove figure if it was created
    if os.path.isfile(fig_file):
        os.remove(fig_file)
    os.unlink(filt_pairs.name)


def read_pairs(path):
    """Header lines and data lines of a pairs file."""
    with open(path) as pairs:
        lines = pairs.readlines()
    return [x for x in lines if x.startswith("#")], [x for x in lines if not x.startswith("#")]


def write_pairs(path, header, lines):
    with open(path, "w") as pairs:
        pairs.writelines(header + list(lines))
    return str(path)


def by_coord(line, chrom_order):
    """Sort key of a pairs line: chromosome (in the given order), then position."""
    fields = line.split("\t")
    return chrom_order.index(fields[1]), int(fields[2])


def test_sample_pairs_whole_small_file():
    """A file of at most the sample size is used whole, as before."""
    _, lines = read_pairs(PAIRS)
    assert list(hcf._sample_pairs(PAIRS, len(lines))) == lines
    assert list(hcf._sample_pairs(PAIRS, 1_000_000)) == lines


def test_sample_pairs_ignore_order(tmp_path):
    """The same pairs are drawn whatever the order of the file, from all of it."""
    header, lines = read_pairs(PAIRS)
    rng = np.random.default_rng(0)
    orders = {
        "seq1_first": sorted(lines, key=lambda x: by_coord(x, ["seq1", "seq2"])),
        "seq2_first": sorted(lines, key=lambda x: by_coord(x, ["seq2", "seq1"])),
        "shuffled": [lines[i] for i in rng.permutation(len(lines))],
    }
    samples = {
        name: sorted(hcf._sample_pairs(write_pairs(tmp_path / name, header, ordered), 1000))
        for name, ordered in orders.items()
    }
    sample = samples["shuffled"]
    assert samples["seq1_first"] == sample
    assert samples["seq2_first"] == sample
    # About 1000 of the 10000 pairs (standard deviation 30)...
    assert 900 <= len(sample) <= 1100

    # ...spread over both chromosomes as the whole library is, whereas the
    # first 1000 pairs of a sorted file are all on its first chromosome.
    def seq2_share(pairs):
        return np.mean([x.split("\t")[1] == "seq2" for x in pairs])

    assert abs(seq2_share(sample) - seq2_share(lines)) < 0.05


def event_library(rng):
    """Pairs on two chromosomes, a and b. All four orientations are equally
    frequent at every distance (40 pairs per orientation and number of
    restriction sites), but uncuts (+-) and loops (-+) are 10 times as
    frequent up to 7 sites apart on a, and only up to 2 sites apart on b."""
    events = []
    for chrom, reach in (("a", 7), ("b", 2)):
        for nsites in range(50):
            for strands in ("++", "--", "+-", "-+"):
                count = 400 if strands in ("+-", "-+") and nsites <= reach else 40
                for frag in rng.integers(0, 10_000, size=count):
                    events.append((chrom, int(frag), strands, nsites))
    names = rng.permutation(len(events))
    return [
        f"read{name}\t{chrom}\t{frag * 100 + 10}\t{chrom}\t{(frag + nsites) * 100 + 60}"
        f"\t{strands[0]}\t{strands[1]}\t{frag}\t{frag + nsites}\n"
        for name, (chrom, frag, strands, nsites) in zip(names, events)
    ]


def test_get_thresholds_ignore_order(tmp_path):
    """Thresholds estimated from a sample do not depend on the order of the
    pairs, and reflect the whole library."""
    header = [
        "## pairs format v1.0\n",
        "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 frag1 frag2\n",
    ]
    rng = np.random.default_rng(1)
    lines = event_library(rng)
    orders = {
        "a_first": sorted(lines, key=lambda x: by_coord(x, ["a", "b"])),
        "b_first": sorted(lines, key=lambda x: by_coord(x, ["b", "a"])),
        "shuffled": [lines[i] for i in rng.permutation(len(lines))],
    }
    paths = {
        name: write_pairs(tmp_path / name, header, ordered) for name, ordered in orders.items()
    }
    thresholds = {name: hcf.get_thresholds(path, sample_size=5000) for name, path in paths.items()}
    assert thresholds["a_first"] == thresholds["shuffled"]
    assert thresholds["b_first"] == thresholds["shuffled"]
    # The uncuts and loops of chromosome a are seen even when b comes first.
    # The first 5000 pairs of that file, all on b, would give lower thresholds.
    uncut, loop = thresholds["b_first"]
    assert uncut > 7 and loop > 7
    head = write_pairs(tmp_path / "head", header, orders["b_first"][:5000])
    head_uncut, head_loop = hcf.get_thresholds(head)
    assert head_uncut < uncut and head_loop < loop
