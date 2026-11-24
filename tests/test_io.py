# test input/output functions from hicstuff
# 20190402
import bz2
import filecmp
import gzip
import logging
import os
import random
import re
import shlex
import shutil
import subprocess
import zipfile
from pathlib import Path
from tempfile import NamedTemporaryFile

import cooler
import numpy as np
import pandas as pd
import pytest
from Bio import SeqIO

import hicstuff.io as hio

GENOME = "test_data/genome/seq.fa"
MAT_GRAAL = hio.load_sparse_matrix("test_data/abs_fragments_contacts_weighted.txt")
FRAGS_GRAAL = pd.read_csv("test_data/fragments_list.txt", delimiter="\t")
PAIRS_FILE = "test_data/valid_idx_filtered.pairs"
PAIRS_SORT_KEYS = ["chr1", "chr2", "pos1", "pos2"]


def test_compress():
    """Test reading and checking of compressed files"""

    # Generate temp files and store path
    f = NamedTemporaryFile(delete=False, mode="w")
    fgz = NamedTemporaryFile(delete=False, mode="wb")
    fbz = NamedTemporaryFile(delete=False)
    fz = NamedTemporaryFile(delete=False)
    fgz.close()
    fbz.close()
    fz.close()
    # Fill with some text
    f.write("xyz")
    f.close()
    # Write text to compressed file using different compression types
    raw = open(f.name, mode="rb").read()
    gz = gzip.open(fgz.name, mode="wb")
    gz.write(raw)
    gz.close()
    bz = bz2.BZ2File(fbz.name, "wb")
    bz.write(raw)
    bz.close()
    z = zipfile.ZipFile(fz.name, "w")
    z.write(f.name)
    z.close()
    exp_compress = {
        f.name: False,
        fgz.name: True,
        fz.name: True,
        fbz.name: True,
    }
    for fh in [f, fgz, fbz, fz]:
        content = hio.read_compressed(fh.name).read()
        # Check reading
        assert content == "xyz"
        # Check guessing compression state
        assert hio.is_compressed(fh.name) == exp_compress[fh.name]
        # Clean files
        os.unlink(fh.name)


@pytest.mark.parametrize("compression", ["none", "gzip", "bz2", "zip"])
def test_read_compressed_binary(tmp_path, compression):
    """read_compressed must return the decompressed bytes in "rb" mode."""
    plain = tmp_path / "content.txt"
    plain.write_bytes(b"xyz\n")
    path = plain if compression == "none" else _compress(plain, compression)
    with hio.read_compressed(str(path), "rb") as f:
        assert f.read() == b"xyz\n"


def _compress(plain, compression):
    """Compress a file with gzip, bgzip, bz2 or zip and return the new path."""
    if compression in ("gzip", "bgzip"):
        out = plain.with_name(plain.name + ".gz")
        if compression == "bgzip":
            with open(out, "wb") as f:
                subprocess.run(["bgzip", "-c", str(plain)], stdout=f, check=True)
        else:
            with gzip.open(out, "wb") as f:
                f.write(plain.read_bytes())
    elif compression == "bz2":
        out = plain.with_name(plain.name + ".bz2")
        with bz2.open(out, "wb") as f:
            f.write(plain.read_bytes())
    else:
        out = plain.with_name(plain.name + ".zip")
        with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
            z.write(plain, arcname=plain.name)
    return out


def _write_shuffled_pairs(path, n_copies=1):
    """Write the test pairs, data lines shuffled, and return (header, data)."""
    with open(PAIRS_FILE) as f:
        lines = f.read().splitlines()
    header = [line for line in lines if line.startswith("#")]
    data = [line for line in lines if not line.startswith("#")] * n_copies
    random.Random(42).shuffle(data)
    path.write_text("\n".join(header + data) + "\n")
    return header, data


@pytest.mark.parametrize("compression", ["gzip", "bgzip", "bz2", "zip"])
def test_sort_compressed_pairs(tmp_path, compression):
    """Sorting compressed pairs must give the same output as plain pairs."""
    if compression == "bgzip" and shutil.which("bgzip") is None:
        pytest.skip("bgzip is not installed")
    plain = tmp_path / "unsorted.pairs"
    header, data = _write_shuffled_pairs(plain)
    compressed = _compress(plain, compression)
    assert hio.is_compressed(str(compressed))
    assert hio.get_pairs_header(str(compressed)) == header

    sorted_plain = tmp_path / "sorted_plain.pairs"
    sorted_comp = tmp_path / "sorted_compressed.pairs"
    hio.sort_pairs(str(plain), str(sorted_plain), keys=PAIRS_SORT_KEYS)
    hio.sort_pairs(str(compressed), str(sorted_comp), keys=PAIRS_SORT_KEYS)

    # The output is plain text, identical to the one from the plain input
    assert sorted_comp.read_bytes() == sorted_plain.read_bytes()
    lines = sorted_comp.read_text().splitlines()
    # The header is kept, with the new sort order
    assert "#sorted: readID" in header
    expected_header = [
        "#sorted: chr1-chr2-pos1-pos2" if line.startswith("#sorted") else line for line in header
    ]
    assert lines[: len(header)] == expected_header
    sorted_data = lines[len(header) :]
    assert sorted(sorted_data) == sorted(data)
    keys = [(f[1], f[3], int(f[2]), int(f[4])) for f in (line.split("\t") for line in sorted_data)]
    assert keys == sorted(keys)


def _sort_pairs_lines(tmp_path, lines, keys):
    """Sort shuffled pairs lines with sort_pairs and return the sorted rows."""
    lines = list(lines)
    random.Random(0).shuffle(lines)
    unsorted, out = tmp_path / "unsorted.pairs", tmp_path / "sorted.pairs"
    header = "## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n"
    unsorted.write_text(header + "\n".join(lines) + "\n")
    hio.sort_pairs(str(unsorted), str(out), keys=keys)
    rows = [line.split("\t") for line in out.read_text().splitlines() if not line.startswith("#")]
    assert sorted("\t".join(row) for row in rows) == sorted(lines)
    return rows


def _blocks(values):
    """The values in order, with each run of equal values counted once."""
    return [v for i, v in enumerate(values) if i == 0 or v != values[i - 1]]


#: Contig names whose byte and version orders differ (contig_2 and contig_10),
#: and names that sort -V (coreutils >= 9.1) finds equal although they differ:
#: contig_01 and contig_1, a and a0, k141_10 and k141_010.
TIED_CONTIGS = ["contig_01", "contig_1", "contig_2", "contig_10", "a", "a0", "k141_10", "k141_010"]


@pytest.mark.parametrize(
    "keys", [["chr1", "chr2", "pos1", "pos2"], ["chr1", "pos1", "chr2", "pos2"]], ids="-".join
)
def test_sort_pairs_one_block_per_contig(tmp_path, keys):
    """Each contig forms one block, in byte order, also when version order would tie it."""
    rng = random.Random(1)
    lines = [
        f"r{i}\t{chr1}\t{rng.randint(1, 50)}\t{chr2}\t{rng.randint(1, 50)}\t+\t-"
        for i, (chr1, chr2) in enumerate(
            (chr1, chr2) for chr1 in TIED_CONTIGS for chr2 in TIED_CONTIGS for _ in range(4)
        )
    ]
    rows = _sort_pairs_lines(tmp_path, lines, keys)

    # One block per contig, in byte order (contig_10 before contig_2), as
    # pairtools sort
    assert _blocks([row[1] for row in rows]) == sorted(TIED_CONTIGS)
    if keys[1] == "chr2":
        # One block per chromosome pair, sorted by positions: what pairix
        # needs to index the file
        groups = [(row[1], row[3]) for row in rows]
        within = [(int(row[2]), int(row[4])) for row in rows]
    else:
        groups = [row[1] for row in rows]
        within = [int(row[2]) for row in rows]
    assert _blocks(groups) == sorted(set(groups))
    for group in set(groups):
        values = [value for g, value in zip(groups, within) if g == group]
        assert values == sorted(values)


def test_sort_pairs_one_block_per_read(tmp_path):
    """Each read forms one block, also when dictionary order ties it with another."""
    # Dictionary order ignores punctuation: r.1, r:1 and r1 compare equal, and
    # so do these two Illumina read names (tile 1101, x 111 and 11)
    names = ["r.1", "r1", "r:1", "A1:1101:111:1000", "A1:1101:11:11000"]
    lines = [f"{name}\t{chrom}\t1\t{chrom}\t2\t+\t+" for name in names for chrom in "ABC"]
    rows = _sort_pairs_lines(tmp_path, lines, keys=["readID", "chr1"])
    # Tied names come in byte order
    assert _blocks([row[0] for row in rows]) == [
        "A1:1101:111:1000",
        "A1:1101:11:11000",
        "r.1",
        "r1",
        "r:1",
    ]


@pytest.mark.needs("sort")
def test_sort_pairs_runs_sort_in_c_locale(tmp_path, monkeypatch):
    """sort runs in the C locale whatever the caller's, so ties are broken byte by byte."""
    monkeypatch.setenv("LC_ALL", "en_US.UTF-8")
    real_popen, envs = hio.sp.Popen, []

    def popen(cmd, *args, **kwargs):
        if cmd[0] == "sort" and "--version" not in cmd:
            envs.append(kwargs.get("env"))
        return real_popen(cmd, *args, **kwargs)

    monkeypatch.setattr(hio.sp, "Popen", popen)
    _sort_pairs_lines(tmp_path, ["r1\tA\t1\tA\t2\t+\t+"], keys=["readID"])
    assert [env and env.get("LC_ALL") for env in envs] == ["C"]


def test_sort_pairs_duplicates_adjacent(tmp_path):
    """Pairs with the same coordinates are adjacent, as filter_pcr_dup needs."""
    lines = [
        f"r{i}\tA\t10\t{chr2}\t20\t+\t+" for i, chr2 in enumerate(["contig_01", "contig_1"] * 3)
    ]
    rows = _sort_pairs_lines(tmp_path, lines, keys=["chr1", "pos1", "chr2", "pos2"])
    assert _blocks([row[3] for row in rows]) == ["contig_01", "contig_1"]


@pytest.mark.parametrize("key", ["strand1", "strand2"])
def test_sort_pairs_by_strand(tmp_path, key):
    """Strands sort: dictionary order alone ignores + and -."""
    lines = [f"r{i}\tA\t{i}\tA\t{i}\t{'+-'[i % 2]}\t{'-+'[i % 3 % 2]}" for i in range(12)]
    rows = _sort_pairs_lines(tmp_path, lines, keys=[key])
    strands = [row[5 if key == "strand1" else 6] for row in rows]
    assert strands == sorted(strands)
    assert _blocks(strands) == ["+", "-"]


def test_sort_pairs_unknown_key(tmp_path):
    """An unknown sort key is refused before the output is created."""
    out = tmp_path / "sorted.pairs"
    with pytest.raises(ValueError, match="chrom1"):
        hio.sort_pairs(PAIRS_FILE, str(out), keys=["chrom1", "pos1"])
    assert not out.exists()


#: What sort --version prints on stdout, how sort_pairs identifies it, and the
#: warning it gives when it does not pass --parallel. GNU 9.4 is Ubuntu 24.04's
#: sort, captured; the uutils line is that of Ubuntu 26.04's /usr/bin/sort
#: (rust-coreutils 0.10.0); the BSD lines are VERSION in FreeBSD's and Apple's
#: sort.h; busybox sort has no --version and only prints its usage on stderr.
SORT_VERSIONS = {
    "GNU-9.4": (
        "sort (GNU coreutils) 9.4\n"
        "Copyright (C) 2023 Free Software Foundation, Inc.\n"
        "License GPLv3+: GNU GPL version 3 or later <https://gnu.org/licenses/gpl.html>.\n"
        "This is free software: you are free to change and redistribute it.\n"
        "There is NO WARRANTY, to the extent permitted by law.\n"
        "\n"
        "Written by Mike Haertel and Paul Eggert.\n",
        ("GNU", (9, 4)),
        None,
    ),
    "GNU-8.23": (
        "sort (GNU coreutils) 8.23\nCopyright (C) 2014 Free Software Foundation, Inc.\n",
        ("GNU", (8, 23)),
        None,
    ),
    "GNU-8.22": (
        "sort (GNU coreutils) 8.22\nCopyright (C) 2013 Free Software Foundation, Inc.\n",
        ("GNU", (8, 22)),
        "GNU sort version is 8.22 but 8.23 or later is required",
    ),
    "GNU-8.4": (
        "sort (GNU coreutils) 8.4\nCopyright (C) 2010 Free Software Foundation, Inc.\n",
        ("GNU", (8, 4)),
        "GNU sort version is 8.4 but 8.23 or later is required",
    ),
    "uutils": ("sort (uutils coreutils) 0.10.0\n", ("uutils", (0, 10)), None),
    "FreeBSD": ("2.3-FreeBSD\n", ("BSD", None), "Using BSD sort"),
    "macOS": ("2.3-Apple (199)\n", ("BSD", None), "Using BSD sort"),
    "busybox": (
        "",
        ("unknown", None),
        "Unknown sort (sort --version printed ''): sorting pairs in Python",
    ),
}


@pytest.mark.parametrize(
    ("text", "implementation"), [v[:2] for v in SORT_VERSIONS.values()], ids=list(SORT_VERSIONS)
)
def test_sort_implementation(text, implementation):
    """sort --version tells GNU, uutils and BSD sort apart, with their version."""
    assert hio._sort_implementation(text) == implementation


@pytest.mark.parametrize(
    ("text", "warning"), [(v[0], v[2]) for v in SORT_VERSIONS.values()], ids=list(SORT_VERSIONS)
)
def test_sort_pairs_parallel(tmp_path, monkeypatch, caplog, capfd, text, warning):
    """--parallel goes to GNU sort from 8.23 and to uutils sort, as Ubuntu 26.04's.

    A sort that prints text for --version stands first on PATH, and runs the
    real sort. An unknown sort, as busybox's, is not run: the pairs are
    sorted in Python.
    """
    if hio._sort_implementation(text)[0] != "unknown" and shutil.which("sort") is None:
        pytest.skip("not on PATH: sort")
    fake = tmp_path / "bin" / "sort"
    fake.parent.mkdir()
    fake.with_suffix(".version").write_text(text)
    fake.write_text(
        "#!/bin/sh\n"
        'if [ "$1" = --version ]; then\n'
        '    cat "$0.version"\n'
        '    [ -s "$0.version" ] || { echo "sort: unrecognized option --version" >&2; exit 1; }\n'
        "    exit 0\n"
        "fi\n"
        'echo "$@" >> "$0.args"\n'
        f'exec {shlex.quote(shutil.which("sort"))} "$@"\n'
    )
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", f"{fake.parent}{os.pathsep}{os.environ['PATH']}")
    out = tmp_path / "sorted.pairs"
    hio.sort_pairs(PAIRS_FILE, str(out), keys=PAIRS_SORT_KEYS, threads=2)

    unknown = hio._sort_implementation(text)[0] == "unknown"
    args = [] if unknown else fake.with_suffix(".args").read_text().split()
    assert fake.with_suffix(".args").exists() is not unknown
    assert ("--parallel=2" in args) is (warning is None)
    if warning is None:
        assert not [r for r in caplog.records if r.levelno >= logging.WARNING]
    else:
        assert warning in caplog.text
    # busybox's usage is not shown
    assert "unrecognized option" not in capfd.readouterr().err
    # The pairs are all there
    rows = [line for line in out.read_text().splitlines() if not line.startswith("#")]
    pairs = Path(PAIRS_FILE).read_text().splitlines()
    assert sorted(rows) == sorted(line for line in pairs if not line.startswith("#"))


#: Pairs data lines that sort orders with care: blanks inside, before,
#: between and after fields, missing and extra fields, signs, decimals,
#: leading zeros and other text where numbers go, punctuation, both cases,
#: bytes above 127, a carriage return, and duplicates.
_NAMES = ["r.1", "r1", "r:1", "A1:1101:111:1000", "A1:1101:11:11000", "r1/1", "R1", "r 1", ""]
_CHROMS = TIED_CONTIGS + ["chr1", "chr10", "chr2", "A", "chré", "chrX ", "c-1", "c.1", "c1", ""]
_NUMBERS = ["0", "1", "10", "2", "007", "7", "-1", "-0", "-", "1.5", "1.50", ".5", "5.", "NA"]
_NUMBERS += ["1e5", "+5", " 12", "12abc", "", "99999999999999999999999", "-.5", "3\r"]
_STRANDS = ["+", "-", ".", "", "*", "+-", "- "]


def _tricky_pairs_lines(n_lines, seed):
    rng = random.Random(seed)
    lines = []
    for _ in range(n_lines):
        fields = [rng.choice(_NAMES), rng.choice(_CHROMS), rng.choice(_NUMBERS)]
        fields += [rng.choice(_CHROMS), rng.choice(_NUMBERS), rng.choice(_STRANDS)]
        fields += [rng.choice(_STRANDS), rng.choice(_NUMBERS), rng.choice(_NUMBERS), "extra"]
        fields = fields[: rng.choice([9] * 8 + [3, 5, 7, 10])]
        line = fields[0]
        for field in fields[1:]:
            line += rng.choice(["\t"] * 20 + ["\t\t", " ", "\t "]) + field
        line = rng.choice([""] * 18 + ["\t"]) + line + rng.choice([""] * 18 + ["\t", "\r"])
        lines.append(line)
    return lines + lines[: n_lines // 10]


def _gnu_sort():
    """Whether the sort on PATH is GNU sort."""
    if shutil.which("sort") is None:
        return False
    version = subprocess.run(["sort", "--version"], capture_output=True, text=True).stdout
    return hio._sort_implementation(version)[0] == "GNU"


@pytest.mark.parametrize(
    "keys",
    [[column] for column in ("readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2")]
    + [[column] for column in ("frag1", "frag2")]
    + [PAIRS_SORT_KEYS, ["frag1", "frag2"], ["chr1", "pos1", "chr2", "pos2", "strand1", "strand2"]]
    + [["readID", "pos1", "strand2", "frag2"]],
    ids="-".join,
)
def test_python_sort_orders_as_gnu_sort(keys, tmp_path, hide_tools):
    """Without sort, the pairs are sorted in Python, in the order of GNU sort,
    in one chunk or merging many temporary files, which are removed."""
    if not _gnu_sort():
        pytest.skip("GNU sort is needed to compare with")
    unsorted = tmp_path / "unsorted.pairs"
    header = "## pairs format v1.0\n#sorted: readID\n#columns: readID chr1 pos1 chr2 pos2\n"
    unsorted.write_text(header + "\n".join(_tricky_pairs_lines(500, seed=len(keys))))
    expected = tmp_path / "gnu.pairs"
    hio.sort_pairs(str(unsorted), str(expected), keys=keys)
    assert expected.read_bytes().endswith(b"\n")
    hide_tools("sort")
    tmp_dir = tmp_path / "tmp"
    tmp_dir.mkdir()
    # 2 lines in each temporary file with 1K, merged in two rounds
    for buffer in ("2G", "1K"):
        out = tmp_path / f"python_{buffer}.pairs"
        hio.sort_pairs(str(unsorted), str(out), keys=keys, tmp_dir=str(tmp_dir), buffer=buffer)
        assert out.read_bytes() == expected.read_bytes()
        assert not list(tmp_dir.iterdir())


#: Pairs data lines in the order that LC_ALL=C GNU sort gives them for keys of
#: sort_pairs: numbers with signs, decimals, leading zeros or text; names that
#: dictionary order ties; blanks around and in fields, and missing fields;
#: strands that dictionary order ignores. Recorded with GNU sort 9.4, so that
#: no sort program is needed to check them.
GNU_ORDERS = {
    "pos1": [
        "r2\tchr1\t-1",
        "r10\tchr1\t-.5",
        "r11\tchr1\tNA",
        "r12\tchr1\t",
        "r14\tchr1\t+5",
        "r18\tchr1\t-",
        "r3\tchr1\t-0",
        "r4\tchr1\t0",
        "r9\tchr1\t.5",
        "r13\tchr1\t1e5",
        "r7\tchr1\t1.5",
        "r8\tchr1\t1.50",
        "r17\tchr1\t5.",
        "r5\tchr1\t007",
        "r6\tchr1\t7",
        "r1\tchr1\t9",
        "r0\tchr1\t10",
        "r15\tchr1\t 12",
        "r16\tchr1\t12abc",
        "r19\tchr1\t99999999999999999999",
    ],
    "readID": [
        "A1:1101:111:1000\tchr1\t5",
        "A1:1101:11:11000\tchr1\t5",
        "R1\tchr1\t5",
        "r b\tchr1\t5",
        "r.1\tchr1\t5",
        "r1\tchr1\t5",
        "r:1\tchr1\t5",
        "r_1\tchr1\t5",
        "r1/1\tchr1\t5",
        "r-b\tchr1\t5",
        "rb\tchr1\t5",
    ],
    "chr1-pos1": [
        "r7",
        "r2\t\tchr1\t4",
        "r9\t chr1\t7",
        "r4\tchr1\t2\t",
        "\tr3\tchr1\t3",
        "r1\tchr1\t5",
        "r11\tchr1\t5",
        "r8\tchr1 \t6",
        "r10\tchr10\t8",
        "r6\tchr2",
        "r5 x\tchr1\t1",
    ],
    "strand1-strand2": [
        "r4\tA\t1\tA\t2\t*\t*",
        "r0\tA\t1\tA\t2\t+\t-",
        "r5\tA\t1\tA\t2\t+-\t+",
        "r1\tA\t1\tA\t2\t-\t+",
        "r2\tA\t1\tA\t2\t.\t+",
        "r7\tA\t1\tA\t2\t~\ta",
        "r3\tA\t1\tA\t2\t\t-",
        "r6\tA\t1\tA\t2\ta\t~",
    ],
}


@pytest.mark.parametrize("keys", list(GNU_ORDERS))
def test_python_sort_orders_lines_as_gnu_sort_did(keys, tmp_path, hide_tools):
    """Without sort, the pairs are sorted in Python, in the order of GNU sort."""
    expected = [line.encode() for line in GNU_ORDERS[keys]]
    lines = list(expected)
    random.Random(0).shuffle(lines)
    unsorted, out = tmp_path / "unsorted.pairs", tmp_path / "sorted.pairs"
    unsorted.write_bytes(b"## pairs format v1.0\n#sorted: readID\n" + b"\n".join(lines) + b"\n")
    hide_tools("sort")
    hio.sort_pairs(str(unsorted), str(out), keys=keys.split("-"))
    assert out.read_bytes().split(b"\n")[2:-1] == expected


def test_python_sort_merges_temporary_files(tmp_path, hide_tools):
    """Pairs sorted in chunks of two lines, then merged from hundreds of
    temporary files in two rounds, come in the order of one chunk, and no
    temporary file is left. The last chunk holds one line, without newline."""
    unsorted = tmp_path / "unsorted.pairs"
    _write_shuffled_pairs(unsorted, n_copies=3)
    unsorted.write_bytes(unsorted.read_bytes().rstrip(b"\n").rsplit(b"\n", 1)[0])
    keys = ["chr1", "pos1", "chr2", "pos2", "strand1", "strand2"]
    hide_tools("sort")
    tmp_dir = tmp_path / "tmp"
    tmp_dir.mkdir()
    one, many = tmp_path / "one.pairs", tmp_path / "many.pairs"
    hio.sort_pairs(str(unsorted), str(one), keys=keys, tmp_dir=str(tmp_dir))
    hio.sort_pairs(str(unsorted), str(many), keys=keys, tmp_dir=str(tmp_dir), buffer="1K")
    assert many.read_bytes() == one.read_bytes()
    assert not list(tmp_dir.iterdir())


def test_sort_env_sorts_in_python(tmp_path, monkeypatch, caplog):
    """HICSTUFF_SORT=python sorts in Python although a sort is on PATH."""
    caplog.set_level(logging.INFO)
    unsorted = tmp_path / "unsorted.pairs"
    _write_shuffled_pairs(unsorted)
    expected = tmp_path / "expected.pairs"
    hio.sort_pairs(str(unsorted), str(expected), keys=PAIRS_SORT_KEYS)
    # A sort that fails, first on PATH
    fake = tmp_path / "bin" / "sort"
    fake.parent.mkdir()
    fake.write_text("#!/bin/sh\nexit 2\n")
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", f"{fake.parent}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("HICSTUFF_SORT", "python")
    out = tmp_path / "sorted.pairs"
    hio.sort_pairs(str(unsorted), str(out), keys=PAIRS_SORT_KEYS)
    assert out.read_bytes() == expected.read_bytes()
    assert "Sorting pairs in Python (HICSTUFF_SORT=python)." in caplog.text


def test_sort_that_does_not_run_sorts_in_python(tmp_path, monkeypatch, caplog):
    """A sort on PATH that cannot be run is not used: the pairs are sorted in Python."""
    caplog.set_level(logging.INFO)
    unsorted = tmp_path / "unsorted.pairs"
    _write_shuffled_pairs(unsorted)
    expected = tmp_path / "expected.pairs"
    hio.sort_pairs(str(unsorted), str(expected), keys=PAIRS_SORT_KEYS)
    # Its interpreter does not exist, and it is the only sort on PATH: a run
    # looks for another one further on PATH when it cannot run the first
    fake = tmp_path / "bin" / "sort"
    fake.parent.mkdir()
    fake.write_text("#!/nonexistent/interpreter\n")
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", str(fake.parent))
    out = tmp_path / "sorted.pairs"
    hio.sort_pairs(str(unsorted), str(out), keys=PAIRS_SORT_KEYS)
    assert out.read_bytes() == expected.read_bytes()
    assert "Unknown sort (sort --version printed ''): sorting pairs in Python." in caplog.text


def test_python_sort_removes_temporary_files_on_error(tmp_path, hide_tools):
    """A truncated input raises, and leaves no temporary file behind."""
    plain = tmp_path / "unsorted.pairs"
    _write_shuffled_pairs(plain, n_copies=20)
    raw = _compress(plain, "gzip").read_bytes()
    truncated = tmp_path / "truncated.pairs.gz"
    truncated.write_bytes(raw[: len(raw) // 2])
    tmp_dir = tmp_path / "tmp"
    tmp_dir.mkdir()
    hide_tools("sort")
    with pytest.raises(EOFError):
        hio.sort_pairs(
            str(truncated),
            str(tmp_path / "out.pairs"),
            keys=PAIRS_SORT_KEYS,
            tmp_dir=str(tmp_dir),
            buffer="10K",
        )
    assert not list(tmp_dir.iterdir())


@pytest.mark.parametrize(
    ("buffer", "size"),
    [("2G", 2 << 30), ("500M", 500 << 20), ("64k", 64 << 10), ("1024", 1 << 20), ("300b", 300)],
)
def test_sort_buffer_bytes(buffer, size):
    """Sizes are read as sort -S reads them, in kibibytes by default."""
    assert hio._sort_buffer_bytes(buffer) == size
    assert 0 < hio._sort_buffer_bytes("1%") < hio._sort_buffer_bytes("2%")


def test_python_sort_refuses_invalid_buffer(tmp_path, hide_tools):
    """An invalid buffer size is refused before the output is created."""
    hide_tools("sort")
    out = tmp_path / "sorted.pairs"
    with pytest.raises(ValueError, match="Invalid sort buffer size 'not_a_size'"):
        hio.sort_pairs(PAIRS_FILE, str(out), keys=PAIRS_SORT_KEYS, buffer="not_a_size")
    assert not out.exists()


def test_sort_compressed_pairs_header_only(tmp_path):
    """A compressed pairs file without any pair gives a header-only output."""
    plain = tmp_path / "empty.pairs"
    plain.write_text("## pairs format v1.0\n#sorted: readID\n#columns: readID chr1 pos1\n")
    out = tmp_path / "sorted.pairs"
    hio.sort_pairs(str(_compress(plain, "gzip")), str(out), keys=["chr1", "pos1"])
    assert (
        out.read_text() == "## pairs format v1.0\n#sorted: chr1-pos1\n#columns: readID chr1 pos1\n"
    )


@pytest.mark.needs("sort")
def test_sort_compressed_pairs_failure_raises(tmp_path):
    """A failing UNIX sort must raise, also when it stops reading early."""
    plain = tmp_path / "unsorted.pairs"
    # Larger than the pipe buffers, so that sort exits while pairs are still sent
    _write_shuffled_pairs(plain, n_copies=100)
    assert plain.stat().st_size > 2 * 1024 * 1024
    with pytest.raises(RuntimeError, match="Sorting pairs"):
        hio.sort_pairs(
            str(_compress(plain, "gzip")),
            str(tmp_path / "out.pairs"),
            keys=PAIRS_SORT_KEYS,
            buffer="not_a_size",
        )


def test_sort_truncated_pairs_raises(tmp_path):
    """A truncated compressed file must raise instead of being sorted partially."""
    plain = tmp_path / "unsorted.pairs"
    _write_shuffled_pairs(plain, n_copies=100)
    compressed = _compress(plain, "gzip")
    truncated = tmp_path / "truncated.pairs.gz"
    raw = compressed.read_bytes()
    truncated.write_bytes(raw[: len(raw) // 2])
    with pytest.raises(EOFError):
        hio.sort_pairs(str(truncated), str(tmp_path / "out.pairs"), keys=PAIRS_SORT_KEYS)


def test_save_bedgraph2d():
    """Test saving 2D bedgraph files"""
    # Create temp file and write a 2D bedgraph matrix inside
    f = NamedTemporaryFile("w", delete=False)
    f.close()
    # Load GRAAL Matrix from test_data
    hio.save_bedgraph2d(MAT_GRAAL, FRAGS_GRAAL, f.name)
    # Check if the file created is identical to the 2D bedgraph matrix
    # containing the same dataset as the GRAAL matrix used as input
    assert filecmp.cmp(f.name, "test_data/mat.bg2")
    os.unlink(f.name)


def test_load_bedgraph2d():
    """Test loading sparse matrices from 2D bedgraph files"""
    mat_bg = hio.load_bedgraph2d(
        "test_data/mat.bg2", fragments_file="test_data/fragments_list.txt"
    )[0]
    assert np.allclose(MAT_GRAAL.todense(), mat_bg.todense())

    # Load using fixed bin sizes
    mat_bg = hio.load_bedgraph2d("test_data/mat_5kb.bg2", bin_size=5000)[0]
    mat_graal = hio.load_sparse_matrix("test_data/mat_5kb.tsv")
    assert mat_bg.shape == mat_graal.shape


def test_cooler_io():
    """Test input output operations on cool files"""
    f = NamedTemporaryFile("w", delete=False)
    f.close()
    # Write cool from GRAAL objects
    hio.save_cool(f.name, MAT_GRAAL, FRAGS_GRAAL)
    c = cooler.Cooler(f.name)
    # Just copy the start column and write it as a new one
    dummycol = c.bins()[:].loc[:, "start"]  # pylint: disable=no-member
    hio.add_cool_column(c, dummycol, "dummy", dtype=np.int64)
    # Read custom cool into GRAAL objects
    mat, frags, chroms = hio.load_cool(f.name)
    os.unlink(f.name)


def test_hic_format():
    assert hio.get_hic_format("test_data/mat_5kb.bg2") == "bg2"
    assert hio.get_hic_format("test_data/abs_fragments_contacts_weighted.txt") == "graal"
    with pytest.raises(ValueError):
        assert hio.get_hic_format("test_data/valid.pairs")


#: Files of an index after its prefix, as each aligner (and layout) writes them
INDEX_FILES = {
    "bowtie2": [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"],
    "bowtie2-large": [".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"],
    "bwa": [".amb", ".ann", ".bwt", ".pac", ".sa"],
    "bwa-64": [".64.amb", ".64.ann", ".64.bwt", ".64.pac", ".64.sa"],
}


def _touch_index(prefix, suffixes):
    for suffix in suffixes:
        Path(f"{prefix}{suffix}").touch()


@pytest.mark.parametrize("layout", INDEX_FILES)
@pytest.mark.parametrize("name", ["genome", "genome.v2", "genome.v2.fa"])
def test_check_fasta_index(tmp_path, monkeypatch, layout, name):
    """An index is found from its prefix written in any form, and that prefix
    is returned as given: relative with ./ or ../, absolute, dots in names."""
    mode = layout.split("-")[0]
    other_mode = {"bowtie2": "bwa", "bwa": "bowtie2"}[mode]
    # A hidden directory, with a dot in its name, holding the index next to
    # the index of another version of the genome.
    index_dir = tmp_path / ".refs" / "v2.0"
    index_dir.mkdir(parents=True)
    _touch_index(index_dir / name, INDEX_FILES[layout])
    _touch_index(index_dir / "genome.v1.fa", INDEX_FILES[layout])
    (tmp_path / "work").mkdir()
    monkeypatch.chdir(tmp_path / "work")

    rel = f"../.refs/v2.0/{name}"
    for ref in (rel, f"./{rel}", str(index_dir / name)):
        assert hio.check_fasta_index(ref, mode=mode) == ref
        assert hio.check_fasta_index(ref, mode=other_mode) is None
    assert hio.check_fasta_index(Path(rel), mode=mode) == rel
    # Neither a partial prefix, nor a prefix with an extra extension, nor the
    # prefix of an index that does not exist.
    for ref in ("../.refs/v2.0/gen", f"{rel}.fa", "../.refs/v2.0/genome.v3.fa"):
        assert hio.check_fasta_index(ref, mode=mode) is None

    monkeypatch.chdir(index_dir)
    for ref in (name, f"./{name}"):
        assert hio.check_fasta_index(ref, mode=mode) == ref


@pytest.mark.parametrize("layout", INDEX_FILES)
@pytest.mark.parametrize(
    ("fasta", "prefixes"),
    [
        ("genome.fa", ["genome.fa", "genome"]),
        ("genome.v2.FASTA", ["genome.v2.FASTA", "genome.v2"]),
        ("genome.fna.gz", ["genome.fna.gz", "genome.fna", "genome"]),
        ("genome.gz", ["genome.gz", "genome"]),
    ],
)
def test_check_fasta_index_without_extension(tmp_path, monkeypatch, layout, fasta, prefixes):
    """The index of a FASTA file is found named after it with or without its
    extensions (bowtie2-build genome.fa genome), the FASTA path as given
    first, and then without its compression extension."""
    mode = layout.split("-")[0]
    ref_dir = tmp_path / "ref.v2"
    ref_dir.mkdir()
    monkeypatch.chdir(ref_dir)
    (ref_dir / fasta).write_text(">chr1\nACGT\n")
    # Each prefix, from the last one tried to the first, is found once its
    # index exists, and returned as written: relative with ./ or ../, or
    # absolute.
    for prefix in reversed(prefixes):
        _touch_index(ref_dir / prefix, INDEX_FILES[layout])
        for ref, root in ((fasta, ""), (f"./{fasta}", "./"), (f"../ref.v2/{fasta}", "../ref.v2/")):
            assert hio.check_fasta_index(ref, mode=mode) == root + prefix
        assert hio.check_fasta_index(ref_dir / fasta, mode=mode) == str(ref_dir / prefix)


def test_check_fasta_index_keeps_other_extensions(tmp_path):
    """Only FASTA and compression extensions are left out, from a file: an
    index prefix is never cut, nor the path of a missing FASTA file."""
    _touch_index(tmp_path / "genome", INDEX_FILES["bowtie2"])
    _touch_index(tmp_path / "genome.v2", INDEX_FILES["bowtie2"])
    (tmp_path / "genome.txt").write_text(">chr1\nACGT\n")
    (tmp_path / "genome.v2.fa").mkdir()
    for ref in ("genome.txt", "genome.fa", "genome.v2.fa", "genome.v2.1"):
        assert hio.check_fasta_index(tmp_path / ref, mode="bowtie2") is None


def test_check_fasta_index_incomplete(tmp_path):
    """An index missing a file is no index, nor are bowtie2 files of both sizes."""
    prefix = tmp_path / "genome.fa"
    _touch_index(prefix, INDEX_FILES["bowtie2"][:3] + INDEX_FILES["bowtie2-large"][3:])
    _touch_index(prefix, INDEX_FILES["bwa"][:-1] + INDEX_FILES["bwa-64"][-1:])
    for mode in ("bowtie2", "bwa"):
        assert hio.check_fasta_index(prefix, mode=mode) is None
    # The missing file completes the index
    _touch_index(prefix, INDEX_FILES["bwa"][-1:])
    assert hio.check_fasta_index(prefix, mode="bwa") == str(prefix)


@pytest.mark.parametrize(
    "ref", ["../ref/genome.v2.fa", "./genome.v2.fa.gz", "../genome.v2.mmi", "/data/genome.fasta"]
)
def test_check_fasta_index_minimap2(ref):
    """minimap2 aligns on the FASTA or .mmi file itself: its path is kept as given."""
    assert hio.check_fasta_index(ref, mode="minimap2") == ref


def test_renamed_fasta_next_to_input(tmp_path, monkeypatch):
    """rename_genome and reorder_fasta write by default <input without its
    extension>_renamed.fa, next to the input."""
    ref_dir = tmp_path / "ref.v2"
    ref_dir.mkdir()
    shutil.copy(GENOME, ref_dir / "genome.v2.fa")
    expected = [rec.id for rec in SeqIO.parse(GENOME, "fasta")]
    renamed = ref_dir / "genome.v2_renamed.fa"
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.chdir(work)

    hio.rename_genome("../ref.v2/genome.v2.fa")
    assert [rec.id for rec in SeqIO.parse(renamed, "fasta")] == expected
    renamed.unlink()
    hio.reorder_fasta(Path("./../ref.v2/genome.v2.fa"), threshold=0)
    assert sorted(rec.id for rec in SeqIO.parse(renamed, "fasta")) == sorted(expected)
    assert list(work.iterdir()) == []


def test_check_is_fasta(tmp_path):
    """FASTA input, plain or compressed, is recognised; anything else is not."""
    assert hio.check_is_fasta(GENOME)
    fasta_gz = tmp_path / "seq.fa.gz"
    with open(GENOME, "rb") as src, gzip.open(fasta_gz, "wb") as dst:
        shutil.copyfileobj(src, dst)
    assert hio.check_is_fasta(str(fasta_gz))
    # A text file that is not FASTA, a binary file, and a bowtie2 index prefix
    assert not hio.check_is_fasta(PAIRS_FILE)
    binary = tmp_path / "index.1.bt2"
    binary.write_bytes(bytes(range(256)) * 4)
    assert not hio.check_is_fasta(str(binary))
    assert not hio.check_is_fasta(str(tmp_path / "index"))


def test_gc_bins():
    """Test function to compute GC content per bin in genome"""
    bins = FRAGS_GRAAL
    gc = hio.gc_bins(GENOME, bins)
    assert gc.shape[0] == bins.shape[0]
    assert gc.min() >= 0.0
    assert gc.max() <= 1.0

    # Check values correctness on first chromosome
    first_record = next(SeqIO.parse(GENOME, "fasta"))
    seq = str(first_record.seq)
    starts = bins.start_pos[bins.chrom == first_record.id].values
    ends = bins.end_pos[bins.chrom == first_record.id].values
    for i, (start, end) in enumerate(zip(starts, ends)):
        bin_seq = seq[start:end]
        gc_count = len(re.findall(r"[GgCc]", bin_seq))
        gc_prop = gc_count / len(bin_seq)
        assert np.isclose(gc_prop, gc[i])


def test_get_pos_col():
    """Test function to guess chrom, start, end columns"""
    dummy = np.random.random((10, 5))  # pylint: disable=no-member
    chrom, start, end = "CHROMOSOME", "sTaRt", "EnD"
    df = pd.DataFrame(dummy, columns=["something", "different", start, end, chrom])
    cols = hio.get_pos_cols(df)
    assert len(cols) == 3
    assert cols == (chrom, start, end)

    bad_df = pd.DataFrame(dummy, columns=[el for el in "abcde"])
    with pytest.raises(ValueError):
        hio.get_pos_cols(bad_df)


def test_check_fastq_entries():
    """Test function to check number of entries in fastq files"""
    filen = "test_data/sample.reads_for.fastq.gz"
    nreads = hio.check_fastq_entries(filen)
    assert nreads == 10000
    nreads = hio.check_fastq_entries(filen, threads=2)
    assert nreads == 10000


@pytest.mark.parametrize(
    "name, stripped",
    [
        ("SRR1.7/1", "SRR1.7"),
        ("SRR1.7/2", "SRR1.7"),
        ("read", "read"),
        # Only a trailing /1 or /2 is a mate suffix
        ("read/3", "read/3"),
        ("read/12", "read/12"),
        ("read_1", "read_1"),
        ("read/1/2", "read/1"),
        ("read/1x", "read/1x"),
    ],
)
def test_strip_mate_suffix(name, stripped):
    assert hio.strip_mate_suffix(name) == stripped
