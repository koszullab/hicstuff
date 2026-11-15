# test input/output functions from hicstuff
# 20190402
import bz2
import filecmp
import gzip
import os
import random
import re
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


def test_sort_compressed_pairs_header_only(tmp_path):
    """A compressed pairs file without any pair gives a header-only output."""
    plain = tmp_path / "empty.pairs"
    plain.write_text("## pairs format v1.0\n#sorted: readID\n#columns: readID chr1 pos1\n")
    out = tmp_path / "sorted.pairs"
    hio.sort_pairs(str(_compress(plain, "gzip")), str(out), keys=["chr1", "pos1"])
    assert (
        out.read_text() == "## pairs format v1.0\n#sorted: chr1-pos1\n#columns: readID chr1 pos1\n"
    )


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
