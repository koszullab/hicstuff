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


def test_check_fasta_index():
    f = NamedTemporaryFile("w", delete=False)
    f.close()
    assert hio.check_fasta_index(f.name, mode="minimap2") == f.name
    for i in range(6):
        Path(f.name + f".{i}.bt2").touch()
    assert hio.check_fasta_index(f.name, mode="bowtie2") == f.name
    assert hio.check_fasta_index(f.name, mode="bwa") is None
    os.unlink(f.name)
    for i in range(6):
        os.unlink(f.name + f".{i}.bt2")


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
