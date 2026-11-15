"""Iterative alignment writes every read exactly once.

bwa mem and minimap2 split a read that spans a ligation junction into a
primary record and a supplementary (0x800) one. Iterative alignment must
decide whether a read is aligned from its primary record alone: a read whose
primary record passes is written and never aligned again, whatever other
records it has.

The library below makes this happen before the last round. The 5' end of
each ``chimeric`` read is a 45 bp repeat (MAPQ 0), so it fails every
truncation up to 80 bp; at 100 or 120 bp, its unique 3' end becomes the
primary record, with a MAPQ that passes, and the repeat a supplementary
record.
"""

import collections
import subprocess as sp

import numpy as np
import pysam
import pytest

import hicstuff.iteralign as hci
import hicstuff.pipeline as hpi

CHROM_LEN = 10_000
READ_LEN = 120
#: Distance between a repeat and its copy on the same chromosome.
COPY_SHIFT = 5_000
#: Read pair classes, in read ID order, with the number of pairs of each.
CLASSES = {"chimeric": 8, "unique": 4, "repeat": 4, "unmapped": 4}
#: Length of the repeat at the 5' end of chimeric and repeat forward mates.
REPEAT_LEN = {"chimeric": 45, "repeat": 80}
#: Mapping quality threshold of the iterative alignment.
MIN_QUAL = 30

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def _write_fastq(path, reads):
    with open(path, "w") as fq:
        for name, seq in reads:
            fq.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


@pytest.fixture(scope="module")
def chimeric_lib(tmp_path_factory, require_tools):
    """Two random chromosomes and 120 bp read pairs of four classes:

    * ``chimeric``: both mates are a 45 bp repeat followed by 75 bp of the
      other chromosome (the forward mate on chr1 then chr2, the reverse mate
      the reverse complement of chr2 then chr1);
    * ``unique``: both mates read unique sequence;
    * ``repeat``: the forward mate is an 80 bp repeat followed by 40 bp of
      chr2, so its primary record is the repeat even at full length;
    * ``unmapped``: the forward mate is random sequence.

    Returns the genome (indexed for all three aligners), both FASTQ files,
    the info_contigs file and, per read ID, the class and the 0-based
    position the unique 3' end of a chimeric forward mate starts at.
    """
    root = tmp_path_factory.mktemp("chimeric")
    rng = np.random.default_rng(1984)
    chroms = {c: list(rng.choice(list("ACGT"), CHROM_LEN)) for c in ("chr1", "chr2")}
    kinds = [kind for kind, n in CLASSES.items() for _ in range(n)]
    slots = [300 + 220 * i for i in range(len(kinds))]
    # Copy the repeats first, so that the reads see the final sequence
    for kind, p in zip(kinds, slots):
        if kind in REPEAT_LEN:
            n = REPEAT_LEN[kind]
            for chrom in ("chr1", "chr2") if kind == "chimeric" else ("chr1",):
                seq = chroms[chrom]
                seq[p + COPY_SHIFT : p + COPY_SHIFT + n] = seq[p : p + n]
    chr1, chr2 = ("".join(chroms[c]) for c in ("chr1", "chr2"))
    reads1, reads2, truth = [], [], {}
    for i, (kind, p) in enumerate(zip(kinds, slots)):
        name = f"read{i:04d}"
        if kind in REPEAT_LEN:
            # Repeat, then the unique end, up to READ_LEN
            n, u = REPEAT_LEN[kind], p + 100
            read1 = chr1[p : p + n] + chr2[u : u + READ_LEN - n]
        elif kind == "unique":
            read1 = chr1[p : p + READ_LEN]
        else:
            read1 = "".join(rng.choice(list("ACGT"), READ_LEN))
        if kind == "chimeric":
            read2 = _revcomp(chr2[p : p + n]) + _revcomp(chr1[u : u + READ_LEN - n])
        else:
            read2 = _revcomp(chr2[p : p + READ_LEN])
        reads1.append((name, read1))
        reads2.append((name, read2))
        truth[name] = (kind, p + 100)

    genome = root / "genome.fa"
    genome.write_text(f">chr1\n{chr1}\n>chr2\n{chr2}\n")
    sp.run(["bowtie2-build", "-q", str(genome), str(genome)], check=True)
    sp.run(["bwa", "index", str(genome)], check=True, stderr=sp.DEVNULL)
    _write_fastq(root / "reads1.fq", reads1)
    _write_fastq(root / "reads2.fq", reads2)
    info_contigs = root / "info_contigs.txt"
    info_contigs.write_text(f"contig\tlength\nchr1\t{CHROM_LEN}\nchr2\t{CHROM_LEN}\n")
    return {
        "genome": str(genome),
        "reads1": str(root / "reads1.fq"),
        "reads2": str(root / "reads2.fq"),
        "info_contigs": str(info_contigs),
        "truth": truth,
    }


def _records(bam):
    with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as handle:
        return list(handle)


@pytest.mark.parametrize("aligner", ["bwa", "minimap2"])
def test_library_has_supplementary_records(chimeric_lib, aligner, tmp_path):
    """Premise of the tests below: at full length, the aligner reports the
    unique 3' end of chimeric reads as their primary record, with a MAPQ that
    passes, and their repeat 5' end as a supplementary record."""
    if aligner == "bwa":
        cmd = ["bwa", "mem", "-v", "1", chimeric_lib["genome"], chimeric_lib["reads1"]]
    else:
        cmd = ["minimap2", "-x", "sr", "-a", chimeric_lib["genome"], chimeric_lib["reads1"]]
    sam = tmp_path / "full.sam"
    with open(sam, "w") as out:
        sp.run(cmd, check=True, stdout=out, stderr=sp.DEVNULL)
    with pysam.AlignmentFile(str(sam), "r") as handle:
        records = list(handle)
    for name, (kind, _) in chimeric_lib["truth"].items():
        if kind != "chimeric":
            continue
        primary = [r for r in records if r.query_name == name and not r.flag & 0x900]
        supplementary = [r for r in records if r.query_name == name and r.is_supplementary]
        assert len(primary) == 1 and len(supplementary) == 1
        assert primary[0].reference_name == "chr2"
        assert primary[0].mapping_quality >= MIN_QUAL
        assert supplementary[0].reference_name == "chr1"


@pytest.mark.parametrize("aligner", ["bowtie2", "bwa", "minimap2"])
def test_iterative_align_writes_each_read_once(chimeric_lib, aligner, tmp_path):
    """Each read is written once, with its primary record: from the round its
    primary record passed in, or from the unaligned reads of the last round."""
    out = tmp_path / "reads1.bam"
    hci.iterative_align(
        chimeric_lib["reads1"],
        str(tmp_path),
        chimeric_lib["genome"],
        1,
        str(out),
        aligner=aligner,
        min_qual=MIN_QUAL,
    )
    records = _records(out)
    counts = collections.Counter(r.query_name for r in records)
    assert {name: n for name, n in counts.items() if n > 1} == {}
    assert sorted(counts) == sorted(chimeric_lib["truth"])
    # Only primary records: one record per read, as non-iterative alignment
    assert not [r.query_name for r in records if r.flag & 0x900]
    for r in records:
        kind, start = chimeric_lib["truth"][r.query_name]
        # bowtie2 reports no supplementary records and does not pass chimeric
        # reads here: its local MAPQ for the unique end stays below 30 (the
        # repeat scores close behind), and its last, end-to-end round leaves
        # them unmapped.
        if kind == "chimeric" and aligner != "bowtie2":
            # The read is kept with the primary record that passed: its unique 3' end
            assert r.reference_name == "chr2" and abs(r.reference_start - start) <= 5
            assert r.mapping_quality >= MIN_QUAL
        elif kind == "repeat":
            # Never passed: written once with the primary record of the last round
            assert r.mapping_quality < MIN_QUAL
        elif kind == "unmapped":
            assert r.is_unmapped


@pytest.mark.parametrize("aligner", ["bwa", "minimap2"])
def test_iterative_pairs_written_once(chimeric_lib, aligner, tmp_path, caplog):
    """The pipeline writes one pair per read pair whose mates both pass."""
    bams = []
    for mate in ("reads1", "reads2"):
        bam = str(tmp_path / f"{mate}.bam")
        hpi.align_reads(
            chimeric_lib[mate],
            chimeric_lib["genome"],
            bam,
            tmp_dir=str(tmp_path),
            aligner=aligner,
            iterative=True,
            min_qual=MIN_QUAL,
        )
        bams.append(bam)
    pairs = str(tmp_path / "valid.pairs")
    hpi.bam2pairs(*bams, pairs, chimeric_lib["info_contigs"], min_qual=MIN_QUAL)
    with open(pairs) as handle:
        ids = [line.split("\t", 1)[0] for line in handle if not line.startswith("#")]
    expected = sorted(
        name for name, (kind, _) in chimeric_lib["truth"].items() if kind in ("chimeric", "unique")
    )
    assert sorted(ids) == expected
    assert "only present in one BAM file" not in caplog.text


# ---------------------------------------------------------------------------
# _filter_bamfile, on hand-written records
# ---------------------------------------------------------------------------


def _write_bam(path, records):
    """Write a BAM file of single-end records given as
    (name, flag, reference_id, position, mapping quality)."""
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for name, flag, ref, pos, mapq in records:
            r = pysam.AlignedSegment(bam.header)
            r.query_name = name
            r.flag = flag
            r.reference_id = ref
            r.reference_start = pos
            r.mapping_quality = mapq
            if ref >= 0:
                r.cigarstring = "20M"
            r.query_sequence = "A" * 20
            r.query_qualities = pysam.qualitystring_to_array("I" * 20)
            bam.write(r)


def test_filter_bamfile_classifies_reads_by_primary_record(tmp_path):
    """Secondary (0x100) and supplementary (0x800) records neither keep a
    read nor send it to the next round, and are not written."""
    records = [
        # Passing primary record, with a supplementary and a secondary one
        ("a", 0, 0, 10, 60),
        ("a", 2048, 0, 500, 0),
        ("a", 256, 0, 700, 0),
        # Passing reverse primary record, at exactly the threshold
        ("b", 16, 0, 20, 30),
        ("b", 2064, 0, 600, 30),
        # Primary record below the threshold, with a passing supplementary
        ("c", 0, 0, 30, 10),
        ("c", 2048, 0, 300, 60),
        # Unmapped
        ("d", 4, -1, -1, 0),
        # Primary record below the threshold, with a passing secondary
        ("e", 16, 0, 40, 0),
        ("e", 256, 0, 400, 60),
    ]
    bam_in, bam_out = tmp_path / "in.bam", tmp_path / "out.bam"
    _write_bam(bam_in, records)
    remaining = hci._filter_bamfile(str(bam_in), str(bam_out), min_qual=30)
    assert remaining == {"c", "d", "e"}
    kept = [(r.query_name, r.flag, r.reference_start) for r in _records(bam_out)]
    assert kept == [("a", 0, 10), ("b", 16, 20)]
