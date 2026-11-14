"""Deterministic synthetic Hi-C library for the test-suite.

The generator writes, into a directory of the caller's choosing:

* ``genome.fa`` (and ``genome.fa.gz``): three chromosomes of 10-30 kb whose
  base composition is drawn per chromosome (GC 40-60 %). Every natural DpnII
  (``GATC``) and HindIII (``AAGCTT``) site is erased first, then sites are
  written back at known positions, so both restriction maps are exactly known
  and every DpnII fragment is long enough to hold a whole read.
* ``reads_R1.fq.gz`` / ``reads_R2.fq.gz``: a DpnII Hi-C library. Each read pair
  comes from the two sides of a ligation junction: R1 is read from the far end
  of one ligated fragment towards the junction, R2 from the far end of the
  other one. Intra-chromosomal contacts decay as ``s ** -alpha`` with the
  genomic distance ``s``; a known fraction of contacts is inter-chromosomal.
  PCR duplicates, uncut fragments, self-ligations (loops), junction-spanning
  (chimeric) R1 reads and pairs with an unmappable mate are injected at known
  rates.
* ``truth.tsv``: one row per read pair with what an aligner should report for
  both mates (chromosome, strand, leftmost position, restriction fragment) and
  how the pair was generated.

Everything is drawn from a single seeded :class:`numpy.random.Generator`, so
the same parameters always produce byte-identical files.

Physical model
--------------
DpnII cuts ``^GATC``. After fill-in, a fragment spanning ``[start, end)``
covers ``[start, end + 4)`` and both of its ends carry a ``GATC``; joining
two ends yields the ``GATCGATC`` junction. A read pointing towards the
junction from fragment ``X`` maps on the forward strand when the junction is
at ``X``'s right end, and on the reverse strand when it is at its left end.
With hicstuff's conventions (pairs ordered by position), this gives:

* ``+-`` for uncut fragments (reads facing each other),
* ``-+`` for self-ligated fragments (reads pointing away from each other),
* all four orientations, equally likely, for genuine contacts.
"""

from __future__ import annotations

import gzip
import json
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

#: Recognition site and 0-based offset of the top-strand cut inside it.
ENZYMES = {"DpnII": ("GATC", 0), "HindIII": ("AAGCTT", 1)}
#: Enzyme used to build the simulated library.
LIBRARY_ENZYME = "DpnII"
#: Ligation junction of two filled-in DpnII ends.
JUNCTION = "GATCGATC"
#: Length of the DpnII 5' overhang (filled in before ligation).
OVERHANG = 4

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
_BASES = np.array(list("ACGT"))
_KINDS = ("contact", "uncut", "loop", "duplicate", "unmapped")


def revcomp(seq: str) -> str:
    """Reverse complement of a DNA string."""
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class SimParams:
    """Parameters of a synthetic library. Rates are fractions of ``n_pairs``."""

    seed: int = 1984
    n_chroms: int = 3
    chrom_len_range: Tuple[int, int] = (10_000, 30_000)
    gc_range: Tuple[float, float] = (0.40, 0.60)
    #: Distance between consecutive DpnII sites (= DpnII fragment length).
    dpnii_spacing: Tuple[int, int] = (150, 600)
    #: Distance between consecutive HindIII sites.
    hindiii_spacing: Tuple[int, int] = (2_500, 6_000)
    n_pairs: int = 6_000
    read_len_range: Tuple[int, int] = (50, 100)
    #: Maximum distance between a read's far end and the ligation junction.
    max_half_insert: int = 400
    #: Exponent of the contact probability decay, P(s) ~ s ** -alpha.
    alpha: float = 1.2
    #: Shortest genomic distance of a simulated cis contact.
    s_min: int = 1_000
    #: Fraction of genuine contacts that are inter-chromosomal.
    trans_fraction: float = 0.10
    #: Fraction of genuine contacts whose R1 reads through the junction.
    chimeric_fraction: float = 0.10
    uncut_rate: float = 0.06
    loop_rate: float = 0.04
    dup_rate: float = 0.05
    unmapped_rate: float = 0.02


@dataclass
class SyntheticLibrary:
    """Paths to a generated library and its ground truth."""

    root: Path
    params: SimParams
    genome: Path
    genome_gz: Path
    reads1: Path
    reads2: Path
    truth_path: Path
    chroms: Dict[str, str]
    #: 0-based positions of the first base of each recognition site.
    sites: Dict[str, Dict[str, np.ndarray]]
    truth: pd.DataFrame = field(repr=False)

    @property
    def chrom_names(self) -> List[str]:
        return list(self.chroms)

    @property
    def chrom_sizes(self) -> Dict[str, int]:
        return {name: len(seq) for name, seq in self.chroms.items()}

    def boundaries(self, enzymes="DpnII") -> Dict[str, np.ndarray]:
        """Fragment boundaries (0, cut positions..., chrom length) per chromosome.

        ``enzymes`` is an enzyme name or a list of names (cut sites are merged).
        """
        if isinstance(enzymes, str):
            enzymes = [enzymes]
        out = {}
        for chrom, seq in self.chroms.items():
            cuts = set()
            for enz in enzymes:
                offset = ENZYMES[enz][1]
                cuts.update(int(p) + offset for p in self.sites[enz][chrom])
            out[chrom] = np.array(sorted({0, len(seq)} | cuts))
        return out

    def fragments(self, enzymes="DpnII") -> pd.DataFrame:
        """Expected ``fragments_list.txt`` content (without gc_content)."""
        rows = []
        for chrom, bounds in self.boundaries(enzymes).items():
            for i, (start, end) in enumerate(zip(bounds[:-1], bounds[1:])):
                rows.append((i + 1, chrom, int(start), int(end), int(end - start)))
        return pd.DataFrame(rows, columns=["id", "chrom", "start_pos", "end_pos", "size"])

    def frag_offsets(self, enzymes="DpnII") -> Dict[str, int]:
        """Global index of the first fragment of each chromosome."""
        offsets, total = {}, 0
        for chrom, bounds in self.boundaries(enzymes).items():
            offsets[chrom] = total
            total += len(bounds) - 1
        return offsets

    def mapped_truth(self) -> pd.DataFrame:
        """Truth rows of pairs whose two mates are expected to align."""
        return self.truth[self.truth.kind != "unmapped"]


# ---------------------------------------------------------------------------
# Genome
# ---------------------------------------------------------------------------


def _random_seq(rng: np.random.Generator, length: int, gc: float) -> List[str]:
    probs = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return list(rng.choice(_BASES, size=length, p=probs))


def _erase_motifs(rng: np.random.Generator, seq: List[str], motifs: List[str]) -> None:
    """Mutate ``seq`` in place until none of ``motifs`` occurs, including
    occurrences wrapping around the end (Bio.Restriction searches circularly)."""
    pattern = re.compile("|".join(motifs))
    span = max(len(m) for m in motifs)
    n = len(seq)
    while True:
        text = "".join(seq)
        circular = text + text[: span - 1]
        hits = [m.start() for m in pattern.finditer(circular)]
        if not hits:
            return
        for start in hits:
            pos = (start + int(rng.integers(0, 2))) % n
            seq[pos] = str(rng.choice([b for b in "ACGT" if b != seq[pos]]))


def _place_sites(
    rng: np.random.Generator, length: int, spacing: Tuple[int, int], margin: int
) -> List[int]:
    """Positions separated by a random distance in ``spacing``."""
    positions = []
    pos = int(rng.integers(spacing[0], spacing[1] + 1))
    while pos < length - margin:
        positions.append(pos)
        pos += int(rng.integers(spacing[0], spacing[1] + 1))
    return positions


def _make_genome(rng: np.random.Generator, params: SimParams):
    chroms, sites = {}, {enz: {} for enz in ENZYMES}
    motifs = [site for site, _ in ENZYMES.values()]
    lo, hi = params.dpnii_spacing
    for i in range(params.n_chroms):
        name = f"chr{i + 1}"
        length = int(rng.integers(params.chrom_len_range[0], params.chrom_len_range[1] + 1))
        gc = float(rng.uniform(*params.gc_range))
        seq = _random_seq(rng, length, gc)
        _erase_motifs(rng, seq, motifs)
        # HindIII sites first, then DpnII sites at a safe distance from them.
        hind = _place_sites(rng, length, params.hindiii_spacing, params.hindiii_spacing[0])
        dpn = [
            p
            for p in _place_sites(rng, length, params.dpnii_spacing, lo)
            if all(abs(p - h) > 10 for h in hind)
        ]
        # Dropping DpnII sites next to HindIII ones may leave long fragments:
        # that is fine, fragments are only required to be at least `lo` long.
        for p in hind:
            seq[p : p + 6] = list(ENZYMES["HindIII"][0])
        for p in dpn:
            seq[p : p + 4] = list(ENZYMES["DpnII"][0])
        text = "".join(seq)
        # The written sites must be the only ones (no accidental site was
        # created at the edge of a written one).
        for enz, expected in (("HindIII", hind), ("DpnII", dpn)):
            motif = ENZYMES[enz][0]
            circular = text + text[: len(motif) - 1]
            found = [m.start() for m in re.finditer(f"(?={motif})", circular)]
            assert found == expected, f"unexpected {enz} sites in {name}"
        chroms[name] = text
        sites["HindIII"][name] = np.array(hind, dtype=int)
        sites["DpnII"][name] = np.array(dpn, dtype=int)
    return chroms, sites


# ---------------------------------------------------------------------------
# Library
# ---------------------------------------------------------------------------


class _Sampler:
    """Draws read pairs from a genome with a DpnII restriction map."""

    def __init__(self, rng, params, chroms, bounds, offsets):
        self.rng = rng
        self.p = params
        self.chroms = chroms
        self.names = list(chroms)
        self.lengths = np.array([len(chroms[c]) for c in self.names], dtype=float)
        self.bounds = bounds
        self.offsets = offsets

    # -- helpers -----------------------------------------------------------

    def _chrom(self, exclude: str | None = None) -> str:
        weights = self.lengths.copy()
        if exclude is not None:
            weights[self.names.index(exclude)] = 0
        return self.names[int(self.rng.choice(len(self.names), p=weights / weights.sum()))]

    def _frag(self, chrom: str, pos: int) -> int:
        """Local index of the fragment containing 0-based ``pos``."""
        return int(np.searchsorted(self.bounds[chrom], pos, side="right") - 1)

    def _frag_span(self, chrom: str, idx: int) -> Tuple[int, int]:
        b = self.bounds[chrom]
        return int(b[idx]), int(b[idx + 1])

    def _read_len(self, lo: int | None = None) -> int:
        lo = self.p.read_len_range[0] if lo is None else lo
        return int(self.rng.integers(lo, self.p.read_len_range[1] + 1))

    def _ends(self, chrom: str, idx: int) -> List[str]:
        """Ligatable ends of a fragment (chromosome ends are not cut sites)."""
        ends = []
        if idx > 0:
            ends.append("left")
        if idx < len(self.bounds[chrom]) - 2:
            ends.append("right")
        return ends

    def _read(self, chrom, strand, start, read_len, aligned=None):
        """Truth record for one mate. ``start`` is the 0-based leftmost base of
        the aligned part; ``five`` the 0-based position of the 5' end."""
        aligned = read_len if aligned is None else aligned
        five = start if strand == "+" else start + aligned - 1
        local = self._frag(chrom, start)
        return {
            "chr": chrom,
            "strand": strand,
            "start": start,
            "five": five,
            "len": read_len,
            "alen": aligned,
            "frag": self.offsets[chrom] + local,
            "lfrag": local,
        }

    def _side(self, chrom, idx, end, m):
        """Molecule segment of length ``m`` ending at the ligated ``end`` of a
        fragment, oriented 5'->3' towards the junction, and its genomic
        interval and strand."""
        start, stop = self._frag_span(chrom, idx)
        seq = self.chroms[chrom]
        if end == "right":
            lo = stop + OVERHANG - m
            return seq[lo : stop + OVERHANG], lo, "+"
        lo = start
        return revcomp(seq[start : start + m]), lo, "-"

    # -- pair kinds -------------------------------------------------------

    def _contact_loci(self):
        rng, p = self.rng, self.p
        if rng.random() < p.trans_fraction:
            c1 = self._chrom()
            c2 = self._chrom(exclude=c1)
            x1 = int(rng.integers(0, len(self.chroms[c1])))
            x2 = int(rng.integers(0, len(self.chroms[c2])))
            return c1, x1, c2, x2
        chrom = self._chrom()
        length = len(self.chroms[chrom])
        a, b, k = float(p.s_min), float(length), 1.0 - p.alpha
        while True:
            # Rejection on x1 + s keeps P(s) ~ s^-alpha per pair of loci.
            u = rng.random()
            s = (a**k + u * (b**k - a**k)) ** (1.0 / k)
            x1 = int(rng.integers(0, length))
            x2 = x1 + int(round(s)) * (1 if rng.random() < 0.5 else -1)
            if 0 <= x2 < length:
                return chrom, x1, chrom, x2

    def contact(self, chimeric: bool):
        rng, p = self.rng, self.p
        while True:
            c1, x1, c2, x2 = self._contact_loci()
            f1, f2 = self._frag(c1, x1), self._frag(c2, x2)
            if c1 == c2 and f1 == f2:
                continue
            ends1, ends2 = self._ends(c1, f1), self._ends(c2, f2)
            if chimeric:
                # Chimeric R1 reads are forward reads crossing the junction
                if "right" not in ends1:
                    continue
                e1 = "right"
            else:
                e1 = str(rng.choice(ends1))
            e2 = str(rng.choice(ends2))
            break
        len1 = self._read_len(lo=80 if chimeric else None)
        len2 = self._read_len()
        size1 = np.diff(self._frag_span(c1, f1))[0]
        size2 = np.diff(self._frag_span(c2, f2))[0]
        if chimeric:
            m1 = int(rng.integers(int(np.ceil(0.6 * len1)), len1 - 25 + 1))
        else:
            m1 = int(rng.integers(len1, min(size1, p.max_half_insert) + 1))
        m2 = int(rng.integers(len2, min(size2, p.max_half_insert) + 1))
        side1, lo1, strand1 = self._side(c1, f1, e1, m1)
        side2, lo2, strand2 = self._side(c2, f2, e2, m2)
        molecule = side1 + revcomp(side2)
        seq1 = molecule[:len1]
        seq2 = revcomp(molecule)[:len2]
        if strand1 == "+":
            read1 = self._read(c1, "+", lo1, len1, aligned=min(m1, len1))
        else:
            read1 = self._read(c1, "-", lo1 + m1 - len1, len1)
        if strand2 == "+":
            read2 = self._read(c2, "+", lo2, len2)
        else:
            read2 = self._read(c2, "-", lo2 + m2 - len2, len2)
        lig1 = self._frag_span(c1, f1)[1 if e1 == "right" else 0]
        lig2 = self._frag_span(c2, f2)[1 if e2 == "right" else 0]
        info = {
            "kind": "contact",
            "chimeric": chimeric,
            "chimeric_split": m1 if chimeric else -1,
            "lig_end1": e1,
            "lig_end2": e2,
            "lig_dist": abs(lig1 - lig2) if c1 == c2 else -1,
        }
        return seq1, seq2, read1, read2, info

    def uncut(self):
        """Undigested molecule: forward read at one end, reverse at the other,
        spanning at most two restriction sites."""
        rng = self.rng
        while True:
            chrom = self._chrom()
            length = len(self.chroms[chrom])
            la, lb = self._read_len(), self._read_len()
            insert = int(rng.integers(max(la, lb) + 50, 501))
            s = int(rng.integers(0, length - insert))
            fa, fa_end = self._frag(chrom, s), self._frag(chrom, s + la - 1)
            fb, fb_end = self._frag(chrom, s + insert - lb), self._frag(chrom, s + insert - 1)
            if fa != fa_end or fb != fb_end or fb - fa > 2:
                continue
            break
        seq = self.chroms[chrom]
        plus = (seq[s : s + la], self._read(chrom, "+", s, la))
        minus = (
            revcomp(seq[s + insert - lb : s + insert]),
            self._read(chrom, "-", s + insert - lb, lb),
        )
        first, second = (plus, minus) if rng.random() < 0.5 else (minus, plus)
        return first[0], second[0], first[1], second[1], {"kind": "uncut"}

    def loop(self):
        """Self-circularised fragment: both mates in one fragment, pointing
        away from each other."""
        rng, p = self.rng, self.p
        while True:
            chrom = self._chrom()
            idx = int(rng.integers(1, len(self.bounds[chrom]) - 2))
            start, stop = self._frag_span(chrom, idx)
            la, lb = self._read_len(), self._read_len()
            if stop - start < la + lb + 10:
                continue
            break
        m1 = int(rng.integers(la, min(stop - start - lb, p.max_half_insert) + 1))
        m2 = int(rng.integers(lb, min(stop - start - m1, p.max_half_insert) + 1))
        seq = self.chroms[chrom]
        lo_plus = stop + OVERHANG - m1
        plus = (seq[lo_plus : lo_plus + la], self._read(chrom, "+", lo_plus, la))
        lo_minus = start + m2 - lb
        minus = (revcomp(seq[lo_minus : lo_minus + lb]), self._read(chrom, "-", lo_minus, lb))
        first, second = (plus, minus) if rng.random() < 0.5 else (minus, plus)
        return first[0], second[0], first[1], second[1], {"kind": "loop"}

    def unmapped(self):
        """R1 is random sequence that cannot align; R2 is a genuine read."""
        rng = self.rng
        la = self._read_len()
        while True:
            junk = "".join(rng.choice(_BASES, size=la))
            if not re.search(JUNCTION, junk):
                break
        _, seq2, _, read2, _ = self.contact(chimeric=False)
        return junk, seq2, None, read2, {"kind": "unmapped"}


def _flatten(read_id, read1, read2, info, dup_of=""):
    row = {"read_id": read_id, "dup_of": dup_of, "chimeric": False, "chimeric_split": -1}
    row.update({"lig_end1": "", "lig_end2": "", "lig_dist": -1})
    row.update(info)
    for mate, read in (("1", read1), ("2", read2)):
        keys = ("chr", "strand", "start", "five", "len", "alen", "frag", "lfrag")
        for key in keys:
            row[key + mate] = None if read is None else read[key]
    return row


def _write_fastq(path: Path, records) -> None:
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as fq:
        for name, seq in records:
            fq.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n".encode())


def _write_fasta(path: Path, chroms: Dict[str, str], width: int = 60) -> None:
    with open(path, "w") as fa:
        for name, seq in chroms.items():
            fa.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fa.write(seq[i : i + width] + "\n")


def generate(out_dir, params: SimParams | None = None) -> SyntheticLibrary:
    """Write a synthetic genome, Hi-C library and truth table to ``out_dir``."""
    params = SimParams() if params is None else params
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(params.seed)

    chroms, sites = _make_genome(rng, params)
    lib = SyntheticLibrary(
        root=out_dir,
        params=params,
        genome=out_dir / "genome.fa",
        genome_gz=out_dir / "genome.fa.gz",
        reads1=out_dir / "reads_R1.fq.gz",
        reads2=out_dir / "reads_R2.fq.gz",
        truth_path=out_dir / "truth.tsv",
        chroms=chroms,
        sites=sites,
        truth=pd.DataFrame(),
    )
    sampler = _Sampler(rng, params, chroms, lib.boundaries(LIBRARY_ENZYME), lib.frag_offsets())

    n = params.n_pairs
    counts = {
        kind: int(round(rate * n))
        for kind, rate in (
            ("uncut", params.uncut_rate),
            ("loop", params.loop_rate),
            ("duplicate", params.dup_rate),
            ("unmapped", params.unmapped_rate),
        )
    }
    counts["contact"] = n - sum(counts.values())

    pairs = []  # (seq1, seq2, read1, read2, info)
    for _ in range(counts["contact"]):
        pairs.append(sampler.contact(chimeric=bool(rng.random() < params.chimeric_fraction)))
    for kind in ("uncut", "loop", "unmapped"):
        pairs.extend(getattr(sampler, kind)() for _ in range(counts[kind]))
    # PCR duplicates are exact copies of genuine contacts.
    originals = rng.choice(counts["contact"], size=counts["duplicate"], replace=True)

    # Shuffle pairs so that kinds are interleaved; read names follow file order.
    n_total = len(pairs) + len(originals)
    slot = np.empty(n_total, dtype=int)
    slot[rng.permutation(n_total)] = np.arange(n_total)
    names = [f"sim{i:06d}" for i in range(n_total)]
    rows, fq1, fq2 = [None] * n_total, [None] * n_total, [None] * n_total
    for k, (seq1, seq2, read1, read2, info) in enumerate(pairs):
        i = int(slot[k])
        rows[i] = _flatten(names[i], read1, read2, info)
        fq1[i], fq2[i] = (names[i], seq1), (names[i], seq2)
    for j, orig in enumerate(originals):
        i = int(slot[len(pairs) + j])
        seq1, seq2, read1, read2, info = pairs[int(orig)]
        info = dict(info, kind="duplicate")
        rows[i] = _flatten(names[i], read1, read2, info, dup_of=names[int(slot[int(orig)])])
        fq1[i], fq2[i] = (names[i], seq1), (names[i], seq2)

    truth = pd.DataFrame(rows)
    for col in ("start1", "five1", "len1", "alen1", "frag1", "lfrag1"):
        truth[col] = truth[col].astype("Int64")
    lib.truth = truth

    _write_fasta(lib.genome, chroms)
    with open(lib.genome, "rb") as src, open(lib.genome_gz, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as dst:
            dst.write(src.read())
    _write_fastq(lib.reads1, fq1)
    _write_fastq(lib.reads2, fq2)
    truth.to_csv(lib.truth_path, sep="\t", index=False)
    with open(out_dir / "params.json", "w") as handle:
        json.dump(asdict(params), handle, indent=2)
    return lib


def read_pairs(path) -> pd.DataFrame:
    """Load a hicstuff pairs file (plain or gzipped) into a DataFrame."""
    names = ["readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2", "frag1", "frag2"]
    df = pd.read_csv(path, sep="\t", comment="#", header=None, dtype={0: str, 1: str, 3: str})
    df.columns = names[: df.shape[1]]
    return df
