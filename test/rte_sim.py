"""Tiny insertion simulator for the tools/rte tests.

Builds one insertion at a synthetic site and returns everything annotate needs, in the SPEC
formats: the combined.txt.gz junction strings (lower clip / UPPER reference), the
insertions.evidence.tsv rows (JunctionEvidence) and the pooled reads (reference-forward,
mates oriented like the allele), plus the discovery genome as a FastaGenome.

    allele = ref[:R] + INSERT + ref[L:]      TSD = ref[L:R] (R > L), deletion if R < L
    title  = "chrS:L-R"
"""
from __future__ import annotations

import os
import random
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from tools.rte.genome import FastaGenome  # noqa: E402
from tools.rte.inputs import InsertionEvidence, JunctionEvidence, EvidenceRead  # noqa: E402
from tools.rte.annotator import InsertionInput  # noqa: E402
from tools.rte.sequtil import rc  # noqa: E402

FIX = os.path.join(REPO, "test", "fixtures", "rte_library")


def rnd(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def build(insert_sense, strand=1, tsd=15, motif="TTTTTAA", seed=1, read_len=150, frag=350,
          step=25, clip_len=100, polya_len=None, ref_homopolymer=0, exclude=(), samples=("S1",),
          contig="chrS", glen=8000, site=4000, flank_len=80):
    """insert_sense: inserted sequence in ELEMENT sense (5'..3' incl. poly-A).
    exclude: insert-sense coordinates; reads covering any of them (+-20 bp) are dropped."""
    rng = random.Random(seed)
    ref = list(rnd(glen, rng))
    L = site
    R = site + tsd
    if motif:
        if strand > 0:
            ref[L - 2:L + 5] = list(rc(motif))
        else:
            ref[R - 5:R + 2] = list(motif)
    if ref_homopolymer:
        # reference A-run (T-run for - strand) touching the 3' breakpoint: slippage context
        if strand > 0:
            ref[L:L + ref_homopolymer] = ["A"] * ref_homopolymer
        else:
            ref[R - ref_homopolymer:R] = ["T"] * ref_homopolymer
    ref = "".join(ref)
    ins = insert_sense if strand > 0 else rc(insert_sense)
    allele = ref[:R] + ins + ref[L:]
    a0, a1 = R, R + len(ins)        # insert span on the allele
    cl = min(clip_len, len(ins))
    right_str = ref[R - flank_len:R].upper() + ins[:cl].lower()
    left_str = ins[-cl:].lower() + ref[L:L + flank_len].upper()
    title = f"{contig}:{L}-{R}"

    # excluded allele intervals
    excl = []
    for p in exclude:
        q = p if strand > 0 else len(ins) - p
        excl.append((a0 + q - 20, a0 + q + 20))

    def keep(s, e):
        return not any(s < xe and e > xs for xs, xe in excl)

    reads = []
    nclip = {"LEFT": set(), "RIGHT": set()}
    fi = 0
    for si, sample in enumerate(samples):
        for f in range(a0 - frag + 30 + si * 7, a1 - 30, step):
            r1 = (f, f + read_len)
            r2 = (f + frag - read_len, f + frag)
            if r1[0] < 0 or r2[1] > len(allele):
                continue
            fi += 1
            fid = f"{fi:x}"
            recs = []
            for r12, (s, e) in ((1, r1), (2, r2)):
                if not keep(s, e):
                    recs = None
                    break
                in_ins_s, in_ins_e = max(s, a0), min(e, a1)
                if s < a0 - 20 and e > a0 + 20 and e - a0 >= 20:
                    recs.append(("RIGHT", "CLIP", r12, s, e))
                elif s < a1 - 20 and e > a1 + 20:
                    recs.append(("LEFT", "CLIP", r12, s, e))
                elif in_ins_e - in_ins_s >= read_len - 2:
                    recs.append(("", "MATE", r12, s, e))
                else:
                    recs.append(("", "REFONLY", r12, s, e))
            if not recs or all(r[1] == "REFONLY" for r in recs):
                continue
            anchor_side = "RIGHT" if r1[0] < a0 else "LEFT"
            for side, role, r12, s, e in recs:
                if role == "REFONLY":
                    role, side = "MATE", anchor_side
                if role == "MATE" and not side:
                    side = anchor_side
                if role == "CLIP":
                    nclip[side].add((sample, fid))
                reads.append(EvidenceRead(side, role, sample, fid, str(r12), allele[s:e]))
    ev = InsertionEvidence(title, reads=reads)
    pl = polya_len if polya_len is not None else _tail_len(insert_sense)
    for side, js in (("LEFT", left_str), ("RIGHT", right_str)):
        three = (side == "LEFT") == (strand > 0)
        n = len(nclip[side])
        ev.junctions[side] = JunctionEvidence(
            side=side, n_reads=n, n_fragments=n, n_independent=n, n_samples=len(samples),
            supported=int(n >= 2), clip_consensus=js,
            polya_len_median=float(pl) if three else 0.0)
    genome = FastaGenome(records={contig: ref})
    inp = InsertionInput(title, left_str, right_str)
    return inp, ev, genome, {"L": L, "R": R, "ref": ref, "allele": allele, "insert": ins}


def _tail_len(s):
    n = 0
    for ch in reversed(s.upper()):
        if ch != "A":
            break
        n += 1
    return n


def library():
    from tools.rte.library import RteLibrary
    return RteLibrary(FIX)


def fasta(name):
    import mappy
    return {n: s.upper() for n, s, _ in mappy.fastx_read(os.path.join(FIX, name))}


EVIDENCE_COLS = ["insertion_id", "side", "n_reads", "n_fragments", "n_independent", "n_samples",
                 "n_mates", "supported", "clip_consensus", "consensus_depth", "polya_len_median",
                 "polya_len_range", "beyond_polya", "beyond_polya_support"]


def write_sidecars(prefix, evs):
    """Write `<prefix>.insertions.evidence.tsv.gz` + `<prefix>.insertions.reads.fa.gz` (SPEC
    formats) for a list of InsertionEvidence."""
    import gzip
    with gzip.open(prefix + ".insertions.evidence.tsv.gz", "wt") as fh:
        fh.write("\t".join(EVIDENCE_COLS) + "\n")
        for ev in evs:
            for side, j in ev.junctions.items():
                row = [ev.insertion_id, side, j.n_reads, j.n_fragments, j.n_independent, j.n_samples,
                       j.n_mates, j.supported, j.clip_consensus or ".", ".", j.polya_len_median,
                       j.polya_len_range or ".", j.beyond_polya or ".", j.beyond_polya_support]
                fh.write("\t".join(str(x) for x in row) + "\n")
    with gzip.open(prefix + ".insertions.reads.fa.gz", "wt") as fh:
        for ev in evs:
            for r in ev.reads:
                fh.write(f">{ev.insertion_id}|{r.side}|{r.role}|{r.sample}|{r.frag}|{r.r12}\n{r.seq}\n")


def fastq_record(title, side, seq):
    return f"@{title}:{side}\n{seq}\n+\n{'I' * len(seq)}\n"
