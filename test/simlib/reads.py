"""Fragment / read-pair sampler with Illumina-like noise, plus alignment-by-construction.

Noise model
-----------
* Fragment (insert) length ~ N(330, 90), >= read_len + 10 (calibrated: val1 CRAM median 319).
* Homopolymer (SBS slippage) jitter: every read crossing a run of >= 8 identical bases sees
  the run length changed by round(N(0, sd)), sd = scale * min(0.3 + 0.045*(n-8), 3.5) —
  ~+-1 for a 20-mer, ~+-3 for a 70-mer poly-A. Drawn independently per read (so PCR
  duplicates and the two mates of one fragment disagree, as on a real flow cell). When a run
  straddles inserted and reference sequence (a poly-A abutting an A-rich TSD) the jitter is
  placed on the inserted part.
* Substitution errors at `error_rate`, elevated (x5) and lower quality for 40 cycles after a
  homopolymer >= 20 (phasing loss after long poly-A).
* PCR duplicates WITHOUT the 0x400 flag: with probability `pcr_dup_frac` a fragment gets
  1 + Geom(0.5) extra copies (<= 6) with IDENTICAL outer coordinates on both mates and
  independently drawn sequencing noise — the read-based independence rule must collapse them.

Alignment by construction (val1, genome-free)
---------------------------------------------
Each read base keeps its haplotype index; haplotype segments say whether a base is reference
(and where / which strand) or inserted. The primary alignment is the longest reference piece;
everything else is soft-clipped (bwa-like local alignment); a second reference piece >= 30 bp
becomes an `SA` tag + supplementary record; homopolymer jitter inside an aligned reference
piece becomes I/D unless it is < 8 bases from the piece end, where it is soft-clipped
(bwa's end bonus); post-homopolymer phasing junk (see add_noise) is soft-clipped — together
the organic source of poly-A slippage clips.
"""
import random
from bisect import bisect_right

from .seqs import homopolymer_runs, revcomp

QUAL_HI, QUAL_MID, QUAL_LO = "F", ":", ","   # 37 / 25 / 11 (NovaSeq binned)


class Hap:
    """A haplotype window: sequence + segment map (val1) + per-run jitter scale."""

    def __init__(self, seq, segs=None, name="hap", jitter_boost=None):
        self.seq = seq
        self.segs = segs or []
        self.name = name
        self.runs = homopolymer_runs(seq, 8)
        self._run_starts = [r[0] for r in self.runs]
        self._seg_starts = [s.a0 for s in self.segs]
        # optional {(run_start): multiplier} for artefact loci (exaggerated slippage)
        self.jitter_boost = jitter_boost or {}

    def seg_index(self, a):
        if not self.segs:
            return -1
        return bisect_right(self._seg_starts, a) - 1

    def runs_in(self, a, b):
        i = max(0, bisect_right(self._run_starts, a) - 1)
        out = []
        while i < len(self.runs) and self.runs[i][0] < b:
            r = self.runs[i]
            if r[1] > a:
                out.append(r)
            i += 1
        return out


def jitter_sd(n, scale):
    if n < 8:
        return 0.0
    return scale * min(0.3 + 0.045 * (n - 8), 3.5)


def read_entries(rng, hap, start, end, forward, read_len, jitter_scale, slip_min_delta=0):
    """Sequence one read from hap[start:end] (a fragment). forward=True reads from `start`
    rightwards; else from `end` leftwards. Returns a list of [base, hap_idx, extra] in HAP
    orientation (for a reverse read this is the reverse complement of what was sequenced)."""
    span = read_len + 60
    a, b = (start, min(end, start + span)) if forward else (max(start, end - span), end)
    pos = list(range(a, b))
    extra = [False] * len(pos)
    if jitter_scale > 0:
        for (r0, r1, base) in hap.runs_in(a, b):
            n = r1 - r0
            boost = hap.jitter_boost.get(r0, 1.0)
            sd = jitter_sd(n, jitter_scale * boost)
            if sd <= 0:
                continue
            d = int(round(rng.gauss(0, sd)))
            if boost > 1.0 and slip_min_delta and abs(d) < slip_min_delta and rng.random() < 0.6:
                d = slip_min_delta * (1 if rng.random() < 0.5 else -1)
            if d == 0:
                continue
            # candidate positions inside the run and inside this read window; prefer INS
            cand = [i for i, p in enumerate(pos) if r0 <= p < r1 and not extra[i]]
            if not cand:
                continue
            if hap.segs:
                ins_c = [i for i in cand if hap.segs[hap.seg_index(pos[i])].kind == "INS"]
                if ins_c:
                    cand = ins_c
            if d < 0:
                k = min(-d, len(cand) - 1)
                drop = set(cand[-k:] if forward else cand[:k]) if k > 0 else set()
                pos = [p for i, p in enumerate(pos) if i not in drop]
                extra = [e for i, e in enumerate(extra) if i not in drop]
            else:
                at = cand[len(cand) // 2]
                pos = pos[:at] + [pos[at]] * d + pos[at:]
                extra = extra[:at] + [True] * d + extra[at:]
    if forward:
        pos, extra = pos[:read_len], extra[:read_len]
    else:
        pos, extra = pos[-read_len:], extra[-read_len:]
    seq = hap.seq
    return [[seq[p], p, e] for p, e in zip(pos, extra)]


def add_noise(rng, entries, forward, error_rate, burst_scale=1.0):
    """Substitution errors + qualities, in SEQUENCING order; returns (seq, qual) in hap
    orientation. After a homopolymer >= 15 the read may lose phasing for the rest of the
    cycles (probability min(0.6, 0.01*(n-12)) * burst_scale): the remaining bases become
    low-quality, mostly homopolymer-base junk (hap index set to -1, so an aligner clips them).
    This is the read-level origin of poly-A slippage clips at reference A-tracts and of
    reads that do NOT carry the 3'-beyond-poly-A sequence."""
    n = len(entries)
    order = list(range(n)) if forward else list(range(n - 1, -1, -1))
    quals = [QUAL_HI] * n
    run_base, run_len, phase_left, junk = None, 0, 0, None
    for k, i in enumerate(order):
        b = entries[i][0]
        if junk is not None:
            entries[i][0] = junk if rng.random() < 0.7 else rng.choice("ACGT")
            entries[i][1] = -1
            quals[i] = rng.choice("#',")
            continue
        if b == run_base:
            run_len += 1
        else:
            if run_len >= 15 and burst_scale > 0 and \
                    rng.random() < min(0.6, 0.01 * (run_len - 12)) * burst_scale:
                junk = run_base
                entries[i][0] = junk if rng.random() < 0.7 else rng.choice("ACGT")
                entries[i][1] = -1
                quals[i] = "#"
                continue
            if run_len >= 20:
                phase_left = 40
            run_base, run_len = b, 1
        er = error_rate * (5 if phase_left > 0 else 1)
        if phase_left > 0:
            quals[i] = QUAL_MID
            phase_left -= 1
        if rng.random() < er:
            entries[i][0] = rng.choice([x for x in "ACGT" if x != b])
            quals[i] = QUAL_LO
    return "".join(e[0] for e in entries), "".join(quals)


class FragmentSampler:
    def __init__(self, rng, read_len=151, flen_mu=330, flen_sd=90, jitter_scale=1.0,
                 error_rate=0.002, pcr_dup_frac=0.0, slip_min_delta=0, burst_scale=1.0):
        self.rng = rng
        self.read_len = read_len
        self.flen_mu, self.flen_sd = flen_mu, flen_sd
        self.jitter_scale = jitter_scale
        self.error_rate = error_rate
        self.pcr_dup_frac = pcr_dup_frac
        self.slip_min_delta = slip_min_delta
        self.burst_scale = burst_scale

    def flen(self):
        return max(self.read_len + 10, int(self.rng.gauss(self.flen_mu, self.flen_sd)))

    def n_fragments(self, hap_len, depth, weight=1.0):
        lam = depth * hap_len / (2.0 * self.read_len) * weight
        # Poisson via normal approx for large lam
        if lam > 50:
            return max(0, int(round(self.rng.gauss(lam, lam ** 0.5))))
        L, k, p = pow(2.718281828, -lam), 0, 1.0
        while True:
            p *= self.rng.random()
            if p < L:
                return k
            k += 1

    def copies(self):
        if self.pcr_dup_frac <= 0 or self.rng.random() >= self.pcr_dup_frac:
            return 0
        k = 1
        while self.rng.random() < 0.5 and k < 6:
            k += 1
        return k

    def pair(self, hap, start, end):
        """Return (read_fwd, read_rev) each = (entries, seq, qual, forward). Which one is
        read1 is decided by the caller (coin flip)."""
        out = []
        for fwd in (True, False):
            e = read_entries(self.rng, hap, start, end, fwd, self.read_len, self.jitter_scale,
                             self.slip_min_delta)
            seq, qual = add_noise(self.rng, e, fwd, self.error_rate, self.burst_scale)
            out.append((e, seq, qual, fwd))
        return out


# ------------------------------------------------------------------ alignment by construction
def _cigar_str(ops):
    out, last, n = [], None, 0
    for op in ops:
        if op == last:
            n += 1
        else:
            if last is not None:
                out.append(f"{n}{last}")
            last, n = op, 1
    if last is not None:
        out.append(f"{n}{last}")
    return "".join(out)


def align_entries(hap, entries, min_aln=20, min_supp=30, end_indel_clip=8):
    """Return dict(primary=..., supp=[...]) or None if unmapped. Each alignment:
    dict(ref_start (window coords), cigar ops list (in REF-forward orientation), reverse
    (bool, relative to hap orientation), q0, q1 (read slice in hap orientation))."""
    pieces = []          # (seg_idx, i0, i1)
    cur, i0 = None, 0
    for i, (b, p, ex) in enumerate(entries):
        si = hap.seg_index(p) if p >= 0 else -2
        if si != cur:
            if cur is not None:
                pieces.append((cur, i0, i))
            cur, i0 = si, i
    if cur is not None:
        pieces.append((cur, i0, len(entries)))
    ref_pieces = []
    for si, a, b in pieces:
        if si < 0 or hap.segs[si].kind != "REF":
            continue
        nm = sum(1 for k in range(a, b) if not entries[k][2])
        ref_pieces.append((nm, si, a, b))
    if not ref_pieces:
        return None
    ref_pieces.sort(key=lambda x: (-x[0], x[2]))
    alns = []
    for rank, (nm, si, a, b) in enumerate(ref_pieces):
        if nm < (min_aln if rank == 0 else min_supp):
            continue
        seg = hap.segs[si]
        ops, rpos = [], []
        prev_r = None
        for k in range(a, b):
            base, p, ex = entries[k]
            if ex:
                ops.append("I"); rpos.append(None)
                continue
            r = seg.refpos(p)
            if prev_r is not None:
                gap = (r - prev_r - 1) if seg.strand == "+" else (prev_r - r - 1)
                ops.extend(["D"] * gap); rpos.extend([None] * gap)
            ops.append("M"); rpos.append(r)
            prev_r = r
        # bwa end-bonus behaviour: an indel within `end_indel_clip` M of a piece end -> clip
        q0, q1 = a, b
        def m_before(idx):
            return sum(1 for o in ops[:idx] if o == "M")
        changed = True
        while changed:
            changed = False
            idxs = [i for i, o in enumerate(ops) if o in "ID"]
            if not idxs:
                break
            first, last = idxs[0], idxs[-1]
            if m_before(first) < end_indel_clip:
                # clip everything up to and including the indel block
                j = first
                while j < len(ops) and ops[j] in "ID":
                    j += 1
                q0 += sum(1 for o in ops[:j] if o != "D")
                ops, rpos = ops[j:], rpos[j:]
                changed = True
                continue
            m_after = sum(1 for o in ops[last + 1:] if o == "M")
            if m_after < end_indel_clip:
                j = last
                while j >= 0 and ops[j] in "ID":
                    j -= 1
                q1 -= sum(1 for o in ops[j + 1:] if o != "D")
                ops, rpos = ops[:j + 1], rpos[:j + 1]
                changed = True
        mcount = sum(1 for o in ops if o == "M")
        if mcount < (min_aln if rank == 0 else min_supp):
            continue
        refs = [r for r in rpos if r is not None]
        if seg.strand == "+":
            ref_start, rev = refs[0], False
        else:
            ref_start, rev = refs[-1], True
            ops = ops[::-1]
        alns.append({"ref_start": ref_start, "ops": ops, "reverse": rev, "q0": q0, "q1": q1,
                     "seg": si, "ref_end": max(refs) + 1})
        if len(alns) == 2:
            break
    if not alns:
        return None
    return {"primary": alns[0], "supp": alns[1:]}


def full_cigar(aln, read_len, hard=False):
    """CIGAR string in ref-forward orientation including clips."""
    clip = "H" if hard else "S"
    lead, trail = aln["q0"], read_len - aln["q1"]
    if aln["reverse"]:
        lead, trail = trail, lead
    core = _cigar_str(aln["ops"])
    s = (f"{lead}{clip}" if lead else "") + core + (f"{trail}{clip}" if trail else "")
    return s


def entries_span(entries):
    ps = [p for _, p, _ in entries if p >= 0]
    if not ps:
        return 0, 0
    return min(ps), max(ps) + 1
