"""Render catalogue events as aligned BAM records without a genome or aligner (val1).

Each event gets its own random reference window placed at a contig offset; reads are
sampled from the alt / ref haplotypes (VAF mixing per sample) and aligned by construction
(reads.align_entries). Mates are real records with correct flags, TLEN, MQ and MC tags;
reads that fall wholly inside inserted sequence are unmapped and placed at their mate.
"""
import pysam

from .models import Seg, apply_event, plant_en_motif
from .reads import Hap, align_entries, entries_span, full_cigar
from .seqs import degenerate_en_motif, revcomp, rnd_seq


def _rec(hdr, name, flag, tid, pos, mapq, cigar, seq, qual, tags=None):
    a = pysam.AlignedSegment(hdr)
    a.query_name = name
    a.flag = flag
    a.reference_id = tid
    a.reference_start = pos
    a.mapping_quality = mapq
    a.cigarstring = cigar
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array(qual)
    if tags:
        a.set_tags(tags)
    return a


class PairWriter:
    """Turns sampled read pairs into pysam records (primary + supplementary + unmapped)."""

    def __init__(self, hdr, tid, contig, base, read_len, mapq_fn=None):
        self.hdr, self.tid, self.contig, self.base = hdr, tid, contig, base
        self.read_len = read_len
        self.mapq_fn = mapq_fn or (lambda: 60)
        self.records = []

    def _one(self, hap, read):
        entries, seq, qual, fwd = read
        aln = align_entries(hap, entries)
        return {"entries": entries, "seq": seq, "qual": qual, "fwd": fwd, "aln": aln,
                "mapq": self.mapq_fn() if aln else 0}

    def write_pair(self, hap, pair, name, r1_is_fwd):
        a, b = self._one(hap, pair[0]), self._one(hap, pair[1])
        r1, r2 = (a, b) if r1_is_fwd else (b, a)
        for r, me, mate in ((r1, 0x40, r2), (r2, 0x80, r1)):
            r["first"] = me
        info = []
        for r in (r1, r2):
            if r["aln"]:
                p = r["aln"]["primary"]
                rev = (not r["fwd"]) ^ p["reverse"]
                r["pos"] = self.base + p["ref_start"]
                r["rev"] = rev
                r["cigar"] = full_cigar(p, len(r["seq"]))
                r["end"] = self.base + p["ref_end"]
            info.append(r)
        both = r1["aln"] and r2["aln"]
        tlen = 0
        if both:
            lo = min(r1["pos"], r2["pos"]); hi = max(r1["end"], r2["end"])
            tlen = hi - lo
        proper = bool(both and r1["rev"] != r2["rev"] and tlen <= 1000)
        for r, mate in ((r1, r2), (r2, r1)):
            flag = 0x1 | r["first"]
            if proper:
                flag |= 0x2
            if not r["aln"]:
                flag |= 0x4
            if not mate["aln"]:
                flag |= 0x8
            tags = []
            if r["aln"]:
                seq = r["seq"] if not r["aln"]["primary"]["reverse"] else revcomp(r["seq"])
                qual = r["qual"] if not r["aln"]["primary"]["reverse"] else r["qual"][::-1]
                if r["rev"]:
                    flag |= 0x10
                pos, cig, mq = r["pos"], r["cigar"], r["mapq"]
            else:
                # unmapped: original read orientation, placed at the mate
                seq = r["seq"] if r["fwd"] else revcomp(r["seq"])
                qual = r["qual"] if r["fwd"] else r["qual"][::-1]
                pos = mate["pos"] if mate["aln"] else -1
                cig, mq = None, 0
            if mate["aln"]:
                if mate["rev"]:
                    flag |= 0x20
                tags += [("MQ", mate["mapq"], "i"), ("MC", mate["cigar"], "Z")]
            supp = r["aln"]["supp"] if r["aln"] else []
            if supp:
                s = supp[0]
                srev = (not r["fwd"]) ^ s["reverse"]
                tags.append(("SA", f"{self.contig},{self.base + s['ref_start'] + 1},"
                                   f"{'-' if srev else '+'},{full_cigar(s, len(r['seq']))},{r['mapq']},0;", "Z"))
            if not r["aln"] and not mate["aln"]:
                return False
            rec = _rec(self.hdr, name, flag, self.tid, pos, mq, cig, seq, qual, tags)
            if mate["aln"]:
                rec.next_reference_id = self.tid
                rec.next_reference_start = mate["pos"]
            elif r["aln"]:
                rec.next_reference_id = self.tid
                rec.next_reference_start = pos
            if both:
                rec.template_length = tlen if r["pos"] <= mate["pos"] else -tlen
            self.records.append(rec)
            if supp:
                s = supp[0]
                srev = (not r["fwd"]) ^ s["reverse"]
                q0, q1 = s["q0"], s["q1"]
                sseq = r["seq"][q0:q1]; squal = r["qual"][q0:q1]
                if s["reverse"]:
                    sseq, squal = revcomp(sseq), squal[::-1]
                p = r["aln"]["primary"]
                prev = (not r["fwd"]) ^ p["reverse"]
                sa_back = (f"{self.contig},{r['pos'] + 1},{'-' if prev else '+'},{r['cigar']},{r['mapq']},0;")
                sflag = (flag & ~0x4) | 0x800
                sflag = (sflag & ~0x10) | (0x10 if srev else 0)
                srec = _rec(self.hdr, name, sflag, self.tid, self.base + s["ref_start"], r["mapq"],
                            full_cigar(s, len(r["seq"]), hard=True), sseq, squal,
                            [("SA", sa_back, "Z")])
                srec.next_reference_id = rec.next_reference_id
                srec.next_reference_start = rec.next_reference_start
                self.records.append(srec)
        return True


def _crosses(entries, b, k=10):
    lo, hi = entries_span(entries)
    return lo <= b - k and hi >= b + k


def window_for(ev, W=1500):
    ext = max(ev.target_len, 0)
    for e in ev.extras5:
        if e[0] == "REFCOPY":
            ext = max(ext, abs(e[2]) + e[3] + 50)
    N = 2 * W + 2 * ext + 200
    return N, W + ext + 100




class Prepared:
    """Sample-independent part of an event: the reference window (with planted EN motif /
    A-tract / reference element) and its haplotypes + truth breakpoints."""

    def __init__(self, ref, nick, strand, ev, haps, junctions, tr):
        self.ref, self.nick, self.strand, self.ev = ref, nick, strand, ev
        self.haps, self.junctions, self.tr = haps, junctions, tr


def _ref_hap(seq, boost=None):
    return Hap(seq, [Seg(0, len(seq), "REF", 0, len(seq), "+", "ref")], "ref", boost)


def prepare_event(rref, ev, strand, lib=None, W=1500):
    """Build the window + haplotypes for `ev` with the window RNG `rref` (deterministic, so
    every sample of a patient sees the same locus)."""
    N, nick = window_for(ev, W)
    ref = rnd_seq(rref, N)
    tr = {"strand": strand}
    if ev.plant_en:
        ref = plant_en_motif(ref, nick, strand, degenerate_en_motif(rref, ev.en_mm))
    blank = dict(ins_len=0, polya_side="-", tsd_seq="", en_motif=".", en_mm=".", x_seq="")
    if ev.render == "slippage":
        n = ev.info["tract_len"]
        ref = ref[:nick] + "A" * n + ref[nick + n:]
        ref = ref[:nick - 1] + ("C" if ref[nick - 1] == "A" else ref[nick - 1]) + ref[nick:]
        ref = ref[:nick + n] + ("G" if ref[nick + n] == "A" else ref[nick + n]) + ref[nick + n + 1:]
        tr.update(left=nick + n, right=nick, tsd=-n, **blank)
        return Prepared(ref, nick, strand, ev, {"ref": _ref_hap(ref, {nick: 3.0})}, {}, tr)
    if ev.render == "mismap":
        src = next(s for s in lib.sources["L1"] if s.id == ev.source_id)
        L = ev.info["ref_elem_len"]
        e_ref = src.element_seq[-L:] + src.tail       # reference L1 3' end ending at `nick`
        p = nick - len(e_ref)
        ref = ref[:p] + e_ref + ref[nick:]
        para = list(e_ref)                            # paralog copy, ~0.5% divergent
        for i in range(len(para)):
            if rref.random() < 0.005:
                para[i] = rref.choice([x for x in "ACGT" if x != para[i]])
        hseq = ref[:p] + "".join(para) + src.flank3[:1500]
        hap_m = Hap(hseq, [Seg(0, nick, "REF", 0, nick, "+", "ref"),
                           Seg(nick, len(hseq), "INS", label="paralog_flank")], "mis")
        tr.update(left=nick, right=nick, tsd=0, **blank)
        return Prepared(ref, nick, strand, ev, {"ref": _ref_hap(ref), "alt": hap_m}, {"R": nick}, tr)
    if ev.render == "foldback_reads":
        j = nick
        hseq = ref[:j] + revcomp(ref[j - 400:j])
        hap_f = Hap(hseq, [Seg(0, j, "REF", 0, j, "+", "ref"),
                           Seg(j, len(hseq), "REF", j - 400, j, "-", "foldback")], "fb")
        tr.update(left=j, right=j, tsd=0, **blank)
        return Prepared(ref, nick, strand, ev, {"ref": _ref_hap(ref), "alt": hap_f}, {"R": j}, tr)
    alt = apply_event(ref, nick, strand, ev)
    x0, x1 = alt.x_span
    tr.update(left=alt.left, right=alt.right, tsd=alt.right - alt.left, ins_len=len(alt.x_seq),
              polya_side=alt.polya_side, tsd_seq=alt.tsd_seq or "", en_motif=alt.en_motif,
              en_mm=alt.en_mm, x_seq=alt.x_seq, mh_seq=alt.mh_seq, parts=alt.parts)
    return Prepared(ref, nick, strand, ev, {"ref": _ref_hap(ref), "alt": Hap(alt.alt, alt.segs, "alt")},
                    {"R": x0, "L": x1}, tr)


def render_sample(rng, prep, writer, sampler, depth, vaf, present, name):
    """Sample one colony/sample's reads for a prepared event into writer.records.
    Returns per-junction support counts {R_frags, R_reads, L_frags, L_reads}."""
    ev = prep.ev
    counts = {"R_reads": 0, "R_frags": 0, "L_reads": 0, "L_frags": 0}

    def sample_into(hap, weight, junctions=None, n_override=None, starts=None, copies_fixed=None):
        n = n_override if n_override is not None else sampler.n_fragments(len(hap.seq), depth, weight)
        for i in range(n):
            fl = sampler.flen()
            st = starts(fl) if starts is not None else rng.randint(0, max(0, len(hap.seq) - fl))
            en = min(len(hap.seq), st + fl)
            ncopy = sampler.copies() if copies_fixed is None else copies_fixed
            r1_fwd = rng.random() < 0.5
            for c in range(ncopy + 1):
                pair = sampler.pair(hap, st, en)
                writer.write_pair(hap, pair, f"{name}_{hap.name}{i}" + (f"_d{c}" if c else ""), r1_fwd)
                for side, b in (junctions or {}).items():
                    if any(_crosses(r[0], b) for r in pair):
                        counts[side + "_reads"] += 1
                        if c == 0:
                            counts[side + "_frags"] += 1

    if ev.render == "slippage":
        old = (sampler.slip_min_delta, sampler.burst_scale)
        if present:      # exaggerated slippage + post-homopolymer phasing loss at this tract
            sampler.slip_min_delta, sampler.burst_scale = 3, max(2.5, 2.5 * sampler.burst_scale)
        hap = prep.haps["ref"] if present else _ref_hap(prep.ref)
        sample_into(hap, 1.0)
        sampler.slip_min_delta, sampler.burst_scale = old
        return counts
    if not present or (ev.role == "TP" and vaf <= 0):
        sample_into(prep.haps["ref"], 1.0)
        return counts
    if ev.render == "mismap":
        sample_into(prep.haps["ref"], 1.0)
        mq_old = writer.mapq_fn
        writer.mapq_fn = lambda: 0 if rng.random() < 0.7 else rng.choice([23, 37, 60])
        sample_into(prep.haps["alt"], 0.5, prep.junctions)
        writer.mapq_fn = mq_old
        return counts
    if ev.render == "foldback_reads":
        j = prep.junctions["R"]
        sample_into(prep.haps["ref"], 1.0)
        sample_into(prep.haps["alt"], 0, prep.junctions, n_override=ev.info["n_molecules"],
                    starts=lambda fl: rng.randint(max(0, j - fl + 40), j - 25),
                    copies_fixed=ev.info["pcr_copies"])
        return counts
    if ev.render == "single_fragment":
        sample_into(prep.haps["ref"], 1.0)
        sides = ev.info.get("sides", "POLYA")
        chosen = (["L" if prep.tr["polya_side"] == "LEFT" else "R"] if sides == "POLYA" else ["R", "L"])
        for side in chosen:
            b = prep.junctions[side]
            sample_into(prep.haps["alt"], 0, prep.junctions, n_override=1,
                        starts=lambda fl, b=b: rng.randint(max(0, b - fl + 30), max(0, b - 30)),
                        copies_fixed=ev.info.get("pcr_copies", 0))
        return counts
    # normal TP: VAF mixing between the reference and the alt haplotype
    sample_into(prep.haps["ref"], 1.0 - vaf)
    sample_into(prep.haps["alt"], vaf, prep.junctions)
    return counts
