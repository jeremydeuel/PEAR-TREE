#!/usr/bin/env python3
"""VAL-1 synthetic spike-in simulator (Phase 3+).

Emits a coordinate-sorted BAM of labelled retrotransposon insertions plus a truth
TSV, so discovery recall/precision can be measured across element classes, their
structural variants, TSD lengths and variant-allele fractions (VAF) — and, in the
same run, so *specificity* can be measured against the artefact catalogue of
`RTE_detection_review/05_sequencing_artefacts.md`.

Each insertion is rendered as the read signals discovery actually keys on — no
reference genome or aligner required, exactly as the differential generators do:

  * a LEFT breakpoint  = left-soft-clipped reads  (CIGAR `<clip>S<anchor>M`),
    breakpoint = reference_start = L  (the 3' / poly-A junction)
  * a RIGHT breakpoint = right-soft-clipped reads (CIGAR `<anchor>M<clip>S`),
    breakpoint = reference_end   = R = L + TSD  (the 5' junction)
  * output pairs the two when tsd_min <= R - L <= tsd_max, emitting `contig:L-R`

VAF is modelled by the alt-read count at each junction versus the count of
fully-aligned reference reads spanning the locus.

Geometry is calibrated to a real GRCm38 (mm10) mouse WGS CRAM (45521#13, PCR-free
NovaSeq, 1.8M reads profiled): 151 bp reads, soft-clip lengths drawn from the measured
distribution (p05/p50/p95 = 15/59/116 bp) so junction reads span each breakpoint at
realistic, varying offsets, ~320 bp insert size, and artefact-class prevalences noted
against each artefact block below (SMS 0.36%, cruciform 2.47%, maps-elsewhere 3.2%).

The truth `class` column encodes `<ELEMENT>_<variant>` (e.g. `L1HS_5p_inversion`,
`HERVK113_full_provirus`) so score.py reports recall per element *and* per variant.

What this simulator covers
--------------------------
1. **Element classes and every structural variant** (`RTE_detection_review` §2.3):
   LINE-1 (L1Hs), Alu (AluY) and SVA as TPRT/poly-A elements, each with full-length,
   5' truncation, 5' inversion (twin priming), partnered and orphan 3' transduction,
   3' (poly-A) deletion, and internal inv/del; ERV/LTR elements (generic ERVK/IAP-like)
   with full provirus, solo-LTR and 5'-inversion variants. The 5'/3' junction *clip
   content* is drawn from per-element consensus termini so a downstream `annotate`
   pass would classify them.
2. **Named human endogenous retroviruses HERV-K113 and HERV-K117** as full-length
   provirus insertions (6 bp TSD, LTR-consensus junctions, no poly-A). These are
   positive controls: HML-2 proviruses are not expected to retrotranspose de novo,
   but if one ever did, this is the signal PEAR-TREE would see. ⚠ The built-in LTR
   sequences are *structural stand-ins* (correct length, GC, TG..CA LTR boundary),
   not the real LTR5_Hs consensus — supply `--element-fasta` with the true
   HERV-K(HML-2) LTR / provirus termini to make them Dfam-annotatable.
3. **The artefact catalogue** (§5), each independently tuneable:
     - structure-specific chimeras: both-ends-soft-clipped (SMS) reads and self-fold
       palindromes (§5.1/5.2) — artefacts PEAR-TREE does *not* yet reject (review
       R15), so at baseline they surface as false positives a future SMS/self-fold
       filter must remove;
     - cruciform / inverted-repeat reads with an SA tag to the same contig within
       1 kb (§5.1) — PEAR-TREE *does* reject these (cluster poisoning), so they must
       stay absent from the calls;
     - poly-G / dark-cycle (§5.4), adapter read-through (§5.5), homopolymer /
       low-complexity anchors (§5.6), PCR duplicates and non-reproducible PCR chimeras
       (§5.7), and segdup/mismap reads that map fully elsewhere via an XA tag (§5.9) —
       all of which existing filters must remove (0 false positives expected);
     - high-coverage pile-up regions (§4.2.6) for the SPEC-3/4 coverage gates;
     - poly-A dropout (§5.3): a genuine L1 whose poly-A junction reads are lost to
       library prep — a real insertion baseline necessarily *misses*.
4. **Feature A/B toggles**: discordant-mate rescue of one-sided junctions (RTE-origin vs
   non-RTE) and processed-pseudogene mate signatures, with the matching RepeatMasker `.out`
   and exon annotation tracks.
5. **Novel processed pseudogenes** (`--n-novel-pseudogene`): de-novo L1-mediated retrocopies
   of a spliced mRNA. Unlike the Feature-B signature check above (which reuses an ERV clip and
   only exercises the mate-spanning-exons signal), these render a biologically faithful
   retrocopy — the inserted / clipped sequence *is* the spliced transcript (parent-gene exons
   concatenated, introns skipped) carrying the L1 TPRT scar (poly-A tail + TSD + EN motif). So
   each is a real insertion (in truth, scored for recall/VAF) whose junction reads *also* carry
   the splice hallmark (mates spanning >=2 parent-gene exons), exercising the whole
   discovery -> splice-annotate path on realistic sequence. The parent-gene exons are written to
   the `--out-exons` track (gene `PG1`).

*** SYNTHETIC-ONLY LIMITATION ***
The reads are still clean and the junctions exact; a synthetic artefact exercises the
*code path* a real artefact would, but cannot reproduce real enzymatic-prep messiness.
Green VAL-1 remains necessary, not sufficient: orthogonal (IGV / long-read / PCR)
confirmation of real calls is still required.
"""
import argparse
import os
import random
import sys

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from simlib import library as _simlib_library  # noqa: E402
from simlib import models as _simlib_models  # noqa: E402
from simlib import truth as _simlib_truth  # noqa: E402
from simlib.reads import FragmentSampler  # noqa: E402
from simlib.seqs import rnd_seq as _rnd_seq  # noqa: E402
from simlib.val1render import PairWriter, prepare_event, render_sample  # noqa: E402

# a fixed "LTR consensus" shared by the SENS/Feature ERV toggles below — models the
# biology that ERV clips are the element LTR consensus (so many insertions share one
# clip). Kept verbatim so those validated toggle targets are byte-stable.
LTR_CONSENSUS = "TGCTAGGCAACCGTATTCAGGTACCGATTCAGGCATTGACCGATTCAGGTACCGATTCAGG"
ANCHOR_POOL = "ACGTTGCACCGATTACGGATCCGTTAAGGCACTTGACCGTATCGGATCACGTTAGCCGATA"

BASES = "ACGT"
_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def rnd_seq(rng, n):
    return "".join(rng.choice(BASES) for _ in range(n))


def clip_from(seq, off, length, rc=False):
    """A `length`-bp fragment of `seq` starting at `off` (cycling if short)."""
    if not seq:
        seq = "A"
    reps = (off + length) // len(seq) + 2
    s = seq * reps
    frag = s[off:off + length]
    if len(frag) < length:
        frag = (frag + s)[:length]
    return revcomp(frag) if rc else frag


# ---------------------------------------------------------------------------
# Element consensus library.
#
# For a genome-free, clip-signal simulator, the *breakpoint geometry* (TSD, poly-A
# presence, split junction) drives discovery recall, while the *clip sequence* matters
# only for a downstream annotate/positive-control step. The sequences below are
# deterministic structural stand-ins (fixed-seed RNG) with the right hallmarks — LTRs
# open `TG` and close `CA`; TPRT elements carry a poly-A tail and an L1 endonuclease
# `TTAAAA` motif in the target-site flank. Replace any of them with real consensus via
# --element-fasta (records `<ELEMENT>_end5`, `<ELEMENT>_end3`, `<ELEMENT>_body`, or
# `<ELEMENT>_ltr`).
# ---------------------------------------------------------------------------
_lib = random.Random(0xE1E5)


def _b(n):
    return "".join(_lib.choice(BASES) for _ in range(n))


def _ltr(n=150):
    return "TG" + _b(n - 4) + "CA"          # retroviral LTR hallmark: 5'-TG..CA-3'


def _tprt(end5_len=150, end3_len=130, body_len=420):
    return {"end5": _b(end5_len), "end3": _b(end3_len), "body": _b(body_len), "polya": True}


def _ltr_elem(ltr_len=150):
    ltr = _ltr(ltr_len)
    return {"end5": ltr, "end3": ltr, "body": ltr + _b(120) + ltr, "polya": False}


def default_elements():
    return {
        # human TPRT poly-A elements
        "L1HS": {**_tprt(), "klass": "LINE", "tsd": (5, 20)},
        "ALUY": {**_tprt(140, 110, 280), "klass": "SINE", "tsd": (5, 20)},
        "SVA":  {**_tprt(150, 130, 360), "klass": "SVA", "tsd": (5, 20)},
        # LTR / ERV elements (integrase-mediated: short TSD, no poly-A)
        "ERVK": {**_ltr_elem(150), "klass": "ERV", "tsd": (4, 6)},
        # named human HML-2 proviruses — positive controls (not expected to jump)
        "HERVK113": {**_ltr_elem(160), "klass": "ERV", "tsd": (6, 6)},
        "HERVK117": {**_ltr_elem(160), "klass": "ERV", "tsd": (6, 6)},
    }


def read_fasta(path):
    recs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    recs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        recs[name] = "".join(buf).upper()
    return recs


def apply_element_fasta(elements, path):
    """Override element termini from a FASTA (`<ELEMENT>_end5/_end3/_body/_ltr`)."""
    for key, seq in read_fasta(path).items():
        if "_" not in key:
            continue
        elem, part = key.rsplit("_", 1)
        if elem not in elements:
            continue
        if part == "ltr":
            elements[elem]["end5"] = elements[elem]["end3"] = seq
            elements[elem]["body"] = seq + seq
        elif part in ("end5", "end3", "body"):
            elements[elem][part] = seq
    return elements


# TPRT and LTR structural variants. Each maps to a pair of junction clips in
# build_junctions(); the truth `class` is `<ELEMENT>_<variant>`.
L1_VARIANTS = ["5p_truncated", "5p_inversion", "3p_transduction_partnered",
               "3p_transduction_orphan", "5p_transduction", "3p_deletion", "internal_inv_del"]
ALU_VARIANTS = ["full", "5p_truncated", "5p_inversion",
                "3p_transduction_partnered", "3p_transduction_orphan", "3p_deletion"]
SVA_VARIANTS = ["full", "5p_truncated", "5p_inversion", "3p_deletion"]
ERV_VARIANTS = ["solo_ltr", "5p_inversion"]

# Novel processed-pseudogene parent-gene exon sequences. A de-novo retrocopy carries the
# *spliced transcript* — its exons concatenated, introns skipped — so the inserted (and
# hence clipped) sequence is this mRNA. Fixed-seed so that mRNA (and thus the parent-gene
# identity a downstream annotate would assign) is byte-stable across runs, independent of
# --seed. Four exons; the spliced mRNA is their concatenation. Exon length exceeds the read
# length so a mate landing in an exon is a full within-exon alignment.
_pg_rng = random.Random(0x9E0D)
PG_EXON_LEN = 300
PG_N_EXONS = 4
PG_EXON_SEQ = ["".join(_pg_rng.choice(BASES) for _ in range(PG_EXON_LEN)) for _ in range(PG_N_EXONS)]


# Insert-terminus length made available for clipping. Reads span a breakpoint at
# varying offsets (see sample_clip_len), so a read may clip up to this many bp of
# inserted sequence; must exceed the max sampled clip length.
INS_LEN = 140

# Empirical soft-clip-length inverse-CDF, measured on a real GRCm38 mouse WGS CRAM
# (45521#13, 1.8M reads: usable one-sided clips p05/p25/p50/p75/p95 = 15/31/59/92/116).
_CLIP_ICDF = [(0.0, 15), (0.05, 15), (0.25, 31), (0.50, 59), (0.75, 92), (0.95, 116), (1.0, 130)]


def sample_clip_len(rng, lo=15, hi=110):
    """Draw a soft-clip length from the real distribution, clamped so every junction
    read stays a valid evidence read (clip >= min_clip_len, anchor = read_len - clip
    >= the 40 bp unclipped floor)."""
    u = rng.random()
    for (p0, v0), (p1, v1) in zip(_CLIP_ICDF, _CLIP_ICDF[1:]):
        if u <= p1:
            v = v0 + (v1 - v0) * (u - p0) / (p1 - p0) if p1 > p0 else v0
            return int(min(max(v, lo), hi))
    return hi


def build_junctions(rng, elem, variant, tail, ins_len=INS_LEN):
    """Return (ins5, ins3): the inserted sequence available at the 5' (RIGHT/R) and
    3' (LEFT/L) junctions. ins5[0] abuts R and extends into the insert; ins3[-1] abuts
    L (so for a TPRT element the poly-A tail is the L-adjacent suffix, exactly where a
    read barely crossing the 3' junction would clip it). Each is `ins_len` bp."""
    end5, end3, body = elem["end5"], elem["end3"], elem["body"]
    tailA = "A" * tail

    if not elem["polya"]:                                   # LTR / ERV: both ends LTR
        if variant == "5p_inversion":
            return clip_from(body, 20, ins_len, rc=True), clip_from(end3, 0, ins_len)
        return clip_from(end5, 0, ins_len), clip_from(end3, 0, ins_len)

    # TPRT poly-A element. ins3 ends in the poly-A tail (abuts the 3' junction L).
    i3 = clip_from(end3, 0, ins_len - tail) + tailA
    if variant == "full":
        return clip_from(end5, 0, ins_len), i3
    if variant == "5p_truncated":                           # ~51% of somatic L1s
        return clip_from(body, 60, ins_len), i3
    if variant == "5p_inversion":                           # twin priming (~30%)
        return clip_from(body, 45, ins_len, rc=True), i3
    if variant == "3p_transduction_partnered":              # body + unique tag + polyA
        i3 = clip_from(end3, 0, 20) + rnd_seq(rng, ins_len - 20 - tail) + tailA
        return clip_from(end5, 0, ins_len), i3
    if variant == "3p_transduction_orphan":                 # body lost: unique tag + polyA
        return rnd_seq(rng, ins_len), rnd_seq(rng, ins_len - tail) + tailA
    if variant == "5p_transduction":                        # unique tag upstream of 5' end
        return rnd_seq(rng, ins_len), i3
    if variant == "3p_deletion":                            # poly-A end lost
        return clip_from(end5, 0, ins_len), clip_from(body, 90, ins_len)
    if variant == "internal_inv_del":                       # canonical ends
        return clip_from(end5, 0, ins_len), i3
    return clip_from(end5, 0, ins_len), i3


def make_read(hdr, tid, name, seq, start, cigar, mapq, flag=0,
              next_tid=-1, next_start=-1, tlen=0, tags=None):
    a = pysam.AlignedSegment(hdr)
    a.query_name = name
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    a.flag = flag
    a.reference_id = tid
    a.reference_start = start
    a.mapping_quality = mapq
    a.cigarstring = cigar
    # mate/discordant fields (Feature A/B): PNEXT/RNEXT + template length. A read is
    # "discordant" when it is paired (0x1) but not proper (no 0x2); its mate is mapped
    # (no 0x8) either on another contig or far away on the same one.
    a.next_reference_id = next_tid
    a.next_reference_start = next_start
    a.template_length = tlen
    if tags:
        a.set_tags(tags)
    return a


PAIRED_R1 = 0x1 | 0x40


def simulate(args):
    rng = random.Random(args.seed)
    contigs = ["1", "2", "3"]
    hdr = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": args.contig_len} for c in contigs],
    })
    tid = {c: i for i, c in enumerate(contigs)}

    elements = default_elements()
    if args.element_fasta:
        apply_element_fasta(elements, args.element_fasta)

    anchor_m = args.anchor_len   # aligned bases for the fixed-geometry probe reads (SENS/artefacts)
    clip_s = args.clip_len       # clipped bases for the fixed-geometry probe reads
    read_len = args.read_len     # 151 in the calibrating mouse CRAM
    ref_m = read_len             # fully-aligned reference read length

    records = []
    truth = []
    # space insertions far apart so their breakpoints never cross-cluster
    step = args.step
    slot = [10_000, 10_000, 10_000]  # next free position per contig
    idx = 0

    # shared ERV clip + anchors for the SENS/Feature toggle targets (byte-stable)
    erv5 = LTR_CONSENSUS[:clip_s].ljust(clip_s, "A")
    anchor_l = (ANCHOR_POOL * 2)[:anchor_m]
    anchor_r = (ANCHOR_POOL[::-1] * 2)[:anchor_m]

    def sample_insert():
        # proper-pair insert size, calibrated to the CRAM (median 319, p05-p95 196-535)
        return max(read_len + 20, int(rng.gauss(330, 90)))

    def coverage_reads(tag, ti, L, n=None):
        if n is None:
            n = rng.randint(args.ref_min, args.ref_max)
        for k in range(n):
            ins = sample_insert()
            st = max(L - ref_m // 2 + (k % 7), 0)
            records.append(make_read(hdr, ti, f"{tag}_ref{k}", rnd_seq(rng, ref_m),
                                     st, f"{ref_m}M", 60,
                                     flag=0x1 | 0x2 | 0x40, next_tid=ti,
                                     next_start=st + ins - ref_m, tlen=ins))
        return n

    def make_anchor(en_motif=False):
        a = rnd_seq(rng, anchor_m)
        if en_motif:                       # L1 endonuclease target motif in the flank
            a = "TTAAAA" + rnd_seq(rng, anchor_m - 6)
        return a

    def make_flank(en_motif=False):
        a = rnd_seq(rng, read_len)
        if en_motif:                       # L1 endonuclease target motif abutting L
            a = "TTAAAA" + rnd_seq(rng, read_len - 6)
        return a

    def junction_mapq():
        # most junction reads are MAPQ 60; a real fraction anchor in a repeat and drop
        # to 0 (12% in the CRAM). Off by default so the recall gate stays deterministic.
        if args.junction_lowmapq_frac and rng.random() < args.junction_lowmapq_frac:
            return 0
        return 60

    def emit_insertion(element, variant):
        """One real insertion: 5' (RIGHT@R) + 3' (LEFT@L) junctions + coverage. Junction
        reads span the breakpoint at varying offsets, so each clips a different amount of
        inserted sequence (clip length ~ the empirical CRAM distribution)."""
        nonlocal idx
        elem = elements[element]
        contig = contigs[idx % len(contigs)]
        ti = tid[contig]
        L = slot[ti]; slot[ti] += step
        lo, hi = elem["tsd"]
        tsd = rng.randint(lo, hi)
        R = L + tsd
        alt = rng.randint(args.alt_min, args.alt_max)
        ref = rng.randint(args.ref_min, args.ref_max)
        vaf = round(alt / (alt + ref), 3) if (alt + ref) else 0.0

        # true poly-A tails vary (CRAM proximal A/T runs: median 7, but those are mostly
        # background; a *detectable* insertion tail must clear the 12 bp POLYA_CUTOFF, so
        # sample 13..polya_len). Off-tail (ERV) insertions ignore this.
        t = rng.randint(13, max(14, args.polya_len))
        ins5, ins3 = build_junctions(rng, elem, variant, t)
        g5 = make_flank()                          # genomic flank 5' of R
        g3 = make_flank(en_motif=elem["polya"])    # genomic flank 3' of L (poly-A side)
        tag = f"{element}_{variant}_{contig}_{L}"

        # A transduced non-RTE "tag" is co-mobilised from a unique SOURCE locus; reads
        # crossing the tag-bearing junction have their mates there (the barcode that
        # traces the event to its master element — Tubio 2014). This is exactly what a
        # transduction shares with a translocation-at-RTE; the poly-A + TSD + EN motif are
        # what set the true retrotransposition apart, so a transduction/translocation
        # discriminator must key on the hallmarks, not merely on the cross-locus pointer.
        tag_side = {"3p_transduction_partnered": "L", "3p_transduction_orphan": "R",
                    "5p_transduction": "R"}.get(variant)
        src = ((ti + 1) % len(contigs), 14_000_000 + idx * 500) if tag_side else (-1, -1)
        mR = src if tag_side == "R" else (-1, -1)
        mL = src if tag_side == "L" else (-1, -1)

        # RIGHT breakpoint at R (5' junction): reads = [genome ...R][insert], clip = insert
        for k in range(alt):
            clip = sample_clip_len(rng); a = read_len - clip
            records.append(make_read(hdr, ti, f"{tag}_R{k}", g5[-a:] + ins5[:clip],
                                     R - a, f"{a}M{clip}S", junction_mapq(),
                                     flag=PAIRED_R1, next_tid=mR[0], next_start=mR[1]))
        # LEFT breakpoint at L (3' junction): reads = [insert][genome L...], clip = insert
        for k in range(alt):
            clip = sample_clip_len(rng); a = read_len - clip
            records.append(make_read(hdr, ti, f"{tag}_L{k}", ins3[-clip:] + g3[:a],
                                     L, f"{clip}S{a}M", junction_mapq(),
                                     flag=PAIRED_R1, next_tid=mL[0], next_start=mL[1]))
        coverage_reads(tag, ti, L, ref)
        truth.append((contig, L, R, f"{element}_{variant}", tsd, alt, ref, vaf))
        idx += 1

    # --- element plan (backward-compatible --n-l1 / --n-erv drive canonical) --------
    plan = [("L1HS", "full", args.n_l1), ("ERVK", "full_provirus", args.n_erv)]
    plan += [("L1HS", v, args.n_l1_variant) for v in L1_VARIANTS]
    plan += [("ALUY", v, args.n_alu) for v in ALU_VARIANTS]
    plan += [("SVA", v, args.n_sva) for v in SVA_VARIANTS]
    plan += [("ERVK", v, args.n_erv_variant) for v in ERV_VARIANTS]
    plan += [("HERVK113", "full_provirus", args.n_hervk113),
             ("HERVK117", "full_provirus", args.n_hervk117)]
    for element, variant, n in plan:
        for _ in range(n):
            emit_insertion(element, variant)

    # SENS-1 targets: ERV insertions whose junction breakpoints wobble over a few bp,
    # so the exact modal count is 1 (missed) but the windowed count clears the floor.
    t_right = anchor_r + erv5   # [anchor | element], boundary slides for RIGHT clips
    t_left = erv5 + anchor_l    # [element | anchor], boundary slides for LEFT clips
    for _ in range(args.n_wobble):
        contig = contigs[idx % len(contigs)]; ti = tid[contig]
        L = slot[ti]; slot[ti] += step
        tsd = rng.randint(args.tsd_min, args.tsd_max); R = L + tsd
        alt = max(args.alt_min, 3)
        tag = f"WOB_{contig}_{L}"
        for k in range(alt):
            off = k - alt // 2  # spread breakpoints across +/- a couple bp
            records.append(make_read(hdr, ti, f"{tag}_R{k}", t_right,
                                     R - anchor_m, f"{anchor_m + off}M{clip_s - off}S", 60, flag=PAIRED_R1))
            records.append(make_read(hdr, ti, f"{tag}_L{k}", t_left,
                                     L + off, f"{clip_s + off}S{anchor_m - off}M", 60, flag=PAIRED_R1))
        coverage_reads(tag, ti, L)
        truth.append((contig, L, R, "WOBBLE", tsd, alt, 0, 0.0))
        idx += 1

    # SENS-2 targets: ERV insertions with one high-MAPQ + one low-MAPQ read per junction.
    for _ in range(args.n_lowmapq):
        contig = contigs[idx % len(contigs)]; ti = tid[contig]
        L = slot[ti]; slot[ti] += step
        tsd = rng.randint(args.tsd_min, args.tsd_max); R = L + tsd
        tag = f"LMQ_{contig}_{L}"
        for mq in (60, 30):
            records.append(make_read(hdr, ti, f"{tag}_R{mq}", anchor_r + erv5,
                                     R - anchor_m, f"{anchor_m}M{clip_s}S", mq, flag=PAIRED_R1))
            records.append(make_read(hdr, ti, f"{tag}_L{mq}", erv5 + anchor_l,
                                     L, f"{clip_s}S{anchor_m}M", mq, flag=PAIRED_R1))
        coverage_reads(tag, ti, L)
        truth.append((contig, L, R, "LOWMAPQ", tsd, 2, 0, 0.0))
        idx += 1

    # SENS-8 targets: an insertion whose RIGHT junction clip is a short pure poly-A
    # tail (>= min_good_bases but <= the 12 bp floor), so baseline rejects that side.
    short = "A" * args.shortpolya_len
    for _ in range(args.n_shortpolya):
        contig = contigs[idx % len(contigs)]; ti = tid[contig]
        L = slot[ti]; slot[ti] += step
        tsd = rng.randint(args.tsd_min, min(args.tsd_max, 20)); R = L + tsd
        alt = max(args.alt_min, 2)
        tag = f"SPA_{contig}_{L}"
        for k in range(alt):
            records.append(make_read(hdr, ti, f"{tag}_R{k}", anchor_r + short,
                                     R - anchor_m, f"{anchor_m}M{len(short)}S", 60, flag=PAIRED_R1))
            records.append(make_read(hdr, ti, f"{tag}_L{k}", erv5 + anchor_l,
                                     L, f"{clip_s}S{anchor_m}M", 60, flag=PAIRED_R1))
        coverage_reads(tag, ti, L)
        truth.append((contig, L, R, "SHORTPOLYA", tsd, alt, 0, 0.0))
        idx += 1

    # §5.3 poly-A dropout: a genuine L1 whose 3' poly-A junction reads are lost to the
    # library prep, leaving only the 5' side. Baseline cannot pair it -> MISS.
    l1 = elements["L1HS"]
    for _ in range(args.n_polya_dropout):
        contig = contigs[idx % len(contigs)]; ti = tid[contig]
        L = slot[ti]; slot[ti] += step
        tsd = rng.randint(l1["tsd"][0], l1["tsd"][1]); R = L + tsd
        alt = rng.randint(args.alt_min, args.alt_max)
        clip5 = clip_from(l1["end5"], 0, clip_s)
        tag = f"PAD_{contig}_{L}"
        for k in range(alt):
            records.append(make_read(hdr, ti, f"{tag}_R{k}", make_anchor() + clip5,
                                     R - anchor_m, f"{anchor_m}M{clip_s}S", 60, flag=PAIRED_R1))
        coverage_reads(tag, ti, L)
        truth.append((contig, L, R, "L1HS_polya_dropout", tsd, alt, 0, 0.0))
        idx += 1

    # MATE-RESCUE targets: a real L1 insertion into a low-mapability-but-mate-unique
    # flank. The junction clipped reads anchor in a repeat, so bwa gives their OWN
    # alignment a MAPQ below the floor (materescue_mapq), but each is a proper pair whose
    # mate maps uniquely (MQ tag = 60). Baseline discovery drops them (mapq<min_mapq);
    # `mate_anchor_rescue` keeps them because the unique mate pins the locus. This is the
    # fullstack id=32 scenario (L1HS 3'-transduction at a low-mapability hg38 homolog).
    for _ in range(args.n_materescue):
        contig = contigs[idx % len(contigs)]; ti = tid[contig]
        L = slot[ti]; slot[ti] += step
        tsd = rng.randint(args.tsd_min, min(args.tsd_max, 20)); R = L + tsd
        alt = rng.randint(max(args.alt_min, 2), args.alt_max)
        t = rng.randint(13, max(14, args.polya_len))
        ins5, ins3 = build_junctions(rng, elements["L1HS"], "full", t)
        g5 = make_flank(); g3 = make_flank(en_motif=True)
        tag = f"MRESC_{contig}_{L}"
        mq_tag = [("MQ", 60, "i")]                 # mate maps uniquely
        for k in range(alt):                       # RIGHT junction @R, own MAPQ low
            clip = sample_clip_len(rng); a = read_len - clip
            ins = sample_insert(); mst = max(R - ins, 0)
            records.append(make_read(hdr, ti, f"{tag}_R{k}", g5[-a:] + ins5[:clip],
                                     R - a, f"{a}M{clip}S", args.materescue_mapq,
                                     flag=0x1 | 0x2 | 0x40, next_tid=ti, next_start=mst,
                                     tlen=ins, tags=mq_tag))
        for k in range(alt):                       # LEFT junction @L, own MAPQ low
            clip = sample_clip_len(rng); a = read_len - clip
            ins = sample_insert(); mst = L + a + ins - ref_m
            records.append(make_read(hdr, ti, f"{tag}_L{k}", ins3[-clip:] + g3[:a],
                                     L, f"{clip}S{a}M", args.materescue_mapq,
                                     flag=0x1 | 0x2 | 0x40 | 0x20, next_tid=ti,
                                     next_start=mst, tlen=ins, tags=mq_tag))
        coverage_reads(tag, ti, L)
        truth.append((contig, L, R, "L1HS_materescue", tsd, alt, 0, 0.0))
        idx += 1

    # ------------------------------------------------------------------ Feature A
    # DISCORDANT targets (A1): a one-sided junction — a real LEFT soft-clip at L, no
    # RIGHT soft-clip — rescued by a cluster of reverse discordant anchor reads that
    # start at ~R (reference_start = R => RIGHT-role partner) whose mates map into an
    # element locus. `discordant` variants send mates to the RTE band (contig 3,
    # covered by the emitted rmsk track => RTE-origin true); `disc_artefact` variants
    # send mates to a non-RTE band. With `discordant_anchor` on (no gate) BOTH are
    # rescued; under `discordant_rte_only` only the RTE-origin ones survive.
    elem_tid = tid[contigs[2]]  # contig "3" holds the synthetic element / gene loci
    RTE_BAND = 15_000_000       # mates here fall inside the emitted rmsk track (young)
    NONRTE_BAND = 16_000_000    # mates here are random genome (no rmsk entry)
    rmsk_rows = []  # (contig, begin, end, div) young RepeatMasker element intervals

    def discordant_block(n, mate_band, cls, in_truth):
        nonlocal idx
        for j in range(n):
            contig = contigs[idx % len(contigs)]
            if contig == contigs[2]:
                contig = contigs[0]  # keep the insertion off the element contig
            ti = tid[contig]
            L = slot[ti]; slot[ti] += step
            tsd = rng.randint(args.tsd_min, min(args.tsd_max, 20)); R = L + tsd
            alt = max(args.alt_min, 2)
            ndisc = max(args.discordant_min_reads, 3)
            tag = f"{cls}_{contig}_{L}"
            # real LEFT breakpoint at L
            for k in range(alt):
                records.append(make_read(hdr, ti, f"{tag}_L{k}", erv5 + anchor_l,
                                         L, f"{clip_s}S{anchor_m}M", 60, flag=PAIRED_R1))
            # reverse discordant anchors at R, mate -> element band on contig 3
            for k in range(ndisc):
                mpos = mate_band + j * 2000 + k * 10
                records.append(make_read(hdr, ti, f"{tag}_D{k}", rnd_seq(rng, ref_m),
                                         R, f"{ref_m}M", 60,
                                         flag=0x1 | 0x10 | 0x80,  # paired, reverse, read2, NOT proper
                                         next_tid=elem_tid, next_start=mpos))
            coverage_reads(tag, ti, L)
            if in_truth:
                truth.append((contig, L, R, cls.upper(), tsd, alt, 0, 0.0))
            idx += 1

    discordant_block(args.n_discordant, RTE_BAND, "discordant", in_truth=True)
    discordant_block(args.n_disc_artefact, NONRTE_BAND, "disc_artefact", in_truth=False)

    # Breakpoint in the unsequenced insert GAP (between read1 and read2): NO read spans
    # either junction with a soft-clip, so the insertion is supported ONLY by discordant
    # read pairs — a genomic read on each flank whose mate maps into the element. This is
    # PEAR-TREE's architectural blind spot (§7.7#5: no discordant-pair discovery), so it is
    # MISSED at baseline AND by the one-sided-clip Feature-A rescue (which still needs a
    # clip). In truth, opt-in; the target for validating a future discordant-pair caller.
    for _ in range(args.n_discordant_only):
        contig = contigs[idx % len(contigs)]
        if contig == contigs[2]:
            contig = contigs[0]                    # keep the site off the element contig
        ti = tid[contig]
        L = slot[ti]; slot[ti] += step
        tsd = rng.randint(args.tsd_min, min(args.tsd_max, 20)); R = L + tsd
        ndisc = max(args.discordant_min_reads, 3)
        tag = f"discordant_only_{contig}_{L}"
        for k in range(ndisc):
            # LEFT flank read ends at L (no clip); mate maps into the element band
            fl = max(L - ref_m - k * 13, 0); ml = RTE_BAND + k * 40
            records.append(make_read(hdr, ti, f"{tag}_FL{k}", rnd_seq(rng, ref_m), fl,
                                     f"{ref_m}M", 60, flag=0x1 | 0x40, next_tid=elem_tid, next_start=ml))
            records.append(make_read(hdr, elem_tid, f"{tag}_FL{k}", rnd_seq(rng, ref_m), ml,
                                     f"{ref_m}M", 60, flag=0x1 | 0x80 | 0x10, next_tid=ti, next_start=fl))
            # RIGHT flank read starts at R (no clip); mate maps into the element band
            fr = R + k * 13; mr = RTE_BAND + 2000 + k * 40
            records.append(make_read(hdr, ti, f"{tag}_FR{k}", rnd_seq(rng, ref_m), fr,
                                     f"{ref_m}M", 60, flag=0x1 | 0x40 | 0x10, next_tid=elem_tid, next_start=mr))
            records.append(make_read(hdr, elem_tid, f"{tag}_FR{k}", rnd_seq(rng, ref_m), mr,
                                     f"{ref_m}M", 60, flag=0x1 | 0x80, next_tid=ti, next_start=fr))
        coverage_reads(tag, ti, L)
        truth.append((contig, L, R, "DISCORDANT_ONLY", tsd, ndisc, 0, 0.0))
        idx += 1

    if args.n_discordant or args.n_disc_artefact or args.n_discordant_only:
        # one young-element interval spanning the whole RTE band the mates map into
        rmsk_rows.append((contigs[2], RTE_BAND, RTE_BAND + 1_000_000, 3.0))

    # ------------------------------------------------------------------ Feature B
    # PSEUDOGENE targets (B3): a normal insertion whose RIGHT-junction reads are read1
    # of pairs whose mates map into >=2 exons of a single synthetic gene (introns
    # skipped) — the processed-pseudogene signature. A single-exon variant must NOT be
    # flagged. Exons are emitted to the companion annotation file.
    GENE_START = 17_000_000
    EXON_LEN = 400
    EXON_GAP = 10_000  # intron between exons
    exon_rows = []  # (contig, begin, end, gene_id)
    n_gene_exons = 3
    for e in range(n_gene_exons):
        b = GENE_START + e * EXON_GAP
        exon_rows.append((contigs[2], b, b + EXON_LEN, "G1"))

    def pseudogene_block(n, n_exons_hit, cls, in_truth):
        nonlocal idx
        for _ in range(n):
            contig = contigs[idx % len(contigs)]
            if contig == contigs[2]:
                contig = contigs[0]
            ti = tid[contig]
            L = slot[ti]; slot[ti] += step
            tsd = rng.randint(args.tsd_min, min(args.tsd_max, 20)); R = L + tsd
            alt = max(args.alt_min, n_exons_hit)
            tag = f"{cls}_{contig}_{L}"
            # LEFT breakpoint (plain) at L
            for k in range(alt):
                records.append(make_read(hdr, ti, f"{tag}_L{k}", erv5 + anchor_l,
                                         L, f"{clip_s}S{anchor_m}M", 60, flag=PAIRED_R1))
            # RIGHT breakpoint at R: forward read1 => has_mate; mate maps into exon (k mod n_exons_hit)
            for k in range(alt):
                exon = k % n_exons_hit
                mpos = GENE_START + exon * EXON_GAP + 50
                nm = f"{tag}_R{k}"
                records.append(make_read(hdr, ti, nm, anchor_r + erv5,
                                         R - anchor_m, f"{anchor_m}M{clip_s}S", 60,
                                         flag=PAIRED_R1,  # paired, read1, forward
                                         next_tid=elem_tid, next_start=mpos))
                # the mate as a real record inside the exon (read2)
                records.append(make_read(hdr, elem_tid, nm, rnd_seq(rng, ref_m),
                                         mpos, f"{ref_m}M", 60,
                                         flag=0x1 | 0x80, next_tid=ti, next_start=R - anchor_m))
            coverage_reads(tag, ti, L)
            if in_truth:
                truth.append((contig, L, R, cls.upper(), tsd, alt, 0, 0.0))
            idx += 1

    pseudogene_block(args.n_pseudogene, n_gene_exons, "pseudogene", in_truth=True)
    pseudogene_block(args.n_single_exon, 1, "single_exon", in_truth=True)

    # NOVEL PROCESSED PSEUDOGENE (de-novo retrocopy). Where --n-pseudogene above reuses an
    # ERV clip and only exercises the mate-spanning-exons *signal*, this renders a
    # biologically faithful retrocopy: the inserted (hence clipped) sequence IS the spliced
    # transcript — the parent gene's exons concatenated, introns skipped — and it carries the
    # L1 TPRT scar (poly-A tail + TSD + EN-motif flank). The review's processed pseudogene:
    # "the same poly-A/TSD hallmarks but carry spliced exonic sequence" (§2.3). So it is a
    # REAL insertion (in truth, scored for recall/VAF like any TPRT element) whose junction
    # reads ALSO carry the splice hallmark — their mates map into >=2 exons of the parent gene
    # (introns skipped) — exercising the whole discovery -> splice-annotate path on realistic
    # sequence rather than an ERV stand-in. Parent-gene exons go into the --out-exons track.
    PG_GENE_START = 18_000_000                 # parent gene: on the annotated contig 3, clear of G1 (17 Mb)
    PG_EXON_GAP = 8_000                        # intron between exons (>> exon length => intron-skip test passes)
    for e in range(PG_N_EXONS):
        b = PG_GENE_START + e * PG_EXON_GAP
        exon_rows.append((contigs[2], b, b + PG_EXON_LEN, "PG1"))
    pg_mrna = "".join(PG_EXON_SEQ)
    # a TPRT-style element whose termini are the mRNA ends (5' = transcript start, 3' = the
    # bases abutting the poly-A tail), so build_junctions renders exonic clips + a poly-A tail.
    pg_elem = {"end5": pg_mrna[:INS_LEN], "end3": pg_mrna[-INS_LEN:], "body": pg_mrna, "polya": True}
    for _ in range(args.n_novel_pseudogene):
        contig = contigs[idx % len(contigs)]
        if contig == contigs[2]:
            contig = contigs[0]                # keep the insertion off the parent-gene contig
        ti = tid[contig]
        L = slot[ti]; slot[ti] += step
        tsd = rng.randint(5, 20); R = L + tsd
        # enough alt reads to span every exon (mate exon = k % PG_N_EXONS), and to clear the
        # >= splice_min_exons floor even at the low end.
        alt = max(rng.randint(args.alt_min, args.alt_max), PG_N_EXONS)
        ref = rng.randint(args.ref_min, args.ref_max)
        vaf = round(alt / (alt + ref), 3) if (alt + ref) else 0.0
        t = rng.randint(13, max(14, args.polya_len))
        ins5, ins3 = build_junctions(rng, pg_elem, "full", t)
        g5 = make_flank(); g3 = make_flank(en_motif=True)   # L1 EN motif on the poly-A (3') flank
        tag = f"NPSG_{contig}_{L}"
        # RIGHT (5') junction @R: read1 forward, clip = spliced-mRNA 5' end; mate -> a parent exon
        for k in range(alt):
            clip = sample_clip_len(rng); a = read_len - clip
            exon = k % PG_N_EXONS
            mpos = PG_GENE_START + exon * PG_EXON_GAP + 40
            nm = f"{tag}_R{k}"
            records.append(make_read(hdr, ti, nm, g5[-a:] + ins5[:clip],
                                     R - a, f"{a}M{clip}S", junction_mapq(),
                                     flag=PAIRED_R1, next_tid=elem_tid, next_start=mpos))
            # the mate as a real record inside the parent-gene exon (read2): its sequence is a
            # slice of the exon, so the retrocopy clip and the parent mate share the transcript.
            records.append(make_read(hdr, elem_tid, nm, PG_EXON_SEQ[exon][40:40 + ref_m],
                                     mpos, f"{ref_m}M", 60,
                                     flag=0x1 | 0x80, next_tid=ti, next_start=R - a))
        # LEFT (3') junction @L: read1 forward, clip = spliced-mRNA 3' end + poly-A tail
        for k in range(alt):
            clip = sample_clip_len(rng); a = read_len - clip
            records.append(make_read(hdr, ti, f"{tag}_L{k}", ins3[-clip:] + g3[:a],
                                     L, f"{clip}S{a}M", junction_mapq(), flag=PAIRED_R1))
        coverage_reads(tag, ti, L, ref)
        truth.append((contig, L, R, "PSEUDOGENE_novel", tsd, alt, ref, vaf))
        idx += 1

    # =====================================================================
    # TPRT-hallmark insertion-type catalogue (--types / --n-per-type). Read-level,
    # multi-sample, literature-calibrated (see test/simlib). Off by default, so the legacy
    # output above is byte-identical when --types is not given. Uses its own RNG.
    # =====================================================================
    new_rows = []                                   # (truth8 tuple, label dict, ins_seq)
    sample_records = [records] + [[] for _ in range(max(1, args.samples) - 1)]
    type_keys = _simlib_models.parse_types(args.types)
    if type_keys and args.n_per_type > 0:
        trng = random.Random(args.seed * 7919 + 17)
        lib = _simlib_library.load_library(args.rte_library, args.hs1_2bit, args.hs1_rmsk,
                                           args.library_cache)
        print(lib.summary(), file=sys.stderr)
        genes = _val1_genes(args, trng)
        ctx = _simlib_models.Ctx(lib, genes, polya_scale=args.polya_scale,
                                 max_del=args.max_l1_deletion, max_dup=args.max_l1_duplication)
        sampler = FragmentSampler(trng, read_len=read_len, jitter_scale=args.polya_jitter,
                                  error_rate=args.error_rate,
                                  pcr_dup_frac=args.pcr_dup_unflagged_frac,
                                  burst_scale=args.phasing_burst)
        nsamp = max(1, args.samples)
        for key in type_keys:
            for rep in range(args.n_per_type):
                ev = _simlib_models.build_event(key, trng, ctx)
                strand = trng.choice("+-")
                rref = random.Random(trng.getrandbits(64))
                prep = prepare_event(rref, ev, strand, lib)
                contig = contigs[idx % len(contigs)]; ti = tid[contig]
                base = slot[ti]
                slot[ti] += max(step, len(prep.ref) + 20_000)
                idx += 1
                # sample presence + per-sample VAF (clonal het 0.5 or subclonal)
                if ev.role == "ARTEFACT":
                    present = {trng.randrange(nsamp)}  # library artefact: one sample only
                elif nsamp == 1 or trng.random() < args.present_all_frac:
                    present = set(range(nsamp))
                else:
                    k = trng.randint(1, nsamp)
                    present = set(trng.sample(range(nsamp), k))
                vafs, counts = [], []
                for si in range(nsamp):
                    v = 0.0
                    if si in present:
                        v = 0.5 if trng.random() < args.vaf_clonal_frac else round(trng.uniform(0.1, 0.5), 3)
                    vafs.append(v if ev.role == "TP" else (1.0 if si in present else 0.0))
                    w = PairWriter(hdr, ti, contig, base, read_len)
                    c = render_sample(trng, prep, w, sampler, args.depth, v, si in present,
                                      f"{key}_{contig}_{base}_S{si + 1}")
                    sample_records[si].extend(w.records)
                    counts.append(c)
                tr = prep.tr
                left, right = base + tr["left"], base + tr["right"]
                alt_reads = sum(c["R_reads"] + c["L_reads"] for c in counts)
                vaf_mean = round(sum(vafs) / nsamp, 3)
                row8 = (contig, left, right, f"{ev.element}_{key}", right - left, alt_reads,
                        0, vaf_mean)
                labels = _simlib_truth.event_labels(ev, tr, sorted(s + 1 for s in present), vafs,
                                                    counts)
                new_rows.append((row8, labels, tr.get("x_seq", "")))

    # =====================================================================
    # False-positive artefacts — NOT in truth. A call at one of these loci is a false
    # positive. They live in a dedicated coordinate band (>=10 Mb) clear of the real
    # insertions (<~5 Mb) and the Feature A/B bands (15-18 Mb on contig 3).
    # =====================================================================
    # start past the real insertions (default band is 10 Mb) but before the Feature
    # A/B bands (15-18 Mb, contig 3); shifts up automatically if counts push slots higher
    astate = [max(s, 10_000_000) for s in slot]
    astep = 20_000
    acount = [0]

    def next_art():
        contig = contigs[acount[0] % len(contigs)]
        ti = tid[contig]
        pos = astate[ti]; astate[ti] += astep
        acount[0] += 1
        return contig, ti, pos

    def fp_pair(prefix, clip5, clip3, a5, a3, mapq=60, flag=PAIRED_R1,
                r_tags=None, l_tags=None, n=2):
        contig, ti, A = next_art()
        tsd = rng.randint(args.tsd_min, args.tsd_max)
        R = A + tsd
        tag = f"{prefix}_{contig}_{A}"
        for k in range(n):
            records.append(make_read(hdr, ti, f"{tag}_R{k}", a5 + clip5,
                                     R - anchor_m, f"{anchor_m}M{len(clip5)}S", mapq,
                                     flag=flag, tags=r_tags))
            records.append(make_read(hdr, ti, f"{tag}_L{k}", clip3 + a3,
                                     A, f"{len(clip3)}S{anchor_m}M", mapq,
                                     flag=flag, tags=l_tags))
        coverage_reads(tag, ti, A)
        return contig, ti, A

    # ADAPTER0 is the exact adapter confirmed as a top recurrent clip in the CRAM
    # (AGATCGGAAGAGC... — adapter read-through is 5.5% of real clips, the commonest class).
    ADAPTER0 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

    # §5.5 adapter read-through (CRAM: 5.5% of clips): clip is the sequencing adapter
    # -> is_adapter drops it.
    for _ in range(args.n_adapter):
        fp_pair("ADAPT", ADAPTER0, ADAPTER0, make_anchor(), make_anchor())

    # §5.4 poly-G / dark cycle (CRAM: 2.8% of clips): clip is a poly-G run -> is_adapter drops it.
    for _ in range(args.n_polyg):
        fp_pair("POLYG", "G" * clip_s, "G" * clip_s, make_anchor(), make_anchor())

    # §5.6 homopolymer / low-complexity (CRAM: 2.6% of clips): the *anchor* is a 2-bp
    # microsatellite, so the unclipped-consensus n-polymer filter rejects the breakpoint.
    for _ in range(args.n_lowcomplexity):
        lc = ("CA" * (anchor_m // 2 + 1))[:anchor_m]
        fp_pair("LOWCX", clip_from(l1["body"], 10, clip_s),
                clip_from(l1["body"], 80, clip_s), lc, lc)

    # §5.7 PCR duplicates: every supporting read carries the duplicate flag, so
    # discovery skips them and no evidence remains.
    for _ in range(args.n_pcr_dup):
        fp_pair("PCRDUP", clip_from(l1["end5"], 0, clip_s), clip_from(l1["end3"], 0, clip_s),
                make_anchor(), make_anchor(), flag=PAIRED_R1 | 0x400)

    # §5.7 PCR chimera / template switch: reads at a junction disagree base-to-base,
    # so no consensus forms.
    for _ in range(args.n_pcr_chimera):
        contig, ti, A = next_art()
        tsd = rng.randint(args.tsd_min, args.tsd_max); R = A + tsd
        tag = f"CHIM_{contig}_{A}"
        for k in range(3):
            records.append(make_read(hdr, ti, f"{tag}_R{k}", rnd_seq(rng, anchor_m) + rnd_seq(rng, clip_s),
                                     R - anchor_m, f"{anchor_m}M{clip_s}S", 60, flag=PAIRED_R1))
            records.append(make_read(hdr, ti, f"{tag}_L{k}", rnd_seq(rng, clip_s) + rnd_seq(rng, anchor_m),
                                     A, f"{clip_s}S{anchor_m}M", 60, flag=PAIRED_R1))
        coverage_reads(tag, ti, A)

    # §5.9 segdup / mismapping: clipped reads carry an XA tag showing the whole read
    # maps contiguously elsewhere -> reject_fully_mapping_reads drops them.
    for _ in range(args.n_mismap):
        xa = [("XA", "2,+18000000,150M,0;", "Z")]
        fp_pair("MISMAP", clip_from(l1["end5"], 0, clip_s), clip_from(l1["end3"], 0, clip_s),
                make_anchor(), make_anchor(), r_tags=xa, l_tags=xa)

    # mapping-ambiguity FP — models the full-stack pericentromere / telomere / segdup /
    # assembly-discordance calls (the 66 FPs the hs1->hg38 run produced). Its ANCHOR
    # reads are mostly ambiguously placed: low MAPQ, XS within a hair of AS, and a
    # PARTIAL (not full-length) XA alt so `reject_fully_mapping_reads` does NOT fire. A
    # minority of reads look locally unique (MAPQ 60, no XS) — enough to clear the
    # >=2-per-side call floor, so it is a FALSE POSITIVE at baseline. This is the target
    # class for the anchor-uniqueness filter (#2): a true insertion has all-unique
    # anchors, this has few, so a "keep only if >=K unique-anchor reads" rule removes it.
    # uniqfrac sets the unique/ambiguous mix (0.25 ~ the full-stack FP profile).
    for _ in range(args.n_mismap_ambiguous):
        contig, ti, A = next_art()
        tsd = rng.randint(max(args.tsd_min, 3), args.tsd_max); R = A + tsd
        tag = f"MMAMB_{contig}_{A}"
        nreads = args.mismap_ambiguous_reads
        nuniq = max(2, round(args.mismap_ambiguous_uniqfrac * nreads))
        a5 = make_anchor(); a3 = make_anchor()
        clip5 = clip_from(l1["end5"], 0, clip_s); clip3 = clip_from(l1["end3"], 0, clip_s)
        part_xa = [("XA", f"2,+18000000,{anchor_m // 2}M{anchor_m // 2 + clip_s}S,2;", "Z")]
        amb_tags = part_xa + [("AS", anchor_m, "i"), ("XS", anchor_m - 2, "i")]
        for k in range(nreads):
            uniq = k < nuniq
            mq = 60 if uniq else rng.choice([0, 25, 35])
            tg = None if uniq else amb_tags
            records.append(make_read(hdr, ti, f"{tag}_R{k}", a5 + clip5,
                                     R - anchor_m, f"{anchor_m}M{clip_s}S", mq,
                                     flag=PAIRED_R1, tags=tg))
            records.append(make_read(hdr, ti, f"{tag}_L{k}", clip3 + a3,
                                     A, f"{clip_s}S{anchor_m}M", mq,
                                     flag=PAIRED_R1, tags=tg))
        coverage_reads(tag, ti, A)

    # §5.1 cruciform / inverted repeat: supplementary reads with an SA tag to the same
    # contig within 1 kb -> cluster poisoning removes them.
    for _ in range(args.n_cruciform):
        contig, ti, A = next_art()
        tsd = rng.randint(args.tsd_min, args.tsd_max); R = A + tsd
        tag = f"CRUC_{contig}_{A}"
        for k in range(2):
            sa_r = [("SA", f"{contig},{R - anchor_m + 200},+,{anchor_m}M{clip_s}S,60,0;", "Z")]
            sa_l = [("SA", f"{contig},{A + 200},+,{clip_s}S{anchor_m}M,60,0;", "Z")]
            records.append(make_read(hdr, ti, f"{tag}_R{k}", rnd_seq(rng, anchor_m) + erv5,
                                     R - anchor_m, f"{anchor_m}M{clip_s}S", 60,
                                     flag=PAIRED_R1 | 0x800, tags=sa_r))
            records.append(make_read(hdr, ti, f"{tag}_L{k}", erv5 + rnd_seq(rng, anchor_m),
                                     A, f"{clip_s}S{anchor_m}M", 60,
                                     flag=PAIRED_R1 | 0x800, tags=sa_l))
        coverage_reads(tag, ti, A)

    # §5.1 structure-specific chimera — BOTH-ENDS-CLIPPED (SMS). Unequal clips so the
    # longer side wins and a spurious pair forms. PEAR-TREE has NO SMS reject yet
    # (review R15), so this is expected to be a FALSE POSITIVE at baseline.
    for _ in range(args.n_sms):
        contig, ti, A = next_art()
        tsd = rng.randint(args.tsd_min, args.tsd_max); R = A + tsd
        tag = f"SMS_{contig}_{A}"
        big, small = clip_s, 12
        clipL = rnd_seq(rng, big); mid_l = make_anchor(); tailL = rnd_seq(rng, small)
        clipR = rnd_seq(rng, big); mid_r = make_anchor(); headR = rnd_seq(rng, small)
        for k in range(2):
            records.append(make_read(hdr, ti, f"{tag}_L{k}", clipL + mid_l + tailL,
                                     A, f"{big}S{anchor_m}M{small}S", 60, flag=PAIRED_R1))
            records.append(make_read(hdr, ti, f"{tag}_R{k}", headR + mid_r + clipR,
                                     R - anchor_m, f"{small}S{anchor_m}M{big}S", 60, flag=PAIRED_R1))
        coverage_reads(tag, ti, A)

    # §5.1 self-fold palindrome: the clip is the reverse complement of the adjacent
    # anchor. No genome-free self-palindrome check exists yet (R15) -> FALSE POSITIVE.
    for _ in range(args.n_palindrome):
        contig, ti, A = next_art()
        tsd = rng.randint(args.tsd_min, args.tsd_max); R = A + tsd
        tag = f"PAL_{contig}_{A}"
        a5 = make_anchor(); a3 = make_anchor()
        clip5 = revcomp(a5[-clip_s:]); clip3 = revcomp(a3[:clip_s])
        for k in range(2):
            records.append(make_read(hdr, ti, f"{tag}_R{k}", a5 + clip5,
                                     R - anchor_m, f"{anchor_m}M{clip_s}S", 60, flag=PAIRED_R1))
            records.append(make_read(hdr, ti, f"{tag}_L{k}", clip3 + a3,
                                     A, f"{clip_s}S{anchor_m}M", 60, flag=PAIRED_R1))
        coverage_reads(tag, ti, A)

    # §4.2.9 non-RTE true SV — CHROMOSOMAL TRANSLOCATION at a retrotransposon locus, the
    # most dangerous mimic. Two RTE loci on different contigs are joined with a few bp of
    # breakpoint microhomology, so at each locus a LEFT clip (L) and RIGHT clip (R) appear
    # with L < R — a fake TSD — and the clipped sequence is element (LTR) consensus. This
    # is geometrically and by-sequence identical to an ERV insertion (no poly-A either),
    # so BASELINE discovery CALLS it -> FALSE POSITIVE. No discovery-time filter rejects
    # it: each read's split alignment (SA) is to the partner CONTIG (so the same-contig
    # cruciform check never fires) and is partial (so maps-fully-elsewhere never fires);
    # only the combine-step end-to-end remap (the clip maps fully to the partner locus) or
    # a reciprocal-partner / cross-chromosome-SA filter can remove it. Each translocation
    # emits both reciprocal breakpoints (the partner geometry such a filter would key on),
    # with a discordant mate pointing to the partner locus.
    def transloc_locus(ti, pos, ptid, pcontig, ppos, delta, tag):
        R = pos; L = pos - delta                       # microhomology fake TSD = delta
        aR = make_anchor(); aL = make_anchor()         # genomic flank shared across the reads
        for k in range(3):
            sa_r = [("SA", f"{pcontig},{ppos},+,{clip_s}M{anchor_m}S,60,0;", "Z")]
            records.append(make_read(hdr, ti, f"{tag}_R{k}", aR + erv5,
                                     R - anchor_m, f"{anchor_m}M{clip_s}S", 60,
                                     flag=PAIRED_R1, next_tid=ptid, next_start=ppos, tags=sa_r))
            sa_l = [("SA", f"{pcontig},{ppos},+,{anchor_m}S{clip_s}M,60,0;", "Z")]
            records.append(make_read(hdr, ti, f"{tag}_L{k}", erv5 + aL,
                                     L, f"{clip_s}S{anchor_m}M", 60,
                                     flag=PAIRED_R1, next_tid=ptid, next_start=ppos, tags=sa_l))
        coverage_reads(tag, ti, pos)

    for _ in range(args.n_translocation):
        cA, tiA, A = next_art()
        tiB = (tiA + 1) % len(contigs); cB = contigs[tiB]
        B = astate[tiB]; astate[tiB] += astep          # reserve a partner slot on cB
        delta = rng.randint(max(args.tsd_min, 3), min(args.tsd_max, 20))
        transloc_locus(tiA, A, tiB, cB, B, delta, f"TRANSLOC_{cA}_{A}")
        transloc_locus(tiB, B, tiA, cA, A, delta, f"TRANSLOC_{cB}_{B}")  # reciprocal partner

    # §4.2.6 high-coverage pile-up: a dense depth spike plus a spurious clipped pair.
    # A false positive at baseline; removed by the SPEC-3 mask / SPEC-4 adaptive floor.
    for a in range(args.n_artefacts):
        contig = contigs[a % len(contigs)]
        ti = tid[contig]
        A = 5_000_000 + a * 100_000
        pileup = rng.randint(args.artefact_cov_min, args.artefact_cov_max)
        for k in range(pileup):
            st = A + (k % 400)  # starts within one 500 bp bin, aligned with the breakpoints below
            records.append(make_read(hdr, ti, f"art_{contig}_{A}_p{k}", rnd_seq(rng, ref_m),
                                     st, f"{ref_m}M", 60, flag=PAIRED_R1))
        tsd = rng.randint(args.tsd_min, args.tsd_max)
        seq_r = rnd_seq(rng, anchor_m) + rnd_seq(rng, clip_s)
        seq_l = rnd_seq(rng, clip_s) + rnd_seq(rng, anchor_m)
        for k in range(2):  # spurious RIGHT breakpoint at A+tsd
            records.append(make_read(hdr, ti, f"art_{contig}_{A}_R{k}", seq_r,
                                     (A + tsd) - anchor_m, f"{anchor_m}M{clip_s}S", args.artefact_mapq, flag=PAIRED_R1))
        for k in range(2):  # spurious LEFT breakpoint at A
            records.append(make_read(hdr, ti, f"art_{contig}_{A}_L{k}", seq_l,
                                     A, f"{clip_s}S{anchor_m}M", args.artefact_mapq, flag=PAIRED_R1))

    if max(slot + astate) > args.contig_len:
        sys.exit(f"simulated loci run past --contig-len {args.contig_len} (max {max(slot + astate)}); "
                 f"raise --contig-len or lower --step / counts")
    out_paths = sample_bam_paths(args.out_bam, len(sample_records))
    for recs, path in zip(sample_records, out_paths):
        if len(sample_records) > 1:
            hdr_s = pysam.AlignmentHeader.from_dict({
                **hdr.to_dict(), "RG": [{"ID": os.path.basename(path)[:-4], "SM": os.path.basename(path)[:-4]}]})
        else:
            hdr_s = hdr
        recs.sort(key=lambda a: (a.reference_id, a.reference_start))
        with pysam.AlignmentFile(path, "wb", header=hdr_s) as out:
            for a in recs:
                out.write(a)
        pysam.index(path)

    with open(args.out_truth, "w") as f:
        head = ["contig", "left", "right", "class", "tsd", "alt_reads", "ref_reads", "vaf"]
        if new_rows:
            head += _simlib_truth.LABEL_COLUMNS
        f.write("\t".join(head) + "\n")
        for row in truth:
            cols = [str(x) for x in row]
            if new_rows:
                lab = _simlib_truth.legacy_labels(str(row[3]))
                cols += [lab[c] for c in _simlib_truth.LABEL_COLUMNS]
            f.write("\t".join(cols) + "\n")
        for row8, lab, _ in new_rows:
            f.write("\t".join([str(x) for x in row8] + [lab[c] for c in _simlib_truth.LABEL_COLUMNS]) + "\n")
    if new_rows:
        # inserted sequences (element sense) for annotate / structure checks
        with open(args.out_truth + ".ins.fa", "w") as f:
            for row8, lab, seq in new_rows:
                if seq:
                    f.write(f">{row8[0]}:{row8[1]}-{row8[2]}|{lab['variant']}|{lab['role']}\n{seq}\n")

    # Feature A: emit a matching RepeatMasker .out track (div-gated parser format:
    # col1=%div, col4=contig, col5=begin 1-based, col6=end) covering the RTE band the
    # discordant mates map into, so the D3 config can point `discordant_rte_track` at it.
    if args.out_rmsk and rmsk_rows:
        with open(args.out_rmsk, "w") as f:
            f.write("   SW  perc perc perc  query  begin  end (left) strand repeat class ...\n")
            for contig, begin0, end0, div in rmsk_rows:
                f.write(f"1000 {div} 0.0 0.0 {contig} {begin0 + 1} {end0} (0) + ELEM LINE/L1 1 100 (0) 1\n")

    # Feature B: emit the exon annotation (contig, begin, end, gene_id; 0-based
    # half-open) so the D5 config can point `exon_annotation` at it.
    if args.out_exons and exon_rows:
        with open(args.out_exons, "w") as f:
            for contig, begin, end, gene in exon_rows:
                f.write(f"{contig}\t{begin}\t{end}\t{gene}\n")

    print(f"wrote {sum(len(r) for r in sample_records)} reads, {len(truth) + len(new_rows)} truth rows "
          f"to {', '.join(out_paths)}; truth -> {args.out_truth}", file=sys.stderr)


def sample_bam_paths(out_bam, n):
    """One sample: --out-bam as given. n>1 colonies: <stem>.S1.bam .. <stem>.Sn.bam."""
    if n <= 1:
        return [out_bam]
    stem = out_bam[:-4] if out_bam.endswith(".bam") else out_bam
    return [f"{stem}.S{i + 1}.bam" for i in range(n)]


def _val1_genes(args, rng):
    """Parent genes for the pseudogene / decoy types: a real annotation (--gene-model +
    --hs1-2bit), else gene models built on real hs1 sequence with canonical GT..AG splice
    signals (--hs1-2bit), else on random sequence."""
    G = _simlib_library
    if args.hs1_2bit and os.path.exists(args.hs1_2bit):
        from simlib.seqs import Genome
        g = Genome(args.hs1_2bit)
        if args.gene_model:
            genes = G.load_gene_model(args.gene_model, g)
            if genes:
                return genes
        genes = []
        contigs = [c for c in g.lengths if c in ("chr1", "chr2", "chr3", "chr5", "chr12", "chr17")]
        for i in range(30):
            c = rng.choice(contigs)
            for _ in range(20):
                st = rng.randint(10_000_000, g.lengths[c] - 10_000_000)
                seq = g.fetch(c, st, st + 30_000)
                if seq.count("N") == 0:
                    break
            gene = G.gene_from_sequence(f"hs1gene{i}_{c}_{st}", c, seq, st, rng,
                                        strand=rng.choice("+-"))
            if gene:
                genes.append(gene)
        return genes
    genes = []
    for i in range(10):
        gene = G.gene_from_sequence(f"synGene{i}", "synthetic", _rnd_seq(rng, 30_000), 0, rng)
        if gene:
            genes.append(gene)
    return genes


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--out-bam", required=True)
    p.add_argument("--out-truth", required=True)
    p.add_argument("--out-rmsk", help="Feature A: write a RepeatMasker .out track for the RTE band")
    p.add_argument("--out-exons", help="Feature B: write the exon annotation (contig,begin,end,gene_id)")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--element-fasta", default=None,
                   help="FASTA overriding element termini (records <ELEMENT>_end5/_end3/_body/_ltr)")
    p.add_argument("--contig-len", type=int, default=20_000_000,
                   help="length of each of the 3 synthetic contigs (raise for large --n runs)")
    p.add_argument("--step", type=int, default=100_000,
                   help="spacing between successive insertion sites (lower to pack more insertions per contig)")

    g = p.add_argument_group("element counts (real insertions, in truth)")
    g.add_argument("--n-l1", type=int, default=30, help="canonical full-length L1Hs")
    g.add_argument("--n-erv", type=int, default=20, help="canonical full-length ERV-K provirus")
    g.add_argument("--n-l1-variant", type=int, default=4, help="each non-canonical L1 variant")
    g.add_argument("--n-alu", type=int, default=4, help="each Alu variant")
    g.add_argument("--n-sva", type=int, default=3, help="each SVA variant")
    g.add_argument("--n-erv-variant", type=int, default=4, help="each non-canonical ERV variant")
    g.add_argument("--n-hervk113", type=int, default=1, help="HERV-K113 provirus (positive control)")
    g.add_argument("--n-hervk117", type=int, default=1, help="HERV-K117 provirus (positive control)")

    g = p.add_argument_group("geometry / VAF (defaults calibrated to a real GRCm38 mouse WGS CRAM)")
    g.add_argument("--read-len", type=int, default=151, help="read length (151 in the calibrating CRAM)")
    g.add_argument("--junction-lowmapq-frac", type=float, default=0.0,
                   help="fraction of junction reads dropped to MAPQ 0 (~0.12 in the CRAM; 0 keeps the gate deterministic)")
    g.add_argument("--anchor-len", type=int, default=90, help="anchor length for fixed-geometry probe reads")
    g.add_argument("--clip-len", type=int, default=60, help="clip length for fixed-geometry probe/artefact reads")
    g.add_argument("--polya-len", type=int, default=24,
                   help="max poly-A tail length in TPRT 3' clips (sampled 13..this; CRAM proximal A/T runs median 7)")
    g.add_argument("--tsd-min", type=int, default=2)
    g.add_argument("--tsd-max", type=int, default=40)
    g.add_argument("--alt-min", type=int, default=2, help="min alt (junction) reads per side => low VAF")
    g.add_argument("--alt-max", type=int, default=8)
    g.add_argument("--ref-min", type=int, default=5)
    g.add_argument("--ref-max", type=int, default=40)

    g = p.add_argument_group("sensitivity artefacts / toggles (in truth; may be MISSED at baseline)")
    g.add_argument("--n-wobble", type=int, default=0, help="SENS-1: insertions with wobbled breakpoints")
    g.add_argument("--n-lowmapq", type=int, default=0, help="SENS-2: insertions with mixed-MAPQ junction reads")
    g.add_argument("--n-shortpolya", type=int, default=0, help="SENS-8: insertions with a short poly-A clip")
    g.add_argument("--shortpolya-len", type=int, default=11, help="length of the short poly-A clip (<= 12 floor)")
    g.add_argument("--n-polya-dropout", type=int, default=0,
                   help="§5.3: L1 with its poly-A junction reads dropped (baseline miss)")
    g.add_argument("--n-materescue", type=int, default=0,
                   help="L1 into a low-mapability-but-mate-unique flank: junction clips below the MAPQ floor with a uniquely-mapped mate (MQ=60). Baseline miss; recovered by mate_anchor_rescue")
    g.add_argument("--materescue-mapq", type=int, default=20,
                   help="own MAPQ of the mate-rescue junction clips (< min_mapq 40)")

    g = p.add_argument_group("Feature A/B toggles")
    g.add_argument("--n-discordant", type=int, default=0, help="A1: one-sided junctions rescued by RTE-origin discordant mates")
    g.add_argument("--n-disc-artefact", type=int, default=0, help="A2: one-sided junctions whose discordant mates are non-RTE (not in truth)")
    g.add_argument("--discordant-min-reads", type=int, default=3, help="discordant anchor reads per one-sided junction")
    g.add_argument("--n-pseudogene", type=int, default=0, help="B3: insertions whose mates span >=2 exons of one gene")
    g.add_argument("--n-single-exon", type=int, default=0, help="B3 negative control: mates hit a single exon (must not flag)")
    g.add_argument("--n-novel-pseudogene", type=int, default=0,
                   help="novel de-novo processed pseudogene (retrocopy): a REAL TPRT insertion (in truth) whose clip is the spliced transcript + poly-A and whose mates span >=2 parent-gene exons (splice hallmark)")
    g.add_argument("--n-discordant-only", type=int, default=0,
                   help="breakpoints in the read GAP: no soft-clip anywhere, only discordant pairs (baseline miss — no discordant-pair discovery, §7.7#5)")

    g = p.add_argument_group("false-positive artefacts (NOT in truth; a call here is a FP)")
    g.add_argument("--n-adapter", type=int, default=5, help="§5.5 adapter read-through, 5.5%% of real clips (must reject)")
    g.add_argument("--n-polyg", type=int, default=3, help="§5.4 poly-G / dark cycle (must reject)")
    g.add_argument("--n-lowcomplexity", type=int, default=3, help="§5.6 low-complexity anchor (must reject)")
    g.add_argument("--n-pcr-dup", type=int, default=3, help="§5.7 PCR duplicates (must reject)")
    g.add_argument("--n-pcr-chimera", type=int, default=3, help="§5.7 PCR chimera / no consensus (must reject)")
    g.add_argument("--n-mismap", type=int, default=3, help="§5.9 segdup/mismap XA full-length (must reject)")
    g.add_argument("--n-mismap-ambiguous", type=int, default=0,
                   help="mapping-ambiguity FP mirroring the full-stack pericentromere/segdup calls: mostly ambiguous anchors (low MAPQ, XS~AS, PARTIAL XA) + a unique minority. Passes baseline; the target class for the anchor-uniqueness filter (#2)")
    g.add_argument("--mismap-ambiguous-reads", type=int, default=12,
                   help="clipped reads per side at a mapping-ambiguity FP locus")
    g.add_argument("--mismap-ambiguous-uniqfrac", type=float, default=0.85,
                   help="fraction of those reads that look locally unique (MAPQ 60, no XS). Default 0.85 matches the full-stack hs1->hg38 assembly-discordance FPs, whose shifted hg38 homolog is locally unique so most reads look confident — the hard, realistic case for the anchor-uniqueness filter")
    g.add_argument("--n-cruciform", type=int, default=3, help="§5.1 SA same-contig <1kb (must reject)")
    g.add_argument("--n-sms", type=int, default=3,
                   help="§5.1 both-ends-clipped chimera (baseline FP — no SMS filter yet, R15)")
    g.add_argument("--n-palindrome", type=int, default=3,
                   help="§5.1 self-fold palindrome (baseline FP — no self-fold filter yet, R15)")
    g.add_argument("--n-translocation", type=int, default=2,
                   help="§4.2.9 translocation between two RTE loci w/ fake TSD (baseline FP; each emits 2 reciprocal mimic loci — needs the combine-step remap / reciprocal-partner filter)")
    g.add_argument("--n-artefacts", type=int, default=0, help="§4.2.6 high-coverage pile-up region")
    g.add_argument("--artefact-cov-min", type=int, default=150)
    g.add_argument("--artefact-cov-max", type=int, default=300)
    g.add_argument("--artefact-mapq", type=int, default=60, help="MAPQ of pile-up clipped reads (30 => low-MAPQ pileup)")

    g = p.add_argument_group(
        "TPRT-hallmark insertion-type catalogue (read-level, multi-sample; off unless --types)")
    g.add_argument("--types", default="",
                   help="comma list of catalogue keys and/or literature type ids, or all / tp / artefact. "
                        "Keys: " + ", ".join(_simlib_models.CATALOGUE))
    g.add_argument("--n-per-type", type=int, default=0, help="events per selected type")
    g.add_argument("--samples", type=int, default=1,
                   help="colonies/samples of one patient; >1 writes <out-bam stem>.S1..Sn.bam "
                        "(legacy catalogue goes to S1 only)")
    g.add_argument("--depth", type=float, default=15.0, help="per-sample read depth at catalogue loci")
    g.add_argument("--present-all-frac", type=float, default=0.5,
                   help="fraction of TP events present in ALL samples (rest: random k of n)")
    g.add_argument("--vaf-clonal-frac", type=float, default=0.7,
                   help="per-sample probability an event is clonal heterozygous (VAF 0.5); else U(0.1,0.5)")
    g.add_argument("--rte-library", default=None,
                   help="resources/rte_library dir (real intact elements, transduction sources, flanks)")
    g.add_argument("--hs1-2bit", default=os.environ.get("PEARTREE_HS1_2BIT"),
                   help="hs1 .2bit for the element/flank fallback and real-sequence gene models "
                        "(env PEARTREE_HS1_2BIT)")
    g.add_argument("--hs1-rmsk", default=os.environ.get("PEARTREE_HS1_RMSK"),
                   help="hs1 RepeatMasker .out.gz for the element/flank fallback (env PEARTREE_HS1_RMSK)")
    g.add_argument("--library-cache", default=None,
                   help="cache dir for the hs1 extraction (default <hs1 dir>/simlib_rte_cache)")
    g.add_argument("--gene-model", default=None,
                   help="gene annotation (tools/build_gene_model.py TSV or GTF) for pseudogene parents")
    g.add_argument("--polya-jitter", type=float, default=1.0,
                   help="scale of per-read homopolymer (SBS slippage) length jitter; 0 disables")
    g.add_argument("--phasing-burst", type=float, default=1.0,
                   help="scale of post-homopolymer phasing loss (read turns to low-quality junk after "
                        "a long poly-A, p=min(0.6,0.01*(n-12)) per read); 0 disables")
    g.add_argument("--polya-scale", type=float, default=1.0, help="scale the median poly-A length")
    g.add_argument("--pcr-dup-unflagged-frac", type=float, default=0.05,
                   help="fraction of fragments with PCR/optical duplicate copies (0x400 NOT set)")
    g.add_argument("--error-rate", type=float, default=0.002, help="per-base substitution rate")
    g.add_argument("--max-l1-deletion", type=int, default=20000, help="max L1-mediated deletion (bp)")
    g.add_argument("--max-l1-duplication", type=int, default=5000, help="max L1-mediated duplication (bp)")

    simulate(p.parse_args())


if __name__ == "__main__":
    main()
