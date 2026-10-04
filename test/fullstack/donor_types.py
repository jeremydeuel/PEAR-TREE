#!/usr/bin/env python3
"""Multi-sample, insertion-type-catalogue donor builder (build_donor.py --types ...).

Builds, for ONE patient with n colonies/samples, per-sample haplotype FASTAs of one or more
hs1 windows carrying events from the TPRT-hallmark catalogue (test/simlib/models.py), plus
library-artefact molecules, a junction table for read-support counting, and a truth table in
hs1 coordinates. `simulate_reads.py` turns the haplotypes into reads (poly-A jitter, unflagged
PCR duplicates), `run_multisample.sh` aligns them and runs discovery, `score_types.py` scores.

VAF per sample is realised by haplotype mixing: the reference haplotype has weight 0.5 and
four alt haplotypes weight 0.125 each; an event with VAF k/8 (k = 1..4) is placed in the first
k alt haplotypes (clonal heterozygous = 0.5 = all four). Sample presence: an event is in all
samples with probability --present-all-frac, else in a random k-of-n subset; library
artefacts exist in exactly one sample (they are per-library).

Sites: events that carry an EN motif are placed at REAL hs1 positions whose 6-mer matches
TT|AAAA with the event's sampled mismatch count (0-3; the strand follows the motif); the
genome is never edited except by the insertion itself.

Outputs (in --out-dir):
  ref.hap.fa                       unmodified windows (weight 0.5 in every sample)
  S<i>.hap<k>.fa (k=1..4)          alt haplotypes of sample i
  S<i>.molecules.fa                single-molecule artefacts (+PCR copies) for sample i
  haps.tsv                         sample, hap file, weight
  junctions.tsv                    hap file, event id, side, hap position
  slippage.tsv                     hap file, start, end of artefact A-tracts (boosted slippage)
  truth_types_hs1.tsv              truth (TRUTH_SCHEMA.md), hs1 coordinates
  truth_types.ins.fa               inserted sequences (element sense)
  flanks.fa                        300 bp hs1 flanks per event (bwa-map to lift the truth)
  sources_hs1.fa                   hs1 source loci (element + 15 kb flank) used by transductions
"""
import argparse
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from simlib import library as L  # noqa: E402
from simlib import models as M  # noqa: E402
from simlib import phylo as PH  # noqa: E402
from simlib import truth as T  # noqa: E402
from simlib.seqs import Genome, homopolymer_runs, revcomp, write_fasta  # noqa: E402
from simlib.val1render import window_for  # noqa: E402

N_ALT = 4                 # alt haplotypes per sample, weight 0.125 each
ALT_W = 0.125
UNSUPPORTED = {"ART_SUBFAMILY_MISMAP"}   # organic in fullstack (reads from inserted L1s mismap)


def parse_region(s):
    c, rest = s.split(":")
    a, b = rest.replace(",", "").split("-")
    return c, int(a), int(b)


def quantize_vaf(v):
    k = max(1, min(N_ALT, int(round(v / ALT_W))))
    return k, k * ALT_W


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--types", required=True)
    p.add_argument("--n-per-type", type=int, default=3)
    p.add_argument("--samples", type=int, default=2)
    p.add_argument("--hs1", default="/Users/jeremy/Downloads/hs1.fa", help="hs1 FASTA (indexed) or .2bit")
    p.add_argument("--hs1-rmsk", default=None, help="hs1 RepeatMasker .out.gz (element/flank fallback)")
    p.add_argument("--rte-library", default=None, help="resources/rte_library dir")
    p.add_argument("--library-cache", default=None)
    p.add_argument("--gene-model", default=None, help="hs1 gene annotation (TSV/GTF) for pseudogene parents")
    p.add_argument("--region", action="append", default=None,
                   help="hs1 window(s) to implant into (repeatable), default chr22:26000000-30000000")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--min-site-gap", type=int, default=15000)
    p.add_argument("--present-all-frac", type=float, default=0.5)
    p.add_argument("--vaf-clonal-frac", type=float, default=0.7)
    p.add_argument("--polya-scale", type=float, default=1.0)
    p.add_argument("--max-l1-deletion", type=int, default=8000)
    p.add_argument("--max-l1-duplication", type=int, default=3000)
    # tree-structured presence (test/simlib/phylo.py); off = the legacy random k-of-n subsets
    p.add_argument("--tree", default=None,
                   help="Newick file (tips -> S1..Sn in leaf order) or random:N. Overrides --samples. "
                        "TP events go on one branch (carried by exactly its clade, clonal het x colony "
                        "purity), a fraction become non-clade decoys; writes tree.nwk, samples.tsv, "
                        "phylo_truth.tsv")
    p.add_argument("--tree-root-frac", type=float, default=0.1, help="TP events carried by every colony")
    p.add_argument("--tree-branch-weight", choices=["uniform", "length"], default="uniform")
    p.add_argument("--tree-nonclade-frac", type=float, default=0.15,
                   help="TP events placed on a random NON-clade subset (>= 2 carriers, >= 1 non-carrier)")
    p.add_argument("--tree-germline-frac", type=float, default=0.5,
                   help="of the ROOT events, the fraction that are germline het (VAF 0.5 in every "
                        "colony, contaminating cells included) rather than somatic before the MRCA")
    p.add_argument("--purity-range", default="0.7,1.0", help="tree mode: per-colony purity U(lo,hi)")
    p.add_argument("--low-depth-frac", type=float, default=0.25, help="tree mode: low-depth colonies")
    p.add_argument("--low-depth-range", default="0.25,0.5", help="tree mode: their depth factor U(lo,hi)")
    a = p.parse_args(argv)

    rng = random.Random(a.seed)
    os.makedirs(a.out_dir, exist_ok=True)
    g = Genome(a.hs1)
    lib = L.load_library(a.rte_library, a.hs1, a.hs1_rmsk, a.library_cache)
    print(lib.summary(), file=sys.stderr)
    regions = [parse_region(r) for r in (a.region or ["chr22:26000000-30000000"])]
    nsamp = max(1, a.samples)
    design, trng, placements = None, None, {}
    if a.tree:
        trng = random.Random(a.seed * 7919 + 17)          # separate stream: event building unchanged
        design = PH.build_design(a.tree, trng, tuple(map(float, a.purity_range.split(","))),
                                 a.low_depth_frac, tuple(map(float, a.low_depth_range.split(","))))
        if a.samples not in (design.n, 2):
            print(f"  --tree has {design.n} tips; overriding --samples {a.samples}", file=sys.stderr)
        nsamp = design.n
    keys = [k for k in M.parse_types(a.types) if k not in UNSUPPORTED]

    # ---- per-region sequence, EN sites, A-tracts, gene models ---------------------------
    reg = []
    for (c, rs, re_) in regions:
        seq = g.fetch(c, rs, re_)
        sites = M.find_en_sites(seq, max_mm=3, margin=5000)
        by_mm = {}
        for nick, st, mm in sites:
            by_mm.setdefault(mm, []).append((nick, st))
        tracts = [(r0, r1) for (r0, r1, b) in homopolymer_runs(seq, 18) if b in "AT"
                  and 5000 < r0 < len(seq) - 5000]
        reg.append({"contig": c, "start": rs, "seq": seq, "by_mm": by_mm, "tracts": tracts,
                    "used": []})
    genes = []
    if a.gene_model:
        genes = [x for x in L.load_gene_model(a.gene_model, g)
                 if any(x.contig == r["contig"] and r["start"] <= x.exons[0][0] and
                        x.exons[-1][1] <= r["start"] + len(r["seq"]) for r in reg)]
    if not genes:   # parent genes on REAL window sequence with canonical GT..AG introns
        for i in range(12):
            r = rng.choice(reg)
            off = rng.randint(10_000, len(r["seq"]) - 40_000)
            gene = L.gene_from_sequence(f"hs1gene{i}_{r['contig']}_{r['start'] + off}", r["contig"],
                                        r["seq"][off:off + 30_000], r["start"] + off, rng,
                                        strand=rng.choice("+-"))
            if gene and "N" not in gene.mrna():
                genes.append(gene)
    ctx = M.Ctx(lib, genes, a.polya_scale, a.max_l1_deletion, a.max_l1_duplication)

    def free(r, lo, hi):
        if lo < 5000 or hi > len(r["seq"]) - 5000 or "N" in r["seq"][lo:hi]:
            return False
        return all(hi + a.min_site_gap <= u0 or lo >= u1 + a.min_site_gap for u0, u1 in r["used"])

    # ---- events -------------------------------------------------------------------------
    events = []      # dicts
    eid = 0
    for key in keys:
        for _ in range(a.n_per_type):
            ev = M.build_event(key, rng, ctx)
            placed = None
            for attempt in range(400):
                r = rng.choice(reg)
                if ev.render == "slippage":
                    if not r["tracts"]:
                        continue
                    t0, t1 = rng.choice(r["tracts"])
                    lo, hi = t0 - 400, t1 + 400
                    if free(r, lo, hi):
                        placed = (r, None, None, lo, hi, (t0, t1))
                        break
                    continue
                Nw, nick_off = window_for(ev, W=400)
                if ev.plant_en:
                    pool = r["by_mm"].get(ev.en_mm) or r["by_mm"].get(1) or []
                    if not pool:
                        continue
                    nick, strand = rng.choice(pool)
                else:
                    nick, strand = rng.randint(6000, len(r["seq"]) - 6000), rng.choice("+-")
                lo = nick - nick_off
                hi = lo + Nw
                if free(r, lo, hi):
                    placed = (r, nick, strand, lo, hi, None)
                    break
            if placed is None:
                print(f"  could not place {key}; skipped", file=sys.stderr)
                continue
            r, nick, strand, lo, hi, tract = placed
            r["used"].append((lo, hi))
            e = {"id": eid, "key": key, "ev": ev, "reg": r, "lo": lo, "hi": hi}
            eid += 1
            if tract:
                e.update(strand="+", alt=None, tract=tract,
                         left=tract[1], right=tract[0])
            elif ev.render == "foldback_reads":
                e.update(strand=strand, alt=None, left=nick, right=nick)
            else:
                win = r["seq"][lo:hi]
                alt = M.apply_event(win, nick - lo, strand, ev)
                e.update(strand=strand, alt=alt, left=lo + alt.left, right=lo + alt.right)
            # presence / VAF
            if design is not None:
                if ev.role == "ARTEFACT":
                    s = PH.draw_nonclade(trng, design) if ev.render == "slippage" else None
                    if s is not None:      # systematic artefact: same tract slips in many colonies
                        present, pl = set(s), "NONCLADE_ARTEFACT"
                    else:                  # per-library artefact: one colony
                        present, pl = {trng.randrange(nsamp)}, "ARTEFACT"
                else:
                    s = PH.draw_nonclade(trng, design) if trng.random() < a.tree_nonclade_frac else None
                    if s is not None:
                        present, pl = set(s), "NONCLADE"
                    else:
                        pl, s = PH.draw_branch(trng, design, a.tree_root_frac, a.tree_branch_weight)
                        present = set(s)
                        if pl == "ROOT" and trng.random() < a.tree_germline_frac:
                            pl, e["germline"] = "GERMLINE", True
                placements[e["id"]] = (pl, frozenset(present))
            elif ev.role == "ARTEFACT":
                present = {rng.randrange(nsamp)}
            elif nsamp == 1 or rng.random() < a.present_all_frac:
                present = set(range(nsamp))
            else:
                present = set(rng.sample(range(nsamp), rng.randint(1, nsamp)))
            kv = []
            for si in range(nsamp):
                if si not in present:
                    kv.append((0, 0.0))
                elif ev.role == "ARTEFACT":
                    kv.append((0, 1.0))
                elif design is not None and e.get("germline"):   # in every cell: VAF 0.5
                    kv.append((N_ALT, 0.5))
                elif design is not None:   # clonal het in every carrier; purity dilutes it
                    kv.append((N_ALT, 0.5 * design.purity[si]))
                else:
                    v = 0.5 if rng.random() < a.vaf_clonal_frac else rng.uniform(0.1, 0.5)
                    kv.append(quantize_vaf(v))
            e["present"], e["kv"] = present, kv
            events.append(e)

    # ---- haplotypes ---------------------------------------------------------------------
    haps_rows, junc_rows, slip_rows = [], [], []
    ref_recs = {f"{r['contig']}_{r['start']}": r["seq"] for r in reg}
    write_fasta(os.path.join(a.out_dir, "ref.hap.fa"), ref_recs, width=0)
    has_germ = any(e.get("germline") for e in events)
    for si in range(nsamp):
        pur = design.purity[si] if design is not None else 1.0
        # tree mode: founder-derived cells (fraction pur) carry the colony's somatic events on one
        # haplotype (4 alt haps, 0.125 pur each); germline-het events are on one haplotype of
        # EVERY cell, so the contaminating cells' copy is an extra "germ" hap (0.5 (1 - pur))
        w_germ = 0.5 * (1.0 - pur) if has_germ else 0.0
        haps_rows.append((si + 1, "ref.hap.fa", 1.0 - 0.5 * pur - w_germ))
        hap_specs = [(f"S{si + 1}.hap{k}.fa", ALT_W * pur, (lambda e, k=k: e["kv"][si][0] >= k))
                     for k in range(1, N_ALT + 1)]
        if w_germ > 0:
            hap_specs.append((f"S{si + 1}.germ.fa", w_germ, lambda e: bool(e.get("germline"))))
        for fname, hw, pick in hap_specs:
            recs = {}
            for r in reg:
                rname = f"{r['contig']}_{r['start']}"
                evs = sorted((e for e in events if e["reg"] is r and e["alt"] is not None
                              and e["ev"].render == "normal" and pick(e)),
                             key=lambda e: e["lo"])
                chunks, prev, off = [], 0, 0
                for e in evs:
                    chunks.append(r["seq"][prev:e["lo"]])
                    off += e["lo"] - prev
                    x0, x1 = e["alt"].x_span
                    junc_rows.append((fname, e["id"], "R", rname, off + x0))
                    junc_rows.append((fname, e["id"], "L", rname, off + x1))
                    chunks.append(e["alt"].alt)
                    off += len(e["alt"].alt)
                    prev = e["hi"]
                chunks.append(r["seq"][prev:])
                recs[rname] = "".join(chunks)
                # artefact A-tracts present in this sample: boosted slippage (positions shift
                # by the net length of the insertions placed before them)
                for e in events:
                    if e["reg"] is r and e["ev"].render == "slippage" and si in e["present"]:
                        t0, t1 = e["tract"]
                        shift = sum(len(x["alt"].alt) - (x["hi"] - x["lo"]) for x in evs if x["hi"] <= t0)
                        slip_rows.append((fname, rname, t0 + shift, t1 + shift))
            write_fasta(os.path.join(a.out_dir, fname), recs)
            haps_rows.append((si + 1, fname, hw))
        for e in events:            # slippage also shows on the reference haplotype
            if e["ev"].render == "slippage" and si in e["present"]:
                rname = f"{e['reg']['contig']}_{e['reg']['start']}"
                slip_rows.append((f"S{si + 1}:ref.hap.fa", rname, e["tract"][0], e["tract"][1]))

    # ---- single-molecule artefacts ------------------------------------------------------
    for si in range(nsamp):
        mols = {}
        for e in events:
            ev = e["ev"]
            if si not in e["present"] or ev.render not in ("single_fragment", "foldback_reads"):
                continue
            r = e["reg"]
            copies = ev.info.get("pcr_copies", 0)
            if ev.render == "foldback_reads":
                j = e["left"]
                for m in range(ev.info["n_molecules"]):
                    a1 = rng.randint(120, 300); b1 = rng.randint(20, 200)
                    mol = r["seq"][j - a1:j] + revcomp(r["seq"][j - b1:j])
                    mols[f"e{e['id']}m{m}|copies={copies}"] = mol
                continue
            alt = e["alt"]
            x0, x1 = alt.x_span
            sides = ev.info.get("sides", "POLYA")
            bs = ([x1 if alt.polya_side == "LEFT" else x0] if sides == "POLYA" else [x0, x1])
            for m, b in enumerate(bs):
                fl = max(180, int(rng.gauss(330, 90)))
                st = max(0, b - rng.randint(30, fl - 30))
                mols[f"e{e['id']}m{m}|copies={copies}"] = alt.alt[st:st + fl]
        write_fasta(os.path.join(a.out_dir, f"S{si + 1}.molecules.fa"), mols)

    with open(os.path.join(a.out_dir, "haps.tsv"), "w") as f:
        f.write("sample\thap\tweight\n")
        for row in haps_rows:
            f.write("\t".join(map(str, row)) + "\n")
    with open(os.path.join(a.out_dir, "junctions.tsv"), "w") as f:
        f.write("hap\tevent_id\tside\trecord\tpos\n")
        for row in junc_rows:
            f.write("\t".join(map(str, row)) + "\n")
    with open(os.path.join(a.out_dir, "slippage.tsv"), "w") as f:
        f.write("hap\trecord\tstart\tend\n")
        for row in slip_rows:
            f.write("\t".join(map(str, row)) + "\n")

    # ---- truth --------------------------------------------------------------------------
    cols = ["id", "hs1_contig", "hs1_left", "hs1_right", "tsd", "donor_record"] + T.LABEL_COLUMNS
    flanks, ins_fa, src_used = {}, {}, set()
    with open(os.path.join(a.out_dir, "truth_types_hs1.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for e in events:
            ev, r = e["ev"], e["reg"]
            alt = e["alt"]
            tr = {"strand": e["strand"]}
            if alt is not None:
                tr.update(ins_len=len(alt.x_seq), polya_side=alt.polya_side, tsd_seq=alt.tsd_seq,
                          en_motif=alt.en_motif, en_mm=alt.en_mm, mh_seq=alt.mh_seq, parts=alt.parts)
                ins_fa[f"{e['id']}|{e['key']}|{ev.role}"] = alt.x_seq
            vafs = [v for _, v in e["kv"]]
            counts = [{"R_frags": ".", "L_frags": ".", "R_reads": ".", "L_reads": "."}] * nsamp
            lab = T.event_labels(ev, tr, sorted(s + 1 for s in e["present"]), vafs, counts)
            hl, hr = r["start"] + e["left"], r["start"] + e["right"]
            row = [e["id"], r["contig"], hl, hr, e["right"] - e["left"], f"{r['contig']}_{r['start']}"]
            f.write("\t".join(map(str, row + [lab[c] for c in T.LABEL_COLUMNS])) + "\n")
            lo, hi = sorted((e["left"], e["right"]))
            flanks[f"{e['id']}_L"] = r["seq"][lo - 300:lo]
            flanks[f"{e['id']}_R"] = r["seq"][hi:hi + 300]
            if ev.source_id != ".":
                src_used.add(ev.source_id)
    write_fasta(os.path.join(a.out_dir, "flanks.fa"), flanks)
    write_fasta(os.path.join(a.out_dir, "truth_types.ins.fa"), ins_fa)
    src_recs = {}
    for klass in ("L1", "SVA"):
        for s in lib.sources[klass]:
            if s.id in src_used:
                src_recs[f"hs1src_{s.id}"] = s.flank5 + s.element_seq + s.tail + s.flank3
    write_fasta(os.path.join(a.out_dir, "sources_hs1.fa"), src_recs)
    # the pseudogene-parent gene models actually used (synthetic or from --gene-model), in the
    # tools/build_gene_model.py track format (contig start end gene strand, 0-based half-open,
    # hs1): annotate gets the SAME exons as truth (test/e2e/make_gene_track.py)
    with open(os.path.join(a.out_dir, "genes_hs1.tsv"), "w") as f:
        for gm in genes:
            for s, e in gm.exons:
                f.write(f"{gm.contig}\t{s}\t{e}\t{gm.id}\t{gm.strand}\n")
    if design is not None:
        PH.write_design(design, a.out_dir, placements)
    n_tp =sum(1 for e in events if e["ev"].role == "TP")
    print(f"placed {len(events)} events ({n_tp} TP, {len(events) - n_tp} artefact) in "
          f"{len(reg)} window(s) for {nsamp} sample(s) -> {a.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
