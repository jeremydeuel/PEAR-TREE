#!/usr/bin/env python3
"""Append an accepted (novel) 3' transduction source to resources/rte_library/.

Implements the bookkeeping half of the "Novel source rule" (plans/tprt_hallmarks/SPEC.md,
docs/transduction_sources.html): once annotate has reported `TD3P_SOURCE=novel:<hs1 coords>`
and a curator accepts it, this script
  1. checks the candidate against the acceptance criteria: reference L1 >= 5.5 kb with identity
     to the L1HS consensus >= 0.98 (tier A), or 0.95-0.98 with >= 2 independent daughters
     (tier B); or a cohort-called non-reference L1 insertion point,
  2. refuses duplicates (an existing source on the same strand within 1 kb),
  3. extracts the strand-aware downstream flank from hs1 (soft-masked if --hs1-rmsk is given),
  4. appends one row to transduction_sources.tsv and one record to flanks_3p.fa.gz
     (re-bgzipped + re-indexed) and refreshes manifest.tsv.

Reference element present on hs1 (1-based inclusive element coordinates):
  python tools/rte_library/add_source.py --library resources/rte_library \
      --hs1-2bit /path/hs1.2bit --chrom chr5 --start 1000001 --end 1006030 --strand + \
      --evidence "PEAR-TREE:PD12345 (3 daughters)" --n-daughters 3

Cohort-called non-reference L1 (0-based junction offset on hs1; flank starts there):
  python tools/rte_library/add_source.py --library resources/rte_library \
      --hs1-2bit /path/hs1.2bit --chrom chr5 --junction 1000000 --strand - \
      --evidence "PEAR-TREE cohort call PD12345:chr5:1000000" --n-daughters 2
"""
import argparse
import gzip
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from common import TwoBit, revcomp, read_fasta, cons_identity, write_fasta  # noqa: E402

MIN_LEN = 5500          # >= 5.5 kb: 5'UTR promoter present (transcription needs it)
ID_TIER_A = 0.98        # L1HS-like: 95 % of full-length L1HS >= 0.989, 95 % of published sources >= 0.991
ID_TIER_B = 0.95        # L1PA2 / young L1PA3 range (L1PA2 median 0.974); needs >= 2 daughters
FLANK = 15000           # 99.8 % of published transduction distal ends lie within 15 kb


def load_sources(path):
    with open(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(cols, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]
    return cols, rows


def source_tier(length, ident, n_daughters):
    """Novel-source acceptance (docs/transduction_sources.html). Returns (tier, failures):
    tier 'A' = credible (L1HS-like), 'B' = reasonably similar (accepted with >= 2 independent
    daughters), 'reject'. Identity = common.cons_identity to the L1HS consensus."""
    fails = []
    if length < MIN_LEN:
        fails.append("length %d < %d" % (length, MIN_LEN))
    if ident >= ID_TIER_A:
        tier = "A"
    elif ident >= ID_TIER_B:
        tier = "B"
        if n_daughters < 2:
            fails.append("identity %.4f in tier B (%.2f-%.2f) needs >= 2 independent daughters"
                         % (ident, ID_TIER_B, ID_TIER_A))
    else:
        tier = "reject"
        fails.append("identity to L1HS consensus %.4f < %.2f" % (ident, ID_TIER_B))
    return tier, fails


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--library", required=True)
    ap.add_argument("--hs1-2bit", required=True)
    ap.add_argument("--hs1-rmsk", help="hs1 RepeatMasker .out(.gz) for soft-masking (optional)")
    ap.add_argument("--work", default=None, help="cache dir for the parsed rmsk (default: library/../.cache)")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, help="element start, 1-based (reference element)")
    ap.add_argument("--end", type=int, help="element end, 1-based inclusive (reference element)")
    ap.add_argument("--junction", type=int, help="0-based 3' junction offset (non-reference insertion)")
    ap.add_argument("--strand", required=True, choices=["+", "-"])
    ap.add_argument("--element-class", default="L1", choices=["L1", "SVA"])
    ap.add_argument("--evidence", required=True)
    ap.add_argument("--n-daughters", type=int, default=1)
    ap.add_argument("--flank", type=int, default=FLANK)
    ap.add_argument("--id", default=None)
    ap.add_argument("--force", action="store_true", help="append even if criteria fail")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    lib = a.library
    tsv = os.path.join(lib, "transduction_sources.tsv")
    fa = os.path.join(lib, "flanks_3p.fa.gz")
    cols, rows = load_sources(tsv)
    hs1 = TwoBit(a.hs1_2bit)
    cons = {n: s for n, _, s in read_fasta(os.path.join(lib, "consensus.fa"))}

    if a.junction is None and (a.start is None or a.end is None):
        raise SystemExit("give either --start/--end (reference element) or --junction (non-reference)")
    point = a.junction is not None
    ident = "."
    fails = []
    tier = "cohort_nonref" if point else "."
    if not point:
        g = hs1.seq(a.chrom, a.start - 1, a.end)
        el = g if a.strand == "+" else revcomp(g)
        idt = cons_identity(el, cons["L1HS"]) if a.element_class == "L1" else 1.0
        ident = "%.4f" % idt
        if a.element_class == "L1":
            tier, fails = source_tier(len(el), idt, a.n_daughters)
        anchor = a.end if a.strand == "+" else a.start
    else:
        anchor = a.junction
    for f in fails:
        print("criterion failed: %s" % f, file=sys.stderr)
    if fails and not a.force:
        raise SystemExit("not accepted (use --force to override and say why in --evidence)")

    # duplicate check: same contig + strand, 3' end within 1 kb
    for r in rows:
        if r["hs1_chrom"] != a.chrom or r["hs1_strand"] != a.strand:
            continue
        try:
            rs, re_ = int(r["hs1_start"]), int(r["hs1_end"])
        except ValueError:
            continue
        r_anchor = re_ if a.strand == "+" else rs
        if abs(r_anchor - anchor) <= 1000:
            raise SystemExit("already in the library as %s" % r["id"])

    sl = "f" if a.strand == "+" else "r"
    sid = a.id or "%s_%s_%d_%s" % (a.element_class, a.chrom, (a.junction if point else a.start), sl)
    if any(r["id"] == sid for r in rows):
        raise SystemExit("id %s exists" % sid)

    from sources import flank_window, pas_hexamers
    if point:
        f0, f1 = flank_window(a.junction, a.junction, a.strand, a.flank, "3p", point=True)
    else:
        f0, f1 = flank_window(a.start, a.end, a.strand, a.flank, "3p")
    g = hs1.seq(a.chrom, f0, f1)
    if a.hs1_rmsk:
        import rmsk as rm
        work = a.work or os.path.join(os.path.dirname(os.path.abspath(lib)), ".cache")
        os.makedirs(work, exist_ok=True)
        _, mask = rm.load(a.hs1_rmsk, work, "hs1")
        g = rm.softmask(g, a.chrom, f0, mask)
    seq = g if a.strand == "+" else revcomp(g)
    desc = "hs1:%s:%d-%d(%s) 3p_flank len=%d" % (a.chrom, f0 + 1, f1, a.strand, len(seq))

    row = {c: "." for c in cols}
    row.update(dict(
        id=sid, element_class=a.element_class, subfamily=".", ta_status=".",
        reference="no" if point else "yes",
        hs1_status="insertion_point" if point else "present", hs1_chrom=a.chrom,
        hs1_start=str(a.junction if point else a.start), hs1_end=str(a.junction if point else a.end),
        hs1_strand=a.strand, hg38_start="-1", hg38_end="-1", strand=a.strand,
        strand_source="curator", identity_L1HS=ident, evidence=a.evidence, seed="novel_accepted",
        n_daughters=str(a.n_daughters), daughters_by_study="cohort:%d" % a.n_daughters,
        hotness="strong" if a.n_daughters >= 5 else "active", flank_3p=sid,
        flank_3p_len=str(a.flank), flank_3p_genome="hs1", pas_hexamers_3p=pas_hexamers(seq),
        flank_5p=".", notes="NOVEL_SOURCE;tier=%s" % tier + (";forced: " + "; ".join(fails) if fails else "")))
    if not point:
        from build import ta_status
        t = ta_status(el)
        row["ta_status"] = t
    print("\t".join(row[c] for c in cols))
    if a.dry_run:
        return

    with open(tsv, "a") as fh:
        fh.write("\t".join(row[c] for c in cols) + "\n")
    # append to the bgzipped FASTA: decompress, append, recompress, re-index
    tmp = fa[:-3] + ".tmp"
    with gzip.open(fa, "rt") as src, open(tmp, "w") as dst:
        shutil.copyfileobj(src, dst)
        write_fasta(dst, sid, seq, desc)
    import pysam
    for ext in ("", ".fai", ".gzi"):
        if os.path.exists(fa + ext):
            os.remove(fa + ext)
    pysam.tabix_compress(tmp, fa, force=True)
    os.remove(tmp)
    pysam.faidx(fa)
    from build import manifest
    manifest(lib)
    print("appended %s (%d bp flank)" % (sid, len(seq)), file=sys.stderr)


if __name__ == "__main__":
    main()
