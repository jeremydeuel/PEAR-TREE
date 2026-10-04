#!/usr/bin/env python3
"""Render docs/transduction_sources.html from resources/rte_library/ (self-contained HTML:
inline CSS + SVG + a few lines of inline JS for the table filter; no external resources).

  python tools/rte_library/make_doc.py --library resources/rte_library --out docs/transduction_sources.html
"""
import argparse
import csv
import html
import os


def tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def esc(x):
    return html.escape(str(x))


CSS = """
:root{--bg:#fbfaf7;--fg:#1d1d1f;--muted:#5f6368;--card:#ffffff;--line:#e3e1dc;--accent:#2f6f4e;
--accent2:#9a3b2f;--hot:#b3261e;--strong:#c26a00;--active:#2f6f4e;--none:#6b6b6b;--code:#f2f0eb;
--l1:#2f6f4e;--flank:#c26a00;--polya:#4a5fc1;--tsd:#7a4fa3}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#141414;--fg:#ececec;
--muted:#a0a0a0;--card:#1e1e1e;--line:#333;--accent:#7cc79f;--accent2:#e38b7d;--hot:#ff8a80;
--strong:#ffb74d;--active:#7cc79f;--none:#9e9e9e;--code:#262626;--l1:#7cc79f;--flank:#ffb74d;
--polya:#9fa8ff;--tsd:#c9a3ef}}
:root[data-theme="dark"]{--bg:#141414;--fg:#ececec;--muted:#a0a0a0;--card:#1e1e1e;--line:#333;
--accent:#7cc79f;--accent2:#e38b7d;--hot:#ff8a80;--strong:#ffb74d;--active:#7cc79f;--none:#9e9e9e;
--code:#262626;--l1:#7cc79f;--flank:#ffb74d;--polya:#9fa8ff;--tsd:#c9a3ef}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:16px/1.55 -apple-system,BlinkMacSystemFont,
"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
main{max-width:980px;margin:0 auto;padding:24px 16px 64px}
h1{font-size:1.7rem;margin:.2em 0 .1em}h2{font-size:1.25rem;margin-top:2.2em;border-bottom:1px solid var(--line);padding-bottom:.25em}
h3{font-size:1.05rem;margin-top:1.6em}
p,li{max-width:72ch}.muted{color:var(--muted)}
code{background:var(--code);padding:.1em .35em;border-radius:4px;font-size:.9em;word-break:break-word}
pre{background:var(--code);padding:12px;border-radius:8px;overflow-x:auto;font-size:.85em;line-height:1.45}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px;margin:14px 0}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:10px}
.stat{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:10px 12px}
.stat b{display:block;font-size:1.4rem}
.tablewrap{overflow-x:auto;border:1px solid var(--line);border-radius:10px;background:var(--card)}
table{border-collapse:collapse;width:100%;font-size:.86rem}
th,td{padding:6px 8px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top;white-space:nowrap}
th{position:sticky;top:0;background:var(--card)}
td.wrap{white-space:normal;min-width:180px}
.pill{display:inline-block;padding:0 .5em;border-radius:999px;font-size:.8em;border:1px solid currentColor}
.hot{color:var(--hot)}.strong{color:var(--strong)}.active{color:var(--active)}.none_reported,.candidate{color:var(--none)}
input[type=search]{width:100%;max-width:420px;padding:8px 10px;border:1px solid var(--line);border-radius:8px;
background:var(--card);color:var(--fg);font-size:1rem;margin:8px 0}
svg{max-width:100%;height:auto}
.fig{overflow-x:auto}.fig svg{min-width:600px}
svg text{fill:var(--fg);font-size:12px}
"""

FIG = """
<svg viewBox="0 0 760 250" role="img" aria-label="Source L1, downstream flank window and a daughter insertion carrying a 3' transduction">
  <text x="10" y="18" font-weight="600">Source locus (reference or polymorphic full-length L1, + strand)</text>
  <line x1="10" y1="50" x2="750" y2="50" stroke="var(--muted)" stroke-width="2"/>
  <rect x="60" y="38" width="300" height="24" rx="3" fill="var(--l1)" opacity=".85"/>
  <text x="150" y="55" style="fill:#fff">L1 5'UTR · ORF1 · ORF2 · 3'UTR</text>
  <text x="330" y="80" class="muted" style="font-size:11px">weak pA</text>
  <line x1="360" y1="30" x2="360" y2="70" stroke="var(--fg)" stroke-dasharray="3 3"/>
  <rect x="360" y="40" width="380" height="20" rx="3" fill="var(--flank)" opacity=".25"/>
  <rect x="360" y="40" width="120" height="20" rx="3" fill="var(--flank)" opacity=".9"/>
  <text x="365" y="98">flanks_3p window: 0 – 15 kb downstream, element sense</text>
  <text x="476" y="34" style="font-size:11px">AATAAA</text>
  <line x1="490" y1="38" x2="490" y2="62" stroke="var(--fg)"/>
  <text x="365" y="115" class="muted" style="font-size:11px">readthrough transcript ends at a downstream pA site (endpoint clusters per source)</text>
  <path d="M 380 125 C 380 150, 300 150, 300 168" stroke="var(--muted)" fill="none" marker-end="url(#a)"/>
  <defs><marker id="a" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="var(--muted)"/></marker></defs>
  <text x="10" y="165" font-weight="600">Daughter insertion (TPRT)</text>
  <line x1="10" y1="200" x2="750" y2="200" stroke="var(--muted)" stroke-width="2"/>
  <rect x="120" y="190" width="24" height="20" fill="var(--tsd)" opacity=".8"/>
  <rect x="144" y="190" width="180" height="20" fill="var(--l1)" opacity=".85"/>
  <text x="150" y="205" style="fill:#fff;font-size:11px">5'-truncated L1</text>
  <rect x="324" y="190" width="120" height="20" fill="var(--flank)" opacity=".9"/>
  <text x="330" y="205" style="font-size:11px">transduced flank</text>
  <rect x="444" y="190" width="60" height="20" fill="var(--polya)" opacity=".85"/>
  <text x="452" y="205" style="fill:#fff;font-size:11px">(A)n</text>
  <rect x="504" y="190" width="24" height="20" fill="var(--tsd)" opacity=".8"/>
  <text x="120" y="232" style="font-size:11px">TSD</text><text x="504" y="232" style="font-size:11px">TSD</text>
  <text x="560" y="205" style="font-size:11px">orphan TD = flank + (A)n only</text>
</svg>
"""


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    L = a.library
    src = tsv(os.path.join(L, "transduction_sources.tsv"))
    stats = tsv(os.path.join(L, "transduction_stats.tsv"))
    act = tsv(os.path.join(L, "active.tsv"))
    l1 = tsv(os.path.join(L, "l1_intact.tsv"))

    def cnt(f):
        return sum(1 for r in src if f(r))
    n_pub = cnt(lambda r: "published" in r["seed"])
    n_ref = cnt(lambda r: r["element_class"] == "L1" and r["reference"] == "yes")
    n_non = cnt(lambda r: r["element_class"] == "L1" and r["reference"] == "no")
    n_hs1 = cnt(lambda r: r["reference"] == "hs1_only")
    n_sva = cnt(lambda r: r["element_class"] == "SVA")
    n_hot = cnt(lambda r: r["hotness"] == "hot")
    n_strong = cnt(lambda r: r["hotness"] == "strong")
    n_absent = cnt(lambda r: r["element_class"] == "L1" and r["reference"] == "yes" and r["hs1_status"] == "insertion_point")
    n_in_hs1 = cnt(lambda r: r["hs1_status"] == "present_in_hs1")

    st = {(r["study"], r["metric"]): r for r in stats}
    allr = st[("all", "distal")]
    tdl = st[("all", "td_len")]

    rows_html = []
    pub = sorted([r for r in src if "published" in r["seed"]], key=lambda r: (-int(r["n_daughters"]), r["id"]))
    for r in pub:
        loc = "%s:%s-%s(%s)" % (r["hs1_chrom"], r["hs1_start"], r["hs1_end"], r["hs1_strand"]) if r["hs1_chrom"] != "." else "hg38 %s:%s(%s)" % (r["hg38_chrom"], r["hg38_start"], r["strand"])
        ev = r["evidence"].replace(";", ", ")
        rows_html.append(
            "<tr><td><code>%s</code></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>"
            "<td class=wrap>%s</td><td>%s</td><td><span class='pill %s'>%s</span></td></tr>" % (
                esc(r["id"]), esc(r["band"]), esc(loc), esc(r["subfamily"]), esc(r["ta_status"]),
                esc({"yes": "ref", "no": "non-ref", "hs1_only": "hs1 only"}.get(r["reference"], r["reference"]) +
                    (" (absent in hs1)" if r["reference"] == "yes" and r["hs1_status"] == "insertion_point" else "") +
                    (" (present in hs1)" if r["hs1_status"] == "present_in_hs1" else "")),
                esc(ev), esc(r["n_daughters"]), esc(r["hotness"]), esc(r["hotness"])))

    stat_rows = "".join(
        "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>" % (
            esc(r["study"]), esc({"td_len": "transduced length", "distal": "distal end from source 3' end"}[r["metric"]]),
            esc(r["n"]), esc(r["median"]), esc(r["p95"]), esc(r["p99"]), esc(r["max"]),
            esc("%.1f%%" % (100 * float(r["frac_le_10kb"]))), esc("%.1f%%" % (100 * float(r["frac_le_15kb"]))))
        for r in stats)

    n_ta = sum(1 for r in l1 if r["subfamily_call"] == "L1HS" and r["ta_status"] == "Ta")
    n_preta = sum(1 for r in l1 if r["subfamily_call"] == "L1HS" and r["ta_status"] == "preTa")

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Transduction source library</title>
<meta name="description" content="How PEAR-TREE's 3' transduction source library is built, looked up, and extended with novel sources.">
<style>{CSS}</style></head>
<body><main>
<p class="muted">PEAR-TREE · TPRT-hallmark annotator · <code>resources/rte_library/</code></p>
<h1>Transduction source library</h1>
<p>Which L1 (or SVA) produced an insertion? When the insertion carries a piece of unique genomic
sequence that sat downstream of its progenitor, the answer can be read off directly. This page
documents the library of candidate progenitors (“sources”) that annotate uses for that lookup,
how it was built, how a lookup works, and how a source that is not in the library is accepted
and added.</p>

<div class="grid">
 <div class="stat"><b>{len(src)}</b>sources in total</div>
 <div class="stat"><b>{n_pub}</b>published loci (3 studies + compendium)</div>
 <div class="stat"><b>{n_hot} / {n_strong}</b>hot (≥20) / strong (5–19 daughters)</div>
 <div class="stat"><b>{n_ref} / {n_non}</b>reference / non-reference L1</div>
 <div class="stat"><b>{n_hs1}</b>hs1-only full-length L1HS</div>
 <div class="stat"><b>{n_sva}</b>SVA_E/F sources (3' and 5' flanks)</div>
</div>

<h2>1. Transductions in one paragraph each</h2>
<div class="card fig">{FIG}</div>
<p><b>3' (partnered) transduction.</b> L1's own polyadenylation signal (AATAAA at the very 3' end,
overlapping the poly-A) is weak, so a fraction of transcripts read through into the downstream
genome and stop at the next usable pA site. Reverse transcription of such a transcript by TPRT
inserts L1 sequence (usually 5'-truncated) followed by the unique downstream segment and then the
poly-A (Moran et al. 1999; Goodier et al. 2000). The segment identifies the source.</p>
<p><b>Orphan transduction.</b> The same readthrough transcript, but reverse transcription stops
before reaching L1 sequence (extreme 5' truncation): the insertion is only the transduced segment
plus poly-A, flanked by a TSD — it looks like a “random” templated insertion unless the segment is
recognised as the downstream flank of a source (Tubio et al. 2014). Annotate
calls these <code>ORPHAN_TD</code>.</p>
<p><b>SVA 5' transduction.</b> SVAs can be transcribed from an upstream (external) promoter; the
transcript then starts upstream of the SVA hexamer and the insertion carries upstream genomic
sequence at its 5' end (Damert et al. 2009, MAST2-driven SVA_F1 group; ~8 % of SVAs). Annotate tags
<code>TD5P</code> using <code>flanks_5p_sva.fa.gz</code>. SVAs also make 3' transductions (Xing et
al. 2006), so SVA sources have 3' flanks too.</p>

<h2>2. How the library was built</h2>
<p>Rebuild: <code>bash tools/rte_library/fetch_inputs.sh</code> then <code>python
tools/rte_library/build.py</code> (details and provenance in <code>resources/rte_library/README.md</code>).</p>
<h3>Source lists</h3>
<ul>
<li><b>Nam et al. 2023</b> (Nature 617:540), Supplementary Table 4: 276 sources active in normal
colorectal epithelium / fibroblast clones and tumours, GRCh37, with strand, long-read subfamily
(Ta-0/Ta-1d/Ta-1nd/pre-Ta/PA2) and per-clone transduction counts.</li>
<li><b>Rodriguez-Martin et al. 2020</b> (Nat Genet 52:306, PCAWG), Supplementary Table 5: 124
germline sources with per-tumour counts; strands inferred from Supplementary Table 2 (which side
of the source each transduced segment lies on) — 43/44 agree with RepeatMasker.</li>
<li><b>Gardner et al. 2017</b> (MELT, Genome Res 27:1916), Supplemental Table S9 B (38 sources
with offspring counts, the Tubio 2014 counts and Brouha/Beck cell-culture activity) and S9 C
(literature compendium: Tubio, Beck, Brouha, Scott, Solyom, Helman, Evrony, …).</li>
<li><b>Not obtainable:</b> Tubio et al. 2014 Table S5 (not open access; PMC serves the supplement
only behind a browser proof-of-work check) — its counts enter through Gardner S9; Brouha et al. 2003
gives no coordinates in the open text; Damert et al. 2009 has no coordinate supplement. No
coordinates were invented.</li>
<li><b>Seeds</b> (no publication needed): every full-length (≥5.9 kb) L1HS in hg38 RepeatMasker,
L1Base-intact L1PA2/3, full-length L1HS present only in hs1 (T2T-resolved or present only in the
CHM13 haplotype), near-full-length SVA_E/F. These carry <code>hotness=candidate</code>.</li>
</ul>
<h3>Coordinates</h3>
<p>Published GRCh37 positions are lifted hg19→hg38→hs1 with the UCSC chains. Only same-chromosome
hits are accepted and an element's lift must be bracketed by its ±1 kb anchors, because the
over.chain files contain small paralog chains (we saw an SVA on chr6 “lift” onto an SVA on chr1).
A source is <i>reference</i> when a young L1 (L1HS/L1PA2–8, ≥4 kb) lies within 300 bp of the
published position in hg38; otherwise it is <i>non-reference</i> and the published position is
the insertion point. The flank of a non-reference source is ordinary reference sequence — exactly
what gets transduced — so it is extracted from hs1 at the lifted junction. {n_absent} reference
(hg38) sources are absent from CHM13 and {n_in_hs1} non-reference ones are present in it; both are
handled (<code>hs1_status</code>).</p>
<h3>Flank length</h3>
<p>From the published catalogues (n = {allr['n']} transductions with source and segment
coordinates): the transduced segment is short (median {tdl['median']} bp, 95th percentile
{tdl['p95']} bp), but its distal end — where the readthrough transcript was polyadenylated —
lies a median {allr['median']} bp and up to {allr['max']} bp downstream of the source's 3' end;
{100*float(allr['frac_le_10kb']):.1f}% lie within 10 kb (the window used by TraFiC-mem and MEIGA)
and {100*float(allr['frac_le_15kb']):.1f}% within <b>15 kb</b>, the window used here (Tubio et al.
2014 saw transductions reaching ~12 kb).</p>
<div class="tablewrap"><table><thead><tr><th>study</th><th>metric</th><th>n</th><th>median</th>
<th>p95</th><th>p99</th><th>max</th><th>≤10 kb</th><th>≤15 kb</th></tr></thead><tbody>{stat_rows}</tbody></table></div>
<h3>Transduction endpoints cluster at the source's pA sites</h3>
<p>Zumalave et al. (2024, bioRxiv; long reads) showed that the transductions of one source end at
the same place: all 83 transductions of the 2q24.1 source in tumour PD0270a end exactly 234 bp
downstream, at an alternative pA site — that source lacks the canonical L1 pA signal — whereas a
1p22.3 source with a strong canonical signal makes mostly solo-L1s and a few transductions ending
at an APARENT-predicted downstream site. Two consequences used here: (i) each source row lists the
first AATAAA/ATTAAA hexamers of its flank (<code>pas_hexamers_3p</code>) and whether it keeps the
canonical 3' AATAAA (<code>canonical_pas_3p</code>); a transduced segment ending 10–35 bp past one
of them is extra support; (ii) a source's activity is under-counted by transductions alone when
its canonical signal is strong — <code>hotness</code> is a lower bound.</p>
<h3>Validation</h3>
<p>383 of 400 randomly drawn published transduced segments (Rodriguez-Martin 2020) map onto
<code>flanks_3p.fa.gz</code> in sense orientation (1 antisense, 16 unmapped), median offset 316 bp
into the flank. Ta/pre-Ta typing of the 146 L1Base intact L1s (two 3'UTR diagnostic sites, L1.3
5931 = Boissinot et al. 2000's ACA/G site, and 5712) reproduces the long-read subfamily calls of Nam
2023 for 50/50 Ta and 28/29 pre-Ta sources; result: {n_ta} L1HS-Ta and {n_preta} L1HS-pre-Ta.
</p>

<h2>3. Source table</h2>
<p>Published loci, hottest first ({len(pub)} rows; seeds are in the TSV). Daughters are summed
across independent cohorts. Locus is hs1 (1-based, element; for an insertion point start = end =
0-based junction).</p>
<input type="search" id="q" placeholder="filter: band, id, study, hot…" aria-label="filter table">
<div class="tablewrap" style="max-height:70vh"><table id="t"><thead><tr><th>id</th><th>band (hg38)</th>
<th>hs1 locus</th><th>subfamily</th><th>Ta</th><th>status</th><th>evidence</th><th>daughters</th><th>hotness</th></tr></thead>
<tbody>{''.join(rows_html)}</tbody></table></div>

<h2>4. Lookup algorithm (annotate)</h2>
<ol>
<li><b>Extract the candidate segment.</b> From the insertion consensus (pooled reads, indel-aware):
the sequence between the element's 3' end and the poly-A (partnered), or the whole insert before
the poly-A when no element sequence is present (orphan). Discard ≤ 20 bp segments (too short to
place uniquely; these are handled as <code>TEMPLATED_LOCAL</code> or untemplated).</li>
<li><b>Known sources.</b> Align to <code>flanks_3p.fa.gz</code> (mappy/minimap2 <code>sr</code>, or
bowtie2 <code>--local</code>). Accept a hit if ≥ 30 bp aligned at ≥ 95 % identity, on the
<b>forward</b> strand of the flank (flanks are in element sense, so a genuine transduced segment is
always forward), outside soft-masked repeat for most of its length, and the poly-A follows the
segment's distal end in the insertion. Report <code>TD3P</code> + <code>TD3P_SOURCE=&lt;id&gt;</code>,
the flank offsets (start ≈ distance from the source's 3' end; end = transcript termination) and
whether the end sits 10–35 bp past a <code>pas_hexamers_3p</code> hexamer. Several equally good
sources (segmental duplications) → list all ids separated by <code>|</code>.</li>
<li><b>Element–source concordance.</b> For a partnered transduction the element part should be the
source's class (L1 tag ↔ L1 source). A tag in the flank of an SVA source after L1 sequence is
<code>CHIMERIC_ENDS</code>, not a transduction.</li>
<li><b>Unknown flank → novel source.</b> If nothing in the library matches, map the segment to hs1
(unique, MAPQ ≥ 20). Look upstream of the hit, strand-aware, for a source element: for a hit on
+ at position p, a + strand L1 whose 3' end lies in [p − 15 kb, p]; for a hit on −, a − strand L1
whose 5'-most base lies in [p, p + 15 kb]. Apply the rule in section 5. Without such an element, check
the cohort's own non-reference L1 calls; otherwise the segment is not a transduction
(<code>TEMPLATED_LOCAL</code> if it comes from ≤ 15 bp of the insertion site, else <code>UNKNOWN</code>).</li>
<li><b>SVA 5' transduction.</b> Sequence 5' of an SVA hexamer that aligns forward to
<code>flanks_5p_sva.fa.gz</code> → <code>TD5P</code>.</li>
</ol>

<h2>5. Accepting a novel source</h2>
<div class="card">
<p><b>Rule</b> (refines the SPEC “Novel source rule”). The transduced segment is unique on hs1
(MAPQ ≥ 20), not in <code>flanks_3p.fa.gz</code>, lies 0–15 kb downstream (strand-aware) of a
candidate element, and the insertion carries TPRT hallmarks (poly-A after the segment, TSD or EN
motif). The candidate element is scored with <code>common.cons_identity</code> (best of consensus-
in-element and element-in-consensus infix identity to the L1HS consensus):</p>
<ul>
<li><b>Tier A — credible:</b> reference L1 ≥ 5.5 kb (5'UTR promoter present) and identity ≥ 0.98,
or a full-length L1 insertion called elsewhere in the same cohort (non-reference; no sequence
needed). Report <code>TD3P_SOURCE=novel:&lt;hs1 chr:start-end(strand)&gt;</code>,
<code>NOVEL_SOURCE</code> and the identity.</li>
<li><b>Tier B — reasonably similar:</b> ≥ 5.5 kb and 0.95 ≤ identity &lt; 0.98. Report the same tags
(with the identity, so the score can down-weight it); append to the library only with ≥ 2
independent daughters (different insertion sites, or different patients).</li>
<li><b>Not a source:</b> identity &lt; 0.95 or &lt; 5.5 kb.</li>
</ul>
<p>Why these numbers (identity to our L1HS consensus, 40 random full-length copies per family):
L1HS median 0.993 (5th percentile 0.989); L1PA2 median 0.974 (0.963–0.985); L1PA3 median 0.949;
L1PA4 0.913; L1PA5 0.891. The 56 published reference sources with ≥ 1 daughter: median 0.995, 5th
percentile 0.991, minimum 0.966 (two L1PA2-annotated copies typed pre-Ta by long reads — the
reason tier B exists). Retrotransposition-competent L1s are L1HS (Brouha et al. 2003: ~80–100 per
genome, a handful hot); the youngest L1PA2 can still transduce when reactivated (Nam 2023 lists 11
PA2 sources).</p>
</div>
<h3>Adding it to the library</h3>
<p><code>tools/rte_library/add_source.py</code> re-checks the tier, refuses a duplicate (existing
source on the same strand with its 3' end within 1 kb), extracts the 15 kb downstream flank from
hs1 (soft-masked when the hs1 RepeatMasker track is given), appends one row to
<code>transduction_sources.tsv</code> (<code>seed=novel_accepted</code>, <code>notes=NOVEL_SOURCE;tier=…</code>)
and one record to <code>flanks_3p.fa.gz</code> (re-bgzipped and re-indexed), and refreshes
<code>manifest.tsv</code>. With the inputs fetched by <code>fetch_inputs.sh</code> into
<code>/data/rte_inputs</code>:</p>
<pre>python tools/rte_library/add_source.py --library resources/rte_library \\
    --hs1-2bit /data/rte_inputs/genomes/hs1.2bit \\
    --hs1-rmsk /data/rte_inputs/genomes/hs1.repeatMasker.out.gz \\
    --chrom chr5 --start 1000001 --end 1006030 --strand + \\
    --evidence "PEAR-TREE PD12345: 3 daughters (chr2:..., chr7:..., chr11:...)" --n-daughters 3</pre>
<p>For a cohort-called non-reference L1 use <code>--junction &lt;0-based offset&gt;</code> instead of
<code>--start/--end</code>. Rebuilding with <code>build.py</code> regenerates the library from the
published lists and seeds, so appended sources must be re-applied (keep the add_source commands with
the cohort's results).</p>

<h2>6. Files</h2>
<ul>
<li><code>transduction_sources.tsv</code> — one row per source: id, class, subfamily (rmsk + long-read
call in brackets), Ta status, hg38 band, reference status, hs1 and hg38 coordinates, strand and how
it was determined, published hg19 coordinates, L1Base id, identity to L1HS, canonical pA, evidence,
seed, daughters (total and per study), hotness, alt ids (e.g. LRE3), cell-culture activity, flank
names and pA hexamers.</li>
<li><code>flanks_3p.fa.gz</code> — element-sense downstream flanks; description
<code>hs1:chr:start-end(strand)</code>. Sources with unknown strand ({cnt(lambda r: r['strand'] == '.')} non-reference loci) have
two records, <code>&lt;id&gt;/+</code> and <code>&lt;id&gt;/-</code>.</li>
<li><code>flanks_5p_sva.fa.gz</code> — 5 kb upstream of each SVA source, ending at its 5' end.</li>
<li><code>active.tsv</code> — {len(act)} L1s regarded as active (all L1HS-class intact L1s + published
sources with ≥ 5 daughters) with identity to the consensus; used for <code>nearest_active</code> /
<code>element_identity</code>.</li>
</ul>
<p class="muted">References: Boissinot et al. 2000 Mol Biol Evol 17:915 · Brouha et al. 2003 PNAS
100:5280 · Damert et al. 2009 Genome Res 19:1992 · Gardner et al. 2017 Genome Res 27:1916 ·
Goodier et al. 2000 Hum Mol Genet 9:653 · Moran et al. 1999 Science 283:1530 · Nam et al. 2023
Nature 617:540 · Penzkofer et al. 2017 NAR 45:D68 · Rodriguez-Martin et al. 2020 Nat Genet 52:306 ·
Tubio et al. 2014 Science 345:1251343 · Xing et al. 2006 PNAS
103:17608 · Zumalave et al. 2024 bioRxiv 2024.08.27.</p>
</main>
<script>
(function(){{var q=document.getElementById('q'),rows=document.querySelectorAll('#t tbody tr');
q.addEventListener('input',function(){{var v=q.value.toLowerCase();rows.forEach(function(r){{
r.style.display=r.textContent.toLowerCase().indexOf(v)>=0?'':'none';}});}});}})();
</script>
</body></html>
"""
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write(doc)


if __name__ == "__main__":
    main()
