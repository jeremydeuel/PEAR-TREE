#!/usr/bin/env python3
"""Generate docs/insertion_types.html — the PEAR-TREE retrotransposition insertion-type catalogue.

The page is static and self-contained: every figure is an inline SVG built here from a small
block/read description so all figures share one geometry and one colour legend (CSS custom
properties on :root, with a dark-mode override).  Re-run after editing:

    python3 docs/src/make_insertion_types.py

Vocabulary (element / structure / tags) follows plans/tprt_hallmarks/SPEC.md.
"""
from __future__ import annotations

import html
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "insertion_types.html"))

# ----------------------------------------------------------------------------------------
# Figure engine
# ----------------------------------------------------------------------------------------
W = 800            # SVG viewBox width
X0, X1 = 10, 790   # drawable range
READ_LEN = 62      # px, one short read
ROW_H = 13         # read row pitch
READ_H = 8
ALIGNED = {"flank", "flank2", "tsd", "al"}   # pieces that align at this locus (drawn as such)

# kind -> (css class, bar height)
KIND = {
    "flank": ("k-flank", 10),
    "flank2": ("k-flank2", 10),
    "tsd": ("k-tsd", 22),
    "l1": ("k-l1", 26),
    "l1inv": ("k-l1inv", 26),
    "alu": ("k-alu", 26),
    "sva": ("k-sva", 26),
    "polya": ("k-polya", 16),
    "td": ("k-td", 20),
    "exon": ("k-exon", 20),
    "exonb": ("k-exonb", 20),
    "intron": ("k-intron", 12),
    "local": ("k-local", 18),
    "gap": ("k-gap", 18),
}
ARROW_KINDS = {"l1", "l1inv", "alu", "sva", "td", "local", "exon", "exonb"}


def esc(s: str) -> str:
    return html.escape(s, quote=True)


def text_w(s: str, size: float = 10.5) -> float:
    return 0.56 * size * len(s)


class Block:
    def __init__(self, kind, w, label="", dir=None, inner=None):
        self.kind, self.w, self.label, self.dir, self.inner = kind, w, label, dir, inner
        self.x = 0.0
        self.px = 0.0


class Lane:
    """One horizontal allele (alt allele, reference or derivative chromosome) + its reads."""

    def __init__(self, title, blocks, scale_to=None):
        self.title = title
        self.blocks = [Block(*b) if isinstance(b, tuple) else b for b in blocks]
        total = sum(b.w for b in self.blocks)
        span = (X1 - X0) if scale_to is None else scale_to
        x = X0
        for b in self.blocks:
            b.x = x
            b.px = b.w / total * span
            x += b.px
        self.reads = []   # list of fragments; fragment = list of reads; read = list of pieces
        self.marks = []   # (x, text)

    # geometry helpers -------------------------------------------------------------------
    def bx(self, i):
        return self.blocks[i].x

    def ex(self, i):
        return self.blocks[i].x + self.blocks[i].px

    def mid(self, i):
        return self.blocks[i].x + self.blocks[i].px / 2

    def mark(self, x, text):
        self.marks.append((x, text))

    # read construction ------------------------------------------------------------------
    def _pieces(self, a, b):
        out = []
        for blk in self.blocks:
            lo, hi = max(a, blk.x), min(b, blk.x + blk.px)
            if hi - lo > 0.5 and blk.kind != "gap":   # absent sequence: the read simply skips it
                out.append((lo, hi, blk.kind))
        return out

    def read(self, a, length=READ_LEN, strand="+"):
        a = max(X0, min(a, X1 - length))
        return {"pieces": self._pieces(a, a + length), "strand": strand}

    def frag(self, *reads):
        self.reads.append(list(reads))

    def split(self, i, f=0.5, mate_dx=None):
        """Read straddling the boundary in front of block i (fraction f of read left of it)."""
        a = self.bx(i) - f * READ_LEN
        r = self.read(a)
        if mate_dx is None:
            self.frag(r)
        else:
            m = self.read(a + mate_dx, strand="-" if mate_dx > 0 else "+")
            if mate_dx < 0:
                r["strand"] = "-"
            self.frag(r, m)

    def pair(self, a, b):
        r1 = self.read(a, strand="+")
        r2 = self.read(b, strand="-")
        self.frag(r1, r2)

    def custom(self, *reads):
        """reads: lists of (x0, x1, kind) pieces given explicitly."""
        rs = [{"pieces": list(p), "strand": "+" if k == 0 else "-"} for k, p in enumerate(reads)]
        self.frag(*rs)


def _arrow(x, w, y, h, d):
    tip = min(12.0, w / 3)
    t = h / 2
    if d == "-":
        pts = [(x, y), (x + tip, y - t), (x + w, y - t), (x + w, y + t), (x + tip, y + t)]
    else:
        pts = [(x, y - t), (x + w - tip, y - t), (x + w, y), (x + w - tip, y + t), (x, y + t)]
    return " ".join(f"{px:.1f},{py:.1f}" for px, py in pts)


def render_lane(lane: Lane, y0: float) -> tuple[str, float]:
    out = []
    out.append(f'<text class="lt" x="{X0}" y="{y0 + 11:.1f}">{esc(lane.title)}</text>')
    # marks row(s)
    mark_rows = []
    mk_svg = []
    for x, t in lane.marks:
        tw = text_w(t, 11)
        lo, hi = x - tw / 2, x + tw / 2
        if lo < X0:
            lo, hi = X0, X0 + tw
        if hi > X1:
            lo, hi = X1 - tw, X1
        r = 0
        while r < len(mark_rows) and any(not (hi + 6 < a or lo - 6 > b) for a, b in mark_rows[r]):
            r += 1
        if r == len(mark_rows):
            mark_rows.append([])
        mark_rows[r].append((lo, hi))
        mk_svg.append((x, lo, t, r))
    n_mrows = max(1, len(mark_rows))
    top = y0 + 20 + 13 * n_mrows
    cy = top + 16
    for x, lo, t, r in mk_svg:
        ty = y0 + 29 + 13 * r
        out.append(f'<text class="mk" x="{lo:.1f}" y="{ty:.1f}">{esc(t)}</text>')
        out.append(f'<line class="mkl" x1="{x:.1f}" y1="{ty + 2:.1f}" x2="{x:.1f}" y2="{cy - 14:.1f}"/>')
    # blocks
    for b in lane.blocks:
        cls, h = KIND[b.kind]
        if b.kind == "gap":
            out.append(f'<line class="gapl" x1="{b.x:.1f}" y1="{cy:.1f}" x2="{b.x + b.px:.1f}" y2="{cy:.1f}"/>')
            out.append(f'<rect class="gapr" x="{b.x + 1:.1f}" y="{cy - h / 2:.1f}" width="{max(1, b.px - 2):.1f}" height="{h}" rx="3"/>')
        elif b.kind in ARROW_KINDS and b.dir in ("+", "-"):
            out.append(f'<polygon class="{cls}" points="{_arrow(b.x, b.px, cy, h, b.dir)}"/>')
            if b.px > 70 and b.kind in {"l1", "l1inv", "alu", "sva"}:
                l5, l3 = ("5′", "3′") if b.dir == "+" else ("3′", "5′")
                out.append(f'<text class="it" x="{b.x + 5:.1f}" y="{cy + 3.5:.1f}">{l5}</text>')
                out.append(f'<text class="it" x="{b.x + b.px - 22:.1f}" y="{cy + 3.5:.1f}">{l3}</text>')
        else:
            out.append(f'<rect class="{cls}" x="{b.x:.1f}" y="{cy - h / 2:.1f}" width="{b.px:.1f}" height="{h}"/>')
        inner = b.inner
        if inner is None:
            inner = {"tsd": "TSD", "polya": "A(n)"}.get(b.kind, "")
        if inner and text_w(inner, 10) < b.px - 3:
            c = "itd" if b.kind in {"tsd", "exon", "exonb", "intron"} else "it"
            out.append(f'<text class="{c}" text-anchor="middle" x="{b.x + b.px / 2:.1f}" y="{cy + 3.3:.1f}">{esc(inner)}</text>')
    # labels below, greedy rows
    rows = []
    lab_svg = []
    for b in lane.blocks:
        if not b.label:
            continue
        tw = text_w(b.label, 12)
        c = b.x + b.px / 2
        lo, hi = c - tw / 2, c + tw / 2
        if lo < X0:
            lo, hi = X0, X0 + tw
        if hi > X1:
            lo, hi = X1 - tw, X1
        r = 0
        while r < len(rows) and any(not (hi + 6 < a or lo - 6 > bb) for a, bb in rows[r]):
            r += 1
        if r == len(rows):
            rows.append([])
        rows[r].append((lo, hi))
        lab_svg.append((lo, c, r, b.label))
    lab_top = cy + 17
    for lo, c, r, t in lab_svg:
        ty = lab_top + 13 + 14 * r
        if r > 0:
            out.append(f'<line class="mkl" x1="{c:.1f}" y1="{cy + 14:.1f}" x2="{c:.1f}" y2="{ty - 9:.1f}"/>')
        out.append(f'<text class="bl" x="{lo:.1f}" y="{ty:.1f}">{esc(t)}</text>')
    y = lab_top + 14 * len(rows) + 10
    # reads: pack fragments into rows
    if lane.reads:
        out.append(f'<text class="lt2" x="{X0}" y="{y + 4:.1f}">short reads</text>')
        y += 12
        frows = []
        placed = []
        for fr in sorted(lane.reads, key=lambda f: min(p[0] for r in f for p in r["pieces"])):
            lo = min(p[0] for r in fr for p in r["pieces"])
            hi = max(p[1] for r in fr for p in r["pieces"])
            r = 0
            while r < len(frows) and frows[r] + 6 > lo:
                r += 1
            if r == len(frows):
                frows.append(-1e9)
            frows[r] = hi
            placed.append((fr, r))
        for fr, r in placed:
            ry = y + r * ROW_H
            spans = []
            for rd in fr:
                ps = rd["pieces"]
                a, b = ps[0][0], ps[-1][1]
                spans.append((a, b))
            if len(fr) == 2:
                (a1, b1), (a2, b2) = sorted(spans)
                disc = any(all(p[2] not in ALIGNED for p in rd["pieces"]) for rd in fr)
                if a2 > b1:
                    out.append(f'<line class="{"pld" if disc else "pl"}" x1="{b1:.1f}" y1="{ry + READ_H / 2:.1f}" x2="{a2:.1f}" y2="{ry + READ_H / 2:.1f}"/>')
            for rd in fr:
                for a, b, k in rd["pieces"]:
                    cls = "r-al" if k in ("flank", "flank2", "al") else f"r-{k}"
                    out.append(f'<rect class="{cls}" x="{a:.1f}" y="{ry:.1f}" width="{b - a:.1f}" height="{READ_H}"/>')
                a, b = rd["pieces"][0][0], rd["pieces"][-1][1]
                if rd["strand"] == "+":
                    out.append(f'<polygon class="rh" points="{b:.1f},{ry:.1f} {b + 4:.1f},{ry + READ_H / 2:.1f} {b:.1f},{ry + READ_H:.1f}"/>')
                else:
                    out.append(f'<polygon class="rh" points="{a:.1f},{ry:.1f} {a - 4:.1f},{ry + READ_H / 2:.1f} {a:.1f},{ry + READ_H:.1f}"/>')
        y += len(frows) * ROW_H + 6
    return "\n".join(out), y


def figure(fid, lanes, caption):
    parts = []
    y = 4.0
    for k, lane in enumerate(lanes):
        if k:
            parts.append(f'<line class="sep" x1="{X0}" y1="{y + 2:.1f}" x2="{X1}" y2="{y + 2:.1f}"/>')
            y += 8
        s, y = render_lane(lane, y)
        parts.append(s)
    h = y + 4
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {h:.0f}" role="img" '
           f'aria-labelledby="{fid}-t"><title id="{fid}-t">{esc(caption)}</title>\n' + "\n".join(parts) + "\n</svg>")
    return f'<figure id="{fid}">{svg}<figcaption>{caption}</figcaption></figure>'


# ----------------------------------------------------------------------------------------
# Figure definitions
# ----------------------------------------------------------------------------------------
L = READ_LEN
FR = 175   # px between read starts of a pair (schematic fragment)


def std_reads(lane, i5, i3, polya_idx=None, disc5=True, disc3=True):
    """Split reads at 5' (front of block i5) and 3' (front of block i3) junctions + discordants."""
    lane.split(i5, 0.55)
    lane.split(i5, 0.3, mate_dx=-FR)
    lane.split(i3, 0.5)
    lane.split(i3, 0.7, mate_dx=FR)
    if disc5:
        lane.pair(lane.bx(i5) - L - 40, lane.bx(i5) + 70)
    if disc3:
        lane.pair(lane.bx(i3) - L - 60, lane.bx(i3) + 40)
    if polya_idx is not None:
        lane.pair(lane.bx(polya_idx) + 2, lane.bx(polya_idx) + FR)


def fig_solo_full():
    ln = Lane("insertion allele — full-length L1 (FULL_LENGTH)", [
        ("flank", 14), ("tsd", 4), ("l1", 56, "L1HS ~6.0 kb, intact 5′ UTR", "+"), ("polya", 5),
        ("tsd", 4), ("flank", 14)])
    ln.mark(ln.bx(2), "5′ junction"); ln.mark(ln.ex(3), "3′ junction · EN nick TTTT/AA")
    std_reads(ln, 2, 4, polya_idx=3)
    return ln


def fig_solo_trunc():
    ln = Lane("insertion allele — 5′-truncated L1 (TRUNCATED_5P)", [
        ("flank", 22), ("tsd", 5), ("l1", 30, "L1 3′ end only (median ~0.5 kb)", "+"), ("polya", 7),
        ("tsd", 5), ("flank", 22)])
    ln.mark(ln.bx(2), "5′ junction = truncation point")
    ln.mark(ln.ex(3), "3′ junction (poly-A)")
    std_reads(ln, 2, 4, polya_idx=3)
    return ln


def fig_solo_inv():
    ln = Lane("insertion allele — 5′-inverted L1, twin priming (INVERTED_5P)", [
        ("flank", 14), ("tsd", 4), ("l1inv", 22, "inverted 5′ part (internal primer)", "-"),
        ("gap", 3, "i-del / i-dup"),
        ("l1", 30, "forward part (poly-dT primer)", "+"), ("polya", 5), ("tsd", 4), ("flank", 14)])
    ln.mark(ln.bx(2), "5′ junction (clip = L1, reverse strand)")
    ln.mark(ln.mid(3), "inversion junction (median 5.15 kb on L1)")
    ln.mark(ln.ex(5), "3′ junction")
    ln.split(2, 0.55)
    ln.split(2, 0.3, mate_dx=-FR)
    ln.split(3, 0.5)
    ln.split(5, 0.5)
    ln.split(6, 0.6, mate_dx=FR)
    ln.pair(ln.bx(2) - L - 40, ln.bx(2) + 60)
    ln.pair(ln.bx(6) - L - 50, ln.bx(6) + 40)
    return ln


def fig_td_partnered():
    ln = Lane("insertion allele — partnered 3′ transduction (element L1, tag TD3P)", [
        ("flank", 12), ("tsd", 4), ("l1", 26, "L1 (5′-truncated)", "+"), ("polya", 4, "", None, "A"),
        ("td", 18, "transduced 3′ flank of source L1", "+"), ("polya", 5), ("tsd", 4), ("flank", 12)])
    ln.mark(ln.mid(3), "source L1 poly-A (short / absent)")
    ln.mark(ln.ex(5), "3′ junction: clip = A(n) then source flank")
    std_reads(ln, 2, 6, polya_idx=5)
    ln.pair(ln.bx(4) + 5, ln.bx(4) + FR - 20)
    return ln


def fig_td_orphan():
    ln = Lane("insertion allele — orphan 3′ transduction (element ORPHAN_TD)", [
        ("flank", 20), ("tsd", 5), ("td", 26, "unique 3′ flank of a source L1 (no L1 sequence)", "+"),
        ("polya", 6), ("tsd", 5), ("flank", 20)])
    ln.mark(ln.bx(2), "5′ junction: clip maps uniquely (MAPQ high) to source flank")
    std_reads(ln, 2, 4, polya_idx=3)
    return ln


def fig_alu():
    ln = Lane("insertion allele — Alu (element ALU) ~0.3 kb, shorter than a fragment", [
        ("flank", 26), ("tsd", 5), ("alu", 20, "AluY (left/right monomer)", "+"), ("polya", 6),
        ("tsd", 5), ("flank", 26)])
    std_reads(ln, 2, 4)
    ln.pair(ln.bx(1) - L - 30, ln.ex(4) + 20)   # SPAN-like pair across the whole insertion
    return ln


def fig_sva():
    ln = Lane("insertion allele — SVA with 5′ transduction (element SVA, tag TD5P)", [
        ("flank", 12), ("tsd", 4), ("td", 12, "upstream flank of source SVA", "+"),
        ("sva", 40, "SVA: CCCTCT hexamer · Alu-like · VNTR · SINE-R", "+"), ("polya", 5), ("tsd", 4),
        ("flank", 12)])
    ln.mark(ln.bx(2), "5′ junction: clip = source upstream (often a gene exon 1)")
    std_reads(ln, 2, 5, polya_idx=4)
    return ln


def fig_pseudogene():
    ln = Lane("insertion allele — processed pseudogene (element PSEUDOGENE, tag EXON_JUNCTION)", [
        ("flank", 14), ("tsd", 4), ("exon", 14, "exon 1", "+"), ("exonb", 12, "exon 2", "+"),
        ("exon", 16, "exon 3 (3′ UTR)", "+"), ("polya", 6), ("tsd", 4), ("flank", 14)])
    ln.mark(ln.bx(3), "exon–exon junction (no intron)")
    ln.mark(ln.bx(4), "exon–exon junction")
    std_reads(ln, 2, 6, polya_idx=5)
    ln.split(3, 0.5)
    ln.split(4, 0.45)
    return ln


def fig_polya_only():
    ln = Lane("insertion allele — solitary poly(A/T) (element POLYA_ONLY)", [
        ("flank", 30), ("tsd", 6), ("polya", 10, "", None, "A(n) only"), ("tsd", 6), ("flank", 30)])
    ln.mark(ln.bx(2), "5′ junction: clip = A(n)")
    ln.mark(ln.ex(2), "3′ junction: clip = A(n)")
    ln.split(2, 0.6)
    ln.split(2, 0.75, mate_dx=-FR)
    ln.split(3, 0.4)
    ln.split(3, 0.25, mate_dx=FR)
    return ln


def fig_l1_del():
    ref = Lane("reference", [("flank", 20, "A"), ("gap", 34, "deleted segment 0.5 kb – 53 Mb"),
                             ("flank", 4, "target"), ("flank", 20, "B")])
    alt = Lane("insertion allele — L1-mediated deletion (tag L1_MED_DELETION)", [
        ("flank", 20, "A (upstream DSB)"), ("l1", 30, "L1 (no 5′ TSD)", "+"), ("polya", 6),
        ("flank", 22, "B (original target site)")])
    alt.mark(alt.bx(1), "5′: 1–5 bp microhomology, no TSD")
    alt.mark(alt.ex(2), "3′: EN motif, poly-A")
    std_reads(alt, 1, 3, polya_idx=2)
    return [ref, alt]


def fig_foldback():
    alt = Lane("derivative chromosome — L1-bridged fold-back inversion (tag L1_MED_DUPLICATION)", [
        ("flank", 30, "chr arm →"), ("l1", 18, "short L1", "+"), ("polya", 5),
        ("flank2", 30, "← same arm, inverted copy (sister chromatid)")])
    alt.mark(alt.bx(1), "junction 1")
    alt.mark(alt.ex(2), "junction 2 (poly-A)")
    alt.split(1, 0.5)
    alt.split(1, 0.25, mate_dx=-FR)
    alt.split(3, 0.5)
    alt.split(3, 0.75, mate_dx=FR)
    alt.pair(alt.bx(1) - L - 50, alt.bx(1) + 40)
    alt.pair(alt.ex(2) - 20 - L, alt.ex(2) + 70)
    return [alt]


def fig_rtrg():
    lanes = []
    for title, lab in [("deletion-like [+/−] (tag L1_MED_DELETION)", "B: downstream, intervening DNA lost"),
                       ("duplication-like [−/+] (tag L1_MED_DUPLICATION)", "B: upstream, intervening DNA duplicated"),
                       ("inversion-like [+/+] or [−/−] (out of scope: bridge)", "← B, inverted orientation")]:
        ln = Lane(title, [("flank", 26, "A"), ("l1", 22, "L1 bridge", "+"), ("polya", 5),
                          ("flank2", 26, lab)])
        ln.split(1, 0.5)
        ln.split(3, 0.5)
        ln.pair(ln.bx(1) - L - 40, ln.bx(1) + 40)
        ln.pair(ln.ex(2) - 30 - L, ln.ex(2) + 60)
        lanes.append(ln)
    return lanes


def fig_recip_tl():
    d1 = Lane("der(A) — bridge 1", [("flank", 22, "chr A (5′ side)"), ("tsd", 4, "", None, "a"),
                                    ("l1", 24, "L1 insertion #1", "+"), ("polya", 5), ("tsd", 4, "", None, "b"),
                                    ("flank2", 22, "chr B")])
    d2 = Lane("der(B) — bridge 2", [("flank2", 22, "chr B"), ("tsd", 4, "", None, "b′"),
                                    ("l1", 24, "L1 insertion #2", "+"), ("polya", 5), ("tsd", 4, "", None, "a′"),
                                    ("flank", 22, "chr A (3′ side)")])
    d1.mark(d1.mid(1), "heteroduplicated TSD: copies a / a′ on different derivatives")
    for d in (d1, d2):
        d.split(2, 0.5)
        d.split(4, 0.5)
        d.pair(d.bx(2) - L - 40, d.bx(2) + 50)
        d.pair(d.ex(3) - 20 - L, d.ex(3) + 60)
    return [d1, d2]


def fig_chimeric_bridge():
    d1 = Lane("bridge with both TPRT ends (two poly-A/T tails, two EN motifs)", [
        ("flank", 18, "chr A"), ("polya", 5, "", None, "T(n)"), ("td", 8, "TD #2", "-"),
        ("l1inv", 18, "L1 #2 (−)", "-"), ("l1", 18, "L1 #1", "+"), ("td", 8, "TD #1", "+"), ("polya", 5),
        ("flank2", 18, "chr B")])
    d2 = Lane("partner bridge: remnants of the same two L1s, no poly-A, no EN motif", [
        ("flank2", 22, "chr B′"), ("l1", 16, "L1 #1 piece", "+"), ("l1inv", 16, "L1 #2 piece", "-"),
        ("flank", 22, "chr A′")])
    d1.split(2, 0.5); d1.split(7, 0.5)
    d2.split(1, 0.5); d2.split(3, 0.5); d2.split(2, 0.5)
    return [d1, d2]


def fig_recip_inv():
    b1 = Lane("bridge 1 — two 5′-truncated L1s in opposite orientation, joined at their 5′ ends", [
        ("flank", 18, "chr5 a"), ("polya", 5, "", None, "T(n)"), ("l1inv", 20, "L1 (−)", "-"),
        ("l1", 20, "L1 (+)", "+"), ("polya", 5), ("flank2", 18, "chr5 c (inverted)")])
    b2 = Lane("bridge 2 (Mb away) — interstitial L1 fragment within a heteroduplicated target site", [
        ("flank", 22, "chr5 b"), ("tsd", 4), ("l1", 18, "L1 fragment (other half)", "+"), ("tsd", 4),
        ("flank2", 22, "chr5 d")])
    b1.split(2, 0.5); b1.split(5, 0.5); b1.split(3, 0.5)
    b2.split(2, 0.5); b2.split(4, 0.5)
    return [b1, b2]


def fig_switch():
    a = Lane("twin priming + 5′ switching (INVERTED_5P_SWITCH)", [
        ("flank", 12), ("tsd", 4), ("l1inv", 16, "inverted 5′ part", "-"), ("l1", 10, "internal inversion", "+"),
        ("l1inv", 8, "", "-"), ("l1", 28, "forward part", "+"), ("polya", 5), ("tsd", 4), ("flank", 12)])
    a.mark(a.bx(3), "switch junctions")
    a.mark(a.bx(5), "twin-priming junction")
    a.split(2, 0.5); a.split(3, 0.5); a.split(5, 0.5); a.split(7, 0.5)
    b = Lane("twin priming + 3′ switching (INVERTED_5P_SWITCH)", [
        ("flank", 12), ("tsd", 4), ("l1inv", 18, "inverted 5′ part", "-"),
        ("l1", 12, "non-inverted piece from a distant L1 region", "+"),
        ("l1", 26, "forward part", "+"), ("polya", 5), ("tsd", 4), ("flank", 12)])
    b.mark(b.bx(3), "inversion junction")
    b.mark(b.bx(4), "discontinuity in L1 coordinates")
    b.split(2, 0.5); b.split(3, 0.5); b.split(4, 0.5); b.split(6, 0.5)
    return [a, b]


def fig_templated():
    ln = Lane("insertion allele — L1 with a templated local segment (tag TEMPLATED_LOCAL)", [
        ("flank", 18, "site (template origin nearby)"), ("tsd", 4), ("local", 8, "local copy (< 250 bp)", "-"),
        ("l1", 30, "L1", "+"), ("polya", 5), ("tsd", 4), ("flank", 18)])
    ln.mark(ln.mid(2), "clip maps back next to the site (often inverted)")
    ln.split(2, 0.6)
    ln.split(3, 0.5)
    ln.split(5, 0.5)
    ln.split(5, 0.7, mate_dx=FR)
    ln.pair(ln.bx(2) - L - 40, ln.bx(3) + 30)
    return ln


def fig_premrna():
    ln = Lane("insertion allele — co-inserted local pre-mRNA (tag PREMRNA_COINSERT)", [
        ("flank", 14, "site, inside/near a transcribed gene"), ("tsd", 4), ("exon", 8, "exon", "+"),
        ("intron", 10, "intron", None, "intron"), ("l1", 26, "L1", "+"), ("polya", 5), ("tsd", 4),
        ("flank", 14)])
    ln.mark(ln.bx(4), "RT template switch L1 RNA → nearby pre-mRNA")
    ln.split(2, 0.6); ln.split(3, 0.5); ln.split(4, 0.5); ln.split(6, 0.5)
    ln.pair(ln.bx(2) - L - 30, ln.bx(3) + 10)
    return ln


def fig_foldback_5p():
    ln = Lane("insertion allele — fold-back inverted duplication 5′ of the site (tag FOLDBACK_INVDUP_5P)", [
        ("flank", 14, "upstream segment S →"), ("tsd", 4), ("local", 10, "S′ (inverted copy, 52–220 bp)", "-"),
        ("l1", 26, "L1", "+"), ("polya", 5), ("tsd", 4), ("flank", 14)])
    ln.mark(ln.bx(2), "clip = reverse complement of adjacent reference")
    ln.mark(ln.bx(3), "L1 continues beyond the fold-back")
    ln.split(2, 0.55); ln.split(3, 0.5); ln.split(5, 0.5)
    ln.split(5, 0.7, mate_dx=FR)
    return ln


def fig_tsd_del():
    ref = Lane("reference", [("flank", 30), ("gap", 6, "target (deleted)"), ("flank", 30)])
    alt = Lane("insertion allele — target-site deletion (tag TSD_DELETION)", [
        ("flank", 30), ("l1", 26, "L1", "+"), ("polya", 6), ("flank", 30)])
    alt.mark(alt.bx(1), "no duplicated sequence: left and right flanks abut the insert")
    std_reads(alt, 1, 3)
    return [ref, alt]


def fig_en_indep():
    alt = Lane("insertion allele — EN-independent insertion (tag EN_INDEPENDENT)", [
        ("flank", 30, "pre-existing DSB / dysfunctional telomere"), ("l1", 26, "L1, truncated at 5′ AND 3′ end", "+"),
        ("flank", 30)])
    alt.mark(alt.bx(1), "no TSD")
    alt.mark(alt.ex(1), "no poly-A, no EN motif")
    alt.split(1, 0.5); alt.split(2, 0.5)
    alt.pair(alt.bx(1) - L - 40, alt.bx(1) + 50)
    alt.pair(alt.ex(1) - 40 - L, alt.ex(1) + 50)
    return [alt]


def fig_art_chimera():
    a = Lane("library molecule only — ligation / PCR / MDA chimera", [
        ("flank", 36, "genomic fragment from locus X"), ("l1", 30, "L1 piece from elsewhere (any coordinate)", "+"),
        ])
    a.mark(a.bx(1), "one junction · no TSD · no EN motif")
    a.split(1, 0.5)
    b = Lane("A-tailed ligation chimera", [
        ("flank", 36, "fragment from locus X"), ("polya", 4, "", None, "A"), ("l1", 30, "L1 or other fragment", "+"),
        ])
    b.mark(b.bx(1), "short A-tail at the join mimics a TPRT poly-A")
    b.split(1, 0.6)
    return [a, b]


def fig_art_slippage():
    ref = Lane("reference (no insertion) — Alu with its own A-tract", [
        ("flank", 24), ("alu", 22, "reference AluY", "+"), ("polya", 8, "", None, "ref A(n)"), ("flank", 26)])
    x = ref.ex(2)
    ref.mark(x, "soft-clips pile up here; clip = extra A")
    for c, a in ((18, 44), (26, 36), (14, 48)):
        ref.custom([(x - c, x, "polya"), (x, x + a, "al")])
    return [ref]


def fig_art_mismap():
    ref = Lane("reference — L1HS copy P (paralog of an active, non-reference copy Q)", [
        ("flank", 26, "locus of P"), ("l1", 30, "reference L1HS P", "+"), ("polya", 5), ("flank", 26)])
    x = ref.ex(2)
    ref.custom([(x - 50, x, "al"), (x, x + 12, "clip2")])
    ref.custom([(x - 40, x, "al"), (x, x + 22, "clip2")])
    ref.custom([(x - 34, x, "al"), (x, x + 28, "clip2")])
    ref.mark(x, "clip = flank of Q (reads of Q placed on P)")
    return [ref]


def fig_art_foldback():
    ref = Lane("reference / library — inverted repeat or hairpin", [
        ("flank", 30), ("local", 10, "arm →", "+"), ("flank", 8), ("local", 10, "← arm", "-"), ("flank", 30)])
    x = ref.bx(1)
    ref.custom([(x - 40, x, "flank"), (x, x + 18, "local")])
    ref.custom([(x - 30, x, "flank"), (x, x + 25, "local")])
    ref.mark(x, "clip = reverse complement of adjacent reference; no element, no poly-A")
    return [ref]


# ----------------------------------------------------------------------------------------
# Content
# ----------------------------------------------------------------------------------------
IN, PART, OUT_, ART = "in", "partial", "out", "art"
BADGE = {IN: "in scope", PART: "partly in scope", OUT_: "out of scope this round", ART: "artefact — score negative"}


def pts(*items):
    return "<ul class='pts'>" + "".join(f"<li>{i}</li>" for i in items) + "</ul>"


SECTIONS = []


def sec(sid, num, title, scope, body, fig, detect, freq, vocab, signature):
    SECTIONS.append(dict(sid=sid, num=num, title=title, scope=scope, body=body, fig=fig, detect=detect,
                         freq=freq, vocab=vocab, signature=signature))


sec("t1", "1", "Solo-L1 (full-length, 5′-truncated, 5′-inverted)", IN,
    """<p>An L1 RNA is reverse-transcribed at its own insertion site by target-primed reverse transcription
    (TPRT): the L1 endonuclease nicks the bottom strand at a T-rich motif (consensus 5′-TTTT/AA-3′), the
    exposed 3′-OH primes cDNA synthesis on the RNA poly(A) tail, and a staggered second nick creates the
    target-site duplication (TSD). Because reverse transcription usually stops early, almost all somatic
    copies keep only the 3′ end of L1 (5′ truncation). In <b>twin priming</b> the second DNA end of the
    break anneals internally on the same RNA and primes a second cDNA, so the 5′ part is inverted relative
    to the 3′ part; at the inversion point a short internal deletion (poly-dT cDNA stopped short of the
    internal primer) or duplication (it ran past it) is common.</p>""",
    [figure("f1a", [fig_solo_full()], "Full-length L1: TSD on both sides, poly-A at the 3′ junction, the 5′ clip is L1 position ~1."),
     figure("f1b", [fig_solo_trunc()], "5′-truncated L1: the 5′ clip starts at an internal L1 coordinate (the truncation point)."),
     figure("f1c", [fig_solo_inv()], "5′-inverted L1 (twin priming): the 5′ clip and its mates read L1 on the opposite strand to the 3′ end.")],
    pts("discovery: CLIP reads at both junctions; POLYA reads placed by their mate on the 3′ side; DISC anchors whose mates hit L1",
        "annotate: element <code>L1</code>; <code>structure</code> from the 5′ clip consensus — <code>FULL_LENGTH</code> (clip reaches L1 5′ UTR start), <code>TRUNCATED_5P</code>, <code>INVERTED_5P</code> (5′ clip/mates on the opposite strand), else <code>5P_UNRESOLVED</code>",
        "score: TSD 4–25 bp with identity; poly-A ≥ 10 on the strand-consistent side; EN motif at the nick; 5′/3′ class concordant; identity to an active L1HS; twin-priming junction ≥ 590 bp into L1; ≥ 2 independent fragments per junction pooled over colonies"),
    "Zumalave 2026 (10 high-rate tumours, 6,418 events): 56% of all events (3,611); of 5′-truncated solo-L1s 49% 5′-deletion only and 51% 5′-inverted; full-length only 0.3% (10/3,611). Of inverted copies 66% carry an internal deletion (median 14.5 bp), 17% an internal duplication (median 22 bp). PCAWG: 14,967/19,166 (78%). Nam 2023 normal colon: 89% solo-L1 (1,063/1,198), 29.5% with a short 5′ inversion. Tubio 2014: ~5% of non-inverted cancer L1s full-length.",
    "element=L1; structure=FULL_LENGTH | TRUNCATED_5P | INVERTED_5P | 5P_UNRESOLVED",
    "L1 clip at 5′ and A(n) clip at 3′, TSD between them; 5′ clip strand flips for inverted copies")

sec("t2", "2", "Partnered 3′ transduction", IN,
    """<p>Transcription of a source L1 reads through its weak polyadenylation signal and stops at a downstream
    site, so the RNA carries unique flanking sequence after the L1 3′ end. The new copy therefore reads
    <i>L1 → (source poly-A, often short) → source 3′ flank → poly-A</i>. The transduced segment is a barcode
    for the source element. Zumalave et al. showed that the strength of the canonical vs. downstream
    polyadenylation sites of a source decides whether it makes mostly solo-L1s or mostly transductions.</p>""",
    [figure("f2", [fig_td_partnered()], "Partnered transduction: the 3′ clip contains the poly-A followed by unique sequence whose mates map to the source locus.")],
    pts("discovery: the 3′ CLIP/POLYA reads carry poly-A then non-L1 sequence; DISC mates on the 3′ side map with high MAPQ to another chromosome (translocation-like)",
        "combine: indel-aware consensus must keep the sequence <i>beyond</i> the poly-A (<code>beyond_polya</code>)",
        "annotate: element <code>L1</code>, tag <code>TD3P</code> + <code>TD3P_SOURCE=&lt;id&gt;</code> from <code>flanks_3p.fa</code>, or <code>NOVEL_SOURCE</code> under the novel-source rule",
        "score: covered 3′-beyond-poly-A sequence (extra points with ≥ 2 independent fragments); the 3′ tag being the transduced flank of exactly the element at the 5′ end counts as concordant"),
    "Zumalave 2026: 1,535 partnered (24% of 6,418). Tubio 2014: transductions 24% of somatic L1s (655/2,756), about half partnered. Nam 2023 normal colon: only 1% (11/1,198). PCAWG: 3,669 transductions (partnered + orphan, 19%). Transduced segments reach 12 kb but are mostly short.",
    "element=L1; tags=TD3P,TD3P_SOURCE=…",
    "3′ clip = A(n) + unique source flank; 3′ mates to the source locus")

sec("t3", "3", "Orphan 3′ transduction", IN,
    """<p>The same read-through transcript, but reverse transcription is so 5′-truncated that no L1 sequence
    is copied: only unique downstream sequence of a source plus a poly-A lands at the new site. With
    short reads it looks like a small unique-sequence insertion or a fake translocation to the source
    locus; the poly-A, TSD and EN motif are what identify it as retrotransposition.</p>""",
    [figure("f3", [fig_td_orphan()], "Orphan transduction: no L1 at all; both junction clips and mates are unique sequence from the source flank.")],
    pts("discovery: CLIP reads whose clips map uniquely elsewhere (not to a repeat) — must not be dropped by an element filter; POLYA reads on the 3′ side",
        "annotate: element <code>ORPHAN_TD</code>, <code>TD3P_SOURCE=&lt;id&gt;</code>; the segment must lie within 15 kb downstream (strand-aware) of a source L1, otherwise it is not an orphan transduction",
        "score: poly-A + TSD + EN motif carry the call because there is no element sequence; concordance = 5′ clip and 3′ tag map to the same source flank"),
    "Zumalave 2026: 705 (11% of 6,418). Tubio 2014: 333/655 transductions (≈12% of all somatic L1). Nam 2023 normal colon: 124/1,198 (10%) — orphans far outnumber partnered events there.",
    "element=ORPHAN_TD; tags=TD3P_SOURCE=…",
    "unique-sequence clips + poly-A; mates to a source flank")

sec("t4", "4", "Alu and SVA insertions; SVA 5′ transduction", IN,
    """<p>Alu (~0.3 kb) and SVA (a few kb; hexamer, Alu-like, VNTR and SINE-R modules) are non-autonomous:
    they are mobilised in trans by L1 ORF2p, so they carry the same TPRT hallmarks (TSD, poly-A, EN motif).
    An Alu is shorter than a sequencing fragment, so read pairs often span the whole insert. SVAs driven by
    an upstream cellular promoter can carry <b>5′ transductions</b> — upstream flank (frequently an exon of
    the host gene, e.g. the MAST2-driven SVA_F group) placed in front of the SVA.</p>""",
    [figure("f4a", [fig_alu()], "Alu: both junctions in one fragment; spanning pairs and reads entirely inside the Alu."),
     figure("f4b", [fig_sva()], "SVA with 5′ transduction: the 5′ clip is unique upstream sequence of the source SVA, not SVA.")],
    pts("discovery: as for L1; SPAN pairs (D2 flag) help very short inserts",
        "annotate: element <code>ALU</code> / <code>SVA</code>; tag <code>TD5P</code> when the 5′ clip maps to <code>flanks_5p_sva.fa</code>",
        "score: identity to young AluY / SVA_E-F consensus; old subfamily only scores negative"),
    "Zumalave 2026: Alu 119 (2%), SVA 6 (&lt;0.1%). PCAWG: Alu 130, SVA 23 (together &lt;1%). SVA 5′ transductions: ~8% of all reference SVAs (Damert 2009, germline); no somatic rate available.",
    "element=ALU | SVA; tags=TD5P",
    "Alu/SVA clips + A(n); short insert spanned by pairs")

sec("t5", "5", "Processed pseudogene", IN,
    """<p>L1 ORF2p occasionally reverse-transcribes a spliced cellular mRNA. The new copy has a TSD and a
    poly-A like an L1 insertion, but its body is exons joined without introns. Discordant mates map to the
    exons of the source gene, which alone would look like many small deletions at the source; the
    <b>exon–exon junction</b> read (a read aligned in exon <i>n</i> whose clip is exon <i>n+1</i>) is the proof
    that the inserted sequence is processed mRNA.</p>""",
    [figure("f5", [fig_pseudogene()], "Processed pseudogene: reads across exon boundaries clip exactly at splice sites; insertion junctions carry TSD and poly-A.")],
    pts("discovery: CLIP at both insertion junctions with clips that map to a gene; poly-A reads",
        "annotate: element <code>PSEUDOGENE</code> + tag <code>EXON_JUNCTION</code> when the consensus/reads.fa contain a junction joining two annotated exons at their splice sites (exon track on hs1)",
        "score: exon–exon junction is a positive feature; without it the call stays <code>UNCERTAIN</code>"),
    "Zumalave 2026: 133 (2%), up to 31 in one tumour. PCAWG: 274 (1.4%), 70 of them in a single pancreatic tumour. Germline (Ewing 2013): 48 gene retrocopy insertions in 1000 Genomes data, exon–exon junctions recovered for 39.",
    "element=PSEUDOGENE; tags=EXON_JUNCTION",
    "clips map to a gene; exon–exon junction reads; A(n) + TSD")

sec("t6", "6", "Solitary poly(A/T) insertion", IN,
    """<p>An extreme 5′ truncation: reverse transcription copies only the RNA poly(A) tail before the
    second strand is resolved, leaving a TSD-flanked poly(A) (or poly(T) on the minus strand) with no
    element sequence. It is the hardest true type because a poly-A clip is also the commonest artefact
    signature (see <a href="#a2">poly-A slippage</a>).</p>""",
    [figure("f6", [fig_polya_only()], "Solitary poly(A/T): both junction clips are A(n) and the TSD is the only other hallmark.")],
    pts("discovery: CLIP reads with A(n) clips on both sides; the single-read poly-A rescue passes single fragments through",
        "annotate: element <code>POLYA_ONLY</code>",
        "score: requires a real TSD and EN motif plus ≥ 2 independent fragments at each end; negative if the site is an A-rich reference tract (slippage context)"),
    "Zumalave 2026: 157 (2% of 6,418). Not reported separately in PCAWG or Nam 2023.",
    "element=POLYA_ONLY",
    "A(n) clips at both junctions, TSD, no element")

sec("t7", "7", "L1-mediated deletion", IN,
    """<p>The L1 cDNA, primed at the original target site, pairs with a 3′ overhang of a distant pre-existing
    break upstream instead of the second nick; aberrant repair removes everything in between. The insertion
    keeps a poly-A and EN motif at the 3′ end but has <b>no TSD</b>; the 5′ end joins via 1–5 bp
    microhomology. Short reads show a one-sided cluster (no reciprocal cluster within ~500 bp) plus a
    copy-number loss.</p>""",
    [figure("f7", fig_l1_del(), "L1-mediated deletion: the two junction clusters are far apart on the reference and copy number drops between them.")],
    pts("discovery: two one-sided junctions separated by the deleted segment; may be emitted as two loci",
        "annotate: tag <code>L1_MED_DELETION</code> (with element <code>L1</code> or <code>ORPHAN_TD</code>); 5′ microhomology instead of TSD",
        "score: EN motif + poly-A at the 3′ junction; no TSD bonus; tag carried, not treated as an artefact"),
    "PCAWG: 90 events, 0.5 kb – 53.4 Mb; 5′ microhomology (median 3 bp) in 75% (47/63) resolved junctions; EN motif at the 3′ end in 82% (74/90); 8% (7/90) EN-independent. Zumalave 2026: deletion-like junctions are the largest rearrangement class, 43% (66/152).",
    "tags=L1_MED_DELETION",
    "one-sided clusters, Mb apart, copy-number loss, no TSD")

sec("t8", "8", "L1-mediated duplication, fold-back inversion and BFB initiation", IN,
    """<p>When the L1 cDNA engages a break on the sister chromatid or upstream of the target, the intervening
    segment is duplicated rather than lost. In the fold-back form an L1 bridges the two sister chromatids in
    head-to-head orientation: two clusters with the same orientation sit a few kb apart, copy number rises
    on one side, and the dicentric product can start breakage–fusion–bridge (BFB) cycles that amplify
    oncogenes (CCND1 in two PCAWG tumours).</p>""",
    [figure("f8", fig_foldback(), "Fold-back inversion bridged by an L1: both flanks of the bridge are the same chromosome arm in opposite orientation.")],
    pts("discovery: clusters whose mates hit L1 on both sides with the same orientation",
        "annotate: tag <code>L1_MED_DUPLICATION</code> (copy-number context is outside PEAR-TREE; the tag records the junction geometry)",
        "score: TPRT hallmarks at the poly-A junction score as usual"),
    "PCAWG: individual cases (79.6 Mb 14q duplication; CCND1 BFB amplification in 2 tumours) among 103 L1-mediated rearrangements. Zumalave 2026: duplication-like 9% (14/152) of rearrangements.",
    "tags=L1_MED_DUPLICATION",
    "same-orientation clusters into L1, copy-number gain")

sec("t9", "9", "Retrotransposon-mediated rearrangements (deletion-, duplication-, inversion-like)", PART,
    """<p>Zumalave et al. define an RT-RG as two distant breakpoints joined by a somatic retrotransposition
    bridge and classify intrachromosomal ones by breakpoint orientation: deletion-like [+/−],
    duplication-like [−/+], inversion-like [+/+] or [−/−]; interchromosomal ones are translocation-like.
    The bridge is most often a solo-L1 but can be any type, including transductions and pseudogenes.</p>""",
    [figure("f9", fig_rtrg(), "Three intrachromosomal RT-RG geometries; only the far flank (B) differs.")],
    pts("deletion-like → <code>L1_MED_DELETION</code>, duplication-like → <code>L1_MED_DUPLICATION</code> (in scope as tags)",
        "inversion-like and translocation-like bridges are out of scope this round; they should still come out as single junctions with TPRT hallmarks rather than be scored as artefacts"),
    "Zumalave 2026: 152 RT-RGs = 2% of events; deletion-like 43%, interchromosomal 30%, inversion-like 18%, duplication-like 9%; bridges are solo-L1 in 66%. PCAWG: 103 L1-mediated rearrangements, mostly deletions.",
    "tags=L1_MED_DELETION | L1_MED_DUPLICATION (inversion-/translocation-like: not called)",
    "L1 bridge between two distant breakpoints")

sec("t10", "10", "Reciprocal translocation bridged by two concurrent L1s", OUT_,
    """<p>Two synchronous L1 insertions on different chromosomes exchange partners during second-strand
    resolution, so each derivative chromosome carries one complete L1 insertion (own poly-A, own EN motif)
    and each insertion's two TSD copies end up on different derivatives — a <b>target-site
    heteroduplication</b> (e.g. 3 bp and −10 bp in one case). Short reads see four one-sided junctions that
    look like two translocation partners each bridged by an L1.</p>""",
    [figure("f10", fig_recip_tl(), "Reciprocal translocation: copies a/a′ and b/b′ of each TSD sit on different derivative chromosomes.")],
    pts("out of scope this round: no pairing of junctions across chromosomes",
        "each junction may still be reported as a one-sided L1 junction; it must not be penalised as <code>CHIMERIC_ENDS</code> just because the two ends map to different chromosomes"),
    "Zumalave 2026: 13 reciprocal translocations (20 of 45 interchromosomal bridges); 10/13 from two distinct synchronous insertions, 3/13 from a single event.",
    "— (not called)",
    "four one-sided junctions on two chromosomes; heteroduplicated TSD")

sec("t11", "11", "Chimeric bridges", OUT_,
    """<p>In some reciprocal translocations the two bridges are mosaics of two retrotransposition events
    (e.g. two different partnered transductions) that recombined while the cDNAs were growing. One bridge
    has both TPRT ends (a poly-A, a poly-T and two EN motifs); its partner contains the leftover L1
    pieces with no poly-A and no EN motif. The exchange points coincide with orientation flips, as in twin
    priming.</p>""",
    [figure("f11", fig_chimeric_bridge(), "Chimeric bridges: the poly-A-carrying ends of two insertions end up in one bridge, remnants in the other.")],
    pts("out of scope this round",
        "a hallmark-free bridge looks like an artefact chimera (<code>CHIMERIC_ENDS</code>); that is acceptable while bridges are out of scope"),
    "Zumalave 2026: a subset of the 13 reciprocal translocations (example PD0331a, 22p/12q); no separate rate.",
    "— (not called)",
    "mosaic L1/TD bridges, one without hallmarks")

sec("t12", "12", "Reciprocal inversions and complex events", OUT_,
    """<p>The same in-trans second-strand synthesis between two concurrent L1 reactions (one canonical, one
    twin-primed) can produce paracentric inversions where each breakpoint is bridged by part of the same
    L1, or chains of rearrangements with acentric/circular products and large deletions when a single
    twin-primed cDNA repairs two pre-existing breaks.</p>""",
    [figure("f12", fig_recip_inv(), "Reciprocal inversion (PD0307a-like): the two bridges are halves of related L1 insertions, megabases apart.")],
    pts("out of scope this round"),
    "Zumalave 2026: rare — e.g. a 2.67 Mb paracentric inversion (PD0307a) and a chromosome-1 complex event with a 30 Mb loss (PD0331a).",
    "— (not called)",
    "L1 fragments at Mb-distant breakpoints with shared inversion breakpoints")

sec("t13", "13", "Twin priming with 5′ or 3′ switching", IN,
    """<p>Variants of twin priming in which one of the two growing cDNAs changes template on the L1 RNA.
    In <b>5′ switching</b> the internally primed cDNA jumps to a region of the RNA not yet copied; in
    <b>3′ switching</b> the poly-dT-primed cDNA jumps to a region downstream of the internal primer. Both
    leave an extra internal inversion or a coordinate discontinuity at the twin-priming junction, on top of
    the usual 5′ inversion.</p>""",
    [figure("f13", fig_switch(), "Twin priming with switching: more than one orientation change or a jump in L1 coordinates within the 5′ part.")],
    pts("requires the covered-element consensus from <code>insertions.reads.fa</code>; usually only resolved when the insert is short enough for reads to cross the extra junctions",
        "annotate: <code>structure=INVERTED_5P_SWITCH</code>; falls back to <code>INVERTED_5P</code> / <code>5P_UNRESOLVED</code> when the extra junction is not covered"),
    "Zumalave 2026: 20 with 5′ switching and 3 with 3′ switching, out of 1,836 5′-inverted solo-L1s (~1%). Szak 2002 had already noted twice-inverted genomic L1s.",
    "structure=INVERTED_5P_SWITCH",
    "≥ 2 orientation changes / coordinate jump inside the L1 body")

sec("t14", "14", "Templated local insertions (local chimeric L1)", IN,
    """<p>A non-L1 segment is copied into the L1 during integration, giving a chimeric L1 that duplicates a
    small genomic segment. Zumalave et al. report templates typically shorter than 250 bp lying within 15 bp
    of the insertion site, and favour a model in which the growing cDNA leaves the L1 RNA and anneals to
    nearby homologous DNA (akin to synthesis-dependent strand annealing). In reads the extra segment is a clip
    that maps back to the immediate neighbourhood, followed by L1. Where in the insert the segment sits varies;
    the figure draws it at the 5′ end.</p>""",
    [figure("f14", [fig_templated()], "Templated local insertion: the 5′ clip first matches nearby reference sequence, then continues into L1.")],
    pts("annotate: tag <code>TEMPLATED_LOCAL</code> when a non-element segment of the insert (typically < 250 bp) maps within 15 bp of the site",
        "must not trigger the fold-back artefact penalty: the clip continues into element sequence and the 3′ end carries poly-A + TSD"),
    "Zumalave 2026: 22 events among ~3,600 solo-L1s (&lt;1%).",
    "tags=TEMPLATED_LOCAL",
    "5′ clip = nearby reference segment then L1")

sec("t15", "15", "Co-inserted local pre-mRNA", IN,
    """<p>The reverse transcriptase switches from the L1 RNA to a pre-mRNA transcribed near the insertion
    site and co-inserts part of it, introns included. Unlike a processed pseudogene the co-inserted piece is
    unspliced and comes from a gene next to the site.</p>""",
    [figure("f15", [fig_premrna()], "Co-inserted pre-mRNA: the 5′ part of the insert maps to an exon and intron of a nearby gene.")],
    pts("annotate: tag <code>PREMRNA_COINSERT</code> when the 5′ part of the insert maps to a gene within the neighbourhood and spans an exon–intron boundary (unspliced); distinct from <code>EXON_JUNCTION</code>"),
    "Nam 2023: a single event (in a clone from an adenoma) among ~1,200 resolved insertions.",
    "tags=PREMRNA_COINSERT",
    "5′ clip maps to an unspliced nearby transcript")

sec("t16", "16", "Fold-back inverted duplication 5′ of the target site", IN,
    """<p>After integration, extra DNA synthesis at the 5′ end copies a short stretch of upstream flank
    back on itself, leaving an inverted duplication of 52–220 bp between the 5′ TSD and the L1. Its 5′ clip
    is therefore the reverse complement of the adjacent reference — the same as a fold-back artefact —
    but it is followed by L1 and the 3′ end carries a normal poly-A and TSD.</p>""",
    [figure("f16", [fig_foldback_5p()], "Fold-back inverted duplication: reverse-complement local clip at the 5′ junction, then L1.")],
    pts("annotate: tag <code>FOLDBACK_INVDUP_5P</code>; the score's “clip identical to the adjacent reference” penalty applies only when the clip does not continue into element sequence and there is no TPRT 3′ end"),
    "Nam 2023: 3 events (0.3%) plus 1 combined with a 5′ inversion (0.1%) among ~1,200.",
    "tags=FOLDBACK_INVDUP_5P",
    "5′ clip = reverse-complement adjacent reference, then L1")

sec("tA", "A", "Target-site deletion instead of a TSD", IN,
    """<p>If the second nick lies upstream of the first, or the overhang is trimmed, the insertion is
    flanked by a small deletion of the target rather than a duplication. The 5′ and 3′ junction positions
    are then separated by the deleted bases instead of overlapping.</p>""",
    [figure("fA", fig_tsd_del(), "Target-site deletion: junctions do not overlap; the few deleted bases are absent from the insertion allele.")],
    pts("annotate: tag <code>TSD_DELETION</code>, <code>tsd_len</code> negative",
        "score: weak positive for deletions ≤ 20 bp (strong positive is reserved for TSD 4–25 bp with identity)"),
    "Common enough that Nam 2023 plot negative target-site lengths alongside TSDs and Zumalave call it fairly frequent; neither source gives a fraction in the text.",
    "tags=TSD_DELETION",
    "junctions do not overlap; no duplicated bases")

sec("tB", "B", "Endonuclease-independent insertions", IN,
    """<p>L1 can use a pre-existing break (or dysfunctional telomere) instead of its endonuclease nick. The
    hallmarks disappear: no TSD, no EN motif, and the element is often truncated at both ends (no poly-A).</p>""",
    [figure("fB", fig_en_indep(), "EN-independent insertion: L1 sequence with neither TSD nor poly-A.")],
    pts("annotate: tag <code>EN_INDEPENDENT</code> when there is no TSD, no poly-A and both ends are truncated",
        "score: few TPRT points by construction; such calls stay <code>UNCERTAIN</code> unless independent-fragment and colony evidence is strong"),
    "PCAWG: 7/90 (8%) of L1-mediated deletions. Somatic rate for plain insertions is unknown.",
    "tags=EN_INDEPENDENT",
    "L1 clips, no TSD, no poly-A")

TOC_TITLE = {
    "t1": "Solo-L1", "t2": "Partnered 3′ transduction", "t3": "Orphan 3′ transduction", "t4": "Alu, SVA, SVA 5′ TD",
    "t5": "Processed pseudogene", "t6": "Solitary poly(A/T)", "t7": "L1-mediated deletion",
    "t8": "L1-mediated duplication / fold-back", "t9": "RT-mediated rearrangements",
    "t10": "Reciprocal translocation (2 L1s)", "t11": "Chimeric bridges", "t12": "Reciprocal inversion / complex",
    "t13": "Twin priming + switching", "t14": "Templated local insertion", "t15": "Co-inserted pre-mRNA",
    "t16": "Fold-back inv-dup 5′ of site", "tA": "Target-site deletion", "tB": "EN-independent insertion",
}

SHORT_FREQ = {
    "t1": "56% of events (Z); 78% (PCAWG); 89% normal colon (Nam)",
    "t2": "24% (Z); 1% normal colon (Nam)",
    "t3": "11% (Z); 10% normal colon (Nam)",
    "t4": "Alu 2%, SVA &lt;0.1% (Z)",
    "t5": "2% (Z); 1.4% (PCAWG)",
    "t6": "2% (Z)",
    "t7": "90 in PCAWG; 43% of RT-RGs (Z)",
    "t8": "9% of RT-RGs (Z); case reports (PCAWG)",
    "t9": "2% of events (Z)",
    "t10": "13 events in 10 tumours (Z)",
    "t11": "subset of #10 (Z)",
    "t12": "rare, case level (Z)",
    "t13": "23/1,836 inverted L1s (Z)",
    "t14": "22 events (Z)",
    "t15": "1 event (Nam)",
    "t16": "0.4% normal colon (Nam)",
    "tA": "not quantified in sources",
    "tB": "8% of L1-mediated deletions (PCAWG)",
}

ARTS = []


def art(sid, title, body, fig, detect):
    ARTS.append(dict(sid=sid, title=title, body=body, fig=fig, detect=detect))


art("a1", "Ligation, PCR and MDA chimeras",
    """<p>Library ligation, PCR template switching and (for single cells) whole-genome amplification join
    unrelated molecules. A genomic fragment fused to an L1 looks like an insertion junction. Evrony et al.
    showed that 7 of 13 PCR-validated candidates in an earlier single-neuron study were such chimeras: they had
    no TSD, could start at any L1 coordinate (including inactive or truncated copies), and 5′ and 3′
    junctions from <i>different</i> L1s had been paired as one insertion. A-tailing during library prep can
    add a short poly-A at the join.</p>""",
    figure("fa1", fig_art_chimera(), "Chimera: a single junction, supported by one molecule, with no TSD and no EN motif."),
    pts("independence rule: one fragment (plus its PCR duplicates) never makes ≥ 2 independent fragments",
        "<code>CHIMERIC_ENDS</code> when 5′ and 3′ ends come from incompatible elements",
        "no TSD, no EN motif, no poly-A or only a short A-tail → few points; cross-sample identical fragments score negative"))
art("a2", "Poly-A slippage at reference A-tracts",
    """<p>Polymerase slippage and aligner behaviour at long reference A-tracts (typically the tails of
    reference Alus and L1s) produce reads soft-clipped right where the tract ends, with a clip that is
    simply more A. This mimics the 3′ junction of a real insertion and is the main source of discovery
    false positives in PEAR-TREE.</p>""",
    figure("fa2", fig_art_slippage(), "Poly-A slippage: clipped reads pile up at the end of a reference A-tract; no insertion exists."),
    pts("score: poly-A at both junctions or a reference A-rich context is negative",
        "no element sequence beyond the poly-A, no TSD, recurrence across many unrelated loci"))
art("a3", "Mismapping among same-subfamily elements",
    """<p>Reads from a polymorphic or non-reference L1HS (or Alu) copy align to a near-identical reference copy
    and are soft-clipped where the two copies' flanks diverge. The clip is the flank of the other copy, not an
    insertion; a germline polymorphism can thus look somatic at the wrong locus.</p>""",
    figure("fa3", fig_art_mismap(), "Mismapping: clips at the end of a reference L1HS that are the flank of a paralogous copy."),
    pts("clip anchored inside a reference element of the same subfamily; flank clip maps elsewhere",
        "recurrence across colonies/patients and cross-sample identical fragments score negative"))
art("a4", "Fold-back palindromes",
    """<p>Inverted repeats in the reference and hairpin library artefacts produce reads whose clip is the reverse
    complement of the adjacent reference. Without element sequence, poly-A or TSD, this is an artefact; with
    them it may be a real <a href="#t16">fold-back inverted duplication</a>.</p>""",
    figure("fa4", fig_art_foldback(), "Fold-back: clip = reverse complement of the neighbouring reference, nothing else."),
    pts("score: “clip identical to the adjacent reference (fold-back)” is negative unless the clip continues into element sequence and the 3′ end has TPRT hallmarks"))


# ----------------------------------------------------------------------------------------
# Page
# ----------------------------------------------------------------------------------------
LEGEND = [
    ("k-flank", "flank at the site (aligned)"), ("k-flank2", "flank of another locus / chromosome"),
    ("k-tsd", "target-site duplication (TSD)"), ("k-l1", "L1, sense"), ("k-l1inv", "L1, inverted segment"),
    ("k-alu", "Alu"), ("k-sva", "SVA"), ("k-polya", "poly(A/T) tail"), ("k-td", "transduced flank of a source element"),
    ("k-exon", "exon (pseudogene / pre-mRNA)"), ("k-intron", "intron (unspliced)"),
    ("k-local", "local sequence (templated / fold-back copy)"), ("k-gap", "deleted / absent sequence"),
]

REFS = [
    ("Zumalave S, Santamarina M, Espasandín NP, … Rodriguez-Martin B, Tubio JMC. Synchronous L1 retrotransposition events promote chromosomal crossover early in human tumorigenesis. Science (2026), PMID 41747018; preprint bioRxiv", "10.1101/2024.08.27.596794"),
    ("Nam CH, Youk J, Kim JY, … Ju YS. Widespread somatic L1 retrotransposition in normal colorectal epithelium. Nature (2023)", "10.1038/s41586-023-06046-z"),
    ("Rodriguez-Martin B, Alvarez EG, Baez-Ortega A, … Tubio JMC. Pan-cancer analysis of whole genomes identifies driver rearrangements promoted by LINE-1 retrotransposition. Nat Genet (2020)", "10.1038/s41588-019-0562-0"),
    ("Tubio JMC, Li Y, Ju YS, … Campbell PJ. Extensive transduction of nonrepetitive DNA mediated by L1 retrotransposition in cancer genomes. Science 345, 1251343 (2014)", "10.1126/science.1251343"),
    ("Ostertag EM, Kazazian HH Jr. Twin priming: a proposed mechanism for the creation of inversions in L1 retrotransposition. Genome Res 11, 2059–2065 (2001)", "10.1101/gr.205701"),
    ("Flasch DA, Macia Á, Sánchez L, … Moran JV. Genome-wide de novo L1 retrotransposition connects endonuclease activity with replication. Cell (2019)", "10.1016/j.cell.2019.02.050"),
    ("Szak ST, Pickeral OK, Makalowski W, … Boeke JD. Molecular archeology of L1 insertions in the human genome. Genome Biol 3, research0052 (2002)", "10.1186/gb-2002-3-10-research0052"),
    ("Ewing AD, Ballinger TJ, Earl D, … Haussler D. Retrotransposition of gene transcripts leads to structural variation in mammalian genomes. Genome Biol 14, R22 (2013)", "10.1186/gb-2013-14-3-r22"),
    ("Evrony GD, Lee E, Park PJ, Walsh CA. Resolving rates of mutation in the brain using single-neuron genomics. eLife 5, e12966 (2016)", "10.7554/eLife.12966"),
    ("Gardner EJ, Lam VK, Harris DN, … Devine SE. The Mobile Element Locator Tool (MELT): population-scale mobile element discovery and biology. Genome Res (2017)", "10.1101/gr.218032.116"),
    ("Scott EC, Gardner EJ, Masood A, … Devine SE. A hot L1 retrotransposon evades somatic repression and initiates human colorectal cancer. Genome Res (2016)", "10.1101/gr.201814.115"),
    ("Brouha B, Schustak J, Badge RM, … Kazazian HH Jr. Hot L1s account for the bulk of retrotransposition in the human population. PNAS 100, 5280–5285 (2003)", "10.1073/pnas.0831042100"),
    ("Damert A, Raiz J, Horn AV, … Schumann GG. 5′-Transducing SVA retrotransposon groups spread efficiently throughout the human genome. Genome Res 19, 1992–2008 (2009)", "10.1101/gr.093435.109"),
    ("Gilbert N, Lutz-Prigge S, Moran JV. Genomic deletions created upon LINE-1 retrotransposition. Cell 110, 315–325 (2002)", "10.1016/S0092-8674(02)00828-0"),
    ("Morrish TA, Gilbert N, Myers JS, … Moran JV. DNA repair mediated by endonuclease-independent LINE-1 retrotransposition. Nat Genet 31, 159–165 (2002)", "10.1038/ng898"),
]

CSS = r"""
:root{
  --bg:#fbfaf7; --panel:#ffffff; --ink:#1d232b; --muted:#5b6672; --line:#d9dde2; --accent:#2b5fb8;
  --c-flank:#a7b0ba; --c-flank2:#5e6d7f; --c-tsd:#f0a33a; --c-l1:#2f66d0; --c-l1inv:#16a2a2;
  --c-alu:#2f9a58; --c-sva:#c0408f; --c-polya:#d43c3c; --c-td:#8a56e8; --c-exon:#c9971c; --c-exonb:#e3bd55;
  --c-intron:#efdfa8; --c-local:#a65a2c; --c-gap:#9aa3ad; --c-onfill:#ffffff; --c-ondark:#1d232b;
  --b-in:#1f7a46; --b-in-bg:#e2f3e8; --b-part:#8a5a00; --b-part-bg:#fbefd2; --b-out:#8b2b2b; --b-out-bg:#f8e1e1;
  --b-art:#4a4f7a; --b-art-bg:#e6e7f4;
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --bg:#13171c; --panel:#1b2128; --ink:#e4e8ec; --muted:#9aa5b1; --line:#2e3640; --accent:#7aa7f5;
    --c-flank:#6f7a86; --c-flank2:#9fb0c4; --c-tsd:#f3b55c; --c-l1:#5d8ef0; --c-l1inv:#33c4c4;
    --c-alu:#4cc17b; --c-sva:#e06ab2; --c-polya:#f06464; --c-td:#a983ff; --c-exon:#dcae3a; --c-exonb:#efd27e;
    --c-intron:#6e6440; --c-local:#d07d4b; --c-gap:#6f7a86; --c-onfill:#0f1317; --c-ondark:#0f1317;
    --b-in:#7fd8a2; --b-in-bg:#173322; --b-part:#f2c56b; --b-part-bg:#3a2d10; --b-out:#f29a9a; --b-out-bg:#3b1a1a;
    --b-art:#b9bdf0; --b-art-bg:#262a4a;
  }
}
:root[data-theme="dark"]{
    --bg:#13171c; --panel:#1b2128; --ink:#e4e8ec; --muted:#9aa5b1; --line:#2e3640; --accent:#7aa7f5;
    --c-flank:#6f7a86; --c-flank2:#9fb0c4; --c-tsd:#f3b55c; --c-l1:#5d8ef0; --c-l1inv:#33c4c4;
    --c-alu:#4cc17b; --c-sva:#e06ab2; --c-polya:#f06464; --c-td:#a983ff; --c-exon:#dcae3a; --c-exonb:#efd27e;
    --c-intron:#6e6440; --c-local:#d07d4b; --c-gap:#6f7a86; --c-onfill:#0f1317; --c-ondark:#0f1317;
    --b-in:#7fd8a2; --b-in-bg:#173322; --b-part:#f2c56b; --b-part-bg:#3a2d10; --b-out:#f29a9a; --b-out-bg:#3b1a1a;
    --b-art:#b9bdf0; --b-art-bg:#262a4a;
}
*{box-sizing:border-box}
html{scroll-behavior:smooth}
body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.55 "Source Sans 3",system-ui,-apple-system,"Segoe UI",sans-serif}
code{overflow-wrap:anywhere;font-family:"JetBrains Mono",ui-monospace,Menlo,monospace;font-size:.86em;background:color-mix(in srgb,var(--line) 45%,transparent);padding:.05em .3em;border-radius:4px}
a{color:var(--accent)}
.wrap{display:grid;grid-template-columns:minmax(0,1fr);max-width:1240px;margin:0 auto;padding:0 16px}
nav.toc{border-bottom:1px solid var(--line);padding:12px 0;font-size:.9rem}
nav.toc ul{margin:0;padding:0;list-style:none;columns:2;column-gap:1.5em}
nav.toc li{break-inside:avoid;margin:2px 0;padding-left:1.8em;text-indent:-1.8em}
nav.toc .tn{display:inline-block;min-width:1.8em;text-indent:0;color:var(--muted);font-variant-numeric:tabular-nums}
@media (max-width:520px){nav.toc ul{columns:1}}
nav.toc a{text-decoration:none;color:var(--ink)}
nav.toc a:hover{color:var(--accent)}
nav.toc .h{font-weight:700;margin:0 0 6px;color:var(--muted);text-transform:uppercase;letter-spacing:.06em;font-size:.75rem}
main{min-width:0;padding-bottom:64px}
@media (min-width:1020px){
  .wrap{grid-template-columns:250px minmax(0,1fr);gap:36px}
  nav.toc{position:sticky;top:0;align-self:start;max-height:100vh;overflow:auto;border-bottom:0;border-right:1px solid var(--line);padding:24px 12px 24px 0}
  nav.toc ul{columns:1}
}
header h1{font-size:2rem;line-height:1.2;margin:28px 0 6px}
header p.sub{color:var(--muted);margin:0 0 16px}
h2{font-size:1.4rem;margin:48px 0 8px;padding-top:8px;border-top:1px solid var(--line);display:flex;flex-wrap:wrap;align-items:baseline;gap:10px}
h2 .n{color:var(--muted);font-variant-numeric:tabular-nums;min-width:1.6em}
h3{font-size:1rem;margin:18px 0 6px;color:var(--muted);text-transform:uppercase;letter-spacing:.05em}
.badge{font-size:.72rem;font-weight:700;padding:.15em .6em;border-radius:999px;letter-spacing:.03em;white-space:nowrap}
.badge.in{color:var(--b-in);background:var(--b-in-bg)}
.badge.partial{color:var(--b-part);background:var(--b-part-bg)}
.badge.out{color:var(--b-out);background:var(--b-out-bg)}
.badge.art{color:var(--b-art);background:var(--b-art-bg)}
figure{margin:14px 0;background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:10px 12px}
figure{overflow-x:auto}
figure svg{width:100%;min-width:560px;height:auto;display:block}
figcaption{font-size:.86rem;color:var(--muted);margin-top:6px}
.freq{background:var(--panel);border-left:3px solid var(--c-tsd);padding:8px 12px;border-radius:0 8px 8px 0;font-size:.93rem}
.vocab{font-size:.9rem}
ul.pts{margin:4px 0;padding-left:1.2em}
ul.pts li{margin:3px 0}
.legend{display:grid;grid-template-columns:repeat(auto-fill,minmax(220px,1fr));gap:6px 18px;background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:12px 14px;font-size:.88rem}
.legend span.sw{display:inline-block;width:26px;height:12px;border-radius:2px;vertical-align:-1px;margin-right:8px}
.sw.k-gap{background:transparent!important;border:1.5px dashed var(--c-gap)}
.legend .rd{grid-column:1/-1;color:var(--muted);border-top:1px solid var(--line);padding-top:8px;margin-top:4px}
.tbl{overflow-x:auto;border:1px solid var(--line);border-radius:10px;background:var(--panel)}
table{border-collapse:collapse;width:100%;font-size:.86rem;min-width:760px}
th,td{text-align:left;vertical-align:top;padding:7px 10px;border-bottom:1px solid var(--line)}
th{position:sticky;top:0;background:var(--panel);font-weight:700}
ol.refs li{margin:6px 0;font-size:.9rem}
.note{font-size:.9rem;color:var(--muted)}
.cols{display:grid;grid-template-columns:1fr;gap:12px}
@media (min-width:760px){.cols{grid-template-columns:1fr 1fr}}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:10px 14px;font-size:.92rem}
/* SVG classes */
svg text{font-family:"Source Sans 3",system-ui,sans-serif}
.lt{font-size:13px;font-weight:700;fill:var(--ink)}
.lt2{font-size:11px;fill:var(--muted);font-style:italic}
.mk{font-size:11px;fill:var(--muted)}
.mkl{stroke:var(--muted);stroke-width:.7;stroke-dasharray:2 2}
.bl{font-size:12px;fill:var(--ink)}
.it{font-size:10px;fill:var(--c-onfill);font-weight:700}
.itd{font-size:10px;fill:var(--c-ondark);font-weight:700}
.sep{stroke:var(--line);stroke-width:1}
.k-flank,.r-al,.sw.k-flank{fill:var(--c-flank);background:var(--c-flank)}
.k-flank2,.sw.k-flank2{fill:var(--c-flank2);background:var(--c-flank2)}
.k-tsd,.r-tsd,.sw.k-tsd{fill:var(--c-tsd);background:var(--c-tsd)}
.k-l1,.r-l1,.sw.k-l1{fill:var(--c-l1);background:var(--c-l1)}
.k-l1inv,.r-l1inv,.sw.k-l1inv{fill:var(--c-l1inv);background:var(--c-l1inv)}
.k-alu,.r-alu,.sw.k-alu{fill:var(--c-alu);background:var(--c-alu)}
.k-sva,.r-sva,.sw.k-sva{fill:var(--c-sva);background:var(--c-sva)}
.k-polya,.r-polya,.sw.k-polya{fill:var(--c-polya);background:var(--c-polya)}
.k-td,.r-td,.sw.k-td{fill:var(--c-td);background:var(--c-td)}
.k-exon,.r-exon,.sw.k-exon{fill:var(--c-exon);background:var(--c-exon)}
.k-exonb,.r-exonb{fill:var(--c-exonb)}
.k-intron,.r-intron,.sw.k-intron{fill:var(--c-intron);background:var(--c-intron)}
.k-local,.r-local,.sw.k-local{fill:var(--c-local);background:var(--c-local)}
.gapr{fill:none;stroke:var(--c-gap);stroke-width:1.4;stroke-dasharray:4 3}
.gapl{stroke:var(--c-gap);stroke-width:1}
.r-clip2{fill:var(--c-flank2)}
.rh{fill:var(--muted)}
.pl{stroke:var(--muted);stroke-width:1}
.pld{stroke:var(--muted);stroke-width:1;stroke-dasharray:3 2}
"""


def legend_svg():
    """Small key figure for read conventions."""
    ln = Lane("how reads are drawn (insertion-allele coordinates)", [
        ("flank", 30, "flank"), ("tsd", 4), ("l1", 28, "element", "+"), ("polya", 5), ("tsd", 4), ("flank", 30, "flank")])
    ln.split(2, 0.5)
    ln.pair(ln.bx(2) - L - 70, ln.bx(2) + 60)
    ln.split(4, 0.6, mate_dx=FR)
    return figure("fkey", [ln], "Grey read segments align at this locus; coloured segments are soft-clipped or (for a whole read) a mate that maps elsewhere, coloured by what they contain. Solid line = concordant pair, dashed = discordant pair; arrowheads give read strand.")


def build():
    toc = ["<li><a href='#intro'><span class='tn'></span>How to read this page</a></li>",
           "<li><a href='#vocab'><span class='tn'></span>Vocabulary &amp; evidence</a></li>",
           "<li><a href='#structure'><span class='tn'></span>5′ structure classes</a></li>"]
    for s in SECTIONS:
        toc.append(f"<li><a href='#{s['sid']}'><span class='tn'>{s['num']}</span>{esc(TOC_TITLE[s['sid']])}</a></li>")
    toc.append("<li><a href='#artefacts'><span class='tn'>✕</span>Artefacts that mimic MEIs</a></li>")
    toc.append("<li><a href='#summary'><span class='tn'></span>Summary table</a></li><li><a href='#refs'><span class='tn'></span>References</a></li>")

    legend = "".join(f"<div><span class='sw {c}'></span>{esc(t)}</div>" for c, t in LEGEND)

    body = []
    body.append(f"""
<header>
<h1>Retrotransposition insertion types</h1>
<p class="sub">Shared reference for the PEAR-TREE caller, the simulators and annotation — what each insertion
type looks like in the genome and in short reads, how often it occurs, and how PEAR-TREE names and scores it.
Branch <code>tprt-hallmarks</code>; vocabulary per <code>plans/tprt_hallmarks/SPEC.md</code>.</p>
</header>
<section id="intro">
<h2>How to read this page</h2>
<p>Every figure has the same two lanes. The <b>top lane</b> is the insertion allele (not to scale): grey
flanks, orange TSD boxes, the element as an arrow pointing 5′→3′, the poly(A) tail in red, transduced
source flank in violet, inverted segments in teal, deletions as dashed boxes. The <b>bottom lane</b>
shows the short-read evidence drawn on the same coordinates.</p>
<div class="legend">{legend}
<div class="rd">One legend for every figure. Colours are CSS custom properties and switch with the dark theme.</div></div>
{legend_svg()}
<p class="note">Frequencies come from different settings and are not interchangeable: Zumalave 2026 selected 10
tumours with &gt; 100 somatic events each (6,418 events, long + short reads); PCAWG (Rodriguez-Martin 2020) is
2,954 cancer genomes with short reads (19,166 events); Nam 2023 is ~1,200 resolved events in clonally expanded
normal colorectal cells — the setting closest to PEAR-TREE's colony data; Szak 2002 and Ewing 2013 describe
germline/reference copies.</p>
</section>

<section id="vocab">
<h2>Vocabulary and evidence</h2>
<div class="cols">
<div class="card"><b>TPRT hallmarks</b>
<ul class="pts"><li><b>TSD</b> — target-site duplication, typically 4–25 bp (in vivo average ~14 bp), A-rich.</li>
<li><b>EN motif</b> — the L1 endonuclease nick, consensus 5′-TTTT/AA-3′ (a preference, not a fixed sequence: even the most used single 7-mer, TTTTT/AA, accounts for &lt; 10% of insertions).</li>
<li><b>poly(A)</b> at the 3′ end, on the strand consistent with element orientation (poly(T) when the element is on the minus strand).</li>
<li><b>5′ truncation / inversion</b> as products of early termination and twin priming.</li></ul></div>
<div class="card"><b>PEAR-TREE output (SPEC)</b>
<ul class="pts"><li><code>element</code>: L1, ALU, SVA, PSEUDOGENE, POLYA_ONLY, ORPHAN_TD, NON_TPRT, UNKNOWN</li>
<li><code>structure</code>: FULL_LENGTH, TRUNCATED_5P, INVERTED_5P, INVERTED_5P_SWITCH, 5P_UNRESOLVED</li>
<li><code>tags</code>: TD3P, TD5P, TD3P_SOURCE=…, NOVEL_SOURCE, TEMPLATED_LOCAL, PREMRNA_COINSERT, TSD_DELETION, EN_INDEPENDENT, L1_MED_DELETION, L1_MED_DUPLICATION, FOLDBACK_INVDUP_5P, CHIMERIC_ENDS, EXON_JUNCTION</li>
<li>evidence roles: CLIP, POLYA, MATE, DISC, SPAN; a junction is supported with ≥ 2 independent fragments pooled over all colonies of a patient — the poly-A end included; independence is decided from the reads, never from the duplicate flag.</li>
<li><code>tprt_call</code>: TPRT / LIKELY_TPRT / UNCERTAIN / ARTEFACT_LIKE from the additive <code>tprt_score</code>.</li></ul></div>
</div>
</section>

<section id="structure">
<h2>5′ structure classes</h2>
<p>The 3′ end (poly-A junction) is shared by all TPRT insertions; what differs is the 5′ end. <code>structure</code>
is only assigned when the 5′ junction is covered (clip consensus or mates), otherwise <code>5P_UNRESOLVED</code>.</p>
<div class="tbl"><table><thead><tr><th>structure</th><th>5′ junction evidence</th><th>mechanism</th><th>frequency (Zumalave 2026 solo-L1)</th></tr></thead><tbody>
<tr><td><code>FULL_LENGTH</code></td><td>5′ clip reaches the start of the L1 5′ UTR, same strand as the 3′ end</td><td>complete reverse transcription</td><td>0.3% (10/3,611)</td></tr>
<tr><td><code>TRUNCATED_5P</code></td><td>5′ clip starts at an internal L1 coordinate, same strand</td><td>early termination of TPRT</td><td>49% of truncated (1,762/3,598)</td></tr>
<tr><td><code>INVERTED_5P</code></td><td>5′ clip/mates are L1 on the opposite strand; inversion junction inside the insert</td><td>twin priming (± internal del/dup)</td><td>51% of truncated (1,836/3,598)</td></tr>
<tr><td><code>INVERTED_5P_SWITCH</code></td><td>extra orientation change or coordinate jump inside the 5′ part</td><td>twin priming with 5′/3′ switching</td><td>23 of 1,836 inverted</td></tr>
<tr><td><code>5P_UNRESOLVED</code></td><td>no 5′ junction coverage</td><td>—</td><td>—</td></tr>
</tbody></table></div>
<p class="note">Figures for the first four classes are in <a href="#t1">type 1</a> and <a href="#t13">type 13</a>.</p>
</section>
""")

    for s in SECTIONS:
        figs = "".join(s["fig"])
        body.append(f"""
<section id="{s['sid']}">
<h2><span class="n">{s['num']}</span>{esc(s['title'])} <span class="badge {s['scope']}">{BADGE[s['scope']]}</span></h2>
{s['body']}
{figs}
<h3>Reported frequency</h3><p class="freq">{s['freq']}</p>
<h3>How PEAR-TREE detects it</h3>{s['detect']}
<p class="vocab"><b>Vocabulary:</b> <code>{esc(s['vocab'])}</code></p>
</section>""")

    body.append("""<section id="artefacts"><h2><span class="n">✕</span>Artefacts that mimic MEIs <span class="badge art">score negative</span></h2>
<p>The TPRT point system exists to separate these from the types above. None of them produces the full
combination of a TSD with identity, an EN motif at the nick, a strand-consistent poly-A, concordant 5′/3′
element classes and ≥ 2 independent fragments.</p>""")
    for a in ARTS:
        body.append(f"""<section id="{a['sid']}"><h3 style="text-transform:none;letter-spacing:0;color:var(--ink);font-size:1.1rem">{esc(a['title'])}</h3>
{a['body']}{a['fig']}<p><b>Separating it:</b></p>{a['detect']}</section>""")
    body.append("</section>")

    rows = []
    for s in SECTIONS:
        rows.append(f"<tr><td><a href='#{s['sid']}'>{s['num']}</a></td><td>{esc(s['title'])}</td><td>{esc(s['signature'])}</td>"
                    f"<td><code>{esc(s['vocab'])}</code></td><td>{SHORT_FREQ[s['sid']]}</td>"
                    f"<td><span class='badge {s['scope']}'>{BADGE[s['scope']]}</span></td></tr>")
    for a in ARTS:
        rows.append(f"<tr><td><a href='#{a['sid']}'>✕</a></td><td>{esc(a['title'])}</td><td>see section</td>"
                    f"<td><code>NON_TPRT / ARTEFACT_LIKE</code></td><td>—</td><td><span class='badge art'>{BADGE[ART]}</span></td></tr>")
    body.append(f"""<section id="summary"><h2>Summary</h2>
<p class="note">Frequency column abbreviates the sources: Z = Zumalave 2026, PCAWG = Rodriguez-Martin 2020, Nam = Nam 2023; see each section for denominators and the other sources.</p>
<div class="tbl"><table><thead><tr><th>#</th><th>type</th><th>short-read signature</th><th>vocabulary</th><th>frequency</th><th>scope</th></tr></thead>
<tbody>{''.join(rows)}</tbody></table></div></section>""")

    refs = "".join(f"<li>{esc(t)}. <a href='https://doi.org/{d}'>doi:{d}</a></li>" for t, d in REFS)
    body.append(f"""<section id="refs"><h2>References</h2><ol class="refs">{refs}</ol>
<p class="note">Generated by <code>docs/src/make_insertion_types.py</code>; edit the script and re-run, do not edit the HTML by hand.</p></section>""")

    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Retrotransposition insertion types</title>
<meta name="description" content="Catalogue of somatic retrotransposition insertion types with genomic and short-read figures, frequencies and PEAR-TREE vocabulary.">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&amp;family=Source+Sans+3:wght@400;600;700&amp;display=swap" rel="stylesheet">
<style>{CSS}</style>
</head>
<body>
<div class="wrap">
<nav class="toc" aria-label="Contents"><p class="h">Contents</p><ul>{''.join(toc)}</ul></nav>
<main>
{''.join(body)}
</main>
</div>
</body>
</html>
"""
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(page)
    print("wrote", OUT)


if __name__ == "__main__":
    build()
