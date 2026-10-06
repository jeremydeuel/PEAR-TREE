"""Orchestrates the tools/rte annotation of one patient's insertions (annotate_v2 plug-in).

Activation (all under CONFIG['annotate']):
    rte_library        path to resources/rte_library (or a stand-in); REQUIRED to enable
    rte_evidence_file  callable(sample)->path or path; default: insertions_file with
                       `.combined.txt.gz` -> `.insertions.evidence.tsv.gz` (optional file)
    rte_reads_file     same, `.insertions.reads.fa.gz` (optional file)
    genome_2bit        discovery genome (GRCh38/hg19/hs1) for TSD / EN motif / slippage /
                       templated / pre-mRNA (optional; without it TSD falls back to the flanks)
    remap_2bit         remap genome (hs1) for exon-junction cores and novel-source identity
    remap_index        minimap2 index/FASTA of the remap genome (novel-source locator)
    remap_rmsk         RepeatMasker .out / rmsk.txt of the remap genome (novel-source rule)
    exon_annotation    exon track on the clip-remap genome (already used by annotate_v2 for the
                       pseudogene CANDIDATE genes; hs1 in a real run)
    rte_exon_annotation  exon track on `remap_2bit` used for the exon-exon junction cores; default
                       exon_annotation (they differ only in the E2E, whose clip remap is a reduced
                       GRCh38 while remap_2bit is hs1)
    rte_score          {'weights': {...}, 'thresholds': {...}} overrides of score.py
    rte_structure / rte_assembly / rte_transduction / rte_pseudogene   option overrides
    rte_recurrence_max (3) loci sharing one truncation signature before 'recurrence' fires
    rte_premrna_window (1_000_000) bp around the site searched for pre-mRNA templates (gene model)
    rte_wide_window    (10_000) bp either side of the site in which an unexplained insert piece
                       is looked up as a local template (TEMPLATED_LOCAL <= 250 bp away, beyond
                       that a co-inserted local pre-mRNA, PREMRNA_COINSERT)
    rte_max_reads      (400) reads per insertion used for the assembly
"""
from __future__ import annotations

import os
from collections import Counter
from dataclasses import dataclass, field

from .assembly import Assembler, SiteContext
from .genome import open_genome
from .hallmarks import (split_junction, polya_info, locate_site, target_site, en_motif,
                        slippage_context, foldback, parse_locus)
from .inputs import read_evidence_tsv, read_reads_fa, InsertionEvidence
from .library import RteLibrary
from .pseudogene import ExonJunctionIndex, load_exons_by_gene, load_gene_strands
from .record import RteRecord
from .score import ScoreInput, score
from .structure import classify
from .transduction import NovelSourceFinder, MappyLocator


@dataclass
class InsertionInput:
    title: str
    left_seq: str = ""
    right_seq: str = ""
    pseudogene_genes: list = field(default_factory=list)
    legacy_class: str | None = None
    sv: tuple | None = None          # annotate_v2 _sv_subtype() result (rank, desc, contig, pos)

    @classmethod
    def from_legacy(cls, ins, legacy_class=None):
        genes = []
        try:
            pg = ins._pseudogene()
            if pg:
                genes.append(pg[0])
        except Exception:
            pass
        for g, _side, _n in getattr(ins, "splice_hits", []) or []:
            if g not in genes:
                genes.append(g)
        # every gene with an exon under a clip remap is a candidate for the exon-exon junction
        # search: the junction read is the proof, the candidate list only bounds the search
        # (_pseudogene() also wants a poly-A on the clip, which a jittered tail can miss)
        for g, _s, _e in (getattr(ins, "left_exons", []) or []) + (getattr(ins, "right_exons", []) or []):
            if g not in genes:
                genes.append(g)
        sv = None
        try:
            sv = ins._sv_subtype(allow_rte=True)
        except Exception:
            sv = None
        return cls(ins.title, ins.left_seq or "", ins.right_seq or "", genes, legacy_class, sv)


def title_gap(title):
    """R - L from a numeric locus name `contig:L-R` (None for polyA_/oneside_/disc_ tokens)."""
    loc = parse_locus(title)
    return None if loc is None else loc[2] - loc[1]


def _resolve_path(spec, sample):
    if spec is None:
        return None
    return spec(sample) if callable(spec) else spec


def default_sidecars(insertions_file):
    base = insertions_file
    for suf in (".combined.txt.gz", ".txt.gz"):
        if base.endswith(suf):
            base = base[: -len(suf)]
            break
    return base + ".insertions.evidence.tsv.gz", base + ".insertions.reads.fa.gz"


class RteAnnotator:
    def __init__(self, cfg: dict, genome=None, remap_genome=None, locator=None, rmsk=None,
                 exons_by_gene=None, gene_model=None, cohort_l1=None):
        self.cfg = cfg
        self.lib = RteLibrary(cfg["rte_library"], cfg)
        self.assembler = Assembler(self.lib, cfg.get("rte_assembly"))
        self.genome = open_genome(genome if genome is not None else cfg.get("genome_2bit"))
        self.remap = open_genome(remap_genome if remap_genome is not None else cfg.get("remap_2bit"))
        if locator is None and cfg.get("remap_index"):
            locator = MappyLocator(cfg["remap_index"])
        if rmsk is None:
            rmsk = cfg.get("remap_rmsk")
        self.novel = NovelSourceFinder(self.lib, cfg.get("rte_transduction"), rmsk, locator,
                                       self.remap, cohort_l1)
        exon_track = cfg.get("rte_exon_annotation") or cfg.get("exon_annotation")
        if exons_by_gene is None and exon_track and os.path.exists(str(exon_track)):
            exons_by_gene = load_exons_by_gene(exon_track)
        strands = load_gene_strands(exon_track) if exon_track and os.path.exists(str(exon_track)) else {}
        self.exon_index = ExonJunctionIndex(exons_by_gene or {}, self.remap, cfg.get("rte_pseudogene"),
                                            strands)
        self.gene_model = gene_model
        self.evidence = {}

    # ------------------------------------------------------------------ inputs
    def load_evidence(self, evidence_path=None, reads_path=None, wanted=None):
        if evidence_path and os.path.exists(evidence_path):
            read_evidence_tsv(evidence_path, self.evidence)
            print(f"[rte] read junction evidence for {len(self.evidence)} insertions from {evidence_path}")
        if reads_path and os.path.exists(reads_path):
            read_reads_fa(reads_path, self.evidence, wanted)
            print(f"[rte] read pooled reads from {reads_path}")

    # ------------------------------------------------------------------ pre-mRNA helper
    def _premrna_fn(self, site):
        if self.genome is None or self.gene_model is None or site.contig is None:
            return None
        anchor = site.L if site.L is not None else site.R
        if anchor is None:
            return None
        win = int(self.cfg.get("rte_premrna_window", 1_000_000))
        lo = max(0, anchor - win)
        state = {}

        def fn(seq):
            if "al" not in state:
                import mappy
                ref = self.genome.fetch(site.contig, lo, anchor + win)
                state["al"] = mappy.Aligner(seq=ref, preset="sr") if ref else None
            al = state["al"]
            if al is None:
                return None
            for h in al.map(seq):
                if h.mlen / max(1, h.blen) < 0.9:
                    continue
                p = lo + (h.r_st + h.r_en) // 2
                if abs(p - anchor) <= 250:
                    continue      # that is a local template, not a distal pre-mRNA
                cands = self.gene_model._candidates(site.contig, p) or []
                for gs, ge, name, strand, exons in cands:
                    if gs <= p < ge:
                        feat = self.gene_model._genic_feature(exons, strand, p)[1]
                        return f"{name}:{feat}@{site.contig}:{p}"
            return None
        return fn

    # ------------------------------------------------------------------ one insertion
    def annotate(self, inp: InsertionInput, ev: InsertionEvidence | None = None) -> RteRecord:
        rec = RteRecord(inp.title)
        ev = ev or InsertionEvidence(inp.title)
        jl, jr = ev.junctions.get("LEFT"), ev.junctions.get("RIGHT")
        left_str = jl.clip_consensus if (jl and jl.clip_consensus) else inp.left_seq
        right_str = jr.clip_consensus if (jr and jr.clip_consensus) else inp.right_seq
        li, lf, rf, ri = split_junction(left_str, right_str)
        pa = polya_info(li, ri, jl, jr)
        site = locate_site(inp.title, lf, rf, self.genome)
        target_site(site, lf, rf, self.genome)
        ctx = SiteContext(inp.title, site.contig, site.L, site.R, left_flank=lf, right_flank=rf)
        if self.genome is not None and (site.L is not None or site.R is not None):
            w = int(self.assembler.cfg["local_window"])
            pts = [p for p in (site.L, site.R) if p is not None]
            ctx.window_start = max(0, min(pts) - w)
            ctx.window_seq = self.genome.fetch(site.contig, ctx.window_start, max(pts) + w)
            ww = int(self.cfg.get("rte_wide_window", 10000))
            if ww > w:
                ctx.wide_start = max(0, min(pts) - ww)
                ctx.wide_seq = self.genome.fetch(site.contig, ctx.wide_start, max(pts) + ww)
        junction_seqs = {}
        if left_str:
            junction_seqs["LEFT"] = (left_str, (len(li), len(left_str)))
        if right_str:
            junction_seqs["RIGHT"] = (right_str, (0, len(rf)))
        reads = self._cap_reads(ev.reads, int(self.cfg.get("rte_max_reads", 400)))
        hint = (pa.strand, pa.source) if pa.strand else None
        asm = self.assembler.assemble(ctx, junction_seqs, reads, hint)
        strand = asm.strand
        polya_reads = None
        if not pa.strand and asm.strand_source == "polya_reads":
            # the junction strings carry no poly-A but clip reads do (A-run | REF): that tail set
            # the strand. Reported (median run over the voting reads; a lower bound, a read can
            # end inside the tail) but NOT scored: the score's poly-A weights are calibrated on
            # the junction/evidence poly-A, and 178/207 PD37590 loci have none there, so
            # scoring the read tails would shift most of the patient (artefacts included)
            _st, lens = asm._strand_from_polya(asm.raw_layouts)
            if lens and _st == strand:
                polya_reads = sorted(lens)[len(lens) // 2]
        en_motif(site, strand, self.genome)
        slippage_context(site, strand, self.genome)
        # pseudogene proof
        pg_hits = []
        if inp.pseudogene_genes:
            seqs = [(f"junction_{k}", v[0].upper()) for k, v in junction_seqs.items()]
            seqs += [(f"{r.role}|{r.sample}|{r.frag}", r.seq) for r in reads]
            pg_hits = self.exon_index.find(inp.pseudogene_genes, seqs)
        call = classify(asm, self.lib, ctx, self.cfg.get("rte_structure"), self.novel,
                        self._premrna_fn(site),
                        (inp.pseudogene_genes, pg_hits, self._pseudogene_structure_fn(inp.pseudogene_genes)),
                        inp.legacy_class)
        rte_elem = call.element in ("L1", "ALU", "SVA")
        # ---- site-level tags. Discovery pairing modes (SPEC "Pairing modes and locus names")
        # are recoverable from the locus-name geometry gap = R - L: [-30, -1] target-site
        # deletion, < -30 L1-mediated deletion, {0, 1} blunt, > 40 L1-mediated duplication OR a
        # long-TSD chimera (same geometry; decided here from element + poly-A polarity).
        gap = site.tsd_len if site.tsd_len is not None else title_gap(inp.title)
        if gap is not None and site.tsd_len is None:
            site.tsd_len = gap                    # genome-free fallback: the locus geometry
        md = self.cfg.get("rte_max_target_site_deletion", 30)
        if gap is not None and -md <= gap < 0:
            call.add("TSD_DELETION")
        if rte_elem:
            # an intrachromosomal clip partner explained by a co-inserted local template /
            # pre-mRNA is that template, not an L1-mediated rearrangement (E2E: 3/5
            # PREMRNA_COINSERT TPs were tagged L1_MED_DUPLICATION)
            sv_intra = (inp.sv is not None and inp.sv[0] == 2
                        and not {"PREMRNA_COINSERT", "TEMPLATED_LOCAL"} & set(call.tags))
            polarised = pa.length >= 10 and not pa.both_sided and "CHIMERIC_ENDS" not in call.tags
            if (gap is not None and gap < -md) or (sv_intra and (gap or 0) <= 0):
                call.add("L1_MED_DELETION")
            elif (sv_intra and (gap or 0) > 0) or (gap is not None and gap > 40 and polarised):
                call.add("L1_MED_DUPLICATION")
            # EN-independent (Morrish 2002): no poly-A, no TSD (blunt = gap 0/1, SPEC pairing
            # modes), 3' end short of the consensus end
            if (pa.length < 10 and (site.tsd_len is None or site.tsd_len <= 1)
                    and call.three_prime_short and call.structure != "FULL_LENGTH"):
                call.add("EN_INDEPENDENT")
        # ---- beyond poly-A on the 3' side
        side3 = "LEFT" if strand >= 0 else "RIGHT"
        ev3 = ev.junctions.get(side3)
        bseq, bsup = self._beyond_polya(asm)
        if ev3 is not None and ev3.beyond_polya:
            if not bseq:
                bseq = ev3.beyond_polya
            bsup = max(bsup, ev3.beyond_polya_support)
        # ---- support
        supported = sum(1 for j in ev.junctions.values() if j.supported)
        n_samples = max((j.n_samples for j in ev.junctions.values()), default=0)
        csi = any(j.cross_sample_identical for j in ev.junctions.values())
        fb = foldback(li, lf, ri, rf)
        td_match = False
        if call.source is not None and not call.source.novel:
            se = self.lib.source_element(call.source.source_id)
            td_match = se is not None and se == asm.nearest_intact and call.j5_class == "L1"
        # ---- record
        rec.element, rec.structure, rec.tags = call.element, call.structure, list(call.tags)
        rec.detail = dict(call.detail)
        rec.consensus = asm.consensus
        rec.strand = strand
        if asm.covered:
            rec.covered_5p, rec.covered_3p = asm.covered_5p, asm.covered_3p
            rec.covered_intervals = list(asm.covered)
            rec.covered_seqs = list(asm.covered_seqs)
        if asm.nearest_active != ".":
            rec.element_identity = asm.element_identity
            rec.nearest_active = asm.nearest_active
        rec.site = (site.contig, site.L, site.R)
        rec.tsd_seq, rec.tsd_len = site.tsd_seq, site.tsd_len
        rec.en_motif, rec.en_mismatches = site.en_motif, site.en_mismatches
        rec.polya_len = pa.length if pa.length else 0.0
        rec.beyond_polya, rec.beyond_polya_support = bseq, bsup
        if site.slippage:
            rec.detail["slippage"] = site.slippage_detail
        if polya_reads is not None:
            rec.detail["polya_reads"] = polya_reads
        if fb:
            rec.detail["foldback_bp"] = fb
        n_short = sum(j.n_short_used for j in ev.junctions.values())
        if n_short:      # combine counted SHORT-overhang fragments (count_short_overhang)
            rec.detail["short_used"] = n_short
            rec.detail["short_mate_inside"] = sum(j.n_short_mate_inside for j in ev.junctions.values())
        rec.score_input = ScoreInput(
            element=call.element, structure=call.structure, tags=list(call.tags),
            tsd_len=site.tsd_len, tsd_verified=site.tsd_verified, polya_len=rec.polya_len,
            polya_both_sides=pa.both_sided, slippage=site.slippage,
            beyond_polya_len=len(bseq), beyond_polya_support=bsup,
            en_mismatches=site.en_mismatches,
            ends_concordant=bool(call.j5_class) and call.j5_class == call.j3_class,
            td_source_matches_5p=td_match,
            element_identity=asm.element_identity or 0.0,
            inactive_only=bool(asm.consensus) and not self.lib.is_young(asm.consensus),
            inv_p1=call.inv_p1, junctions_supported=supported, n_samples=n_samples,
            foldback=bool(fb), recurrent=False, cross_sample_identical=csi,
            novel_tier=(call.source.tier if call.source is not None else ""))
        self._score(rec)
        return rec

    def _pseudogene_structure_fn(self, genes):
        """FULL_LENGTH / TRUNCATED_5P of a proven pseudogene from the inserted sequence at the
        5' junction (element-sense layouts) against the spliced transcript."""
        if not genes:
            return None
        tol = int((self.cfg.get("rte_structure") or {}).get("pseudogene_full_length_tol", 15))

        def fn(layouts):
            seqs = []
            for lay in sorted(layouts, key=lambda l: l.role != "JUNCTION"):
                segs = lay.segments
                if len(segs) >= 2 and segs[0].kind == "REF":
                    seqs.append(lay.seq[segs[0].q_en:])
            return self.exon_index.structure(genes, seqs, tol=tol)
        return fn

    @staticmethod
    def _cap_reads(reads, cap):
        """Keep every junction read (CLIP/POLYA/SPAN); thin the rest evenly so coverage along
        the element is preserved (deterministic)."""
        if len(reads) <= cap:
            return list(reads)
        keep = [r for r in reads if r.role in ("CLIP", "POLYA", "SPAN")][:cap]
        rest = [r for r in reads if r.role not in ("CLIP", "POLYA", "SPAN")]
        room = cap - len(keep)
        if room > 0 and rest:
            step = len(rest) / room
            keep += [rest[int(i * step)] for i in range(room)]
        return keep

    def _score(self, rec):
        sc = self.cfg.get("rte_score") or {}
        rec.tprt_score, rec.tprt_points, rec.tprt_call = score(
            rec.score_input, sc.get("weights"), sc.get("thresholds"))

    @staticmethod
    def _beyond_polya(asm):
        """Sequence on the far side of the 3' poly-A (element sense: ...X | POLYA | REF) and the
        number of independent fragments (distinct sample+frag) carrying >= 10 bp of it."""
        best = ""
        frags = set()
        for lay in asm.layouts:
            segs = lay.segments
            if len(segs) < 3 or segs[-1].kind != "REF" or segs[-2].kind != "POLYA":
                continue
            x = segs[-3]
            if x.kind in ("POLYA", "REF", "LOCAL") or x.qlen < 10:
                continue
            s = lay.seq[x.q_st:x.q_en]
            if lay.role == "JUNCTION":
                best = s if len(s) > len(best) else best
            else:
                frags.add(lay.frag_key)
                if not best:
                    best = s
        return best, len(frags)

    # ------------------------------------------------------------------ all insertions
    def annotate_all(self, inputs: dict) -> dict:
        out = {}
        for key, inp in inputs.items():
            try:
                out[key] = self.annotate(inp, self.evidence.get(key))
            except Exception as e:      # never let one locus kill the table
                rec = RteRecord(key)
                rec.detail["error"] = f"{type(e).__name__}: {e}"[:200]
                out[key] = rec
        self.cohort_source_pass(inputs, out)
        self.recurrence_pass(out)
        return out

    def cohort_source_pass(self, inputs, records):
        """SPEC novel-source rule, cohort branch: an unsourced 3' transduction may come from an
        L1 insertion called elsewhere in the same cohort. Only possible when the discovery genome
        IS the remap genome (cohort sites and tag placements share coordinates)."""
        g, r = self.cfg.get("genome_2bit"), self.cfg.get("remap_2bit")
        if not g or g != r or self.novel.locator is None:
            return
        l1 = []
        for rec in records.values():
            if rec.element == "L1" and rec.tprt_call in ("TPRT", "LIKELY_TPRT") and rec.site[0]:
                pos = rec.site[1] if rec.strand >= 0 else rec.site[2]
                if pos is not None:
                    l1.append((rec.site[0], pos, "+" if rec.strand >= 0 else "-"))
        redo = [k for k, rec in records.items()
                if "TD3P" in rec.tags and not any(t.startswith("TD3P_SOURCE=") for t in rec.tags)]
        if not l1 or not redo:
            return
        self.novel.cohort_l1 = l1
        for k in redo:
            try:
                records[k] = self.annotate(inputs[k], self.evidence.get(k))
            except Exception:
                pass

    def recurrence_pass(self, records: dict):
        """Same 5' truncation signature (class, structure, junction position +-2 bp) at more
        than `rte_recurrence_max` unrelated loci = one molecule/artefact seen many times (index
        hopping / mapping artefact), not independent insertions. Full-length elements are
        excluded (they legitimately share the consensus 5' end)."""
        mx = int(self.cfg.get("rte_recurrence_max", 3))
        sig = Counter()
        keys = {}
        for k, r in records.items():
            j5 = r.detail.get("j5", r.detail.get("fwd_start"))
            if r.element not in ("L1", "ALU", "SVA") or r.structure in ("FULL_LENGTH", "5P_UNRESOLVED") or j5 is None:
                continue
            s = (r.element, r.structure, int(j5) // 5)
            keys[k] = s
            sig[s] += 1
        for k, s in keys.items():
            if sig[s] > mx and records[k].score_input is not None:
                records[k].score_input.recurrent = True
                records[k].detail["recurrence"] = sig[s]
                self._score(records[k])
