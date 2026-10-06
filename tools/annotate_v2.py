# PEAR-TREE - paired ends of aberrant retrotransposons in phylogenetic trees
#
# Copyright (C) 2025 Jeremy Deuel <jeremy.deuel@usz.ch>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.


# this is an updated version of annotate.py that does not require excel file generation and also is more sophisticated in removing false positives.

DEBUG = True

import gzip
import os
import pysam
from math import floor, log2
import re
import subprocess
import sys
from src.config import CONFIG

class RepeatMasker_Annotation:
    def __init__(self, line):
        self.score = float(line[0])
        self.millidiv = float(line[1]) if line[1] != ' ' else None
        self.millidel = float(line[2]) if line[2] != ' ' else None
        self.milliins = float(line[3]) if line[3] != ' ' else None
        self.seqname = line[4]
        self.start = int(line[5])
        self.end = int(line[6])
        self.strand = line[8]
        if self.strand == 'C': self.strand = "-" #legacy annotation
        self.repName = line[9]
        cl = line[10].split("/", maxsplit=1)
        if len(cl) == 2:
            self.repClass, self.repFamily = cl
        else:
            self.repClass = cl
            self.repFamily = None
        self.repStart = int(line[11]) if line[11][0] != '(' else int(line[11][1:-1])
        self.repEnd = int(line[12]) if line[12] != '.' else None

    def __str__(self):
        return f"{self.repName}|{self.repFamily}|{self.repClass}|{self.strand}"
class Dfam_Annotation:
    def __init__(self, line):
        self.model = line[0]
        self.is_active = re.match(r"^L1(HS|PA)", self.model) or re.match(r"^Alu", self.model) or re.match(r"^SVA",self.model)
        self.accession = line[1]
        self.bits = float(line[3])
        self.evalue = float(line[4])
        self.bias = float(line[5])
        self.hmm_start = int(line[6])
        self.hmm_end = int(line[7])
        self.strand = line[8]
        self.ali_start = int(line[9])
        self.ali_end = int(line[10])
        self.env_start = int(line[11])
        self.env_end = int(line[12])
        self.model_length = int(line[13])
        self.description = line[14]

    def __str__(self):
        return f"{self.model} {self.strand} ({self.hmm_start}-{self.hmm_end}|{self.model_length})"


class GeneModel:
    """Insertion-SITE annotator: given the reference locus where an insertion LANDED (the title
    "contig:start-end"), report whether it sits inside a gene and, if so, which gene and which
    genic feature -- exon, splice donor / acceptor / polypyrimidine tract / branch point (lariat),
    or deep intron -- and whether it disrupts or lies near a gene's promoter (its TSS).

    This is ORTHOGONAL to every existing annotation. The Dfam/rmsk/pseudogene/SV logic answers
    "what was inserted" (element identity) using the clip sequence and where the clip REMAPS (in
    bowtie2_index2); the site annotation answers "where did it land" using the junction's own
    reference locus (the title), which lives on the SAMPLE's discovery/BAM genome. So the gene
    model must be built for that genome (GRCh38 or GRCh37), and its contigs must match the title
    contigs -- lookup is chr-prefix tolerant (chr8 <-> 8) so one track serves both conventions.

    Track format (build_gene_model.py): a (optionally gzipped) TSV `contig<TAB>start<TAB>end<TAB>
    gene<TAB>strand`, 0-based half-open, ONE ROW PER (merged) EXON, sorted. Everything -- the gene
    span, the TSS, and every intron/exon boundary -- is reconstructed from the per-gene exon rows,
    so a single file drives all features (an unstranded exon-only track, e.g. the pseudogene track,
    is NOT enough: promoter needs the TSS = strand + gene 5' end, and donor vs acceptor needs the
    strand). None in config = site annotation off; no other behaviour changes.
    """

    def __init__(self, path, cfg=None):
        cfg = cfg or {}
        # Splice windows, measured in nt INTO the intron from the exon boundary (sense-aware).
        #   donor    : the 5' splice site (GT..) -- intron positions +1..+donor_window from an
        #              exon's 3' end.
        #   acceptor : the 3' splice site (..AG) -- intron positions -1..-acceptor_window before
        #              an exon's 5' end.
        #   ppt      : the polypyrimidine tract, just upstream of the acceptor.
        #   branch   : the branch-point A / lariat, further upstream still.
        # The three 3' bins are nested distance thresholds (acceptor < ppt < branch), so the
        # tightest matching one wins. Defaults follow standard splicing anatomy; override per
        # deployment via CONFIG['annotate'][...].
        self.donor_window    = cfg.get('splice_donor_window', 6)
        self.acceptor_window = cfg.get('splice_acceptor_window', 3)
        self.ppt_window      = cfg.get('splice_ppt_window', 17)
        self.branch_window   = cfg.get('splice_branch_window', 45)
        # Promoter windows around the TSS (sense-aware: *_up is 5' of the TSS, *_down is 3').
        # core = "disrupts promoter"; the wider proximal window = "near promoter".
        self.prom_core_up    = cfg.get('promoter_core_up', 250)
        self.prom_core_down  = cfg.get('promoter_core_down', 250)
        self.prom_up         = cfg.get('promoter_up', 2000)
        self.prom_down       = cfg.get('promoter_down', 500)
        self.genes = {}        # contig -> [(start, end, name, strand, exons), ...] sorted by start
        self._maxspan = {}     # contig -> longest gene on it (bounds the overlap back-scan)
        if path is not None:
            self._load(path)

    # -- contig-name tolerance (title "chr8" vs a numeric-contig GRCh37 track, or vice versa) --
    def _resolve(self, contig):
        if contig in self.genes:
            return contig
        alt = contig[3:] if contig.startswith("chr") else "chr" + contig
        return alt if alt in self.genes else None

    def _load(self, path):
        print(f"reading gene model {path}")
        acc = {}   # (contig, gene) -> [strand, [(s, e), ...]]
        opener = gzip.open if str(path).endswith(".gz") else open
        n = 0
        with opener(path, 'rt') as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line[0] == '#':
                    continue
                f = line.split("\t")
                if len(f) < 5:
                    continue
                contig, start, end, gene, strand = f[0], int(f[1]), int(f[2]), f[3], f[4]
                a = acc.setdefault((contig, gene), [strand, []])
                a[1].append((start, end))
                n += 1
        for (contig, gene), (strand, ivs) in acc.items():
            ivs.sort()
            merged = [list(ivs[0])]
            for s, e in ivs[1:]:                 # defensively merge (builder already merges)
                if s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            exons = tuple((s, e) for s, e in merged)
            gs, ge = exons[0][0], exons[-1][1]
            self.genes.setdefault(contig, []).append((gs, ge, gene, strand, exons))
        for contig in self.genes:
            self.genes[contig].sort()
            self._maxspan[contig] = max(ge - gs for gs, ge, *_ in self.genes[contig])
        print(f"imported {n} exon rows for {sum(len(v) for v in self.genes.values())} genes "
              f"on {len(self.genes)} contigs from {path}")

    def _candidates(self, contig, p):
        """Genes on `contig` whose body comes within prom_up of p (so both body overlaps and
        promoter-proximity are covered). Returns [] for a locus far from any gene, None if the
        contig is absent from the model (can't judge -> caller treats as intergenic)."""
        key = self._resolve(contig)
        if key is None:
            return None
        genes = self.genes[key]
        reach = self.prom_up                      # widest window either side of p we care about
        hi = p + reach
        # binary search: rightmost gene whose start <= hi
        lo_i, hi_i = 0, len(genes)
        while lo_i < hi_i:
            mid = (lo_i + hi_i) // 2
            if genes[mid][0] <= hi:
                lo_i = mid + 1
            else:
                hi_i = mid
        # scan left over the window, bounded by the longest gene (a long gene can start far left
        # and still reach p); collect those whose body/promoter window actually reaches p.
        floor = p - reach - self._maxspan.get(key, 0)
        out = []
        for i in range(lo_i - 1, -1, -1):
            gs, ge, name, strand, exons = genes[i]
            if gs < floor:
                break
            if ge + reach >= p:
                out.append((gs, ge, name, strand, exons))
        return out

    def _splice_class(self, dd, da):
        """Classify an intronic point from its distance (nt) to the DONOR boundary (dd) and the
        ACCEPTOR boundary (da). Tightest match wins."""
        if dd <= self.donor_window and dd <= da:
            return (1, 'splice_donor', 'splice donor')
        if da <= self.acceptor_window:
            return (2, 'splice_acceptor', 'splice acceptor')
        if da <= self.ppt_window:
            return (3, 'polypyrimidine_tract', 'polypyrimidine tract')
        if da <= self.branch_window:
            return (4, 'branch_point', 'branch point (lariat)')
        if dd <= self.donor_window:
            return (1, 'splice_donor', 'splice donor')
        return (5, 'intron', 'intron')

    def _genic_feature(self, exons, strand, p):
        """(rank, keyword, label) for point p inside a gene body. rank orders disruptiveness so
        the most salient feature wins when genes overlap (exon < donor < acceptor < ppt < branch
        < intron)."""
        for s, e in exons:
            if s <= p < e:
                return (0, 'exon', 'exonic')
        for i in range(len(exons) - 1):
            le = exons[i][1]                       # first intron base is le (0-based)
            rs = exons[i + 1][0]                   # last intron base is rs-1
            if le <= p < rs:
                d_left = p - le + 1                # nt into the intron from the left boundary
                d_right = rs - p                   # nt into the intron from the right boundary
                # + strand: exon 3' end is the genomic-left boundary -> donor on the left.
                # - strand: everything is reversed.
                if strand == '+':
                    dd, da = d_left, d_right
                else:
                    dd, da = d_right, d_left
                return self._splice_class(dd, da)
        return (5, 'intron', 'intron')             # single-exon genes never reach here

    def _promoter(self, gs, ge, strand, p):
        """(abs_dist, keyword, label) if p is within a promoter window of this gene's TSS, else
        None. TSS = gene 5' end (start for +, end for -). Windows are sense-aware."""
        tss = gs if strand == '+' else ge
        off = (p - tss) if strand == '+' else (tss - p)   # <0 = upstream (5') of the TSS
        if -self.prom_core_up <= off <= self.prom_core_down:
            return (abs(off), 'promoter_core', 'disrupts promoter')
        if -self.prom_up <= off <= self.prom_down:
            return (abs(off), 'promoter_proximal', 'near promoter')
        return None

    def annotate(self, contig, lo, hi):
        """Classify the insertion site spanning [lo, hi] on `contig`. Returns
        (region_keyword, gene, strand, human_text). Uses the locus midpoint as the representative
        point (the target-site window is a few bp; every feature window is much wider)."""
        p = (lo + hi) // 2
        cands = self._candidates(contig, p)
        if cands is None or not cands:
            return ('intergenic', '.', '.', 'intergenic')
        best_genic = None   # (rank, kw, label, gene, strand)
        best_prom = None    # (dist, kw, label, gene, strand)
        for gs, ge, name, strand, exons in cands:
            if gs <= p < ge:
                rank, kw, label = self._genic_feature(exons, strand, p)
                if best_genic is None or rank < best_genic[0]:
                    best_genic = (rank, kw, label, name, strand)
            prom = self._promoter(gs, ge, strand, p)
            if prom is not None and (best_prom is None or prom[0] < best_prom[0]):
                best_prom = (prom[0], prom[1], prom[2], name, strand)
        parts = []
        region_kw, gene_out, strand_out = 'intergenic', '.', '.'
        if best_genic is not None:
            _, region_kw, label, gene_out, strand_out = best_genic
            parts.append(f"{label} of {gene_out} ({strand_out})")
        if best_prom is not None:
            _, pkw, plabel, pgene, pstrand = best_prom
            # Suppress a redundant promoter note for the gene we already reported as genic --
            # unless it is a *core* disruption, which is a distinct functional statement worth
            # keeping even when we already know the locus is in that gene (e.g. landed in exon 1).
            same_gene = best_genic is not None and pgene == gene_out
            if not same_gene or pkw == 'promoter_core':
                parts.append(f"{plabel} of {pgene} ({pstrand})")
            if best_genic is None:
                region_kw, gene_out, strand_out = pkw, pgene, pstrand
        if not parts:
            return ('intergenic', '.', '.', 'intergenic')
        return (region_kw, gene_out, strand_out, '; '.join(parts))


class Insertion:
    # Shared insertion-SITE annotator (GeneModel), set once by the container from
    # CONFIG['annotate']['gene_model']. None => no site annotation (default), so conclusion()
    # appends nothing and every existing test that asserts an exact conclusion string is
    # unaffected. See GeneModel and _site_annotation().
    gene_model = None
    def __init__(self, title, left_seq, right_seq):
        self.title = title
        self.left_seq = left_seq
        self.right_seq = right_seq
        self.nins = 0
        self.nart = 0
        self.nwt = 0
        self.right_dfams = []
        self.left_dfams = []
        self.right_maps = []
        self.left_maps = []
        # (E) partial placements from the bowtie2 --local pass, populated by read_sam_local()
        # only. Each is (contig, ref_pos, strand, query_coverage_bp, mapq). Empty unless the
        # local-remap channel ran, so every consumer is a no-op on --end-to-end-only data.
        self.right_local_maps = []
        self.left_local_maps = []
        # exon hits of the mapped clips, for processed-pseudogene detection:
        # each is (gene_id, exon_start, exon_end). Empty unless an exon track is configured.
        self.right_exons = []
        self.left_exons = []
        # discovery splice-hallmark evidence (Feature B), aggregated by combine_insertions:
        # each is (gene_id, side, n_exons) — mates span >= n_exons exons of one gene.
        self.splice_hits = []
        # reciprocal translocation partner locus ("contig:pos"), set by the container's
        # link_reciprocal_translocations() pass when this junction and another point back at
        # each other (a balanced/reciprocal translocation). None until/unless that pass runs.
        self.reciprocal_partner = None
        # processed-pseudogene proof from tools/rte (an exon-exon junction covered by a read):
        # None = not evaluated (rte off), True/False once the rte pass ran. Only consulted by
        # _element_conclusion() when CONFIG['annotate']['pseudogene_require_exon_junction'] is set.
        self.exon_junction_proven = None
        # extract inserted sequences
    # Sequence case encodes the junction: UPPER = aligned to the reference, lower = clipped.
    # A tail is therefore a homopolymer run in the *clipped* part, flush against the aligned
    # part: `a{N}[ACGT]` at the end of a left clip, `[ACGT]t{N}` at the start of a right clip.
    #
    # N was hardcoded at 6, which is far too permissive: a real TPRT poly-A tail is 15-40 bp,
    # while an incidental 6-mer A/T run occurs throughout ordinary sequence -- A-rich
    # Low_complexity tracts and the A-rich 3' ends of L1/MIR/L2 all carry one -- so junctions
    # with no tail at all scored poly-A positive. Sweep on the fp10k FP-stress harness
    # (true insertions vs hallmark-free decoys), fraction scored poly-A positive:
    #     min_run:   6      12     14     20
    #     true:    97.1%  97.0%  97.0%  93.7%
    #     decoys:  13.3%   7.3%   5.3%   2.4%
    # 12 (the default) removes ~45% of the false poly-A at ~0.1% cost. The true-call curve is
    # flat out to 14 partly because that harness's shortest simulated tail is 15 bp, so 12
    # keeps margin for genuinely short/truncated tails in real data.
    #
    # NB this only changes what gets *called* once `require_polya_hallmark` is on: with the
    # homology fallbacks open (the default), conclusion() accepts without consulting the
    # poly-A verdict at all, so tightening N alone moves ~1 call on the harness.
    # Override with CONFIG['annotate']['polya_min_len'] (src/config.py is deployment-local).
    @staticmethod
    def _polya_min_len():
        return CONFIG['annotate'].get('polya_min_len', 12)

    def has_right_polyA(self):
        return re.search(r"[ACGT]t{%d}" % self._polya_min_len(), self.right_seq)
    def has_left_polyA(self):
        return re.search(r"a{%d}[ACGT]" % self._polya_min_len(), self.left_seq)

    # SVA is a composite element (CCCTCT hexamer + Alu-like region + VNTR + SINE-R + poly-A).
    # Its diagnostic 5' motif is the (CCCTCT)n hexamer (revcomp (AGAGGG)n), which Alu and L1
    # lack entirely -- so a hexamer tandem in a clip is an SVA-specific fingerprint. A single
    # 6-mer occurs by chance in ~5% of clips, so it is only ever corroborating, never a sole
    # trigger; the tandem (>= 2 copies, 12 bp) is specific enough to act as a hallmark.
    _SVA_HEX_TANDEM = re.compile(r"(?:CCCTCT){2,}|(?:AGAGGG){2,}", re.I)
    _SVA_HEX_ANY = re.compile(r"CCCTCT|AGAGGG", re.I)

    # Minimum SVA-model bit score for a clip's SVA hit to count as the *dominant* element
    # identity (over a competing Alu on the same clip). The Alu-like body of an SVA scores an
    # Alu model, so a trace SVA hit (< a few bits) on a clip that is really Alu must not win;
    # a genuine SVA VNTR/SINE-R clip scores the SVA models in the tens of bits. Override with
    # CONFIG['annotate']['sva_min_bits'].
    @staticmethod
    def _sva_min_bits():
        return CONFIG['annotate'].get('sva_min_bits', 5)

    @staticmethod
    def _best_bits(dfams, prefix):
        b = [m.bits for m in dfams if m.model.startswith(prefix)]
        return max(b) if b else None

    def _sva_conclusion(self):
        """SVA-specific classifier, run before the Alu/L1 logic. SVA's Alu-like region makes a
        real SVA score an Alu model (so it was mislabelled ALU), and its poly-A tail can sit on
        either junction regardless of which strand the SVA HMM hits (so the strand-gated Alu
        acceptance rejected it) -- and when the SVA hit is on the LEFT clip with a bare poly-A
        on the right, the Alu/L1 block (nested under `len(right_dfams)>0`) never even ran. This
        method resolves the composite structure directly.

        A call needs an SVA *identity* signal and a TPRT/composite *hallmark*:
          identity: an SVA Dfam hit that dominates any Alu on its clip (sva_bits >= sva_min_bits
                    and >= that clip's best Alu), OR a CCCTCT/AGAGGG hexamer tandem (Alu/L1 lack it).
          hallmark: a poly-A tail on either clip, OR an SVA hit on BOTH clips (the element spans
                    both breakpoints), OR -- corroborating a dominant SVA hit -- a lone hexamer.
        Returns the conclusion string, or None to fall through to the generic logic. The
        homology-only paths are gated by _strict_hallmark() exactly like the Alu/L1 branches."""
        sva_L = self._best_bits(self.left_dfams, 'SVA')
        sva_R = self._best_bits(self.right_dfams, 'SVA')
        alu_L = self._best_bits(self.left_dfams, 'Alu')
        alu_R = self._best_bits(self.right_dfams, 'Alu')
        floor = self._sva_min_bits()
        dom_side = None
        for sv, al, side in ((sva_L, alu_L, 'left'), (sva_R, alu_R, 'right')):
            if sv is not None and sv >= floor and (al is None or sv >= al):
                dom_side = side
                break
        sva_dominant = dom_side is not None
        sva_both = sva_L is not None and sva_L > 0 and sva_R is not None and sva_R > 0
        hex_tandem = bool(self._SVA_HEX_TANDEM.search(self.left_seq)
                          or self._SVA_HEX_TANDEM.search(self.right_seq))
        hex_any = bool(self._SVA_HEX_ANY.search(self.left_seq)
                       or self._SVA_HEX_ANY.search(self.right_seq))
        polyA = bool(self.has_left_polyA() or self.has_right_polyA())
        hallmark = polyA or sva_both
        if self._strict_hallmark():
            # strict: the poly-A tail itself must be present; homology (dominance/hexamer) alone
            # is not accepted, mirroring the Alu/L1 strict paths.
            fire = polyA and (hex_tandem or sva_dominant)
        else:
            fire = (hallmark and (hex_tandem or sva_dominant)) or (sva_dominant and hex_any)
        if not fire:
            return None
        best = max([b for b in (sva_L, sva_R) if b is not None], default=0.0)
        ev = []
        if sva_dominant:
            ev.append(f"{dom_side} SVA {best:.0f}b")
        if hex_tandem:
            ev.append("hexamer")
        if polyA:
            ev.append("polyA")
        elif sva_both:
            ev.append("both junctions")
        DEBUG and print(f"    - SVA composite call ({', '.join(ev)}) -> accepted")
        return f"SVA (composite: {', '.join(ev)})"

    @staticmethod
    def _strict_hallmark():
        """CONFIG['annotate']['require_polya_hallmark'] (default False): when set, an RTE family
        call needs the poly-A hallmark itself -- the homology-only fallbacks are skipped, since
        sequence homology alone cannot tell a pasted-in element fragment from a genuine TPRT
        insertion. conclusion() otherwise accepts on three conditions: (1) a poly-A on the other
        side, (2) the same/an L1 dfam model on the other side, (3) an rmsk annotation of the same
        family near the other side's mapping -- (2) and (3) accept a junction carrying no TPRT
        hallmark at all. The multi-exon processed-pseudogene path is deliberately NOT gated:
        intron-skipping is its own retrotransposition hallmark, not mere homology.

        Measured on the fp10k FP-stress harness (10k true insertions + hallmark-free decoys):
            gate  polyA_min   true called   FP leaked
            off       6          92.4%        1349
            off      12          92.4%        1348
            ON        6          89.8%         286
            ON       12          89.7%         164
        i.e. the gate is the lever (-79% FP), and polya_min_len only pays off behind it.
        Default off: the -2.7pt recall cost is measured against decoys that differ from true
        insertions ONLY by the missing hallmark, which is a deliberately harsh case -- validate
        on your own data before enabling."""
        return bool(CONFIG['annotate'].get('require_polya_hallmark', False))

    # RepeatMasker classes that are retrotransposons (the only things retrotransposition
    # produces). A clip mapping to any of these is RTE-consistent; a clip mapping only to
    # unique / non-RTE sequence is not.
    _RTE_REPCLASSES = ("LINE", "SINE", "LTR", "Retroposon")

    def _pseudogene(self):
        """Processed-pseudogene signature: the inserted mRNA's clips map into exon(s) of a
        single gene (L1 machinery retrotransposes a spliced, poly-adenylated transcript).
        Returns (gene_id, n_distinct_exons, has_polyA) or None. Two distinct exons of one
        gene (introns skipped) is the strong splice signal; a single exon needs a poly-A
        tail to be called (else a lone exon overlap is not specific). Empty (no call) unless
        an exon track is configured. The multi-exon *mate-splice* test stays in discovery
        (Feature B / splice_hallmark); this is the annotate-stage, clip-level complement."""
        genes = {}
        for g, s, e in self.left_exons + self.right_exons:
            genes.setdefault(g, set()).add((s, e))
        if not genes:
            return None
        g = max(genes, key=lambda k: len(genes[k]))
        nexon = len(genes[g])
        has_polya = bool(self.has_left_polyA() or self.has_right_polyA())
        if nexon >= 2:
            return (g, nexon, has_polya)
        if nexon == 1 and has_polya:
            return (g, 1, True)
        return None

    def _maps_uniquely_to_nonrte(self, maps, min_mapq=30):
        """True if a clip maps uniquely (MAPQ >= min_mapq) to a locus carrying no
        retrotransposon annotation — i.e. a unique genomic partner locus, the signature of a
        rearrangement (e.g. a chromosomal translocation) rather than a dispersed RTE. A
        single-copy RTE maps to its RTE-annotated source, so it does NOT trip this."""
        for pos, qual, rmsks, strand in maps:
            if qual is None or qual < min_mapq:
                continue
            if not any(getattr(r, 'repClass', None) in self._RTE_REPCLASSES for r in (rmsks or [])):
                return True
        return False

    # ------------------------------------------------------------------ structural variants
    # A rearrangement (inversion, translocation, large deletion/duplication) presents to
    # PEAR-TREE as a "junction" whose clipped side, instead of being an inserted mobile
    # element, is *reference sequence from the partner breakpoint*. It is therefore an
    # annotation ORTHOGONAL to the element identity: the junction may fall inside a
    # retrotransposon, or the rearrangement may itself be MEI-caused, so an SV note is
    # *added to* (never substituted for) the Alu/L1/SVA call. See conclusion().
    #
    # The locus itself is encoded in the title as "contig:start-end". Comparing it to where
    # a clip uniquely remaps gives the subtype:
    #   * other contig                       -> translocation junction
    #   * same contig, distal, reverse clip  -> inversion (the inverted segment reads
    #                                           reverse-complemented, so its clip aligns '-';
    #                                           a co-linear deletion partner would align '+')
    #   * same contig, distal, forward clip  -> intrachromosomal SV (deletion/duplication)
    _LOCUS_RE = re.compile(r'^(.+):(\d+)-(\d+)$')
    _MAP_RE = re.compile(r'^(.+):(\d+)[+-]?$')

    @staticmethod
    def _fmt_dist(d):
        if d >= 1_000_000:
            return f"{d/1e6:.1f} Mb"
        return f"{d/1000:.0f} kb" if d >= 1000 else f"{d} bp"

    def _parse_locus(self):
        m = self._LOCUS_RE.match(self.title)
        return (m.group(1), int(m.group(2)), int(m.group(3))) if m else None

    def _sv_partner_maps(self, allow_rte):
        """Clip remaps that can be trusted as a rearrangement partner: MAPQ >= sv_min_mapq
        (unique — a dispersed-repeat clip multi-maps at MAPQ 0 and is excluded, which is what
        keeps the whole Alu family from reading as translocations). allow_rte=False also drops
        partners carrying any RTE annotation (the conservative pure-SV gate, unchanged from the
        legacy flag); allow_rte=True keeps them, so a repeat-mediated / MEI-caused SV whose
        partner happens to sit in a young uniquely-mapping element is still seen.
        Yields (side, partner_contig, partner_pos, strand, partner_is_rte)."""
        floor = CONFIG['annotate'].get('sv_min_mapq', 30)
        out = []
        for side, maps in (('left', self.left_maps), ('right', self.right_maps)):
            # (D) clip-map trust guard: a short / AT-rich / low-entropy clip yields chance or
            # paralogous 'unique' hits that masquerade as translocation partners. Drop this
            # side's maps entirely if the clip that produced them is not credibly unique.
            if not self._clip_trustworthy_for_sv(side):
                continue
            for pos, qual, rmsks, strand in maps:
                if qual is None or qual < floor:
                    continue
                is_rte = any(getattr(r, 'repClass', None) in self._RTE_REPCLASSES
                             for r in (rmsks or []))
                if not allow_rte and is_rte:
                    continue
                m = self._MAP_RE.match(pos)
                if not m:
                    continue
                out.append((side, m.group(1), int(m.group(2)), strand, is_rte))
        return out

    def _sv_subtype(self, allow_rte, inversion_only=False):
        """Best SV descriptor for this junction, or None. Returns a tuple
        (rank, description, partner_contig, partner_pos); lower rank wins (inversion 0 <
        translocation 1 < deletion/duplication 2). `inversion_only` keeps only the inversion
        signature (used as the transduction guard on poly-A-bearing RTE calls)."""
        loc = self._parse_locus()
        if loc is None:
            return None
        contig, s, e = loc
        mind = CONFIG['annotate'].get('sv_min_distance', 1000)
        best = None
        for side, pc, pp, ps, is_rte in self._sv_partner_maps(allow_rte):
            if pc == contig:
                dist = min(abs(pp - s), abs(pp - e))
                if dist < mind:
                    continue                       # local micro-context, not a rearrangement
                if ps == '-':
                    rank = 0
                    desc = f"inversion junction -> {pc}:{pp} ({self._fmt_dist(dist)})"
                elif inversion_only:
                    continue
                else:
                    rank = 2
                    desc = (f"intrachromosomal SV (deletion/duplication) -> {pc}:{pp} "
                            f"({self._fmt_dist(dist)})")
            elif inversion_only:
                continue
            else:
                rank = 1
                if self.reciprocal_partner is not None:
                    desc = f"balanced translocation (reciprocal {self.reciprocal_partner})"
                else:
                    desc = f"translocation junction -> {pc}:{pp}"
            if is_rte:
                desc += " (repeat-mediated)"
            if best is None or rank < best[0]:
                best = (rank, desc, pc, pp)
        return best

    def get_fasta(self) -> str:
        """
        This function returns a FASTA chunk with the inserted sequences for the insertion.
        """
        # one-sided records (poly-A / discordant-anchored end) carry an empty junction string;
        # the min() over an empty list used to crash here.
        self.left_ins_seq, self.right_ins_seq = (s.upper() for s in self._insert_clips())
        return f'>{self.title}:R\n{self.right_ins_seq}\n>{self.title}:L\n{self.left_ins_seq}\n'

    # ------------------------------------------------------ uncharacterised complex insertion
    # A junction with no Dfam hit, no clip remap and no poly-A used to collapse -- whatever its
    # clips contained -- into 'artefact', the same bin as genotyping noise and poly-A slippage.
    # But a junction whose two breakpoints carry SUBSTANTIAL, HIGH-COMPLEXITY inserted sequence
    # that simply matched nothing is not noise: it is a real but uncharacterised insertion -- a
    # non-MEI / complex / templated insertion whose inserted sequence is novel or chimeric, so
    # bowtie2's --end-to-end remap (which must align the WHOLE clip to one reference block)
    # drops it and no retrotransposon model scores it. Separate the two so these do not
    # masquerade as artefacts. Thresholds are deployment-local (defaulted in code).
    @staticmethod
    def _shannon(seq):
        """Per-base Shannon entropy (bits) of a sequence; 0 for empty/homopolymer, ~2 for a
        balanced 4-letter mix. Distinguishes complex inserted sequence from low-complexity
        homopolymer/STR slippage."""
        seq = seq.upper()
        if not seq:
            return 0.0
        counts = {}
        for b in seq:
            counts[b] = counts.get(b, 0) + 1
        n = len(seq)
        return -sum((c / n) * log2(c / n) for c in counts.values())

    def _insert_clips(self):
        """The inserted (clipped, lower-case) sequence on each junction -- the left junction's
        leading run before the first reference (upper-case) base and the right junction's
        trailing run from the first inserted (lower-case) base. These are exactly the
        substrings get_fasta() submits to Dfam / bowtie2. Returns (left_insert, right_insert)."""
        ls, rs = self.left_seq, self.right_seq
        li_end = next((i for i, ch in enumerate(ls) if ch.isupper()), len(ls))
        ri_start = next((i for i, ch in enumerate(rs) if ch.islower()), len(rs))
        return ls[:li_end], rs[ri_start:]

    def _is_complex_insertion(self):
        """True when BOTH breakpoints carry substantial (>= complex_ins_min_len bp),
        high-complexity (>= complex_ins_min_entropy bits on at least one side) inserted
        sequence -- the signature of a real but uncharacterised insertion rather than noise."""
        min_len = CONFIG['annotate'].get('complex_ins_min_len', 20)
        min_ent = CONFIG['annotate'].get('complex_ins_min_entropy', 1.6)
        left_ins, right_ins = self._insert_clips()
        if len(left_ins) < min_len or len(right_ins) < min_len:
            return False
        return max(self._shannon(left_ins), self._shannon(right_ins)) >= min_ent

    # ------------------------------------------------ reference-free local-rearrangement subtypes
    # Ten worker analyses of TP "unknown" loci (2026-07-18) found they are overwhelmingly REAL,
    # non-MEI local events whose clips are NOT novel sequence but *copies of this locus's own
    # reference flanks*, joined across the breakpoint. Because such a clip is chimeric (part = one
    # flank, part = the other / a novel seam) it fails bowtie2 --end-to-end by construction, and it
    # carries no Dfam/poly-A signal -- so it fell into "unknown". They are cheaply nameable from
    # the two junction strings ALONE (no reference lookup): the inserted (lower-case) clip of one
    # junction reappears as the aligned (upper-case) flank of the other. Forward -> tandem /
    # segmental duplication; reverse-complement -> inverted duplication; a short-period tandem
    # repeat flush to both breakpoints -> microsatellite length change.
    _STR_RUN = re.compile(r'(([ACGT]{1,3})\2{4,})')

    @staticmethod
    def _rc(seq):
        return seq.upper().translate(str.maketrans('ACGT', 'TGCA'))[::-1]

    @staticmethod
    def _longest_submatch(a, b):
        """Length of the longest exact substring of `a` that occurs in `b` (both upper-cased).
        O(len(a)*best); fine for the short (<300 bp) clips/flanks here."""
        a = a.upper(); b = b.upper()
        best = 0
        la = len(a)
        for i in range(la):
            k = best + 1
            while i + k <= la and a[i:i + k] in b:
                k += 1
            if k - 1 > best:
                best = k - 1
        return best

    def _flank_uppers(self):
        """The reference-aligned (upper-case) genomic flank of each junction."""
        lf = ''.join(c for c in self.left_seq if c.isupper())
        rf = ''.join(c for c in self.right_seq if c.isupper())
        return lf, rf

    @staticmethod
    def _norm_unit(unit):
        """Rotation-invariant canonical form of a repeat unit (its lexicographically minimal
        rotation), so (CA)n and (AC)n compare equal."""
        return min(unit[i:] + unit[:i] for i in range(len(unit))) if unit else unit

    def _microsatellite_subtype(self, li, ri, lf, rf):
        """A polymorphic microsatellite / STR length change (not an insertion): the SAME tandem
        repeat unit (rotation-invariant), long and flush to BOTH breakpoints, with each soft-clip
        also sharing a stretch with the opposite genomic flank (same locus). Returns a label or
        None. Requiring the *same unit* at both breakpoints -- not merely a repeat near each -- is
        what stops an incidental short repeat in a complex insertion from being mis-called STR (so
        those stay on the duplication / complex-insertion path). Runs before the (entropy-gated)
        duplication test so a low-complexity (CA)n tract is called STR, not a segmental dup."""
        ls, rs = self.left_seq.upper(), self.right_seq.upper()
        lb = len(li)                 # left  breakpoint = insert|flank boundary
        rb = len(rs) - len(ri)       # right breakpoint = flank|insert boundary
        tol = 3   # the repeat need only be flush to the breakpoint, not straddle it
        min_run = CONFIG['annotate'].get('str_min_run', 14)
        def crossing_unit(seq, bnd):
            best = None
            for m in self._STR_RUN.finditer(seq):
                if m.start() <= bnd + tol and m.end() >= bnd - tol and (m.end() - m.start()) >= min_run:
                    run = m.group(1)
                    p = next((q for q in (1, 2, 3)
                              if run == (run[:q] * (len(run) // q + 1))[:len(run)]), len(m.group(2)))
                    cand = (m.end() - m.start(), self._norm_unit(run[:p]))
                    if best is None or cand[0] > best[0]:
                        best = cand
            return best
        uL = crossing_unit(ls, lb)
        uR = crossing_unit(rs, rb)
        if uL is None or uR is None or uL[1] != uR[1]:   # same canonical repeat unit both sides
            return None
        share = (max(self._longest_submatch(li, rf), self._longest_submatch(li, lf)) >= 10
                 and max(self._longest_submatch(ri, lf), self._longest_submatch(ri, rf)) >= 10)
        if not share:
            return None
        return f"microsatellite ({uL[1]})n length change"

    def _reciprocal_dup_subtype(self):
        """Reference-free local-duplication / STR subtype for a junction the element and SV logic
        could not explain, or None. Fires only on the *same-locus* signature (clip = copy of the
        opposite flank), so it is safe to consult only in the unknown/artefact fallbacks -- an
        accepted Alu/L1/SVA call never reaches here, so its (TSD-driven) flank match is moot."""
        li, ri = (s.upper() for s in self._insert_clips())
        if not li or not ri:
            return None
        lf, rf = (s.upper() for s in self._flank_uppers())
        mind = CONFIG['annotate'].get('dup_min_match', 18)
        ment = CONFIG['annotate'].get('dup_min_entropy', 1.7)
        # microsatellite first (low-entropy, own detector)
        micro = self._microsatellite_subtype(li, ri, lf, rf)
        if micro is not None:
            return micro
        # tandem / segmental duplication: BOTH inserts carry a >= mind exact copy of a flank,
        # forward orientation, and both inserts are complex (entropy gate excludes STR/homopolymer).
        if min(self._shannon(li), self._shannon(ri)) >= ment:
            mL = max(self._longest_submatch(li, rf), self._longest_submatch(li, lf))
            mR = max(self._longest_submatch(ri, lf), self._longest_submatch(ri, rf))
            if mL >= mind and mR >= mind:
                return f"tandem/segmental duplication (dup unit >= {min(mL, mR)} bp)"
            mLrc = max(self._longest_submatch(li, self._rc(rf)), self._longest_submatch(li, self._rc(lf)))
            mRrc = max(self._longest_submatch(ri, self._rc(lf)), self._longest_submatch(ri, self._rc(rf)))
            if mLrc >= mind and mRrc >= mind:
                return f"inverted duplication (>= {min(mLrc, mRrc)} bp)"
        return None

    # ------------------------------------------------------------------- flank-leak Dfam demotion
    def _dfam_is_flank_leak(self, dfam, side):
        """True when a Dfam hit sits on clip sequence that is actually genomic FLANK that leaked
        into the soft-clip (a short insert forces the aligner to clip contiguous reference past
        the breakpoint, and that reference tail can carry an old-repeat HMM hit describing the
        SITE, not the inserted element). Detected by: the clip sub-sequence under the hit is a
        long exact copy of one of this locus's reference flanks. Purely local; no reference lookup."""
        need = CONFIG['annotate'].get('flankleak_min', 30)
        li, ri = (s.upper() for s in self._insert_clips())
        clip = li if side == 'left' else ri
        s = max(0, dfam.ali_start - 1)
        e = min(len(clip), dfam.ali_end)
        sub = clip[s:e]
        if len(sub) < need:
            return False
        lf, rf = (s.upper() for s in self._flank_uppers())
        want = min(len(sub), need)
        return self._longest_submatch(sub, lf) >= want or self._longest_submatch(sub, rf) >= want

    # --------------------------------------------------------------- clip-map trust guard (for SV)
    def _clip_trustworthy_for_sv(self, side):
        """A clip remap is only a credible rearrangement partner if the clip that produced it is
        long enough and complex enough to map uniquely on merit. A short, AT-rich / low-entropy
        clip yields chance / paralogous 'unique' hits that masquerade as translocation partners
        (seen repeatedly in the worker analyses). Gate them out of the SV partner set."""
        li, ri = (s.upper() for s in self._insert_clips())
        clip = li if side == 'left' else ri
        if len(clip) < CONFIG['annotate'].get('sv_clip_min_len', 25):
            return False
        if self._shannon(clip) < CONFIG['annotate'].get('sv_clip_min_entropy', 1.9):
            return False
        at = (clip.count('A') + clip.count('T')) / len(clip)
        if at > CONFIG['annotate'].get('sv_clip_max_at', 0.72):
            return False
        return True

    # ---------------------------------------------------------- (E) split / local-remap subtype
    def _local_remap_subtype(self):
        """Feature E: when the whole-clip --end-to-end remap failed but a bowtie2 --local pass
        placed part of a clip, use that partial placement to name the event. A clip whose novel
        core maps uniquely to a DISTAL locus is a templated / complex insertion sourced there; a
        clip that maps back within this locus corroborates a local duplication. Returns a label or
        None. Inert unless the local-remap channel ran (left_local_maps / right_local_maps
        populated by read_sam_local); it is a no-op on data scored without the --local pass."""
        loc = self._parse_locus()
        if loc is None:
            return None
        contig, s, e = loc
        # bowtie2 --local caps MAPQ lower than --end-to-end (the soft-clipped read is shorter),
        # so a unique local placement scores ~20-24; use a channel-specific floor, not sv_min_mapq.
        floor = CONFIG['annotate'].get('local_min_mapq', 20)
        near = CONFIG['annotate'].get('sv_min_distance', 1000)
        best = None
        for side, maps in (('left', self.left_local_maps), ('right', self.right_local_maps)):
            for pc, pp, strand, qcov, mapq in maps:
                if mapq is None or mapq < floor:
                    continue
                if qcov < CONFIG['annotate'].get('local_min_qcov', 20):
                    continue
                if pc == contig and min(abs(pp - s), abs(pp - e)) < near:
                    cand = (2, "local duplication (confirmed by split remap)")
                else:
                    cand = (1, f"templated/complex insertion (source {pc}:{pp})")
                if best is None or cand[0] < best[0]:
                    best = cand
        return best[1] if best else None

    def site(self):
        """Insertion-site annotation for this junction's own reference locus (the title), or
        None when no gene model is configured or the title carries no parseable locus. Returns
        the GeneModel.annotate() tuple (region_keyword, gene, strand, human_text). Orthogonal to
        the element identity: it says WHERE the insertion landed, not WHAT landed."""
        if self.gene_model is None:
            return None
        loc = self._parse_locus()
        if loc is None:
            return None
        contig, s, e = loc
        return self.gene_model.annotate(contig, s, e)

    def conclusion(self) -> str:
        """Full annotation: the element identity from _element_conclusion(), plus — because a
        structural variant and a mobile element are NOT mutually exclusive (a rearrangement
        junction may fall inside a retrotransposon, or the SV may itself be MEI-caused) — an
        additive `[SV: ...]` note when a clip uniquely maps to a distal / other-contig
        rearrangement partner, and — orthogonally to both — an additive `[site: ...]` note
        describing where the insertion LANDED (gene / exon / splice site / intron / promoter),
        from GeneModel via the title locus. Neither note ever changes the element class and
        nothing is dropped; they only annotate. Pure structural variants (no element call) are
        handled inside _element_conclusion() and carry the subtype directly."""
        out = self._conclusion_no_site()
        site = self.site()
        if site is not None:
            out = f"{out} [site: {site[3]}]"
        return out

    def _conclusion_no_site(self) -> str:
        base = self._element_conclusion()
        cls = VariantAnnotationContainer.element_class(base)
        if cls in ('ALU', 'LINE1', 'SVA', 'RTE_other'):
            # An element call gets an additive SV note only for a *trustworthy* partner: a clip
            # that maps uniquely to NON-RTE sequence (the genuine MEI-caused-SV signal — e.g. an
            # Alu-mediated deletion whose far breakpoint is a unique locus). A clip matching only
            # a paralogous repeat copy is NOT reported: it is far more likely the element itself
            # aligning to a family member than a real rearrangement, and would decorate hundreds
            # of ordinary Alu/L1 insertions with spurious translocations. A repeat partner is
            # admitted only when reciprocally confirmed (two junctions pointing back at each
            # other — corroboration a lone paralog match cannot fake).
            allow_rte = self.reciprocal_partner is not None
            # Transduction guard: a 3' transduction drags a unique, poly-adenylated genomic
            # flank that mimics a co-oriented SV partner. Under a poly-A hallmark accept only
            # the inversion signature (a reverse-strand partner), which a transduction — always
            # co-oriented — cannot produce.
            polya = bool(self.has_left_polyA() or self.has_right_polyA())
            sub = self._sv_subtype(allow_rte=allow_rte, inversion_only=polya)
            if sub is not None:
                return f"{base} [SV: {sub[1]}]"
        return base

    def _element_conclusion(self) -> str:
        """
        This function aggregates all information available to come to a conclusion
        """
        # Discovery splice-hallmark (Feature B): mate reads span >= 2 exons of a single gene
        # with the introns skipped — the definitive processed-pseudogene signal, which the
        # clip-level exon check below cannot see (annotate's clips are the terminal 5' exon +
        # 3' poly-A). Reported first, even when the terminal clips did not map.
        if self.splice_hits:
            gene, side, nex = max(self.splice_hits, key=lambda x: x[2])
            return f"processed pseudogene of {gene} (splice: {nex} exons, discovery mates)"
        if len(self.right_dfams)==0 and len(self.left_dfams)==0 and len(self.right_maps)==0 and len(self.left_maps)==0:
            # Nothing identified this junction: no Dfam hit, no clip remap, no poly-A. If BOTH
            # clips nonetheless carry substantial high-complexity inserted sequence it is a real
            # but uncharacterised insertion (a non-MEI / complex / templated insertion whose novel
            # or chimeric clip fails the --end-to-end remap and matches no retrotransposon model),
            # not genotyping noise -- keep those out of the artefact bin.
            # First try to NAME it: most such loci are local tandem/segmental duplications or
            # microsatellite length changes whose clips are copies of this locus's own flanks
            # (reference-free signature), or -- if the --local pass ran -- a templated insertion.
            dup = self._reciprocal_dup_subtype() or self._local_remap_subtype()
            if dup is not None:
                return dup
            if self._is_complex_insertion():
                return 'unknown (unmapped complex insertion)'
            return 'artefact'
        # SVA is a composite element: its Alu-like body scores an Alu model and its poly-A can be
        # on either junction, so the strand-gated Alu/L1 logic below both mislabels and (when the
        # SVA hit is on the left clip with a bare poly-A right, since that block is nested under
        # `len(right_dfams)>0`) skips real SVAs. Give the SVA-specific composite test first claim.
        sva = self._sva_conclusion()
        if sva is not None:
            return sva
        if len(self.right_dfams)>0:
            for m in self.right_dfams:
                if True or m.is_active:
                    if "L1" in m.model:
                        DEBUG and print(f"    - detected RIGHT L1 dfam model {m.model}...")
                        # ignore strand for LINEs, since these can sometimes reverse a part of the LINE at the 5' insertion site.
                        # accept this in three conditions
                        # 1) polyA on the other side
                        if self.has_left_polyA():
                            DEBUG and print(f"      has left polyA -> accepted")
                            return f"polyA < {m}"
                        else:
                            DEBUG and print(f"      does not have left polyA -> rejected, checking left dfams")
                        if self._strict_hallmark():
                            DEBUG and print(f"      strict hallmark: no left polyA -> no call")
                        else:
                            # 2) L1 element on the other side
                            for m2 in self.left_dfams:
                                if "L1" in m2.model:
                                    DEBUG and print(f"      found another L1 ({m2.model}) in left defams -> accepted.")
                                    return f"{m2} <- {m}"
                            DEBUG and print(f"      does not have a suitable left dfam, checking mappings.")
                            # 3) match near an L1 element on the other side.
                            for pos, qual, rmsks, strand in self.left_maps:
                                for r in rmsks:
                                    if r.repFamily == "L1":
                                        DEBUG and print(f"      found a suitable rmsk annotation ({r}) in left mapping -> accepted.")
                                        return f"{r}({strand}) <- {m}"
                    else:
                        DEBUG and print(f"    - detected RIGHT dfam model {m.model} on strand {m.strand}, checking...")
                        # dont ignore strands for Alus and since these should have the same orientation
                        # accept this in three conditions
                        # 1) the other side has a polyA, but only if the RTE is on reverse
                        if m.strand == "+" and self.has_left_polyA():
                            DEBUG and print(
                                f"      found a suitable polyA in the left mapping -> accepted.")
                            return f"polyA <- {m}"
                        if self._strict_hallmark():
                            DEBUG and print(f"      strict hallmark: no left polyA -> no call")
                        else:
                            # 2) same element on the other side, oriented in the same direction
                            for m2 in self.left_dfams:
                                if m2.is_active and m2.model[:3] == m.model[:3] and m2.strand == m.strand:
                                    DEBUG and print(
                                        f"      found a suitable model in the left dfams ({m2.model} on strand {m2.strand}) in the left mapping -> accepted.")
                                    return f"{m2} <- {m}"
                            # 3) match near an L1 element on the other side.
                            for pos, qual, rmsks, strand in self.left_maps:
                                for r in rmsks:
                                    if r.repName[:3] == m.model[:3]:
                                        # check same strand
                                        if (strand == "+") ^ (r.strand == "+") == (m.strand == "-"):
                                            DEBUG and print(
                                                f"      found a suitable mapping in the left dfams ({r} on strand {r.strand}, mapping is on strand {strand})-> accepted.")
                                            return f"{r}({strand}) <- {m}"
            for m in self.left_dfams:
                if True or m.is_active:
                    if "L1" in m.model:
                        DEBUG and print(f"    - detected left L1 dfam model {m.model}...")
                        # ignore strand for LINEs, since these can sometimes reverse a part of the LINE at the 5' insertion site.
                        # accept this in three conditions
                        # 1) polyA on the other side
                        if self.has_right_polyA():
                            DEBUG and print(f"      has right polyA -> accepted")
                            return f"{m} -> polyA"
                        # 2) L1 element on the other side
                        # ignore, already covered above.
                        if self._strict_hallmark():
                            DEBUG and print(f"      strict hallmark: no right polyA -> no call")
                        else:
                            # 3) match near an L1 element on the other side.
                            for pos, qual, rmsks, strand in self.right_maps:
                                for r in rmsks:
                                    if r.repFamily == "L1":
                                        DEBUG and print(
                                            f"      found a suitable rmsk annotation ({r}) in right mapping -> accepted.")
                                        return f"{m} -> {r}({strand})"
                    else:
                        # dont ignore strands for Alus and since these should have the same orientation
                        # accept this in three conditions
                        # 1) the other side has a polyA, but only if the RTE is on reverse
                        DEBUG and print(f"    - detected left dfam model {m.model} on strand {m.strand}, checking...")
                        if m.strand == "-" and self.has_right_polyA():
                            return f"{m} -> polyA"
                        # 2) same element on the other side, oriented in the same direction
                        # ignore, already covered above.
                        if self._strict_hallmark():
                            DEBUG and print(f"      strict hallmark: no right polyA -> no call")
                        else:
                            # 3) match near an L1 element on the other side.
                            for pos, qual, rmsks, strand in self.right_maps:
                                for r in rmsks:
                                    if r.repName[:3] == m.model[:3]:
                                        # check same strand
                                        if (strand == "+") ^ (r.strand == "+") == (m.strand == "-"):
                                            return f"{m} -> {r}({strand})"
        # check elements only mapped, but not identified by dfam
        l1_map_check = False
        for pos, qual, rmsks, strand in self.left_maps:
            for r in rmsks:
                if 'L1' in r.repName[:2] and not l1_map_check:
                    if self.has_right_polyA(): return r.repName
                    l1_map_check = True #dont check this twice
                    if self._strict_hallmark():
                        continue          # homology on both sides is not a TPRT hallmark
                    for pos2, qual2, rmsks2, strand2 in self.right_maps:
                        for r2 in rmsks2:
                            if 'L1' in r2.repName[:2]:
                                return r.repName
        if not l1_map_check:
            # check left side for L1 match and right for polyA
            for pos, qual, rmsks, strand in self.right_maps:
                if l1_map_check: break
                for r in rmsks:
                    if 'L1' in r.repName[:2] and not l1_map_check:
                        if self.has_left_polyA(): return r.repName
                        l1_map_check = True
                        break
        plus_elements = set()
        minus_elements = set()
        for pos, qual, rmsks, strand in self.left_maps:
            for r in rmsks:
                if (strand == "+") ^ (r.strand == "+"):
                    if self.has_right_polyA(): return r.repName
                    minus_elements.add(r.repName[:3])
                else:
                    plus_elements.add(r.repName[:3])
        for pos, qual, rmsks, strand in self.right_maps:
            for r in rmsks:
                if (strand == "+") ^ (r.strand == "+"):
                    if not self._strict_hallmark() and r.repName[:3] in minus_elements:
                        return r.repName
                else:
                    if self.has_left_polyA(): return r.repName
                    if not self._strict_hallmark() and r.repName[:3] in plus_elements:
                        return r.repName
        # --- processed-pseudogene annotation (non-gating) ---
        # A retrotransposed spliced mRNA: not an RTE (no Dfam/rmsk RTE hit above), but its
        # clips map into exon(s) of a gene, with a poly-A tail. Checked before the SV flag
        # because a pseudogene clip is a unique non-RTE map too — this is the more specific
        # explanation.
        pg = self._pseudogene()
        if pg and CONFIG['annotate'].get('pseudogene_require_exon_junction') \
                and self.exon_junction_proven is False:
            # tools/rte looked for an exon-exon junction read and found none: candidate only
            DEBUG and print(f"    - pseudogene of {pg[0]} lacks an exon-exon junction read -> not called")
            pg = None
        if pg:
            gene, nexon, has_polya = pg
            detail = f"{nexon} exons" if nexon >= 2 else "1 exon + polyA"
            return f"processed pseudogene of {gene} ({detail})"
        # --- reference-free local duplication / microsatellite (precedes the distal-SV flag) ---
        # A clip that is a copy of THIS locus's own flank is a local duplication, not a distal
        # rearrangement — so this test takes precedence over the translocation flag below, which
        # otherwise fires on a paralogous / spurious cross-locus map (worker analyses 7/9/10).
        dup = self._reciprocal_dup_subtype()
        if dup is not None:
            return dup
        # --- annotation-only non-RTE structural-variant flag (non-gating; nothing removed) ---
        # This call could not be explained as a retrotransposition above. If NEITHER junction
        # carries a poly-A tail AND a clip maps uniquely to a locus with no retrotransposon
        # annotation, it looks like a rearrangement partner (a chromosomal translocation
        # joining two loci), not an insertion. It is only FLAGGED for review, never dropped:
        # the cross-locus pointer alone cannot separate a translocation from a single-copy RTE
        # or a transduction, so we lean on the RTE hallmarks instead — a single-copy RTE maps
        # to its RTE-annotated source and a transduction keeps its poly-A, so both are resolved
        # above and never reach here.
        # (D) the SV flag only fires for a partner from a *trustworthy* clip (long / complex /
        # not AT-rich); a short low-complexity clip's 'unique' map is a chance/paralogous hit, not
        # a rearrangement partner, so it must not even raise the generic flag.
        if (not self.has_left_polyA() and not self.has_right_polyA()
                and ((self._clip_trustworthy_for_sv('left')
                      and self._maps_uniquely_to_nonrte(self.left_maps))
                     or (self._clip_trustworthy_for_sv('right')
                         and self._maps_uniquely_to_nonrte(self.right_maps)))):
            # A pure (non-RTE) rearrangement: subtype it (translocation / inversion /
            # deletion-duplication) when the partner resolves, else keep the legacy generic
            # flag. Same firing condition as before, so no locus that used to be flagged is
            # lost — the change is only that the label is now specific where it can be.
            sub = self._sv_subtype(allow_rte=False)
            if sub is not None:
                return sub[1]
            return 'unknown (possible non-RTE SV, e.g. translocation)'
        # (E) last resort: a split/--local placement of a chimeric clip, if that pass ran.
        loc_sub = self._local_remap_subtype()
        if loc_sub is not None:
            return loc_sub
        return 'unknown'


class VariantAnnotationContainer:
    def __init__(self, sample, output):
        self.sample = sample
        self.insertions_file = CONFIG['annotate']['insertions_file'](sample)
        self.genotyping_file = CONFIG['annotate']['genotyping_file'](sample)
        self.dfam_file = CONFIG['annotate']['tmp']('dfam')(sample)
        self.sam_file = CONFIG['annotate']['tmp']('sam')(sample)
        self.local_sam_file = CONFIG['annotate']['tmp']('local.sam')(sample)
        self.fasta_file = CONFIG['annotate']['tmp']('fa.gz')(sample)
        self.output = output
        self.insertions = {}
        if not os.path.exists(self.insertions_file):
            raise FileNotFoundError(f"Insertions file {self.insertions_file} not found")
        self.read_insertions()
        self.read_splice()
        if not os.path.exists(self.genotyping_file):
            raise FileNotFoundError(f"Genotyping file {self.genotyping_file} not found")
        self.read_genotyping()
        if not os.path.exists(self.fasta_file) or os.path.getsize(self.fasta_file) == 0:
            self.generate_fasta_file()
        if not os.path.exists(self.dfam_file) or os.path.getsize(self.dfam_file) == 0:
            self.generate_dfam_file()
        self.read_dfam()
        # (C) drop Dfam hits that sit on genomic flank leaked into the soft-clip, so a short
        # insert's flanking old-repeat context cannot masquerade as a mobile-element hallmark.
        self.demote_flank_leak_dfams()
        if not os.path.exists(self.sam_file) or os.path.getsize(self.sam_file) == 0:
            self.generate_sam_file()
        self.read_sam()
        # (E) optional bowtie2 --local pass: places the mappable part of a chimeric clip that the
        # --end-to-end pass dropped, letting a templated/complex insertion be sourced and a local
        # duplication corroborated. Guarded — a no-op unless the aligner+index are present and
        # the channel is enabled, so --end-to-end-only runs are unaffected.
        if CONFIG['annotate'].get('sv_local_remap', True):
            if not os.path.exists(self.local_sam_file) or os.path.getsize(self.local_sam_file) == 0:
                self.generate_sam_local_file()
            if os.path.exists(self.local_sam_file) and os.path.getsize(self.local_sam_file) > 0:
                self.read_sam_local()
        self.link_reciprocal_translocations()
        self.read_gene_model()
        # TPRT-hallmark annotation (tools/rte): additive columns, only when configured
        self.rte_records = self.run_rte()

    def run_rte(self):
        """tools/rte plug-in: element / structure / tags / TSD / EN / poly-A / TPRT score per
        insertion, from the combine sidecars (insertions.evidence.tsv.gz + insertions.reads.fa.gz)
        and the reference RTE library. Enabled by CONFIG['annotate']['rte_library']; returns {}
        (no new columns, legacy output byte-identical) otherwise. The sidecars default to the
        insertions file's siblings and may be absent (the junction strings are used alone)."""
        cfg = CONFIG['annotate']
        if not cfg.get('rte_library'):
            return {}
        try:
            from tools.rte.annotator import RteAnnotator, InsertionInput, default_sidecars
        except ImportError:                       # run as `python tools/annotate_v2.py`
            from rte.annotator import RteAnnotator, InsertionInput, default_sidecars
        ev_default, rd_default = default_sidecars(self.insertions_file)
        ev_path = cfg.get('rte_evidence_file') or ev_default
        rd_path = cfg.get('rte_reads_file') or rd_default
        ev_path = ev_path(self.sample) if callable(ev_path) else ev_path
        rd_path = rd_path(self.sample) if callable(rd_path) else rd_path
        ann = RteAnnotator(cfg, gene_model=Insertion.gene_model)
        self.rte_lib = ann.lib
        ann.load_evidence(ev_path, rd_path, wanted=set(self.insertions))
        inputs = {k: InsertionInput.from_legacy(ins, self.element_class(ins.conclusion()))
                  for k, ins in self.insertions.items()}
        records = ann.annotate_all(inputs)
        for k, rec in records.items():
            ins = self.insertions.get(k)
            if ins is not None and inputs[k].pseudogene_genes:
                ins.exon_junction_proven = 'EXON_JUNCTION' in rec.tags
        print(f"[rte] annotated {len(records)} insertions "
              f"({sum(1 for r in records.values() if r.tprt_call == 'TPRT')} TPRT)")
        return records

    def read_gene_model(self):
        """Load the insertion-SITE gene model (CONFIG['annotate']['gene_model']) once and share it
        across all insertions via the Insertion.gene_model class attribute. Optional: None (the
        default) leaves site annotation off. Window parameters are read from the same config
        block (splice_*_window, promoter_*), so a deployment can tune them without code changes."""
        path = CONFIG['annotate'].get('gene_model')
        if not path:
            Insertion.gene_model = None
            return
        Insertion.gene_model = GeneModel(path, CONFIG['annotate'])

    def link_reciprocal_translocations(self, window=100000):
        """Pair reciprocal (balanced) translocation junctions. A single junction only shows a
        one-way pointer locus_A -> partner_B; it is a *balanced* translocation when another
        junction near B points back near A. This pass finds those pairs and records the partner
        locus on each Insertion (`reciprocal_partner`), which _sv_subtype() then reports as
        'balanced translocation (reciprocal ...)'. `window` is the coordinate tolerance for
        matching a junction to a partner breakpoint. Interchromosomal junctions only (few
        hundred typically), so the pairwise scan is cheap."""
        trans = []   # (key, locus_contig, locus_mid, partner_contig, partner_pos)
        for key, ins in self.insertions.items():
            loc = ins._parse_locus()
            if loc is None:
                continue
            sub = ins._sv_subtype(allow_rte=True)          # reciprocal_partner still None here
            if sub is None or sub[0] != 1:                 # rank 1 == interchromosomal
                continue
            trans.append((key, loc[0], (loc[1] + loc[2]) // 2, sub[2], sub[3]))
        n = 0
        for ki, lci, lmi, pci, ppi in trans:
            for kj, lcj, lmj, pcj, ppj in trans:
                if ki == kj:
                    continue
                # kj's locus sits at ki's partner, and kj's partner points back at ki's locus.
                if (lcj == pci and abs(lmj - ppi) <= window
                        and pcj == lci and abs(ppj - lmi) <= window):
                    self.insertions[ki].reciprocal_partner = f"{lcj}:{lmj}"
                    n += 1
                    break
        print(f"linked {n} reciprocal (balanced) translocation junction(s).")

    def read_splice(self):
        """Attach discovery splice-hallmark (Feature B) evidence to the insertions. The
        sidecar `<combined>.splice.tsv` is written by combine_insertions (which re-keys the
        per-discovery-file splice rows onto the combined insertion names). Optional: absent
        when splice_hallmark was off during discovery, so this is a no-op then."""
        path = (self.insertions_file[:-7] if self.insertions_file.endswith(".txt.gz")
                else self.insertions_file) + ".splice.tsv"
        if not os.path.exists(path):
            return
        n = 0
        with open(path) as fh:
            fh.readline()  # header: insertion gene side n_exons intron_bp span_bp
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 4:
                    continue
                name, gene, side = p[0], p[1], p[2]
                try:
                    nex = int(p[3])
                except ValueError:
                    continue
                if name in self.insertions:
                    self.insertions[name].splice_hits.append((gene, side, nex))
                    n += 1
        print(f"attached {n} discovery splice-hallmark rows from {path}")

    def read_insertions(self):
        """
        Read the insertions.combined.txt.gz FASTQ (records `@<title>:L` / `@<title>:R`) into
        Insertion objects. Records are paired by title, so one-sided insertions (a poly-A- or
        discordant-anchored end has no junction reads; combine writes only the other side) get
        an empty string for the missing side instead of crashing the L-then-R assertion.
        """
        print(f"reading insertions file {self.insertions_file}...")
        sides = {}
        order = []
        with gzip.open(self.insertions_file, 'rt') as ifh:
            while True:
                header = ifh.readline()
                if not header:
                    break
                header = header.strip()
                if not header:
                    continue
                assert header[0] == '@', f"malformed record header {header!r}"
                seq = ifh.readline().strip()
                assert ifh.readline().strip().startswith('+')
                ifh.readline()                       # quality
                title, side = header[1:-2], header[-1]
                assert side in ('L', 'R'), f"unknown junction side in {header!r}"
                if title not in sides:
                    sides[title] = {}
                    order.append(title)
                sides[title][side] = seq
        for title in order:
            self.insertions[title] = Insertion(title, sides[title].get('L', ''), sides[title].get('R', ''))
        n1 = sum(1 for t in order if len(sides[t]) == 1)
        print(f"done reading insertions file {self.insertions_file}, read {len(self.insertions)} insertions "
              f"({n1} one-sided).")

    def read_genotyping(self):
        """
        This function read sthe genotyping file, calculates the number of tips with insertions, artefacts and wild-type calls, updates the Insertion object with these numbers and filteres insertions with at least one heterozygous or homozygous call.
        """
        print(f"reading genotyping file {self.genotyping_file}...")
        present = self.present_calls()
        titles = None
        with gzip.open(self.genotyping_file, 'rt') as ifh:
            for line in ifh:
                line = line.strip()
                if titles is None:
                    titles = line.split(";")
                    assert len(titles) > 1
                    continue
                line = line.split(";", maxsplit=len(titles))
                title = line[0]
                if title in self.insertions.keys():
                    self.insertions[title].nins = sum([1 for gt in line if gt in present])
                    self.insertions[title].nart = sum([1 for gt in line if gt == "artefact"])
                    self.insertions[title].nwt = sum([1 for gt in line if gt == "wild-type"])
        self.insertions = {key: value for key, value in self.insertions.items() if value.nins>0}
        print(f"imported genotypes, found {len(self.insertions)} insertions with one or more tips containing insertions.")

    @staticmethod
    def present_calls():
        """Genotype labels that make a colony a carrier. Legacy: heterozygous / homozygous only.
        The genotyper's `insertion` call (presence certain, zygosity indeterminate -- what the
        low-coverage one-sided loci and far L1DEL/L1DUP pairs typically get with the .tprt
        genotyping keys) counts too in the TPRT pipeline mode (CONFIG['annotate']['rte_library']
        set) or when CONFIG['annotate']['count_insertion_call'] says so; otherwise those loci
        had nins == 0 and were dropped here. Default (no rte keys) stays byte-identical."""
        a = CONFIG['annotate']
        use = a.get('count_insertion_call')
        if use is None:
            use = bool(a.get('rte_library'))
        return ('heterozygous', 'homozygous', 'insertion') if use else ('heterozygous', 'homozygous')

    def generate_fasta_file(self):
        """
        This function generates the fasta file necessary for dfam and samtools
        """
        with gzip.open(self.fasta_file, 'wt') as fh:
            for insertion in self.insertions.values():
                fh.write(insertion.get_fasta())

    def generate_dfam_file(self):
        """
        Run the nucleotide HMM scan of the inserted-sequence clips and write a Dfam-format
        hit table (the format read_dfam parses).

        Two back-ends, selected by config:
          * dfamscan.pl (CONFIG['annotate']['dfamscan'] set and present) — the production
            wrapper that applies the models' GA thresholds and per-chromosome filtering.
          * nhmmscan --dfamtblout (default fallback) — needs only HMMER on PATH. nhmmscan
            emits the identical Dfam table, so no external Perl script is required. Used for
            the test harness and any host with a hmmpress'd HMM library but no dfamscan.pl.
        """
        print(f"running DFAM on {self.sample}")
        assert os.path.exists(self.fasta_file)
        hmm = CONFIG['annotate']['hmm']
        assert os.path.exists(hmm), f"HMM library {hmm} not found (run hmmpress on it first)"
        hmmer = CONFIG['annotate'].get('hmmer')
        if hmmer:
            # prepend the configured HMMER bin dir so nhmmscan / dfamscan.pl resolve.
            os.environ["PATH"] = hmmer + ":" + os.environ["PATH"]
        assert 0 == os.system("nhmmscan -h > /dev/null 2>&1"), "nhmmscan not found on PATH"
        dfamscan = CONFIG['annotate'].get('dfamscan')
        # cores this job owns (LSF), not the whole node's
        cpu = int(os.environ.get("LSB_DJOB_NUMPROC") or os.cpu_count() or 1)
        # nhmmscan cannot read gzip ("Sequence file ... is empty or misformatted"), whether called
        # directly or by dfamscan.pl -> decompress the clip fasta to a temp file for both
        fa = self.fasta_file
        tmp_fa = None
        if fa.endswith(".gz"):
            tmp_fa = self.dfam_file + ".query.fa"
            with gzip.open(fa, "rt") as i, open(tmp_fa, "w") as o:
                o.write(i.read())
            fa = tmp_fa
        try:
            if dfamscan and os.path.exists(dfamscan):
                assert os.access(dfamscan, os.X_OK)
                # a PERL5LIB/PERL_* inherited from the submitting shell's modules points this perl
                # at another perl's XS modules ("ListUtil.c: loadable library and perl binaries
                # are mismatched", PD37449 farm run) -> run dfamscan.pl without perl variables
                env = {k: v for k, v in os.environ.items() if not k.startswith("PERL")}
                rc = subprocess.call([dfamscan, "--fastafile", fa, "--hmmfile", hmm,
                                      "--cpu", str(cpu), "--dfam_outfile", self.dfam_file], env=env)
                assert rc == 0, f"dfamscan.pl failed (exit {rc})"
            else:
                rc = os.system(
                    f"nhmmscan --cpu {cpu} --dfamtblout {self.dfam_file} {hmm} {fa} "
                    f"> /dev/null 2>&1")
                assert rc == 0, "nhmmscan failed"
        finally:
            if tmp_fa and os.path.exists(tmp_fa):
                os.remove(tmp_fa)
        assert os.path.exists(self.dfam_file)

    def generate_sam_file(self):
        """
        This function generates the SAM alignment file with bowtie2
        """
        print(f"running bowtie2 on {self.sample}")
        assert os.path.exists(self.fasta_file)
        assert os.path.exists(CONFIG['combine_insertions']['bowtie2_executable'])
        assert os.access(CONFIG['combine_insertions']['bowtie2_executable'], os.X_OK)
        assert 0 == os.system(f"{CONFIG['combine_insertions']['bowtie2_executable']} -x {CONFIG['combine_insertions']['bowtie2_index2']} --end-to-end -f {self.fasta_file} > {self.sam_file}")
        assert os.path.exists(self.sam_file)
        assert os.path.getsize(self.sam_file)>0

    def generate_sam_local_file(self):
        """(E) bowtie2 --local pass: unlike --end-to-end (which drops a clip unless the WHOLE clip
        aligns to one contiguous block), --local soft-clips the unmatched ends and reports the
        best-matching *sub*-sequence's placement. That places the mappable portion of a chimeric,
        junction-spanning clip. Guarded: silently skipped (leaving the local channel empty, so
        every consumer is a no-op) unless the bowtie2 executable AND its index are actually
        present -- deployment-local cluster paths are absent off-cluster, so this never fires on
        the local re-score, only on a real run."""
        exe = CONFIG['combine_insertions']['bowtie2_executable']
        idx = CONFIG['combine_insertions']['bowtie2_index2']
        if not (os.path.exists(exe) and os.access(exe, os.X_OK) and os.path.exists(idx + '.1.bt2')):
            print(f"bowtie2 --local pass skipped (executable or index {idx} absent); "
                  f"local-remap channel inert for {self.sample}")
            return
        print(f"running bowtie2 --local on {self.sample}")
        rc = os.system(f"{exe} -x {idx} --local -f {self.fasta_file} > {self.local_sam_file}")
        if rc != 0:
            print(f"bowtie2 --local returned {rc}; local-remap channel left empty")

    def demote_flank_leak_dfams(self):
        """(C) Remove Dfam hits whose clip span is genomic flank leaked into the soft-clip (see
        Insertion._dfam_is_flank_leak). Such a hit describes the insertion SITE's old-repeat
        context, not the inserted element, and would otherwise seed a false Alu/L1/SVA hallmark."""
        n = 0
        for ins in self.insertions.values():
            for side, attr in (('left', 'left_dfams'), ('right', 'right_dfams')):
                kept = [d for d in getattr(ins, attr) if not ins._dfam_is_flank_leak(d, side)]
                n += len(getattr(ins, attr)) - len(kept)
                setattr(ins, attr, kept)
        print(f"demoted {n} flank-leak Dfam hit(s) (genomic flank in the soft-clip)")

    def read_dfam(self):
        """
        This function reads the output of DFAM and decorates the insertion object with it.
        """
        n_annot_left = 0
        n_annot_right = 0
        with open(self.dfam_file, 'r') as dfam:
            for line in dfam:
                line = line.strip()
                if not line: continue
                if line[0] == '#':
                    continue
                line = line.split(None, 14)
                insertion = line[2][:-2]
                if not insertion in self.insertions:
                    print(f"insertion {insertion} is not in the insertion list")
                    continue
                if line[2][-1] == 'L':
                    self.insertions[insertion].left_dfams.append(Dfam_Annotation(line))
                    n_annot_right += 1
                elif line[2][-1] == 'R':
                    self.insertions[insertion].right_dfams.append(Dfam_Annotation(line))
                    n_annot_left += 1
                else:
                    raise ValueError(f"Unknown insertion side {line[2]}, expected R or L.")
        print(f"read dfam for sample {self.sample}, found {n_annot_left} dfam annotations left and {n_annot_right} dfam annotations right.")

    def read_rmsk(self, rmsk_library_path):
        """
        reads the repeatmasker library from UCSC
        """
        print(f"reading rmsk file {rmsk_library_path}")
        rmsk_library = {}
        n_rmsk = 0
        with gzip.open(rmsk_library_path, 'rt') as rmsk:
            units = None
            titles = None
            for line in rmsk:
                line = line.strip()
                if units is None:
                    units = line.split(None)
                    continue
                if titles is None:
                    titles = line.split(None)
                    continue
                if line == '': continue
                line = line.split(None, maxsplit=15)
                r = RepeatMasker_Annotation(line)
                if r.repClass in ('SINE','LINE','LTR') and r.strand in ('+','-'):
                    if line[4] not in rmsk_library.keys():
                        rmsk_library[line[4]] = []
                    n_rmsk += 1
                    rmsk_library[line[4]].append(RepeatMasker_Annotation(line))
        print(f"imported {n_rmsk} rmsk entries for LINE, SINE and LTRs from {rmsk_library_path}")
        return rmsk_library

    def read_exons(self, exon_path):
        """Read a gene-exon BED (contig<TAB>start<TAB>end<TAB>gene_id, 0-based half-open;
        plain or gzipped) into {contig: [(start, end, gene_id), ...]} sorted by start, for
        processed-pseudogene detection. Same format as discovery Feature B's exon track."""
        print(f"reading exon annotation {exon_path}")
        lib = {}
        n = 0
        opener = gzip.open if str(exon_path).endswith(".gz") else open
        with opener(exon_path, 'rt') as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line[0] == '#':
                    continue
                f = line.split("\t")
                if len(f) < 4:
                    continue
                contig, start, end, gene = f[0], int(f[1]), int(f[2]), f[3]
                lib.setdefault(contig, []).append((start, end, gene))
                n += 1
        for contig in lib:
            lib[contig].sort()
        print(f"imported {n} exons for {len(lib)} contigs from {exon_path}")
        return lib

    def get_exon(self, exon_library: dict, pos: tuple, pad: int = 100) -> list:
        """Exons overlapping the mapped clip position (+/- pad). Returns (gene, start, end)."""
        seqname, p = pos
        exons = exon_library.get(seqname)
        if not exons:
            return []
        lo, hi = 0, len(exons)
        while lo < hi - 1:                       # binary search to the neighbourhood by start
            mid = (lo + hi) // 2
            if exons[mid][0] < p - 500:
                lo = mid
            else:
                hi = mid
        out = []
        for i in range(lo, len(exons)):
            s, e, gene = exons[i]
            if s > p + pad:
                break
            if e >= p - pad:                     # exon interval overlaps [p-pad, p+pad]
                out.append((gene, s, e))
        return out

    def get_rmsk(self, rmsk_library: dict[str, RepeatMasker_Annotation], pos: tuple[str, int, str]) -> list[RepeatMasker_Annotation]:
        seqname, start = pos
        if not seqname in rmsk_library.keys():
            print(f"seqname {seqname} not in rmsk library")
            return []
        search_pos = start - 500
        lpos = 0
        upos = len(rmsk_library[seqname])
        while lpos < upos-1:
            spos = lpos + floor((upos-lpos)/2)
            if rmsk_library[seqname][spos].start < search_pos:
                lpos = spos
            else:
                upos = spos
        output = []
        for lpos in range(lpos, len(rmsk_library[seqname])):
            if rmsk_library[seqname][lpos].start > start+8500: continue
            delta_start = rmsk_library[seqname][lpos].start - start
            if delta_start > 100: continue
            delta_end = rmsk_library[seqname][lpos].end - start
            if delta_end < -100: continue
            if delta_end > -100 and delta_start < 100:
                output.append(rmsk_library[seqname][lpos])
        return output




    def read_sam(self):
        """
        this function reads the output of bowtie2 and decorates the insertion object with it.
        positions are lifted over using the chainfile, if one is provided. This is usefull if re-mapping is done to another genome version, e.g. to hs1 if originally mapped to hg38.
        """
        #if CONFIG['combine_insertions']['bowtie2_index2_lo'] is not None:
        #    print(f"reading CONFIG['combine_insertions']['bowtie2_index2_lo'] {CONFIG['combine_insertions']['bowtie2_index2_lo']}, this might take a while...")
        #    lo = pyliftover.LiftOver(CONFIG['combine_insertions']['bowtie2_index2_lo'])
        #    print(f"done reading chainfile {CONFIG['combine_insertions']['bowtie2_index2_lo']}")
        #else:
        #    lo = None
        lo = None
        rmsk_library = self.read_rmsk(CONFIG['annotate']['rmsk'])
        exon_path = CONFIG['annotate'].get('exon_annotation')
        exon_library = self.read_exons(exon_path) if exon_path else {}

        rightn = 0
        leftn = 0
        with pysam.AlignmentFile(self.sam_file) as sam:
            for read in sam:
                if read.is_qcfail: continue
                #if read.is_secondary: continue
                #if read.is_supplementary: continue
                if read.is_unmapped: continue
                local_rmsks = None
                local_exons = []
                insertion = read.query_name[:-2]
                if not insertion in self.insertions.keys(): continue
                if (read.query_name[-1] == "R") ^ read.is_forward:
                    local_rmsks = self.get_rmsk(rmsk_library, (read.reference_name, read.reference_start))
                    if exon_library:
                        local_exons = self.get_exon(exon_library, (read.reference_name, read.reference_start))
                    if lo is not None:
                        co = lo.convert_coordinate(read.reference_name, read.reference_start,
                                                   '+' if read.is_forward else '-')
                    else:
                        co = [(read.reference_name, read.reference_start, '+' if read.is_forward else '-')]
                else:
                    local_rmsks = self.get_rmsk(rmsk_library, (read.reference_name, read.reference_end))
                    if exon_library:
                        local_exons = self.get_exon(exon_library, (read.reference_name, read.reference_end))
                    if lo is not None:
                        co = lo.convert_coordinate(read.reference_name, read.reference_end,
                                                   '+' if read.is_forward else '-')
                    else:
                        co = [(read.reference_name, read.reference_end, '+' if read.is_forward else '-')]
                if co:
                    if read.query_name[-1] == "R":
                        rightn += 1
                        self.insertions[insertion].right_maps.append((f"{co[0][0]}:{co[0][1]}{co[0][2]}", read.mapping_quality, local_rmsks, co[0][2]))
                        # Exon evidence is recorded regardless of MAPQ: a processed pseudogene
                        # of a MULTI-COPY source gene has clips that legitimately multi-map
                        # (many genomic copies), so a MAPQ floor here would silently blind us to
                        # exactly those genes. Specificity comes from _pseudogene() instead — it
                        # requires two distinct exons of ONE gene (or one exon + a poly-A tail),
                        # which a chance multimapper does not satisfy — and RTE clips are resolved
                        # by the Dfam/rmsk branches before _pseudogene() is ever consulted.
                        self.insertions[insertion].right_exons.extend(local_exons)
                    elif read.query_name[-1] == "L":
                        leftn += 1
                        self.insertions[insertion].left_maps.append((f"{co[0][0]}:{co[0][1]}{co[0][2]}", read.mapping_quality, local_rmsks,
                                                                        co[0][2]))
                        self.insertions[insertion].left_exons.extend(local_exons)
                    else:
                        raise ValueError(f"Unknown insertion side {read.query_name}, expected R or L.")
        print(f"imported {rightn} right mappings and {leftn} left mappings.")

    def read_sam_local(self):
        """(E) Read the bowtie2 --local SAM and record each clip's partial placement:
        (contig, ref_start, strand, query_coverage_bp, mapq) on the insertion's
        left_local_maps / right_local_maps. query_coverage_bp is how much of the clip aligned
        (the rest was soft-clipped) -- a chimeric clip aligns only its mappable half. Coordinates
        are NOT lifted over (the local pass targets the same index as the main remap)."""
        rightn = leftn = 0
        with pysam.AlignmentFile(self.local_sam_file) as sam:
            for read in sam:
                if read.is_qcfail or read.is_unmapped:
                    continue
                name = read.query_name
                insertion = name[:-2]
                if insertion not in self.insertions:
                    continue
                qcov = read.query_alignment_length or 0
                rec = (read.reference_name, read.reference_start,
                       '+' if read.is_forward else '-', qcov, read.mapping_quality)
                if name[-1] == 'R':
                    self.insertions[insertion].right_local_maps.append(rec); rightn += 1
                elif name[-1] == 'L':
                    self.insertions[insertion].left_local_maps.append(rec); leftn += 1
        print(f"imported {rightn} right and {leftn} left --local (split) placements.")

    # Coarse element class from a conclusion() string, for the flat table. The order
    # matters: the non-RTE SV flag and the pseudogene calls also read as "unknown"/contain
    # element-like tokens, so the specific cases are tested before the generic ones.
    @staticmethod
    def element_class(conclusion: str) -> str:
        # The element identity wins over an additive `[SV: ...]` suffix: a rearrangement
        # junction that falls inside a retrotransposon is still classed by that element. Strip
        # the suffix and classify the element head, so even an RTE_other element (which carries
        # no big-three token) keeps its class. A pure SV has no element head — the whole string
        # is the SV descriptor — so it still falls through to the non_RTE_SV umbrella below.
        head = conclusion.split(' [SV:', 1)[0].split(' [site:', 1)[0]
        u = head.upper()
        if 'PSEUDOGENE' in u:
            return 'processed_pseudogene'
        if head == 'artefact':
            return 'artefact'
        if 'SVA' in u:
            return 'SVA'
        if 'ALU' in u:
            return 'ALU'
        if 'L1' in u or 'LINE' in u:
            return 'LINE1'
        if 'MICROSATELLITE' in u:
            return 'microsatellite'
        if 'TEMPLATED' in u:
            return 'templated_insertion'   # (E) distal-sourced; kept separate from the clean SV set
        if any(t in u for t in ('TRANSLOCATION', 'INVERSION', 'INTRACHROMOSOMAL SV',
                                'DELETION', 'DUPLICATION', 'NON-RTE SV')):
            return 'non_RTE_SV'
        if u.startswith('UNKNOWN'):
            return 'unknown'
        return 'RTE_other'   # a mapped/dfam RTE that is none of the big three (HERV/LTR, MIR, ...)

    _TEMPLATED_SOURCE = re.compile(r"templated/complex insertion \(source ([^:\s()]+):(\d+)\)")

    def locus_class(self, key, concl):
        """(class, conclusion) for the table: element_class() of the legacy conclusion, overridden
        by a confident tools/rte verdict when the legacy class is only a fallback (unknown /
        artefact / templated_insertion) -- precedence documented in tools/rte/locus_class.py.
        Without tools/rte (no rte_library) this is exactly element_class() + the conclusion."""
        cls = self.element_class(concl)
        rte = getattr(self, 'rte_records', None) or {}
        if not rte:
            return cls, concl
        try:
            from tools.rte.locus_class import locus_class
        except ImportError:                       # run as `python tools/annotate_v2.py`
            from rte.locus_class import locus_class
        m = self._TEMPLATED_SOURCE.search(concl)
        tsrc = (m.group(1), int(m.group(2))) if m else None
        new, note = locus_class(cls, rte.get(key), getattr(self, 'rte_lib', None), tsrc,
                                CONFIG['annotate'].get('rte_remap_assembly', 'hs1'))
        if new != cls and note:
            concl = f"{concl}; RTE: {note}"
        return new, concl

    def write_table(self, path: str):
        """Flat, machine-readable annotation table: one row per called locus. Complements
        print() (the verbose per-junction report). Gzipped when `path` ends in .gz."""
        opener = gzip.open if path.endswith('.gz') else open
        cols = ['locus', 'class', 'conclusion', 'n_ins', 'n_wt', 'n_art',
                'left_polyA', 'right_polyA', 'left_dfam', 'right_dfam', 'left_map', 'right_map',
                'site_region', 'site_gene', 'site_strand']
        rte = getattr(self, 'rte_records', None) or {}
        if rte:
            # tools/rte columns (plans/tprt_hallmarks/SPEC.md), from the structured RteRecord
            rte_cols = next(iter(rte.values())).COLUMNS
            cols = cols + rte_cols
        n = 0
        with opener(path, 'wt') as fh:
            fh.write('\t'.join(cols) + '\n')
            for key, ins in self.insertions.items():
                cls, concl = self.locus_class(key, ins.conclusion())
                ld = ','.join(sorted({m.model for m in ins.left_dfams})) or '.'
                rd = ','.join(sorted({m.model for m in ins.right_dfams})) or '.'
                lm = (','.join(sorted({p for p, _, _, _ in ins.left_maps})) or '.')[:80]
                rm = (','.join(sorted({p for p, _, _, _ in ins.right_maps})) or '.')[:80]
                site = ins.site()
                sregion, sgene, sstrand = (site[0], site[1], site[2]) if site else ('.', '.', '.')
                row = [key, cls, concl.replace('\t', ' ').replace('\n', ' '),
                       str(ins.nins), str(ins.nwt), str(ins.nart),
                       'Y' if ins.has_left_polyA() else 'N',
                       'Y' if ins.has_right_polyA() else 'N',
                       ld, rd, lm, rm, sregion, sgene, sstrand]
                if rte:
                    rec = rte.get(key)
                    row += rec.row() if rec is not None else ['.'] * len(rte_cols)
                fh.write('\t'.join(row) + '\n')
                n += 1
        print(f"wrote annotation table ({n} loci) to {path}")

    def print(self):
        for key, insertion in self.insertions.items():
            print(f"> insertion {key} found in {insertion.nins} tip{'s' if insertion.nins!=1 else ''} ({insertion.nwt}=wt, {insertion.nart}=art)")
            print(f" {insertion.conclusion()}")
            site = insertion.site()
            if site is not None:
                print(f"  SITE: {site[3]}")
            rec = (getattr(self, 'rte_records', None) or {}).get(key)
            if rec is not None:
                print(f"  RTE: {rec.element} {rec.structure} [{','.join(rec.tags) or '-'}] "
                      f"{rec.tprt_call} score={rec.tprt_score} ({rec.tprt_points})")
            print(f"  RIGHT INSERTION: {insertion.right_seq}")
            for dfam in insertion.right_dfams:
                print(f"    {str(dfam)}")
            else:
                print(f"    - no dfam entries found")
            for pos, mq, rmsks, strand in insertion.right_maps:
                print(f"    {pos} {[str(r) for r in rmsks]}")
            if insertion.has_right_polyA():
                print(f"    is polyA")
            print(f"  LEFT INSERTION: {insertion.left_seq}")
            for dfam in insertion.left_dfams:
                print(f"    {str(dfam)}")
            else:
                print(f"    - no dfam entries found")
            for pos, mq, rmsks, strand in insertion.left_maps:
                print(f"    {pos} {[str(r) for r in rmsks]}")
            if insertion.has_left_polyA():
                print(f"    is polyA")


if __name__== '__main__':
    if not len(sys.argv) > 1:
        raise NotImplementedError(f"this script takes one or two arguments.\nUsage: annotate_v2.py [sample_name] ([output_path])")
    sample = sys.argv[1]
    if len(sys.argv)>2:
        output_path = sys.argv[2]
    else:
        output_path = f"{sample}.out"
    f = VariantAnnotationContainer(sample, output_path)
    # output_path is the flat table (pipeline.sh expects <PATIENT_ID>.annotated.csv.gz);
    # the verbose per-junction report still goes to stdout (captured in the job log).
    f.write_table(output_path)
    f.print()




