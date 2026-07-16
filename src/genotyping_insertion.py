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


import gzip
import pysam
from genotyping_evidence_read import EvidenceRead
from typing import Iterator, List
from config import CONFIG

# genotype call vocabulary. These strings are the on-disk contract consumed by
# combine_genotypes.py, so keep the two in sync.
GT_ARTEFACT = 'artefact'
GT_WILDTYPE = 'wild-type'
GT_HETEROZYGOUS = 'heterozygous'
GT_HOMOZYGOUS = 'homozygous'
# low-confidence calls (evidence points one way but below the confidence gate,
# see CONFIG['genotyping']['min_score_for_call'] / 'min_supporting_reads').
GT_INSERTION_UNCERTAIN = 'insertion?'
GT_WILDTYPE_UNCERTAIN = 'wild-type?'
# not-assessable calls, counted as NA by combine_genotypes.
GT_NO_COVERAGE = 'no-coverage'
GT_HIGH_COVERAGE = 'high-coverage'
GT_ERROR = 'error'


class Insertion:
    def __init__(self, name: str):
        self.name = name
        # A locus name is "<contig>:<left>-<right>". Split from the right so that
        # contigs whose own names contain ':' (hg38 ALT/HLA e.g. HLA-A*01:01:01)
        # or '-' are parsed correctly instead of crashing the unpack.
        self.chr, pos = name.rsplit(":", 1)
        left, right = pos.rsplit("-", 1)
        self.left_pos, self.right_pos = int(left), int(right)
        self.left_clipped = None
        self.right_clipped = None
        self.left_ref = None
        self.right_ref = None
        self.evidence_reads: List[EvidenceRead] = []
    def __str__(self):
        return self.name

    @staticmethod
    def import_file(path: str) -> Iterator['Insertion']:
        insertion = None
        with gzip.open(path, 'rt') as f:
            for line in f:
                line = line.strip()
                if not len(line): continue
                if line[0] == ">":
                    if insertion: yield insertion
                    insertion = Insertion(line[1:])
                    status = None
                    continue
                if line[0] == "@":
                    status = line[1:]
                    continue
                if status == 'RIGHT_INSERTION':
                    insertion.right_clipped = line
                elif status == 'LEFT_INSERTION':
                    insertion.left_clipped = line
                elif status == 'RIGHT_REFERENCE':
                    insertion.right_ref = line
                elif status == 'LEFT_REFERENCE':
                    insertion.left_ref = line
        if insertion: yield insertion
    def genotype(self, bam: pysam.AlignmentFile):
        start = max(0, min(self.left_pos, self.right_pos) - 1)
        end = max(self.left_pos, self.right_pos) + 1
        min_mapq = CONFIG['genotyping'].get('min_mapq', 40)
        seen_qnames = set()
        for read in bam.fetch(self.chr, start, end):
            if read.mapping_quality < min_mapq: continue
            if read.is_secondary: continue
            if read.is_supplementary: continue
            if read.is_unmapped: continue
            if read.is_qcfail: continue
            if read.is_duplicate: continue
            if read.reference_name != self.chr: continue
            # A read without a stored sequence/qualities (SEQ '*') cannot be scored.
            if read.query_sequence is None or read.query_qualities is None: continue
            # Deduplicate by fragment: the two overlapping mates of one fragment,
            # or PCR/optical duplicates sharing a qname, must count as one piece of
            # evidence, not several (the old "remove qname duplicates" comment
            # described behaviour that did not exist).
            if read.query_name in seen_qnames: continue
            seen_qnames.add(read.query_name)
            evi_read = EvidenceRead(read)
            # left_pos/right_pos are genomic coordinates; test "is not None" rather
            # than truthiness so a legitimate coordinate 0 (contig start) is scored.
            if self.left_pos is not None:
                evi_read.qleft(self.left_pos, self.left_ref, self.left_clipped)
            if self.right_pos is not None:
                evi_read.qright(self.right_pos, self.right_ref, self.right_clipped)
            self.evidence_reads.append(evi_read)
        return

    @staticmethod
    def _side_call(ref: float, alt: float, art: float, art_min: float):
        """Classify one breakpoint side of one read as ref/alt/art.

        Returns (kind, margin) where kind is 'ref' | 'alt' | 'art' | None and
        margin is the (non-negative) confidence of that side. None means the side
        did not cover the breakpoint or was an exact tie (uninformative).
        """
        if ref == 0 and alt == 0 and art == 0:
            return None, 0.0
        if art > ref and art > alt and art >= art_min:
            return 'art', float(art)
        if alt > ref:
            return 'alt', float(alt - ref)
        if ref > alt:
            return 'ref', float(ref - alt)
        return None, 0.0

    def summarise_evidence(self):
        """Genotype the locus from per-read evidence using an allele-fraction model.

        Each spanning read casts ONE vote — ref (wild-type), alt (carries the
        inserted junction) or art (artefact) — and the call is driven by the
        variant allele fraction VAF = n_alt / (n_alt + n_ref), not by raw
        quality-sum deltas. This makes het/hom/wild-type a function of allele
        balance (het ~0.5, hom ~1.0 in clonal material) and lets a few
        contaminating alt reads read as wild-type rather than a spurious het.
        Thresholds are configurable under CONFIG['genotyping'].
        """
        cfg = CONFIG['genotyping']
        art_min = cfg.get('art_min_score', 60)
        double_alt_is_artefact = cfg.get('double_alt_is_artefact', True)
        vaf_wildtype_max = cfg.get('vaf_wildtype_max', 0.10)
        vaf_het_min = cfg.get('vaf_het_min', 0.30)
        vaf_hom_min = cfg.get('vaf_hom_min', 0.85)
        artefact_read_fraction = cfg.get('artefact_read_fraction', 0.5)
        min_artefact_reads = cfg.get('min_artefact_reads', 2)
        min_supporting_reads = cfg.get('min_supporting_reads', 2)
        min_score_for_call = cfg.get('min_score_for_call', 6)

        if not len(self.evidence_reads):
            return GT_NO_COVERAGE, 0, 0

        n_ref = n_alt = n_art = 0
        ref_score = alt_score = art_score = 0.0
        for er in self.evidence_reads:
            left = self._side_call(*er.left_genotype, art_min)
            right = self._side_call(*er.right_genotype, art_min)
            covered = [s for s in (left, right) if s[0] is not None]
            if not covered:
                continue
            kinds = {kind for kind, _ in covered}
            if 'art' in kinds:
                n_art += 1
                art_score += sum(m for kind, m in covered if kind == 'art')
                continue
            if double_alt_is_artefact and kinds == {'alt'} and len(covered) == 2:
                # Both junctions of a single short read match the inserted element:
                # geometrically impossible for a real (long) insertion whose two
                # ends are hundreds of bp apart, so this is a chimeric artefact.
                # NOTE (element-class caveat): combine_insertions already drops loci
                # whose alt consensus resembles the local reference, so a genuine
                # wild-type read cannot land here by flank homology alone; set
                # double_alt_is_artefact=False for element families (e.g. mouse
                # ERV/LTR) where that guarantee is weaker.
                n_art += 1
                art_score += sum(m for _, m in covered)
                continue
            if 'alt' in kinds:
                n_alt += 1
                alt_score += sum(m for kind, m in covered if kind == 'alt')
            else:
                n_ref += 1
                ref_score += sum(m for kind, m in covered if kind == 'ref')

        informative = n_ref + n_alt
        total = informative + n_art

        if total == 0:
            return GT_NO_COVERAGE, 0, 0
        # artefact-dominated locus: enough artefact reads AND a dominant fraction.
        if n_art >= min_artefact_reads and n_art >= artefact_read_fraction * total:
            return GT_ARTEFACT, int(art_score), int(max(ref_score, alt_score))
        if informative == 0:
            # only a handful of artefact reads, below the artefact threshold:
            # not assessable rather than a (false) wild-type.
            return GT_NO_COVERAGE, 0, 0

        vaf = n_alt / informative
        confident_ins = n_alt >= min_supporting_reads and alt_score >= min_score_for_call
        if vaf >= vaf_hom_min:
            gt = GT_HOMOZYGOUS if confident_ins else GT_INSERTION_UNCERTAIN
            return gt, int(alt_score), int(ref_score)
        if vaf >= vaf_het_min:
            gt = GT_HETEROZYGOUS if confident_ins else GT_INSERTION_UNCERTAIN
            return gt, int(alt_score), int(ref_score)
        if vaf > vaf_wildtype_max:
            # some alt evidence but too little for a clean het — likely a subclone,
            # index hopping or residual artefact; flag as uncertain wild-type.
            return GT_WILDTYPE_UNCERTAIN, int(ref_score), int(alt_score)
        confident_wt = n_ref >= min_supporting_reads and ref_score >= min_score_for_call
        gt = GT_WILDTYPE if confident_wt else GT_WILDTYPE_UNCERTAIN
        return gt, int(ref_score), int(alt_score)
