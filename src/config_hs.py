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


# configuration file for PEAR-TREE

CONFIG = {

    # settings for step 1: discovery phase:
    'discovery': {
        'min_mapq': 60, #minimal mapq of an alignment to be considered in discovery phase
        'min_clip_len': 12, #minimal length of clipped bases to be considered
        'min_evidence_reads_per_breakpoint': 2, #minimal number of independent evidence reads for each breakpoint
        'max_bp_window': 40, #max number of allowed bases deleted or duplicated between two breakpoints
        'max_homopolymer_len': 6, #maximal length of homopolymer in mapped part immediately adjacent to the breakpoint
        'min_adapterlen_for_clip': 4, #minimal adapter length to be clipped from an already clipped part.
        'min_good_bases': 10, #minimal number of bases of a clipped read to have a good quality*
        'min_consensus_score_for_good_base': 2, #minimal score for a good base (delta best hit vs. second best hit)
        'min_breakpoints_aggregated_during_first_step': 2, #minimal number of breakpoints aggregated during first step in order to even start filtering
        'max_read_count': 120, #maximum numbers of reads in the span of a breakpoint allowed (exclude high-coverage artefact-rich regions)
        'exclude_same_contig_supplementary': 1000,  # minimum distance between a supplementary read to not be excluded (not interested in micro indels)
        'reject_fully_mapping_reads': True, #drop clipped reads whose XA/SA shows the whole read maps contiguously elsewhere (not a real junction)
        # SPEC-3 pileup gate — recommended ON for human WGS. Drops breakpoints whose local
        # depth exceeds coverage_mask_multiplier x the genome-wide median, removing the
        # pericentromere/telomere classical-satellite mismap pileups that stack to tens of x
        # median at MAPQ 60 (so neither the MAPQ floor nor the combine remap catches them).
        # Validated on test/fullstack/scale10k: 10k implants -> 0 genuine FP, recall unchanged.
        # NB consumed by the *rust* discovery via `--config discovery_hs.config`; the legacy
        # Python discovery predates the SPEC-3 gate and ignores these keys.
        'coverage_mask': True,
        'coverage_mask_multiplier': 5.0,
        'adaptive_evidence': False, #SPEC-4 also cuts satellite FPs but loses low-VAF TPs -> leave off
        # mate_anchor_rescue — recommended ON for human WGS. Accepts a soft-clipped read
        # below min_mapq when its mate maps uniquely (MQ tag >= min_mapq); recovers
        # insertions into low-mapability-but-mate-unique flanks. On test/fullstack/scale10k
        # it added +62 true insertions (recall 95.4 -> 96.0%) with 0 added FP after combine —
        # the sweep's clean sensitivity win. Needs the MQ tag from `samtools fixmate`.
        # NB: for the rust discovery the operative settings (incl. min_mapq lowered to 40)
        # live in test/fullstack/scale10k/discovery_hs.config; pass that with `--config`.
        'mate_anchor_rescue': True,
    },

    #define adapter sequences
    'adapters': ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',  # fwd NebNext Adapter
                'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',  # rev NebNext Adapter
                'AGATCGGAAAGCACACGTCTGAACTCCAGTCA',  # common sequencing error of fwd adapter
                'AGATCGGAAAGCGTCGTGTAGGGAAAGAGTGT',  # common sequencing error of rev adapter
                ],
    'genotyping': {
        'max_bases': 12, #max bases (per side) used for genotyping (applied in combine_insertions)
        'min_mapq': 40, #minimal mapq of a spanning read to be used as genotyping evidence
        'min_score_for_call': 6, #minimal aggregate quality-margin for a confident (non-uncertain) call
        'min_supporting_reads': 2, #minimal number of allele-supporting reads for a confident het/hom/wt call
        'min_reads_for_zygosity': 6, #min informative reads before an alt-dominant locus is called homozygous rather than 'insertion' (zygosity unclear)
        'recover_low_coverage_presence': True, #promote a single strong alt read at a known locus to 'insertion' when the reference allele is not confidently present (n_ref<min_supporting_reads); disable where index hopping is uncontrolled
        'reads_for_high_coverage': 60, #read count above which a locus is flagged high-coverage (counted as NA)
        'art_min_score': 60, #per-side quality above which a read side counts as artefact (matches neither ref nor alt)
        'vaf_wildtype_max': 0.10, #VAF at or below which a locus is called wild-type
        'vaf_het_min': 0.30, #minimal VAF for a heterozygous call
        'vaf_hom_min': 0.85, #minimal VAF for a homozygous call
        'artefact_read_fraction': 0.5, #fraction of artefact reads (of all evidence) needed to call the locus an artefact
        'min_artefact_reads': 2, #minimal number of artefact reads to call the locus an artefact
        'double_alt_is_artefact': True, #a single read matching the inserted element on BOTH junctions is a chimeric artefact
    },
    'combine_genotypes': {
        'min_wild-types': 20, #minimal number of wild-type colonies. For enriched/targeted (non-WGS) data set this to 0 and raise max_na.
        'min_insertions': 1, #minimal number of colonies with insertion
        'max_artefact': 24, #maximal number of colonies with artefacts
        'max_na': 24, #maximal number of colonies with NA genotype (high coverage, no coverage or error)
        'min_best_score': 800, #minimal best het/hom support score across samples for a locus to pass
        # CLONALITY gate — recommended ON for any multi-colony tree study. Every other
        # gate here counts CALLS, and a call is a thresholded allele fraction, so none of
        # them can separate a locus that is genuinely present in some colonies from one
        # whose alt reads are a constant per-locus error rate that the VAF bands slice at
        # random. This tests the read counts directly: a real clonal event must be
        # over-dispersed relative to a single binomial (~4-5 at 10 informative reads),
        # whereas a constant-rate artefact sits at ~1.
        # NB `min_wild-types` alone is a BAND-PASS on the allele fraction, not a germline
        # filter: on PD44579 (174 colonies) loci with p_hat >= 0.4 passed at 0.0-0.1%
        # while p_hat 0.05-0.30 passed at 36-45% — i.e. it removed clean germline and
        # kept exactly the loci too ambiguous to call. It reported 3,792 loci of which
        # 3,641 were statistically indistinguishable from random colony sets; with this
        # gate the same data yields 17.
        # COST (simulated at PD44579 depths, 174 colonies): at 10 informative reads it
        # keeps 86% of true private events and ~100% of k>=2; at 7 reads only 62% of
        # private. It also drops sub-clonal events (carrier VAF ~0.2 -> <30% kept), so
        # raise/lower with depth and lower it for impure colonies. Set to 0/None to
        # disable. Needs n_alt/n_ref in the genotype files.
        'min_dispersion': 3.0,
    },
    'combine_insertions': {
        'genome_2bit': '/Users/jeremy/Documents/genomes/hg38.2bit', #path to genome, has to be 2bit file
        'exclude_files_with_many_insertions': 1_000_000, #exclude files with more than this number of insertions. No single-leaf insertions of these files can be identified.
        'samtools_executable': '/opt/homebrew/bin/samtools', #path to bowtie2 executable
        'bowtie2_executable': '/opt/homebrew/bin/bowtie2', #path to bowtie2 executable
        'bowtie2_index': '/Users/jeremy/Documents/genomes/bowtie2_indices/hs1' #path to bowtie2 index
    },
    'version': '1.1'
}


