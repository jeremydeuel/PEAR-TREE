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


# configuration file for PEAR-TREE, specific for mus musculus

CONFIG = {

    # settings for step 1: discovery phase:
    'discovery': {
        'min_mapq': 40, #minimal mapq of an alignment to be considered in discovery phase
        'min_clip_len': 12, #minimal length of clipped bases to be considered
        'min_evidence_reads_per_breakpoint': 2, #minimal number of independent evidence reads for each breakpoint
        'max_bp_window': 40, #max number of allowed bases deleted or duplicated between two breakpoints
        'max_homopolymer_len': 6, #maximal length of homopolymer in mapped part immediately adjacent to the breakpoint
        'min_adapterlen_for_clip': 4, #minimal adapter length to be clipped from an already clipped part.
        'min_good_bases': 10, #minimal number of bases of a clipped read to have a good quality*
        'min_consensus_score_for_good_base': 2, #minimal score for a good base (delta best hit vs. second best hit)
        'min_breakpoints_aggregated_during_first_step': 2, #minimal number of breakpoints aggregated during first step in order to even start filtering
        'max_read_count': 120, #maximum numbers of reads in the span of a breakpoint allowed (exclude high-coverage artefact-rich regions)
        'exclude_same_contig_supplementary': 1000, #minimum distance between a supplementary read to not be excluded (not interested in micro indels)
        'reject_fully_mapping_reads': True, #drop clipped reads whose XA/SA shows the whole read maps contiguously elsewhere (not a real junction)
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
        'min_wild-types': 20, #minimal number of wild-type colonies
        'min_insertions': 1, #minimal number of colonies with insertion
        'max_artefact': 0, #maximal number of colonies with artefacts
        'max_na': 0, #maximal number of colonies with NA genotype (high coverage, no coverage or error)
        'min_best_score': 800, #minimal best het/hom support score across samples for a locus to pass
    },
    'combine_insertions': {
        'genome_2bit': '/lustre/scratch125/casm/teams/team273/users/jd43/pt_hu_trees/hs1/hs1.2bit', #path to genome, has to be 2bit file
        'exclude_files_with_many_insertions': 1_000_000, #exclude files with more than this number of insertions. No single-leaf insertions of these files can be identified.
        'samtools_executable': '/software/spack_environments/default/00/opt/spack/linux-ubuntu22.04-x86_64_v3/gcc-13.1.0/samtools-1.19-ufgcnbyuj24ufumlyimozp6habconpvy/bin/samtools', #path to bowtie2 executable
        'bowtie2_executable': '/nfs/users/nfs_j/jd43/software/bowtie2-2.5.4-linux-x86_64/bowtie2', #path to bowtie2 executable
        'bowtie2_index': '/lustre/scratch125/casm/teams/team273/users/jd43/pt_hu_trees/hs1', #path to bowtie2 index used to align the entire string of aligned read and clipped read, use a t2t genome if available.
        'bowtie2_index2': '/lustre/scratch125/casm/teams/team273/users/jd43/pt_hu_trees/hs1', #path to bowtie2 index used to align only the clipped part of the read, use the latest annotated genome here.
        'bowtie2_index2_lo': '' # path to chainfile linking bowtie2_index to the index used for alignment in the bam file.

    },
    'version': '1.0'
}


