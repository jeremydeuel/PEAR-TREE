# configuration file for GEMINATE

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
        'max_read_count': 60, #maximum numbers of reads in the span of a breakpoint allowed (exclude high-coverage artefact-rich regions)
    },

    #define adapter sequences
    'adapters': ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',  # fwd NebNext Adapter
                'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',  # rev NebNext Adapter
                'AGATCGGAAAGCACACGTCTGAACTCCAGTCA',  # common sequencing error of fwd adapter
                'AGATCGGAAAGCGTCGTGTAGGGAAAGAGTGT',  # common sequencing error of rev adapter
                ],
    'genotyping': {
        'max_bases': 12, #max bases used for genotyping
        'min_score_for_call': 6, #minimal score to call
        'reads_for_high_coverage': 60, #number of reads required to call high coverage
    },
    'combine_genotypes': {
        'min_wild-types': 20, #minimal number of wild-type colonies
        'min_insertions': 1, #minimal number of colonies with insertion
        'max_artefact': 200, #maximal number of colonies with artefacts
        'max_na': 18, #maximal number of colonies with NA genotype (high coverage or no coverage)
    },
    'combine_insertions': {
        'genome_2bit': '/Users/jeremy/Desktop/isa_all/mm39.2bit', #path to genome, has to be 2bit file
        'exclude_files_with_many_insertions': 1_000_000, #exclude files with more than this number of insertions. No single-leaf insertions of these files can be identified.
        'samtools_executable': '/opt/homebrew/bin/samtools', #path to bowtie2 executable
        'bowtie2_executable': '/opt/homebrew/bin/bowtie2', #path to bowtie2 executable
        'bowtie2_index': '/Users/jeremy/Desktop/isa_all/mm39/mm39' #path to bowtie2 index
    },
    'version': '1.0'
}


