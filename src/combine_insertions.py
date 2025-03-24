import gzip
import sys, os
from revcomp import revcomp
import pysam
from config import CONFIG
from combine_insertions_insertion import Insertion
from combine_insertions_intersect_insertions import intersect_insertions
from combine_insertions_get_sequence import get_sequence


def combine_insertions(input_files, insertions_genotyping_file, combined_insertions, insertions_fasta, insertion_bam, threads):
    all_insertions = {}
    for f in input_files:
        fb = os.path.basename(f)
        all_insertions[fb] = {i.name: i for i in Insertion.parseFile(f) if len(i.chr) < 6 and i.chr != "MT" and i.chr !="chrM"}
        print(f"\033[37mFile {fb}: imported {len(all_insertions[fb])} insertions\033[0m")
        if len(all_insertions[fb]) > CONFIG['combine_insertions']['exclude_files_with_many_insertions']:
            del all_insertions[fb]
            print(f"\033[31mremoved file {f} since it contains too many insertions.\033[0m")
    print(f"intersecting insertions from {len(all_insertions)} files...")
    insertions = intersect_insertions(all_insertions)
    print(f"writing summarised insertions fasta file {insertions_fasta}")
    with gzip.open(insertions_fasta, 'wt') as f:
        f.writelines(
            [f'>{i.name}:L\n{i.left_consensus}\n\n>{i.name}:R\n{i.right_consensus}\n\n' for i in insertions])
    print(f"running bowtie2 {CONFIG['combine_insertions']['bowtie2_executable']} with index {CONFIG['combine_insertions']['bowtie2_index']}")
    cmd = f"{CONFIG['combine_insertions']['bowtie2_executable']} -f {insertions_fasta} -x {CONFIG['combine_insertions']['bowtie2_index']} --end-to-end --sensitive --threads {threads} --qc-filter | {CONFIG['combine_insertions']['samtools_executable']} view -F 4 -b -o {insertion_bam}"
    os.system(cmd)
    filter_reads = set()
    with pysam.AlignmentFile(insertion_bam) as f:
        for read in f:
            filter_reads.add(read.query_name[:-2])
    print(
        f"detected {len(filter_reads)} insertions where at least one end maps entirely to the reference genome, removing these (since they can not be chimeric)...")
    with gzip.open(combined_insertions, 'wt') as f:
        f.writelines([str(i) for i in insertions if i.name not in filter_reads])
    print(f'wrote {len([i for i in insertions if i.name not in filter_reads])} insertions to insertions.txt.gz')
    n_excluded = 0
    n_included = 0
    with gzip.open(insertions_genotyping_file, 'wt') as f:
        for i in insertions:
            if i.name in filter_reads:
                continue
            right = i.right_clipped.upper()[:CONFIG['genotyping']['max_bases']]
            left = revcomp(i.left_clipped[:CONFIG['genotyping']['max_bases']]).upper()
            if len(i.chr) > 5:
                print(f"Excluding {i.name} since the breakpoint is on {i.chr}.")
                n_excluded += 1
                continue
            right_ref = get_sequence(i.chr, i.right_pos,
                                     i.right_pos + CONFIG['genotyping']['max_bases']).upper()
            left_ref = get_sequence(i.chr, i.left_pos - CONFIG['genotyping']['max_bases'],
                                    i.left_pos).upper()
            consensus = get_sequence(i.chr, i.right_pos, i.left_pos).upper()
            if not len(right_ref) or not len(left_ref):
                print(f"Excluding {i.name} due to missing reference")
                n_excluded += 1
                continue
            if 'N' in right_ref:
                print(f"Excluding {i.name} due to Ns in right reference")
                n_excluded += 1
                continue
            if 'N' in left_ref:
                print(f"Excluding {i.name} due to Ns in left reference")
                n_excluded += 1
                continue
            if len(consensus) and 'N' in consensus:
                print(f"Excluding {i.name} due to Ns in the consensus sequence")
                n_excluded += 1
                continue
            n_included += 1
            f.write(f'>{i.name}\n')
            f.write(f'@RIGHT_INSERTION\n{right}\n')
            f.write(f'@LEFT_INSERTION\n{left}\n')
            f.write(f'@RIGHT_REFERENCE\n{right_ref}\n')
            f.write(f'@LEFT_REFERENCE\n{left_ref}\n')
    print(f"wrote {n_included} of {n_included + n_excluded}, excluded {n_excluded}")
