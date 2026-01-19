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
import os
from revcomp import revcomp
import pysam
import pyliftover
from config import CONFIG
from combine_insertions_insertion import Insertion, TYPE_LEFT_POLYA, TYPE_RIGHT_POLYA, TYPE_FULL_INFO
from combine_insertions_intersect_insertions import intersect_insertions
from combine_insertions_get_sequence import get_sequence
from sequence_checks import sequence_matching_score
from collections import Counter

def combine_insertions(input_files, insertions_genotyping_file, combined_insertions, insertions_fasta, insertion_bam, threads):

    all_insertions = []
    j_cutoff = 0
    for f in input_files:
        j_cutoff += 1
        #if j_cutoff > 5: break
        fb = os.path.basename(f)
        file_insertions = [i for i in Insertion.parseFile(f) if len(i.reference_name) < 6 and i.reference_name != "MT" and i.reference_name !="chrM"]
        print(f"\033[37mFile {fb}: imported {len(file_insertions)} insertions\033[0m")
        if len(file_insertions) > CONFIG['combine_insertions']['exclude_files_with_many_insertions']:
            print(f"\033[31mremoved file {f} since it contains too many insertions.\033[0m")
        else:
            all_insertions += file_insertions
    print(f"intersecting insertions from {len(input_files)} files...")
    insertions = intersect_insertions(all_insertions)
    #remove insertions in regions with far too high count
    bins = []
    bin_range = 100
    ins_cutoff = 4
    for i in insertions:
        bins.append((i.reference_name,int(i.right_pos/bin_range)*bin_range if i.right_pos is not None else int(i.left_pos/bin_range)*bin_range))
    c = Counter(bins)
    removed = 0
    regions = 0
    for (reference_name, bin), n in c.items():
        clean_ins = []
        if n>=ins_cutoff:
            regions += 1
            for i in insertions:
                pos = i.right_pos if i.right_pos is not None else i.left_pos
                if i.reference_name == reference_name and pos > bin - bin_range/2 and pos < bin + bin_range * 1.5:
                    removed += 1
                else:
                    clean_ins.append(i)
            insertions = clean_ins
    print(f"filtering regions with very high insertion rate of {ins_cutoff} or higher per {bin_range} bases ,removed {removed} insertions in {regions} regions, {len(insertions)} insertions are remaining")
    print(f"writing summarised insertions fasta file {insertions_fasta}")
    with gzip.open(insertions_fasta, 'wt') as f:
        f.writelines(
            [f'{i.left_consensus.fastq(f"{i.name}:L")}{i.right_consensus.fastq(f"{i.name}:R")}' for i in insertions if i.type is TYPE_FULL_INFO])
        f.writelines(
            [f'{i.left_consensus.fastq(f"{i.name}:L")}' for i in insertions if i.type is TYPE_RIGHT_POLYA])
        f.writelines(
            [f'{i.right_consensus.fastq(f"{i.name}:R")}' for i in insertions if i.type is TYPE_LEFT_POLYA])
    if not os.path.exists(insertion_bam):
        print(f"running bowtie2 {CONFIG['combine_insertions']['bowtie2_executable']} with index {CONFIG['combine_insertions']['bowtie2_index']}")
        cmd = f"{CONFIG['combine_insertions']['bowtie2_executable']} {insertions_fasta} -x {CONFIG['combine_insertions']['bowtie2_index']} --end-to-end --sensitive --threads {threads} --qc-filter | {CONFIG['combine_insertions']['samtools_executable']} view -F 4 -b -o {insertion_bam}"
        os.system(cmd)
    filter_reads = set()
    with pysam.AlignmentFile(insertion_bam) as f:
        for read in f:
            filter_reads.add(read.query_name[:-2])
    print(
        f"detected {len(filter_reads)} insertions where at least one end maps entirely to the reference genome, removing these (since they can not be chimeric)...")
    insertions = [i for i in insertions if i.name not in filter_reads]

    print(f"now re-mapping in local mode all clipped parts of reads")
    clipped_bam = insertion_bam[:-4] + ".insertionsonly.bam"
    if not os.path.exists(clipped_bam):
        with gzip.open(insertions_fasta, 'wt') as f:
            f.writelines(
                [f'{i.left_clipped.fastq(f"{i.name}:L")}{i.right_clipped.fastq(f"{i.name}:R")}' for i in insertions])
        print(f"running bowtie2 {CONFIG['combine_insertions']['bowtie2_executable']} with index {CONFIG['combine_insertions']['bowtie2_index2']}")
        cmd = f"{CONFIG['combine_insertions']['bowtie2_executable']} {insertions_fasta} -k 1000 -x {CONFIG['combine_insertions']['bowtie2_index2']} --local --very-fast --threads {threads} --qc-filter | {CONFIG['combine_insertions']['samtools_executable']} view -F 4 -b -o {clipped_bam}"
        os.system(cmd)
    filter_reads = set()
    print(f"removing all reads where one of the clipped ends maps within 1000bp of the breakpoint.")
    lo = pyliftover.LiftOver(CONFIG['combine_insertions']['bowtie2_index2_lo'])
    filter_reads = set()
    delta_sampler = []
    with pysam.AlignmentFile(clipped_bam) as f:
        for read in f:
            if read.is_mapped:
                reference_name, pos, side = read.query_name.split(":")
                if len(reference_name) < 3 or reference_name[:3] != 'chr':
                    reference_name = f"chr{reference_name}"
                left_pos, right_pos = pos.split("-")
                if side == 'R':
                    pos = int(right_pos)
                else:
                    pos = int(left_pos)
                map = lo.convert_coordinate(read.reference_name, read.reference_start if read.is_forward else read.reference_end)
                if map is not None:
                    for map_rn, map_coord, map_strand, map_len in map:
                        if reference_name == map_rn:
                            if abs(pos-map_coord) < 1000:
                                delta_sampler.append(abs(pos-map_coord))
                                filter_reads.add(read.query_name[:-2])
    #print(Counter(delta_sampler))
    print(f"detected {len(filter_reads)} insertions where the clipped part maps near the breakpoint. Removing these")
    insertions = [i for i in insertions if i.name not in filter_reads]
    with gzip.open(combined_insertions, 'wt') as f:
        f.writelines(
            [f'{i.left_consensus.fastq(f"{i.name}:L")}{i.right_consensus.fastq(f"{i.name}:R")}' for i in insertions if
             i.type is TYPE_FULL_INFO])
        f.writelines(
            [f'{i.left_consensus.fastq(f"{i.name}:L")}' for i in insertions if i.type is TYPE_RIGHT_POLYA])
        f.writelines(
            [f'{i.right_consensus.fastq(f"{i.name}:R")}' for i in insertions if i.type is TYPE_LEFT_POLYA])
    print(f'wrote {len([i for i in insertions if i.name not in filter_reads])} insertions to insertions.txt.gz')
    n_excluded = 0
    n_included = 0
    with gzip.open(insertions_genotyping_file, 'wt') as f:
        for i in insertions:
            if i.name in filter_reads:
                continue
            right = i.right_clipped.upper()[:CONFIG['genotyping']['max_bases']]
            left = i.left_clipped[:CONFIG['genotyping']['max_bases']].upper().revcomp()
            if len(i.reference_name) > 5:
                print(f"Excluding {i.name} since the breakpoint is on {i.reference_name}.")
                n_excluded += 1
                continue
            if i.type is not TYPE_RIGHT_POLYA:
                right_ref = get_sequence(i.reference_name, i.right_pos,
                                         i.right_pos + CONFIG['genotyping']['max_bases']).upper()
            else:
                right_ref = None
            if i.type is not TYPE_LEFT_POLYA:
                left_ref = get_sequence(i.reference_name, i.left_pos - CONFIG['genotyping']['max_bases'],
                                    i.left_pos).upper()
            else:
                left_ref = None
            if i.type is not TYPE_RIGHT_POLYA and not len(right_ref):
                print(f"Excluding {i.name} due to missing reference")
                n_excluded += 1
                continue
            if i.type is not TYPE_LEFT_POLYA and not len(left_ref):
                print(f"Excluding {i.name} due to missing reference")
                n_excluded += 1
                continue
            if right_ref is not None and 'N' in right_ref:
                print(f"Excluding {i.name} due to Ns in right reference")
                n_excluded += 1
                continue
            if left_ref is not None and 'N' in left_ref:
                print(f"Excluding {i.name} due to Ns in left reference")
                n_excluded += 1
                continue
            if i.type is TYPE_RIGHT_POLYA  is not None and sequence_matching_score([left_ref, str(left)]) > 0:
                print(f"Excluding {i.name} due to similar left sequence between ref and alt {left_ref} vs {left}")
                n_excluded += 1
                continue
            if i.type is TYPE_LEFT_POLYA and sequence_matching_score([right_ref, str(right)]) > 0:
                print(f"Excluding {i.name} due to similar right sequence between ref and alt: {right_ref} vs {right}")
                n_excluded += 1
                continue
            if i.type is TYPE_FULL_INFO and sequence_matching_score([right_ref, str(right)]) > 0 and sequence_matching_score([left_ref, str(left)]) > 0:
                print(f"Excluding {i.name} due to similar right and left sequence between ref and alt: {right_ref} vs {right} (right) and {left_ref} vs {left} (left)")
                n_excluded += 1
                continue
            if i.type is TYPE_FULL_INFO and (left_ref == revcomp(right_ref) or left_ref == right_ref):
                print(f"Excluding {i.name} due to identical clipped sequences: {left_ref} vs {right_ref}")
                n_excluded += 1
                continue
            n_included += 1
            f.write(f'>{i.name}\n')
            if i.type is not TYPE_RIGHT_POLYA:
                f.write(f'@RIGHT_INSERTION\n{right}\n')
                f.write(f'@RIGHT_REFERENCE\n{right_ref}\n')
            if i.type is not TYPE_LEFT_POLYA:
                f.write(f'@LEFT_INSERTION\n{left}\n')
                f.write(f'@LEFT_REFERENCE\n{left_ref}\n')
    print(f"wrote {n_included} of {n_included + n_excluded}, excluded {n_excluded}")
