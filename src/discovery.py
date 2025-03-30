# Geminate - Genomic events of mutation through insertional alterations by transposable elements
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


# this script detects potential insertion points in a bam file
# an insertion point is defined as a pair of clipped reads facing each other within MAX_BP_WINDOW bases
# === -> matches to reference genome
# **** -> does not match reference genome
# all directions are relative to the (+) strand (bam convention)
# breakpoints:
# R is the position of the first base of the right-clipped read
# L is the position of the first base of the mapped part of a left-clipped read
#                R                                  L
#      ==========*****************....**************============
#        RIGHT_CLIPPED                    LEFT_CLIPPED
# a breakpoint is named as seqname:R-L
# for instance, if a 6 bp target site duplication of an ltr exists, this will be the situation relative to the ref genome
#                  R
#     =============********
#     *******=================
#            L
#            0123456
# thus L<R and L-R = -6
#
# this script outputs a special .txt.gz. This is a custom text based file format similar to a fasta file
# layout
# >[seqname]:[R]-[L]
# @RIGHT_CONSENSUS
# [RIGHT_MAPPED_SEQ][right_clipped_seq]
# @RIGHT_MATES
# one or several mates, separated by newline
# @LEFT_CONSENSUS
# [left_clipped_seq][LEFT_MAPPED_SEQ]
# @LEFT_MATES
# one or several mates, separated by newline

# Mate extraction:
#         +-strand 5'  ===================--------------------------------======================= 3'   => insertion
#                 R1++++++++++.......------------R2                                                 -> not interesting, mate gives no additional information
#                 R2++++++++++.......------------R1
#                                  R1++++++++++.......------------R2                                -> interesting, mate gives additional information and originally maps to the - strand
#                                  R2++++++++++.......------------R1
#                                                R1++++++++++.......-----------R2                   -> interesting, mate gives additional information and originally maps to the + strand
#                                                R2++++++++++.......-----------R1
#                                                                R1++++++++++.......-----------R2    -> not interesting, mate gives no additional information
#                                                                R2++++++++++.......-----------R1    -> not interesting, mate gives no additional information


import pysam
from typing import Tuple, Iterable, List
from revcomp import revcomp
import gzip
import sys
from adapter import is_adapter
from consensus import find_consensus, is_good_consensus
from sequence_checks import is_low_complexity, has_well_defined_breakpoint, clean_clipped_seq
from config import CONFIG

DEBUG = True #print debug messages

class Discovery:
    # constants
    CLIP_LEFT = 1
    CLIP_RIGHT = 2

    def __init__(self):
        # dicts holding temprary breakpoints per seqname
        self.bp_left = {}
        self.bp_right = {}
        # dict holding mates to be searched.
        self.mates = {}
        # dict holding the final breakpoints
        self.breakpoints = {}
        #current reference:
        self.reference_name = None
        #alignment file
        self.filepath = None
    def cleanup(self) -> None:
        """
        Perform cleanup of the current reference name
        """
        if self.reference_name is None: return
        assert self.filepath is not None
        if not len(self.bp_left) or not len(self.bp_right):
            DEBUG and print(f"seqname {self.reference_name} has {len(self.bp_left)} left bps and {len(self.bp_right)} right bps, skipping")
            return

        #dictionaries with "clean" breakpoints
        clean_left = {}
        clean_right = {}
        mates_left = {}
        mates_right = {}

        #cleaning statistics
        stats = {
            'removed_single': {
                'left': 0,
                'right': 0
            },
            'short_indel': {
                'left': 0,
                'right': 0
            },
            'removed_coverage': {
                'left': 0,
                'right': 0
            },
            'removed_no_consensus': {
                'left': 0,
                'right': 0
            },
            'removed_illdefined_breakpoint': {
                'left': 0,
                'right': 0
            },
        }

        #clean left breakpoints
        for breakpoint, bps in self.bp_left.items():
            if len(bps) >= CONFIG['discovery']['min_breakpoints_aggregated_during_first_step']:
                mates_left[breakpoint] = []
                # aggregate clipped sequences
                clipped_seqs = {}
                unclipped_seqs = {}
                is_excluded = False
                for qname, clipped_part, unclipped_part, is_read1, is_forward, exclude_flag in bps:
                    # clean clipped sequence from polyG and adapter sequences
                    if qname not in clipped_seqs.keys():
                        clipped_seqs[qname] = []
                        unclipped_seqs[qname] = []
                    clipped_seqs[qname].append(clean_clipped_seq(revcomp(clipped_part)))
                    unclipped_seqs[qname].append(unclipped_part)
                    if not is_forward:  # only reads originally mapping to the - strand can have a mate within the insertion
                        mates_left[breakpoint].append((self.reference_name, breakpoint, not is_read1, qname))
                    is_excluded = is_excluded or exclude_flag
                for qname, seqs in clipped_seqs.items():
                    clipped_seqs[qname] = sorted(seqs, key=len, reverse=True)[0]  # choose only longest seq
                for qname, seqs in unclipped_seqs.items():
                    unclipped_seqs[qname] = sorted(seqs, key=len, reverse=True)[0]
                # remove seqs that do not reach minimal length
                clipped_seqs = [clipped_seq for clipped_seq in clipped_seqs.values() if len(clipped_seq) > CONFIG['discovery']['min_clip_len']]
                unclipped_seqs = [unclipped_seq for unclipped_seq in unclipped_seqs.values() if
                                  len(unclipped_seq) > CONFIG['discovery']['min_clip_len']]
                if is_excluded:
                    stats['short_indel']['left'] += 1
                    continue
                # remove breakpoints with less than 2 independent reads
                if len(set([len(clipped_seq) for clipped_seq in clipped_seqs])) < 2:
                    stats['removed_single']['left'] += 1
                    continue  # ignore insertions with not at least three different length evidence reads
                if len(set([len(unclipped_seq) for unclipped_seq in unclipped_seqs])) < 2:
                    stats['removed_single']['left'] += 1
                    continue  # ignore insertions with not at least three different length evidence reads
                # ignore high-coverage regions
                if len(clipped_seqs) > 50:
                    stats['removed_coverage']['left'] += 1
                    continue
                if len(unclipped_seqs) > 50:
                    stats['removed_coverage']['left'] += 1
                    continue

                # find consensus for clipped part
                clipped_n, clipped_score, clipped_seq = find_consensus(clipped_seqs)
                unclipped_n, unclipped_score, unclipped_seq = find_consensus(unclipped_seqs)
                if not is_good_consensus(clipped_n, clipped_score):
                    stats['removed_no_consensus']['left'] += 1
                    continue
                if not is_good_consensus(unclipped_n, unclipped_score):
                    stats['removed_no_consensus']['left'] += 1
                    continue
                if not has_well_defined_breakpoint(clipped_seq, unclipped_seq):
                    stats['removed_illdefined_breakpoint']['left'] += 1
                    continue
                # print(f"L {breakpoint}: {revcomp(clipped_seq)} {unclipped_seq}")
                clean_left[breakpoint] = (revcomp(clipped_seq).lower(), unclipped_seq)
            else:
                stats['removed_single']['left'] += 1

        for breakpoint, bps in self.bp_right.items():
            mates_right[breakpoint] = []
            if len(bps) >= CONFIG['discovery']['min_breakpoints_aggregated_during_first_step']:
                # aggregate clipped sequences
                clipped_seqs = {}
                unclipped_seqs = {}
                is_excluded = False
                for qname, clipped_part, unclipped_part, is_read1, is_forward, exclude_flag in bps:
                    # clean clipped sequence from polyG and adapter sequences
                    if qname not in clipped_seqs.keys():
                        clipped_seqs[qname] = []
                        unclipped_seqs[qname] = []
                    clipped_seqs[qname].append(clean_clipped_seq(clipped_part))
                    unclipped_seqs[qname].append(revcomp(unclipped_part))
                    if is_forward:  # only reads originally mapping to the + strand can have a mate within the insertion
                        mates_right[breakpoint].append((self.reference_name, breakpoint, not is_read1, qname))
                    is_excluded = is_excluded or exclude_flag
                for qname, seqs in clipped_seqs.items():
                    clipped_seqs[qname] = sorted(seqs, key=len, reverse=True)[0]  # choose only longest seq
                for qname, seqs in unclipped_seqs.items():
                    unclipped_seqs[qname] = sorted(seqs, key=len, reverse=True)[0]
                # remove seqs that do not reach minimal length
                clipped_seqs = [clipped_seq for clipped_seq in clipped_seqs.values() if len(clipped_seq) > CONFIG['discovery']['min_clip_len']]
                unclipped_seqs = [unclipped_seq for unclipped_seq in unclipped_seqs.values() if
                                  len(unclipped_seq) > CONFIG['discovery']['min_clip_len']]
                if is_excluded:
                    stats['short_indel']['right'] += 1
                    continue
                # remove breakpoints with less than 2 independent reads
                if len(set([len(clipped_seq) for clipped_seq in clipped_seqs])) < 2:
                    stats['removed_single']['right'] += 1
                    continue  # ignore insertions with not at least three different length evidence reads
                if len(set([len(unclipped_seq) for unclipped_seq in unclipped_seqs])) < 2:
                    stats['removed_single']['right'] += 1
                    continue  # ignore insertions with not at least three different length evidence reads

                # ignore high-coverage regions
                if len(clipped_seqs) > 50:
                    stats['removed_coverage']['right'] += 1
                    continue
                if len(unclipped_seqs) > 50:
                    stats['removed_coverage']['right'] += 1
                    continue

                # find consensus for clipped part
                clipped_n, clipped_score, clipped_seq = find_consensus(clipped_seqs)
                unclipped_n, unclipped_score, unclipped_seq = find_consensus(unclipped_seqs)
                if not is_good_consensus(clipped_n, clipped_score):
                    stats['removed_no_consensus']['right'] += 1
                    continue
                if not is_good_consensus(unclipped_n, unclipped_score):
                    stats['removed_no_consensus']['right'] += 1
                    continue
                if not has_well_defined_breakpoint(unclipped_seq, clipped_seq):
                    stats['removed_illdefined_breakpoint']['right'] += 1
                    continue
                # print(f"R {breakpoint}: {revcomp(unclipped_seq)} {clipped_seq.lower()}")
                clean_right[breakpoint] = (clipped_seq.lower(), revcomp(unclipped_seq))
            else:
                stats['removed_single']['right'] += 1

        DEBUG and print(f"cleanup stats")
        DEBUG and print(
            f"right: kept={len(clean_right)}, removed single={stats['removed_single']['right']} high_cov={stats['removed_coverage']['right']}, no_consensus={stats['removed_no_consensus']['right']}, illdefined_bp={stats['removed_illdefined_breakpoint']['right']} short_indel={stats['short_indel']['right']}")
        DEBUG and print(
            f"left: kept={len(clean_left)}, removed single={stats['removed_single']['left']} high_cov={stats['removed_coverage']['left']}, no_consensus={stats['removed_no_consensus']['left']}, illdefined_bp={stats['removed_illdefined_breakpoint']['left']}  short_indel={stats['short_indel']['left']}")

        left_bps = sorted(clean_left.keys())
        right_bps = sorted(clean_right.keys())

        left_i = 0
        right_i = 0

        pairs = 0
        removed_nomate_left = 0
        removed_nomate_right = 0
        multibp_in_window = 0
        sequence_contained_in_eachother = 0
        high_bp_coverage = 0
        with (pysam.AlignmentFile(self.filepath) as f):
            while left_i < len(clean_left) and right_i < len(clean_right):
                delta = left_bps[left_i] - right_bps[right_i]
                if abs(delta) <= CONFIG['discovery']['max_bp_window']:

                    # make sure there is no other possibility in the vincinity
                    if left_i < len(clean_left) - 1:
                        if abs(left_bps[left_i + 1] - right_bps[right_i]) <= CONFIG['discovery']['max_bp_window']:
                            multibp_in_window += 1
                            left_i += 2
                            continue
                    if right_i < len(clean_right) - 1:
                        if abs(left_bps[left_i] - right_bps[right_i + 1]) <= CONFIG['discovery']['max_bp_window']:
                            multibp_in_window += 1
                            right_i += 2
                            continue

                    left_seq = clean_left[left_bps[left_i]][0]
                    left_mapped_seq = clean_left[left_bps[left_i]][1]
                    right_seq = clean_right[right_bps[right_i]][0]
                    right_mapped_seq = clean_right[right_bps[right_i]][1]
                    # check if clipped sequences are contained within each other (this is a sign for a microinsertion)
                    if left_seq[:6] in right_mapped_seq or right_seq[-6:] in left_mapped_seq:
                        left_i += 1
                        right_i += 1
                        sequence_contained_in_eachother += 1
                        continue

                    start = right_bps[right_i]
                    end = left_bps[left_i]
                    start, end = min(start, end), max(start, end)
                    coverage = f.count(self.reference_name, start, end + 1)
                    if coverage > CONFIG['discovery']['max_read_count']:
                        right_i += 1
                        left_i += 1
                        high_bp_coverage += 1
                        continue
                    self.breakpoints[f'{self.reference_name}:{right_bps[right_i]}-{left_bps[left_i]}'] = (
                        right_seq,
                        left_seq,
                        clean_right[right_bps[right_i]][1],
                        clean_left[left_bps[left_i]][1],
                        [],  # rescue fragments right
                        []  # rescue fragments left
                    )
                    for ref, bp, read1, qname in mates_right[right_bps[right_i]]:
                        if f"{qname}:{1 if read1 else 0}" in self.mates.keys():
                            if self.mates[
                                    f"{qname}:{1 if read1 else 2}"] == \
                                    f'{self.reference_name}:{right_bps[right_i]}-{left_bps[left_i]}:R':
                                continue
                            else:
                                self.mates[f"{qname}:{1 if read1 else 2}"] = 'invalid'
                        else:
                            self.mates[
                                f"{qname}:{1 if read1 else 2}"] = \
                                f'{self.reference_name}:{right_bps[right_i]}-{left_bps[left_i]}:R'

                    for ref, bp, read1, qname in mates_left[left_bps[left_i]]:
                        if f"{qname}:{1 if read1 else 0}" in self.mates.keys():
                            if self.mates[
                                    f"{qname}:{1 if read1 else 2}"] == f'{self.reference_name}:{right_bps[right_i]}-{left_bps[left_i]}:L':
                                continue
                            else:
                                self.mates[f"{qname}:{1 if read1 else 2}"] = 'invalid'
                        else:
                            self.mates[
                                f"{qname}:{1 if read1 else 2}"] = f'{self.reference_name}:{right_bps[right_i]}-{left_bps[left_i]}:L'
                    DEBUG and print(f"found novel insertion number {len(self.breakpoints)}")
                    DEBUG and print(f">{self.reference_name}:{right_bps[right_i]}-{left_bps[left_i]}:R")
                    DEBUG and print(right_seq)
                    DEBUG and print(f">{self.reference_name}:{right_bps[right_i]}-{left_bps[left_i]}:L")
                    DEBUG and print(left_seq)
                    DEBUG and print("")
                    right_i += 1
                    left_i += 1
                elif delta > CONFIG['discovery']['max_bp_window']:
                    right_i += 1
                elif delta < -CONFIG['discovery']['max_bp_window']:
                    left_i += 1
        DEBUG and print(
            f"paired stats: kept={len(self.breakpoints)}, nomate_left={removed_nomate_left}, nomate_right={removed_nomate_right}, coverage={high_bp_coverage}, multibp={multibp_in_window}, containing_seq={sequence_contained_in_eachother}")

        self.bp_left = {}
        self.bp_right = {}
        DEBUG and print(f"... currently holding {len(self.mates)} rescue fragments.")


    def add_breakpoint(self, side: int, breakpoint: int, qname: str, clipped_seq: str, unclipped_seq: str, is_read1: bool,
                       is_forward: bool, exclude_flag: bool) -> None:  # this function has side effects
        """
        side = self.CLIP_LEFT (clipped part on left hand side, usually 3' end of insertion)
             = self.CLIP_RIGHT (right hand side, usually 5' end of insertion)
        breakpoint: position of the first un-clipped base for self.CLIP_LEFT or the last clipped base of self.CLIP_RIGHT.
        qname: Name of the fragment
        clipped_seq: clipped sequence
        unclipped_seq: unclipped sequence.
            NOTE: The original sequence is:
                self.CLIP_LEFT: [cipped_seq][unclipped_seq]
                self.CLIP_RIGHT: [unclipped_seq][clipped_seq]
        is_read1: boolean indicating of this is a read 1 fragment or a read 2 fragment
        is_forward: boolean indicating if this read ORIGINALLY (not presently) mapped to the fwd strand.

        """
        clipped_seq = clipped_seq.upper()
        unclipped_seq = unclipped_seq.upper()
        # check if sequence is too short
        if len(clipped_seq) < CONFIG['discovery']['min_clip_len']:
            return
        if side is self.CLIP_LEFT:
            # check if unclipped part has homopolymer at breakpoint
            if 'T' * CONFIG['discovery']['max_homopolymer_len'] in unclipped_seq[:CONFIG['discovery']['max_homopolymer_len']] or \
                    'A' * CONFIG['discovery']['max_homopolymer_len'] in unclipped_seq[:CONFIG['discovery']['max_homopolymer_len']] or \
                    'C' * CONFIG['discovery']['max_homopolymer_len'] in unclipped_seq[:CONFIG['discovery']['max_homopolymer_len']] or \
                    'G' * CONFIG['discovery']['max_homopolymer_len'] in unclipped_seq[:CONFIG['discovery']['max_homopolymer_len']]:
                return
            if breakpoint not in self.bp_left.keys():
                self.bp_left[breakpoint] = []
            self.bp_left[breakpoint].append((qname, clipped_seq, unclipped_seq, is_read1, is_forward, exclude_flag))
        else:
            # check if unclipped part has homopolymer at breakpoint
            if 'T' * CONFIG['discovery']['max_homopolymer_len'] in unclipped_seq[-CONFIG['discovery']['max_homopolymer_len']:] or \
                    'A' * CONFIG['discovery']['max_homopolymer_len'] in unclipped_seq[-CONFIG['discovery']['max_homopolymer_len']:] or \
                    'C' * CONFIG['discovery']['max_homopolymer_len'] in unclipped_seq[-CONFIG['discovery']['max_homopolymer_len']:] or \
                    'G' * CONFIG['discovery']['max_homopolymer_len'] in unclipped_seq[-CONFIG['discovery']['max_homopolymer_len']:]:
                return
            if breakpoint not in self.bp_right.keys():
                self.bp_right[breakpoint] = []
            self.bp_right[breakpoint].append((qname, clipped_seq, unclipped_seq, is_read1, is_forward, exclude_flag))


    def extract_chimeric(self) -> None:  # this function has side effects!
        assert self.filepath is not None
        with pysam.AlignmentFile(self.filepath) as f:
            self.reference_name = None
            for read in f:
                # determine if this read is a clipped read
                if read.mapq < CONFIG['discovery']['min_mapq']: continue
                if read.is_secondary: continue
                if read.is_qcfail: continue
                if read.is_duplicate: continue
                if read.is_unmapped: continue
                if read.reference_name != self.reference_name:
                    self.cleanup()
                    self.reference_name = read.reference_name
                if len(read.reference_name) > 5:
                    continue  # do not process weird contigs
                if read.reference_name == "MT" or read.reference_name == "chrM":
                    continue  # ignore mitochondrial reads
                try:
                    left_class, left_len = read.cigartuples[0]
                    right_class, right_len = read.cigartuples[-1]
                except TypeError as e:
                    print(e)
                clip = None
                # find cruciform dna artefacts here

                if left_class == pysam.CSOFT_CLIP and not right_class == pysam.CSOFT_CLIP:
                    clip = self.CLIP_LEFT
                elif right_class == pysam.CSOFT_CLIP and not left_class == pysam.CSOFT_CLIP:
                    clip = self.CLIP_RIGHT
                elif right_class == pysam.CSOFT_CLIP and left_class == pysam.CSOFT_CLIP:
                    if left_len > right_len:
                        clip = self.CLIP_RIGHT
                    elif right_len > left_len:
                        clip = self.CLIP_LEFT
                if clip:
                    exclude_flag = False
                    # check if this is a cruciform DNA artefact or a short indel. If yes, remove
                    if read.is_supplementary:
                        for sa in read.get_tag('SA').split(";"):
                            if not len(sa): continue
                            seqname, start, strand, cigar, mapq, _ = sa.split(",")
                            if seqname==read.reference_name:
                                if abs(read.reference_start-int(start))<CONFIG['discovery']['exclude_same_contig_supplementary']:
                                    exclude_flag = True
                    if clip is self.CLIP_LEFT and left_len >= CONFIG['discovery']['min_clip_len']:
                        if is_adapter(revcomp(read.seq[:left_len][-CONFIG['discovery']['min_clip_len']:])):
                            continue
                        # remove right-clipped part of mapped read, if necessary
                        if right_class is pysam.CSOFT_CLIP:
                            self.add_breakpoint(clip, read.reference_start, read.query_name, read.seq[:left_len],
                                           read.seq[left_len:-right_len], read.is_read1, read.is_forward, exclude_flag)
                        else:
                            self.add_breakpoint(clip, read.reference_start, read.query_name, read.seq[:left_len],
                                           read.seq[left_len:], read.is_read1, read.is_forward, exclude_flag)
                    elif clip is self.CLIP_RIGHT and right_len >= CONFIG['discovery']['min_clip_len']:
                        if is_adapter(read.seq[-right_len:][:CONFIG['discovery']['min_clip_len']]):
                            continue
                        # remove left-clipped part of mapped read, if necessary
                        if left_class is pysam.CSOFT_CLIP:
                            self.add_breakpoint(clip, read.reference_end, read.query_name, read.seq[-right_len:],
                                           read.seq[left_len:-right_len], read.is_read1, read.is_forward, exclude_flag)
                        else:
                            # add_clip_stats(read.seq[-right_len:][:CONFIG['discovery']['min_clip_len']])
                            self.add_breakpoint(clip, read.reference_end, read.query_name, read.seq[-right_len:],
                                           read.seq[:-right_len], read.is_read1, read.is_forward, exclude_flag)
            # final cleanup
            self.cleanup()


    def cleanup_mates(self) -> None:  # this function has side effects!
        fragments_to_be_deleted = set()
        for f in self.mates.keys():
            other_side = '2' if f[-1:] == '1' else '1'
            if f[:-1] + other_side in self.mates.keys():
                fragments_to_be_deleted.add(f)
                fragments_to_be_deleted.add(f[:-1] + other_side)
        DEBUG and print(
            f"removing {len(fragments_to_be_deleted)} mates that have already been used to initially identify the breakpoint.")
        for i in fragments_to_be_deleted:
            del self.mates[i]
        DEBUG and print(f"left with {len(self.mates)} mates that have to be extracted.")


    def find_mates(self) -> None:  # this function has side effects!
        """
        Find mates listed in self.mates
        """
        with pysam.AlignmentFile(self.filepath) as f:
            for read in f:
                if read.is_secondary: continue
                if read.is_qcfail: continue
                if read.is_duplicate: continue
                if read.is_supplementary: continue
                qname = f"{read.query_name}:{1 if read.is_read1 else 2}"
                if qname in self.mates.keys():
                    breakpoint = self.mates[qname]
                    if breakpoint == 'invalid': continue
                    side = self.CLIP_LEFT if breakpoint[-1] == 'L' else self.CLIP_RIGHT
                    breakpoint = breakpoint[:-2]
                    seq = read.seq
                    # check if the original sequence is inverted in the bam file.
                    # this is only the case if the sequence is mapped to the forward strand xor is read 2.
                    seq = clean_clipped_seq(seq)
                    if side is self.CLIP_RIGHT:
                        # the mate of a read to self.CLIP_RIGHT originally always maps to the - strand. Thus is has to be inverted, except if the bam file has already done this for us.
                        if read.is_forward:
                            self.breakpoints[breakpoint][4].append(revcomp(
                                seq))  # , f";revcomp=T, side=R, mapped={'no' if read.is_unmapped else 'yes'} dir={'fwd' if read.is_forward else 'rev'}, read={'1' if read.is_read1 else '2'}, align={read.reference_name}:{read.reference_start}-{read.reference_end} mapq={read.mapq} cigar={read.cigarstring} qname={read.query_name}"))
                        else:
                            self.breakpoints[breakpoint][4].append(
                                seq)  # , f";revcomp=F, side=R, mapped={'no' if read.is_unmapped else 'yes'} dir={'fwd' if read.is_forward else 'rev'}, read={'1' if read.is_read1 else '2'}, align={read.reference_name}:{read.reference_start}-{read.reference_end} mapq={read.mapq} cigar={read.cigarstring} qname={read.query_name}"))
                    else:  # side is self.CLIP_LEFT
                        # the mate of a read to self.CLIP_LEFT originally always maps to the + strand. Thus is has to be inverted, if the bam file has done this for us.
                        if read.is_forward:
                            self.breakpoints[breakpoint][5].append(
                                seq)  # , f";revcomp=F, side=L, mapped={'no' if read.is_unmapped else 'yes'} dir={'fwd' if read.is_forward else 'rev'}, read={'1' if read.is_read1 else '2'}, align={read.reference_name}:{read.reference_start}-{read.reference_end} mapq={read.mapq} cigar={read.cigarstring} qname={read.query_name}"))
                        else:
                            self.breakpoints[breakpoint][5].append(revcomp(
                                seq))  # , f";revcomp=T, side=L, mapped={'no' if read.is_unmapped else 'yes'} dir={'fwd' if read.is_forward else 'rev'}, read={'1' if read.is_read1 else '2'}, align={read.reference_name}:{read.reference_start}-{read.reference_end} mapq={read.mapq} cigar={read.cigarstring} qname={read.query_name}"))


    def extend_mates(self) -> None:  # this function has side effects
        for breakpoint, (
                right_clipped, left_clipped, right_matched, left_matched, right_rescues,
                left_rescues) in self.breakpoints.items():
            # extend right
            while len(right_rescues):
                search_seq = right_clipped[-min(len(right_clipped), 10):].upper()
                carry_over = []
                extension = []
                for r in right_rescues:
                    count = r.count(search_seq)
                    if count == 0:
                        carry_over.append(r)
                    elif count == 1:
                        end_base = r.find(search_seq) + len(search_seq)
                        extension.append(r[end_base:])
                    else:  # multi-hit: repetitive
                        extension = None
                        break
                if extension is None:
                    break
                right_rescues = carry_over
                if not len(extension):
                    break
                extension = sorted(extension, key=len, reverse=True)[0]
                if not len(extension):
                    break
                right_clipped += extension.lower()
                #DEBUG and print(f"{breakpoint}: extended right by {extension.lower()} to {right_clipped}")

            while len(left_rescues):
                search_seq = left_clipped[:min(len(left_clipped), 10)].upper()
                carry_over = []
                extension = []
                for r in left_rescues:
                    count = r.count(search_seq)
                    if count == 0:
                        carry_over.append(r)
                    elif count == 1:
                        start_base = r.find(search_seq)
                        extension.append(r[:start_base])
                    else:  # multi-hit: repetitive
                        extension = None
                        break
                if extension is None:
                    break
                left_rescues = carry_over
                if not len(extension):
                    break
                extension = sorted(extension, key=len, reverse=True)[0]
                if not len(extension):
                    break
                left_clipped = extension.lower() + left_clipped
                #DEBUG and print(f"{breakpoint}: extended left by {extension.lower()} to {left_clipped}")

            self.breakpoints[breakpoint] = (
                right_clipped, left_clipped, right_matched, left_matched, right_rescues, left_rescues)

    def discovery(self, bamfile_name: str):
        self.filepath = bamfile_name
        print(f"starting discovery in file {bamfile_name}")
        self.extract_chimeric()
        print(f"finding mates of reads covering breakpoints in file {bamfile_name}")
        self.cleanup_mates()
        self.find_mates()
        print("extend clipped sequences using the mates sequences")
        self.extend_mates()
        print(f"disovery phase completed for file {bamfile_name}, found {len(self.breakpoints)} breakpoints.")

    def output(self, outfile_path: str):
        print(f"writing result to {outfile_path}")
        with gzip.open(outfile_path, 'wt') as outfile:
            for breakpoint, (
                    right_clipped, left_clipped, right_matched, left_matched, right_rescues,
                    left_rescues) in self.breakpoints.items():
                outfile.write(f">{breakpoint}\n")
                outfile.write(f"@RIGHT_CONSENSUS\n")
                outfile.write(f"{right_matched}{right_clipped}\n")
                outfile.write(f"@RIGHT_MATES\n")
                outfile.writelines([f"{s}\n" for s in right_rescues])
                outfile.write(f"@LEFT_CONSENSUS\n")
                outfile.write(f"{left_clipped}{left_matched}\n")
                outfile.write("@LEFT_MATES\n")
                outfile.writelines([f"{s}\n" for s in left_rescues])
                outfile.write("\n")
        print("done.")