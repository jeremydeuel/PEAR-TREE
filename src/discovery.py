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
from revcomp import revcomp
import gzip
from adapter import is_adapter
from consensus import find_consensus, is_good_consensus
from sequence_checks import has_well_defined_breakpoint, clean_clipped_seq
from config import CONFIG
from quality_seq import QualitySeq
from breakpoint import Breakpoint, CLIP_LEFT, CLIP_RIGHT
from typing import Tuple, List, Dict
from polyABreakpoint import PolyABreakpoint
DEBUG = True #print debug messages

class Discovery:
    def __init__(self):
        # dicts holding temprary breakpoints per seqname
        self.temporary_breakpoints : List[Breakpoint] = []
        self.final_left_breakpoints : List[Breakpoint] = []
        self.final_right_breakpoints: List[Breakpoint] = []
        self.polyA : List[PolyABreakpoint] = []
        # dict holding mates to be searched.
        self.mates : Dict[str, PolyABreakpoint | Breakpoint ] =  {}
        # dict holding the final breakpoints
        #current reference:
        self.reference_name : str | None = None
        #alignment file
        self.filepath : str | None = None
    def cleanup(self) -> None:
        """
        Perform cleanup of the current reference name
        """
        if self.reference_name is None: return
        assert self.filepath is not None
        if not len(self.temporary_breakpoints):
            DEBUG and print(f"seqname {self.reference_name} has no breakpoints, skipping")
            return
        left_bps = [b for b in self.temporary_breakpoints if b.side is CLIP_LEFT]
        right_bps = [b for b in self.temporary_breakpoints if b.side is CLIP_RIGHT]
        print(f"doing cleanup with {len(self.temporary_breakpoints)} breakpoints, thereof {len(left_bps)} left and {len(right_bps)} right.")
        left_bps = sorted(left_bps, key=lambda x: x.breakpoint)
        right_bps = sorted(right_bps, key=lambda x: x.breakpoint)
        left_clean_bps = []
        right_clean_bps = []
        # join bps together
        last_bps = []
        for bp in left_bps:
            if len(last_bps):
                if abs(bp.breakpoint-last_bps[-1].breakpoint) < 6:
                    last_bps.append(bp)
                    continue
                else:
                    if b := Breakpoint.join(last_bps): left_clean_bps.append(b)
            last_bps = [bp]
        if len(last_bps):
            if b := Breakpoint.join(last_bps): left_clean_bps.append(b)
        last_bps = []
        for bp in right_bps:
            if len(last_bps):
                if abs(bp.breakpoint-last_bps[-1].breakpoint) < 6:
                    last_bps.append(bp)
                    continue
                else:
                    # cleanup last_bp
                    if b:= Breakpoint.join(last_bps): right_clean_bps.append(b)
            last_bps = [bp]
        if len(last_bps):
            if b:= Breakpoint.join(last_bps): right_clean_bps.append(b)
        print(f"left={len(left_clean_bps)}, right={len(right_clean_bps)}")
        self.temporary_breakpoints = []
        self.final_left_breakpoints += left_clean_bps
        self.final_right_breakpoints += right_clean_bps
        print(Breakpoint.stats)

    def add_breakpoint(self, side: int, breakpoint: int, qname: str, clipped_seq: str, unclipped_seq: str, is_read1: bool,
                       is_forward: bool, exclude_flag: bool, is_precise: bool = True) -> None:  # this function has side effects
        """
        side = CLIP_LEFT (clipped part on left hand side, usually 3' end of insertion)
             = CLIP_RIGHT (right hand side, usually 5' end of insertion)
        breakpoint: position of the first un-clipped base for CLIP_LEFT or the last clipped base of CLIP_RIGHT.
        qname: Name of the fragment
        clipped_seq: clipped sequence
        unclipped_seq: unclipped sequence.
            NOTE: The original sequence is:
                CLIP_LEFT: [cipped_seq][unclipped_seq]
                CLIP_RIGHT: [unclipped_seq][clipped_seq]
        is_read1: boolean indicating of this is a read 1 fragment or a read 2 fragment
        is_forward: boolean indicating if this read ORIGINALLY (not presently) mapped to the fwd strand.

        """
        #if len(clipped_seq) < CONFIG['discovery']['min_clip_len']:
        #    return
        bp = Breakpoint(side, self.reference_name, breakpoint, qname, clipped_seq.upper(), unclipped_seq.upper(), is_read1, is_forward, exclude_flag)
        if not is_precise: #discordant mate
            bp.bp_precise = False
        self.temporary_breakpoints.append(bp)
        return


    def extract_chimeric(self) -> None:  # this function has side effects!
        assert self.filepath is not None
        with pysam.AlignmentFile(self.filepath) as f:
            self.reference_name = None
            for read in f:
                # determine if this read is a clipped read
                if read.mapping_quality < CONFIG['discovery']['min_mapq']:
                    if b := PolyABreakpoint.findPolyA(read):
                        self.polyA.append(b)
                        #print(f"found valid polyA in mapq={read.mapq} with sequence {str(b.clipped)}")
                    continue
                if read.is_secondary: continue
                if read.is_qcfail: continue
                if read.is_duplicate: continue
                if read.reference_name != self.reference_name:
                    print(f'Cleanup for {self.reference_name}')
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
                    clip = CLIP_LEFT
                elif right_class == pysam.CSOFT_CLIP and not left_class == pysam.CSOFT_CLIP:
                    clip = CLIP_RIGHT
                elif right_class == pysam.CSOFT_CLIP and left_class == pysam.CSOFT_CLIP:
                    if left_len > right_len:
                        clip = CLIP_LEFT
                    elif right_len > left_len:
                        clip = CLIP_RIGHT
                if clip:
                    exclude_flag = False
                    # check if this is a cruciform DNA artefact or a short indel. If yes, remove
                    if read.is_supplementary:
                        for sa in read.get_tag('SA').split(";"):
                            if not len(sa): continue
                            reference_name, start, strand, cigar, mapq, _ = sa.split(",")
                            if reference_name == read.reference_name:
                                if abs(read.reference_start-int(start))<CONFIG['discovery']['exclude_same_contig_supplementary']:
                                    exclude_flag = True
                    if clip is CLIP_LEFT:
                        if is_adapter(revcomp(read.query_sequence[:left_len][-CONFIG['discovery']['min_clip_len']:])):
                            continue
                        # remove right-clipped part of mapped read, if necessary

                        if right_class is pysam.CSOFT_CLIP:
                            self.add_breakpoint(clip, read.reference_start, read.query_name, \
                                                QualitySeq(read.query_sequence[:left_len], read.query_qualities[:left_len]), \
                                                QualitySeq(read.query_sequence[left_len:-right_len],read.query_qualities[left_len:-right_len]), \
                                                read.is_read1, read.is_forward, exclude_flag)
                        else:
                            self.add_breakpoint(clip, read.reference_start, read.query_name, \
                                                QualitySeq(read.query_sequence[:left_len],read.query_qualities[:left_len]), \
                                           QualitySeq(read.query_sequence[left_len:],read.query_qualities[left_len:]),\
                                                read.is_read1, read.is_forward, exclude_flag)
                    elif clip is CLIP_RIGHT:
                        if is_adapter(read.query_sequence[-right_len:][:CONFIG['discovery']['min_clip_len']]):
                            continue
                        # remove left-clipped part of mapped read, if necessary
                        if left_class is pysam.CSOFT_CLIP:
                            self.add_breakpoint(clip, read.reference_end, read.query_name, \
                                                QualitySeq(read.query_sequence[-right_len:],read.query_qualities[-right_len:]), \
                                           QualitySeq(read.query_sequence[left_len:-right_len], read.query_qualities[left_len:-right_len]), \
                                                read.is_read1, read.is_forward, exclude_flag)
                        else:
                            # add_clip_stats(read.seq[-right_len:][:CONFIG['discovery']['min_clip_len']])
                            self.add_breakpoint(clip, read.reference_end, read.query_name, \
                                                QualitySeq(read.query_sequence[-right_len:],read.query_qualities[-right_len:]),
                                           QualitySeq(read.query_sequence[:-right_len],read.query_qualities[:-right_len]), \
                                                read.is_read1, read.is_forward, exclude_flag)
                else: #check if there is a discordant mate
                    continue #dont do any of this
                    exclude_flag = False #no clipped bases, thus also not showing evidence of a cruciform dna piece
                    if not read.is_proper_pair:
                        if read.is_forward:
                            if read.is_read1:
                                self.add_breakpoint(CLIP_RIGHT, read.reference_end, read.query_name,
                                                    QualitySeq('',[]),
                                                    QualitySeq(read.seq, read.query_qualities),
                                                    read.is_read1, read.is_forward, exclude_flag, False)
                            else:
                                self.add_breakpoint(CLIP_LEFT, read.reference_start, read.query_name,
                                                    QualitySeq(read.seq, read.query_qualities),
                                                    QualitySeq('', []),
                                                    read.is_read1, read.is_forward, exclude_flag, False)
                        else:
                            if read.is_read1:
                                self.add_breakpoint(CLIP_LEFT, read.reference_end, read.query_name,
                                                    QualitySeq(read.seq, read.query_qualities),
                                                    QualitySeq('', []),
                                                    read.is_read1, read.is_forward, exclude_flag, False)
                            else:
                                self.add_breakpoint(CLIP_RIGHT, read.reference_start, read.query_name,
                                                    QualitySeq('', []),
                                                    QualitySeq(read.seq, read.query_qualities),
                                                    read.is_read1, read.is_forward, exclude_flag, False)
            # final cleanup
            self.cleanup()


    def get_mates(self) -> Tuple[set, set, dict]:
        """extract a list of mates"""
        read1_bp = set()
        read2_bp = set()
        qmap = dict()
        for breakpoint in self.final_left_breakpoints:
            for read_1, qname in breakpoint.mates:
                if read_1:
                    read1_bp.add(qname)
                else:
                    read2_bp.add(qname)
                qmap[qname] = breakpoint
        for breakpoint in self.final_right_breakpoints:
            for read_1, qname in breakpoint.mates:
                if read_1:
                    read1_bp.add(qname)
                else:
                    read2_bp.add(qname)
                qmap[qname] = breakpoint
        for polyA in self.polyA:
            if polyA.qname in qmap.keys():
                continue #real breakpoints have priority
            qmap[polyA.qname] = polyA
            if polyA.is_read1:
                read2_bp.add(polyA.qname)
            else:
                read1_bp.add(polyA.qname)
        bp_union = read1_bp & read2_bp
        read1_bp = read1_bp - bp_union
        read2_bp = read2_bp - bp_union
        DEBUG and print(f"found {len(read1_bp)+len(read2_bp)} mates to be extracted, excluded {len(bp_union)} mates that themselves were already used as chimeric reads.")
        return read1_bp, read2_bp, qmap


    def find_mates(self) -> None:  # this function has side effects!
        """
        Find mates listed in self.mates
        """

        read1_mates, read2_mates, qmap = self.get_mates()

        with pysam.AlignmentFile(self.filepath) as f:
            for read in f:
                if read.is_secondary: continue
                if read.is_qcfail: continue
                if read.is_duplicate: continue
                if read.is_read1:
                    if not read.qname in read1_mates:
                        continue
                else:
                    if not read.qname in read2_mates:
                        continue
                bp = qmap[read.qname]
                seq = clean_clipped_seq(QualitySeq(read.seq, read.query_qualities))
                if isinstance(bp, PolyABreakpoint):
                    bp.setMate(read)
                elif isinstance(bp, Breakpoint):
                    if bp.side is CLIP_RIGHT:
                        # the mate of a read to CLIP_RIGHT originally always maps to the + strand. Thus is has to be inverted, except if the bam file has already done this for us.
                        if read.is_forward:
                            bp.mate_seqs.append(seq)  # , f";revcomp=T, side=R, mapped={'no' if read.is_unmapped else 'yes'} dir={'fwd' if read.is_forward else 'rev'}, read={'1' if read.is_read1 else '2'}, align={read.reference_name}:{read.reference_start}-{read.reference_end} mapq={read.mapq} cigar={read.cigarstring} qname={read.query_name}"))
                        else:
                            bp.mate_seqs.append(seq.revcomp())  # , f";revcomp=F, side=R, mapped={'no' if read.is_unmapped else 'yes'} dir={'fwd' if read.is_forward else 'rev'}, read={'1' if read.is_read1 else '2'}, align={read.reference_name}:{read.reference_start}-{read.reference_end} mapq={read.mapq} cigar={read.cigarstring} qname={read.query_name}"))
                    else:  # side is CLIP_LEFT
                        # the mate of a read to CLIP_LEFT originally always maps to the - strand. Thus is has to be inverted, if the bam file has done this for us.
                        # update: Since the left clipped sequence is additionally inverted to make sure the first bases is nearest to the breakpoint, revcomp has to be flipped again, making this kind of redundant.
                        if read.is_forward:
                            bp.mate_seqs.append(seq)  # , f";revcomp=F, side=L, mapped={'no' if read.is_unmapped else 'yes'} dir={'fwd' if read.is_forward else 'rev'}, read={'1' if read.is_read1 else '2'}, align={read.reference_name}:{read.reference_start}-{read.reference_end} mapq={read.mapq} cigar={read.cigarstring} qname={read.query_name}"))
                        else:
                            bp.mate_seqs.append(seq.revcomp())  # , f";revcomp=T, side=L, mapped={'no' if read.is_unmapped else 'yes'} dir={'fwd' if read.is_forward else 'rev'}, read={'1' if read.is_read1 else '2'}, align={read.reference_name}:{read.reference_start}-{read.reference_end} mapq={read.mapq} cigar={read.cigarstring} qname={read.query_name}"))


    def extend_mates(self) -> None:  # this function has side effects
        for bp in self.temporary_breakpoints:
            mates = bp.mate_seqs
            while len(mates):
                carry_over = []
                extension = []
                search_seq = bp.clipped[-min(len(bp.clipped), 10):].upper() #last 10 bases
                for r in mates:
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
                    mates = carry_over
                    extension = sorted(extension, key=len, reverse=True)[0]
                    if not len(extension):
                        break
                    extension = find_consensus(extension)
                    bp.clipped += extension

    def discovery(self, bamfile_name: str):
        self.filepath = bamfile_name
        print(f"starting discovery in file {bamfile_name}")
        self.extract_chimeric()
        print(f"finding mates of reads covering breakpoints in file {bamfile_name}")
        self.find_mates()
        print("extend clipped sequences using the mates sequences")
        self.extend_mates()
        self.polyA = [pA for pA in self.polyA if self.reference_name is not None] #mate found.
        print(f"disovery phase completed for file {bamfile_name}, found {len(self.temporary_breakpoints)} breakpoints.")
    def printOutput(self, filehandle, left: Breakpoint | PolyABreakpoint, right: Breakpoint | PolyABreakpoint):
        reference_name = None
        left_bp = None
        right_bp = None
        if isinstance(left, PolyABreakpoint):
            left_bp = f"polyA_{left.breakpoint}"
            right_bp = right.breakpoint
            reference_name = left.reference_name
        elif isinstance(right, PolyABreakpoint):
            right_bp = f"polyA_{right.breakpoint}"
            left_bp = left.breakpoint
            reference_name = right.reference_name
        else:
            left_bp = left.breakpoint
            right_bp = right.breakpoint
            reference_name = left.reference_name
        bp_name = f"{reference_name}:{left_bp}-{right_bp}"
        if isinstance(left, PolyABreakpoint):
            filehandle.write(left.clipped.revcomp().fastq(f"{bp_name}:LEFT:CLIPPED_POLYA"))
        else:
            filehandle.write(left.clipped.revcomp().fastq(f"{bp_name}:LEFT:CLIPPED"))
            filehandle.write(left.unclipped.fastq(f"{bp_name}:LEFT:ALIGNED"))
            mate_i = 0
            for mate in left.mate_seqs:
                filehandle.write(mate.revcomp().fastq(f"{bp_name}:LEFT:MATE{mate_i}"))
                mate_i += 1
        if isinstance(right, PolyABreakpoint):
            filehandle.write(right.clipped.fastq(f"{bp_name}:RIGHT:CLIPPED_POLYA"))
        else:
            filehandle.write(right.unclipped.revcomp().fastq(f"{bp_name}:RIGHT:ALIGNED"))
            filehandle.write(right.clipped.fastq(f"{bp_name}:RIGHT:CLIPPED"))
            mate_i = 0
            for mate in right.mate_seqs:
                filehandle.write(mate.revcomp().fastq(f"{bp_name}:RIGHT:MATE{mate_i}"))
                mate_i += 1

        pass
    def output(self, outfile_path: str):
        polyA = {}
        left_bp = {}
        right_bp = {}
        for pA in self.polyA:
            if pA.reference_name is None: continue
            if pA.breakpoint is None: continue
            if pA.reference_name not in polyA.keys():
                polyA[pA.reference_name] = []
            polyA[pA.reference_name].append(pA)
        for bp in self.final_left_breakpoints:
            if bp.reference_name is None: continue
            if not bp.reference_name in left_bp.keys():
                left_bp[bp.reference_name] = []
            left_bp[bp.reference_name].append(bp)
        for bp in self.final_right_breakpoints:
            if bp.reference_name is None: continue
            if not bp.reference_name in right_bp.keys():
                right_bp[bp.reference_name] = []
            right_bp[bp.reference_name].append(bp)
        for reference_name in sorted(set(polyA.keys()) | set(left_bp.keys()) | set(right_bp.keys())):
            if len(reference_name) > 5: continue #ignore weird contigs
            if not reference_name in polyA.keys():
                polyA[reference_name] = []
            if not reference_name in left_bp.keys():
                left_bp[reference_name] = []
            if not reference_name in right_bp.keys():
                right_bp[reference_name] = []
        print(f"writing result to {outfile_path}")
        with gzip.open(outfile_path, 'wt') as outfile:
            for reference_name in polyA.keys():
                print(f"processing seqname {reference_name}")
                if len(reference_name) > 5: continue #same as above
                l : List[Breakpoint] = sorted(left_bp[reference_name], key=lambda x: x.breakpoint)
                r : List[Breakpoint] = sorted(right_bp[reference_name], key=lambda x: x.breakpoint)
                p : List[polyA] = sorted(polyA[reference_name], key=lambda x: x.breakpoint)
                print(f"have {len(r)} right, {len(l)} left and {len(p)} polyA breakpoints")
                iL: int = 0
                iR: int = 0
                iP: int = 0
                nP: int = 0
                nO: int = 0
                while iL < len(l) and iR < len(r):
                    tsd: int = r[iR].breakpoint-l[iL].breakpoint
                    if tsd < 2:
                        while iP < len(p) and p[iP].breakpoint - l[iL].breakpoint < 12:
                            iP += 1
                        if iP != len(p) and p[iP].breakpoint - l[iL].breakpoint < 120 and p[iP].clip is CLIP_RIGHT:
                            self.printOutput(outfile, l[iL], p[iP])
                            nP += 1
                            iL += 1
                            continue
                        iR += 1
                        continue
                    elif tsd > 40:
                        while iP < len(p) and r[iR].breakpoint - p[iP].breakpoint > 120:
                            iP += 1
                        if iP != len(p) and r[iR].breakpoint - p[iP].breakpoint > 12 and p[iP].clip is CLIP_LEFT:
                            self.printOutput(outfile, p[iP], r[iR])
                            nP += 1
                            iR += 1
                            continue
                        iL += 1
                        continue
                    else:
                        self.printOutput(outfile, l[iL],r[iR])
                        nO += 1
                        iL += 1
                print(f"found {nO} ordinary and {nP} polyA rescued insertions in {reference_name}.")
        print("done.")