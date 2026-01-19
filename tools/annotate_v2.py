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

import pandas as pd
import gzip
import os
import pysam
import pyliftover
from math import floor
import re
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


class Insertion:
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
        # extract inserted sequences
    def has_right_polyA(self):
        return re.search(r"[ACGT]t{6}", self.right_seq)
    def has_left_polyA(self):
        return re.search(r"a{6}[ACGT]", self.left_seq)

    def get_fasta(self) -> str:
        """
        This function returns a FASTA chunk with the inserted sequences for the insertion.
        """
        self.right_ins_seq = self.right_seq[min([self.right_seq.find(b) for b in 'acgt' if b in self.right_seq]):].upper()
        self.left_ins_seq = self.left_seq[:min([self.left_seq.find(b) for b in 'ACGT' if b in self.left_seq])].upper()
        return f'>{self.title}:R\n{self.left_ins_seq}\n>{self.title}:L\n{self.right_ins_seq}\n'

    def conclusion(self) -> str:
        """
        This function aggregates all information available to come to a conclusion
        """
        if len(self.right_dfams)==0 and len(self.left_dfams)==0 and len(self.right_maps)==0 and len(self.left_maps)==0:
            return 'artefact'
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
                            return "L1"
                        else:
                            DEBUG and print(f"      does not have left polyA -> rejected, checking left dfams")
                        # 2) L1 element on the other side
                        for m2 in self.left_dfams:
                            if "L1" in m2.model:
                                DEBUG and print(f"      found another L1 ({m2.model}) in left defams -> accepted.")
                                return m.model
                        DEBUG and print(f"      does not have a suitable left dfam, checking mappings.")
                        # 3) match near an L1 element on the other side.
                        for pos, qual, rmsks, strand in self.left_maps:
                            for r in rmsks:
                                if r.repFamily == "L1":
                                    DEBUG and print(f"      found a suitable rmsk annotation ({r}) in left mapping -> accepted.")
                                    return m.model
                    else:
                        DEBUG and print(f"    - detected RIGHT dfam model {m.model} on strand {m.strand}, checking...")
                        # dont ignore strands for Alus and since these should have the same orientation
                        # accept this in three conditions
                        # 1) the other side has a polyA, but only if the RTE is on reverse
                        if m.strand == "+" and self.has_left_polyA():
                            DEBUG and print(
                                f"      found a suitable polyA in the left mapping -> accepted.")
                            return m.model
                        # 2) same element on the other side, oriented in the same direction
                        for m2 in self.left_dfams:
                            if m2.is_active and m2.model[:3] == m.model[:3] and m2.strand == m.strand:
                                DEBUG and print(
                                    f"      found a suitable model in the left dfams ({m2.model} on strand {m2.strand}) in the left mapping -> accepted.")
                                return m.model
                        # 3) match near an L1 element on the other side.
                        for pos, qual, rmsks, strand in self.left_maps:
                            for r in rmsks:
                                if r.repName[:3] == m.model[:3]:
                                    # check same strand
                                    if (strand == "+") ^ (r.strand == "+") == (m.strand == "-"):
                                        DEBUG and print(
                                            f"      found a suitable mapping in the left dfams ({r} on strand {r.strand}, mapping is on strand {strand})-> accepted.")
                                        return m.model
            for m in self.left_dfams:
                if True or m.is_active:
                    if "L1" in m.model:
                        DEBUG and print(f"    - detected left L1 dfam model {m.model}...")
                        # ignore strand for LINEs, since these can sometimes reverse a part of the LINE at the 5' insertion site.
                        # accept this in three conditions
                        # 1) polyA on the other side
                        if self.has_right_polyA():
                            DEBUG and print(f"      has right polyA -> accepted")
                            return "L1"
                        # 2) L1 element on the other side
                        # ignore, already covered above.
                        # 3) match near an L1 element on the other side.
                        for pos, qual, rmsks, strand in self.right_maps:
                            for r in rmsks:
                                if r.repFamily == "L1":
                                    DEBUG and print(
                                        f"      found a suitable rmsk annotation ({r}) in right mapping -> accepted.")
                                    return m.model
                    else:
                        # dont ignore strands for Alus and since these should have the same orientation
                        # accept this in three conditions
                        # 1) the other side has a polyA, but only if the RTE is on reverse
                        DEBUG and print(f"    - detected left dfam model {m.model} on strand {m.strand}, checking...")
                        if m.strand == "-" and self.has_right_polyA():
                            return m.model
                        # 2) same element on the other side, oriented in the same direction
                        # ignore, already covered above.
                        # 3) match near an L1 element on the other side.
                        for pos, qual, rmsk, strand in self.right_maps:
                            for r in rmsks:
                                if r.repName[:3] == m.model[:3]:
                                    # check same strand
                                    if (strand == "+") ^ (r.strand == "+") == (m.strand == "-"):
                                        return m.model
        # check elements only mapped, but not identified by dfam
        l1_map_check = False
        for pos, qual, rmsks, strand in self.left_maps:
            for r in rmsks:
                if 'L1' in r.repName[:2] and not l1_map_check:
                    if self.has_right_polyA(): return r.repName
                    l1_map_check = True #dont check this twice
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
                    if r.repName[:3] in minus_elements:
                        return r.repName
                else:
                    if self.has_left_polyA(): return r.repName
                    if r.repName[:3] in plus_elements:
                        return r.repName
        return 'unknown'


class VariantAnnotationContainer:
    def __init__(self, sample, output):
        self.sample = sample
        self.insertions_file = CONFIG['annotate']['insertions_file'](sample)
        self.genotyping_file = CONFIG['annotate']['genotyping_file'](sample)
        self.dfam_file = CONFIG['annotate']['tmp']('dfam')(sample)
        self.sam_file = CONFIG['annotate']['tmp']('sam')(sample)
        self.fasta_file = CONFIG['annotate']['tmp']('fa.gz')(sample)
        self.output = output
        self.insertions = {}
        if not os.path.exists(self.insertions_file):
            raise FileNotFoundError(f"Insertions file {self.insertions_file} not found")
        self.read_insertions()
        if not os.path.exists(self.genotyping_file):
            raise FileNotFoundError(f"Genotyping file {self.genotyping_file} not found")
        self.read_genotyping()
        if not os.path.exists(self.fasta_file) or os.path.getsize(self.fasta_file) == 0:
            self.generate_fasta_file()
        if not os.path.exists(self.dfam_file) or os.path.getsize(self.dfam_file) == 0:
            self.generate_dfam_file()
        self.read_dfam()
        if not os.path.exists(self.sam_file) or os.path.getsize(self.sam_file) == 0:
            self.generate_sam_file()
        self.read_sam()

    def read_insertions(self):
        """
        This fucntion reads the insertions.combined.txt.gz file and generates new Insertions objects including sequences and populates the self.insertions dictionary with these.
        """
        print(f"reading insertions file {self.insertions_file}...")
        with gzip.open(self.insertions_file, 'rt') as ifh:
            for line in ifh:
                line = line.strip()
                if not line: continue
                assert line[0] == '@'
                assert line[-1] == 'L'
                left_title = line[1:-2]
                left_seq = ifh.readline().strip()
                assert ifh.readline().strip() == '+'
                left_qual = ifh.readline().strip()
                right_title = ifh.readline().strip()
                assert right_title[0] == '@'
                assert right_title[-1] == 'R'
                right_title = right_title[1:-2]
                assert right_title == left_title
                right_seq = ifh.readline().strip()
                assert ifh.readline().strip() == '+'
                right_qual = ifh.readline().strip()
                self.insertions[left_title] = Insertion(left_title, left_seq, right_seq)
        print(f"done reading insertions file {self.insertions_file}, read {len(self.insertions)} insertions.")

    def read_genotyping(self):
        """
        This function read sthe genotyping file, calculates the number of tips with insertions, artefacts and wild-type calls, updates the Insertion object with these numbers and filteres insertions with at least one heterozygous or homozygous call.
        """
        print(f"reading genotyping file {self.genotyping_file}...")
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
                    self.insertions[title].nins = sum([1 for gt in line if gt in ('heterozygous','homozygous')])
                    self.insertions[title].nart = sum([1 for gt in line if gt == "artefact"])
                    self.insertions[title].nwt = sum([1 for gt in line if gt == "wild-type"])
        self.insertions = {key: value for key, value in self.insertions.items() if value.nins>0}
        print(f"imported genotypes, found {len(self.insertions)} insertions with one or more tips containing insertions.")

    def generate_fasta_file(self):
        """
        This function generates the fasta file necessary for dfam and samtools
        """
        with gzip.open(self.fasta_file, 'wt') as fh:
            for insertion in self.insertions.values():
                fh.write(insertion.get_fasta())

    def generate_dfam_file(self):
        """
        This function generates the dfam file using dfamscan.pl
        """
        print(f"running DFAM on {self.sample}")
        assert os.path.exists(self.fasta_file)
        assert os.path.exists(CONFIG['annotate']['dfamscan'])
        assert os.path.exists(CONFIG['annotate']['hmm'])
        assert os.path.exists(CONFIG['annotate']['dfamscan'])
        assert os.access(CONFIG['annotate']['dfamscan'], os.X_OK)
        # set PATH variable to include hmmer scripts.
        os.environ["PATH"] = CONFIG['annotate']['hmmer'] + ":" + os.environ["PATH"]
        assert 0 == os.system("nhmmscan -h > /dev/null 2>&1") # check that nhmmscan exists in path, otherwise dfamscan.pl will not work.
        assert 0 == os.system(
            f"{CONFIG['annotate']['dfamscan']} --fastafile {self.fasta_file} --hmmfile {CONFIG['annotate']['hmm']} --cpu {os.cpu_count()} --dfam_outfile {self.dfam_file}")
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
                    self.insertions[insertion].right_dfams.append(Dfam_Annotation(line))
                    n_annot_right += 1
                elif line[2][-1] == 'R':
                    self.insertions[insertion].left_dfams.append(Dfam_Annotation(line))
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
            if delta_start > 500: continue
            delta_end = rmsk_library[seqname][lpos].end - start
            if delta_end < -500: continue
            if delta_end > -500 and delta_start < 500:
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

        rightn = 0
        leftn = 0
        with pysam.AlignmentFile(self.sam_file) as sam:
            for read in sam:
                if read.is_qcfail: continue
                #if read.is_secondary: continue
                #if read.is_supplementary: continue
                if read.is_unmapped: continue
                local_rmsks = None
                insertion = read.query_name[:-2]
                if not insertion in self.insertions.keys(): continue
                if (read.query_name[-1] == "R") ^ read.is_forward:
                    local_rmsks = self.get_rmsk(rmsk_library, (read.reference_name, read.reference_start))
                    if lo is not None:
                        co = lo.convert_coordinate(read.reference_name, read.reference_start,
                                                   '+' if read.is_forward else '-')
                    else:
                        co = [(read.reference_name, read.reference_start, '+' if read.is_forward else '-')]
                else:
                    local_rmsks = self.get_rmsk(rmsk_library, (read.reference_name, read.reference_end))
                    if lo is not None:
                        co = lo.convert_coordinate(read.reference_name, read.reference_end,
                                                   '+' if read.is_forward else '-')
                    else:
                        co = [(read.reference_name, read.reference_end, '+' if read.is_forward else '-')]
                if co:
                    if read.query_name[-1] == "L":
                        rightn += 1
                        self.insertions[insertion].right_maps.append((f"{co[0][0]}:{co[0][1]}{co[0][2]}", read.mapping_quality, local_rmsks, co[0][2]))
                    elif read.query_name[-1] == "R":
                        leftn += 1
                        self.insertions[insertion].left_maps.append((f"{co[0][0]}:{co[0][1]}{co[0][2]}", read.mapping_quality, local_rmsks,
                                                                        co[0][2]))
                    else:
                        raise ValueError(f"Unknown insertion side {read.query_name}, expected R or L.")
        print(f"imported {rightn} right mappings and {leftn} left mappings.")

    def print(self):
        for key, insertion in self.insertions.items():
            print(f"> insertion {key} found in {insertion.nins} tip{'s' if insertion.nins!=1 else ''} ({insertion.nwt}=wt, {insertion.nart}=art)")
            print(f" {insertion.conclusion()}")
            print(f"  RIGHT INSERTION: {insertion.right_seq}")
            for dfam in insertion.right_dfams:
                print(f"    {str(dfam)}")
            for pos, mq, rmsks, strand in insertion.right_maps:
                print(f"    {pos} {[str(r) for r in rmsks]}")
            if re.match(r"[ACGT]t{6,}", insertion.right_seq):
                print(f"    is polyA")
            print(f"  LEFT INSERTION: {insertion.left_seq}")
            for dfam in insertion.left_dfams:
                print(f"    {str(dfam)}")
            for pos, mq, rmsks, strand in insertion.left_maps:
                print(f"    {pos} {[str(r) for r in rmsks]}")
            if re.match(r"a{6,}[ACGT]",insertion.left_seq):
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
    f.print()




