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

import pysam
import gzip
import sys
from genotyping_insertion import Insertion

from multiprocessing import Process, Queue, Pipe
from queue import Empty
from config import CONFIG

DEBUG = True


TERMINATION_SIGNAL = 1


def mc_count_insertion(insertion, bam_file,  connection):
    start = min(insertion.left_pos, insertion.right_pos)
    end = max(insertion.left_pos, insertion.right_pos, start + 1)
    with pysam.AlignmentFile(bam_file) as bam:
        read_count = bam.count(insertion.chr, start, end)
    connection.send(read_count)

def mc_process_insertion(input_queue: Queue, output_queue: Queue, bam_file: str) -> None:
    count_out, count_in = Pipe(duplex=False)
    with pysam.AlignmentFile(bam_file) as bam:
        while ins := input_queue.get():
            if ins is TERMINATION_SIGNAL:
                return
            name, left_clipped, right_clipped, left_ref, right_ref = ins
            i = Insertion(name)
            i.left_clipped = left_clipped
            i.right_clipped = right_clipped
            i.right_ref = right_ref
            i.left_ref = left_ref
            read_count = bam.count(i.chr, min(i.left_pos, i.right_pos), max(i.left_pos, i.right_pos)+1)
            if read_count > CONFIG['genotyping']['reads_for_high_coverage']:
                i.evidence['hc'] = read_count
            else:
                i.genotype(bam)
            i.summarise_evidence()
            output_queue.put((
                i.name,
                i.genotype_string,
                i.evidence
            ))

def mc_output(output_queue: Queue, output_file: str, output_list_queue: Queue):
    print(f"started output process")
    output_buffer = {}
    output_list = []
    with gzip.open(output_file, 'wt') as out:
        out.write('insertion\tgenotype\tright_ins\tleft_ins\twt\tart_left\tart_right\tart_both\tins_both\tno_cov\thigh_cov\n')
        while o := output_queue.get():
            if o is TERMINATION_SIGNAL:
                return
            try:
                while ol := output_list_queue.get(False):
                    output_list.append(ol)
            except Empty:
                pass
            name, genotype_string, evidence = o
            #print(f"output got {name}")
            output_buffer[name] = (genotype_string, evidence)
            while len(output_list):
                if output_list[0] not in output_buffer.keys():
                    break
                current_name = output_list.pop(0)
                genotype_string, evidence = output_buffer[current_name]
                out.write(f'{current_name}\t{genotype_string}\t{evidence["ir"]}\t{evidence["il"]}\t{evidence["wt"]}\t{evidence["ar"]}\t{evidence["al"]}\t{evidence["ad"]}\t{evidence["id"]}\t{evidence["nc"]}\t{evidence["hc"]}\n')
                del output_buffer[current_name]
                #print(f"outputted {current_name}, current length = {len(output_buffer)}")

def genotype(insertion_file: str, bam_file: str, output_file: str, threads=1):
    input_queue = Queue()
    output_queue = Queue()
    output_list = Queue()
    pool = [Process(target=mc_process_insertion, args=(input_queue, output_queue, bam_file)) for _ in range(threads)]
    output_process = Process(target=mc_output, args=(output_queue, output_file, output_list))
    [process.start() for process in pool]
    output_process.start()
    for i in Insertion.import_file(insertion_file):
        input_queue.put((
            i.name,
            i.left_clipped,
            i.right_clipped,
            i.left_ref,
            i.right_ref
        ))
        output_list.put(i.name)
    print(f"submitted everything")
    for i in range(threads):
        input_queue.put(TERMINATION_SIGNAL)  # send termination signal
    [process.join() for process in pool]
    print(f"terminated input")
    output_queue.put(TERMINATION_SIGNAL)
    output_process.join()
    print(f"terminated output")
    print("done.")
