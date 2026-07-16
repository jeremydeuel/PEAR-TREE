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
from genotyping_insertion import Insertion, GT_HIGH_COVERAGE, GT_ERROR

from multiprocessing import Process, Queue
from queue import Empty
from config import CONFIG
import os

DEBUG = os.environ.get("PEARTREE_DEBUG", "") not in ("", "0", "false", "False")


TERMINATION_SIGNAL = 1


def mc_process_insertion(input_queue: Queue, output_queue: Queue, bam_file: str) -> None:
    high_cov = CONFIG['genotyping']['reads_for_high_coverage']
    with pysam.AlignmentFile(bam_file) as bam:
        while (ins := input_queue.get()) != TERMINATION_SIGNAL:
            name, left_clipped, right_clipped, left_ref, right_ref = ins
            # Never let one bad locus kill a worker: a dead worker leaves its
            # name in the ordered writer's queue with no matching result, which
            # would stall the writer and silently truncate the output tail. Emit
            # an explicit GT_ERROR result instead so every submitted name is
            # accounted for.
            try:
                i = Insertion(name)
                i.left_clipped = left_clipped
                i.right_clipped = right_clipped
                i.right_ref = right_ref
                i.left_ref = left_ref
                start = max(0, min(i.left_pos, i.right_pos))
                end = max(i.left_pos, i.right_pos) + 1
                read_count = bam.count(i.chr, start, end)
                if read_count > high_cov:
                    # coverage spikes are a rich artefact source; do not attempt a
                    # genotype here, flag as NA-like high coverage (counted by
                    # combine_genotypes' max_na).
                    gt, score_gt, score_other = GT_HIGH_COVERAGE, read_count, 0
                else:
                    i.genotype(bam)
                    gt, score_gt, score_other = i.summarise_evidence()
            except Exception as exc:  # noqa: BLE001 - must not crash the worker
                sys.stderr.write(f"genotyping failed for {name!r}: {exc!r}\n")
                gt, score_gt, score_other = GT_ERROR, 0, 0
            output_queue.put((name, gt, score_gt, score_other))


def mc_output(output_queue: Queue, output_file: str, output_list_queue: Queue):
    print("started output process")
    output_buffer = {}      # name -> (genotype, score_gt, score_other), awaiting its turn
    pending = []            # names in submission order, not yet written
    list_done = False       # seen the output-list sentinel

    def drain_names(block: bool):
        # Pull submission-order names off the list queue. Non-blocking during the
        # run (streaming fast path); blocking at the end until the sentinel, so
        # no name is lost to feeder-thread timing.
        nonlocal list_done
        while not list_done:
            try:
                n = output_list_queue.get(block)
            except Empty:
                return
            if n == TERMINATION_SIGNAL:
                list_done = True
                return
            pending.append(n)

    def flush(out):
        while pending and pending[0] in output_buffer:
            name = pending.pop(0)
            gt, score_gt, score_other = output_buffer.pop(name)
            out.write(f'{name}\t{gt}\t{int(score_gt)}\t{int(score_other)}\n')

    with gzip.open(output_file, 'wt') as out:
        out.write('insertion\tgenotype\tscore_genotype\tscore_alternative\n')
        while (o := output_queue.get()) != TERMINATION_SIGNAL:
            name, gt, score_gt, score_other = o
            output_buffer[name] = (gt, score_gt, score_other)
            drain_names(block=False)
            flush(out)
        # All workers are done; block-drain every remaining name, then flush.
        drain_names(block=True)
        flush(out)
        if output_buffer or pending:
            # Should not happen: a submitted name never got a result. Write what we
            # have out of order rather than dropping it, and fail loudly.
            sys.stderr.write(
                f"WARNING: {len(pending)} insertion(s) missing a genotype result and "
                f"{len(output_buffer)} buffered result(s) had no matching name; "
                f"output may be out of order.\n")
            for name in list(pending):
                if name in output_buffer:
                    gt, score_gt, score_other = output_buffer.pop(name)
                    out.write(f'{name}\t{gt}\t{int(score_gt)}\t{int(score_other)}\n')


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
    output_list.put(TERMINATION_SIGNAL)  # marks the end of the submission-order names
    print("submitted everything")
    for _ in range(threads):
        input_queue.put(TERMINATION_SIGNAL)  # send termination signal
    [process.join() for process in pool]
    print("terminated input")
    output_queue.put(TERMINATION_SIGNAL)
    output_process.join()
    print("terminated output")
    print("done.")
