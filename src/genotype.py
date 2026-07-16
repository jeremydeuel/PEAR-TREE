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


def score_locus(bam: pysam.AlignmentFile, payload, high_cov: int):
    """Genotype one locus against an already-open BAM.

    ``payload`` is the ``(name, left_clipped, right_clipped, left_ref, right_ref)``
    tuple carried on the work queue. Returns the 7-tuple
    ``(genotype, score_gt, score_other, coverage, n_alt, n_ref, n_art)``.

    This is the single per-locus code path shared by the single-sample driver
    (:func:`genotype`) and the many-sample batch driver (:func:`genotype_batch`),
    so both emit byte-identical rows for the same (locus, BAM).

    Never let one bad locus kill a worker: a dead worker leaves its name in the
    ordered writer's queue with no matching result, which would stall the writer
    and silently truncate the output tail. Emit an explicit GT_ERROR result
    instead so every submitted name is accounted for. The total spanning depth +
    per-read vote counts are emitted alongside the call so downstream (e.g. a
    phylogeny-aware single-read vs barcode-switch likelihood) has the raw
    evidence. Defaults cover the high-coverage / error paths, which do not run
    summarise_evidence.
    """
    name, left_clipped, right_clipped, left_ref, right_ref = payload
    coverage = n_alt = n_ref = n_art = 0
    try:
        i = Insertion(name)
        i.left_clipped = left_clipped
        i.right_clipped = right_clipped
        i.right_ref = right_ref
        i.left_ref = left_ref
        start = max(0, min(i.left_pos, i.right_pos))
        end = max(i.left_pos, i.right_pos) + 1
        read_count = bam.count(i.chr, start, end)
        coverage = read_count
        if read_count > high_cov:
            # coverage spikes are a rich artefact source; do not attempt a
            # genotype here, flag as NA-like high coverage (counted by
            # combine_genotypes' max_na).
            gt, score_gt, score_other = GT_HIGH_COVERAGE, read_count, 0
        else:
            i.genotype(bam)
            gt, score_gt, score_other = i.summarise_evidence()
            n_alt, n_ref, n_art = i.n_alt, i.n_ref, i.n_art
    except Exception as exc:  # noqa: BLE001 - must not crash the worker
        sys.stderr.write(f"genotyping failed for {name!r}: {exc!r}\n")
        gt, score_gt, score_other = GT_ERROR, 0, 0
    return gt, score_gt, score_other, coverage, n_alt, n_ref, n_art


# on-disk header + row format, shared by the single-sample and batch writers so
# every genotype output file is byte-identical regardless of which driver wrote it.
OUTPUT_HEADER = ('insertion\tgenotype\tscore_genotype\tscore_alternative\t'
                 'coverage\tn_alt\tn_ref\tn_art\n')


def format_row(name, gt, score_gt, score_other, coverage, n_alt, n_ref, n_art):
    return (f'{name}\t{gt}\t{int(score_gt)}\t{int(score_other)}\t'
            f'{int(coverage)}\t{int(n_alt)}\t{int(n_ref)}\t{int(n_art)}\n')


def mc_process_insertion(input_queue: Queue, output_queue: Queue, bam_file: str) -> None:
    high_cov = CONFIG['genotyping']['reads_for_high_coverage']
    with pysam.AlignmentFile(bam_file) as bam:
        while (ins := input_queue.get()) != TERMINATION_SIGNAL:
            name = ins[0]
            result = score_locus(bam, ins, high_cov)
            output_queue.put((name, *result))


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
            out.write(format_row(name, *output_buffer.pop(name)))

    with gzip.open(output_file, 'wt') as out:
        # score_genotype/score_alternative are quality-margin sums; coverage is the total
        # spanning read depth and n_alt/n_ref/n_art are the per-read vote counts (raw
        # evidence for downstream phylogeny / barcode-switch models).
        out.write(OUTPUT_HEADER)
        while (o := output_queue.get()) != TERMINATION_SIGNAL:
            name, gt, score_gt, score_other, coverage, n_alt, n_ref, n_art = o
            output_buffer[name] = (gt, score_gt, score_other, coverage, n_alt, n_ref, n_art)
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
                    out.write(format_row(name, *output_buffer.pop(name)))


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


# --------------------------------------------------------------------------- #
# Batch genotyping: one insertion contract, many sample BAMs.
#
# A phylogenetic tree of hundreds of colonies genotypes the SAME insertion
# contract against every colony's BAM. Running the single-sample `genotype()`
# once per colony re-parses the contract, re-imports the interpreter/pysam and
# re-spawns a worker pool for each of the hundreds of samples, and each of those
# runs pays its own tail-latency underutilisation. `genotype_batch()` collapses
# all of that into one invocation: the contract is parsed ONCE, a single
# persistent pool drains one (sample, locus) work queue, and every worker stays
# busy across sample boundaries. Each sample's output file is byte-identical to
# what `genotype()` would have written for that (BAM, contract) — the per-locus
# scoring goes through the same `score_locus()` and the same `format_row()`.
# --------------------------------------------------------------------------- #

# how many BAM handles a worker keeps open at once. Tasks are fed sample-major,
# so at any moment a worker touches at most a couple of adjacent samples; a small
# cache bounds open file descriptors (threads * MAX_OPEN_BAMS) well under ulimit
# while never reopening a BAM within a sample.
MAX_OPEN_BAMS = 4


def mc_batch_worker(task_queue: Queue, result_queue: Queue,
                    payloads, bam_files, high_cov: int) -> None:
    """Drain (sample_idx, locus_idx) tasks, scoring each locus on its sample BAM.

    ``payloads`` (the whole parsed contract) and ``bam_files`` are broadcast to
    every worker once at spawn, so a task on the queue is just two small ints.
    Open BAM handles are cached (LRU by open order) to avoid reopening within a
    sample while bounding descriptors.
    """
    open_bams = {}   # sample_idx -> pysam.AlignmentFile
    open_order = []  # sample_idx in open order (oldest first)
    try:
        while (task := task_queue.get()) != TERMINATION_SIGNAL:
            s_idx, l_idx = task
            bam = open_bams.get(s_idx)
            if bam is None:
                bam = pysam.AlignmentFile(bam_files[s_idx])
                open_bams[s_idx] = bam
                open_order.append(s_idx)
                if len(open_order) > MAX_OPEN_BAMS:
                    open_bams.pop(open_order.pop(0)).close()
            result = score_locus(bam, payloads[l_idx], high_cov)
            result_queue.put((s_idx, l_idx, result))
    finally:
        for bam in open_bams.values():
            bam.close()


def _write_sample(output_file: str, names, results) -> None:
    """Write one sample's gzip output: header + one row per locus, in contract order."""
    with gzip.open(output_file, 'wt') as out:
        out.write(OUTPUT_HEADER)
        for name, result in zip(names, results):
            out.write(format_row(name, *result))


def mc_batch_writer(result_queue: Queue, output_files, names, n_loci: int) -> None:
    """Collect per-sample results and flush each sample's file once it is complete.

    Results arrive interleaved across samples. Each sample gets a pre-sized slot
    list indexed by locus_idx (so contract order is preserved with no sorting) and
    a remaining counter; when a sample's last locus lands, its file is written and
    the buffer freed. With sample-major feeding only a few samples are ever
    in-flight, so memory stays bounded.
    """
    print("started batch output process")
    buffers = {}  # sample_idx -> [slots, remaining]
    while (item := result_queue.get()) != TERMINATION_SIGNAL:
        s_idx, l_idx, result = item
        buf = buffers.get(s_idx)
        if buf is None:
            buf = [[None] * n_loci, n_loci]
            buffers[s_idx] = buf
        buf[0][l_idx] = result
        buf[1] -= 1
        if buf[1] == 0:
            _write_sample(output_files[s_idx], names, buf[0])
            del buffers[s_idx]
            print(f"wrote {output_files[s_idx]}")
    if buffers:
        # Should not happen: a sample never received all its loci. Write what we
        # have (missing rows would raise in format_row) and fail loudly.
        sys.stderr.write(
            f"WARNING: {len(buffers)} sample(s) did not receive every locus result; "
            f"their output files were not written.\n")


def genotype_batch(insertion_file: str, bam_files, output_files, threads: int = 1) -> None:
    """Genotype one insertion contract against many sample BAMs in a single run.

    ``bam_files[k]`` is genotyped into ``output_files[k]``. The contract is parsed
    once; a persistent pool of ``threads`` workers drains a (sample, locus) queue
    fed sample-major; a single writer flushes each sample's file as it completes.
    """
    if len(bam_files) != len(output_files):
        raise ValueError("bam_files and output_files must have the same length")
    insertions = list(Insertion.import_file(insertion_file))
    payloads = [(i.name, i.left_clipped, i.right_clipped, i.left_ref, i.right_ref)
                for i in insertions]
    names = [i.name for i in insertions]
    n_loci = len(payloads)
    high_cov = CONFIG['genotyping']['reads_for_high_coverage']
    print(f"batch genotyping {len(bam_files)} sample(s) x {n_loci} loci, threads={threads}")

    task_queue = Queue()
    result_queue = Queue()
    workers = [Process(target=mc_batch_worker,
                       args=(task_queue, result_queue, payloads, bam_files, high_cov))
               for _ in range(threads)]
    writer = Process(target=mc_batch_writer,
                     args=(result_queue, output_files, names, n_loci))
    [w.start() for w in workers]
    writer.start()

    # sample-major task order keeps each worker's open-BAM set tiny.
    for s_idx in range(len(bam_files)):
        for l_idx in range(n_loci):
            task_queue.put((s_idx, l_idx))
    for _ in range(threads):
        task_queue.put(TERMINATION_SIGNAL)
    [w.join() for w in workers]
    result_queue.put(TERMINATION_SIGNAL)
    writer.join()
    print("done.")
