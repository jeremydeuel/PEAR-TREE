"""Shared insertion-type simulation library (TPRT-hallmark overhaul).

Used by `test/val1/simulate.py` (genome-free, alignment-by-construction BAMs) and
`test/fullstack/build_donor.py` (hs1 donor haplotypes -> reads -> bwa).

Modules
-------
seqs      sequence utilities, genome access (.2bit / FASTA), literature-calibrated
          parameter distributions (TSD, poly-A, truncation, twin priming, ...)
library   real element sequences + real 3'/5' source flanks: `resources/rte_library/`
          (--rte-library) or extracted from hs1 + RepeatMasker as a fallback
models    the insertion-type catalogue: builds an `Event` (inserted sequence + target-site
          geometry + SPEC truth labels) and applies it to a reference window
reads     fragment / read-pair sampler with Illumina poly-A (homopolymer) jitter,
          unflagged PCR duplicates, per-read errors, and alignment-by-construction
truth     the truth-TSV schema (see test/simlib/TRUTH_SCHEMA.md)
"""
