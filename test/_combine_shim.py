"""Import shim for the combine_insertions tests.

combine_insertions does `from config import CONFIG` (src/config.py is a gitignored copy of a
cluster/config.py.*) and combine_insertions_get_sequence opens a 2bit genome at import time.
Neither exists off-farm, so we register an in-memory `config` and a fake
`combine_insertions_get_sequence` (backed by a synthetic genome) just long enough to import
the modules under test, then restore sys.modules so other test files are unaffected. The
imported modules keep their reference to TEST_CONFIG, so tests tweak it in place.
"""
import os
import random
import sys
import types

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

_rng = random.Random(20261004)
GENOME = {"chr1": "".join(_rng.choice("ACGT") for _ in range(40000))}

TEST_CONFIG = {
    "version": "test",
    "discovery": {"min_clip_len": 12, "min_good_bases": 10, "min_adapterlen_for_clip": 4},
    "adapters": ["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"],
    "genotyping": {"max_bases": 12},
    "combine_insertions": {
        "genome_2bit": "/nonexistent.2bit",
        "exclude_files_with_many_insertions": 1_000_000,
        "samtools_executable": "samtools",
        "bowtie2_executable": "bowtie2",
        "bowtie2_index": "/nonexistent",
        "bowtie2_index2": "/nonexistent",
        "bowtie2_index2_lo": "/nonexistent.chain.gz",
        "clean_remap_max_insertion": 12,
        "clean_remap_min_as": -15,
        # TPRT keys deliberately absent: the code must default them to legacy behaviour.
    },
}


def _get_sequence(seqname, start, end):
    s = GENOME.get(seqname, "")
    if start >= end:
        return ""
    return s[max(0, start):end]


def _import_under_shim():
    saved = {k: sys.modules.get(k) for k in ("config", "combine_insertions_get_sequence")}
    cfg = types.ModuleType("config")
    cfg.CONFIG = TEST_CONFIG
    gs = types.ModuleType("combine_insertions_get_sequence")
    gs.get_sequence = _get_sequence
    sys.modules["config"] = cfg
    sys.modules["combine_insertions_get_sequence"] = gs
    try:
        import consensus  # noqa: F401  (legacy column-wise consensus, for comparison)
        import combine_insertions  # noqa: F401
        return consensus, combine_insertions
    finally:
        for k, v in saved.items():
            if v is None:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = v


consensus, combine_insertions = _import_under_shim()
