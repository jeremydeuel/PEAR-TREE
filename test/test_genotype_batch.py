# PEAR-TREE - batch genotyping integration test
#
# Locks in the core contract of genotype_batch(): genotyping one insertion
# contract against many sample BAMs in a single run must produce, for each
# sample, output byte-identical to what the single-sample genotype() driver
# writes for that (BAM, contract). This is what lets a phylogenetic-tree run of
# hundreds of colonies use the batch path without changing any downstream call.
#
# Requires pysam and the test_data/ fixtures; skipped automatically otherwise.
#
# Run with:  python test/test_genotype_batch.py     (from the repo root)
#        or:  pytest test/test_genotype_batch.py

import gzip
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

REPO = os.path.join(os.path.dirname(__file__), "..")
BAM = os.path.join(REPO, "test_data", "test.bam")
CONTRACT = os.path.join(REPO, "test_data", "test_step2.txt.genotyping.txt.gz")


def _payload(path):
    """Decompressed text of a gzip output file (ignores the gzip header, whose
    mtime/name bytes legitimately differ between two writes of identical data)."""
    with gzip.open(path, "rt") as f:
        return f.read()


def run():
    try:
        import pysam  # noqa: F401
    except ImportError:
        print("SKIP: pysam not available")
        return True
    if not (os.path.exists(BAM) and os.path.exists(CONTRACT)):
        print("SKIP: test_data fixtures not found")
        return True

    from genotype import genotype, genotype_batch

    with tempfile.TemporaryDirectory() as tmp:
        single = os.path.join(tmp, "single.txt.gz")
        genotype(CONTRACT, BAM, single, threads=1)

        # two samples, both the same BAM: each must reproduce the single-sample output
        outs = [os.path.join(tmp, "batch_A.txt.gz"), os.path.join(tmp, "batch_B.txt.gz")]
        genotype_batch(CONTRACT, [BAM, BAM], outs, threads=2)

        ref = _payload(single)
        assert ref.startswith("insertion\tgenotype\t"), "single-sample output malformed"
        for out in outs:
            got = _payload(out)
            assert got == ref, (
                f"batch output {os.path.basename(out)} differs from single-sample:\n"
                f"--- single ---\n{ref}\n--- batch ---\n{got}"
            )
    print("OK: batch output byte-identical to single-sample (2 samples)")
    return True


# pytest entry point
def test_batch_matches_single_sample():
    assert run()


if __name__ == "__main__":
    ok = run()
    sys.exit(0 if ok else 1)
