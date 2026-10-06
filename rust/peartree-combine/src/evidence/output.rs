//! `<stem>.insertions.evidence.tsv.gz` + `<stem>.insertions.reads.fa.gz`. OWNER: P5.
//!
//! Mirrors `write_evidence_outputs` (src/combine_insertions_evidence.py:1307) and
//! `_evidence_paths` (combine_insertions.py:103). SPEC.md §6.4, §6.5.
//!
//! Streaming: one TSV row (`JunctionRecord::tsv`) and one reads text at a time; detached
//! records' reads come back from the store by `ReadsRef` (contiguous per chunk, so the
//! survivor section reads the scratch file almost sequentially).

use crate::context::Ctx;
use crate::evidence::EvidenceState;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

/// python `EVIDENCE_TSV_COLUMNS` (evidence.py:49).
pub const EVIDENCE_TSV_COLUMNS: [&str; 26] = [
    "insertion_id",
    "side",
    "n_reads",
    "n_fragments",
    "n_independent",
    "n_samples",
    "n_mates",
    "supported",
    "clip_consensus",
    "consensus_depth",
    "polya_len_median",
    "polya_len_range",
    "beyond_polya",
    "beyond_polya_support",
    "polya_end",
    "n_duplicates",
    "n_cross_sample_identical",
    "fail_reason",
    "consensus_stop",
    "n_dup_coord",
    "n_dup_seq",
    "member_loci",
    "n_short_used",
    "n_short_rejected",
    "n_short_mate_inside",
    "n_independent_no_short",
];

/// `_evidence_paths(combined)`: strip `.combined.txt.gz` (else `.txt.gz`) ->
/// (`{stem}.insertions.evidence.tsv.gz`, `{stem}.insertions.reads.fa.gz`).
pub fn evidence_paths(combined: &str) -> (PathBuf, PathBuf) {
    let stem = combined
        .strip_suffix(".combined.txt.gz")
        .or_else(|| combined.strip_suffix(".txt.gz"))
        .unwrap_or(combined);
    (
        PathBuf::from(format!("{stem}.insertions.evidence.tsv.gz")),
        PathBuf::from(format!("{stem}.insertions.reads.fa.gz")),
    )
}

fn gz_create(p: &Path) -> Result<GzEncoder<BufWriter<File>>, String> {
    let f = File::create(p).map_err(|e| format!("cannot create {}: {e}", p.display()))?;
    Ok(GzEncoder::new(BufWriter::with_capacity(1 << 20, f), Compression::new(6)))
}

/// Write the header + one TSV row per record of each name in `names` (in order; a name without
/// records writes nothing; a repeated name is written again), and the matching reads. `names` =
/// surviving insertion names in combined.txt.gz order, then `sorted(failed)` (byte order).
/// Attached records (rows still held, e.g. absorb re-evaluations) are rendered with the output
/// name like python; detached ones were rendered with `insertion_id` at detach time.
pub fn write_evidence_outputs(state: &EvidenceState, names: &[String], tsv: &Path, fa: &Path, ctx: &Ctx) -> Result<(), String> {
    let mut t = gz_create(tsv)?;
    let mut f = gz_create(fa)?;
    let we = |p: &Path, e: std::io::Error| format!("write {}: {e}", p.display());
    t.write_all(EVIDENCE_TSV_COLUMNS.join("\t").as_bytes()).map_err(|e| we(tsv, e))?;
    t.write_all(b"\n").map_err(|e| we(tsv, e))?;
    let (mut n_rows, mut n_reads) = (0u64, 0u64);
    for name in names {
        let Some(&uid) = state.records.get(name) else { continue };
        let Some(recs) = state.recmap.get(uid as usize).and_then(|r| r.as_ref()) else { continue };
        for rec in recs {
            t.write_all(rec.tsv().as_bytes()).map_err(|e| we(tsv, e))?;
            n_rows += 1;
            let txt = match (&rec.rows, rec.reads_ref) {
                (Some(_), _) => rec.render_reads(name, &ctx.files),
                (None, Some(r)) => state.store.read_reads(r),
                (None, None) => return Err(format!("evidence record {name} {} has neither rows nor reads", rec.side.as_str())),
            };
            f.write_all(txt.as_bytes()).map_err(|e| we(fa, e))?;
            n_reads += (txt.bytes().filter(|&b| b == b'\n').count() / 2) as u64;
        }
    }
    t.finish().and_then(|mut w| w.flush()).map_err(|e| we(tsv, e))?;
    f.finish().and_then(|mut w| w.flush()).map_err(|e| we(fa, e))?;
    println!("wrote {n_rows} junction rows -> {}; {n_reads} reads -> {}", tsv.display(), fa.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paths() {
        let (t, f) = evidence_paths("/x/P1.combined.txt.gz");
        assert_eq!(t, PathBuf::from("/x/P1.insertions.evidence.tsv.gz"));
        assert_eq!(f, PathBuf::from("/x/P1.insertions.reads.fa.gz"));
        assert_eq!(evidence_paths("P1.txt.gz").0, PathBuf::from("P1.insertions.evidence.tsv.gz"));
        assert_eq!(evidence_paths("P1").0, PathBuf::from("P1.insertions.evidence.tsv.gz"));
    }

    #[test]
    fn header_matches_python() {
        // python "\t".join(EVIDENCE_TSV_COLUMNS)
        assert_eq!(
            EVIDENCE_TSV_COLUMNS.join("\t"),
            "insertion_id\tside\tn_reads\tn_fragments\tn_independent\tn_samples\tn_mates\tsupported\t\
             clip_consensus\tconsensus_depth\tpolya_len_median\tpolya_len_range\tbeyond_polya\t\
             beyond_polya_support\tpolya_end\tn_duplicates\tn_cross_sample_identical\tfail_reason\t\
             consensus_stop\tn_dup_coord\tn_dup_seq\tmember_loci\tn_short_used\tn_short_rejected\t\
             n_short_mate_inside\tn_independent_no_short"
        );
    }
}
