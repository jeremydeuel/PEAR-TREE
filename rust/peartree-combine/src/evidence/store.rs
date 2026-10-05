//! Sidecar row routing + detached reads storage (python `_ShardStore`). OWNER: P5.
//!
//! Mirrors src/combine_insertions_evidence.py:838-980 in PURPOSE, not in file format: the
//! on-disk layout is private (the python shard files are not an output). Contract:
//!   * one streaming pass per sidecar (files in parallel), keeping only rows whose
//!     (file, parsed locus) is a member locus of some insertion, routed to every chunk that
//!     needs them; within a (member, side) the row order is the sidecar line order;
//!   * `chunk_rows(c)` gives a chunk's rows keyed by (member, side);
//!   * `lookup(member, side)` serves the parent-side re-evaluation in absorb_one_sided;
//!   * evaluated chunks move each record's rendered reads (FASTA text) to disk
//!     (`put_reads` -> `ReadsRef`), so peak memory is bounded by the chunks in flight;
//!   * the scratch directory `<stem>.evidence_shards` is created fresh (rm -rf first) and
//!     removed by `cleanup()` after the evidence outputs are written.
//! Sidecar path rule (python `sidecar_path`): `<txt.gz>.evidence.tsv.gz` (Rust discovery), else
//! `<stem>.evidence.tsv.gz` when only that one exists. SPEC.md §4.0, §7.

use crate::context::Ctx;
use crate::evidence::junction::ReadsRef;
use crate::evidence::row::EvidenceRow;
use crate::model::{FileId, Insertion, Member, Side};
use std::path::{Path, PathBuf};

/// python `sidecar_path(txt_gz)`.
pub fn sidecar_path(txt_gz: &str) -> PathBuf {
    let rust = format!("{txt_gz}.evidence.tsv.gz");
    let stem = txt_gz.strip_suffix(".txt.gz").unwrap_or(txt_gz);
    let alt = format!("{stem}.evidence.tsv.gz");
    if !Path::new(&rust).exists() && Path::new(&alt).exists() {
        PathBuf::from(alt)
    } else {
        PathBuf::from(rust)
    }
}

/// Rows of one chunk, keyed by (member, side), each list in sidecar order.
pub struct ChunkRows {
    // P5
}

impl ChunkRows {
    /// python `lookup(m, side)` -> rows (empty when none).
    pub fn get(&self, m: &Member, side: Side) -> &[EvidenceRow] {
        todo!("P5")
    }
}

pub struct Store {
    // P5
}

impl Store {
    /// Route every wanted sidecar row (python `_ShardStore.build`). `members[k]` = member loci of
    /// `insertions[k]`. Returns the store and the chunks as half-open insertion index ranges
    /// (any chunking is allowed: outputs never depend on it).
    pub fn build(dir: &Path, accepted: &[FileId], insertions: &[Insertion], members: &[Vec<Member>], ctx: &Ctx) -> Result<(Store, Vec<(usize, usize)>), String> {
        todo!("P5")
    }

    /// basenames with an existing sidecar (python `have`), in input order
    pub fn have(&self) -> &[FileId] {
        todo!("P5")
    }

    /// number of routed rows (log only)
    pub fn n_rows(&self) -> u64 {
        todo!("P5")
    }

    pub fn chunk_rows(&self, chunk: usize) -> Result<ChunkRows, String> {
        todo!("P5")
    }

    /// parent-side lookup (absorb_one_sided re-evaluation)
    pub fn lookup(&self, m: &Member, side: Side) -> Vec<EvidenceRow> {
        todo!("P5")
    }

    /// store a record's rendered reads text; thread-safe
    pub fn put_reads(&self, chunk: usize, text: &str) -> ReadsRef {
        todo!("P5")
    }

    pub fn read_reads(&self, r: ReadsRef) -> String {
        todo!("P5")
    }

    pub fn cleanup(&mut self) {
        todo!("P5")
    }
}
