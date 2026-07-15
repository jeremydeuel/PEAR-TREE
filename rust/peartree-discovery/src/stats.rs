//! Reject-counter sidecar (OBS-1). Mirrors `Breakpoint.stats` in src/breakpoint.py
//! (the per-side dict printed at discovery.py:202) with the identical field set, so
//! the JSON sidecar can be compared directly against the Python counters.
//!
//! Every counter is incremented at exactly the return path in `model::join` that the
//! matching `Breakpoint.stats[...][...] += 1` occupies in the Python join.
//!
//! NOTE (SPD-3): single-threaded plain fields for now. When contig-level rayon lands,
//! switch to per-thread accumulation or atomics and sum at the merge.

use crate::config::CLIP_LEFT;

#[derive(Default, Clone)]
pub struct SideStats {
    pub too_few: u64,
    pub rescued_pa: u64,
    pub excluded: u64,
    pub too_few_after_filter: u64,
    pub clipped_failed: u64,
    pub unclipped_failed: u64,
    pub polymer: u64,
    pub passed: u64,
}

#[derive(Default, Clone)]
pub struct Stats {
    pub left: SideStats,
    pub right: SideStats,
    // Feature A discordant-anchor counters (global, not per-side). All 0 unless
    // `discordant_anchor` is on, so the default sidecar gains only zero-valued keys.
    pub disc_obs: u64,
    pub disc_clusters: u64,
    pub disc_paired: u64,
    pub disc_rejected_rte: u64,
}

impl Stats {
    /// Select the per-side counters, keyed like Python's `stats[side]`.
    pub fn side_mut(&mut self, side: i32) -> &mut SideStats {
        if side == CLIP_LEFT {
            &mut self.left
        } else {
            &mut self.right
        }
    }

    /// Sum another Stats into this one (SPD-3: merge per-contig worker counters).
    pub fn merge(&mut self, other: &Stats) {
        for (dst, src) in [(&mut self.left, &other.left), (&mut self.right, &other.right)] {
            dst.too_few += src.too_few;
            dst.rescued_pa += src.rescued_pa;
            dst.excluded += src.excluded;
            dst.too_few_after_filter += src.too_few_after_filter;
            dst.clipped_failed += src.clipped_failed;
            dst.unclipped_failed += src.unclipped_failed;
            dst.polymer += src.polymer;
            dst.passed += src.passed;
        }
        self.disc_obs += other.disc_obs;
        self.disc_clusters += other.disc_clusters;
        self.disc_paired += other.disc_paired;
        self.disc_rejected_rte += other.disc_rejected_rte;
    }

    /// Serialise to JSON. The `left`/`right` field names match the Python keys
    /// exactly (incl. `rescued_pA`); the `discordant` block is a Rust-only addition
    /// for Feature A (Python has no equivalent), zero-valued unless enabled.
    pub fn to_json(&self) -> String {
        format!(
            "{{\n  \"left\": {},\n  \"right\": {},\n  \"discordant\": {}\n}}\n",
            side_json(&self.left),
            side_json(&self.right),
            self.discordant_json()
        )
    }

    fn discordant_json(&self) -> String {
        format!(
            "{{\"obs\": {}, \"clusters\": {}, \"paired\": {}, \"rejected_rte\": {}}}",
            self.disc_obs, self.disc_clusters, self.disc_paired, self.disc_rejected_rte
        )
    }
}

fn side_json(s: &SideStats) -> String {
    format!(
        "{{\"too_few\": {}, \"rescued_pA\": {}, \"excluded\": {}, \"too_few_after_filter\": {}, \
         \"clipped_failed\": {}, \"unclipped_failed\": {}, \"polymer\": {}, \"passed\": {}}}",
        s.too_few,
        s.rescued_pa,
        s.excluded,
        s.too_few_after_filter,
        s.clipped_failed,
        s.unclipped_failed,
        s.polymer,
        s.passed
    )
}
