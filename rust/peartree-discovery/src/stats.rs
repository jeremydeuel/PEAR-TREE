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
    // `polya_rescue_min_fragments` counters (global). Serialised only when the key is on
    // (`polya_floor_on`), so the default stats sidecar is unchanged.
    pub polya_floor_on: bool,
    /// join() poly-A rescues refused: the cluster's poly-A reads are < N fragments
    pub pa_rescue_floor_rejected: u64,
    /// Bp+polyA emissions refused: the poly-A-mate pool is < N fragments
    pub pa_pair_floor_rejected: u64,
    /// Feature A anchor/partner candidates skipped: < N fragments
    pub disc_floor_rejected: u64,
    // `disc_agree_second_fragment` counters (global). Serialised only when the key is active
    // (`disc_agree_on`), so the default stats sidecar is unchanged.
    pub disc_agree_on: bool,
    /// single-molecule junctions that passed every consensus gate and were held pending
    pub disc_agree_pending: u64,
    /// ... dropped at once: no candidate discordant anchor (other molecule, right side/strand)
    pub disc_agree_no_anchor: u64,
    /// ... promoted to 2 fragments: an anchor's inside mate agrees with the clip's insert
    pub disc_agree_promoted: u64,
    /// ... dropped after the mate pass: anchors present, none of their mates agrees
    pub disc_agree_rejected: u64,
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
        self.pa_rescue_floor_rejected += other.pa_rescue_floor_rejected;
        self.pa_pair_floor_rejected += other.pa_pair_floor_rejected;
        self.disc_floor_rejected += other.disc_floor_rejected;
        self.disc_agree_pending += other.disc_agree_pending;
        self.disc_agree_no_anchor += other.disc_agree_no_anchor;
        self.disc_agree_promoted += other.disc_agree_promoted;
        self.disc_agree_rejected += other.disc_agree_rejected;
    }

    /// Serialise to JSON. The `left`/`right` field names match the Python keys
    /// exactly (incl. `rescued_pA`); the `discordant` block is a Rust-only addition
    /// for Feature A (Python has no equivalent), zero-valued unless enabled.
    pub fn to_json(&self) -> String {
        let floor = if self.polya_floor_on {
            format!(
                ",\n  \"polya_rescue_floor\": {{\"rescue_rejected\": {}, \"pair_rejected\": {}, \"disc_rejected\": {}}}",
                self.pa_rescue_floor_rejected, self.pa_pair_floor_rejected, self.disc_floor_rejected
            )
        } else {
            String::new()
        };
        let agree = if self.disc_agree_on {
            format!(
                ",\n  \"disc_agree\": {{\"pending\": {}, \"no_anchor\": {}, \"promoted\": {}, \"rejected\": {}}}",
                self.disc_agree_pending, self.disc_agree_no_anchor, self.disc_agree_promoted, self.disc_agree_rejected
            )
        } else {
            String::new()
        };
        format!(
            "{{\n  \"left\": {},\n  \"right\": {},\n  \"discordant\": {}{}{}\n}}\n",
            side_json(&self.left),
            side_json(&self.right),
            self.discordant_json(),
            floor,
            agree
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
