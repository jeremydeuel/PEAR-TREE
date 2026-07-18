//! Debug-gated memory instrumentation. Everything here is a no-op unless the
//! environment variable `PEARTREE_MEM_DEBUG` is set, and it is only ever called
//! at phase boundaries (per-contig cleanup, end of extract, after find_mates,
//! before output) — never in the per-read hot loop — so the cost when enabled is
//! also negligible. Its job is to answer "where does the memory peak come from":
//! it attributes resident bytes to the individual accumulators.

use crate::model::Breakpoint;
use crate::polya::PolyABreakpoint;

/// True when `PEARTREE_MEM_DEBUG` is set (any value). Cheap enough to call at a
/// phase boundary; do not call per read.
pub fn enabled() -> bool {
    std::env::var_os("PEARTREE_MEM_DEBUG").is_some()
}

/// Current resident set size in bytes. Linux only (reads /proc/self/statm); the
/// production cluster is Linux, which is where the OOM happens. Returns None
/// elsewhere (on macOS, use `/usr/bin/time -l` for the peak instead).
pub fn rss_bytes() -> Option<u64> {
    #[cfg(target_os = "linux")]
    {
        let s = std::fs::read_to_string("/proc/self/statm").ok()?;
        let resident_pages: u64 = s.split_whitespace().nth(1)?.parse().ok()?;
        // statm counts in pages; 4096 is the near-universal page size on x86-64 Linux.
        Some(resident_pages * 4096)
    }
    #[cfg(not(target_os = "linux"))]
    {
        None
    }
}

/// Approximate heap footprint of one raw breakpoint: the struct itself plus the
/// four seq/qual buffers, the query name, and any mate buffers.
pub fn breakpoint_bytes(bp: &Breakpoint) -> usize {
    let qs = |s: &crate::qseq::QualitySeq| s.seq.capacity() + s.qual.capacity();
    std::mem::size_of::<Breakpoint>()
        + qs(&bp.clipped)
        + qs(&bp.unclipped)
        + bp.query_name.as_ref().map_or(0, |s| s.capacity())
        + bp.mates.capacity() * std::mem::size_of::<(bool, String)>()
        + bp.mate_seqs.iter().map(qs).sum::<usize>()
        + bp.mate_dests.capacity() * std::mem::size_of::<(Option<usize>, i64)>()
}

/// Sum of `breakpoint_bytes` over a slice (len headers included per element).
pub fn breakpoints_bytes(v: &[Breakpoint]) -> usize {
    v.iter().map(breakpoint_bytes).sum()
}

/// Approximate heap footprint of the polyA accumulator (struct + qname + the
/// optional clipped seq/qual + reference name).
pub fn polya_bytes(v: &[PolyABreakpoint]) -> usize {
    v.iter()
        .map(|p| {
            std::mem::size_of::<PolyABreakpoint>()
                + p.qname.capacity()
                + p.reference_name.as_ref().map_or(0, |s| s.capacity())
                + p.clipped.as_ref().map_or(0, |s| s.seq.capacity() + s.qual.capacity())
        })
        .sum()
}

fn mib(bytes: usize) -> f64 {
    bytes as f64 / (1024.0 * 1024.0)
}

/// Emit a bare RSS reading tagged with a phase name. Used at the top-level phase
/// boundaries in main (after discovery / output / rescue / splice) that the
/// per-structure `report` does not reach — this is where the fragmentation
/// high-water mark actually lands.
pub fn phase(tag: &str) {
    if !enabled() {
        return;
    }
    let rss = rss_bytes()
        .map(|b| format!("{:.0}MiB", mib(b as usize)))
        .unwrap_or_else(|| "n/a".to_string());
    eprintln!("[mem] PHASE {tag}: rss={rss}");
}

/// Emit one phase-boundary line: `tag`, the RSS if available, and the estimated
/// bytes of each accumulator. `tmp` is the peak `temporary_breakpoints` for this
/// tag (pass it before it is drained to see the real peak).
#[allow(clippy::too_many_arguments)]
pub fn report(
    tag: &str,
    tmp: &[Breakpoint],
    final_left: &[Breakpoint],
    final_right: &[Breakpoint],
    polya: &[PolyABreakpoint],
    discordant_obs_len: usize,
    discordant_obs_bytes: usize,
) {
    let rss = rss_bytes()
        .map(|b| format!("rss={:.0}MiB ", mib(b as usize)))
        .unwrap_or_default();
    eprintln!(
        "[mem] {tag}: {rss}\
         tmp_bp={} ({:.1}MiB)  final_L={} final_R={} ({:.1}MiB)  \
         polyA={} ({:.1}MiB)  disc_obs={} ({:.1}MiB)",
        tmp.len(),
        mib(breakpoints_bytes(tmp)),
        final_left.len(),
        final_right.len(),
        mib(breakpoints_bytes(final_left) + breakpoints_bytes(final_right)),
        polya.len(),
        mib(polya_bytes(polya)),
        discordant_obs_len,
        mib(discordant_obs_bytes),
    );
}
