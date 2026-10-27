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

/// Per-category heap bytes of a set of breakpoints. Every field that owns heap memory
/// is attributed somewhere, so `total()` is the real (allocator-rounding aside) live
/// footprint — the earlier estimate ignored `ev`, `agree` and the mate qname strings,
/// which is exactly where a PENDING disc_agree junction keeps most of its bytes.
#[derive(Default, Clone, Copy)]
pub struct BpBytes {
    pub n: usize,
    /// breakpoints held PENDING by `disc_agree_second_fragment` (`agree.is_some()`)
    pub n_pending: usize,
    /// `size_of::<Breakpoint>()` per element (the Vec's own storage, len only)
    pub structs: usize,
    /// clipped + unclipped seq/qual buffers
    pub consensus: usize,
    /// query_name + reference_name strings
    pub names: usize,
    /// `mates` (qname strings) + `mate_seqs` + `mate_dests`
    pub mates: usize,
    /// `ev`: EvExtra box + clip_lite/disc_lite/short_lite + captured EvRecs
    pub ev: usize,
    pub ev_clip_lite: usize,
    pub ev_disc_lite: usize,
    pub ev_short_lite: usize,
    /// `agree`: AgreePending box + inserts + clips + anchors
    pub agree: usize,
    pub agree_anchors: usize,
}

impl BpBytes {
    pub fn total(&self) -> usize {
        self.structs + self.consensus + self.names + self.mates + self.ev + self.agree
    }
    pub fn add(&mut self, o: &BpBytes) {
        self.n += o.n;
        self.n_pending += o.n_pending;
        self.structs += o.structs;
        self.consensus += o.consensus;
        self.names += o.names;
        self.mates += o.mates;
        self.ev += o.ev;
        self.ev_clip_lite += o.ev_clip_lite;
        self.ev_disc_lite += o.ev_disc_lite;
        self.ev_short_lite += o.ev_short_lite;
        self.agree += o.agree;
        self.agree_anchors += o.agree_anchors;
    }
}

fn qs_bytes(s: &crate::qseq::QualitySeq) -> usize {
    s.seq.capacity() + s.qual.capacity()
}

/// Full heap footprint of one breakpoint, by category.
pub fn breakpoint_breakdown(bp: &Breakpoint) -> BpBytes {
    use std::mem::size_of;
    let mut b = BpBytes { n: 1, ..Default::default() };
    b.structs = size_of::<Breakpoint>();
    b.consensus = qs_bytes(&bp.clipped) + qs_bytes(&bp.unclipped);
    b.names = bp.query_name.as_ref().map_or(0, |s| s.capacity()) + bp.reference_name.capacity();
    b.mates = bp.mates.capacity() * size_of::<(bool, String)>()
        + bp.mates.iter().map(|(_, q)| q.capacity()).sum::<usize>()
        + bp.mate_seqs.capacity() * size_of::<crate::qseq::QualitySeq>()
        + bp.mate_seqs.iter().map(qs_bytes).sum::<usize>()
        + bp.mate_dests.capacity() * size_of::<(Option<usize>, i64)>();
    if let Some(ev) = bp.ev.as_deref() {
        let rec = |v: &[crate::evidence::EvRec]| -> usize {
            v.iter().map(|r| r.cigar.len() + r.seq.len() + r.qual.len()).sum::<usize>()
        };
        b.ev_clip_lite = ev.clip_lite.capacity() * size_of::<crate::evidence::ClipLite>();
        b.ev_disc_lite = ev.disc_lite.capacity() * size_of::<crate::evidence::DiscLite>();
        b.ev_short_lite = ev.short_lite.capacity() * size_of::<crate::evidence::ShortLite>();
        b.ev = size_of::<crate::evidence::EvExtra>() + b.ev_clip_lite + b.ev_disc_lite + b.ev_short_lite;
        if let Some(rc) = ev.recs.as_deref() {
            b.ev += size_of::<crate::evidence::EvRecs>()
                + (rc.clip.capacity() + rc.disc.capacity() + rc.short.capacity() + rc.mates.capacity())
                    * size_of::<crate::evidence::EvRec>()
                + rec(&rc.clip)
                + rec(&rc.disc)
                + rec(&rc.short)
                + rec(&rc.mates);
        }
    }
    if let Some(ag) = bp.agree.as_deref() {
        b.n_pending = 1;
        b.agree_anchors = ag.anchors.capacity() * size_of::<crate::evidence::DiscLite>();
        b.agree = size_of::<crate::model::AgreePending>()
            + ag.inserts.capacity() * size_of::<Vec<u8>>()
            + ag.inserts.iter().map(|v| v.capacity()).sum::<usize>()
            + ag.clips.capacity() * size_of::<crate::model::ClipMol>()
            + b.agree_anchors;
    }
    b
}

/// Per-category sum over a slice (the Vec's own len storage included via `structs`).
pub fn breakdown(v: &[Breakpoint]) -> BpBytes {
    let mut t = BpBytes::default();
    for bp in v {
        t.add(&breakpoint_breakdown(bp));
    }
    t
}

/// Sum of `breakpoint_bytes` over a slice.
pub fn breakpoints_bytes(v: &[Breakpoint]) -> usize {
    breakdown(v).total()
}

/// Unused capacity of a `Vec<Breakpoint>` in bytes (doubling growth leaves up to 1x len).
pub fn slack_bytes(v: &Vec<Breakpoint>) -> usize {
    (v.capacity() - v.len()) * std::mem::size_of::<Breakpoint>()
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
///
/// `collectors` = (sidecar short-overhang kept, buffered, bytes; sc_disc_tmp len) for the
/// per-contig collectors that live only until `cleanup` (0 when off / already drained).
#[allow(clippy::too_many_arguments)]
pub fn report(
    tag: &str,
    tmp: &[Breakpoint],
    final_left: &Vec<Breakpoint>,
    final_right: &Vec<Breakpoint>,
    polya: &[PolyABreakpoint],
    discordant_obs_len: usize,
    discordant_obs_bytes: usize,
    collectors: (usize, usize, usize, usize),
) {
    let rss = rss_bytes()
        .map(|b| format!("rss={:.0}MiB ", mib(b as usize)))
        .unwrap_or_default();
    let mut fin = breakdown(final_left);
    fin.add(&breakdown(final_right));
    let slack = slack_bytes(final_left) + slack_bytes(final_right);
    eprintln!(
        "[mem] {tag}: {rss}\
         tmp_bp={} ({:.1}MiB)  final_L={} final_R={} ({:.1}MiB)  \
         polyA={} ({:.1}MiB)  disc_obs={} ({:.1}MiB)",
        tmp.len(),
        mib(breakpoints_bytes(tmp)),
        final_left.len(),
        final_right.len(),
        mib(fin.total()),
        polya.len(),
        mib(polya_bytes(polya)),
        discordant_obs_len,
        mib(discordant_obs_bytes),
    );
    // Second line: where the final-breakpoint bytes sit, so a blow-up can be attributed to
    // a field rather than to "breakpoints". Per-breakpoint averages are over all finals.
    let per = |b: usize| if fin.n == 0 { 0.0 } else { b as f64 / fin.n as f64 };
    let (sk, sb, sbytes, dl) = collectors;
    eprintln!(
        "[mem]   final breakdown: pending={} ({:.0}%)  struct={:.1}MiB  consensus={:.1}MiB  names={:.1}MiB  \
         mates={:.1}MiB  ev={:.1}MiB (clip_lite {:.1} / disc_lite {:.1} / short_lite {:.1})  \
         agree={:.1}MiB (anchors {:.1})  avg={:.0}B/bp  vec_slack={:.1}MiB  \
         | collectors: short kept={} buf={} ({:.1}MiB) disc_tmp={} ({:.1}MiB)",
        fin.n_pending,
        if fin.n == 0 { 0.0 } else { 100.0 * fin.n_pending as f64 / fin.n as f64 },
        mib(fin.structs),
        mib(fin.consensus),
        mib(fin.names),
        mib(fin.mates),
        mib(fin.ev),
        mib(fin.ev_clip_lite),
        mib(fin.ev_disc_lite),
        mib(fin.ev_short_lite),
        mib(fin.agree),
        mib(fin.agree_anchors),
        per(fin.total()),
        mib(slack),
        sk,
        sb,
        mib(sbytes),
        dl,
        mib(dl * std::mem::size_of::<crate::evidence::DiscLite>()),
    );
}

#[cfg(test)]
mod size_probe {
    #[test]
    fn print_struct_sizes() {
        use std::mem::size_of;
        eprintln!(
            "SIZE Breakpoint={} QualitySeq={} EvExtra={} AgreePending={} ClipMol={} DiscLite={} ClipLite={} ShortLite={} EvRec={} PolyABreakpoint={}",
            size_of::<crate::model::Breakpoint>(),
            size_of::<crate::qseq::QualitySeq>(),
            size_of::<crate::evidence::EvExtra>(),
            size_of::<crate::model::AgreePending>(),
            size_of::<crate::model::ClipMol>(),
            size_of::<crate::evidence::DiscLite>(),
            size_of::<crate::evidence::ClipLite>(),
            size_of::<crate::evidence::ShortLite>(),
            size_of::<crate::evidence::EvRec>(),
            size_of::<crate::polya::PolyABreakpoint>()
        );
    }
}
