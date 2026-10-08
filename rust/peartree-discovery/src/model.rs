//! Port of src/breakpoint.py (Breakpoint + Breakpoint.join).

use crate::config::*;
use crate::evidence::{frag_hash, select_lowest, ClipLite, DiscLite, EvExtra};
use crate::filters::{clean_clipped_seq, find_consensus};
use crate::qseq::QualitySeq;
use crate::stats::Stats;

#[derive(Clone)]
pub struct Breakpoint {
    pub exclude: bool,
    pub is_forward: Option<bool>,
    pub is_read1: Option<bool>,
    pub query_name: Option<String>,
    pub breakpoint: i64,
    pub clipped: QualitySeq,
    pub unclipped: QualitySeq,
    pub reference_name: String,
    pub side: i32,
    pub has_mate: bool,
    pub bp_precise: bool,
    /// mapping quality of the source read (SENS-2 low-MAPQ-clip-ratio guard); 0 for
    /// a synthesised consensus breakpoint.
    pub mapq: u8,
    /// (mate_is_read1, qname)
    pub mates: Vec<(bool, String)>,
    pub mate_seqs: Vec<QualitySeq>,
    /// Feature B: mate landing sites (ref_id, 0-based pos), filled in find_mates when
    /// `splice_hallmark`/`discordant_anchor` is on. Used by the splice/pseudogene
    /// annotation; does not affect the main output.
    pub mate_dests: Vec<(Option<usize>, i64)>,
    pub n_reads: usize,
    /// distinct fragments (qname hashes) among the reads supporting the consensus
    /// position; 1 for a per-read or rescued breakpoint. Used by the one-sided floor.
    pub n_frags: usize,
    /// TPRT sidecar payload (`evidence_sidecar`); None otherwise (8 bytes).
    pub ev: Option<Box<EvExtra>>,
    /// raw SAM flag of the source read (0 for a synthesised consensus breakpoint); the
    /// sidecar uses it to re-find the record in the mate pass.
    pub flag: u16,
    /// raw RNEXT/PNEXT (0-based) of the source read: where its mate record sits
    /// (sidecar indexed fetch); -1 when unset or for a synthesised breakpoint.
    pub mref: i32,
    pub mpos: i64,
    /// distinct fragments whose clip spans the whole poly-A tail into structured sequence
    /// (`spans_polya`); computed only when `one_sided_min_spanning_fragments` > 0
    pub n_frags_span_polya: usize,
    /// `disc_agree_second_fragment`: Some = a single-molecule junction held PENDING until the
    /// mate pass decides whether a discordant pair agrees on its insert (see `AgreePending`).
    /// Always None after `Discovery::resolve_disc_agree` and when the key is off.
    pub agree: Option<Box<AgreePending>>,
}

impl Breakpoint {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        side: i32,
        reference_name: String,
        breakpoint: i64,
        query_name: Option<String>,
        clipped: QualitySeq,
        unclipped: QualitySeq,
        is_read1: Option<bool>,
        is_forward: Option<bool>,
        exclude: bool,
        mapq: u8,
    ) -> Self {
        Breakpoint {
            exclude,
            is_forward,
            is_read1,
            query_name,
            breakpoint,
            clipped,
            unclipped,
            reference_name,
            side,
            has_mate: false,
            bp_precise: true,
            mapq,
            mates: Vec::new(),
            mate_seqs: Vec::new(),
            mate_dests: Vec::new(),
            n_reads: 1,
            n_frags: 1,
            ev: None,
            flag: 0,
            mref: -1,
            mpos: -1,
            n_frags_span_polya: 0,
            agree: None,
        }
    }
}

/// `disc_agree_second_fragment` agreement thresholds -- identical to cluster/somatic_table.py
/// AGREE_MIN_BP / AGREE_MIN_ID (and its k = 12 seed).
pub const AGREE_MIN_BP: usize = 25;
pub const AGREE_MIN_ID: f64 = 0.9;
const AGREE_SEED_K: usize = 12;
/// raw clip sequences kept per pending junction (distinct, longest first)
const AGREE_MAX_INSERTS: usize = 4;
/// clip reads kept per pending junction for the anchor duplicate check
const AGREE_MAX_CLIPS: usize = 32;

/// A clip read of a pending junction, for the "different molecule" check of a discordant
/// anchor: fragment (qname hash), the read's 5' outer (unclipped) end and its mate placement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClipMol {
    pub frag: u64,
    pub outer5: i64,
    pub mref: i32,
    pub mpos: i64,
}

/// Payload of a PENDING single-molecule junction (`disc_agree_second_fragment`).
#[derive(Clone, Debug, Default)]
pub struct AgreePending {
    /// raw inserted (soft-clipped) sequences of the supporting clip reads, site-forward
    /// (= BAM SEQ orientation of a read aligned at the site), each >= AGREE_MIN_BP
    pub inserts: Vec<Vec<u8>>,
    /// every clip read of the cluster (capped), for the duplicate check
    pub clips: Vec<ClipMol>,
    /// candidate discordant anchors (filled in `Discovery::cleanup`)
    pub anchors: Vec<DiscLite>,
    /// an anchor's inside mate agreed with an insert (set in the mate pass)
    pub promoted: bool,
}

impl AgreePending {
    /// True when the anchor `d` is (by the lenient rule) the same molecule as one of the clip
    /// reads: same qname, or (tol > 0) 5' outer end AND mate start each within `tol` bp
    /// (both mates unplaced: the outer end alone decides) -- the `count_frag_keys` rule.
    pub fn same_molecule(&self, d: &DiscLite, tol: i64) -> bool {
        let outer5 = if d.is_reverse() { d.end } else { d.start };
        self.clips.iter().any(|c| {
            c.frag == d.frag
                || (tol > 0
                    && (c.outer5 - outer5).abs() <= tol
                    && if c.mref >= 0 && d.mref >= 0 {
                        c.mref == d.mref && (c.mpos - d.mpos).abs() <= tol
                    } else {
                        c.mref < 0 && d.mref < 0
                    })
        })
    }
}

/// `mate` and `insert` share >= AGREE_MIN_BP bases on one diagonal at >= AGREE_MIN_ID identity
/// over their whole overlap on that diagonal (both site-forward): they agree on what was
/// inserted. Port of cluster/somatic_table.py `agrees` (k = 12 seeds, each diagonal tried once).
pub fn agrees(mate: &[u8], insert: &[u8]) -> bool {
    let k = AGREE_SEED_K;
    if mate.is_empty() || insert.len() < AGREE_MIN_BP || insert.len() < k || mate.len() < k {
        return false;
    }
    let mut seeds: rustc_hash::FxHashMap<&[u8], Vec<usize>> = rustc_hash::FxHashMap::default();
    for i in 0..=(insert.len() - k) {
        seeds.entry(&insert[i..i + k]).or_default().push(i);
    }
    let mut tried: rustc_hash::FxHashSet<i64> = rustc_hash::FxHashSet::default();
    for j in 0..=(mate.len() - k) {
        let Some(is) = seeds.get(&mate[j..j + k]) else { continue };
        for &i in is {
            let d = j as i64 - i as i64;
            if !tried.insert(d) {
                continue;
            }
            let lo = 0i64.max(-d);
            let hi = (insert.len() as i64).min(mate.len() as i64 - d);
            let n = hi - lo;
            if n >= AGREE_MIN_BP as i64 {
                let same = (lo..hi).filter(|&x| insert[x as usize] == mate[(x + d) as usize]).count();
                if same as f64 >= AGREE_MIN_ID * n as f64 {
                    return true;
                }
            }
        }
    }
    false
}

/// Site-forward sequence of a discordant anchor's inside mate: in an FR pair the mate lies
/// opposite to its anchor, so its stored (BAM) SEQ is site-forward when (mate 0x10) !=
/// (anchor 0x10), else the reverse complement (= combine's `allele_forward_seq`).
pub fn mate_site_forward(stored: &[u8], mate_reverse: bool, anchor_reverse: bool) -> Vec<u8> {
    if mate_reverse != anchor_reverse {
        stored.to_vec()
    } else {
        crate::qseq::revcomp_bytes(stored)
    }
}

/// SENS-8: a clipped consensus is a "pure" poly-A/T terminus when >=90% of its
/// bases are A (or >=90% T). Used to relax the clip-length floor for genuine
/// poly-A tails only.
fn is_pure_polya(qs: &QualitySeq) -> bool {
    let n = qs.len();
    if n == 0 {
        return false;
    }
    let a = qs.seq.iter().filter(|&&b| b == b'A' || b == b'a').count();
    let t = qs.seq.iter().filter(|&&b| b == b'T' || b == b't').count();
    (a.max(t) as f64 / n as f64) >= 0.9
}

/// Counter.most_common()[0]: most frequent value, ties broken by first
/// insertion order. Returns (value, count).
fn most_common_first(values: &[i64]) -> (i64, usize) {
    // preserve first-seen order
    let mut order: Vec<i64> = Vec::new();
    let mut counts: Vec<usize> = Vec::new();
    for &v in values {
        if let Some(idx) = order.iter().position(|&x| x == v) {
            counts[idx] += 1;
        } else {
            order.push(v);
            counts.push(1);
        }
    }
    let mut best = 0usize;
    for i in 1..order.len() {
        if counts[i] > counts[best] {
            best = i;
        }
    }
    (order[best], counts[best])
}

/// Outcome of a polyA-tail rescue attempt. Distinguishes the three Python paths so
/// OBS-1 can count them exactly: a rescued breakpoint (`rescued_pA`), a non-polyA
/// solo (`too_few`), and a polyA whose cleaned clip fell below min_clip_len — which
/// Python returns as None **without** incrementing any counter.
enum Rescue {
    Rescued(Breakpoint),
    TooShort,
    NotPolyA,
}

/// polyA-tail rescue for a single (solo or low-support) breakpoint. Mutates and
/// returns the breakpoint if it qualifies. Mirrors the two rescue blocks in
/// Breakpoint.join. Counter bookkeeping is left to the caller (see `join`).
fn rescue_polya(mut bp: Breakpoint) -> Rescue {
    if bp.side == CLIP_LEFT && bp.clipped.pyslice(Some(-8), None).eq_bytes(b"AAAAAAAA") {
        bp.clipped = clean_clipped_seq(&bp.clipped.revcomp());
        if bp.clipped.len() < MIN_CLIP_LEN {
            return Rescue::TooShort;
        }
        if bp.is_forward == Some(false) {
            bp.has_mate = true;
            bp.mates = vec![(!(bp.is_read1.unwrap_or(false)), bp.query_name.clone().unwrap_or_default())];
        }
        Rescue::Rescued(bp)
    } else if bp.side == CLIP_RIGHT && bp.clipped.pyslice(None, Some(8)).eq_bytes(b"TTTTTTTT") {
        bp.clipped = clean_clipped_seq(&bp.clipped);
        if bp.clipped.len() < MIN_CLIP_LEN {
            return Rescue::TooShort;
        }
        if bp.is_forward == Some(true) {
            bp.has_mate = true;
            bp.mates = vec![(!(bp.is_read1.unwrap_or(false)), bp.query_name.clone().unwrap_or_default())];
        }
        Rescue::Rescued(bp)
    } else {
        Rescue::NotPolyA
    }
}

/// Fragment id of a per-read breakpoint (qname hash). Mates and primary +
/// supplementary records of one template share it.
#[inline]
/// True when an oriented (junction-outward) clip SPANS a poly-A tail: it starts with a T run of
/// >= `min_polya` (>= 80 % T, sequencing errors tolerated), the run ENDS inside the read (a
/// 5-base window with <= 1 T), and >= `beyond` bases follow that are structured sequence (no
/// homopolymer >= 8, >= 3 distinct bases) -- the element's 3' end beyond the tail. A reference
/// A-tract slippage clip is poly-A to the end of the read and fails.
pub fn spans_polya(seq: &[u8], min_polya: usize, beyond: usize) -> bool {
    let up = |b: u8| b.to_ascii_uppercase();
    let n = min_polya.max(1);
    if seq.len() < n + beyond || seq[..n].iter().filter(|&&b| up(b) == b'T').count() * 5 < n * 4 {
        return false;
    }
    // end of the run: first position >= n whose next 5 bases hold <= 1 T
    let mut end = None;
    let mut j = n;
    while j + 5 <= seq.len() {
        if seq[j..j + 5].iter().filter(|&&b| up(b) == b'T').count() <= 1 {
            end = Some(j);
            break;
        }
        j += 1;
    }
    let Some(e) = end else { return false };
    let rest = &seq[e..];
    if rest.len() < beyond {
        return false;
    }
    let pre = &rest[..beyond];
    let mut seen = [false; 4];
    let (mut run, mut best, mut last) = (0usize, 0usize, 0u8);
    for &b in pre {
        let b = up(b);
        match b {
            b'A' => seen[0] = true,
            b'C' => seen[1] = true,
            b'G' => seen[2] = true,
            b'T' => seen[3] = true,
            _ => {}
        }
        run = if b == last { run + 1 } else { 1 };
        last = b;
        best = best.max(run);
    }
    best < 8 && seen.iter().filter(|&&x| x).count() >= 3
}

fn bp_frag(bp: &Breakpoint) -> u64 {
    frag_hash(bp.query_name.as_deref().unwrap_or("").as_bytes())
}

/// A supporting clip read's fragment identity for the lenient duplicate collapse:
/// qname hash, clip-side outer end (breakpoint -/+ RAW clip length, taken before cleaning)
/// and the mate's placement (RNEXT/PNEXT; -1 when unplaced).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FragKey {
    pub frag: u64,
    pub end: i64,
    pub mref: i32,
    pub mpos: i64,
}

fn frag_key(bp: &Breakpoint) -> FragKey {
    let clen = bp.clipped.len() as i64;
    let end = if bp.side == CLIP_LEFT { bp.breakpoint - clen } else { bp.breakpoint + clen };
    FragKey { frag: bp_frag(bp), end, mref: bp.mref, mpos: bp.mpos }
}

/// Independent fragments among `keys`: distinct qnames, then (tol > 0) templates whose
/// clip-side end and mate start both agree within `tol` bp are merged (single linkage) as
/// one molecule -- PCR/optical duplicates markdup missed (shifted ends, R1/R2 swapped).
/// tol <= 0 = distinct qnames only (legacy).
pub fn count_frag_keys(keys: &[FragKey], tol: i64) -> usize {
    let mut v: Vec<FragKey> = Vec::with_capacity(keys.len());
    for k in keys {
        if !v.iter().any(|x| x.frag == k.frag) {
            v.push(*k);
        }
    }
    if tol <= 0 || v.len() < 2 {
        return v.len();
    }
    let n = v.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn root(p: &mut [usize], mut i: usize) -> usize {
        while p[i] != i {
            p[i] = p[p[i]];
            i = p[i];
        }
        i
    }
    let dup = |a: &FragKey, b: &FragKey| {
        (a.end - b.end).abs() <= tol
            && if a.mref >= 0 && b.mref >= 0 {
                a.mref == b.mref && (a.mpos - b.mpos).abs() <= tol
            } else {
                a.mref < 0 && b.mref < 0
            }
    };
    for i in 0..n {
        for j in (i + 1)..n {
            if dup(&v[i], &v[j]) {
                let (ri, rj) = (root(&mut parent, i), root(&mut parent, j));
                parent[ri.max(rj)] = ri.min(rj);
            }
        }
    }
    (0..n).filter(|&i| root(&mut parent, i) == i).count()
}

/// Number of distinct fragments in `frags`.
pub fn count_fragments(frags: &[u64]) -> usize {
    let mut v = frags.to_vec();
    v.sort_unstable();
    v.dedup();
    v.len()
}

/// The sidecar identities (qname hash, flag, breakpoint coordinate) of every read of a
/// cluster, capped deterministically at `max_evidence_reads_per_breakpoint` (lowest
/// frag hash, then flag). None when the sidecar is off.
fn merge_ev(breakpoints: &[Breakpoint], cfg: &DiscoveryConfig) -> Option<Box<EvExtra>> {
    if !cfg.evidence_sidecar {
        return None;
    }
    let lite: Vec<ClipLite> = breakpoints
        .iter()
        .map(|bp| ClipLite { frag: bp_frag(bp), flag: bp.flag, pos: bp.breakpoint, mref: bp.mref, mpos: bp.mpos })
        .collect();
    let clip_lite = select_lowest(lite, cfg.max_evidence_reads_per_breakpoint, |r| (r.frag, r.flag));
    Some(Box::new(EvExtra { clip_lite, ..Default::default() }))
}

/// Port of Breakpoint.join. Consumes the group of breakpoints in a <6bp window
/// and returns a single consensus breakpoint, or None if filtered out.
pub fn join(mut breakpoints: Vec<Breakpoint>, cfg: &DiscoveryConfig, evidence_floor: usize, stats: &mut Stats) -> Option<Breakpoint> {
    // TPRT fragment mode: the floor counts distinct fragments, and a single-fragment
    // cluster takes the normal consensus path when the floor admits it.
    let frag_mode = cfg.min_evidence_fragments_per_sample.is_some();
    // disc_agree_second_fragment: a junction with ONE clip molecule against a floor of 2 may be
    // held pending (see `AgreePending`)
    let agree_ok = frag_mode && cfg.disc_agree_second_fragment && evidence_floor == 2;
    // a lone clip read takes the consensus path (and may go pending) only when the legacy
    // lone-read path would emit nothing and its clip is long enough to compare
    let lone_agree = agree_ok
        && breakpoints.len() == 1
        && breakpoints[0].clipped.len() >= AGREE_MIN_BP
        && !(cfg.polya_rescue_min_fragments <= 1 && matches!(rescue_polya(breakpoints[0].clone()), Rescue::Rescued(_)));
    if breakpoints.len() < 2 && !(frag_mode && evidence_floor <= 1) && !lone_agree {
        let solo_ev = merge_ev(&breakpoints, cfg);
        let bp = breakpoints.pop().unwrap();
        let side = bp.side;
        return match rescue_polya(bp) {
            // poly-A rescue fragment floor: a lone read is ONE fragment
            Rescue::Rescued(_) if cfg.polya_rescue_min_fragments > 1 => {
                stats.pa_rescue_floor_rejected += 1;
                None
            }
            Rescue::Rescued(mut b) => {
                stats.side_mut(side).rescued_pa += 1;
                b.ev = solo_ev;
                Some(b)
            }
            // polyA candidate whose clip cleaned below min_clip_len: Python counts nothing
            Rescue::TooShort => None,
            Rescue::NotPolyA => {
                stats.side_mut(side).too_few += 1;
                None
            }
        };
    }

    // if any breakpoint has the exclude flag, poison the whole group
    if breakpoints.iter().any(|b| b.exclude) {
        stats.side_mut(breakpoints[0].side).excluded += 1;
        return None;
    }

    let side = breakpoints[0].side;
    let reference_name = breakpoints[0].reference_name.clone();

    // SENS-2 guard: with a lowered min_mapq, drop a locus dominated by low-MAPQ
    // clipped reads (a low-quality pileup, not a real junction). Ships with the
    // relaxation so recall gains near repeats do not come with low-MAPQ artefacts.
    if let Some(max_ratio) = cfg.max_lowq_clip_ratio {
        let total = breakpoints.len();
        let lowq = breakpoints.iter().filter(|b| b.mapq < cfg.lowq_mapq_threshold).count();
        if total > 0 && (lowq as f64 / total as f64) > max_ratio {
            stats.side_mut(side).excluded += 1;
            return None;
        }
    }

    // fragment identities from the RAW clips (before cleaning shortens them)
    let tol = cfg.dedup_coord_tolerance;
    let raw_keys: Vec<FragKey> = breakpoints.iter().map(frag_key).collect();
    // disc_agree: keep the raw (site-forward) clips and clip-read identities
    let agree_raw: Vec<(Vec<u8>, ClipMol)> = if agree_ok {
        breakpoints
            .iter()
            .map(|bp| {
                let (clen, ulen) = (bp.clipped.len() as i64, bp.unclipped.len() as i64);
                let fwd = bp.is_forward != Some(false);
                let outer5 = match (bp.side == CLIP_LEFT, fwd) {
                    (true, true) => bp.breakpoint - clen,
                    (true, false) => bp.breakpoint + ulen,
                    (false, true) => bp.breakpoint - ulen,
                    (false, false) => bp.breakpoint + clen,
                };
                (bp.clipped.seq.clone(), ClipMol { frag: bp_frag(bp), outer5, mref: bp.mref, mpos: bp.mpos })
            })
            .collect()
    } else {
        Vec::new()
    };
    // clean/orient clipped sequences and collect qc-passing precise positions
    let mut bps: Vec<i64> = Vec::new();
    let mut frags: Vec<FragKey> = Vec::new();
    // index (into `breakpoints`) of each entry of `bps`
    let mut bps_idx: Vec<usize> = Vec::new();
    for (bi, bp) in breakpoints.iter_mut().enumerate() {
        if side == CLIP_LEFT {
            // left clipped is reverse complemented from now on
            bp.clipped = clean_clipped_seq(&bp.clipped.revcomp());
        } else {
            bp.clipped = clean_clipped_seq(&bp.clipped);
        }
        if side == CLIP_LEFT && bp.is_forward == Some(false) {
            bp.has_mate = true;
        }
        if side == CLIP_RIGHT && bp.is_forward == Some(true) {
            bp.has_mate = true;
        }
        if bp.bp_precise && bp.clipped.len() >= cfg.min_good_bases {
            bps.push(bp.breakpoint);
            frags.push(raw_keys[bi]);
            bps_idx.push(bi);
        }
    }
    if bps.is_empty() {
        stats.side_mut(side).too_few_after_filter += 1;
        return None;
    }

    let (best_bp, exact_n) = most_common_first(&bps);
    // SENS-1/OBS-3: count support within +/- evidence_window of the mode, not only
    // at the exact modal position. 0 = exact (legacy, byte-identical).
    let n = if frag_mode {
        // distinct fragments among the supporting reads (window or exact mode)
        let support: Vec<FragKey> = bps
            .iter()
            .zip(&frags)
            .filter(|(&b, _)| if cfg.evidence_window > 0 { (b - best_bp).abs() <= cfg.evidence_window } else { b == best_bp })
            .map(|(_, &f)| f)
            .collect();
        count_frag_keys(&support, tol)
    } else if cfg.evidence_window > 0 {
        bps.iter().filter(|&&b| (b - best_bp).abs() <= cfg.evidence_window).count()
    } else {
        exact_n
    };

    let in_support = |b: i64| if cfg.evidence_window > 0 { (b - best_bp).abs() <= cfg.evidence_window } else { b == best_bp };
    // disc_agree: one clip molecule, floor 2, and a supporting clip long enough to compare
    let mut pending = agree_ok
        && n == 1
        && bps.iter().zip(&bps_idx).any(|(&b, &bi)| in_support(b) && agree_raw[bi].0.len() >= AGREE_MIN_BP);
    let is_pa = |bp: &Breakpoint| {
        (bp.side == CLIP_LEFT && bp.clipped.pyslice(Some(-8), None).eq_bytes(b"AAAAAAAA"))
            || (bp.side == CLIP_RIGHT && bp.clipped.pyslice(None, Some(8)).eq_bytes(b"TTTTTTTT"))
    };
    // poly-A rescue fragment floor: distinct fragments among the cluster's poly-A reads
    // (the reads the rescued breakpoint stands for)
    let pa_frags_ok = n >= evidence_floor || cfg.polya_rescue_min_fragments == 0 || {
        let f: Vec<FragKey> =
            breakpoints.iter().zip(&raw_keys).filter(|(bp, _)| is_pa(bp)).map(|(_, &k)| k).collect();
        count_frag_keys(&f, tol) >= cfg.polya_rescue_min_fragments
    };
    if pending {
        // the legacy outcome wins when it emits something: a successful poly-A rescue
        if let Some(bp) = breakpoints.iter().find(|bp| is_pa(bp)) {
            if pa_frags_ok && matches!(rescue_polya(bp.clone()), Rescue::Rescued(_)) {
                pending = false;
            }
        }
    }
    if n < evidence_floor && !pending {
        // the sidecar keeps every read of the cluster on the rescued breakpoint
        let group_ev = merge_ev(&breakpoints, cfg);
        // try polyA rescue on the individual breakpoints; first hit wins
        for bp in breakpoints.into_iter() {
            if is_pa(&bp) {
                return match rescue_polya(bp) {
                    Rescue::Rescued(_) if !pa_frags_ok => {
                        stats.pa_rescue_floor_rejected += 1;
                        None
                    }
                    Rescue::Rescued(mut b) => {
                        stats.side_mut(side).rescued_pa += 1;
                        b.ev = group_ev;
                        Some(b)
                    }
                    // guard above guarantees polyA, so NotPolyA is unreachable here
                    Rescue::TooShort | Rescue::NotPolyA => None,
                };
            }
        }
        stats.side_mut(side).too_few_after_filter += 1;
        return None;
    }

    // build clipped/unclipped consensus inputs
    let mut mates: Vec<(bool, String)> = Vec::new();
    let mut clipped: Vec<QualitySeq> = Vec::new();
    let mut unclipped: Vec<QualitySeq> = Vec::new();
    // reads actually contributing to the consensus (LEFT delta 0 is pushed twice into
    // `clipped`, so `clipped.len()` over-counts; used only in fragment mode)
    let mut n_used: usize = 0;
    for bp in breakpoints.iter() {
        if bp.has_mate {
            mates.push((!(bp.is_read1.unwrap_or(false)), bp.query_name.clone().unwrap_or_default()));
        }
        if !bp.bp_precise {
            continue;
        }
        let delta_bp = bp.breakpoint - best_bp;
        let clen = bp.clipped.len() as i64;
        let ulen = bp.unclipped.len() as i64;
        let before = clipped.len();
        if side == CLIP_RIGHT {
            if delta_bp == 0 {
                clipped.push(bp.clipped.clone());
                unclipped.push(bp.unclipped.revcomp());
            } else if delta_bp > 0 && delta_bp < clen {
                clipped.push(bp.unclipped.pyslice(Some(-delta_bp as isize), None).concat(&bp.clipped));
                unclipped.push(bp.unclipped.pyslice(None, Some(-delta_bp as isize)).revcomp());
            } else if delta_bp < 0 && delta_bp > -ulen {
                clipped.push(bp.clipped.pyslice(Some(-delta_bp as isize), None));
                unclipped.push(
                    bp.clipped
                        .pyslice(None, Some(-delta_bp as isize))
                        .revcomp()
                        .concat(&bp.unclipped.revcomp()),
                );
            }
        } else if side == CLIP_LEFT {
            clipped.push(bp.clipped.clone());
            unclipped.push(bp.unclipped.clone());
            if delta_bp == 0 {
                clipped.push(bp.clipped.clone());
                unclipped.push(bp.unclipped.clone());
            } else if delta_bp > 0 && delta_bp < ulen {
                clipped.push(bp.clipped.pyslice(Some(delta_bp as isize), None));
                unclipped.push(
                    bp.clipped
                        .pyslice(None, Some(delta_bp as isize))
                        .revcomp()
                        .concat(&bp.unclipped),
                );
            } else if delta_bp < 0 && delta_bp > -clen {
                clipped.push(
                    bp.unclipped
                        .pyslice(None, Some(-delta_bp as isize))
                        .revcomp()
                        .concat(&bp.clipped),
                );
                unclipped.push(bp.unclipped.pyslice(Some(-delta_bp as isize), None));
            }
        }
        if clipped.len() > before {
            n_used += 1;
        }
    }

    let clipped_cons = find_consensus(&clipped, cfg.consensus_tolerant);
    let unclipped_cons = find_consensus(&unclipped, cfg.consensus_tolerant);
    // SENS-8: keep the 12 bp floor, but allow shorter clips when the clipped
    // consensus is a pure poly-A/T terminus (a genuine tail). Off = legacy floor.
    let clip_floor = if cfg.short_polya_clip && is_pure_polya(&clipped_cons) {
        cfg.short_polya_min_clip
    } else {
        MIN_CLIP_LEN
    };
    if clipped_cons.len() <= clip_floor {
        stats.side_mut(side).clipped_failed += 1;
        return None;
    }
    if unclipped_cons.len() <= 40 {
        stats.side_mut(side).unclipped_failed += 1;
        return None;
    }
    // reject n-polymer at the breakpoint (n = 1..=4)
    for nn in 1..5usize {
        let base = unclipped_cons.pyslice(None, Some(nn as isize));
        let reps = 24 / nn;
        let query: Vec<u8> = base.seq.iter().cycle().take(base.seq.len() * reps).copied().collect();
        if unclipped_cons.pyslice(None, Some(24)).eq_bytes(&query) {
            stats.side_mut(side).polymer += 1;
            return None;
        }
    }

    let mut b = Breakpoint::new(
        side,
        reference_name,
        best_bp,
        None,
        clipped_cons,
        unclipped_cons,
        None,
        None,
        false,
        0, // synthesised consensus breakpoint: mapq unused
    );
    b.mates = mates;
    b.n_frags = {
        let support: Vec<FragKey> = bps
            .iter()
            .zip(&frags)
            .filter(|(&x, _)| if cfg.evidence_window > 0 { (x - best_bp).abs() <= cfg.evidence_window } else { x == best_bp })
            .map(|(_, &f)| f)
            .collect();
        count_frag_keys(&support, tol)
    };
    if cfg.one_sided_loci && cfg.one_sided_min_spanning_fragments > 0 {
        // one-sided gate input: fragments whose (oriented, junction-outward) clip runs through
        // the entire poly-A tail into the element, among the reads supporting the consensus
        let span: Vec<FragKey> = breakpoints
            .iter()
            .zip(&raw_keys)
            .filter(|(bp, _)| {
                bp.bp_precise
                    && (if cfg.evidence_window > 0 { (bp.breakpoint - best_bp).abs() <= cfg.evidence_window } else { bp.breakpoint == best_bp })
                    && spans_polya(&bp.clipped.seq, cfg.one_sided_min_polya, cfg.one_sided_span_beyond)
            })
            .map(|(_, &k)| k)
            .collect();
        b.n_frags_span_polya = count_frag_keys(&span, tol);
    }
    // legacy n_reads double-counts LEFT reads at delta 0; it only feeds the Feature A
    // rescue (off by default), so the legacy value is kept unless in fragment mode.
    b.n_reads = if frag_mode { n_used } else { clipped.len() };
    b.ev = merge_ev(&breakpoints, cfg);
    if pending {
        // held until the mate pass; `passed` is counted on promotion
        let mut inserts: Vec<Vec<u8>> = Vec::new();
        for (&x, &bi) in bps.iter().zip(&bps_idx) {
            let ins = &agree_raw[bi].0;
            if in_support(x) && ins.len() >= AGREE_MIN_BP && !inserts.contains(ins) {
                inserts.push(ins.clone());
            }
        }
        inserts.sort_by(|a, c| c.len().cmp(&a.len()).then_with(|| a.cmp(c)));
        inserts.truncate(AGREE_MAX_INSERTS);
        let clips: Vec<ClipMol> = agree_raw.iter().take(AGREE_MAX_CLIPS).map(|(_, m)| *m).collect();
        b.agree = Some(Box::new(AgreePending { inserts, clips, ..Default::default() }));
        stats.disc_agree_pending += 1;
        return Some(b);
    }
    stats.side_mut(side).passed += 1;
    Some(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// deterministic pseudo-random DNA (LCG), so consensus filters are not tripped
    fn dna(seed: u64, n: usize) -> Vec<u8> {
        let mut x = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect()
    }

    fn bp(side: i32, pos: i64, qname: &str, read1: bool, forward: bool) -> Breakpoint {
        let clip = dna(side as u64, 30);
        let flank = dna(10 + side as u64, 80);
        Breakpoint::new(
            side,
            "chr1".into(),
            pos,
            Some(qname.into()),
            QualitySeq::new(clip, vec![30; 30]),
            QualitySeq::new(flank, vec![30; 80]),
            Some(read1),
            Some(forward),
            false,
            60,
        )
    }

    /// RIGHT-clipped read at position 1000 sharing one junction (same clip/flank).
    fn right_bp(qname: &str, read1: bool) -> Breakpoint {
        bp(CLIP_RIGHT, 1000, qname, read1, true)
    }

    fn cfg(frag: Option<usize>) -> DiscoveryConfig {
        DiscoveryConfig { min_evidence_fragments_per_sample: frag, ..DiscoveryConfig::default() }
    }

    /// LEFT-clipped read at `pos` with a clip of `clen` bases and its mate at chr0:`mpos`.
    fn left_read(qname: &str, read1: bool, pos: i64, clen: usize, mpos: i64) -> Breakpoint {
        let clip = dna(7, 60)[60 - clen..].to_vec();
        let mut b = Breakpoint::new(
            CLIP_LEFT,
            "chr1".into(),
            pos,
            Some(qname.into()),
            QualitySeq::new(clip, vec![30; clen]),
            QualitySeq::new(dna(11, 120), vec![30; 120]),
            Some(read1),
            Some(true),
            false,
            60,
        );
        b.mref = 0;
        b.mpos = mpos;
        b
    }

    fn cfg_dedup(frag: usize, tol: i64) -> DiscoveryConfig {
        DiscoveryConfig { min_evidence_fragments_per_sample: Some(frag), dedup_coord_tolerance: tol, ..DiscoveryConfig::default() }
    }

    /// PD51635 chr16:72393518 LEFT junction: two templates, clip 12 vs 15 bp (outer end 3 bp
    /// apart), mates 1 bp apart, R1/R2 swapped -- one PCR molecule markdup missed.
    #[test]
    fn shifted_strand_swapped_duplicate_counts_once() {
        let mut st = Stats::default();
        let g = vec![left_read("A", false, 1000, 12, 1140), left_read("B", true, 1000, 15, 1141)];
        assert!(join(g.clone(), &cfg_dedup(2, 0), 2, &mut st).is_some(), "legacy qname count passes");
        assert!(join(g.clone(), &cfg_dedup(2, 5), 2, &mut st).is_none(), "lenient dedup rejects");
        assert_eq!(join(g, &cfg_dedup(1, 5), 1, &mut st).unwrap().n_frags, 1);
    }

    #[test]
    fn distinct_molecules_still_count_twice() {
        let mut st = Stats::default();
        // same junction, different outer ends (clip 12 vs 40) and mates 200 bp apart
        let g = vec![left_read("A", true, 1000, 12, 1140), left_read("B", true, 1000, 40, 1340)];
        assert_eq!(join(g.clone(), &cfg_dedup(2, 5), 2, &mut st).unwrap().n_frags, 2);
        // identical clip end but mates far apart: two molecules
        let g = vec![left_read("A", true, 1000, 30, 1140), left_read("B", true, 1000, 30, 1300)];
        assert_eq!(join(g, &cfg_dedup(2, 5), 2, &mut st).unwrap().n_frags, 2);
    }

    #[test]
    fn count_frag_keys_rules() {
        let k = |frag, end, mref, mpos| FragKey { frag, end, mref, mpos };
        // same qname (both mates of one template) counts once even when far apart
        assert_eq!(count_frag_keys(&[k(1, 0, 0, 100), k(1, 500, 0, 9)], 5), 1);
        // chain A~B~C within tol each -> one molecule (single linkage)
        assert_eq!(count_frag_keys(&[k(1, 0, 0, 100), k(2, 4, 0, 104), k(3, 8, 0, 108)], 5), 1);
        // mates on different contigs: never merged
        assert_eq!(count_frag_keys(&[k(1, 0, 0, 100), k(2, 0, 3, 100)], 5), 2);
        // both mates unplaced: the clip end decides; one placed, one not: kept apart
        assert_eq!(count_frag_keys(&[k(1, 0, -1, -1), k(2, 2, -1, -1)], 5), 1);
        assert_eq!(count_frag_keys(&[k(1, 0, -1, -1), k(2, 2, 0, 100)], 5), 2);
        // tol 0 = legacy distinct qnames
        assert_eq!(count_frag_keys(&[k(1, 0, 0, 100), k(2, 0, 0, 100)], 0), 2);
    }

    #[test]
    fn count_fragments_dedups() {
        assert_eq!(count_fragments(&[5, 3, 5, 5, 9]), 3);
        assert_eq!(count_fragments(&[]), 0);
    }

    #[test]
    fn both_mates_of_one_fragment_count_once() {
        let mut st = Stats::default();
        // legacy read floor: 2 reads (the two mates) pass
        let g = vec![right_bp("fragA", true), right_bp("fragA", false)];
        assert!(join(g.clone(), &cfg(None), 2, &mut st).is_some());
        // fragment floor 2: one template -> rejected
        assert!(join(g.clone(), &cfg(Some(2)), 2, &mut st).is_none());
        // fragment floor 1: passes, and n_reads counts the 2 distinct reads
        let b = join(g, &cfg(Some(1)), 1, &mut st).unwrap();
        assert_eq!(b.n_reads, 2);
    }

    #[test]
    fn primary_and_supplementary_count_once() {
        let mut st = Stats::default();
        // read1 primary + read1 supplementary (same qname) + an independent fragment
        let g = vec![right_bp("fragA", true), right_bp("fragA", true), right_bp("fragB", true)];
        assert!(join(g.clone(), &cfg(Some(2)), 2, &mut st).is_some());
        assert!(join(g.clone(), &cfg(Some(3)), 3, &mut st).is_none());
        // legacy counts 3 reads
        assert!(join(g, &cfg(None), 3, &mut st).is_some());
    }

    #[test]
    fn single_fragment_cluster_passes_only_in_fragment_mode() {
        let mut st = Stats::default();
        let g = vec![right_bp("solo", true)];
        assert!(join(g.clone(), &cfg(None), 2, &mut st).is_none()); // legacy: too_few
        let b = join(g.clone(), &cfg(Some(1)), 1, &mut st).expect("floor 1 admits one fragment");
        assert_eq!(b.n_reads, 1);
        assert!(join(g, &cfg(Some(2)), 2, &mut st).is_none());
    }

    /// RIGHT-clipped poly-T (poly-A tail) read at 2000: clip = 10 T + structured tail
    fn polyt_bp(qname: &str, read1: bool) -> Breakpoint {
        let mut clip = vec![b'T'; 10];
        clip.extend(dna(77, 20));
        let mut b = bp(CLIP_RIGHT, 2000, qname, read1, true);
        b.clipped = QualitySeq::new(clip, vec![30; 30]);
        b
    }

    fn cfg_floor(frag: Option<usize>, n: usize) -> DiscoveryConfig {
        DiscoveryConfig { polya_rescue_min_fragments: n, ..cfg(frag) }
    }

    #[test]
    fn polya_rescue_floor_rejects_a_lone_read() {
        let mut st = Stats::default();
        let g = vec![polyt_bp("solo", true)];
        // legacy: a lone poly-A read is rescued
        assert!(join(g.clone(), &cfg_floor(None, 0), 2, &mut st).is_some());
        assert_eq!((st.right.rescued_pa, st.pa_rescue_floor_rejected), (1, 0));
        // floor 2: one fragment -> refused
        assert!(join(g.clone(), &cfg_floor(None, 2), 2, &mut st).is_none());
        assert_eq!((st.right.rescued_pa, st.pa_rescue_floor_rejected), (1, 1));
        // floor 1 admits one fragment
        assert!(join(g, &cfg_floor(None, 1), 2, &mut st).is_some());
    }

    #[test]
    fn polya_rescue_floor_counts_distinct_fragments_of_the_subfloor_cluster() {
        let mut st = Stats::default();
        // fragment floor 3: two poly-A reads form a sub-floor cluster -> rescue path
        let two = vec![polyt_bp("fragA", true), polyt_bp("fragB", true)];
        assert!(join(two, &cfg_floor(Some(3), 2), 3, &mut st).is_some());
        // two records of ONE template (mate / supplementary) are one fragment -> refused
        let same = vec![polyt_bp("fragA", true), polyt_bp("fragA", false)];
        assert!(join(same.clone(), &cfg_floor(Some(3), 2), 3, &mut st).is_none());
        assert_eq!(st.pa_rescue_floor_rejected, 1);
        // key off: legacy rescue
        assert!(join(same, &cfg_floor(Some(3), 0), 3, &mut st).is_some());
    }

    #[test]
    fn left_n_reads_double_count_fixed_only_in_fragment_mode() {
        let mut st = Stats::default();
        let g = vec![bp(CLIP_LEFT, 500, "a", true, false), bp(CLIP_LEFT, 500, "b", true, false)];
        assert_eq!(join(g.clone(), &cfg(None), 2, &mut st).unwrap().n_reads, 4); // legacy x2
        assert_eq!(join(g, &cfg(Some(2)), 2, &mut st).unwrap().n_reads, 2);
    }

    fn cfg_agree(on: bool) -> DiscoveryConfig {
        DiscoveryConfig { disc_agree_second_fragment: on, ..cfg_dedup(2, 5) }
    }

    #[test]
    fn agrees_mirrors_somatic_table() {
        let ins = dna(5, 60);
        // exact overlap of 30 bp (mate starts inside the insert)
        let mut mate = ins[30..].to_vec();
        mate.extend(dna(6, 60));
        assert!(agrees(&mate, &ins));
        // 24 bp overlap: too short
        let mut mate = ins[36..].to_vec();
        mate.extend(dna(6, 60));
        assert!(!agrees(&mate, &ins));
        // 2 mismatches over a 30 bp overlap (93 %) agree; 4 (87 %) do not
        let mut mate = ins[30..].to_vec();
        for i in [3usize, 20] {
            mate[i] = if mate[i] == b'A' { b'C' } else { b'A' };
        }
        assert!(agrees(&mate, &ins));
        for i in [8usize, 27] {
            mate[i] = if mate[i] == b'A' { b'C' } else { b'A' };
        }
        assert!(!agrees(&mate, &ins));
        // reverse complement: different orientation, no agreement
        assert!(!agrees(&crate::qseq::revcomp_bytes(&ins), &ins));
        // an insert shorter than 25 bp never agrees
        assert!(!agrees(&ins, &ins[..24]));
        assert!(agrees(&ins, &ins[..25]));
    }

    #[test]
    fn mate_site_forward_follows_the_fr_pair() {
        let s = b"AACCGT".to_vec();
        assert_eq!(mate_site_forward(&s, true, false), s);
        assert_eq!(mate_site_forward(&s, false, true), s);
        assert_eq!(mate_site_forward(&s, true, true), b"ACGGTT".to_vec());
        assert_eq!(mate_site_forward(&s, false, false), b"ACGGTT".to_vec());
    }

    #[test]
    fn one_clip_molecule_goes_pending_only_with_the_key() {
        let mut st = Stats::default();
        // one molecule = one template with a duplicate shifted 3/1 bp (lenient collapse)
        let g = vec![left_read("A", false, 1000, 30, 1140), left_read("B", true, 1000, 33, 1141)];
        assert!(join(g.clone(), &cfg_agree(false), 2, &mut st).is_none());
        let b = join(g.clone(), &cfg_agree(true), 2, &mut st).expect("pending");
        let a = b.agree.as_ref().unwrap();
        assert_eq!(b.n_frags, 1);
        // raw clips, site-forward (NOT the reverse-complemented LEFT consensus), distinct
        assert_eq!(a.inserts.len(), 2);
        assert!(a.inserts.iter().all(|i| dna(7, 60).ends_with(i)));
        assert_eq!(a.clips.len(), 2);
        assert_eq!((st.disc_agree_pending, st.left.passed), (1, 0));
        // a lone clip read (single record) goes pending as well
        let b = join(vec![left_read("A", true, 1000, 30, 1140)], &cfg_agree(true), 2, &mut st).unwrap();
        assert!(b.agree.is_some());
        // a clip too short to ever agree (< 25 bp) is not held
        assert!(join(vec![left_read("A", true, 1000, 20, 1140)], &cfg_agree(true), 2, &mut st).is_none());
        // only against a floor of exactly 2
        assert!(join(g, &cfg_agree(true), 3, &mut st).is_none());
    }

    #[test]
    fn two_clip_molecules_unaffected_by_the_key() {
        let g = vec![left_read("A", true, 1000, 12, 1140), left_read("B", true, 1000, 40, 1340)];
        let (mut s0, mut s1) = (Stats::default(), Stats::default());
        let off = join(g.clone(), &cfg_agree(false), 2, &mut s0).unwrap();
        let on = join(g, &cfg_agree(true), 2, &mut s1).unwrap();
        assert!(on.agree.is_none());
        assert_eq!((off.n_frags, on.n_frags), (2, 2));
        assert_eq!((off.breakpoint, off.clipped.seq.clone(), off.unclipped.seq.clone()), (on.breakpoint, on.clipped.seq, on.unclipped.seq));
        assert_eq!((s0.left.passed, s1.left.passed, s1.disc_agree_pending), (1, 1, 0));
    }

    #[test]
    fn anchor_duplicate_rule() {
        let a = AgreePending {
            clips: vec![ClipMol { frag: 1, outer5: 1110, mref: 0, mpos: 50_000 }],
            ..Default::default()
        };
        let d = |frag, flag: u16, start, end, mref, mpos| DiscLite { frag, flag, ref_id: 0, start, end, mref, mpos };
        assert!(a.same_molecule(&d(1, 0x10, 5000, 5100, 3, 9), 5), "same qname");
        assert!(a.same_molecule(&d(2, 0x10, 1000, 1113, 0, 49_996), 5), "reverse: end + mate within 5");
        assert!(!a.same_molecule(&d(2, 0x10, 1000, 1113, 0, 49_990), 5), "mate 10 bp away");
        assert!(!a.same_molecule(&d(2, 0x10, 1000, 1120, 0, 50_000), 5), "outer end 10 bp away");
        assert!(a.same_molecule(&d(2, 0x0, 1110, 1200, 0, 50_000), 5), "forward: start is the 5' end -> 1110 matches");
    }
}

#[cfg(test)]
mod span_tests {
    use super::spans_polya;

    #[test]
    fn spans_polya_needs_the_tail_end_and_structured_sequence_beyond() {
        // tail of 15 T, then element 3' end
        assert!(spans_polya(b"TTTTTTTTTTTTTTTGCATTGACCTAGGC", 10, 10));
        // a sequencing error inside the tail is tolerated
        assert!(spans_polya(b"TTTTTTTCTTTTTTTGCATTGACCTAGGC", 10, 10));
        // slippage: poly-T to the end of the read
        assert!(!spans_polya(b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTT", 10, 10));
        // the run ends but too little follows
        assert!(!spans_polya(b"TTTTTTTTTTTTTTTGCATTG", 10, 10));
        // what follows is low complexity (a second homopolymer), not an element
        assert!(!spans_polya(b"TTTTTTTTTTTTTTTGGGGGGGGGGGG", 10, 10));
        // no poly-T start
        assert!(!spans_polya(b"GCATTGACCTAGGCTTTTTTTTTTTTT", 10, 10));
    }

}
