//! Quality-aware affine-gap local aligner with read soft-clipping (SPEC "Read likelihood"). Owner: B.
//!
//! Model (natural-log likelihood units):
//! * the read is GLOBAL: every read base is either aligned (match / mismatch / inserted) or
//!   soft-clipped at one of the two ends at `ln(clip_prob)` per clipped base;
//! * the segment is LOCAL: the aligned part may start and end anywhere in it (free ends);
//! * match `ln(1-e)`, mismatch `ln(e/3)`, `e = min(0.75, e_read + e_seg)`, `e_x = 10^(-q/10)`;
//!   read q is clamped to `[base_q_min, base_q_max]`, segment q only from below (`base_q_min`) so
//!   genome bases at `GENOME_Q` = 60 contribute ~1e-6 (i.e. nothing); `N` on either side: `ln(1/4)`;
//! * a gap of length k costs `open + (k-1)·ext` (Gotoh, three states M / I / D; no I<->D moves);
//!   inserted read bases carry no extra emission term (pair-HMM convention);
//! * homopolymer rule: a DELETION opening at segment base x uses `homopolymer_gap_open_phred`
//!   when x lies in a segment homopolymer run of length >= `homopolymer_min_len`; an INSERTION
//!   opening between segment bases x-1 and x does so when the first inserted read base copies the
//!   base of such a run adjacent to the gap (seg[x-1] or seg[x]) — poly-A length jitter is cheap;
//! * an alignment starts and ends in the match state.
//!
//! Banding: cells with `|j - i - diagonal| <= band_halfwidth` (j = segment base index, i = read base
//! index) are stored in a compact row buffer of width `2w+1` indexed by `k = j - i - diagonal + w`:
//! the diagonal predecessor of `(i, k)` is `(i-1, k)`, the vertical (insertion) one `(i-1, k+1)` and
//! the horizontal (deletion) one `(i, k-1)`. The clipped-prefix start is offered in every cell, so
//! a clipped prefix of length k enters the band at read index k on its own diagonal. The unbanded
//! DP is the same code with a band wide enough to cover the whole matrix (cells outside the
//! segment are skipped, not computed). A 1-byte traceback per band cell recovers the aligned
//! read span and the segment columns.

use crate::config::Config;
use std::cell::RefCell;

/// Result of aligning one read against one segment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Alignment {
    /// natural-log likelihood of the read given the segment (best path)
    pub ll: f64,
    /// read bases aligned (not clipped)
    pub aligned_bases: usize,
    /// segment columns [start, end) covered by the aligned part
    pub seg_start: usize,
    pub seg_end: usize,
}

/// Align `seq`/`qual` (phred) against `seg_seq`/`seg_qual`. `diagonal` = expected segment index
/// of read base 0 (None -> unbanded). Implements the banded DP with full fallback.
#[allow(dead_code)] // single-diagonal convenience API (tests); the scorer uses align_multi_opts
pub fn align(seq: &[u8], qual: &[u8], seg_seq: &[u8], seg_qual: &[u8], diagonal: Option<i64>, cfg: &Config) -> Alignment {
    match diagonal {
        Some(d) => align_multi(seq, qual, seg_seq, seg_qual, &[d], cfg),
        None => align_multi(seq, qual, seg_seq, seg_qual, &[], cfg),
    }
}

/// Like [`align`] with several candidate diagonals (e.g. a read placed inside a target-site
/// duplication matches two genome anchors of one segment): one banded pass per distinct
/// diagonal, the best kept, then ONE unbanded fallback if the best is still below
/// `perfect_ll - fallback_slack_nats`. An empty `diagonals` slice means unbanded.
#[allow(dead_code)] // kept for the single-segment callers / tests
pub(crate) fn align_multi(seq: &[u8], qual: &[u8], seg_seq: &[u8], seg_qual: &[u8], diagonals: &[i64], cfg: &Config) -> Alignment {
    align_multi_opts(seq, qual, seg_seq, seg_qual, diagonals, cfg, cfg.realign_fallback_full).0
}

/// [`align_multi`] with the per-call fallback switchable; also returns `perfect_ll` so the
/// caller can decide on a fallback across SEVERAL segments (readlik.rs: a reference read scored
/// against an alt segment is legitimately far below perfect, that is not a band miss).
pub(crate) fn align_multi_opts(seq: &[u8], qual: &[u8], seg_seq: &[u8], seg_qual: &[u8], diagonals: &[i64], cfg: &Config, fallback: bool) -> (Alignment, f64) {
    WS.with(|ws| {
        let ws = &mut *ws.borrow_mut();
        let p = Params::new(cfg);
        ws.ensure_table(cfg);
        let perfect = ws.prepare_read(seq, qual, cfg);
        if diagonals.is_empty() {
            return (ws.run(seg_seq, seg_qual, None, &p), perfect);
        }
        let mut best: Option<Alignment> = None;
        for (i, &d) in diagonals.iter().enumerate() {
            if diagonals[..i].contains(&d) {
                continue;
            }
            let a = ws.run(seg_seq, seg_qual, Some((d, cfg.band_halfwidth as i64)), &p);
            if best.is_none_or(|b| a.ll > b.ll) {
                best = Some(a);
            }
        }
        let banded = best.expect("non-empty diagonals");
        if fallback && banded.ll < perfect - cfg.fallback_slack_nats {
            let full = ws.run(seg_seq, seg_qual, None, &p);
            if full.ll >= banded.ll {
                return (full, perfect);
            }
        }
        (banded, perfect)
    })
}

const NEG: f32 = f32::NEG_INFINITY;
/// highest phred index of the (read q, segment q) score table
const QMAX_IDX: usize = 93;
const QDIM: usize = QMAX_IDX + 1;

/// ACGT -> 0..3, everything else (N, IUPAC, '=') -> 4. Case-insensitive.
static CODE: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b'a' as usize] = 0;
    t[b'c' as usize] = 1;
    t[b'g' as usize] = 2;
    t[b't' as usize] = 3;
    t
};

#[inline]
fn phred_to_ln(phred: f64) -> f64 {
    -phred / 10.0 * std::f64::consts::LN_10
}

/// Per-call scalar parameters (cheap to derive from the config).
struct Params {
    clip: f32,
    lnq: f32,
    open: f32,
    ext: f32,
    hp_open: f32,
    hp_min: usize,
}

impl Params {
    fn new(cfg: &Config) -> Params {
        Params {
            clip: cfg.clip_prob.ln() as f32,
            lnq: 0.25f64.ln() as f32,
            open: phred_to_ln(cfg.gap_open_phred) as f32,
            ext: phred_to_ln(cfg.gap_ext_phred) as f32,
            hp_open: phred_to_ln(cfg.homopolymer_gap_open_phred) as f32,
            hp_min: cfg.homopolymer_min_len.max(1),
        }
    }
}

/// Reusable buffers (one per thread).
struct Workspace {
    /// `table[rq * QDIM + sq] = [ln(1-e), ln(e/3)]`
    table: Vec<[f32; 2]>,
    table_key: Option<(u8, u8)>,
    /// `ln(1 - 10^(-q/10))` per read phred index (for `perfect_ll`)
    ln_perfect: Vec<f64>,
    qmin: u8,
    qmax: u8,
    // read
    rcode: Vec<u8>,
    rrow: Vec<u32>,
    // columns j in [jmin, jmax] (j = 1-based segment prefix length; column c = j - jmin)
    col_code: Vec<u8>,
    col_q: Vec<u32>,
    col_odel: Vec<f32>,
    col_oins: [Vec<f32>; 5],
    hp: Vec<bool>,
    // DP rows (index p = k + 1; sentinels at 0 and W + 1)
    mp: Vec<f32>,
    ip: Vec<f32>,
    dp: Vec<f32>,
    mc: Vec<f32>,
    ic: Vec<f32>,
    dc: Vec<f32>,
    /// traceback, row-major `tb[(i - 1) * W + k]`: bits 0-1 M predecessor (0 M, 1 I, 2 D,
    /// 3 clipped-prefix start), bit 2 I came from I, bit 3 D came from D
    tb: Vec<u8>,
}

thread_local! {
    static WS: RefCell<Workspace> = RefCell::new(Workspace::new());
}

impl Workspace {
    fn new() -> Workspace {
        Workspace {
            table: Vec::new(),
            table_key: None,
            ln_perfect: Vec::new(),
            qmin: 0,
            qmax: 0,
            rcode: Vec::new(),
            rrow: Vec::new(),
            col_code: Vec::new(),
            col_q: Vec::new(),
            col_odel: Vec::new(),
            col_oins: Default::default(),
            hp: Vec::new(),
            mp: Vec::new(),
            ip: Vec::new(),
            dp: Vec::new(),
            mc: Vec::new(),
            ic: Vec::new(),
            dc: Vec::new(),
            tb: Vec::new(),
        }
    }

    fn ensure_table(&mut self, cfg: &Config) {
        let qmin = cfg.base_q_min.min(QMAX_IDX as u8);
        let qmax = cfg.base_q_max.clamp(qmin, QMAX_IDX as u8);
        if self.table_key == Some((qmin, qmax)) {
            return;
        }
        self.qmin = qmin;
        self.qmax = qmax;
        self.table.clear();
        self.table.resize(QDIM * QDIM, [0.0, 0.0]);
        for rq in 0..QDIM {
            for sq in 0..QDIM {
                let e = (10f64.powf(-(rq as f64) / 10.0) + 10f64.powf(-(sq as f64) / 10.0)).min(0.75);
                self.table[rq * QDIM + sq] = [(1.0 - e).ln() as f32, (e / 3.0).ln() as f32];
            }
        }
        self.ln_perfect = (0..QDIM).map(|q| (1.0 - 10f64.powf(-(q as f64) / 10.0)).max(1e-300).ln()).collect();
        self.table_key = Some((qmin, qmax));
    }

    #[inline]
    fn read_q(&self, q: u8) -> usize {
        q.clamp(self.qmin, self.qmax) as usize
    }

    #[inline]
    fn seg_q(&self, q: u8) -> usize {
        q.clamp(self.qmin, QMAX_IDX as u8) as usize
    }

    /// Encode the read; returns `perfect_ll` = Σ ln(1 - e_read).
    fn prepare_read(&mut self, seq: &[u8], qual: &[u8], _cfg: &Config) -> f64 {
        self.rcode.clear();
        self.rrow.clear();
        let mut perfect = 0.0f64;
        for (i, &b) in seq.iter().enumerate() {
            let rq = self.read_q(qual.get(i).copied().unwrap_or(0));
            self.rcode.push(CODE[b as usize]);
            self.rrow.push((rq * QDIM) as u32);
            perfect += self.ln_perfect[rq];
        }
        perfect
    }

    /// Precompute the per-column arrays for j in [jmin, jmax] (1-based prefix lengths).
    fn prepare_columns(&mut self, seg: &[u8], squal: &[u8], jmin: usize, jmax: usize, p: &Params) {
        let m = seg.len();
        // segment base indices needed: [jmin-1, min(jmax, m-1)]
        let a = jmin - 1;
        let b = jmax.min(m - 1);
        self.hp.clear();
        self.hp.resize(b - a + 1, false);
        let mut pos = a;
        while pos > 0 && seg[pos - 1] == seg[a] {
            pos -= 1;
        }
        while pos <= b {
            let base = seg[pos];
            let mut e = pos + 1;
            while e < m && seg[e] == base {
                e += 1;
            }
            let is_hp = e - pos >= p.hp_min && CODE[base as usize] < 4;
            if is_hp {
                for x in pos.max(a)..e.min(b + 1) {
                    self.hp[x - a] = true;
                }
            }
            pos = e;
        }
        let ncol = jmax - jmin + 1;
        self.col_code.clear();
        self.col_q.clear();
        self.col_odel.clear();
        for v in self.col_oins.iter_mut() {
            v.clear();
        }
        for j in jmin..=jmax {
            let x = j - 1; // segment base consumed by an M/D move into column j
            let code = CODE[seg[x] as usize];
            let sq = self.seg_q(squal.get(x).copied().unwrap_or(0));
            self.col_code.push(code);
            self.col_q.push(sq as u32);
            let hp_left = self.hp[x - a];
            self.col_odel.push(if hp_left { p.hp_open } else { p.open });
            // insertion at the boundary between seg[j-1] and seg[j]
            let left_base = if hp_left { code } else { 255 };
            let right_base = if j < m && self.hp[j - a] { CODE[seg[j] as usize] } else { 255 };
            for (rb, v) in self.col_oins.iter_mut().enumerate() {
                let rb = rb as u8;
                v.push(if rb < 4 && (rb == left_base || rb == right_base) { p.hp_open } else { p.open });
            }
        }
        debug_assert_eq!(self.col_code.len(), ncol);
    }

    fn all_clipped(n: usize, p: &Params) -> Alignment {
        Alignment { ll: (n as f32 * p.clip) as f64, aligned_bases: 0, seg_start: 0, seg_end: 0 }
    }

    /// One DP pass. `band` = Some((diagonal, halfwidth)) or None (whole matrix).
    fn run(&mut self, seg: &[u8], squal: &[u8], band: Option<(i64, i64)>, p: &Params) -> Alignment {
        let n = self.rcode.len();
        let m = seg.len();
        if n == 0 {
            return Alignment { ll: 0.0, aligned_bases: 0, seg_start: 0, seg_end: 0 };
        }
        if m == 0 {
            return Self::all_clipped(n, p);
        }
        let (d, w) = match band {
            Some((d, w)) => (d, w.max(0)),
            None => {
                let w = ((m + n) / 2 + 1) as i64;
                (w + 1 - n as i64, w)
            }
        };
        let width = (2 * w + 1) as usize;
        let jmin = (1 + d - w).max(1);
        let jmax = (n as i64 + d + w).min(m as i64);
        if jmin > jmax {
            return Self::all_clipped(n, p);
        }
        let (jmin, jmax) = (jmin as usize, jmax as usize);
        self.prepare_columns(seg, squal, jmin, jmax, p);

        for v in [&mut self.mp, &mut self.ip, &mut self.dp, &mut self.mc, &mut self.ic, &mut self.dc] {
            v.clear();
            v.resize(width + 2, NEG);
        }
        if self.tb.len() < n * width {
            self.tb.resize(n * width, 0);
        }

        let clip = p.clip;
        let lnq = p.lnq;
        let ext = p.ext;
        let mut best = n as f32 * clip;
        let mut best_cell: Option<(usize, usize)> = None; // (i, p)
        let mut prev_empty = false;

        for i in 1..=n {
            let off = i as i64 + d - w; // j = k + off
            let jlo = off.max(1);
            let jhi = (off + width as i64 - 1).min(m as i64);
            if jlo > jhi {
                if !prev_empty {
                    for v in [&mut self.mp, &mut self.ip, &mut self.dp] {
                        v.fill(NEG);
                    }
                    prev_empty = true;
                }
                continue;
            }
            prev_empty = false;
            let klo = (jlo - off) as usize;
            let len = (jhi - jlo + 1) as usize;
            let plo = klo + 1;
            let clo = jlo as usize - jmin;

            let rc = self.rcode[i - 1];
            let row = &self.table[self.rrow[i - 1] as usize..self.rrow[i - 1] as usize + QDIM];
            let start = (i - 1) as f32 * clip;

            let ccode = &self.col_code[clo..clo + len];
            let cq = &self.col_q[clo..clo + len];
            let codel = &self.col_odel[clo..clo + len];
            let coins = &self.col_oins[rc as usize][clo..clo + len];
            let pm_d = &self.mp[plo..plo + len];
            let pi_d = &self.ip[plo..plo + len];
            let pd_d = &self.dp[plo..plo + len];
            let pm_v = &self.mp[plo + 1..plo + 1 + len];
            let pi_v = &self.ip[plo + 1..plo + 1 + len];
            let (mc, ic, dc) = (&mut self.mc, &mut self.ic, &mut self.dc);
            mc[plo - 1] = NEG;
            ic[plo - 1] = NEG;
            dc[plo - 1] = NEG;
            mc[plo + len] = NEG;
            ic[plo + len] = NEG;
            dc[plo + len] = NEG;
            let mco = &mut mc[plo..plo + len];
            let ico = &mut ic[plo..plo + len];
            let dco = &mut dc[plo..plo + len];
            let tbo = &mut self.tb[(i - 1) * width + klo..(i - 1) * width + klo + len];

            let mut m_left = NEG;
            let mut d_left = NEG;
            let mut row_best = NEG;
            let mut row_best_t = 0usize;
            for t in 0..len {
                let sc = ccode[t];
                let st = row[cq[t] as usize];
                let mut s = if rc == sc { st[0] } else { st[1] };
                if (rc | sc) & 4 != 0 {
                    s = lnq;
                }
                // M
                let mut bm = start;
                let mut tbm = 3u8;
                let pm = pm_d[t];
                if pm >= bm {
                    bm = pm;
                    tbm = 0;
                }
                let pi = pi_d[t];
                if pi > bm {
                    bm = pi;
                    tbm = 1;
                }
                let pd = pd_d[t];
                if pd > bm {
                    bm = pd;
                    tbm = 2;
                }
                let mv = bm + s;
                // I (vertical)
                let ia = pm_v[t] + coins[t];
                let ib = pi_v[t] + ext;
                let (iv, tbi) = if ib > ia { (ib, 4u8) } else { (ia, 0u8) };
                // D (horizontal)
                let da = m_left + codel[t];
                let db = d_left + ext;
                let (dv, tbd) = if db > da { (db, 8u8) } else { (da, 0u8) };
                mco[t] = mv;
                ico[t] = iv;
                dco[t] = dv;
                tbo[t] = tbm | tbi | tbd;
                m_left = mv;
                d_left = dv;
                if mv > row_best {
                    row_best = mv;
                    row_best_t = t;
                }
            }
            let end = row_best + (n - i) as f32 * clip;
            if end > best {
                best = end;
                best_cell = Some((i, plo + row_best_t));
            }
            std::mem::swap(&mut self.mp, &mut self.mc);
            std::mem::swap(&mut self.ip, &mut self.ic);
            std::mem::swap(&mut self.dp, &mut self.dc);
        }

        let Some((i_end, p_end)) = best_cell else {
            return Self::all_clipped(n, p);
        };
        // traceback to the clipped-prefix start
        let mut i = i_end;
        let mut k = p_end - 1;
        let j_end = (k as i64 + i as i64 + d - w) as usize;
        let mut state = 0u8; // 0 M, 1 I, 2 D
        let (start_i, start_j) = loop {
            let byte = self.tb[(i - 1) * width + k];
            match state {
                0 => match byte & 3 {
                    3 => {
                        let j = (k as i64 + i as i64 + d - w) as usize;
                        break (i - 1, j - 1);
                    }
                    s => {
                        state = s;
                        i -= 1;
                    }
                },
                1 => {
                    state = if byte & 4 != 0 { 1 } else { 0 };
                    i -= 1;
                    k += 1;
                }
                _ => {
                    state = if byte & 8 != 0 { 2 } else { 0 };
                    k -= 1;
                }
            }
        };
        Alignment { ll: best as f64, aligned_bases: i_end - start_i, seg_start: start_j, seg_end: j_end }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Simple LCG for deterministic pseudo-random tests.
    struct Lcg(u64);
    impl Lcg {
        pub fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        pub fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
        pub fn seq(&mut self, n: usize) -> Vec<u8> {
            (0..n).map(|_| b"ACGT"[self.below(4) as usize]).collect()
        }
    }

    fn ln_match(q: u8) -> f64 {
        let e = 10f64.powf(-(q as f64) / 10.0) + 1e-6;
        (1.0 - e).ln()
    }

    #[test]
    fn perfect_read_scores_sum_of_matches() {
        let cfg = Config::default();
        let mut r = Lcg(7);
        let seg = r.seq(400);
        let sq = vec![60u8; 400];
        let read = seg[100..251].to_vec();
        let q = vec![30u8; 151];
        let a = align(&read, &q, &seg, &sq, Some(100), &cfg);
        assert!((a.ll - 151.0 * ln_match(30)).abs() < 1e-2, "{a:?}");
        assert_eq!(a.aligned_bases, 151);
        assert_eq!((a.seg_start, a.seg_end), (100, 251));
        let b = align(&read, &q, &seg, &sq, None, &cfg);
        assert!((a.ll - b.ll).abs() < 1e-3);
        assert_eq!((b.seg_start, b.seg_end, b.aligned_bases), (100, 251, 151));
    }

    #[test]
    fn read_overhanging_segment_end_is_clipped() {
        let cfg = Config::default();
        let mut r = Lcg(11);
        let seg = r.seq(200);
        let sq = vec![60u8; 200];
        // last 50 segment bases then 30 foreign bases
        let mut read = seg[150..].to_vec();
        let mut tail = r.seq(30);
        tail[0] = if seg[199] == b'A' { b'C' } else { b'A' };
        read.extend_from_slice(&tail);
        let q = vec![30u8; read.len()];
        let a = align(&read, &q, &seg, &sq, None, &cfg);
        assert_eq!(a.seg_end, 200);
        assert!(a.aligned_bases >= 50 && a.aligned_bases <= 52, "{a:?}");
        let banded = align(&read, &q, &seg, &sq, Some(150), &cfg);
        assert!((a.ll - banded.ll).abs() < 1e-3);
    }

    #[test]
    fn deletion_and_insertion_costs() {
        let cfg = Config::default();
        let open = phred_to_ln(cfg.gap_open_phred);
        let ext = phred_to_ln(cfg.gap_ext_phred);
        // non-homopolymer context: alternate pattern without runs
        let mut r = Lcg(3);
        let mut seg = r.seq(300);
        // break any runs >= 3 to keep the context non-homopolymer
        for i in 2..seg.len() {
            if seg[i] == seg[i - 1] && seg[i] == seg[i - 2] {
                seg[i] = if seg[i] == b'A' { b'C' } else { b'A' };
            }
        }
        let sq = vec![60u8; 300];
        let q = vec![30u8; 150];
        // 3-bp deletion in the read at segment 150
        let mut read = seg[75..150].to_vec();
        read.extend_from_slice(&seg[153..228]);
        let a = align(&read, &q, &seg, &sq, Some(75), &cfg);
        let expect = 150.0 * ln_match(30) + open + 2.0 * ext;
        assert!((a.ll - expect).abs() < 0.05, "{} vs {}", a.ll, expect);
        assert_eq!(a.aligned_bases, 150);
    }
}
