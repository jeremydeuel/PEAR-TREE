//! Indel-aware (homopolymer-compressed star alignment) clip consensus. OWNER: P2.
//!
//! Mirrors src/indel_consensus.py line by line. SPEC.md §4.4 -- the numerics section there is
//! normative: dict-insertion-order maps for column votes, python `max()` first-wins ties,
//! Neumaier `py_sum` for `others`, `py_round` (half-even) for the margin, `py_median`,
//! edlib through `crate::align` (SHW for anchored reads / HW for floating reads, wildcard N pad).
//!
//! Representation notes (none of them changes iteration or float-evaluation order):
//! * `_Col.w` / `.groups` / `.runs` / `.censored` share one insertion-ordered entry per key
//!   (a `runs`/`censored` key is only ever created by the same M op that creates/updates the
//!   `w` key, so the key sets coincide); the gap key is the byte `b'-'`.
//! * group sets only need cardinality: group ids are remapped to dense indices per call and
//!   kept as bitsets.
//! * `_Col.noins_groups` is never read in Python and is not kept.
//! * the orientation-1 RLE of a floating read (`rle(revcomp(seq), qual[::-1])`) is computed once
//!   per call instead of once per vote (it is deterministic).
//! * `cols` (a dict keyed by int column) is a Vec indexed by column with a `present` flag; its
//!   length is `max(cols) + 1`.

use crate::align::{self, Mode, OP_DEL, OP_INS, OP_MATCH, OP_MISMATCH, WILDCARD_N};
use crate::pyfmt::{py_median, py_round, py_sum};

/// python `ClipRead` (indel_consensus.py:48): one read, oriented OUTWARD from the junction.
#[derive(Clone, Debug)]
pub struct ClipRead {
    /// bases as given (rle() uppercases)
    pub seq: Vec<u8>,
    /// one quality per base (phred ints)
    pub qual: Vec<u8>,
    /// independent-fragment cluster id (python `gi`): reads of one cluster count once for depth
    pub group: u32,
    /// `1.0 / len(cluster)`
    pub weight: f64,
    /// false = floating read (a mate) that may start anywhere and either strand
    pub anchored: bool,
}

/// python `ConsensusResult` (indel_consensus.py:61). `Default` == `ConsensusResult()`
/// (empty seq, `stop_reason = "empty"`, no poly-A, beyond_* empty / -1 / 0).
#[derive(Clone, Debug, PartialEq)]
pub struct ConsensusResult {
    /// outward orientation, uppercase
    pub seq: Vec<u8>,
    /// independent fragments per base
    pub depth: Vec<u32>,
    /// `min(93, margin)` per base
    pub score: Vec<u8>,
    /// "empty" | "end" | "disagreement" | "depth"
    pub stop_reason: &'static str,
    /// b'A' / b'T'
    pub polya_base: Option<u8>,
    pub polya_start: i64,
    pub polya_len_median: Option<i64>,
    pub polya_len_min: Option<i64>,
    pub polya_len_max: Option<i64>,
    pub beyond_polya: Vec<u8>,
    pub beyond_start: i64,
    pub beyond_polya_support: i64,
}

impl Default for ConsensusResult {
    fn default() -> Self {
        ConsensusResult {
            seq: Vec::new(),
            depth: Vec::new(),
            score: Vec::new(),
            stop_reason: "empty",
            polya_base: None,
            polya_start: -1,
            polya_len_median: None,
            polya_len_min: None,
            polya_len_max: None,
            beyond_polya: Vec::new(),
            beyond_start: -1,
            beyond_polya_support: 0,
        }
    }
}

impl ConsensusResult {
    /// `polya_len_range`: "" when no poly-A, else "{min}-{max}".
    pub fn polya_len_range(&self) -> String {
        match (self.polya_len_min, self.polya_len_max) {
            (Some(a), Some(b)) => format!("{a}-{b}"),
            _ => String::new(),
        }
    }
}

const ITERATIONS: usize = 4;
const MATE_MIN_OVERLAP: u64 = 20;
const GAP: u8 = b'-';

// ----------------------------------------------------------------------------- RLE

/// `rle(seq, qual)` (indel_consensus.py:84): uppercased run bases, run lengths, mean run
/// quality (`sum/len` as f64).
pub fn rle(seq: &[u8], qual: &[u8]) -> (Vec<u8>, Vec<u32>, Vec<f64>) {
    let r = Rle::new(seq.iter().copied(), qual.iter().copied());
    (r.b, r.l, r.q)
}

struct Rle {
    b: Vec<u8>,
    l: Vec<u32>,
    q: Vec<f64>,
}

impl Rle {
    /// zip(seq.upper(), qual) semantics: stops at the shorter of the two.
    fn new(seq: impl Iterator<Item = u8>, qual: impl Iterator<Item = u8>) -> Rle {
        let mut b: Vec<u8> = Vec::new();
        let mut l: Vec<u32> = Vec::new();
        let mut qs: Vec<u64> = Vec::new();
        let mut prev: Option<u8> = None;
        for (c, q) in seq.zip(qual) {
            let c = c.to_ascii_uppercase();
            if Some(c) == prev {
                *l.last_mut().unwrap() += 1;
                *qs.last_mut().unwrap() += q as u64;
            } else {
                b.push(c);
                l.push(1);
                qs.push(q as u64);
                prev = Some(c);
            }
        }
        // python `q / n`: int / int true division, correctly rounded (== f64 division here)
        let q = qs.iter().zip(&l).map(|(&s, &n)| s as f64 / n as f64).collect();
        Rle { b, l, q }
    }

    fn raw_len(&self) -> u64 {
        self.l.iter().map(|&n| n as u64).sum()
    }
}

/// `revcomp`'s translation table (ACGTNacgtn; anything else unchanged)
fn complement(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'N' => b'N',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        b'n' => b'n',
        x => x,
    }
}

/// `_qw`: quality -> vote weight; Q40 == 1.0, floor 2.
#[inline]
#[allow(clippy::manual_clamp)] // literal min(max(q, 2.0), 40.0)
fn qw(q: f64) -> f64 {
    q.max(2.0).min(40.0) / 40.0
}

// ----------------------------------------------------------------------- alignment

#[derive(Clone, Copy, PartialEq, Eq)]
enum OpKind {
    M,
    I,
    D,
}

#[derive(Clone, Copy)]
struct Op {
    op: OpKind,
    ri: usize,
    col: usize,
}

/// `_align`: RLE read vs RLE consensus padded with `"N" * len(read)`. `false` == (None, None).
fn align_ops(read: &[u8], cons: &[u8], anchored: bool, target: &mut Vec<u8>, out: &mut Vec<Op>) -> bool {
    out.clear();
    target.clear();
    target.extend_from_slice(cons);
    target.resize(cons.len() + read.len(), b'N');
    let mode = if anchored { Mode::Shw } else { Mode::Hw };
    let Some(res) = align::path(read, target, mode, WILDCARD_N) else {
        return false;
    };
    let mut col = res.start as usize;
    let mut ri = 0usize;
    for &o in &res.ops {
        match o {
            OP_MATCH | OP_MISMATCH => {
                out.push(Op { op: OpKind::M, ri, col });
                ri += 1;
                col += 1;
            }
            OP_INS => {
                out.push(Op { op: OpKind::I, ri, col });
                ri += 1;
            }
            OP_DEL => {
                out.push(Op { op: OpKind::D, ri, col });
                col += 1;
            }
            _ => unreachable!("unknown edlib op {o}"),
        }
    }
    true
}

/// Set of dense group indices (python `set(group)`; only the cardinality is ever read).
#[derive(Clone, Default)]
struct GroupSet {
    bits: Vec<u64>,
    n: u32,
}

impl GroupSet {
    #[inline]
    fn insert(&mut self, g: usize, nwords: usize) {
        if self.bits.is_empty() {
            self.bits.resize(nwords, 0);
        }
        let (w, m) = (g >> 6, 1u64 << (g & 63));
        if self.bits[w] & m == 0 {
            self.bits[w] |= m;
            self.n += 1;
        }
    }
}

/// one key of `_Col.w` (+ the same key of `.groups`, `.runs`, `.censored`)
struct BaseVote {
    key: u8,
    w: f64,
    groups: GroupSet,
    runs: Vec<u32>,
    /// 0 == key absent from `censored` (run lengths are >= 1)
    censored: u32,
}

/// one key of `_Col.ins`: `[w, groups, [lens lists]]`, lens flattened with stride `key.len()`
struct InsVote {
    key: Vec<u8>,
    w: f64,
    groups: GroupSet,
    lens: Vec<u32>,
}

#[derive(Default)]
struct Col {
    present: bool,
    w: Vec<BaseVote>,
    ins: Vec<InsVote>,
    noins_w: f64,
}

impl Col {
    fn base(&mut self, key: u8) -> &mut BaseVote {
        let i = match self.w.iter().position(|e| e.key == key) {
            Some(i) => i,
            None => {
                self.w.push(BaseVote { key, w: 0.0, groups: GroupSet::default(), runs: Vec::new(), censored: 0 });
                self.w.len() - 1
            }
        };
        &mut self.w[i]
    }
}

fn col_of(cols: &mut Vec<Col>, c: usize) -> &mut Col {
    if cols.len() <= c {
        cols.resize_with(c + 1, Col::default);
    }
    let x = &mut cols[c];
    x.present = true;
    x
}

/// `_mate_overlap_ok` -> (ok, bp)
fn mate_overlap_ok(ops: &[Op], rb: &[u8], rl: &[u32], cons: &[u8], n_real: usize, min_overlap: u64) -> (bool, u64) {
    let (mut matched, mut mism) = (0u64, 0u64);
    let mut bp = 0u64;
    // M columns are strictly increasing, so python's len(distinct) == matched
    for o in ops {
        if o.op == OpKind::M && o.col < n_real {
            if rb[o.ri] == cons[o.col] {
                matched += 1;
                bp += rl[o.ri] as u64;
            } else {
                mism += 1;
            }
        }
    }
    if bp < min_overlap || matched < 6 {
        return (false, bp);
    }
    (matched as f64 / (matched + mism).max(1) as f64 >= 0.85, bp)
}

/// a read prepared for voting: python `(r, rle(r.seq, r.qual))` plus the precomputed
/// orientation-1 RLE of a floating read
struct Enc<'a> {
    r: &'a ClipRead,
    g: usize,
    fwd: Rle,
    rev: Option<Rle>,
}

struct Scratch {
    target: Vec<u8>,
    ops: Vec<Op>,
    tmp: Vec<Op>,
    good: Vec<bool>,
    nwords: usize,
}

/// `_vote`
fn vote(cons: &[u8], reads: &[Enc], min_overlap: u64, s: &mut Scratch) -> Vec<Col> {
    let n_real = cons.len();
    let mut cols: Vec<Col> = Vec::new();
    let nwords = s.nwords;
    for e in reads {
        let r = e.r;
        if e.fwd.b.is_empty() {
            continue;
        }
        let rle: &Rle;
        if r.anchored {
            if !align_ops(&e.fwd.b, cons, true, &mut s.target, &mut s.ops) {
                continue;
            }
            rle = &e.fwd;
        } else {
            let mut best: Option<(u64, &Rle)> = None;
            for orient in 0..2 {
                let cand = if orient == 0 { &e.fwd } else { e.rev.as_ref().unwrap() };
                if !align_ops(&cand.b, cons, false, &mut s.target, &mut s.tmp) {
                    continue;
                }
                let (ok, bp) = mate_overlap_ok(&s.tmp, &cand.b, &cand.l, cons, n_real, min_overlap);
                if ok && best.is_none_or(|(b0, _)| bp > b0) {
                    best = Some((bp, cand));
                    std::mem::swap(&mut s.ops, &mut s.tmp);
                }
            }
            match best {
                None => continue,
                Some((_, b)) => rle = b,
            }
        }
        let (rb, rl, rq) = (&rle.b, &rle.l, &rle.q);
        let ops = &s.ops;
        let last_ri = rb.len() - 1;
        let n = ops.len();
        // a run length is trusted only if every op within +-2 is a clean match
        s.good.clear();
        s.good.extend(ops.iter().map(|o| o.op == OpKind::M && (o.col >= n_real || rb[o.ri] == cons[o.col])));
        let good = &s.good;
        let mut prev_col_aligned: Option<usize> = None;
        for k in 0..n {
            let Op { op, ri, col: c } = ops[k];
            match op {
                OpKind::M => {
                    let b = rb[ri];
                    let w = r.weight * qw(rq[ri]);
                    let x = col_of(&mut cols, c);
                    let clean = good[k.saturating_sub(2)..(k + 3).min(n)].iter().all(|&g| g);
                    let complete = ri != last_ri && (r.anchored || ri != 0);
                    {
                        let bv = x.base(b);
                        bv.w += w;
                        bv.groups.insert(e.g, nwords);
                        if complete && clean {
                            bv.runs.push(rl[ri]);
                        } else {
                            bv.censored = bv.censored.max(rl[ri]);
                        }
                    }
                    if c >= 1 && prev_col_aligned == Some(c - 1) {
                        x.noins_w += w;
                    }
                    prev_col_aligned = Some(c);
                }
                OpKind::D => {
                    let qn = if ri <= last_ri { rq[ri] } else { rq[last_ri] };
                    let qp = if ri > 0 { rq[ri - 1] } else { qn };
                    let w = r.weight * qw((qn + qp) / 2.0);
                    let x = col_of(&mut cols, c);
                    {
                        let bv = x.base(GAP);
                        bv.w += w;
                        bv.groups.insert(e.g, nwords);
                    }
                    if c >= 1 && prev_col_aligned == Some(c - 1) {
                        x.noins_w += w;
                    }
                    prev_col_aligned = Some(c);
                }
                OpKind::I => {
                    // collect the whole consecutive insertion at its first op
                    if k > 0 && ops[k - 1].op == OpKind::I {
                        continue;
                    }
                    let mut j = k;
                    let mut wsum = 0.0f64;
                    while j < n && ops[j].op == OpKind::I {
                        wsum += qw(rq[ops[j].ri]);
                        j += 1;
                    }
                    // a dangling leading/trailing insertion is not evidence
                    if k == 0 || j >= n || c == 0 {
                        continue;
                    }
                    let nins = j - k;
                    let x = col_of(&mut cols, c);
                    let pos = x
                        .ins
                        .iter()
                        .position(|v| v.key.len() == nins && v.key.iter().zip(&ops[k..j]).all(|(&kb, o)| kb == rb[o.ri]));
                    let iv = match pos {
                        Some(p) => &mut x.ins[p],
                        None => {
                            x.ins.push(InsVote {
                                key: ops[k..j].iter().map(|o| rb[o.ri]).collect(),
                                w: 0.0,
                                groups: GroupSet::default(),
                                lens: Vec::new(),
                            });
                            x.ins.last_mut().unwrap()
                        }
                    };
                    // python: e[0] += r.weight * wsum / len(ins_b)  ==  (w * s) / n
                    iv.w += r.weight * wsum / nins as f64;
                    iv.groups.insert(e.g, nwords);
                    iv.lens.extend(ops[k..j].iter().map(|o| rl[o.ri]));
                    prev_col_aligned = None;
                }
            }
        }
    }
    cols
}

/// one decided column: (base, runlen, depth, margin, complete_runs)
struct Dec {
    b: u8,
    rl: i64,
    d: u32,
    m: i64,
    runs: Vec<i64>,
}

/// python `max(seq, key=f)`: first element with the maximal key
fn first_max<T>(xs: &[T], f: impl Fn(&T) -> f64) -> usize {
    let mut bi = 0;
    let mut bv = f(&xs[0]);
    for (i, x) in xs.iter().enumerate().skip(1) {
        let v = f(x);
        if v > bv {
            bi = i;
            bv = v;
        }
    }
    bi
}

/// `_decide`
fn decide(cols: &[Col], final_: bool, min_depth: usize) -> (Vec<Dec>, &'static str) {
    let mut out: Vec<Dec> = Vec::new();
    let mut reason = "end";
    if !cols.iter().any(|c| c.present) {
        return (out, "empty");
    }
    let mut scratch: Vec<i64> = Vec::new();
    for x in cols {
        if !x.present {
            reason = "end";
            break;
        }
        // insertion before this column (not in the final pass)
        if !final_ && !x.ins.is_empty() {
            let iv = &x.ins[first_max(&x.ins, |v| v.w)];
            if iv.w > x.noins_w && iv.groups.n >= 1 && !out.is_empty() {
                let stride = iv.key.len();
                for (t, &b) in iv.key.iter().enumerate() {
                    scratch.clear();
                    scratch.extend(iv.lens.iter().skip(t).step_by(stride).map(|&l| l as i64));
                    let rl = py_median(&scratch);
                    out.push(Dec { b, rl, d: iv.groups.n, m: 1, runs: Vec::new() });
                }
            }
        }
        if x.w.is_empty() {
            reason = "end";
            break;
        }
        let bi = first_max(&x.w, |v| v.w);
        let bv = &x.w[bi];
        let wb = bv.w;
        let others = py_sum(x.w.iter().enumerate().filter(|&(i, _)| i != bi).map(|(_, v)| v.w));
        let depth = bv.groups.n;
        if final_ {
            if bv.key == GAP {
                if wb > others {
                    continue;
                }
                reason = "disagreement";
                break;
            }
            if wb <= others {
                reason = "disagreement";
                break;
            }
            if (depth as usize) < min_depth {
                reason = "depth";
                break;
            }
        } else if bv.key == GAP {
            continue;
        }
        let cens = if bv.censored == 0 { 1 } else { bv.censored as i64 };
        let (rl, runs) = if !bv.runs.is_empty() {
            let runs: Vec<i64> = bv.runs.iter().map(|&v| v as i64).collect();
            (py_median(&runs), runs)
        } else {
            (cens.max(1), vec![cens])
        };
        let margin = py_round((wb - others) * 40.0).max(1);
        out.push(Dec { b: bv.key, rl, d: depth, m: margin, runs });
    }
    (out, reason)
}

/// `_expand`: merge adjacent same-base columns (a merged run keeps no clean stats)
fn expand(dec: Vec<Dec>) -> Vec<Dec> {
    let mut merged: Vec<Dec> = Vec::with_capacity(dec.len());
    for d in dec {
        match merged.last_mut() {
            Some(p) if p.b == d.b => {
                p.rl += d.rl;
                p.d = p.d.min(d.d);
                p.m = p.m.min(d.m);
                p.runs.clear();
            }
            _ => merged.push(d),
        }
    }
    merged
}

/// `_expand_bases`: re-RLE of the expanded decided bases (bases only)
fn expand_bases(dec: &[Dec]) -> Vec<u8> {
    let mut out: Vec<u8> = Vec::with_capacity(dec.len());
    for d in dec {
        if d.rl <= 0 {
            continue;
        }
        let b = d.b.to_ascii_uppercase();
        if out.last() != Some(&b) {
            out.push(b);
        }
    }
    out
}

/// `_pick_seed`
fn pick_seed(encoded: &[&Rle]) -> Vec<u8> {
    const N_CANDIDATES: usize = 8;
    const N_PROBE: usize = 40;
    let enc: Vec<&Rle> = encoded.iter().copied().filter(|e| !e.b.is_empty()).collect();
    let mut order: Vec<usize> = (0..enc.len()).collect();
    // sorted(..., key=raw length, reverse=True) is stable: ties keep input order
    order.sort_by_key(|&i| std::cmp::Reverse(enc[i].raw_len()));
    let mut cands: Vec<&[u8]> = Vec::new();
    for &i in &order {
        let b = enc[i].b.as_slice();
        if !cands.contains(&b) {
            cands.push(b);
        }
        if cands.len() == N_CANDIDATES {
            break;
        }
    }
    if cands.len() == 1 {
        return cands[0].to_vec();
    }
    assert!(!cands.is_empty(), "_pick_seed: no non-empty anchored read (python raises TypeError here)");
    let bases: Vec<&[u8]> = enc.iter().map(|e| e.b.as_slice()).collect();
    let probe: Vec<&[u8]> = if bases.len() <= N_PROBE {
        bases
    } else {
        let step = (bases.len() / N_PROBE).max(1);
        bases.into_iter().step_by(step).take(N_PROBE).collect()
    };
    let maxp = probe.iter().map(|p| p.len()).max().unwrap();
    let mut target: Vec<u8> = Vec::new();
    let mut best: Option<((i64, i64), &[u8])> = None;
    for &c in &cands {
        target.clear();
        target.extend_from_slice(c);
        target.resize(c.len() + maxp, b'N');
        let mut tot = 0i64;
        for p in &probe {
            let q = &p[..p.len().min(c.len())];
            tot += align::distance(q, &target, Mode::Shw, -1, WILDCARD_N) as i64;
        }
        let key = (tot, -(c.len() as i64));
        if best.is_none_or(|(k, _)| key < k) {
            best = Some((key, c));
        }
    }
    best.unwrap().1.to_vec()
}

/// `_polya`
fn polya(res: &mut ConsensusResult, merged: &[Dec], col_off: &[usize], polya_min_len: usize) {
    let mut best: Option<usize> = None;
    for (k, d) in merged.iter().enumerate() {
        if (d.b == b'A' || d.b == b'T') && d.rl >= polya_min_len as i64 && best.is_none_or(|bk| d.rl > merged[bk].rl) {
            best = Some(k);
        }
    }
    let Some(bk) = best else { return };
    let d = &merged[bk];
    let rl = d.rl;
    res.polya_base = Some(d.b);
    res.polya_start = col_off[bk] as i64;
    res.polya_len_median = Some(rl);
    res.polya_len_min = Some(d.runs.iter().copied().min().unwrap_or(rl));
    res.polya_len_max = Some(d.runs.iter().copied().max().unwrap_or(rl));
    let dep: &[u32];
    if d.b == b'A' {
        let start = col_off[bk] + rl as usize; // <= len(seq): the run is part of seq
        res.beyond_polya = res.seq[start..].to_vec();
        res.beyond_start = start as i64;
        dep = &res.depth[start..];
    } else {
        res.beyond_polya = res.seq[..col_off[bk]].to_vec();
        res.beyond_start = 0;
        dep = &res.depth[..col_off[bk]];
    }
    res.beyond_polya_support = dep.iter().copied().min().map_or(0, |v| v as i64);
}

/// `indel_aware_consensus(reads, min_depth, iterations=4, polya_min_len, mate_min_overlap=20)`
/// (indel_consensus.py:336). Reads are consumed in the given order -- the order is significant
/// (seed choice, vote-map insertion order, float summation order). Pure function.
pub fn indel_aware_consensus(reads: &[ClipRead], min_depth: usize, polya_min_len: usize) -> ConsensusResult {
    if !reads.iter().any(|r| r.anchored && !r.seq.is_empty()) {
        return ConsensusResult::default();
    }
    // dense group indices (only set cardinalities are observable)
    let mut gids: Vec<u32> = reads.iter().filter(|r| !r.seq.is_empty()).map(|r| r.group).collect();
    gids.sort_unstable();
    gids.dedup();
    let enc: Vec<Enc> = reads
        .iter()
        .filter(|r| !r.seq.is_empty())
        .map(|r| Enc {
            r,
            g: gids.binary_search(&r.group).unwrap(),
            fwd: Rle::new(r.seq.iter().copied(), r.qual.iter().copied()),
            rev: if r.anchored {
                None
            } else {
                Some(Rle::new(r.seq.iter().rev().map(|&c| complement(c)), r.qual.iter().rev().copied()))
            },
        })
        .collect();
    let mut s = Scratch {
        target: Vec::new(),
        ops: Vec::new(),
        tmp: Vec::new(),
        good: Vec::new(),
        nwords: gids.len().div_ceil(64),
    };
    let seeds: Vec<&Rle> = enc.iter().filter(|e| e.r.anchored).map(|e| &e.fwd).collect();
    let mut cons = pick_seed(&seeds);
    for _ in 0..ITERATIONS {
        let cols = vote(&cons, &enc, MATE_MIN_OVERLAP, &mut s);
        let (dec, _r) = decide(&cols, false, min_depth);
        let new = expand_bases(&dec);
        if new == cons {
            break;
        }
        cons = new;
    }
    let cols = vote(&cons, &enc, MATE_MIN_OVERLAP, &mut s);
    let (dec, reason) = decide(&cols, true, min_depth);
    drop(cols);
    let merged = expand(dec);
    let mut res = ConsensusResult { stop_reason: reason, ..ConsensusResult::default() };
    let total: usize = merged.iter().map(|d| d.rl as usize).sum();
    res.seq.reserve(total);
    res.depth.reserve(total);
    res.score.reserve(total);
    let mut col_off: Vec<usize> = Vec::with_capacity(merged.len());
    let mut off = 0usize;
    for d in &merged {
        col_off.push(off);
        let rl = d.rl.max(0) as usize;
        off += rl;
        res.seq.extend(std::iter::repeat_n(d.b, rl));
        res.depth.extend(std::iter::repeat_n(d.d, rl));
        res.score.extend(std::iter::repeat_n(d.m.min(93) as u8, rl));
    }
    polya(&mut res, &merged, &col_off, polya_min_len);
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cr(seq: &str, q: u8, group: u32) -> ClipRead {
        ClipRead { seq: seq.as_bytes().to_vec(), qual: vec![q; seq.len()], group, weight: 1.0, anchored: true }
    }

    #[test]
    fn rle_literal() {
        // test_indel_consensus.test_rle
        let (b, l, q) = rle(b"AAACG", &[10, 20, 30, 40, 40]);
        assert_eq!(b, b"ACG");
        assert_eq!(l, vec![3, 1, 1]);
        assert_eq!(q, vec![20.0, 40.0, 40.0]);
    }

    #[test]
    fn empty_and_single() {
        // test_indel_consensus.test_empty_and_single
        assert_eq!(indel_aware_consensus(&[], 2, 8), ConsensusResult::default());
        let r = indel_aware_consensus(&[cr("ACGTACGTAC", 30, 1)], 2, 8);
        assert!(r.seq.is_empty());
        assert_eq!(r.stop_reason, "depth");
    }

    #[test]
    fn genuine_disagreement_stops() {
        // test_indel_consensus.test_genuine_disagreement_stops
        const ELEM: &str = "GGCTCACGCCTGTAATCCCG";
        let a = format!("{ELEM}TTTTGACCA");
        let b = format!("{ELEM}CCGAGGTCA");
        let reads = [cr(&a, 30, 1), cr(&a, 30, 2), cr(&b, 30, 3), cr(&b, 30, 4)];
        let r = indel_aware_consensus(&reads, 2, 8);
        assert_eq!(r.seq, ELEM.as_bytes());
        assert_eq!(r.stop_reason, "disagreement");
    }

    // ---------------------------------------------------------------- python fixtures
    //
    // tests/fixtures/consensus_cases.json.gz is written by
    // tests/fixtures/gen_consensus_fixture.py with the reference venv: the deterministic
    // cases of test/test_indel_consensus.py plus randomised read sets, each with the python
    // ConsensusResult. Set CONSENSUS_FIXTURE=<file.json[.gz]> to diff a larger generated set.

    use serde_json::Value;
    use std::io::Read;

    fn load(path: &str) -> Value {
        let raw = std::fs::read(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let text = if path.ends_with(".gz") {
            let mut s = String::new();
            flate2::read::GzDecoder::new(&raw[..]).read_to_string(&mut s).unwrap();
            s
        } else {
            String::from_utf8(raw).unwrap()
        };
        serde_json::from_str(&text).unwrap()
    }

    fn check_case(c: &Value) -> Result<(), String> {
        let reads: Vec<ClipRead> = c["reads"]
            .as_array()
            .unwrap()
            .iter()
            .map(|r| ClipRead {
                seq: r[0].as_str().unwrap().as_bytes().to_vec(),
                qual: r[1].as_str().unwrap().bytes().map(|b| b - 33).collect(),
                group: r[2].as_u64().unwrap() as u32,
                weight: r[3].as_f64().unwrap(),
                anchored: r[4].as_bool().unwrap(),
            })
            .collect();
        let md = c["min_depth"].as_u64().unwrap() as usize;
        let pml = c["polya_min_len"].as_u64().unwrap() as usize;
        let got = indel_aware_consensus(&reads, md, pml);
        let e = &c["result"];
        let want = ConsensusResult {
            seq: e["seq"].as_str().unwrap().as_bytes().to_vec(),
            depth: e["depth"].as_array().unwrap().iter().map(|v| v.as_u64().unwrap() as u32).collect(),
            score: e["score"].as_array().unwrap().iter().map(|v| v.as_u64().unwrap() as u8).collect(),
            stop_reason: match e["stop_reason"].as_str().unwrap() {
                "empty" => "empty",
                "end" => "end",
                "disagreement" => "disagreement",
                "depth" => "depth",
                x => return Err(format!("unknown stop_reason {x}")),
            },
            polya_base: e["polya_base"].as_str().map(|s| s.as_bytes()[0]),
            polya_start: e["polya_start"].as_i64().unwrap(),
            polya_len_median: e["polya_len_median"].as_i64(),
            polya_len_min: e["polya_len_min"].as_i64(),
            polya_len_max: e["polya_len_max"].as_i64(),
            beyond_polya: e["beyond_polya"].as_str().unwrap().as_bytes().to_vec(),
            beyond_start: e["beyond_start"].as_i64().unwrap(),
            beyond_polya_support: e["beyond_polya_support"].as_i64().unwrap(),
        };
        if got != want {
            return Err(format!("case {}:\n got  {:?}\n want {:?}", c["name"], got, want));
        }
        Ok(())
    }

    fn run_fixture(path: &str) -> usize {
        let v = load(path);
        let cases = v.as_array().unwrap();
        let bad: Vec<String> = cases.iter().filter_map(|c| check_case(c).err()).collect();
        for m in bad.iter().take(5) {
            eprintln!("{m}");
        }
        assert!(bad.is_empty(), "{} / {} cases differ from python", bad.len(), cases.len());
        cases.len()
    }

    #[test]
    fn python_fixture() {
        let n = run_fixture(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/consensus_cases.json.gz"));
        assert!(n > 1000);
    }

    #[test]
    fn python_fixture_external() {
        if let Ok(p) = std::env::var("CONSENSUS_FIXTURE") {
            let n = run_fixture(&p);
            eprintln!("CONSENSUS_FIXTURE: {n} cases identical");
        }
    }

    #[test]
    #[should_panic(expected = "_pick_seed")]
    fn empty_rle_seed_panics() {
        // python: _pick_seed -> TypeError ('NoneType' object is not subscriptable) when every
        // anchored read has a sequence but no qualities (zip -> empty RLE)
        let mut a = cr("ACGT", 30, 1);
        a.qual.clear();
        let mut b = cr("ACGTT", 30, 2);
        b.qual.clear();
        indel_aware_consensus(&[a, b], 2, 8);
    }

    #[test]
    fn is_send_sync() {
        fn f<T: Send + Sync>() {}
        f::<ClipRead>();
        f::<ConsensusResult>();
    }
}
