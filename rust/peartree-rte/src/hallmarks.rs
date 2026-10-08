//! tools/rte/hallmarks.py -- TSD, EN motif, poly-A, slippage, fold-back.
//!
//! WP-HALL: literal port of hallmarks.py at eaa2718 (split_junction, parse_locus, edge_run,
//! polya_info, _find_near, locate_site, tsd_from_flanks, target_site, en_motif, en_bin,
//! slippage_context, foldback, _longest_prefix_match).
//!
//! Golden: events `hallmark.<function name>` (in: the python arguments, genome as a genome id;
//! out: the return value / the mutated SiteInfo).

use crate::genome::Genome;
use crate::inputs::JunctionEvidence;
use crate::sequtil::{edlib_best, hamming, homopolymer_at, rc};

/// (left_insert, left_flank, right_flank, right_insert), all upper case.
pub type SplitJunction = (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>);

/// `split_junction(left_seq, right_seq)`.
pub fn split_junction(left_seq: &[u8], right_seq: &[u8]) -> SplitJunction {
    let li_end = left_seq.iter().position(|c| c.is_ascii_uppercase()).unwrap_or(left_seq.len());
    let ri_start = right_seq.iter().position(|c| c.is_ascii_lowercase()).unwrap_or(right_seq.len());
    (
        left_seq[..li_end].to_ascii_uppercase(),
        left_seq[li_end..].to_ascii_uppercase(),
        right_seq[..ri_start].to_ascii_uppercase(),
        right_seq[ri_start..].to_ascii_uppercase(),
    )
}

/// `-?\d+-(-?\d+)` matching all of `rest` -> (a, b).
fn split_range(rest: &str) -> Option<(i64, i64)> {
    let b = rest.as_bytes();
    let mut i = usize::from(b.first() == Some(&b'-'));
    let ds = i;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    if i == ds || i >= b.len() || b[i] != b'-' {
        return None;
    }
    let a = rest[..i].parse::<i64>().ok()?;
    let r = &rest[i + 1..];
    let rb = r.as_bytes();
    let d0 = usize::from(rb.first() == Some(&b'-'));
    if rb.len() <= d0 || !rb[d0..].iter().all(|c| c.is_ascii_digit()) {
        return None;
    }
    Some((a, r.parse::<i64>().ok()?))
}

/// `parse_locus(title)` -> (contig, a, b) for `contig:a-b` (a, b may be negative).
/// Python `^(.+):(-?\d+)-(-?\d+)$`: greedy contig (backtracks over every ':' from the right);
/// `$` also matches before one trailing newline; `.` does not match a newline.
pub fn parse_locus(title: &str) -> Option<(String, i64, i64)> {
    let t = title.strip_suffix('\n').unwrap_or(title);
    let bytes = t.as_bytes();
    for colon in (1..bytes.len()).rev() {
        if bytes[colon] != b':' {
            continue;
        }
        if bytes[..colon].contains(&b'\n') {
            continue;
        }
        if let Some((a, b)) = split_range(&t[colon + 1..]) {
            return Some((t[..colon].to_string(), a, b));
        }
    }
    None
}

/// `edge_run(insert, at_end)` -> (base (0 = none), length).
pub fn edge_run(insert: &[u8], at_end: bool) -> (u8, i64) {
    let mut s: Vec<u8> = insert.to_ascii_uppercase();
    if s.is_empty() {
        return (0, 0);
    }
    if at_end {
        s.reverse();
    }
    let mut b = s[0];
    if b != b'A' && b != b'T' {
        // allow one non-A/T base right at the junction (ligation / microhomology)
        if s.len() > 1 && (s[1] == b'A' || s[1] == b'T') {
            s.remove(0);
            b = s[0];
        } else {
            return (0, 0);
        }
    }
    let mut n = 0i64;
    let mut mm = 0;
    for i in 0..s.len() {
        if s[i] == b {
            n = i as i64 + 1;
        } else if mm == 0 && n >= 4 && i + 1 < s.len() && s[i + 1] == b {
            mm = 1;
        } else {
            break;
        }
    }
    (b, n)
}

/// hallmarks.PolyAInfo
#[derive(Clone, Debug, PartialEq)]
pub struct PolyAInfo {
    /// +1 poly-A at LEFT, -1 poly-T at RIGHT, 0 unknown
    pub strand: i32,
    pub source: String,
    /// REF-adjacent homopolymer of the LEFT insert (base 0 = none)
    pub left_run: (u8, i64),
    pub right_run: (u8, i64),
    pub both_sided: bool,
    pub length: f64,
}

impl Default for PolyAInfo {
    fn default() -> Self {
        PolyAInfo { strand: 0, source: "none".into(), left_run: (0, 0), right_run: (0, 0), both_sided: false, length: 0.0 }
    }
}

/// `polya_info(left_insert, right_insert, ev_left, ev_right, min_len=10)`.
pub fn polya_info(
    left_insert: &[u8],
    right_insert: &[u8],
    ev_left: Option<&JunctionEvidence>,
    ev_right: Option<&JunctionEvidence>,
    min_len: f64,
) -> PolyAInfo {
    let lr = edge_run(left_insert, true);
    let rr = edge_run(right_insert, false);
    let mut la: f64 = if lr.0 == b'A' { lr.1 as f64 } else { 0.0 };
    let mut rt: f64 = if rr.0 == b'T' { rr.1 as f64 } else { 0.0 };
    // python: la = max(la, median) -> the median only when strictly greater
    if let Some(e) = ev_left {
        if e.polya_len_median != 0.0 && lr.0 == b'A' && e.polya_len_median > la {
            la = e.polya_len_median;
        }
    }
    if let Some(e) = ev_right {
        if e.polya_len_median != 0.0 && rr.0 == b'T' && e.polya_len_median > rt {
            rt = e.polya_len_median;
        }
    }
    let mut info = PolyAInfo { left_run: lr, right_run: rr, ..Default::default() };
    info.both_sided = la >= min_len && rt >= min_len;
    if la >= 5.0 || rt >= 5.0 {
        if la >= rt {
            info.strand = 1;
            info.source = "polyA_left".into();
            info.length = la;
        } else {
            info.strand = -1;
            info.source = "polyT_right".into();
            info.length = rt;
        }
    }
    info
}

/// hallmarks.SiteInfo
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SiteInfo {
    pub contig: Option<String>,
    pub l: Option<i64>,
    pub r: Option<i64>,
    pub located: bool,
    pub tsd_seq: String,
    /// > 0 duplication, < 0 deletion, 0 blunt, None unknown
    pub tsd_len: Option<i64>,
    pub tsd_verified: bool,
    pub en_motif: String,
    pub en_mismatches: Option<i64>,
    pub slippage: bool,
    pub slippage_detail: String,
}

/// `_find_near`: position of `probe` in genome[lo:hi] (exact, else edlib <= 2 edits); the
/// exact hit nearest to `hint` (first on ties), over ALL overlapping occurrences.
fn find_near(genome: &dyn Genome, contig: &str, probe: &[u8], lo: i64, hi: i64, want_end: bool, hint: i64) -> Option<i64> {
    let win = genome.fetch(contig, lo, hi);
    if win.is_empty() || probe.len() < 12 {
        return None;
    }
    if win.len() >= probe.len() {
        let mut best: Option<(i64, i64)> = None; // (distance, offset)
        for p in 0..=(win.len() - probe.len()) {
            if &win[p..p + probe.len()] == probe {
                let d = (lo + p as i64 - hint).abs();
                if best.is_none_or(|(bd, _)| d < bd) {
                    best = Some((d, p as i64));
                }
            }
        }
        if let Some((_, p)) = best {
            return Some(lo + p + if want_end { probe.len() as i64 } else { 0 });
        }
    }
    let (_, ts, te, _) = edlib_best(probe, &win, 2.0 / probe.len() as f64, false)?;
    Some(lo + if want_end { te } else { ts })
}

/// `locate_site(title, left_flank, right_flank, genome, probe_len=30, slack=1500)`.
pub fn locate_site(title: &str, left_flank: &[u8], right_flank: &[u8], genome: Option<&dyn Genome>) -> SiteInfo {
    let probe_len = 30usize;
    let slack = 1500i64;
    let mut si = SiteInfo::default();
    let loc = parse_locus(title);
    if let Some((c, _, _)) = &loc {
        si.contig = Some(c.clone());
    }
    let (Some(genome), Some((contig, a, b))) = (genome, loc) else { return si };
    let lo = a.min(b) - slack;
    let hi = a.max(b) + slack;
    if !right_flank.is_empty() {
        let probe = &right_flank[right_flank.len().saturating_sub(probe_len)..];
        si.r = find_near(genome, &contig, probe, lo, hi, true, b);
    }
    if !left_flank.is_empty() {
        let probe = &left_flank[..probe_len.min(left_flank.len())];
        si.l = find_near(genome, &contig, probe, lo, hi, false, a);
    }
    si.located = si.l.is_some() && si.r.is_some();
    si
}

/// `tsd_from_flanks(left_flank, right_flank, max_len=80, min_len=4, max_mm=1)` ->
/// (seq, len or None, verified).
pub fn tsd_from_flanks(left_flank: &[u8], right_flank: &[u8]) -> (String, Option<i64>, bool) {
    let (max_len, min_len, max_mm) = (80usize, 4usize, 1usize);
    let mut t = max_len.min(left_flank.len()).min(right_flank.len());
    while t >= min_len {
        if hamming(&right_flank[right_flank.len() - t..], &left_flank[..t]) <= max_mm {
            return (String::from_utf8_lossy(&left_flank[..t]).into_owned(), Some(t as i64), true);
        }
        t -= 1;
    }
    (String::new(), None, false)
}

/// `target_site(si, left_flank, right_flank, genome)` (fills tsd_* in place).
pub fn target_site(si: &mut SiteInfo, left_flank: &[u8], right_flank: &[u8], genome: Option<&dyn Genome>) {
    if let (true, Some(genome)) = (si.located, genome) {
        let (l, r) = (si.l.unwrap(), si.r.unwrap());
        let contig = si.contig.clone().unwrap_or_default();
        let t = r - l;
        si.tsd_len = Some(t);
        if t > 0 {
            si.tsd_seq = String::from_utf8_lossy(&genome.fetch(&contig, l, r)).into_owned();
            let tu = t as usize;
            // verify against the read-derived flanks (<= 1 mismatch)
            if right_flank.len() >= tu && left_flank.len() >= tu {
                si.tsd_verified = hamming(&right_flank[right_flank.len() - tu..], &left_flank[..tu]) <= 1;
            } else {
                si.tsd_verified = true;
            }
        } else {
            si.tsd_seq = String::new();
            si.tsd_verified = true;
        }
        return;
    }
    let (seq, t, ok) = tsd_from_flanks(left_flank, right_flank);
    if t.is_some() {
        si.tsd_seq = seq;
        si.tsd_len = t;
        si.tsd_verified = ok;
    }
}

/// `en_motif(si, strand, genome)` (in place).
pub fn en_motif(si: &mut SiteInfo, strand: i32, genome: Option<&dyn Genome>) {
    let Some(genome) = genome else { return };
    if !si.located || strand == 0 {
        return;
    }
    let contig = si.contig.clone().unwrap_or_default();
    let m = if strand > 0 {
        let l = si.l.unwrap();
        rc(&genome.fetch(&contig, l - 2, l + 5))
    } else {
        let r = si.r.unwrap();
        genome.fetch(&contig, r - 5, r + 2)
    };
    if m.len() != 7 {
        return;
    }
    si.en_motif = format!("{}/{}", String::from_utf8_lossy(&m[..5]), String::from_utf8_lossy(&m[5..]));
    let a = m[1..5].iter().filter(|&&c| c != b'T').count() as i64;
    si.en_mismatches = Some(a + if m[5] == b'A' || m[5] == b'G' { 0 } else { 1 });
}

/// `en_bin(mm)`.
pub fn en_bin(mm: Option<i64>) -> &'static str {
    match mm {
        None => ".",
        Some(m) if m <= 1 => "0-1",
        Some(2) => "2",
        Some(3) => "3",
        Some(_) => "4-5",
    }
}

/// `slippage_context(si, strand, genome, min_run=10)` (in place).
pub fn slippage_context(si: &mut SiteInfo, strand: i32, genome: Option<&dyn Genome>) {
    let min_run = 10usize;
    let Some(genome) = genome else { return };
    if !si.located || strand == 0 {
        return;
    }
    let contig = si.contig.clone().unwrap_or_default();
    let (p, base) = if strand > 0 { (si.l.unwrap(), b'A') } else { (si.r.unwrap(), b'T') };
    let ctx = genome.fetch(&contig, p - 30, p + 30);
    if ctx.len() < 60 {
        return;
    }
    let c = 30i64;
    let other = if base == b'T' { b'A' } else { b'T' };
    for b in [base, other] {
        let right = homopolymer_at(&ctx, c, 1);
        let left = homopolymer_at(&ctx, c - 1, -1);
        for ((hb, n), wh) in [(right, "after"), (left, "before")] {
            if hb == b && n >= min_run {
                si.slippage = true;
                si.slippage_detail = format!("ref {}{} {} 3' breakpoint", b as char, n, wh);
                return;
            }
        }
    }
}

/// `foldback(left_insert, left_flank, right_insert, right_flank, min_len=15, max_mm=2)`.
pub fn foldback(left_insert: &[u8], left_flank: &[u8], right_insert: &[u8], right_flank: &[u8]) -> i64 {
    let (min_len, max_mm) = (15i64, 2usize);
    let mut best = 0i64;
    if !right_insert.is_empty() && !right_flank.is_empty() {
        let r = rc(right_flank);
        best = best.max(longest_prefix_match(right_insert, &r, max_mm));
    }
    if !left_insert.is_empty() && !left_flank.is_empty() {
        let r = rc(left_flank);
        let a: Vec<u8> = left_insert.iter().rev().copied().collect();
        let b: Vec<u8> = r.iter().rev().copied().collect();
        best = best.max(longest_prefix_match(&a, &b, max_mm));
    }
    if best >= min_len {
        best
    } else {
        0
    }
}

/// `_longest_prefix_match(a, b, max_mm)`
fn longest_prefix_match(a: &[u8], b: &[u8], max_mm: usize) -> i64 {
    let mut mm = 0;
    let mut n = 0usize;
    for (i, (x, y)) in a.iter().zip(b).enumerate() {
        if x != y {
            mm += 1;
            if mm > max_mm {
                break;
            }
        }
        n = i + 1;
    }
    // do not count trailing mismatches
    while n > 0 && a[n - 1] != b[n - 1] {
        n -= 1;
    }
    n as i64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn locus_parsing() {
        assert_eq!(parse_locus("chr1:10-20"), Some(("chr1".into(), 10, 20)));
        assert_eq!(parse_locus("a:b:-5--3"), Some(("a:b".into(), -5, -3)));
        assert_eq!(parse_locus("chr1:10"), None);
        assert_eq!(parse_locus(""), None);
    }
}
