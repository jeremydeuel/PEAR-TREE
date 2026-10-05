//! Per-junction evaluation and the evidence.tsv row. OWNER: P3.
//!
//! Mirrors `JunctionRecord`, `evaluate_junction`, `_short_overhang_check`, `_mate_inside`,
//! `_clips_start_with_polya` of src/combine_insertions_evidence.py:523-755 and `_render_reads` (:1500). SPEC.md §4.3, §6.4.

use crate::align::{self, Mode};
use crate::config::Config;
use crate::consensus::{indel_aware_consensus, ClipRead, ConsensusResult};
use crate::evidence::dedup::{collapse_fragments, independent_clusters, DedupParams, Fragment};
use crate::evidence::row::{allele_forward_seq, cigar_ref_len, ref_pos_at, EvidenceRow, Role};
use crate::genome::RefFetch;
use crate::model::{InputFile, Interner, Side};
use crate::seq::revcomp;
use rustc_hash::FxHashSet;
use std::fmt::Write as _;

/// `supported` column: "NA" until apply_evidence decides, then 0 / 1.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Supported {
    Na,
    No,
    Yes,
}

/// Where a detached record's reads live (shard store; python `fa_ref`). Opaque to everyone
/// but evidence::store.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadsRef {
    pub chunk: u32,
    pub offset: u64,
    pub len: u32,
}

/// python `JunctionRecord` -- one evidence.tsv row (+ its reads for reads.fa.gz).
#[derive(Clone, Debug)]
pub struct JunctionRecord {
    pub insertion_id: String,
    pub side: Side,
    /// kept rows (unused SHORT fragments removed). None once detached to the store.
    pub rows: Option<Vec<EvidenceRow>>,
    pub reads_ref: Option<ReadsRef>,
    pub n_reads: u32,
    pub n_fragments: u32,
    pub n_independent: u32,
    pub n_samples: u32,
    pub n_mates: u32,
    pub n_duplicates: u32,
    pub n_cross: u32,
    pub n_dup_coord: u32,
    pub n_dup_seq: u32,
    /// "" or the comma-joined sorted member locus ids ("." in the TSV when empty)
    pub member_loci: String,
    pub n_short_used: u32,
    pub n_short_rejected: u32,
    pub n_short_mate_inside: u32,
    pub n_independent_no_short: u32,
    /// SHORT rejection reasons (log only); insertion-ordered (reason, count)
    pub short_reasons: Vec<(&'static str, u32)>,
    pub supported: Supported,
    /// junction reads + mates (TSV)
    pub consensus: ConsensusResult,
    /// junction reads only (combined.txt.gz); `ConsensusResult::default()` unless
    /// `indel_aware_consensus`
    pub combined_consensus: ConsensusResult,
    pub polya_end: bool,
    pub fail_reason: String,
    /// reference part in clip_consensus convention, uppercase (`_aligned_part`)
    pub aligned: Vec<u8>,
}

impl JunctionRecord {
    /// `JunctionRecord.tsv()` (evidence.py:546): the 26 EVIDENCE_TSV_COLUMNS joined by "\t" +
    /// "\n". Byte-exact format in SPEC.md §6.4.
    pub fn tsv(&self) -> String {
        let c = &self.consensus;
        let clip = c.seq.to_ascii_lowercase();
        let (clip_cons, depth, beyond): (Vec<u8>, Vec<u32>, Vec<u8>) = match self.side {
            Side::Left => {
                let mut cc = revcomp(&clip);
                cc.extend_from_slice(&self.aligned);
                (cc, c.depth.iter().rev().copied().collect(), revcomp(&c.beyond_polya))
            }
            Side::Right => {
                let mut cc = self.aligned.clone();
                cc.extend_from_slice(&clip);
                (cc, c.depth.clone(), c.beyond_polya.clone())
            }
        };
        let txt = |b: &[u8]| String::from_utf8_lossy(b).into_owned();
        let mut o = String::with_capacity(256 + 4 * clip_cons.len());
        let _ = write!(
            o,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            self.insertion_id,
            self.side.as_str(),
            self.n_reads,
            self.n_fragments,
            self.n_independent,
            self.n_samples,
            self.n_mates,
            match self.supported {
                Supported::Na => "NA",
                Supported::No => "0",
                Supported::Yes => "1",
            },
            txt(&clip_cons)
        );
        for (i, d) in depth.iter().enumerate() {
            if i > 0 {
                o.push(',');
            }
            let _ = write!(o, "{d}");
        }
        o.push('\t');
        if let Some(m) = c.polya_len_median {
            let _ = write!(o, "{m}");
        }
        let _ = write!(
            o,
            "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            c.polya_len_range(),
            txt(&beyond.to_ascii_lowercase()),
            if beyond.is_empty() { 0 } else { c.beyond_polya_support },
            self.polya_end as u8,
            self.n_duplicates,
            self.n_cross,
            self.fail_reason,
            c.stop_reason,
            self.n_dup_coord,
            self.n_dup_seq,
            if self.member_loci.is_empty() { "." } else { &self.member_loci },
            self.n_short_used,
            self.n_short_rejected,
            self.n_short_mate_inside,
            self.n_independent_no_short
        );
        o
    }

    /// FASTA text of this record's rows (`_render_reads(name, rec)`, evidence.py:1500):
    /// `>{name}|{side}|{role}|{sample}|{frag}|{r12}\n{allele_forward_seq}\n` per row.
    /// Panics if detached (rows None).
    pub fn render_reads(&self, name: &str, files: &[InputFile]) -> String {
        let rows = self.rows.as_ref().expect("render_reads on a detached JunctionRecord");
        let mut o = String::new();
        for r in rows {
            let _ = write!(
                o,
                ">{}|{}|{}|{}|{}|{}\n",
                name,
                self.side.as_str(),
                r.role.as_str(),
                files[r.file as usize].sample,
                r.frag,
                r.r12
            );
            o.push_str(&String::from_utf8_lossy(&allele_forward_seq(r)));
            o.push('\n');
        }
        o
    }
}

/// `evaluate_junction(insertion_id, side, rows, cfg, ref_fetch)` (evidence.py:569).
/// `rows` = the pooled (and re-anchored) rows, in pooling order; they are moved into the record
/// (minus dropped SHORT fragments). `ref_fetch` is only used with `count_short_overhang`.
/// `rec.aligned` / `rec.member_loci` are left empty (the caller sets them).
///
/// Pure function of its arguments (runs on rayon workers). The pooled >=2-fragment GATE is not
/// applied here or anywhere (SPEC.md §0); `n_independent` is only counted.
pub fn evaluate_junction(
    insertion_id: &str,
    side: Side,
    rows: Vec<EvidenceRow>,
    cfg: &Config,
    ref_fetch: Option<&dyn RefFetch>,
    contigs: &Interner,
    files: &[InputFile],
) -> JunctionRecord {
    evaluate_junction_with(insertion_id, side, rows, cfg, ref_fetch, contigs, files, &mut |r, d, p| {
        indel_aware_consensus(r, d, p)
    })
}

/// `evaluate_junction` with the consensus function injected (`consensus(reads, min_depth,
/// polya_min_len)`; called in python's order: TSV consensus, combined consensus, SHORT
/// validation consensus). Unit tests replay python's results through it.
pub(crate) fn evaluate_junction_with(
    insertion_id: &str,
    side: Side,
    rows: Vec<EvidenceRow>,
    cfg: &Config,
    ref_fetch: Option<&dyn RefFetch>,
    contigs: &Interner,
    files: &[InputFile],
    consensus: &mut dyn FnMut(&[ClipRead], usize, usize) -> ConsensusResult,
) -> JunctionRecord {
    let dedup = DedupParams::from_cfg(cfg);
    let min_ind = cfg.min_independent_fragments;
    let polya_min = cfg.polya_min_len;
    let use_short = cfg.count_short_overhang;

    let frags = collapse_fragments(&rows, files);
    // SHORT-only fragments never feed the consensus; validated afterwards, may only ADD support
    let (short, main): (Vec<Fragment>, Vec<Fragment>) =
        frags.into_iter().partition(|f| rows[f.primary].role == Role::Short);
    let mut cl = independent_clusters(&main, &rows, &dedup, contigs);

    // consensus input: every read, weighted 1/|cluster|; group = cluster
    let mut reads: Vec<ClipRead> = Vec::new();
    for (gi, c) in cl.clusters.iter().enumerate() {
        let w = 1.0 / c.len() as f64;
        for &fi in c {
            for &ri in &main[fi].rows {
                let r = &rows[ri];
                match r.role {
                    Role::Short => continue,
                    Role::Clip => {
                        let (s, q) = r.outward_clip();
                        reads.push(ClipRead { seq: s, qual: q, group: gi as u32, weight: w, anchored: true });
                    }
                    _ if !r.seq.is_empty() => {
                        reads.push(ClipRead { seq: r.seq.to_vec(), qual: r.quals(), group: gi as u32, weight: w, anchored: false });
                    }
                    _ => {}
                }
            }
        }
    }
    let cons = consensus(&reads, min_ind, polya_min);
    let all_anchored = reads.iter().all(|r| r.anchored);
    let anchored = || -> Vec<ClipRead> { reads.iter().filter(|r| r.anchored).cloned().collect() };
    let combined = if cfg.indel_aware_consensus {
        if all_anchored {
            cons.clone()
        } else {
            consensus(&anchored(), min_ind, polya_min)
        }
    } else {
        ConsensusResult::default()
    };
    let n_independent_no_short = cl.clusters.len() as u32;

    let mut used: Vec<Fragment> = Vec::new();
    let mut short_reasons: Vec<(&'static str, u32)> = Vec::new();
    let mut n_short_mate_inside = 0u32;
    let (mut n_short_used, mut n_short_rejected) = (0u32, 0u32);
    if use_short && !short.is_empty() {
        let has_clip = main.iter().any(|f| rows[f.primary].role == Role::Clip);
        let anch = anchored();
        let vcons = if anch.is_empty() { Vec::new() } else { consensus(&anch, 1, polya_min).seq };
        for f in &short {
            let reason = short_overhang_check(&rows[f.primary], side, &vcons, has_clip, cfg, ref_fetch, contigs);
            if !reason.is_empty() {
                match short_reasons.iter_mut().find(|(k, _)| *k == reason) {
                    Some((_, n)) => *n += 1,
                    None => short_reasons.push((reason, 1)),
                }
            } else {
                used.push(f.clone());
                n_short_mate_inside += mate_inside(f, &rows, cfg, contigs) as u32;
            }
        }
        n_short_used = used.len() as u32;
        n_short_rejected = (short.len() - used.len()) as u32;
    }
    let mut kept = main;
    kept.extend(used.iter().cloned());
    if !used.is_empty() {
        cl = independent_clusters(&kept, &rows, &dedup, contigs);
    }

    // rows of unused SHORT-only fragments (and their mates) are dropped from every output
    let used_keys: FxHashSet<(u32, &str)> = used.iter().map(|f| (f.file, &*f.frag)).collect();
    let drop: FxHashSet<(u32, &str)> =
        short.iter().map(|f| (f.file, &*f.frag)).filter(|k| !used_keys.contains(k)).collect();
    let keep: Vec<bool> = rows
        .iter()
        .map(|r| {
            let k = (r.file, &*r.frag);
            !drop.contains(&k) && (r.role != Role::Short || used_keys.contains(&k))
        })
        .collect();
    let n_fragments = kept.len() as u32;
    let n_samples = {
        let mut s: Vec<u32> = kept.iter().map(|f| f.file).collect();
        s.sort_unstable();
        s.dedup();
        s.len() as u32
    };
    // polya_end: over the ORIGINAL pooled rows (incl. dropped SHORT)
    let polya_end = cons.polya_base.is_some()
        || rows.iter().any(|r| r.role == Role::PolyA)
        || clips_start_with_polya(rows.iter().filter(|r| r.role == Role::Clip), polya_min);
    let kept_rows: Vec<EvidenceRow> = rows.into_iter().zip(keep).filter_map(|(r, k)| k.then_some(r)).collect();
    let n_mates = kept_rows.iter().filter(|r| r.role == Role::Mate).count() as u32;
    let n_reads = kept_rows.len() as u32 - n_mates;

    JunctionRecord {
        insertion_id: insertion_id.to_string(),
        side,
        rows: Some(kept_rows),
        reads_ref: None,
        n_reads,
        n_fragments,
        n_independent: cl.clusters.len() as u32,
        n_samples,
        n_mates,
        n_duplicates: cl.n_dup,
        n_cross: cl.n_cross,
        n_dup_coord: cl.stats.n_dup_coord,
        n_dup_seq: cl.stats.n_dup_seq,
        member_loci: String::new(),
        n_short_used,
        n_short_rejected,
        n_short_mate_inside,
        n_independent_no_short,
        short_reasons,
        supported: Supported::Na,
        consensus: cons,
        combined_consensus: combined,
        polya_end,
        fail_reason: String::new(),
        aligned: Vec::new(),
    }
}

/// `_short_overhang_check(r, side, cons, has_clip, cfg, ref_fetch)` (evidence.py:665):
/// "" if the SHORT read may count as a fragment, else the rejection reason.
pub(crate) fn short_overhang_check(
    r: &EvidenceRow,
    side: Side,
    cons: &[u8],
    has_clip: bool,
    cfg: &Config,
    ref_fetch: Option<&dyn RefFetch>,
    contigs: &Interner,
) -> &'static str {
    let min_b = cfg.short_overhang_min_bases as i64;
    let min_mm = cfg.short_overhang_min_ref_mismatch;
    if !has_clip {
        return "no_clip_fragment";
    }
    if cons.is_empty() {
        return "no_consensus";
    }
    let seq = r.seq.to_ascii_uppercase();
    let at = r.clip_at;
    if !(0 < at && at < seq.len() as i64) || !r.mapped(contigs) {
        return "bad_record";
    }
    let atu = at as usize;
    let j = ref_pos_at(r.pos, &r.cigar, at);
    let mut over = match side {
        Side::Right => seq[atu..].to_vec(),
        Side::Left => revcomp(&seq[..atu]),
    };
    let n = over.len().min(cons.len());
    if (n as i64) < min_b {
        return "overhang_too_short";
    }
    over.truncate(n);
    let Some(rf) = ref_fetch else { return "no_reference" };
    let ref_len = cigar_ref_len(&r.cigar);
    let pad: i64 = 8;
    let ni = n as i64;
    let contig = contigs.name(r.ref_);
    let (ref_win, ref_in) = match side {
        Side::Right => {
            let (lo, hi) = (j, (j + ni).max(r.pos + ref_len) + pad);
            let win = rf.fetch(contig, lo, hi).to_ascii_uppercase();
            let mut inn = rf.fetch(contig, (j - 6).max(0), j).to_ascii_uppercase();
            inn.reverse();
            (win, inn)
        }
        Side::Left => {
            let (lo, hi) = ((r.pos.min(j - ni) - pad).max(0), j);
            let win = revcomp(&rf.fetch(contig, lo, hi).to_ascii_uppercase());
            let mut inn = revcomp(&rf.fetch(contig, j, j + 6).to_ascii_uppercase());
            inn.reverse();
            (win, inn)
        }
    };
    if ref_win.len() < n {
        return "no_reference";
    }
    let ref_out = &ref_win[..n];
    let c = cons[..(n + pad as usize).min(cons.len())].to_ascii_uppercase();
    let ed_ref = align::distance(&over, &ref_win, Mode::Hw, -1, &[]) as i64;
    let ed_cons = align::distance(&over, &c, Mode::Shw, -1, &[]) as i64;
    let m_cons = ni - ed_cons;
    if m_cons < min_b {
        return "consensus_mismatch";
    }
    if ed_ref < min_mm || ed_cons >= ed_ref {
        return "matches_reference";
    }
    let count = |s: &[u8], b: u8| s.iter().filter(|&&x| x == b).count();
    let mut top = b'A';
    for b in [b'C', b'G', b'T'] {
        if count(&over, b) > count(&over, top) {
            top = b;
        }
    }
    if count(&over, top) as f64 >= 0.8 * n as f64 {
        let head6 = |s: &[u8]| count(&s[..s.len().min(6)], top);
        if head6(&ref_in) >= 4 || head6(ref_out) >= 4 {
            return "ref_homopolymer";
        }
    }
    ""
}

/// `_mate_inside(f, cfg)` (evidence.py:735): the mate of a SHORT fragment lies inside the
/// element (unmapped, other contig, far, or MAPQ < short_mate_min_mapq).
fn mate_inside(f: &Fragment, rows: &[EvidenceRow], cfg: &Config, contigs: &Interner) -> bool {
    let p = &rows[f.primary];
    let Some(mi) = f.mate else { return !p.mate_mapped(contigs) };
    let m = &rows[mi];
    if !m.mapped(contigs) {
        return true;
    }
    if m.ref_ != p.ref_ || (m.pos - p.pos).abs() > cfg.short_mate_max_dist {
        return true;
    }
    m.mapq < cfg.short_mate_min_mapq
}

/// `_clips_start_with_polya(clip_rows, polya_min, within=5)` (evidence.py:749).
fn clips_start_with_polya<'r>(clip_rows: impl Iterator<Item = &'r EvidenceRow>, polya_min: usize) -> bool {
    const WITHIN: usize = 5;
    let (mut n, mut hits) = (0usize, 0usize);
    for r in clip_rows {
        n += 1;
        let mut head = r.outward_clip_seq();
        head.truncate(WITHIN + polya_min);
        head.make_ascii_uppercase();
        let hit = [b'A', b'T'].iter().any(|&c| find_run(&head, c, polya_min).is_some_and(|i| i <= WITHIN));
        hits += hit as usize;
    }
    if n == 0 {
        return false;
    }
    2 * hits >= n
}

/// `s.find(c * k)` (python: lowest index, 0 for k == 0).
fn find_run(s: &[u8], c: u8, k: usize) -> Option<usize> {
    if k == 0 {
        return Some(0);
    }
    let mut run = 0usize;
    for (i, &x) in s.iter().enumerate() {
        run = if x == c { run + 1 } else { 0 };
        if run == k {
            return Some(i + 1 - k);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::row::p3_fixture::{self, Fixture};
    use serde_json::Value;
    use rustc_hash::FxHashMap;

    /// python's `ref_fetch` calls, replayed (a fetch python did not make panics).
    struct Replay(FxHashMap<String, Vec<u8>>);
    impl RefFetch for Replay {
        fn fetch(&self, seqname: &str, start: i64, end: i64) -> Vec<u8> {
            let k = format!("{seqname}:{start}:{end}");
            self.0.get(&k).unwrap_or_else(|| panic!("fetch {k} not made by python")).clone()
        }
    }

    fn fnv(h: &mut u64, b: &[u8]) {
        for &x in b {
            *h = (*h ^ x as u64).wrapping_mul(0x100000001b3);
        }
    }

    /// The generator's digest of a `reads` list (seq|quals|group|weight bits|anchored;)*.
    fn digest(reads: &[ClipRead]) -> String {
        let mut h: u64 = 0xcbf29ce484222325;
        for r in reads {
            let q: Vec<String> = r.qual.iter().map(|x| x.to_string()).collect();
            let s = format!(
                "{}|{}|{}|{}|{};",
                String::from_utf8_lossy(&r.seq),
                q.join(","),
                r.group,
                r.weight.to_bits(),
                r.anchored as u8
            );
            fnv(&mut h, s.as_bytes());
        }
        format!("{h:016x}")
    }

    fn cres(v: &Value) -> ConsensusResult {
        let ints = |k: &str| -> Vec<i64> { v[k].as_array().unwrap().iter().map(|x| x.as_i64().unwrap()).collect() };
        ConsensusResult {
            seq: v["seq"].as_str().unwrap().as_bytes().to_vec(),
            depth: ints("depth").into_iter().map(|x| x as u32).collect(),
            score: ints("score").into_iter().map(|x| x as u8).collect(),
            stop_reason: match v["stop_reason"].as_str().unwrap() {
                "empty" => "empty",
                "end" => "end",
                "disagreement" => "disagreement",
                "depth" => "depth",
                o => Box::leak(o.to_string().into_boxed_str()),
            },
            polya_base: v["polya_base"].as_str().map(|s| s.as_bytes()[0]),
            polya_start: v["polya_start"].as_i64().unwrap(),
            polya_len_median: v["polya_len_median"].as_i64(),
            polya_len_min: v["polya_len_min"].as_i64(),
            polya_len_max: v["polya_len_max"].as_i64(),
            beyond_polya: v["beyond_polya"].as_str().unwrap().as_bytes().to_vec(),
            beyond_start: v["beyond_start"].as_i64().unwrap(),
            beyond_polya_support: v["beyond_polya_support"].as_i64().unwrap(),
        }
    }

    /// Evaluate every fixture junction under configs A/B/C and compare with python: the TSV row,
    /// the reads.fa text, the SHORT rejection reasons and the combined consensus. `real` = use
    /// `consensus::indel_aware_consensus` (needs P2); otherwise python's consensus results are
    /// replayed, after checking that the reads handed to the consensus are python's exactly.
    fn check_fixture(fx: &Fixture, real: bool) -> (usize, u32) {
        let fetch = Replay(
            fx.doc["fetch"].as_object().unwrap().iter().map(|(k, v)| (k.clone(), v.as_str().unwrap().as_bytes().to_vec())).collect(),
        );
        let (mut n, mut short_used) = (0usize, 0u32);
        for (j, (locus, side, rows)) in fx.doc["junctions"].as_array().unwrap().iter().zip(&fx.junctions) {
            for c in ["A", "B", "C"] {
                let ctx = format!("{locus} {} cfg {c}", side.as_str());
                let cfg = p3_fixture::config(&fx.doc["cfgs"][c]);
                let e = &j[c];
                let calls = e["calls"].as_array().unwrap();
                let mut k = 0usize;
                let mut replay = |reads: &[ClipRead], min_depth: usize, polya_min: usize| -> ConsensusResult {
                    let call = calls.get(k).unwrap_or_else(|| panic!("{ctx}: extra consensus call"));
                    k += 1;
                    assert_eq!(reads.len() as u64, call["n"].as_u64().unwrap(), "{ctx} call {k}");
                    assert_eq!(digest(reads), call["digest"].as_str().unwrap(), "{ctx} call {k}");
                    assert_eq!(min_depth as u64, call["min_depth"].as_u64().unwrap(), "{ctx} call {k}");
                    assert_eq!(polya_min as u64, call["polya_min_len"].as_u64().unwrap(), "{ctx} call {k}");
                    if real {
                        let got = indel_aware_consensus(reads, min_depth, polya_min);
                        assert_eq!(got, cres(&call["result"]), "{ctx} call {k}: consensus");
                        got
                    } else {
                        cres(&call["result"])
                    }
                };
                let mut rec = evaluate_junction_with(locus, *side, rows.clone(), &cfg, Some(&fetch), &fx.contigs, &fx.files, &mut replay);
                assert_eq!(k, calls.len(), "{ctx}: consensus call count");
                rec.aligned = e["aligned"].as_str().unwrap().as_bytes().to_vec();
                assert_eq!(rec.tsv(), e["tsv"].as_str().unwrap(), "{ctx}");
                assert_eq!(rec.render_reads(locus, &fx.files), e["reads"].as_str().unwrap(), "{ctx}");
                let sr: Vec<(String, u32)> = rec.short_reasons.iter().map(|(a, b)| (a.to_string(), *b)).collect();
                let esr: Vec<(String, u32)> = serde_json::from_value(e["short_reasons"].clone()).unwrap();
                assert_eq!(sr, esr, "{ctx}");
                assert_eq!(rec.combined_consensus, cres(&e["combined"]), "{ctx}");
                assert_eq!(rec.supported, Supported::Na);
                n += 1;
                short_used += rec.n_short_used;
            }
        }
        (n, short_used)
    }

    #[test]
    fn evaluate_junction_matches_python_replayed_consensus() {
        let (n, used) = check_fixture(&p3_fixture::checked_in(), false);
        assert!(n >= 120 && used > 0, "{n} {used}");
    }

    /// Full fixture, replayed consensus (`P3_FULL_DIR`).
    #[test]
    #[ignore]
    fn evaluate_junction_matches_python_replayed_consensus_full() {
        let Some(fx) = p3_fixture::full() else { return eprintln!("P3_FULL_DIR not set: skipped") };
        let (n, used) = check_fixture(&fx, false);
        eprintln!("full fixture: {n} evaluations, {used} SHORT fragments used");
    }

    /// End to end with the real indel-aware consensus (P2); the full fixture too when
    /// `P3_FULL_DIR` is set.
    #[test]
    fn evaluate_junction_matches_python_end_to_end() {
        check_fixture(&p3_fixture::checked_in(), true);
        if let Some(fx) = p3_fixture::full() {
            check_fixture(&fx, true);
        }
    }

    #[test]
    fn short_overhang_reasons_without_reference_or_clip() {
        let fx = p3_fixture::checked_in();
        let cfg = p3_fixture::config(&fx.doc["cfgs"]["B"]);
        // any SHORT row of the fixture
        let r = fx.junctions.iter().flat_map(|j| j.2.iter()).find(|r| r.role == Role::Short && r.clip_at > 0).unwrap();
        assert_eq!(short_overhang_check(r, r.side, b"ACGTACGTAC", false, &cfg, None, &fx.contigs), "no_clip_fragment");
        assert_eq!(short_overhang_check(r, r.side, b"", true, &cfg, None, &fx.contigs), "no_consensus");
        let long = vec![b'A'; 200];
        let why = short_overhang_check(r, r.side, &long, true, &cfg, None, &fx.contigs);
        assert!(why == "no_reference" || why == "overhang_too_short" || why == "bad_record", "{why}");
    }

    /// An in-memory reference (coordinates only; the contig name is ignored).
    struct Mem(Vec<u8>);
    impl RefFetch for Mem {
        fn fetch(&self, _seqname: &str, start: i64, end: i64) -> Vec<u8> {
            let (s, e) = (start.max(0) as usize, (end.max(0) as usize).min(self.0.len()));
            if s >= e { Vec::new() } else { self.0[s..e].to_vec() }
        }
    }

    /// Step 4 of `_short_overhang_check`: a homopolymer overhang continuing a reference
    /// homopolymer at the junction is rejected; the same overhang next to a non-A reference passes.
    /// Expected values from python `_short_overhang_check` with the same inputs.
    #[test]
    fn short_overhang_ref_homopolymer() {
        let fx = p3_fixture::checked_in();
        let cfg = p3_fixture::config(&serde_json::json!({}));
        let mut r = fx.junctions.iter().flat_map(|j| j.2.iter()).find(|r| r.role == Role::Short && r.mapped(&fx.contigs)).unwrap().clone();
        r.pos = 0;
        r.cigar = "20M8S".into();
        r.clip_at = 20;
        r.seq = b"CGTGCGTGCGTGCGTGCGTGAAAAAAAA".to_vec().into_boxed_slice();
        let cons = b"AAAAAAAAAA";
        let tail = b"GCGTCCGTGCGTCCGTGCGTCCGT";
        // RIGHT side: ref_in = reverse(ref[14:20]); an A-tract there -> ref_homopolymer
        let mut a = b"CGTGCGTGCGTGCG".to_vec();
        a.extend_from_slice(b"AAAAAA");
        a.extend_from_slice(tail);
        assert_eq!(short_overhang_check(&r, Side::Right, cons, true, &cfg, Some(&Mem(a)), &fx.contigs), "ref_homopolymer");
        let mut b = b"CGTGCGTGCGTGCGTGCGTG".to_vec();
        b.extend_from_slice(tail);
        assert_eq!(short_overhang_check(&r, Side::Right, cons, true, &cfg, Some(&Mem(b)), &fx.contigs), "");
    }

    #[test]
    fn tsv_empty_record() {
        let rec = JunctionRecord {
            insertion_id: "chr1:5-9".into(),
            side: Side::Left,
            rows: Some(Vec::new()),
            reads_ref: None,
            n_reads: 0,
            n_fragments: 0,
            n_independent: 0,
            n_samples: 0,
            n_mates: 0,
            n_duplicates: 0,
            n_cross: 0,
            n_dup_coord: 0,
            n_dup_seq: 0,
            member_loci: String::new(),
            n_short_used: 0,
            n_short_rejected: 0,
            n_short_mate_inside: 0,
            n_independent_no_short: 0,
            short_reasons: Vec::new(),
            supported: Supported::Yes,
            consensus: ConsensusResult::default(),
            combined_consensus: ConsensusResult::default(),
            polya_end: false,
            fail_reason: "x".into(),
            aligned: b"ACG".to_vec(),
        };
        // python: JunctionRecord("chr1:5-9", "LEFT") with supported=1, fail_reason="x", aligned="ACG"
        assert_eq!(rec.tsv(), "chr1:5-9\tLEFT\t0\t0\t0\t0\t0\t1\tACG\t\t\t\t\t0\t0\t0\t0\tx\tempty\t0\t0\t.\t0\t0\t0\t0\n");
        assert_eq!(rec.render_reads("n", &[]), "");
    }
}
