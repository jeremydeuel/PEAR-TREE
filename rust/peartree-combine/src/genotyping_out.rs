//! `<stem>.genotyping.txt.gz` (the genotyping contract). OWNER: P6.
//!
//! Mirrors combine_insertions.py:298-365. SPEC.md §6.3.

use crate::context::Ctx;
use crate::genome::RefFetch;
use crate::model::{InsType, Insertion, Interner};
use crate::seq::{revcomp, sequence_matching_score, ScoreSeq};
use rustc_hash::FxHashSet;

/// Returns the uncompressed file text (and logs `wrote N of M, excluded K`).
/// `filter_names` = the clipped-remap filter set (python `filter_reads`, re-checked here).
pub fn genotyping_text(insertions: &[Insertion], filter_names: &FxHashSet<String>, ctx: &Ctx) -> Vec<u8> {
    genotyping_text_with(insertions, filter_names, &ctx.contigs, &ctx.genome, ctx.cfg.genotyping_max_bases, &sequence_matching_score)
}

/// [`genotyping_text`] with its dependencies spelled out (reference, contig names, `max_bases`,
/// `sequence_matching_score`) so it can be driven without a [`Ctx`].
pub fn genotyping_text_with(
    insertions: &[Insertion],
    filter_names: &FxHashSet<String>,
    contigs: &Interner,
    reference: &dyn RefFetch,
    max_bases: usize,
    sms: &dyn Fn(&[ScoreSeq]) -> f64,
) -> Vec<u8> {
    let plain_score = |a: &[u8], b: &[u8]| sms(&[ScoreSeq::Plain(a), ScoreSeq::Plain(b)]);
    let mut out: Vec<u8> = Vec::new();
    let (mut n_excluded, mut n_included) = (0usize, 0usize);
    for ins in insertions {
        let name = ins.name(contigs);
        if filter_names.contains(&name) {
            continue;
        }
        // Feature A: discordant-anchored calls carry one real side and are not genotyped
        if ins.ty == InsType::LeftDisc || ins.ty == InsType::RightDisc {
            n_excluded += 1;
            continue;
        }
        // right = right_clipped.upper()[:max]; left = left_clipped[:max].upper().revcomp()
        // (upper/slice/revcomp commute on the base alphabet, so work on the raw bytes)
        let rc = ins.right_clipped.as_ref().expect("genotyping: insertion without right_clipped");
        let lc = ins.left_clipped.as_ref().expect("genotyping: insertion without left_clipped");
        let right: Vec<u8> = rc.seq.iter().take(max_bases).map(|c| c.to_ascii_uppercase()).collect();
        let left_fwd: Vec<u8> = lc.seq.iter().take(max_bases).map(|c| c.to_ascii_uppercase()).collect();
        let left = revcomp(&left_fwd);
        let contig = contigs.name(ins.contig);
        if contig.len() > 5 {
            println!("Excluding {name} since the breakpoint is on {contig}.");
            n_excluded += 1;
            continue;
        }
        let right_ref: Option<Vec<u8>> = if ins.ty != InsType::RightPolyA {
            let p = ins.right_pos.expect("genotyping: insertion without right_pos");
            Some(upper(reference.fetch(contig, p, p + max_bases as i64)))
        } else {
            None
        };
        let left_ref: Option<Vec<u8>> = if ins.ty != InsType::LeftPolyA {
            let p = ins.left_pos.expect("genotyping: insertion without left_pos");
            Some(upper(reference.fetch(contig, p - max_bases as i64, p)))
        } else {
            None
        };
        if ins.ty != InsType::RightPolyA && right_ref.as_ref().is_some_and(|r| r.is_empty()) {
            println!("Excluding {name} due to missing reference");
            n_excluded += 1;
            continue;
        }
        if ins.ty != InsType::LeftPolyA && left_ref.as_ref().is_some_and(|r| r.is_empty()) {
            println!("Excluding {name} due to missing reference");
            n_excluded += 1;
            continue;
        }
        if right_ref.as_ref().is_some_and(|r| r.contains(&b'N')) {
            println!("Excluding {name} due to Ns in right reference");
            n_excluded += 1;
            continue;
        }
        if left_ref.as_ref().is_some_and(|r| r.contains(&b'N')) {
            println!("Excluding {name} due to Ns in left reference");
            n_excluded += 1;
            continue;
        }
        if ins.ty == InsType::RightPolyA && plain_score(left_ref.as_deref().unwrap(), &left) > 0.0 {
            println!("Excluding {name} due to similar left sequence between ref and alt");
            n_excluded += 1;
            continue;
        }
        if ins.ty == InsType::LeftPolyA && plain_score(right_ref.as_deref().unwrap(), &right) > 0.0 {
            println!("Excluding {name} due to similar right sequence between ref and alt");
            n_excluded += 1;
            continue;
        }
        if ins.ty == InsType::FullInfo
            && plain_score(right_ref.as_deref().unwrap(), &right) > 0.0
            && plain_score(left_ref.as_deref().unwrap(), &left) > 0.0
        {
            println!("Excluding {name} due to similar right and left sequence between ref and alt");
            n_excluded += 1;
            continue;
        }
        if ins.ty == InsType::FullInfo {
            let (l, r) = (left_ref.as_ref().unwrap(), right_ref.as_ref().unwrap());
            if *l == revcomp(r) || l == r {
                println!("Excluding {name} due to identical clipped sequences");
                n_excluded += 1;
                continue;
            }
        }
        n_included += 1;
        out.extend_from_slice(b">");
        out.extend_from_slice(name.as_bytes());
        out.push(b'\n');
        if ins.ty != InsType::RightPolyA {
            out.extend_from_slice(b"@RIGHT_INSERTION\n");
            out.extend_from_slice(&right);
            out.extend_from_slice(b"\n@RIGHT_REFERENCE\n");
            out.extend_from_slice(right_ref.as_ref().unwrap());
            out.push(b'\n');
        }
        if ins.ty != InsType::LeftPolyA {
            out.extend_from_slice(b"@LEFT_INSERTION\n");
            out.extend_from_slice(&left);
            out.extend_from_slice(b"\n@LEFT_REFERENCE\n");
            out.extend_from_slice(left_ref.as_ref().unwrap());
            out.push(b'\n');
        }
    }
    println!("wrote {n_included} of {}, excluded {n_excluded}", n_included + n_excluded);
    out
}

fn upper(mut v: Vec<u8>) -> Vec<u8> {
    v.make_ascii_uppercase();
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Tok, TokKind};
    use crate::seq::QualSeq;

    /// Reference with py2bit/get_sequence semantics (end clamped; start < 0 or start >= end -> "").
    struct MapRef(std::collections::HashMap<String, Vec<u8>>);
    impl RefFetch for MapRef {
        fn fetch(&self, seqname: &str, start: i64, end: i64) -> Vec<u8> {
            let Some(s) = self.0.get(seqname) else { return vec![] };
            if start >= end {
                return vec![];
            }
            let end = end.min(s.len() as i64);
            if start < 0 || start >= end {
                return vec![];
            }
            s[start as usize..end as usize].to_vec()
        }
    }

    /// Output vs the python loop (combine_insertions.py:298-365 executed verbatim on the same
    /// synthetic insertions; generator `tests/data/genotyping_cases_gen.py`): 399 insertions
    /// covering every exclusion branch.
    #[test]
    fn genotyping_matches_python() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/genotyping_cases.json");
        let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
        let max_bases = v["max_bases"].as_u64().unwrap() as usize;
        let reference = MapRef(v["genome"].as_object().unwrap().iter().map(|(k, s)| (k.clone(), s.as_str().unwrap().as_bytes().to_vec())).collect());
        let contigs = Interner::new();
        let qs = |s: &str| QualSeq::new(s.as_bytes().to_vec(), vec![30; s.len()]);
        let mut ins = Vec::new();
        for (uid, x) in v["insertions"].as_array().unwrap().iter().enumerate() {
            let name = x["name"].as_str().unwrap();
            let contig = x["contig"].as_str().unwrap();
            let (a, b) = name.rsplit_once(':').unwrap().1.split_once('-').unwrap();
            let ty = match x["type"].as_i64().unwrap() {
                1 => InsType::RightPolyA,
                2 => InsType::LeftPolyA,
                3 => InsType::FullInfo,
                4 => InsType::RightDisc,
                _ => InsType::LeftDisc,
            };
            ins.push(Insertion {
                uid: uid as u32,
                contig: contigs.intern(contig),
                name_start: Tok { kind: TokKind::Pos, pos: a.parse().unwrap() },
                name_end: Tok { kind: TokKind::Pos, pos: b.parse().unwrap() },
                ty,
                open_side: None,
                left_clipped: Some(qs(x["left_clipped"].as_str().unwrap())),
                left_aligned: None,
                left_pos: x["left_pos"].as_i64(),
                right_clipped: Some(qs(x["right_clipped"].as_str().unwrap())),
                right_aligned: None,
                right_pos: x["right_pos"].as_i64(),
                files: vec![0],
                member_loci: vec![],
                member_sides: None,
            });
        }
        let filter: FxHashSet<String> = v["filter"].as_array().unwrap().iter().map(|s| s.as_str().unwrap().to_string()).collect();
        let got = genotyping_text_with(&ins, &filter, &contigs, &reference, max_bases, &sequence_matching_score);
        let want = v["expected"].as_str().unwrap();
        assert_eq!(String::from_utf8(got).unwrap(), want);
    }
}
