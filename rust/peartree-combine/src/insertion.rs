//! Discovery-file parser and the `+=` merge. OWNER: P1.
//!
//! Mirrors src/combine_insertions_insertion.py (`Insertion.__init__`, `parseFile`,
//! `__iadd__`) and the per-file import loop of combine_insertions.py:116-129. SPEC.md §2.1, §3.1.

use crate::model::{FileId, InsType, Insertion, Interner, LocusKey, Side, Tok, TokKind};
use crate::seq::QualSeq;
use flate2::bufread::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Per-file import result (combine_insertions.py:119-129).
pub struct FileImport {
    pub file: FileId,
    /// records after the contig filter (`len(contig) < 6 and contig not in ("MT", "chrM")`),
    /// in file order. `uid` is unset (0) here.
    pub records: Vec<Insertion>,
}

/// `Insertion.parseFile(path)` + the contig filter, streaming the gzip FASTQ.
///
/// Format: 4-line FASTQ records titled `@<contig>:<start>-<end>:<SIDE>:<FIELD>` (SIDE =
/// LEFT/RIGHT; FIELD = CLIPPED / ALIGNED / CLIPPED_POLYA / MATE<n> / other). Blank lines
/// (after `strip()`) between records are skipped; seq/plus/qual lines are `strip()`ped; qual =
/// `ord(c) - 33`. Consecutive records with the same (contig, start, end) form one Insertion;
/// a repeated non-consecutive id yields a second Insertion (python behaviour). A repeated FIELD
/// within one id overwrites (python dict). MATE records are validated and DISCARDED (dead data,
/// SPEC.md "Dead data"). Type/side decoding exactly as `Insertion.__init__` (SPEC.md §2.1):
/// RIGHT real iff RIGHT:ALIGNED and RIGHT:CLIPPED present (`right_aligned = ALIGNED.revcomp()`,
/// `right_pos = int(end)`), elif end token `oneside_` (type 4, open_side RIGHT), elif `disc_`
/// (type 4), else poly-A (`right_clipped = RIGHT:CLIPPED_POLYA`, type 1); LEFT analogous
/// (`left_clipped = LEFT:CLIPPED.revcomp()`, poly-A: `LEFT:CLIPPED_POLYA.revcomp().lower()`,
/// type 2 / 5); LEFT decoding runs after RIGHT and overwrites `type` (python order).
/// `files = [file]`, `member_loci = [(file, own locus)]`.
///
/// Errors (python would raise): malformed title, missing '+', a poly-A end without its
/// CLIPPED_POLYA record, unparsable coordinate.
///
/// Memory: no per-record String; sequences are boxed slices. Called in parallel (one rayon task
/// per file) by the driver; `contigs` is the shared interner.
pub fn parse_discovery_file(path: &Path, file: FileId, contigs: &Interner) -> Result<FileImport, String> {
    let ctx = |m: String| format!("{}: {m}", path.display());
    let f = File::open(path).map_err(|e| ctx(format!("cannot open: {e}")))?;
    let mut rd = BufReader::with_capacity(1 << 20, MultiGzDecoder::new(BufReader::with_capacity(1 << 16, f)));
    let mut records: Vec<Insertion> = Vec::new();

    // current group
    let mut cur_key: Vec<u8> = Vec::new(); // raw `contig:start-end`
    let mut have_cur = false;
    let mut cur_keep = false; // contig passes the filter
    let mut slots: [Option<QualSeq>; 6] = Default::default();

    let mut line: Vec<u8> = Vec::with_capacity(256);
    let mut seq_buf: Vec<u8> = Vec::with_capacity(256);
    let mut qual_buf: Vec<u8> = Vec::with_capacity(256);
    let mut lineno: u64 = 0;

    macro_rules! flush {
        () => {
            if have_cur && cur_keep {
                let ins = build_insertion(&cur_key, &mut slots, file, contigs).map_err(|e| ctx(e))?;
                records.push(ins);
            }
        };
    }

    loop {
        line.clear();
        let n = rd.read_until(b'\n', &mut line).map_err(|e| ctx(format!("read error: {e}")))?;
        if n == 0 {
            break;
        }
        lineno += 1;
        let title = trim(&line);
        if title.is_empty() {
            continue;
        }
        if title[0] != b'@' {
            return Err(ctx(format!("line {lineno}: expected label field, got {}", String::from_utf8_lossy(title))));
        }
        // `l[1:].split(":")` must give exactly 4 parts; positions.split("-") exactly 2
        let body = &title[1..];
        let mut parts = body.split(|&b| b == b':');
        let (c, pos, side, field) = match (parts.next(), parts.next(), parts.next(), parts.next(), parts.next()) {
            (Some(c), Some(p), Some(s), Some(fl), None) => (c, p, s, fl),
            _ => return Err(ctx(format!("line {lineno}: malformed title {}", String::from_utf8_lossy(title)))),
        };
        if pos.iter().filter(|&&b| b == b'-').count() != 1 {
            return Err(ctx(format!("line {lineno}: malformed positions in {}", String::from_utf8_lossy(title))));
        }
        // group key = title up to the end of the positions field
        let key_len = c.len() + 1 + pos.len();
        let key = &body[..key_len];
        if !have_cur || cur_key != key {
            flush!();
            cur_key.clear();
            cur_key.extend_from_slice(key);
            have_cur = true;
            cur_keep = contig_passes(c);
            slots = Default::default();
        }

        // the three payload lines
        let mut next_line = |rd: &mut BufReader<MultiGzDecoder<BufReader<File>>>, buf: &mut Vec<u8>| -> Result<(), String> {
            buf.clear();
            rd.read_until(b'\n', buf).map_err(|e| format!("read error: {e}"))?;
            Ok(())
        };
        let mut tmp: Vec<u8> = std::mem::take(&mut seq_buf);
        next_line(&mut rd, &mut tmp).map_err(&ctx)?;
        seq_buf = tmp;
        let mut plus: Vec<u8> = std::mem::take(&mut qual_buf);
        next_line(&mut rd, &mut plus).map_err(&ctx)?;
        if trim(&plus) != b"+" {
            return Err(ctx(format!("line {lineno}: expected + separator, got {}", String::from_utf8_lossy(trim(&plus)))));
        }
        next_line(&mut rd, &mut plus).map_err(&ctx)?;
        qual_buf = plus;
        lineno += 3;
        let seq = trim(&seq_buf);
        let qual = trim(&qual_buf);
        if seq.len() != qual.len() {
            return Err(ctx(format!("line {lineno}: sequence/quality length mismatch")));
        }
        if !cur_keep || contains(field, b"MATE") {
            continue; // mates are dead data; filtered contigs are never built
        }
        let slot = match (side, field) {
            (b"RIGHT", b"ALIGNED") => 0,
            (b"RIGHT", b"CLIPPED") => 1,
            (b"RIGHT", b"CLIPPED_POLYA") => 2,
            (b"LEFT", b"ALIGNED") => 3,
            (b"LEFT", b"CLIPPED") => 4,
            (b"LEFT", b"CLIPPED_POLYA") => 5,
            _ => continue,
        };
        let mut q = Vec::with_capacity(qual.len());
        for &c in qual {
            q.push(c.checked_sub(33).ok_or_else(|| ctx(format!("line {lineno}: quality char below '!'")))?);
        }
        slots[slot] = Some(QualSeq::new(seq.to_vec(), q));
    }
    flush!();
    records.shrink_to_fit();
    Ok(FileImport { file, records })
}

fn trim(l: &[u8]) -> &[u8] {
    let mut a = 0;
    let mut b = l.len();
    while a < b && l[a].is_ascii_whitespace() {
        a += 1;
    }
    while b > a && l[b - 1].is_ascii_whitespace() {
        b -= 1;
    }
    &l[a..b]
}

fn contains(h: &[u8], n: &[u8]) -> bool {
    h.windows(n.len()).any(|w| w == n)
}

/// combine_insertions.py:121 `len(i.reference_name) < 6 and i.reference_name not in ("MT", "chrM")`
fn contig_passes(c: &[u8]) -> bool {
    c.len() < 6 && c != b"MT" && c != b"chrM"
}

/// `Insertion.__init__` from one record group. slots: RIGHT ALIGNED/CLIPPED/CLIPPED_POLYA,
/// LEFT ALIGNED/CLIPPED/CLIPPED_POLYA.
fn build_insertion(key: &[u8], slots: &mut [Option<QualSeq>; 6], file: FileId, contigs: &Interner) -> Result<Insertion, String> {
    let key_s = String::from_utf8_lossy(key);
    let (contig, positions) = key_s.split_once(':').ok_or("bad group key")?;
    if contig.is_empty() {
        return Err("reference_name not given".into());
    }
    let (start_s, end_s) = positions.split_once('-').ok_or("bad group key")?;
    let start = Tok::parse(start_s.as_bytes()).ok_or_else(|| format!("unparsable coordinate in {key_s}"))?;
    let end = Tok::parse(end_s.as_bytes()).ok_or_else(|| format!("unparsable coordinate in {key_s}"))?;
    let [r_al, r_cl, r_poly, l_al, l_cl, l_poly] = std::mem::take(slots);

    let mut ty: Option<InsType> = None;
    let mut open_side = None;
    let (right_clipped, right_aligned, right_pos);
    if r_al.is_some() && r_cl.is_some() {
        if end.kind != TokKind::Pos {
            return Err(format!("{key_s}: RIGHT:ALIGNED/CLIPPED present with a non-numeric end token"));
        }
        right_aligned = Some(r_al.unwrap().revcomp());
        right_clipped = r_cl;
        right_pos = Some(end.pos);
    } else if end.kind == TokKind::OneSide {
        right_clipped = None;
        right_aligned = None;
        right_pos = Some(end.pos);
        ty = Some(InsType::RightDisc);
        open_side = Some(Side::Right);
    } else if end.kind == TokKind::Disc {
        right_clipped = None;
        right_aligned = None;
        right_pos = Some(end.pos);
        ty = Some(InsType::RightDisc);
    } else if end.kind == TokKind::PolyA {
        right_clipped = Some(r_poly.ok_or_else(|| format!("{key_s}: poly-A end without RIGHT:CLIPPED_POLYA"))?);
        right_aligned = None;
        right_pos = None;
        ty = Some(InsType::RightPolyA);
    } else {
        return Err(format!("false right poly A detected in {key_s}"));
    }
    let (left_clipped, left_aligned, left_pos);
    if l_al.is_some() && l_cl.is_some() {
        if start.kind != TokKind::Pos {
            return Err(format!("{key_s}: LEFT:ALIGNED/CLIPPED present with a non-numeric start token"));
        }
        left_clipped = Some(l_cl.unwrap().revcomp());
        left_aligned = l_al;
        left_pos = Some(start.pos);
    } else if start.kind == TokKind::OneSide {
        left_clipped = None;
        left_aligned = None;
        left_pos = Some(start.pos);
        ty = Some(InsType::LeftDisc);
        open_side = Some(Side::Left);
    } else if start.kind == TokKind::Disc {
        left_clipped = None;
        left_aligned = None;
        left_pos = Some(start.pos);
        ty = Some(InsType::LeftDisc);
    } else if start.kind == TokKind::PolyA {
        left_clipped = Some(l_poly.ok_or_else(|| format!("{key_s}: poly-A start without LEFT:CLIPPED_POLYA"))?.revcomp().lower());
        left_aligned = None;
        left_pos = None;
        ty = Some(InsType::LeftPolyA);
    } else {
        return Err(format!("false left poly A detected in {key_s}"));
    }
    let contig_id = contigs.intern(contig);
    Ok(Insertion {
        uid: 0,
        contig: contig_id,
        name_start: start,
        name_end: end,
        ty: ty.unwrap_or(InsType::FullInfo),
        open_side,
        left_clipped,
        left_aligned,
        left_pos,
        right_clipped,
        right_aligned,
        right_pos,
        files: vec![file],
        member_loci: vec![(file, LocusKey { contig: contig_id, start, end })],
        member_sides: None,
    })
}

/// `Insertion.__iadd__` (combine_insertions_insertion.py:130). Only the branch reachable from
/// intersect_insertions is exercised (both FULL_INFO -- poly-A merging is dead code since the
/// poly-A loop `continue`s), but port the whole method: for each side, longer clipped replaces,
/// longer aligned replaces (strictly longer); mates are not stored; then
/// `files += other.files`, `member_loci += other.member_loci`. Name/positions unchanged in the
/// reachable branch.
pub fn merge_into(target: &mut Insertion, other: Insertion) {
    // helper: python `if len(a) < len(b): a = b` (None sides cannot occur on reachable paths)
    fn longer(a: &mut Option<QualSeq>, b: Option<QualSeq>) {
        if let (Some(x), Some(y)) = (a.as_ref(), b.as_ref()) {
            if x.len() < y.len() {
                *a = b;
            }
        }
    }
    let Insertion {
        ty: o_ty,
        left_clipped: o_lc,
        left_aligned: o_la,
        left_pos: o_lp,
        right_clipped: o_rc,
        right_aligned: o_ra,
        right_pos: o_rp,
        files: o_files,
        member_loci: o_members,
        ..
    } = other;

    // ---- left
    if target.ty == InsType::LeftPolyA && o_ty != InsType::LeftPolyA {
        target.left_clipped = o_lc;
        target.left_aligned = o_la;
        target.left_pos = o_lp;
        // python: name = f"{ref}:{self.right_pos}-{self.left_pos}"
        if let (Some(r), Some(l)) = (target.right_pos, target.left_pos) {
            target.name_start = Tok::pos(r);
            target.name_end = Tok::pos(l);
        }
        target.ty = InsType::FullInfo;
    } else if target.ty != InsType::LeftPolyA && o_ty == InsType::LeftPolyA {
        // python only appends a mate here (mates are not stored)
    } else {
        longer(&mut target.left_clipped, o_lc);
        if target.ty != InsType::LeftPolyA {
            longer(&mut target.left_aligned, o_la);
        }
    }
    // ---- right
    if target.ty == InsType::RightPolyA && o_ty != InsType::RightPolyA {
        target.right_clipped = o_rc;
        target.right_aligned = o_ra;
        target.right_pos = o_rp;
        if let (Some(r), Some(l)) = (target.right_pos, target.left_pos) {
            target.name_start = Tok::pos(r);
            target.name_end = Tok::pos(l);
        }
        target.ty = InsType::FullInfo;
    } else if target.ty != InsType::RightPolyA && o_ty == InsType::RightPolyA {
    } else {
        longer(&mut target.right_clipped, o_rc);
        if target.ty != InsType::RightPolyA {
            longer(&mut target.right_aligned, o_ra);
        }
    }
    target.files.extend(o_files);
    target.member_loci.extend(o_members);
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{write::GzEncoder, Compression};
    use std::io::Write;

    fn rec(title: &str, seq: &str) -> String {
        format!("@{title}\n{seq}\n+\n{}\n", "I".repeat(seq.len()))
    }

    fn write_gz(name: &str, text: &str) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!("pt_combine_p1_{}_{name}", std::process::id()));
        let mut e = GzEncoder::new(File::create(&p).unwrap(), Compression::fast());
        e.write_all(text.as_bytes()).unwrap();
        e.finish().unwrap();
        p
    }

    #[test]
    fn parses_types_mates_and_filters() {
        let mut t = String::new();
        // full record (+ a MATE that must vanish, a blank line between records, a repeated field)
        t += &rec("chr1:100-120:LEFT:CLIPPED", "AACCGG");
        t += &rec("chr1:100-120:LEFT:ALIGNED", "TTTTGGGA");
        t += "\n";
        t += &rec("chr1:100-120:LEFT:MATE1", "ACGT");
        t += &rec("chr1:100-120:RIGHT:CLIPPED", "CCA");
        t += &rec("chr1:100-120:RIGHT:ALIGNED", "ACGTAC");
        t += &rec("chr1:100-120:RIGHT:ALIGNED", "ACGTACG"); // overwrites
        // oneside right
        t += &rec("chr2:200-oneside_250:LEFT:CLIPPED", "GGGAAA");
        t += &rec("chr2:200-oneside_250:LEFT:ALIGNED", "ACACAC");
        // poly-A right (legacy)
        t += &rec("chr3:300-polyA_333:LEFT:CLIPPED", "CCCCAA");
        t += &rec("chr3:300-polyA_333:LEFT:ALIGNED", "ACACAC");
        t += &rec("chr3:300-polyA_333:RIGHT:CLIPPED_POLYA", "AAAAAAAAAAAA");
        // filtered contigs
        t += &rec("chrM:5-9:LEFT:CLIPPED", "AAAAAA");
        t += &rec("chrUn_KI270742v1:5-9:LEFT:CLIPPED", "AAAAAA");
        let p = write_gz("a.txt.gz", &t);
        let contigs = Interner::new();
        let imp = parse_discovery_file(&p, 7, &contigs).unwrap();
        std::fs::remove_file(&p).ok();
        assert_eq!(imp.records.len(), 3);
        let a = &imp.records[0];
        assert_eq!(a.ty, InsType::FullInfo);
        assert_eq!((a.left_pos, a.right_pos), (Some(100), Some(120)));
        assert_eq!(a.name(&contigs), "chr1:100-120");
        // left_clipped = CLIPPED.revcomp(); right_aligned = ALIGNED.revcomp() of the last one
        assert_eq!(&a.left_clipped.as_ref().unwrap().seq[..], b"CCGGTT");
        assert_eq!(&a.right_aligned.as_ref().unwrap().seq[..], b"CGTACGT");
        assert_eq!(a.files, vec![7]);
        assert_eq!(a.member_loci.len(), 1);
        let b = &imp.records[1];
        assert_eq!(b.ty, InsType::RightDisc);
        assert_eq!(b.open_side, Some(Side::Right));
        assert_eq!(b.right_pos, Some(250));
        assert!(b.right_clipped.is_none() && b.right_aligned.is_none());
        assert_eq!(b.name(&contigs), "chr2:200-oneside_250");
        let c = &imp.records[2];
        assert_eq!(c.ty, InsType::RightPolyA);
        assert_eq!(c.right_pos, None);
        assert_eq!(&c.right_clipped.as_ref().unwrap().seq[..], b"AAAAAAAAAAAA");
        assert!(contigs.get("chrM").is_none() && contigs.get("chrUn_KI270742v1").is_none());
    }

    #[test]
    fn rejects_malformed() {
        let p = write_gz("b.txt.gz", "@chr1:1-2:LEFT\nAC\n+\nII\n");
        assert!(parse_discovery_file(&p, 0, &Interner::new()).is_err());
        let p2 = write_gz("c.txt.gz", "@chr1:1-2:LEFT:CLIPPED\nAC\nX\nII\n");
        assert!(parse_discovery_file(&p2, 0, &Interner::new()).is_err());
        std::fs::remove_file(&p).ok();
        std::fs::remove_file(&p2).ok();
    }

    #[test]
    fn merge_keeps_longer_and_concatenates() {
        let mk = |lc: &str, la: &str, f: FileId| Insertion {
            uid: 0,
            contig: 0,
            name_start: Tok::pos(1),
            name_end: Tok::pos(2),
            ty: InsType::FullInfo,
            open_side: None,
            left_clipped: Some(QualSeq::new(lc.as_bytes().to_vec(), vec![1; lc.len()])),
            left_aligned: Some(QualSeq::new(la.as_bytes().to_vec(), vec![1; la.len()])),
            left_pos: Some(1),
            right_clipped: Some(QualSeq::new(b"AC".to_vec(), vec![1; 2])),
            right_aligned: Some(QualSeq::new(b"AC".to_vec(), vec![1; 2])),
            right_pos: Some(2),
            files: vec![f],
            member_loci: vec![(f, LocusKey { contig: 0, start: Tok::pos(1), end: Tok::pos(2) })],
            member_sides: None,
        };
        let mut a = mk("AAA", "CCCC", 0);
        merge_into(&mut a, mk("GGGG", "TT", 1));
        assert_eq!(&a.left_clipped.as_ref().unwrap().seq[..], b"GGGG");
        assert_eq!(&a.left_aligned.as_ref().unwrap().seq[..], b"CCCC");
        merge_into(&mut a, mk("TTTT", "GGGG", 1)); // equal clip length: not replaced; equal aligned: not replaced
        assert_eq!(&a.left_clipped.as_ref().unwrap().seq[..], b"GGGG");
        assert_eq!(&a.left_aligned.as_ref().unwrap().seq[..], b"CCCC");
        assert_eq!(a.files, vec![0, 1, 1]);
        assert_eq!(a.member_loci.len(), 3);
    }
}
