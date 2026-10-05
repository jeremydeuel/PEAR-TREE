//! Evidence sidecar rows. OWNER: P3.
//!
//! Mirrors `EvidenceRow`, `_int`, `_clip_at_from_cigar`, `_ref_pos_at`, `_junction_pos` cigar
//! helpers and `allele_forward_seq` of src/combine_insertions_evidence.py (lines 80-148, 629-648,
//! 1464-1480). SPEC.md §4.1.

use crate::model::{ContigId, FileId, Interner, LocusKey, Side};
use crate::seq::revcomp;

/// Sidecar `role` column. Unknown roles are kept verbatim (python treats them as priority 5,
/// not evidence, printed as-is in reads.fa).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Role {
    Clip,
    PolyA,
    Disc,
    Span,
    Short,
    Mate,
    Other(Box<str>),
}

impl Role {
    pub fn parse(s: &str) -> Role {
        match s {
            "CLIP" => Role::Clip,
            "POLYA" => Role::PolyA,
            "DISC" => Role::Disc,
            "SPAN" => Role::Span,
            "SHORT" => Role::Short,
            "MATE" => Role::Mate,
            o => Role::Other(o.into()),
        }
    }
    pub fn as_str(&self) -> &str {
        match self {
            Role::Clip => "CLIP",
            Role::PolyA => "POLYA",
            Role::Disc => "DISC",
            Role::Span => "SPAN",
            Role::Short => "SHORT",
            Role::Mate => "MATE",
            Role::Other(s) => s,
        }
    }
    /// `_ROLE_PRIORITY.get(role, 5)`: CLIP 0, POLYA 1, DISC 2, SPAN 3, SHORT 4, MATE 9, other 5.
    pub fn priority(&self) -> u8 {
        match self {
            Role::Clip => 0,
            Role::PolyA => 1,
            Role::Disc => 2,
            Role::Span => 3,
            Role::Short => 4,
            Role::Mate => 9,
            Role::Other(_) => 5,
        }
    }
    /// `role in EVIDENCE_ROLES` (CLIP, POLYA, DISC, SPAN, SHORT).
    pub fn is_evidence(&self) -> bool {
        matches!(self, Role::Clip | Role::PolyA | Role::Disc | Role::Span | Role::Short)
    }
}

/// One sidecar row (python `EvidenceRow`). Missing columns take the python defaults
/// (`r12` 0, `flag` 0, `ref` "*", `pos` -1, `strand` "*", `outer` -1, `mref` "*", `mpos` -1,
/// `mstrand` "*", `tlen` 0, `mapq` 0, `cigar` "*", `clip_at` -1, `seq` "", `qual` ""); integer
/// columns that fail python `int()` take the default too.
#[derive(Clone, Debug)]
pub struct EvidenceRow {
    /// the discovery file the sidecar belongs to (python `sample` = sample_name(file))
    pub file: FileId,
    /// `locus` column, parsed (rows whose locus does not parse are never wanted)
    pub locus: LocusKey,
    pub side: Side,
    pub role: Role,
    pub frag: Box<str>,
    pub r12: i64,
    pub flag: i64,
    /// interned `ref` string ("*" / "" allowed)
    pub ref_: ContigId,
    pub pos: i64,
    /// strand column verbatim (normally "+", "-" or "*")
    pub strand: Box<str>,
    pub outer: i64,
    pub mref: ContigId,
    pub mpos: i64,
    pub mstrand: Box<str>,
    pub tlen: i64,
    pub mapq: i64,
    pub cigar: Box<str>,
    pub clip_at: i64,
    pub seq: Box<[u8]>,
    /// raw quality string (phred+33)
    pub qual: Box<[u8]>,
}

/// Column positions of one sidecar header (`header.index(name)`; None when absent).
#[derive(Clone, Debug, Default)]
pub struct SidecarHeader {
    pub n_cols: usize,
    pub locus: usize,
    pub side: Option<usize>,
    pub role: Option<usize>,
    pub frag: Option<usize>,
    pub r12: Option<usize>,
    pub flag: Option<usize>,
    pub ref_: Option<usize>,
    pub pos: Option<usize>,
    pub strand: Option<usize>,
    pub outer: Option<usize>,
    pub mref: Option<usize>,
    pub mpos: Option<usize>,
    pub mstrand: Option<usize>,
    pub tlen: Option<usize>,
    pub mapq: Option<usize>,
    pub cigar: Option<usize>,
    pub clip_at: Option<usize>,
    pub seq: Option<usize>,
    pub qual: Option<usize>,
}

impl SidecarHeader {
    /// Parse the header line (`line.rstrip("\n").split("\t")`). Error when `locus` is missing
    /// (python `header.index("locus")` raises). NB python also requires `side`, `role`, `frag`
    /// (`d["side"]` KeyError on the first row) -- error lazily, on the first row, like python.
    ///
    /// Duplicate column names: python builds `dict(zip(header, p))`, so the LAST column of a name
    /// wins for every field; only the locus *filter* uses `header.index("locus")` (the first).
    /// `locus` here is the first occurrence and is used for the row's locus too (a sidecar with
    /// two `locus` columns does not exist).
    pub fn parse(line: &str) -> Result<SidecarHeader, String> {
        let line = line.trim_end_matches('\n');
        let mut h = SidecarHeader::default();
        let mut locus = None;
        for (i, name) in line.split('\t').enumerate() {
            h.n_cols = i + 1;
            let slot = match name {
                "locus" => {
                    if locus.is_none() {
                        locus = Some(i);
                    }
                    continue;
                }
                "side" => &mut h.side,
                "role" => &mut h.role,
                "frag" => &mut h.frag,
                "r12" => &mut h.r12,
                "flag" => &mut h.flag,
                "ref" => &mut h.ref_,
                "pos" => &mut h.pos,
                "strand" => &mut h.strand,
                "outer" => &mut h.outer,
                "mref" => &mut h.mref,
                "mpos" => &mut h.mpos,
                "mstrand" => &mut h.mstrand,
                "tlen" => &mut h.tlen,
                "mapq" => &mut h.mapq,
                "cigar" => &mut h.cigar,
                "clip_at" => &mut h.clip_at,
                "seq" => &mut h.seq,
                "qual" => &mut h.qual,
                _ => continue,
            };
            *slot = Some(i);
        }
        h.locus = locus.ok_or_else(|| "evidence sidecar header has no 'locus' column".to_string())?;
        Ok(h)
    }
}

/// python `int(x)` on a str: surrounding whitespace, an optional sign, ASCII digits with single
/// underscores between digits. `None` where python raises ValueError (-> the `_int` default).
/// Values beyond i64 (python big ints) also give `None` (never occur in a sidecar).
pub(crate) fn py_int_str(s: &str) -> Option<i64> {
    let s = s.trim_matches(|c: char| c.is_whitespace() || ('\x1c'..='\x1f').contains(&c));
    let b = s.as_bytes();
    let (neg, digits) = match b.first() {
        Some(b'-') => (true, &b[1..]),
        Some(b'+') => (false, &b[1..]),
        _ => (false, b),
    };
    if digits.is_empty() || !digits[0].is_ascii_digit() || !digits[digits.len() - 1].is_ascii_digit() {
        return None;
    }
    let mut v: i64 = 0;
    let mut prev_us = false;
    for &c in digits {
        if c == b'_' {
            if prev_us {
                return None;
            }
            prev_us = true;
            continue;
        }
        if !c.is_ascii_digit() {
            return None;
        }
        prev_us = false;
        v = v.checked_mul(10)?.checked_add((c - b'0') as i64)?;
    }
    Some(if neg { -v } else { v })
}

/// python slice index normalisation for `s[at:]` / `s[:at]` (negative counts from the end, then
/// clamped to `0..=n`).
#[inline]
pub(crate) fn py_idx(at: i64, n: usize) -> usize {
    if at < 0 {
        (n as i64 + at).max(0) as usize
    } else {
        (at as usize).min(n)
    }
}

impl EvidenceRow {
    /// Parse one data line (without the trailing "\n"). Returns Ok(None) when the column count
    /// differs from the header's (python skips the row) or the locus does not parse.
    ///
    /// A `side` other than LEFT/RIGHT also gives Ok(None): python files such a row under a
    /// `(locus, side)` key that is never looked up. A missing `side`/`role`/`frag` column is an
    /// error (python KeyError).
    pub fn parse(line: &str, h: &SidecarHeader, file: FileId, contigs: &Interner) -> Result<Option<EvidenceRow>, String> {
        let line = line.trim_end_matches('\n');
        let mut p: [&str; 64] = [""; 64];
        let mut n = 0usize;
        for f in line.split('\t') {
            if n < p.len() {
                p[n] = f;
            }
            n += 1;
        }
        if n != h.n_cols {
            return Ok(None);
        }
        let big: Vec<&str>;
        let p: &[&str] = if n <= 64 {
            &p[..n]
        } else {
            big = line.split('\t').collect();
            &big
        };
        let get = |i: Option<usize>| i.map(|i| p[i]);
        let locus_s = p[h.locus];
        let side_s = get(h.side).ok_or("evidence sidecar has no 'side' column (python KeyError)")?;
        let role_s = get(h.role).ok_or("evidence sidecar has no 'role' column (python KeyError)")?;
        let frag_s = get(h.frag).ok_or("evidence sidecar has no 'frag' column (python KeyError)")?;
        let Some(locus) = LocusKey::parse(locus_s, contigs) else { return Ok(None) };
        let Some(side) = Side::parse(side_s) else { return Ok(None) };
        let int = |i: Option<usize>, d: i64| get(i).and_then(py_int_str).unwrap_or(d);
        Ok(Some(EvidenceRow {
            file,
            locus,
            side,
            role: Role::parse(role_s),
            frag: frag_s.into(),
            r12: int(h.r12, 0),
            flag: int(h.flag, 0),
            ref_: contigs.intern(get(h.ref_).unwrap_or("*")),
            pos: int(h.pos, -1),
            strand: get(h.strand).unwrap_or("*").into(),
            outer: int(h.outer, -1),
            mref: contigs.intern(get(h.mref).unwrap_or("*")),
            mpos: int(h.mpos, -1),
            mstrand: get(h.mstrand).unwrap_or("*").into(),
            tlen: int(h.tlen, 0),
            mapq: int(h.mapq, 0),
            cigar: get(h.cigar).unwrap_or("*").into(),
            clip_at: int(h.clip_at, -1),
            seq: get(h.seq).unwrap_or("").as_bytes().into(),
            qual: get(h.qual).unwrap_or("").as_bytes().into(),
        }))
    }

    /// `mapped`: ref not in ("*", "") and pos >= 0 and not flag & 0x4.
    pub fn mapped(&self, contigs: &Interner) -> bool {
        let r = contigs.name(self.ref_);
        r != "*" && !r.is_empty() && self.pos >= 0 && self.flag & 0x4 == 0
    }

    /// `mate_mapped`: mref not in ("*", "") and mpos >= 0 and not flag & 0x8.
    pub fn mate_mapped(&self, contigs: &Interner) -> bool {
        let r = contigs.name(self.mref);
        r != "*" && !r.is_empty() && self.mpos >= 0 && self.flag & 0x8 == 0
    }

    /// `quals()`: `ord(c)-33` per char when len(qual) == len(seq), else `[30] * len(seq)`.
    pub fn quals(&self) -> Vec<u8> {
        if self.qual.len() == self.seq.len() {
            self.qual.iter().map(|&c| c.wrapping_sub(33)).collect()
        } else {
            vec![30; self.seq.len()]
        }
    }

    /// The effective clip offset of `outward_clip` (python slice index, normalised).
    fn outward_at(&self) -> usize {
        let n = self.seq.len();
        let mut at = self.clip_at;
        if !(0 <= at && at <= n as i64) {
            at = clip_at_from_cigar(&self.cigar, self.side, n);
        }
        match self.side {
            Side::Right => py_idx(at, n),
            Side::Left => py_idx(at, n),
        }
    }

    /// `outward_clip()[0]` only (no qualities).
    pub fn outward_clip_seq(&self) -> Vec<u8> {
        let at = self.outward_at();
        match self.side {
            Side::Right => self.seq[at..].to_vec(),
            Side::Left => revcomp(&self.seq[..at]),
        }
    }

    /// `outward_clip()`: (seq, quals) of the junction clip, outward. `at = clip_at` unless not
    /// `0 <= at <= len(seq)`, then `_clip_at_from_cigar`. RIGHT: `seq[at:]`, `q[at:]`;
    /// LEFT: `revcomp(seq[:at])`, `q[:at][::-1]`.
    pub fn outward_clip(&self) -> (Vec<u8>, Vec<u8>) {
        let at = self.outward_at();
        let q = self.quals();
        match self.side {
            Side::Right => (self.seq[at..].to_vec(), q[at..].to_vec()),
            Side::Left => (revcomp(&self.seq[..at]), q[..at].iter().rev().copied().collect()),
        }
    }
}

/// `re.findall(r"(\d+)([MIDNSHP=X])", cigar)` as an iterator of (length, op). A digit run not
/// followed by an op letter yields nothing (shorter runs are followed by a digit, so the regex
/// cannot match inside it either). Lengths beyond i64 saturate (never occur).
pub(crate) fn cigar_ops(cigar: &str) -> impl Iterator<Item = (i64, u8)> + '_ {
    let b = cigar.as_bytes();
    let mut i = 0usize;
    std::iter::from_fn(move || {
        while i < b.len() {
            if !b[i].is_ascii_digit() {
                i += 1;
                continue;
            }
            let s = i;
            while i < b.len() && b[i].is_ascii_digit() {
                i += 1;
            }
            if i < b.len() && b"MIDNSHP=X".contains(&b[i]) {
                let mut v: i64 = 0;
                for &c in &b[s..i] {
                    v = v.saturating_mul(10).saturating_add((c - b'0') as i64);
                }
                let op = b[i];
                i += 1;
                return Some((v, op));
            }
        }
        None
    })
}

/// `_clip_at_from_cigar(cigar, side, n)`: ops = regex `(\d+)([MIDNSHP=X])` findall; none ->
/// 0 (LEFT) / n (RIGHT); LEFT: first op S -> its length else 0; RIGHT: last op S -> n - len else n.
pub fn clip_at_from_cigar(cigar: &str, side: Side, n: usize) -> i64 {
    let mut first = None;
    let mut last = None;
    for op in cigar_ops(cigar) {
        if first.is_none() {
            first = Some(op);
        }
        last = Some(op);
    }
    match side {
        Side::Left => match first {
            Some((l, b'S')) => l,
            _ => 0,
        },
        Side::Right => match last {
            Some((l, b'S')) => n as i64 - l,
            _ => n as i64,
        },
    }
}

/// Sum of `(\d+)([MDN=X])` lengths (python `ref_len` in `_junction_pos` / SHORT check).
/// NB the regex is applied to the raw cigar string with findall -- a malformed cigar still
/// yields whatever digit+op pairs it contains.
///
/// (The narrower regex `(\d+)([MDN=X])` finds exactly the `MDN=X` pairs of the wide one: a
/// digit run followed by I/S/H/P matches neither.)
pub fn cigar_ref_len(cigar: &str) -> i64 {
    cigar_ops(cigar).filter(|&(_, op)| b"MDN=X".contains(&op)).map(|(l, _)| l).sum()
}

/// `_ref_pos_at(pos, cigar, q)` (evidence.py:629).
pub fn ref_pos_at(pos: i64, cigar: &str, q: i64) -> i64 {
    let (mut r, mut qi) = (pos, 0i64);
    for (n, op) in cigar_ops(cigar) {
        match op {
            b'M' | b'=' | b'X' => {
                if q < qi + n {
                    return r + (q - qi);
                }
                r += n;
                qi += n;
            }
            b'I' | b'S' => {
                if q < qi + n {
                    return r;
                }
                qi += n;
            }
            b'D' | b'N' => r += n,
            _ => {}
        }
    }
    r
}

/// `allele_forward_seq(r)` (evidence.py:1464): seq unchanged unless seq non-empty and != "*",
/// role is MATE/POLYA and flag & 0x1; then stored if `bool(flag & 0x10) == (not flag & 0x20)`
/// else revcomp.
pub fn allele_forward_seq(r: &EvidenceRow) -> Vec<u8> {
    if needs_revcomp(r) {
        revcomp(&r.seq)
    } else {
        r.seq.to_vec()
    }
}

/// True when `allele_forward_seq(r)` is the reverse complement of `r.seq`.
#[inline]
pub(crate) fn needs_revcomp(r: &EvidenceRow) -> bool {
    if r.seq.is_empty() || &*r.seq == b"*" || !matches!(r.role, Role::Mate | Role::PolyA) || r.flag & 0x1 == 0 {
        return false;
    }
    let stored_rev = r.flag & 0x10 != 0;
    let partner_fwd = r.flag & 0x20 == 0;
    stored_rev != partner_fwd
}

/// Test fixture shared by the P3 unit tests: sidecar rows of the TPRT E2E fixture (3 colonies)
/// pooled per (locus, side) + the python reference values, generated by the P3 generator
/// script (`gen.py subset`, reference venv, src/ at the port's base commit) and checked in under
/// `src/evidence/testdata/`. `P3_FULL_DIR=<dir>` points the `#[ignore]` full-fixture tests at
/// `gen.py full <dir>` output (all 2872 junctions).
#[cfg(test)]
pub(crate) mod p3_fixture {
    use super::*;
    use crate::model::InputFile;
    use serde_json::Value;
    use std::io::Read;

    pub struct Fixture {
        pub doc: Value,
        pub files: Vec<InputFile>,
        pub contigs: Interner,
        /// (locus id, side, pooled rows) in the order of `doc["junctions"]`
        pub junctions: Vec<(String, Side, Vec<EvidenceRow>)>,
    }

    fn gunzip(path: &std::path::Path) -> String {
        let f = std::fs::File::open(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        let mut s = String::new();
        flate2::read::MultiGzDecoder::new(f).read_to_string(&mut s).unwrap();
        s
    }

    pub fn load(dir: &std::path::Path) -> Fixture {
        let doc: Value = serde_json::from_str(&gunzip(&dir.join("expected.json.gz"))).unwrap();
        let rows_txt = gunzip(&dir.join("rows.tsv.gz"));
        let files: Vec<InputFile> =
            doc["files"].as_array().unwrap().iter().map(|f| InputFile::new(f.as_str().unwrap())).collect();
        let contigs = Interner::new();
        let mut lines = rows_txt.lines();
        let h = SidecarHeader::parse(lines.next().unwrap()).unwrap();
        let mut junctions: Vec<(String, Side, Vec<EvidenceRow>)> = Vec::new();
        let mut ix: rustc_hash::FxHashMap<(LocusKey, Side), usize> = Default::default();
        for l in lines {
            let (fi, line) = l.split_once('\t').unwrap();
            let r = EvidenceRow::parse(line, &h, fi.parse().unwrap(), &contigs).unwrap().unwrap();
            let k = (r.locus, r.side);
            let g = *ix.entry(k).or_insert_with(|| {
                junctions.push((r.locus.name(&contigs), r.side, Vec::new()));
                junctions.len() - 1
            });
            junctions[g].2.push(r);
        }
        let js = doc["junctions"].as_array().unwrap();
        assert_eq!(js.len(), junctions.len());
        for (j, (l, s, _)) in js.iter().zip(&junctions) {
            assert_eq!(j["locus"].as_str().unwrap(), l);
            assert_eq!(j["side"].as_str().unwrap(), s.as_str());
        }
        Fixture { doc, files, contigs, junctions }
    }

    pub fn checked_in() -> Fixture {
        load(&std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src/evidence/testdata").join("p3"))
    }

    /// `P3_FULL_DIR` fixture, or None (test then skips).
    pub fn full() -> Option<Fixture> {
        std::env::var_os("P3_FULL_DIR").map(|d| load(std::path::Path::new(&d)))
    }

    /// python `cfg.get(k, default)` over one generator config dict -> a Config.
    pub fn config(c: &Value) -> crate::config::Config {
        let i = |k: &str, d: i64| c.get(k).and_then(|v| v.as_i64()).unwrap_or(d);
        let b = |k: &str, d: bool| c.get(k).and_then(|v| v.as_bool()).unwrap_or(d);
        let f = |k: &str, d: f64| c.get(k).and_then(|v| v.as_f64()).unwrap_or(d);
        crate::config::Config {
            genome_2bit: String::new(),
            exclude_files_with_many_insertions: 0,
            samtools_executable: String::new(),
            bowtie2_executable: String::new(),
            bowtie2_index: String::new(),
            bowtie2_index2: String::new(),
            bowtie2_index2_lo: String::new(),
            genotyping_max_bases: 30,
            clean_remap_max_insertion: 0,
            clean_remap_min_as: -15,
            trim_far_flank_before_remap: false,
            keep_polya_one_sided: false,
            merge_tolerance_bp: 0,
            polya_aware_clip_agreement: false,
            require_independent_fragments: false,
            min_independent_fragments: i("min_independent_fragments", 2) as usize,
            indel_aware_consensus: b("indel_aware_consensus", false),
            dup_coord_tolerance: i("dup_coord_tolerance", 5),
            dup_max_edit: i("dup_max_edit", 3),
            dup_max_edit_frac: f("dup_max_edit_frac", 0.02),
            polya_min_len: i("polya_min_len", 8) as usize,
            dup_mate_min_mapq: i("dup_mate_min_mapq", 20),
            count_short_overhang: b("count_short_overhang", false),
            short_overhang_min_bases: i("short_overhang_min_bases", 5) as usize,
            short_overhang_min_ref_mismatch: i("short_overhang_min_ref_mismatch", 2),
            short_mate_max_dist: i("short_mate_max_dist", 1000),
            short_mate_min_mapq: i("short_mate_min_mapq", 20),
            slippage_reject: false,
            far_pair_strict: false,
            far_pair_split: true,
            rte_library: String::new(),
            slippage_min_ref_run_combine: 8,
            slippage_min_str_len: 12,
            slippage_min_structured: 10,
            slippage_max_period: 6,
            slippage_junk_frac: 0.5,
            far_pair_max_tsd_deletion: 30,
            far_pair_tsd_max: 40,
            far_pair_min_polya: 10,
            far_pair_allow_antisense: false,
            far_pair_colony_frac: 0.2,
            far_pair_colony_tol: 5,
            config_dir: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn py_int_matches_python() {
        // python int(): whitespace, sign, single underscores between digits
        for (s, v) in [
            ("12", Some(12)),
            (" 13 ", Some(13)),
            ("+3", Some(3)),
            ("-0", Some(0)),
            ("-7", Some(-7)),
            ("1_6", Some(16)),
            ("1__6", None),
            ("_1", None),
            ("1_", None),
            ("", None),
            ("x", None),
            ("1.0", None),
            ("+", None),
            ("- 1", None),
        ] {
            assert_eq!(py_int_str(s), v, "{s:?}");
        }
    }

    #[test]
    fn cigar_helpers_match_python() {
        let fx = p3_fixture::checked_in();
        for c in fx.doc["cigars"].as_array().unwrap() {
            let cig = c[0].as_str().unwrap();
            assert_eq!(clip_at_from_cigar(cig, Side::Left, 40), c[1].as_i64().unwrap(), "{cig}");
            assert_eq!(clip_at_from_cigar(cig, Side::Right, 40), c[2].as_i64().unwrap(), "{cig}");
            assert_eq!(cigar_ref_len(cig), c[3].as_i64().unwrap(), "{cig}");
            let qs = [0, 2, 4, 9, 12, 30];
            for (q, e) in qs.iter().zip(c[4].as_array().unwrap()) {
                assert_eq!(ref_pos_at(100, cig, *q), e.as_i64().unwrap(), "{cig} q={q}");
            }
        }
    }

    #[test]
    fn row_methods_match_python() {
        let fx = p3_fixture::checked_in();
        let h = SidecarHeader::parse(fx.doc["header"].as_str().unwrap()).unwrap();
        let cases = fx.doc["row_cases"].as_array().unwrap();
        assert!(cases.len() >= 50);
        let ints = |v: &serde_json::Value| -> Vec<u8> { v.as_array().unwrap().iter().map(|x| x.as_i64().unwrap() as u8).collect() };
        for c in cases {
            let line = c["line"].as_str().unwrap();
            let r = EvidenceRow::parse(line, &h, 0, &fx.contigs).unwrap().unwrap();
            let (s, q) = r.outward_clip();
            assert_eq!(s, c["outward_seq"].as_str().unwrap().as_bytes(), "{line}");
            assert_eq!(q, ints(&c["outward_q"]), "{line}");
            assert_eq!(r.outward_clip_seq(), s);
            assert_eq!(r.quals(), ints(&c["quals"]), "{line}");
            assert_eq!(allele_forward_seq(&r), c["afs"].as_str().unwrap().as_bytes(), "{line}");
            assert_eq!(r.mapped(&fx.contigs), c["mapped"].as_bool().unwrap(), "{line}");
            assert_eq!(r.mate_mapped(&fx.contigs), c["mate_mapped"].as_bool().unwrap(), "{line}");
            assert_eq!(r.r12, c["r12"].as_i64().unwrap());
            assert_eq!(r.flag, c["flag"].as_i64().unwrap(), "{line}");
            assert_eq!(r.clip_at, c["clip_at"].as_i64().unwrap(), "{line}");
            for (q, e) in [0, 1, 5, 30, 200].iter().zip(c["ref_pos"].as_array().unwrap()) {
                assert_eq!(ref_pos_at(r.pos, &r.cigar, *q), e.as_i64().unwrap());
            }
        }
    }

    #[test]
    fn header_and_row_edge_cases() {
        let contigs = Interner::new();
        assert!(SidecarHeader::parse("side\trole\tfrag").is_err());
        // minimal header: python defaults for every absent column
        let h = SidecarHeader::parse("locus\tside\trole\tfrag\n").unwrap();
        let r = EvidenceRow::parse("chr1:10-20\tLEFT\tCLIP\tf1", &h, 2, &contigs).unwrap().unwrap();
        assert_eq!((r.r12, r.flag, r.pos, r.outer, r.mpos, r.tlen, r.mapq, r.clip_at), (0, 0, -1, -1, -1, 0, 0, -1));
        assert_eq!((contigs.name(r.ref_), &*r.strand, contigs.name(r.mref), &*r.mstrand, &*r.cigar), ("*", "*", "*", "*", "*"));
        assert!(r.seq.is_empty() && r.qual.is_empty() && !r.mapped(&contigs) && !r.mate_mapped(&contigs));
        // wrong field count / unparsable locus / foreign side -> skipped
        assert!(EvidenceRow::parse("chr1:10-20\tLEFT\tCLIP", &h, 0, &contigs).unwrap().is_none());
        assert!(EvidenceRow::parse("chr1:10\tLEFT\tCLIP\tf", &h, 0, &contigs).unwrap().is_none());
        assert!(EvidenceRow::parse("chr1:10-20\tMID\tCLIP\tf", &h, 0, &contigs).unwrap().is_none());
        // a missing required column errors on the first row (python KeyError)
        let h2 = SidecarHeader::parse("locus\tside\tfrag").unwrap();
        assert!(EvidenceRow::parse("chr1:10-20\tLEFT\tf", &h2, 0, &contigs).is_err());
        // unknown role kept verbatim, priority 5
        let r = EvidenceRow::parse("chr1:10-oneside_10\tRIGHT\tODD\tf", &h, 0, &contigs).unwrap().unwrap();
        assert_eq!((r.role.as_str(), r.role.priority(), r.role.is_evidence()), ("ODD", 5, false));
    }
}
