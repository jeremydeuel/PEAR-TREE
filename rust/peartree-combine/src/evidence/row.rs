//! Evidence sidecar rows. OWNER: P3.
//!
//! Mirrors `EvidenceRow`, `_int`, `_clip_at_from_cigar`, `_ref_pos_at`, `_junction_pos` cigar
//! helpers and `allele_forward_seq` of src/combine_insertions_evidence.py (lines 80-148, 245-263,
//! 364-373, 1265-1280). SPEC.md §4.1.

use crate::model::{ContigId, FileId, Interner, LocusKey, Side};

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
    pub fn parse(line: &str) -> Result<SidecarHeader, String> {
        todo!("P3")
    }
}

impl EvidenceRow {
    /// Parse one data line (without the trailing "\n"). Returns Ok(None) when the column count
    /// differs from the header's (python skips the row) or the locus does not parse.
    pub fn parse(line: &str, h: &SidecarHeader, file: FileId, contigs: &Interner) -> Result<Option<EvidenceRow>, String> {
        todo!("P3")
    }

    /// `mapped`: ref not in ("*", "") and pos >= 0 and not flag & 0x4.
    pub fn mapped(&self, contigs: &Interner) -> bool {
        todo!("P3")
    }

    /// `mate_mapped`: mref not in ("*", "") and mpos >= 0 and not flag & 0x8.
    pub fn mate_mapped(&self, contigs: &Interner) -> bool {
        todo!("P3")
    }

    /// `quals()`: `ord(c)-33` per char when len(qual) == len(seq), else `[30] * len(seq)`.
    pub fn quals(&self) -> Vec<u8> {
        todo!("P3")
    }

    /// `outward_clip()`: (seq, quals) of the junction clip, outward. `at = clip_at` unless not
    /// `0 <= at <= len(seq)`, then `_clip_at_from_cigar`. RIGHT: `seq[at:]`, `q[at:]`;
    /// LEFT: `revcomp(seq[:at])`, `q[:at][::-1]`.
    pub fn outward_clip(&self) -> (Vec<u8>, Vec<u8>) {
        todo!("P3")
    }
}

/// `_clip_at_from_cigar(cigar, side, n)`: ops = regex `(\d+)([MIDNSHP=X])` findall; none ->
/// 0 (LEFT) / n (RIGHT); LEFT: first op S -> its length else 0; RIGHT: last op S -> n - len else n.
pub fn clip_at_from_cigar(cigar: &str, side: Side, n: usize) -> i64 {
    todo!("P3")
}

/// Sum of `(\d+)([MDN=X])` lengths (python `ref_len` in `_junction_pos` / SHORT check).
/// NB the regex is applied to the raw cigar string with findall -- a malformed cigar still
/// yields whatever digit+op pairs it contains.
pub fn cigar_ref_len(cigar: &str) -> i64 {
    todo!("P3")
}

/// `_ref_pos_at(pos, cigar, q)` (evidence.py:245).
pub fn ref_pos_at(pos: i64, cigar: &str, q: i64) -> i64 {
    todo!("P3")
}

/// `allele_forward_seq(r)` (evidence.py:1265): seq unchanged unless seq non-empty and != "*",
/// role is MATE/POLYA and flag & 0x1; then stored if `bool(flag & 0x10) == (not flag & 0x20)`
/// else revcomp.
pub fn allele_forward_seq(r: &EvidenceRow) -> Vec<u8> {
    todo!("P3")
}
