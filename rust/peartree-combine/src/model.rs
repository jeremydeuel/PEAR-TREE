//! Shared core types. IMPLEMENTED by the architect except where marked `todo!()`.
//!
//! Python object -> Rust type map (SPEC.md "Core types"):
//!   * `Insertion` (src/combine_insertions_insertion.py:40)            -> [`Insertion`]
//!   * discovery locus id string `contig:start-end` (+ tokens)         -> [`LocusKey`] (+ [`Tok`])
//!   * `(txt.gz basename, locus id)` member tuples                      -> [`Member`]
//!   * txt.gz basenames / sample names                                  -> [`FileId`] into [`InputFile`]s
//!   * contig strings                                                   -> [`ContigId`] via [`Interner`]
//!   * `id(ins)` (evidence by_id / recmap)                              -> [`Insertion::uid`]

use crate::seq::QualSeq;
use rustc_hash::FxHashMap;
use std::fmt::Write as _;
use std::sync::RwLock;

pub type ContigId = u32;
/// Index into the input-file table (`InputFile`), in command-line order. Input basenames
/// must be unique (Python keys files by basename; duplicates are rejected, SPEC.md "CLI").
pub type FileId = u32;

// ------------------------------------------------------------------ interner

/// Thread-safe string interner for contig names (discovery contigs, sidecar `ref`/`mref`
/// values incl. `*`/`=`/``, chain/BAM contig names). Names are leaked: there are at most a few
/// thousand distinct contigs, and `&'static str` keeps the hot paths allocation-free.
#[derive(Default)]
pub struct Interner {
    inner: RwLock<InternerInner>,
}

#[derive(Default)]
struct InternerInner {
    map: FxHashMap<&'static str, ContigId>,
    names: Vec<&'static str>,
}

impl Interner {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn intern(&self, s: &str) -> ContigId {
        if let Some(&id) = self.inner.read().unwrap().map.get(s) {
            return id;
        }
        let mut w = self.inner.write().unwrap();
        if let Some(&id) = w.map.get(s) {
            return id;
        }
        let leaked: &'static str = Box::leak(s.to_owned().into_boxed_str());
        let id = w.names.len() as ContigId;
        w.names.push(leaked);
        w.map.insert(leaked, id);
        id
    }

    pub fn get(&self, s: &str) -> Option<ContigId> {
        self.inner.read().unwrap().map.get(s).copied()
    }

    pub fn name(&self, id: ContigId) -> &'static str {
        self.inner.read().unwrap().names[id as usize]
    }
}

// ------------------------------------------------------------------ locus tokens

/// Kind of one end token of a discovery locus id (`contig:<start>-<end>`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum TokKind {
    /// plain integer coordinate `123`
    Pos,
    /// `polyA_123` (legacy Bp+polyA end; coordinate = poly-A read's mate)
    PolyA,
    /// `disc_123` (Feature A discordant end)
    Disc,
    /// `oneside_123` (TPRT one-sided locus; also produced by `_to_one_sided`)
    OneSide,
}

/// One end token of a locus id. Formatting is injective, so token equality == name equality.
/// NB Python keeps the raw token text in `Insertion.name`; discovery never writes leading
/// zeros / signs, so `Tok::parse` + `Tok::fmt` round-trips every real id (SPEC.md "Names").
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Tok {
    pub kind: TokKind,
    pub pos: i64,
}

impl Tok {
    pub const fn pos(p: i64) -> Tok {
        Tok { kind: TokKind::Pos, pos: p }
    }
    pub const fn one_side(p: i64) -> Tok {
        Tok { kind: TokKind::OneSide, pos: p }
    }

    /// Parse `123` / `polyA_123` / `disc_123` / `oneside_123`.
    pub fn parse(s: &[u8]) -> Option<Tok> {
        let (kind, num) = if let Some(r) = s.strip_prefix(b"oneside_") {
            (TokKind::OneSide, r)
        } else if let Some(r) = s.strip_prefix(b"polyA_") {
            (TokKind::PolyA, r)
        } else if let Some(r) = s.strip_prefix(b"disc_") {
            (TokKind::Disc, r)
        } else {
            (TokKind::Pos, s)
        };
        let pos: i64 = std::str::from_utf8(num).ok()?.parse().ok()?;
        Some(Tok { kind, pos })
    }

    pub fn fmt_into(&self, out: &mut String) {
        let pre = match self.kind {
            TokKind::Pos => "",
            TokKind::PolyA => "polyA_",
            TokKind::Disc => "disc_",
            TokKind::OneSide => "oneside_",
        };
        let _ = write!(out, "{pre}{}", self.pos);
    }
}

/// A discovery locus id `contig:start-end` (one record of one discovery file).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct LocusKey {
    pub contig: ContigId,
    pub start: Tok,
    pub end: Tok,
}

impl LocusKey {
    /// Parse a locus id as written in discovery files / sidecars (`chr1:100-120`,
    /// `chr1:100-oneside_100`, ...). `contig` is everything before the LAST ':'
    /// (python `locus.rsplit(':', 1)`). Returns None when unparsable.
    pub fn parse(s: &str, contigs: &Interner) -> Option<LocusKey> {
        let (c, rest) = s.rsplit_once(':')?;
        let (a, b) = rest.split_once('-')?;
        Some(LocusKey { contig: contigs.intern(c), start: Tok::parse(a.as_bytes())?, end: Tok::parse(b.as_bytes())? })
    }

    pub fn name(&self, contigs: &Interner) -> String {
        let mut s = String::with_capacity(32);
        s.push_str(contigs.name(self.contig));
        s.push(':');
        self.start.fmt_into(&mut s);
        s.push('-');
        self.end.fmt_into(&mut s);
        s
    }

    /// `_locus_junction(locus, side)` (combine_insertions_evidence.py:1245): the coordinate of
    /// `side`'s token -- every token kind carries one.
    pub fn junction(&self, side: Side) -> i64 {
        match side {
            Side::Left => self.start.pos,
            Side::Right => self.end.pos,
        }
    }
}

/// `(txt.gz basename, discovery locus id)` -- one contributing discovery record.
pub type Member = (FileId, LocusKey);

// ------------------------------------------------------------------ sides / types

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Side {
    Left,
    Right,
}

/// `SIDES = ("LEFT", "RIGHT")` -- the Python iteration order.
pub const SIDES: [Side; 2] = [Side::Left, Side::Right];

impl Side {
    pub fn as_str(self) -> &'static str {
        match self {
            Side::Left => "LEFT",
            Side::Right => "RIGHT",
        }
    }
    pub fn other(self) -> Side {
        match self {
            Side::Left => Side::Right,
            Side::Right => Side::Left,
        }
    }
    pub fn parse(s: &str) -> Option<Side> {
        match s {
            "LEFT" => Some(Side::Left),
            "RIGHT" => Some(Side::Right),
            _ => None,
        }
    }
}

/// A set of sides (python `member_sides` values: tuples `("LEFT",)`, `("RIGHT",)`, SIDES).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub struct SideSet(pub u8);

impl SideSet {
    pub const LEFT: SideSet = SideSet(1);
    pub const RIGHT: SideSet = SideSet(2);
    pub const BOTH: SideSet = SideSet(3);
    pub fn only(s: Side) -> SideSet {
        match s {
            Side::Left => Self::LEFT,
            Side::Right => Self::RIGHT,
        }
    }
    pub fn contains(self, s: Side) -> bool {
        self.0 & Self::only(s).0 != 0
    }
}

/// Insertion type; discriminants are the Python ints (TYPE_* in
/// combine_insertions_insertion.py:25-34).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum InsType {
    RightPolyA = 1,
    LeftPolyA = 2,
    FullInfo = 3,
    /// right end has no reads (Feature-A `disc_` or TPRT `oneside_` right token)
    RightDisc = 4,
    /// left end has no reads
    LeftDisc = 5,
}

// ------------------------------------------------------------------ input files

/// One `--discovery_files` entry.
#[derive(Clone, Debug)]
pub struct InputFile {
    /// path as given on the command line
    pub path: String,
    /// `os.path.basename(path)` -- the identity Python uses everywhere (`Insertion.files`)
    pub basename: String,
    /// `sample_name()`: basename without a trailing `.txt.gz`
    pub sample: String,
}

impl InputFile {
    pub fn new(path: &str) -> InputFile {
        let basename = std::path::Path::new(path)
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| path.to_string());
        let sample = basename.strip_suffix(".txt.gz").unwrap_or(&basename).to_string();
        InputFile { path: path.to_string(), basename, sample }
    }
}

// ------------------------------------------------------------------ insertion

/// One (merged) insertion call -- python `Insertion`. Mates are NOT stored: no combine output
/// reads them (SPEC.md "Dead data").
#[derive(Clone, Debug)]
pub struct Insertion {
    /// stable identity (python `id(ins)`): assigned by `intersect_insertions` to every
    /// surviving insertion (0..n in output order) and never changed afterwards.
    pub uid: u32,
    pub contig: ContigId,
    /// the two raw tokens of `self.name` (`contig:name_start-name_end`)
    pub name_start: Tok,
    pub name_end: Tok,
    pub ty: InsType,
    /// TPRT one-sided locus: the end WITHOUT evidence (`open_side` attribute); None otherwise
    pub open_side: Option<Side>,
    /// outward orientation (python stores `LEFT:CLIPPED.revcomp()`)
    pub left_clipped: Option<QualSeq>,
    /// reference-forward, as in the discovery file
    pub left_aligned: Option<QualSeq>,
    pub left_pos: Option<i64>,
    pub right_clipped: Option<QualSeq>,
    /// python stores `RIGHT:ALIGNED.revcomp()` (so outward from the junction)
    pub right_aligned: Option<QualSeq>,
    pub right_pos: Option<i64>,
    /// contributing discovery files (python list of basenames; may contain repeats after `+=`)
    pub files: Vec<FileId>,
    /// (file, locus id) of every record merged into this insertion, in merge order
    pub member_loci: Vec<Member>,
    /// python `member_sides` attribute; an absent attribute and an empty dict behave the same
    /// (`getattr(ins, "member_sides", None) or {}`), both are `None` here. Insertion-ordered
    /// pairs (python dict) -- order is irrelevant for lookups but kept for reproducibility.
    pub member_sides: Option<Vec<(Member, SideSet)>>,
}

impl Insertion {
    /// python `self.name`.
    pub fn name(&self, contigs: &Interner) -> String {
        LocusKey { contig: self.contig, start: self.name_start, end: self.name_end }.name(contigs)
    }

    /// The insertion's own locus id as a key (name == its LocusKey).
    pub fn locus(&self) -> LocusKey {
        LocusKey { contig: self.contig, start: self.name_start, end: self.name_end }
    }

    /// `_open_side(ins)` (combine_insertions_evidence.py:1283): `open_side`, else RIGHT for
    /// type 4 / LEFT for type 5, else None.
    pub fn open_side_eff(&self) -> Option<Side> {
        self.open_side.or(match self.ty {
            InsType::RightDisc => Some(Side::Right),
            InsType::LeftDisc => Some(Side::Left),
            _ => None,
        })
    }

    /// `_ins_junction(ins, side)` = left_pos / right_pos.
    pub fn junction(&self, side: Side) -> Option<i64> {
        match side {
            Side::Left => self.left_pos,
            Side::Right => self.right_pos,
        }
    }

    pub fn clipped(&self, side: Side) -> Option<&QualSeq> {
        match side {
            Side::Left => self.left_clipped.as_ref(),
            Side::Right => self.right_clipped.as_ref(),
        }
    }

    pub fn aligned(&self, side: Side) -> Option<&QualSeq> {
        match side {
            Side::Left => self.left_aligned.as_ref(),
            Side::Right => self.right_aligned.as_ref(),
        }
    }

    /// `member_sides.get(m, SIDES)`.
    pub fn allowed_sides(&self, m: &Member) -> SideSet {
        match &self.member_sides {
            Some(v) => v.iter().find(|(k, _)| k == m).map(|(_, s)| *s).unwrap_or(SideSet::BOTH),
            None => SideSet::BOTH,
        }
    }

    /// python `left_consensus` (combine_insertions_insertion.py:101):
    /// `left_clipped.revcomp().lower() + left_aligned`. Panics on LEFT_POLYA / LEFT_DISC like
    /// the Python assert.
    pub fn left_consensus(&self) -> QualSeq {
        assert!(
            !matches!(self.ty, InsType::LeftPolyA | InsType::LeftDisc),
            "can not extract consensus from left polyA/disc type"
        );
        let lc = self.left_clipped.as_ref().expect("left_consensus: no left_clipped");
        let la = self.left_aligned.as_ref().expect("left_consensus: no left_aligned");
        lc.revcomp().lower().concat(la)
    }

    /// python `right_consensus` (combine_insertions_insertion.py:106):
    /// `right_aligned.revcomp() + right_clipped.lower()`.
    pub fn right_consensus(&self) -> QualSeq {
        assert!(
            !matches!(self.ty, InsType::RightPolyA | InsType::RightDisc),
            "can not extract consensus from right polyA/disc type"
        );
        let rc = self.right_clipped.as_ref().expect("right_consensus: no right_clipped");
        let ra = self.right_aligned.as_ref().expect("right_consensus: no right_aligned");
        ra.revcomp().concat(&rc.lower())
    }
}
