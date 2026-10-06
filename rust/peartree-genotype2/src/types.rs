//! Shared types of peartree-genotype2 — the contract between the modules (plans/genotype_v2/SPEC.md).
//! Every `pub` item here is fixed; change it only together with every user and say so.

/// Locus geometry class from the name (gap = R - L). Mirrors `locus_kind` in
/// src/combine_genotypes.py / tools/phylo/genotype_likelihood.py.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum LocusKind {
    Tsd,
    Blunt,
    TsdDeletion,
    FarDeletion,
    FarDuplication,
    OneSided,
}

impl LocusKind {
    pub fn as_str(self) -> &'static str {
        match self {
            LocusKind::Tsd => "TSD",
            LocusKind::Blunt => "BLUNT",
            LocusKind::TsdDeletion => "TSD_DELETION",
            LocusKind::FarDeletion => "L1_MED_DELETION",
            LocusKind::FarDuplication => "L1_MED_DUPLICATION",
            LocusKind::OneSided => "ONE_SIDED",
        }
    }
}

/// A locus of the genotyping contract. `left_pos` = L, `right_pos` = R, 0-based, so that the
/// alt haplotype is `genome[.., R) ++ INS ++ genome[L, ..)`. For a one-sided locus both
/// positions name the one real breakpoint and `left_open`/`right_open` says which end is missing
/// (`contig:L-oneside_L` -> right_open, real junction at L; `contig:oneside_R-R` -> left_open).
#[derive(Clone, Debug, PartialEq)]
pub struct Locus {
    pub name: String,
    pub chr: String,
    pub left_pos: i64,
    pub right_pos: i64,
    pub left_open: bool,
    pub right_open: bool,
}

impl Locus {
    pub fn is_one_sided(&self) -> bool {
        self.left_open || self.right_open
    }
    /// R - L (0 for one-sided).
    pub fn gap(&self) -> i64 {
        if self.is_one_sided() { 0 } else { self.right_pos - self.left_pos }
    }
    pub fn kind(&self) -> LocusKind {
        if self.is_one_sided() {
            return LocusKind::OneSided;
        }
        match self.gap() {
            g if g < -30 => LocusKind::FarDeletion,
            g if g < 0 => LocusKind::TsdDeletion,
            0 | 1 => LocusKind::Blunt,
            g if g <= 40 => LocusKind::Tsd,
            _ => LocusKind::FarDuplication,
        }
    }
}

/// One junction's consensus as combine wrote it (`combined.txt.gz` record `<locus>:L` / `:R`):
/// `ins_*` = the lowercase (clipped / inserted) part, `flank_*` = the uppercase (reference-aligned)
/// part, both in read orientation on the reference-forward strand. For `:R` the record is
/// flank ++ ins (ins = 5' part of the insertion, downstream of R); for `:L` it is ins ++ flank
/// (ins = 3' part of the insertion, ending right before L). Qualities are phred (0-93), not +33.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct JunctionConsensus {
    pub ins_seq: Vec<u8>,
    pub ins_qual: Vec<u8>,
    pub flank_seq: Vec<u8>,
    pub flank_qual: Vec<u8>,
    /// true when `ins_*` came from the 12-bp contract consensus, not from `--combined`
    pub from_contract_only: bool,
}

/// Which haplotype a segment belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Hyp {
    Ref,
    Alt,
}

/// A local haplotype segment the reads are realigned against.
#[derive(Clone, Debug, PartialEq)]
pub struct Segment {
    pub hyp: Hyp,
    /// human-readable: "REF", "REF_L", "REF_R", "ALT_L", "ALT_R", "ALT_FULL"
    pub label: &'static str,
    /// uppercase ACGTN
    pub seq: Vec<u8>,
    /// phred consensus quality per base (genome bases: `GENOME_Q`)
    pub qual: Vec<u8>,
    /// genome part(s) of the segment: (reference coordinate of segment index `idx`, idx, length).
    /// A read placed by the BAM at reference position p whose leading soft clip is s has its
    /// query base 0 at segment index `idx + (p - ref) - s` when `ref <= p < ref + len`.
    pub anchors: Vec<Anchor>,
    /// segment indices (0-based, between base i-1 and i) where a junction lies — used for
    /// `ReadObs::crosses_junction`.
    pub junction_cols: Vec<usize>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Anchor {
    pub ref_pos: i64,
    pub idx: usize,
    pub len: usize,
}

/// Phred quality assigned to reference-genome bases in a segment.
pub const GENOME_Q: u8 = 60;

/// Everything the aligner and the driver need for one locus.
#[derive(Clone, Debug, PartialEq)]
pub struct LocusModel {
    pub locus: Locus,
    pub kind: LocusKind,
    /// fetch windows, 0-based inclusive breakpoint coordinates `(lo, hi)`; the driver queries
    /// `[lo - 1, hi + 1)` like the legacy genotyper. One window unless the pair is far apart.
    pub windows: Vec<(i64, i64)>,
    /// breakpoints a read must overlap (`[bp - 1, bp + 1)`) to be scored: {L, R} or the real one
    pub breakpoints: Vec<i64>,
    pub segments: Vec<Segment>,
    /// 0.0 normally; 0.5 for a far duplication whose alt haplotype retains both reference
    /// junctions (see SPEC "alt_ref_junction_fraction")
    pub alt_ref_junction_fraction: f64,
    /// true when ALT_FULL (both junction consensuses merged: short insertion) replaced ALT_L/ALT_R
    pub alt_full: bool,
    /// Some(reason) -> the driver writes an `error` row without fetching reads
    pub error: Option<String>,
}

impl LocusModel {
    pub fn error(locus: Locus, reason: impl Into<String>) -> LocusModel {
        let kind = locus.kind();
        LocusModel {
            locus,
            kind,
            windows: Vec::new(),
            breakpoints: Vec::new(),
            segments: Vec::new(),
            alt_ref_junction_fraction: 0.0,
            alt_full: false,
            error: Some(reason.into()),
        }
    }
}

/// Classification of one read after realignment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReadClass {
    Alt,
    Ref,
    Uninformative,
    /// poorly explained by every haplotype (chimera, mismapped read): counted as `n_art`
    Unexplained,
}

/// Per-read realignment result (natural-log likelihoods).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReadObs {
    pub ll_ref: f64,
    pub ll_alt: f64,
    pub class: ReadClass,
    /// fraction of read bases aligned (matched or mismatched, not clipped) under the best hypothesis
    pub explained_frac: f64,
    /// the best alignment covers a junction column
    pub crosses_junction: bool,
}

impl ReadObs {
    pub fn llr(&self) -> f64 {
        self.ll_alt - self.ll_ref
    }
}

/// Minimal decoded alignment record handed to the scorer (fields the driver already decodes).
#[derive(Clone, Debug, PartialEq)]
pub struct ReadInput<'a> {
    pub seq: &'a [u8],
    /// phred (0-93), same length as `seq`
    pub qual: &'a [u8],
    /// BAM CIGAR as (op code 0..=8, len)
    pub cigar: &'a [(u8, usize)],
    /// 0-based leftmost aligned reference position
    pub ref_start: i64,
    /// 0-based exclusive aligned reference end
    pub ref_end: i64,
    pub reverse: bool,
}

/// A finished locus call (owner C fills it; the driver adds coverage/n_disc and formats).
#[derive(Clone, Debug, PartialEq)]
pub struct Call {
    pub genotype: &'static str,
    pub score_genotype: i64,
    pub score_alternative: i64,
    pub n_alt: i64,
    pub n_ref: i64,
    pub n_art: i64,
    pub n_uninf: i64,
    pub vaf: f64,
    pub gq: i32,
    /// PL for dosage 0, 1, 2 (Phred-scaled, min 0)
    pub pl: [i32; 3],
    /// posterior P(dosage) for 0, 1, 2
    pub post: [f64; 3],
}

/// Genotype call vocabulary — the on-disk contract consumed by src/combine_genotypes.py,
/// annotate.py and tools/phylo. Byte-identical to the legacy genotyper.
pub const GT_ARTEFACT: &str = "artefact";
pub const GT_WILDTYPE: &str = "wild-type";
pub const GT_HETEROZYGOUS: &str = "heterozygous";
pub const GT_HOMOZYGOUS: &str = "homozygous";
pub const GT_INSERTION: &str = "insertion";
pub const GT_INSERTION_UNCERTAIN: &str = "insertion?";
pub const GT_WILDTYPE_UNCERTAIN: &str = "wild-type?";
pub const GT_NO_COVERAGE: &str = "no-coverage";
pub const GT_HIGH_COVERAGE: &str = "high-coverage";
pub const GT_ERROR: &str = "error";

/// Output header (gzip TSV). The first 8 columns keep the legacy meaning.
pub const OUTPUT_HEADER: &str = "insertion\tgenotype\tscore_genotype\tscore_alternative\tcoverage\tn_alt\tn_ref\tn_art\tvaf\tgq\tpl_ref\tpl_het\tpl_hom\tn_uninf\tn_disc\n";
