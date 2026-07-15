//! Discovery-relevant configuration, mirroring the `discovery` and `adapters`
//! sections of src/config.py. (Stage 2 will make these load from a file.)

pub const MIN_MAPQ: u8 = 40;
pub const MIN_CLIP_LEN: usize = 12;
pub const MIN_EVIDENCE_READS_PER_BREAKPOINT: usize = 2;
pub const MIN_ADAPTERLEN_FOR_CLIP: usize = 4;
pub const MIN_GOOD_BASES: usize = 10;
pub const EXCLUDE_SAME_CONTIG_SUPPLEMENTARY: i64 = 1000;

pub const POLYA_CUTOFF: usize = 12;

// drop clipped reads whose XA/SA shows the whole read maps contiguously
// elsewhere (not a real junction). On by default; set PEARTREE_KEEP_FULLMAP=1
// to disable (e.g. to reproduce pre-filter output).
pub fn reject_fully_mapping_reads() -> bool {
    !matches!(std::env::var("PEARTREE_KEEP_FULLMAP").as_deref(), Ok("1") | Ok("true") | Ok("True"))
}

// clip side constants (match the Python ints)
pub const CLIP_RIGHT: i32 = 1;
pub const CLIP_LEFT: i32 = 2;

pub const ADAPTERS: [&[u8]; 4] = [
    b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    b"AGATCGGAAAGCACACGTCTGAACTCCAGTCA",
    b"AGATCGGAAAGCGTCGTGTAGGGAAAGAGTGT",
];
