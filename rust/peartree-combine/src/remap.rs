//! Consensus FASTQ writing, the two bowtie2 remaps and their filters. OWNER: P6.
//!
//! Mirrors combine_insertions.py:78-100 (`_far_flank_trimmed`), 160-261. SPEC.md §5.
//! bowtie2 / samtools stay subprocesses; BAMs are read back as SAM text via `samtools view`.

use crate::context::Ctx;
use crate::liftover::LiftOver;
use crate::model::{InsType, Insertion, Interner};
use crate::seq::{revcomp, QualSeq};
use flate2::write::GzEncoder;
use flate2::Compression;
use rustc_hash::FxHashSet;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::{Command, Stdio};

/// The consensus FASTQ text written to `<stem>.fq.gz` (first remap input) and to
/// `<stem>.combined.txt.gz` (final output): five passes over `insertions` in order --
/// FULL_INFO: `{name}:L` left_consensus then `{name}:R` right_consensus; RIGHT_POLYA: `:L`;
/// LEFT_POLYA: `:R`; RIGHT_DISC: `:L`; LEFT_DISC: `:R`. SPEC.md §6.1.
pub fn consensus_fastq(insertions: &[Insertion], ctx: &Ctx) -> Vec<u8> {
    consensus_fastq_with(insertions, &ctx.contigs)
}

/// [`consensus_fastq`] without a [`Ctx`] (only the contig names are needed).
pub fn consensus_fastq_with(insertions: &[Insertion], contigs: &Interner) -> Vec<u8> {
    let mut out: Vec<u8> = Vec::with_capacity(insertions.len() * 600);
    let mut title = String::new();
    // python builds each pass as a list of strings; the passes are emitted in this order
    for pass in [InsType::FullInfo, InsType::RightPolyA, InsType::LeftPolyA, InsType::RightDisc, InsType::LeftDisc] {
        for ins in insertions.iter().filter(|i| i.ty == pass) {
            let name = ins.name(contigs);
            let left = matches!(pass, InsType::FullInfo | InsType::RightPolyA | InsType::RightDisc);
            let right = matches!(pass, InsType::FullInfo | InsType::LeftPolyA | InsType::LeftDisc);
            if left {
                title.clear();
                title.push_str(&name);
                title.push_str(":L");
                ins.left_consensus().fastq_into(&title, &mut out);
            }
            if right {
                title.clear();
                title.push_str(&name);
                title.push_str(":R");
                ins.right_consensus().fastq_into(&title, &mut out);
            }
        }
    }
    out
}

/// Second remap input (overwrites `<stem>.fq.gz`): per insertion in order,
/// `(lc.fastq("{name}:L") if lc else "") + (rc.fastq("{name}:R") if rc else "")` with
/// (lc, rc) = far-flank-trimmed clips when `trim_far_flank_before_remap`, else the raw
/// left_clipped / right_clipped. SPEC.md §6.2.
pub fn clip_fastq(insertions: &[Insertion], ctx: &Ctx) -> Vec<u8> {
    clip_fastq_with(insertions, &ctx.contigs, ctx.cfg.trim_far_flank_before_remap)
}

/// [`clip_fastq`] without a [`Ctx`].
pub fn clip_fastq_with(insertions: &[Insertion], contigs: &Interner, trim_far_flank: bool) -> Vec<u8> {
    let mut out: Vec<u8> = Vec::with_capacity(insertions.len() * 400);
    let mut title = String::new();
    for ins in insertions {
        let name = ins.name(contigs);
        let (lc, rc): (Option<QualSeq>, Option<QualSeq>) = if trim_far_flank {
            (far_flank_trimmed(ins, b'L'), far_flank_trimmed(ins, b'R'))
        } else {
            (ins.left_clipped.clone(), ins.right_clipped.clone())
        };
        if let Some(lc) = &lc {
            title.clear();
            title.push_str(&name);
            title.push_str(":L");
            lc.fastq_into(&title, &mut out);
        }
        if let Some(rc) = &rc {
            title.clear();
            title.push_str(&name);
            title.push_str(":R");
            rc.fastq_into(&title, &mut out);
        }
    }
    out
}

/// `_far_flank_trimmed(ins, side, probe=20, min_insert=10)`. `side` b'L' / b'R'.
pub fn far_flank_trimmed(ins: &Insertion, side: u8) -> Option<QualSeq> {
    const PROBE: usize = 20;
    const MIN_INSERT: usize = 10;
    let clip = if side == b'R' { ins.right_clipped.as_ref() } else { ins.left_clipped.as_ref() }?;
    let far: Vec<u8> = if side == b'R' {
        match &ins.left_aligned {
            // str(left_aligned).upper()[:probe]
            Some(a) => a.seq.iter().take(PROBE).map(|c| c.to_ascii_uppercase()).collect(),
            None => Vec::new(),
        }
    } else {
        match &ins.right_aligned {
            // revcomp(str(right_aligned.revcomp()).upper()[-probe:])
            Some(a) => {
                let rc: Vec<u8> = revcomp(&a.seq).iter().map(|c| c.to_ascii_uppercase()).collect();
                revcomp(&rc[rc.len().saturating_sub(PROBE)..])
            }
            None => Vec::new(),
        }
    };
    if far.len() < PROBE {
        return Some(clip.clone());
    }
    let up: Vec<u8> = clip.seq.iter().map(|c| c.to_ascii_uppercase()).collect();
    match find_sub(&up, &far) {
        Some(k) if k >= MIN_INSERT => Some(QualSeq::new(clip.seq[..k].to_vec(), clip.qual[..k].to_vec())),
        _ => Some(clip.clone()),
    }
}

/// python `str.find` for ASCII bytes (first occurrence).
fn find_sub(hay: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() {
        return Some(0);
    }
    if needle.len() > hay.len() {
        return None;
    }
    hay.windows(needle.len()).position(|w| w == needle)
}

/// gzip-write `data` to `path` (compression level `level`; python uses 1 for `.fq.gz`, the
/// gzip default 9 for outputs -- only the decompressed bytes are specified).
pub fn write_gz(path: &Path, data: &[u8], level: u32) -> Result<(), String> {
    let err = |e: std::io::Error| format!("cannot write {}: {e}", path.display());
    let f = std::fs::File::create(path).map_err(err)?;
    let mut w = GzEncoder::new(BufWriter::with_capacity(1 << 20, f), Compression::new(level));
    w.write_all(data).map_err(err)?;
    let mut inner = w.finish().map_err(err)?;
    inner.flush().map_err(err)
}

/// Run `cmd` through `/bin/sh -c` (python `os.system`), echoing nothing; a non-zero exit is
/// logged, not fatal (python ignores os.system's return value).
pub fn sh(cmd: &str) {
    match Command::new("/bin/sh").arg("-c").arg(cmd).status() {
        Ok(st) if st.success() => {}
        Ok(st) => eprintln!("warning: command exited with {st}: {cmd}"),
        Err(e) => eprintln!("warning: cannot run /bin/sh: {e}"),
    }
}

/// First remap command string (byte-identical to python, combine_insertions.py:178):
/// `{bowtie2} {fq} -x {index} --end-to-end --sensitive --threads {threads} --qc-filter |
/// {samtools} view -F 4 -b -o {bam}`. Skipped by the caller when `bam` exists.
pub fn bowtie2_end_to_end_cmd(ctx: &Ctx, fq: &str, bam: &str) -> String {
    let c = &ctx.cfg;
    format!(
        "{} {} -x {} --end-to-end --sensitive --threads {} --qc-filter | {} view -F 4 -b -o {}",
        c.bowtie2_executable, fq, c.bowtie2_index, ctx.threads, c.samtools_executable, bam
    )
}

/// Second remap command (combine_insertions.py:225): `{bowtie2} {fq} -k 1000 -x {index2}
/// --local --very-fast --threads {threads} --qc-filter | {samtools} view -F 2308 -b -o {bam}`.
pub fn bowtie2_local_cmd(ctx: &Ctx, fq: &str, bam: &str) -> String {
    let c = &ctx.cfg;
    format!(
        "{} {} -k 1000 -x {} --local --very-fast --threads {} --qc-filter | {} view -F 2308 -b -o {}",
        c.bowtie2_executable, fq, c.bowtie2_index2, ctx.threads, c.samtools_executable, bam
    )
}

/// One SAM record as needed by the filters.
#[derive(Clone, Debug)]
pub struct SamRec {
    pub qname: String,
    pub flag: u16,
    pub rname: String,
    /// 0-based leftmost position
    pub pos0: i64,
    /// (len, op) pairs
    pub cigar: Vec<(u32, u8)>,
    /// AS:i tag
    pub as_tag: Option<i64>,
}

impl SamRec {
    pub fn is_unmapped(&self) -> bool {
        self.flag & 4 != 0
    }
    pub fn is_secondary(&self) -> bool {
        self.flag & 0x100 != 0
    }
    pub fn is_supplementary(&self) -> bool {
        self.flag & 0x800 != 0
    }
    pub fn is_reverse(&self) -> bool {
        self.flag & 16 != 0
    }
    /// pysam `reference_end`: start + lengths of M, D, N, =, X.
    pub fn reference_end(&self) -> i64 {
        self.pos0
            + self.cigar.iter().filter(|(_, op)| matches!(op, b'M' | b'D' | b'N' | b'=' | b'X')).map(|(l, _)| *l as i64).sum::<i64>()
    }
}

/// Parse one SAM text line. Header lines (`@`) and blank lines yield `Ok(None)`.
pub fn parse_sam_line(line: &str) -> Result<Option<SamRec>, String> {
    let line = line.trim_end_matches(['\n', '\r']);
    if line.is_empty() || line.starts_with('@') {
        return Ok(None);
    }
    let mut it = line.split('\t');
    let mut next = |what: &str| it.next().ok_or_else(|| format!("malformed SAM line (no {what}): {line}"));
    let qname = next("QNAME")?;
    let flag: u16 = next("FLAG")?.parse().map_err(|_| format!("bad SAM flag: {line}"))?;
    let rname = next("RNAME")?;
    let pos: i64 = next("POS")?.parse().map_err(|_| format!("bad SAM pos: {line}"))?;
    let _mapq = next("MAPQ")?;
    let cigar_s = next("CIGAR")?;
    for w in ["RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"] {
        next(w)?;
    }
    let mut cigar: Vec<(u32, u8)> = Vec::new();
    if cigar_s != "*" {
        let mut n: u32 = 0;
        for c in cigar_s.bytes() {
            if c.is_ascii_digit() {
                n = n * 10 + (c - b'0') as u32;
            } else {
                cigar.push((n, c));
                n = 0;
            }
        }
    }
    let mut as_tag: Option<i64> = None;
    for t in it {
        if let Some(v) = t.strip_prefix("AS:i:") {
            as_tag = v.parse().ok();
        }
    }
    Ok(Some(SamRec { qname: qname.to_string(), flag, rname: rname.to_string(), pos0: pos - 1, cigar, as_tag }))
}

/// Stream `samtools view <bam>`, calling `f` for each record (no whole-file buffering).
pub fn for_each_record<F: FnMut(SamRec) -> Result<(), String>>(samtools: &str, bam: &Path, mut f: F) -> Result<(), String> {
    let mut child = Command::new(samtools)
        .arg("view")
        .arg(bam)
        .stdout(Stdio::piped())
        .spawn()
        .map_err(|e| format!("cannot run {samtools} view {}: {e}", bam.display()))?;
    let out = child.stdout.take().ok_or("no samtools stdout")?;
    let mut rd = BufReader::with_capacity(1 << 20, out);
    let mut buf: Vec<u8> = Vec::new();
    let mut res: Result<(), String> = Ok(());
    loop {
        buf.clear();
        match rd.read_until(b'\n', &mut buf) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => {
                res = Err(format!("read error from samtools: {e}"));
                break;
            }
        }
        match parse_sam_line(&String::from_utf8_lossy(&buf)) {
            Ok(Some(r)) => {
                if let Err(e) = f(r) {
                    res = Err(e);
                    break;
                }
            }
            Ok(None) => {}
            Err(e) => {
                res = Err(e);
                break;
            }
        }
    }
    if res.is_err() {
        let _ = child.kill();
    }
    let st = child.wait().map_err(|e| format!("samtools wait: {e}"))?;
    res?;
    if !st.success() {
        return Err(format!("{samtools} view {} failed ({st})", bam.display()));
    }
    Ok(())
}

/// Stream `samtools view <bam>` and parse each record. SPEC.md §5.1.
pub fn read_bam(ctx: &Ctx, bam: &Path) -> Result<Vec<SamRec>, String> {
    let mut v = Vec::new();
    for_each_record(&ctx.cfg.samtools_executable, bam, |r| {
        v.push(r);
        Ok(())
    })?;
    Ok(v)
}

/// Longest CIGAR `I` (0 if none).
fn max_ins(rec: &SamRec) -> i64 {
    rec.cigar.iter().filter(|(_, op)| *op == b'I').map(|(l, _)| *l as i64).max().unwrap_or(0)
}

/// `read.query_name[:-2]` (python slices characters; names are ASCII).
fn strip_side(qname: &str) -> String {
    let n = qname.chars().count();
    qname.chars().take(n.saturating_sub(2)).collect()
}

fn clean_remap_hit(rec: &SamRec, max_clean_ins: i64, min_clean_as: i64) -> Option<String> {
    if rec.is_unmapped() {
        return None;
    }
    let align_score = rec.as_tag.unwrap_or(-999);
    if max_ins(rec) < max_clean_ins && align_score >= min_clean_as {
        Some(strip_side(&rec.qname))
    } else {
        None
    }
}

/// Clean-remap filter (combine_insertions.py:187-201): names (`qname[:-2]`) of mapped records
/// whose longest CIGAR I < clean_remap_max_insertion and AS (missing -> -999) >= clean_remap_min_as.
pub fn clean_remap_names(recs: &[SamRec], ctx: &Ctx) -> FxHashSet<String> {
    let (mi, ma) = (ctx.cfg.clean_remap_max_insertion, ctx.cfg.clean_remap_min_as);
    recs.iter().filter_map(|r| clean_remap_hit(r, mi, ma)).collect()
}

/// Streaming form of [`clean_remap_names`] straight from the BAM (the driver uses this one).
pub fn clean_remap_names_bam(ctx: &Ctx, bam: &Path) -> Result<FxHashSet<String>, String> {
    let (mi, ma) = (ctx.cfg.clean_remap_max_insertion, ctx.cfg.clean_remap_min_as);
    let mut set = FxHashSet::default();
    for_each_record(&ctx.cfg.samtools_executable, bam, |r| {
        if let Some(n) = clean_remap_hit(&r, mi, ma) {
            set.insert(n);
        }
        Ok(())
    })?;
    Ok(set)
}

/// python `int(str)` for the strings that can occur in locus names / sidecar columns: optional
/// surrounding whitespace, optional sign, digits with single underscores between digits.
pub fn py_int(s: &str) -> Option<i64> {
    let t = s.trim();
    let (neg, digits) = match t.strip_prefix('-') {
        Some(r) => (true, r),
        None => (false, t.strip_prefix('+').unwrap_or(t)),
    };
    if digits.is_empty() || digits.starts_with('_') || digits.ends_with('_') || digits.contains("__") {
        return None;
    }
    if !digits.bytes().all(|c| c.is_ascii_digit() || c == b'_') {
        return None;
    }
    let clean: String = digits.chars().filter(|&c| c != '_').collect();
    let v: i64 = clean.parse().ok()?;
    Some(if neg { -v } else { v })
}

/// Per-record work of the clipped-remap filter; returns the filter name, if any.
fn clipped_remap_hit(r: &SamRec, lo: &LiftOver) -> Result<Option<String>, String> {
    if r.is_secondary() || r.is_supplementary() || r.is_unmapped() {
        return Ok(None);
    }
    let parts: Vec<&str> = r.qname.split(':').collect();
    if parts.len() != 3 {
        return Err(format!("clipped remap: query name {:?} is not contig:L-R:side", r.qname));
    }
    let (contig_raw, pos_s, side) = (parts[0], parts[1], parts[2]);
    let contig =
        if contig_raw.len() < 3 || &contig_raw.as_bytes()[..3] != b"chr" { format!("chr{contig_raw}") } else { contig_raw.to_string() };
    let lr: Vec<&str> = pos_s.split('-').collect();
    if lr.len() != 2 {
        return Err(format!("clipped remap: bad position in query name {:?}", r.qname));
    }
    let tok = if side == "R" { lr[1] } else { lr[0] };
    let pos = py_int(tok).ok_or_else(|| format!("clipped remap: bad integer {tok:?} in {:?}", r.qname))?;
    let coord = if r.is_reverse() { r.reference_end() } else { r.pos0 };
    if let Some(map) = lo.convert_coordinate(&r.rname, coord) {
        for l in map {
            if contig == l.chrom && (pos - l.pos).abs() < 1000 {
                return Ok(Some(strip_side(&r.qname)));
            }
        }
    }
    Ok(None)
}

/// Clipped-remap filter (combine_insertions.py:227-258): primary mapped records only; parse
/// `qname.split(":")` = (contig, "L-R", side), prefix "chr" to the contig unless it starts with
/// "chr" (and has >= 3 chars); pos = R for side 'R' else L (python int()); lift
/// `reference_start` (forward) / `reference_end` (reverse, exclusive end) of the hit; any result
/// on the (prefixed) contig with |pos - lifted| < 1000 -> name `qname[:-2]`. SPEC.md §5.2.
pub fn clipped_remap_names(recs: &[SamRec], lo: &LiftOver) -> Result<FxHashSet<String>, String> {
    let mut set = FxHashSet::default();
    for r in recs {
        if let Some(n) = clipped_remap_hit(r, lo)? {
            set.insert(n);
        }
    }
    Ok(set)
}

/// Streaming form of [`clipped_remap_names`] straight from the BAM (the driver uses this one).
pub fn clipped_remap_names_bam(ctx: &Ctx, bam: &Path, lo: &LiftOver) -> Result<FxHashSet<String>, String> {
    let mut set = FxHashSet::default();
    for_each_record(&ctx.cfg.samtools_executable, bam, |r| {
        if let Some(n) = clipped_remap_hit(&r, lo)? {
            set.insert(n);
        }
        Ok(())
    })?;
    Ok(set)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(line: &str) -> SamRec {
        parse_sam_line(line).unwrap().unwrap()
    }

    #[test]
    fn sam_parsing() {
        let r = rec("chr1:100-120:L\t16\tchr7\t1001\t60\t10S20M5I10M3D5M\t*\t0\t0\tACGT\tIIII\tAS:i:-7\tXS:i:3\tYT:Z:UU");
        assert_eq!(r.pos0, 1000);
        assert_eq!(r.as_tag, Some(-7));
        assert_eq!(max_ins(&r), 5);
        assert!(r.is_reverse());
        assert_eq!(r.reference_end(), 1000 + 20 + 10 + 3 + 5);
        assert!(parse_sam_line("@SQ\tSN:chr1\tLN:5").unwrap().is_none());
        let r = rec("n:1-2:R\t4\t*\t0\t0\t*\t*\t0\t0\tAC\tII");
        assert!(r.is_unmapped() && r.as_tag.is_none() && r.cigar.is_empty());
    }

    #[test]
    fn clean_filter_conditions() {
        // (max_clean_ins=20, min_as=-15)
        let ok = rec("chr1:1-2:L\t0\tchr1\t10\t60\t30M\t*\t0\t0\tA\tI\tAS:i:-3");
        let big_ins = rec("chr1:3-4:L\t0\tchr1\t10\t60\t10M25I10M\t*\t0\t0\tA\tI\tAS:i:-3");
        let low_as = rec("chr1:5-6:L\t0\tchr1\t10\t60\t30M\t*\t0\t0\tA\tI\tAS:i:-16");
        let no_as = rec("chr1:7-8:L\t0\tchr1\t10\t60\t30M\t*\t0\t0\tA\tI");
        let unm = rec("chr1:9-10:L\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\tAS:i:0");
        assert_eq!(clean_remap_hit(&ok, 20, -15).as_deref(), Some("chr1:1-2"));
        assert_eq!(clean_remap_hit(&big_ins, 20, -15), None);
        assert_eq!(clean_remap_hit(&low_as, 20, -15), None);
        assert_eq!(clean_remap_hit(&no_as, 20, -15), None); // -999 < -15
        assert_eq!(clean_remap_hit(&unm, 20, -15), None);
        assert_eq!(clean_remap_hit(&no_as, 20, -1000).as_deref(), Some("chr1:7-8"));
    }

    #[test]
    fn clipped_filter() {
        // chr5 -> chr5 (+); chr6 -> chr6 reverse strand (target size 2000)
        let chain = "chain 1 chr5 5000 + 0 5000 chr5 5000 + 0 5000 1\n5000\n\
                     chain 1 chr6 2000 + 0 2000 chr6 2000 - 0 2000 2\n2000\n";
        let lo = LiftOver::from_reader(chain.as_bytes()).unwrap();
        // forward hit at 1500 near the R junction 1000 on chr5 (|1000-1500| < 1000) -> filtered
        let near = rec("5:900-1000:R\t0\tchr5\t1501\t60\t20M\t*\t0\t0\tA\tI");
        assert_eq!(clipped_remap_hit(&near, &lo).unwrap().as_deref(), Some("5:900-1000"));
        // side L uses the left coordinate (900): |900-2000| = 1100 -> kept
        let far_l = rec("5:900-1000:L\t0\tchr5\t2001\t60\t20M\t*\t0\t0\tA\tI");
        assert_eq!(clipped_remap_hit(&far_l, &lo).unwrap(), None);
        // reverse read uses reference_end: start0=1980 + 20M = 2000 ; R=1990 near, L=900 far
        let rev_r = rec("5:900-1990:R\t16\tchr5\t1981\t60\t20M\t*\t0\t0\tA\tI");
        assert_eq!(clipped_remap_hit(&rev_r, &lo).unwrap().as_deref(), Some("5:900-1990"));
        let rev_l = rec("5:900-1990:L\t16\tchr5\t1981\t60\t20M\t*\t0\t0\tA\tI");
        assert_eq!(clipped_remap_hit(&rev_l, &lo).unwrap(), None);
        // minus-strand target: chr6:100 -> 2000-1-100 = 1899 on chr6
        let minus = rec("chr6:1800-1850:R\t0\tchr6\t101\t60\t10M\t*\t0\t0\tA\tI");
        assert_eq!(clipped_remap_hit(&minus, &lo).unwrap().as_deref(), Some("chr6:1800-1850"));
        // contig mismatch (read maps to chr5, breakpoint contig chr6) -> kept
        let other = rec("chr6:1800-1850:R\t0\tchr5\t101\t60\t10M\t*\t0\t0\tA\tI");
        assert_eq!(clipped_remap_hit(&other, &lo).unwrap(), None);
        // secondary / supplementary / unmapped skipped
        for flag in [256, 2048, 4] {
            let r = rec(&format!("5:900-1000:R\t{flag}\tchr5\t1501\t60\t20M\t*\t0\t0\tA\tI"));
            assert_eq!(clipped_remap_hit(&r, &lo).unwrap(), None);
        }
        // malformed query name -> error (python ValueError)
        assert!(clipped_remap_hit(&rec("a:b:c:d\t0\tchr5\t1\t60\t5M\t*\t0\t0\tA\tI"), &lo).is_err());
        // rname unknown to the chain -> liftover None -> kept
        let unk = rec("5:900-1000:R\t0\tchrUn\t1501\t60\t20M\t*\t0\t0\tA\tI");
        assert_eq!(clipped_remap_hit(&unk, &lo).unwrap(), None);
    }

    #[test]
    fn py_int_forms() {
        assert_eq!(py_int(" 12 "), Some(12));
        assert_eq!(py_int("-3"), Some(-3));
        assert_eq!(py_int("+4"), Some(4));
        assert_eq!(py_int("1_000"), Some(1000));
        assert_eq!(py_int("1__0"), None);
        assert_eq!(py_int("_1"), None);
        assert_eq!(py_int(""), None);
        assert_eq!(py_int("1.5"), None);
    }

    fn qs(s: &str) -> QualSeq {
        QualSeq::new(s.as_bytes().to_vec(), (0..s.len()).map(|i| (i % 40) as u8).collect())
    }

    fn ins_with(left_clip: Option<&str>, left_al: Option<&str>, right_clip: Option<&str>, right_al: Option<&str>) -> Insertion {
        use crate::model::Tok;
        Insertion {
            uid: 0,
            contig: 0,
            name_start: Tok::pos(100),
            name_end: Tok::pos(110),
            ty: InsType::FullInfo,
            open_side: None,
            left_clipped: left_clip.map(qs),
            left_aligned: left_al.map(qs),
            left_pos: Some(100),
            right_clipped: right_clip.map(qs),
            right_aligned: right_al.map(qs),
            right_pos: Some(110),
            files: vec![0],
            member_loci: vec![],
            member_sides: None,
        }
    }

    #[test]
    fn far_flank() {
        let far = "ACGTTGCAACGGTTAACCGG"; // 20
        // R: far = upper(left_aligned)[:20]; the clip holds the far flank after 12 insert bases
        let ins = ins_with(None, Some(&far.to_lowercase()), Some(&format!("tttttttttttt{far}gg")), None);
        let t = far_flank_trimmed(&ins, b'R').unwrap();
        assert_eq!(&*t.seq, b"tttttttttttt");
        assert_eq!(t.qual.len(), 12);
        // k < 10 -> unchanged
        let ins = ins_with(None, Some(far), Some(&format!("ttttttttt{far}")), None);
        assert_eq!(far_flank_trimmed(&ins, b'R').unwrap().len(), 9 + 20);
        // far shorter than probe -> unchanged
        let ins = ins_with(None, Some("ACGT"), Some("ttttttttttttACGT"), None);
        assert_eq!(far_flank_trimmed(&ins, b'R').unwrap().len(), 16);
        // no aligned -> unchanged; no clip -> None
        let ins = ins_with(None, None, Some("ttttttttttttACGT"), None);
        assert_eq!(far_flank_trimmed(&ins, b'R').unwrap().len(), 16);
        assert!(far_flank_trimmed(&ins, b'L').is_none());
        // L: right_aligned `ra` with revcomp(ra) == x; far = revcomp(x) = ra
        let x = "GGGAAACCCTTTGGGAAACC";
        let ra = String::from_utf8(revcomp(x.as_bytes())).unwrap();
        let ins = ins_with(Some(&format!("acgtacgtacgtac{ra}tt")), None, None, Some(&ra));
        let t = far_flank_trimmed(&ins, b'L').unwrap();
        assert_eq!(&*t.seq, b"acgtacgtacgtac");
    }

    /// `samtools view` round trip (skipped when samtools is not on PATH).
    #[test]
    fn reads_bam_via_samtools() {
        let Ok(out) = Command::new("samtools").arg("--version").output() else { return };
        if !out.status.success() {
            return;
        }
        let dir = std::env::temp_dir().join(format!("pt_remap_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let sam = dir.join("t.sam");
        let bam = dir.join("t.bam");
        std::fs::write(
            &sam,
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n\
             1:10-20:L\t0\tchr1\t11\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tAS:i:-2\n\
             1:10-20:R\t16\tchr1\t21\t60\t2M3I3M\t*\t0\t0\tACGTACGT\tIIIIIIII\tAS:i:-9\n",
        )
        .unwrap();
        let st = Command::new("samtools").args(["view", "-b", "-o"]).arg(&bam).arg(&sam).status().unwrap();
        assert!(st.success());
        let mut recs = Vec::new();
        for_each_record("samtools", &bam, |r| {
            recs.push(r);
            Ok(())
        })
        .unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!((recs[0].pos0, recs[0].as_tag, max_ins(&recs[0])), (10, Some(-2), 0));
        assert_eq!((recs[1].pos0, recs[1].reference_end(), max_ins(&recs[1])), (20, 25, 3));
        assert!(for_each_record("samtools", &dir.join("missing.bam"), |_| Ok(())).is_err());
        std::fs::remove_dir_all(&dir).ok();
    }
}
