//! `<stem>.combined.splice.tsv` (Feature B splice-hallmark re-keying). OWNER: P6.
//!
//! Mirrors `_splice_sidecar` / `write_combined_splice` (combine_insertions.py:32-75). SPEC.md §6.6.

use crate::context::Ctx;
use crate::model::{Insertion, Interner};
use crate::remap::py_int;
use rustc_hash::FxHashMap;
use std::io::Write;
use std::path::PathBuf;

/// `_splice_sidecar(combined)`: strip `.txt.gz` then append `.splice.tsv`.
pub fn splice_sidecar(combined: &str) -> PathBuf {
    PathBuf::from(format!("{}.splice.tsv", combined.strip_suffix(".txt.gz").unwrap_or(combined)))
}

/// One discovery splice row: (breakpoint, side, gene, n_exons, intron_bp, span_bp).
type Row = (i64, String, String, i64, i64, i64);

/// Read `<f>.splice.tsv` of EVERY input file (all `--discovery_files`, not only accepted ones);
/// write nothing when none exists. Uncompressed output.
pub fn write_combined_splice(insertions: &[Insertion], combined: &str, ctx: &Ctx, window: i64) -> Result<(), String> {
    let paths: Vec<&str> = ctx.files.iter().map(|f| f.path.as_str()).collect();
    write_combined_splice_with(insertions, combined, &paths, &ctx.contigs, window)
}

/// [`write_combined_splice`] with the input paths and contig names spelled out.
pub fn write_combined_splice_with(
    insertions: &[Insertion],
    combined: &str,
    input_files: &[&str],
    contigs: &Interner,
    window: i64,
) -> Result<(), String> {
    let mut src: FxHashMap<String, Vec<Row>> = FxHashMap::default();
    let mut found = false;
    for f in input_files {
        let sp = format!("{f}.splice.tsv");
        if !std::path::Path::new(&sp).exists() {
            continue;
        }
        found = true;
        let text = std::fs::read(&sp).map_err(|e| format!("cannot read {sp}: {e}"))?;
        let text = String::from_utf8_lossy(&text);
        let mut lines = text.split('\n');
        lines.next(); // header
        for line in lines {
            let line = line.strip_suffix('\r').unwrap_or(line);
            let p: Vec<&str> = line.split('\t').collect();
            if p.len() < 7 {
                continue;
            }
            let (Some(bp), Some(nex), Some(intron), Some(span)) = (py_int(p[1]), py_int(p[4]), py_int(p[5]), py_int(p[6])) else {
                continue;
            };
            src.entry(p[0].to_string()).or_default().push((bp, p[2].to_string(), p[3].to_string(), nex, intron, span));
        }
    }
    if !found {
        return Ok(());
    }
    let out_path = splice_sidecar(combined);
    let mut out: Vec<u8> = Vec::new();
    out.extend_from_slice(b"insertion\tgene\tside\tn_exons\tintron_bp\tspan_bp\n");
    let mut n = 0usize;
    for ins in insertions {
        let Some(rows) = src.get(contigs.name(ins.contig)) else { continue };
        let name = ins.name(contigs);
        for (bp, side, gene, nex, intron, span) in rows {
            let pos = if side == "RIGHT" { ins.right_pos } else { ins.left_pos };
            let Some(pos) = pos else { continue };
            if (bp - pos).abs() <= window {
                writeln!(out, "{name}\t{gene}\t{side}\t{nex}\t{intron}\t{span}").unwrap();
                n += 1;
            }
        }
    }
    std::fs::write(&out_path, &out).map_err(|e| format!("cannot write {}: {e}", out_path.display()))?;
    println!("aggregated {n} discovery splice-hallmark rows -> {}", out_path.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{InsType, Tok};

    fn ins(contigs: &Interner, c: &str, l: i64, r: i64) -> Insertion {
        Insertion {
            uid: 0,
            contig: contigs.intern(c),
            name_start: Tok::pos(l),
            name_end: Tok::pos(r),
            ty: InsType::FullInfo,
            open_side: None,
            left_clipped: None,
            left_aligned: None,
            left_pos: Some(l),
            right_clipped: None,
            right_aligned: None,
            right_pos: Some(r),
            files: vec![0],
            member_loci: vec![],
            member_sides: None,
        }
    }

    #[test]
    fn sidecar_name() {
        assert_eq!(splice_sidecar("/a/x.combined.txt.gz"), PathBuf::from("/a/x.combined.splice.tsv"));
        assert_eq!(splice_sidecar("/a/x"), PathBuf::from("/a/x.splice.tsv"));
    }

    /// Output equals what combine_insertions.write_combined_splice produced for the same input
    /// (generated with the reference venv; window 25, malformed/short rows skipped, rows of two
    /// files concatenated in file order, an insertion with several rows keeps row order).
    #[test]
    fn rekeys_rows() {
        let dir = std::env::temp_dir().join(format!("pt_splice_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let f1 = dir.join("a.txt.gz");
        let f2 = dir.join("b.txt.gz");
        std::fs::write(
            format!("{}.splice.tsv", f1.display()),
            "contig\tbreakpoint\tside\tgene\tn_exons\tintron_bp\tspan_bp\n\
             chr1\t1000\tRIGHT\tGA\t2\t100\t150\n\
             chr1\t1030\tLEFT\tGB\t3\t0_5\t151\n\
             chr1\tx\tLEFT\tBAD\t3\t5\t151\n\
             short\tline\n\
             chr2\t500\tLEFT\tGC\t4\t1\t2\n",
        )
        .unwrap();
        std::fs::write(
            format!("{}.splice.tsv", f2.display()),
            "h\n\
             chr1\t1010\tLEFT\tGD\t2\t7\t8\n\
             chr1\t1100\tRIGHT\tFAR\t2\t7\t8\n",
        )
        .unwrap();
        let contigs = Interner::new();
        let v = vec![ins(&contigs, "chr1", 1004, 1000), ins(&contigs, "chr2", 520, 530), ins(&contigs, "chr3", 5, 6)];
        let combined = dir.join("o.combined.txt.gz");
        let (p1, p2) = (f1.to_string_lossy().to_string(), f2.to_string_lossy().to_string());
        write_combined_splice_with(&v, &combined.to_string_lossy(), &[&p1, &p2], &contigs, 25).unwrap();
        let got = std::fs::read_to_string(dir.join("o.combined.splice.tsv")).unwrap();
        // chr1:1004-1000 -> RIGHT GA (|1000-1000|), LEFT GB (|1030-1004|=26 > 25 -> no), LEFT GD (|1010-1004|)
        assert_eq!(
            got,
            "insertion\tgene\tside\tn_exons\tintron_bp\tspan_bp\n\
             chr1:1004-1000\tGA\tRIGHT\t2\t100\t150\n\
             chr1:1004-1000\tGD\tLEFT\t2\t7\t8\n\
             chr2:520-530\tGC\tLEFT\t4\t1\t2\n"
        );
        // no sidecars -> nothing written
        let none = dir.join("n.combined.txt.gz");
        let missing = dir.join("zzz.txt.gz").to_string_lossy().to_string();
        write_combined_splice_with(&v, &none.to_string_lossy(), &[&missing], &contigs, 25).unwrap();
        assert!(!dir.join("n.combined.splice.tsv").exists());
        std::fs::remove_dir_all(&dir).ok();
    }
}
