//! The annotate_v2 <-> peartree-rte file contract (SPEC.md "Interface"). FOUNDATION (implemented).
//!
//! Input  `--inputs <file>.jsonl[.gz]`: one JSON object per insertion, in annotate_v2's
//!        insertion order (written by tools/rte/rust_bridge.py `write_inputs`):
//!        {"locus", "left_seq", "right_seq", "pseudogene_genes": [..], "legacy_class": str|null,
//!         "sv": [rank, desc, contig, pos]|null}
//! Output `--out <file>.tsv[.gz]`: header + one row per input, same order; columns
//!        [`OUTPUT_COLUMNS`] = locus, the 17 RteRecord.COLUMNS (python `row()` strings),
//!        gt_reads, gt_changed (annotate_v2's optional columns, always written here), then the
//!        machine fields annotate_v2 / locus_class need: consensus, strand, site_contig,
//!        site_L, site_R, rte_detail_json (the typed detail dict).

use crate::annotator::InsertionInput;
use crate::inputs::open_text;
use crate::record::RteRecord;
use flate2::write::GzEncoder;
use flate2::Compression;
use serde_json::Value;
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;

pub const OUTPUT_COLUMNS_EXTRA: [&str; 8] =
    ["gt_reads", "gt_changed", "consensus", "strand", "site_contig", "site_L", "site_R", "rte_detail_json"];

/// The full output header.
pub fn output_columns() -> Vec<&'static str> {
    let mut v = vec!["locus"];
    v.extend(RteRecord::COLUMNS);
    v.extend(OUTPUT_COLUMNS_EXTRA);
    v
}

/// Parse one input JSON line.
pub fn input_from_json(v: &Value) -> Result<InsertionInput, String> {
    let o = v.as_object().ok_or("input line: expected a JSON object")?;
    let s = |k: &str| o.get(k).and_then(|x| x.as_str()).unwrap_or("").to_string();
    let title = o.get("locus").and_then(|x| x.as_str()).ok_or("input line: missing \"locus\"")?.to_string();
    let genes = match o.get("pseudogene_genes") {
        Some(Value::Array(a)) => a.iter().map(|g| g.as_str().map(String::from).ok_or("pseudogene_genes: not a string")).collect::<Result<_, _>>()?,
        None | Some(Value::Null) => Vec::new(),
        Some(x) => return Err(format!("pseudogene_genes: expected a list, got {x}")),
    };
    let legacy_class = o.get("legacy_class").and_then(|x| x.as_str()).map(String::from);
    let sv_rank = match o.get("sv") {
        Some(Value::Array(a)) => a.first().and_then(|r| r.as_i64()),
        _ => None,
    };
    Ok(InsertionInput {
        title,
        left_seq: s("left_seq").into_bytes(),
        right_seq: s("right_seq").into_bytes(),
        pseudogene_genes: genes,
        legacy_class,
        sv_rank,
    })
}

/// Read the inputs JSONL (blank lines skipped). Duplicate loci are an error (annotate_v2 keys
/// insertions by locus).
pub fn read_inputs(path: &Path) -> Result<Vec<InsertionInput>, String> {
    let rdr = open_text(path).map_err(|e| format!("cannot open {}: {e}", path.display()))?;
    let mut out = Vec::new();
    let mut seen = rustc_hash::FxHashSet::default();
    for (i, line) in rdr.lines().enumerate() {
        let line = line.map_err(|e| format!("{}: {e}", path.display()))?;
        if line.trim().is_empty() {
            continue;
        }
        let v: Value = serde_json::from_str(&line).map_err(|e| format!("{}:{}: {e}", path.display(), i + 1))?;
        let inp = input_from_json(&v).map_err(|e| format!("{}:{}: {e}", path.display(), i + 1))?;
        if !seen.insert(inp.title.clone()) {
            return Err(format!("{}:{}: duplicate locus {:?}", path.display(), i + 1, inp.title));
        }
        out.push(inp);
    }
    Ok(out)
}

/// One output row for `locus` (`rec` None -> python `empty_row()` + gt defaults).
pub fn output_row(locus: &str, rec: Option<&RteRecord>) -> Vec<String> {
    let mut row = vec![locus.replace(['\t', '\n'], " ")];
    match rec {
        Some(r) => {
            row.extend(r.row());
            row.push(r.gt_reads.to_string());
            row.push(if r.gt_changed.is_empty() { ".".into() } else { r.gt_changed.replace(['\t', '\n'], " ") });
            row.push(if r.consensus.is_empty() { ".".into() } else { r.consensus.clone() });
            row.push(r.strand.to_string());
            row.push(r.site.0.clone().unwrap_or_else(|| ".".into()));
            row.push(r.site.1.map_or(".".into(), |x| x.to_string()));
            row.push(r.site.2.map_or(".".into(), |x| x.to_string()));
            row.push(r.detail.to_json().to_string());
        }
        None => {
            row.extend(RteRecord::empty_row());
            row.extend(["0", ".", ".", "0", ".", ".", ".", "{}"].map(String::from));
        }
    }
    row
}

/// Write the output table (gzipped when the path ends in .gz).
pub fn write_output(path: &Path, rows: impl IntoIterator<Item = Vec<String>>) -> Result<(), String> {
    let f = std::fs::File::create(path).map_err(|e| format!("cannot create {}: {e}", path.display()))?;
    let we = |e: std::io::Error| format!("write {}: {e}", path.display());
    let mut w: Box<dyn Write> = if path.to_string_lossy().ends_with(".gz") {
        Box::new(BufWriter::new(GzEncoder::new(f, Compression::default())))
    } else {
        Box::new(BufWriter::new(f))
    };
    w.write_all(output_columns().join("\t").as_bytes()).map_err(we)?;
    w.write_all(b"\n").map_err(we)?;
    for r in rows {
        w.write_all(r.join("\t").as_bytes()).map_err(we)?;
        w.write_all(b"\n").map_err(we)?;
    }
    w.flush().map_err(we)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn input_json() {
        let v: Value = serde_json::from_str(
            r#"{"locus": "chr1:5-9|x", "left_seq": "acGT", "right_seq": "", "pseudogene_genes": ["G1", "G2"], "legacy_class": "unknown", "sv": [2, "inv", "chr1", 7]}"#,
        )
        .unwrap();
        let i = input_from_json(&v).unwrap();
        assert_eq!(i.title, "chr1:5-9|x");
        assert_eq!(i.left_seq, b"acGT");
        assert_eq!(i.pseudogene_genes, vec!["G1", "G2"]);
        assert_eq!(i.legacy_class.as_deref(), Some("unknown"));
        assert_eq!(i.sv_rank, Some(2));
        let i = input_from_json(&serde_json::json!({"locus": "a", "legacy_class": null, "sv": null})).unwrap();
        assert_eq!((i.legacy_class, i.sv_rank, i.left_seq.len()), (None, None, 0));
        assert!(input_from_json(&serde_json::json!({"left_seq": "A"})).is_err());
    }

    #[test]
    fn output_shape() {
        assert_eq!(output_columns().len(), 1 + 17 + 8);
        let mut r = RteRecord::new("x");
        r.detail.set("error", "ValueError: bad");
        r.site = (Some("chr1".into()), Some(5), None);
        let row = output_row("x", Some(&r));
        assert_eq!(row.len(), output_columns().len());
        assert_eq!(row[18], "0");
        assert_eq!(row[19], ".");
        assert_eq!(&row[22..25], &["chr1".to_string(), "5".to_string(), ".".to_string()]);
        assert_eq!(row[25], r#"{"error":"ValueError: bad"}"#);
        assert_eq!(output_row("y", None).len(), output_columns().len());
    }
}
