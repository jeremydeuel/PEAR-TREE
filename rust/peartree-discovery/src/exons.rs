//! Feature B: a gene/exon model for the processed-pseudogene (splice) annotation.
//!
//! Loads a simple exon annotation — whitespace-separated `contig begin end gene_id`,
//! 0-based half-open — and, per gene, ranks the exons by genomic position. A mate
//! landing site is mapped to `(gene, exon-rank)`; when a candidate's mates hit
//! several exons of one gene while skipping the introns between them, that is the
//! processed-pseudogene signature. Non-gating: consumed only by the splice sidecar.

use flate2::read::GzDecoder;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{self, BufRead};

#[derive(Clone)]
pub struct Exon {
    pub begin: i64,
    pub end: i64,
    pub gene: String,
    /// 0-based rank of this exon within its gene, ordered by genomic position.
    pub rank: usize,
}

pub struct GeneModel {
    by_contig: FxHashMap<String, Vec<Exon>>,
}

fn open_lines(path: &str) -> io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    if path.ends_with(".gz") {
        Ok(Box::new(io::BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(io::BufReader::new(file)))
    }
}

impl GeneModel {
    /// Load `contig begin end gene_id` records (extra columns ignored; blank / `#` /
    /// `track` lines skipped; malformed coordinate lines skipped, not fatal).
    pub fn load(path: &str) -> io::Result<GeneModel> {
        let reader = open_lines(path)?;
        let mut by_contig: FxHashMap<String, Vec<Exon>> = FxHashMap::default();
        // gene -> its exon spans, to assign ranks after the full file is read
        let mut gene_spans: FxHashMap<String, Vec<(String, i64, i64)>> = FxHashMap::default();
        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') || line.starts_with("track") {
                continue;
            }
            let mut it = line.split_whitespace();
            let (Some(chrom), Some(begin), Some(end), Some(gene)) = (it.next(), it.next(), it.next(), it.next()) else {
                continue;
            };
            let (Ok(begin), Ok(end)) = (begin.parse::<i64>(), end.parse::<i64>()) else {
                continue;
            };
            if end <= begin {
                continue;
            }
            gene_spans.entry(gene.to_string()).or_default().push((chrom.to_string(), begin, end));
        }
        // rank each gene's exons by (contig, begin) so ordering is deterministic
        for (gene, mut spans) in gene_spans {
            spans.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
            for (rank, (chrom, begin, end)) in spans.into_iter().enumerate() {
                by_contig.entry(chrom).or_default().push(Exon { begin, end, gene: gene.clone(), rank });
            }
        }
        for v in by_contig.values_mut() {
            v.sort_by_key(|e| e.begin);
        }
        Ok(GeneModel { by_contig })
    }

    /// The exon containing `pos` on `contig`, if any (exons assumed non-overlapping).
    pub fn lookup(&self, contig: &str, pos: i64) -> Option<&Exon> {
        let v = self.by_contig.get(contig)?;
        let idx = v.partition_point(|e| e.begin <= pos);
        if idx > 0 && pos < v[idx - 1].end {
            Some(&v[idx - 1])
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn ranks_and_lookup() {
        let path = std::env::temp_dir().join("pt_exons_test.bed");
        {
            let mut f = File::create(&path).unwrap();
            // gene G1: three exons with introns between; deliberately out of order
            writeln!(f, "3\t17010000\t17010400\tG1").unwrap();
            writeln!(f, "3\t17000000\t17000400\tG1").unwrap();
            writeln!(f, "3\t17020000\t17020400\tG1").unwrap();
        }
        let gm = GeneModel::load(path.to_str().unwrap()).unwrap();
        let e0 = gm.lookup("3", 17000050).unwrap();
        let e1 = gm.lookup("3", 17010050).unwrap();
        let e2 = gm.lookup("3", 17020050).unwrap();
        assert_eq!((e0.gene.as_str(), e0.rank), ("G1", 0));
        assert_eq!(e1.rank, 1);
        assert_eq!(e2.rank, 2);
        assert!(gm.lookup("3", 17005000).is_none()); // intron
        assert!(gm.lookup("3", 17000400).is_none()); // end exclusive
        std::fs::remove_file(&path).ok();
    }
}
