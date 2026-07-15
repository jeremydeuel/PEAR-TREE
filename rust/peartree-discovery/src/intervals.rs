//! Sorted per-contig interval index for the exclude-BED plumbing (SPEC-5).
//!
//! Positions are 0-based. Intervals are loaded from a BED file (chrom, start, end;
//! 0-based half-open), merged into a non-overlapping sorted set per contig, so a
//! membership query is a single binary search.

use flate2::read::GzDecoder;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{self, BufRead};

pub struct IntervalIndex {
    by_contig: FxHashMap<String, Vec<(i64, i64)>>,
}

/// Open a path as a line reader, transparently gunzipping `.gz`.
fn open_lines(path: &str) -> io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    if path.ends_with(".gz") {
        Ok(Box::new(io::BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(io::BufReader::new(file)))
    }
}

/// Sort + merge overlapping/adjacent intervals per contig into a queryable index.
fn build(mut by_contig: FxHashMap<String, Vec<(i64, i64)>>) -> IntervalIndex {
    for v in by_contig.values_mut() {
        v.sort_unstable();
        let mut merged: Vec<(i64, i64)> = Vec::with_capacity(v.len());
        for &(s, e) in v.iter() {
            match merged.last_mut() {
                Some(last) if s <= last.1 => {
                    if e > last.1 {
                        last.1 = e;
                    }
                }
                _ => merged.push((s, e)),
            }
        }
        *v = merged;
    }
    IntervalIndex { by_contig }
}

impl IntervalIndex {
    /// Load a BED file. Extra columns and blank / `#` / `track` / `browser` lines
    /// are ignored. Malformed coordinate lines are skipped rather than fatal.
    pub fn from_bed(path: &str) -> io::Result<IntervalIndex> {
        let reader = open_lines(path)?;
        let mut by_contig: FxHashMap<String, Vec<(i64, i64)>> = FxHashMap::default();
        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("track")
                || line.starts_with("browser")
            {
                continue;
            }
            let mut it = line.split_whitespace();
            let (Some(chrom), Some(start), Some(end)) = (it.next(), it.next(), it.next()) else {
                continue;
            };
            let (Ok(start), Ok(end)) = (start.parse::<i64>(), end.parse::<i64>()) else {
                continue;
            };
            if end <= start {
                continue;
            }
            by_contig.entry(chrom.to_string()).or_default().push((start, end));
        }
        Ok(build(by_contig))
    }

    /// SPEC-7: load a RepeatMasker `.out[.gz]` track, keeping only *young* copies
    /// (percent divergence <= `div_max`) — the divergence gate, NOT family
    /// membership. Columns: 2 = %div, 5 = query contig, 6 = begin (1-based),
    /// 7 = end. Header/blank lines fail the numeric parse and are skipped.
    pub fn from_repeatmasker(path: &str, div_max: f64) -> io::Result<IntervalIndex> {
        let reader = open_lines(path)?;
        let mut by_contig: FxHashMap<String, Vec<(i64, i64)>> = FxHashMap::default();
        for line in reader.lines() {
            let line = line?;
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.len() < 11 {
                continue;
            }
            let (Ok(div), Ok(start), Ok(end)) = (f[1].parse::<f64>(), f[5].parse::<i64>(), f[6].parse::<i64>()) else {
                continue;
            };
            if div > div_max || end <= start {
                continue;
            }
            // RepeatMasker is 1-based inclusive; store 0-based half-open [start-1, end)
            by_contig.entry(f[4].to_string()).or_default().push((start - 1, end));
        }
        Ok(build(by_contig))
    }

    /// True if `pos` (0-based) falls in any excluded interval on `contig`.
    pub fn contains(&self, contig: &str, pos: i64) -> bool {
        let Some(v) = self.by_contig.get(contig) else {
            return false;
        };
        // last interval whose start <= pos; merged+sorted, so only it can contain pos
        let idx = v.partition_point(|&(s, _)| s <= pos);
        idx > 0 && pos < v[idx - 1].1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn idx(spans: &[(&str, i64, i64)]) -> IntervalIndex {
        let mut by_contig: FxHashMap<String, Vec<(i64, i64)>> = FxHashMap::default();
        for &(c, s, e) in spans {
            by_contig.entry(c.to_string()).or_default().push((s, e));
        }
        for v in by_contig.values_mut() {
            v.sort_unstable();
        }
        IntervalIndex { by_contig }
    }

    #[test]
    fn membership_half_open() {
        let ix = idx(&[("1", 100, 200)]);
        assert!(!ix.contains("1", 99));
        assert!(ix.contains("1", 100)); // start inclusive
        assert!(ix.contains("1", 199));
        assert!(!ix.contains("1", 200)); // end exclusive
        assert!(!ix.contains("2", 150)); // other contig
    }

    #[test]
    fn multiple_sorted_intervals() {
        let ix = idx(&[("1", 100, 200), ("1", 500, 600)]);
        assert!(!ix.contains("1", 50));
        assert!(ix.contains("1", 150));
        assert!(!ix.contains("1", 300));
        assert!(ix.contains("1", 550));
        assert!(!ix.contains("1", 700));
    }

    #[test]
    fn repeatmasker_divergence_gate_and_coords() {
        use std::io::Write;
        let path = std::env::temp_dir().join("pt_rm_parser_test.out");
        {
            let mut f = File::create(&path).unwrap();
            writeln!(f, "   SW perc perc perc query ... header line skipped").unwrap();
            writeln!(f, " 1000  2.0 0.0 0.0 1 101 200 (x) + L1HS LINE/L1 1 100 (0) 1").unwrap(); // young
            writeln!(f, " 1000 30.0 0.0 0.0 1 301 400 (x) + L1MA LINE/L1 1 100 (0) 2").unwrap(); // old
        }
        let ix = IntervalIndex::from_repeatmasker(path.to_str().unwrap(), 5.0).unwrap();
        // young element 101-200 (1-based) -> [100, 200) 0-based half-open
        assert!(ix.contains("1", 100)); // start inclusive
        assert!(ix.contains("1", 199));
        assert!(!ix.contains("1", 200)); // end exclusive
        // old element (div 30 > 5) is not loaded
        assert!(!ix.contains("1", 350));
        std::fs::remove_file(&path).ok();
    }
}
