//! Sorted per-contig interval index for the exclude-BED plumbing (SPEC-5).
//!
//! Positions are 0-based. Intervals are loaded from a BED file (chrom, start, end;
//! 0-based half-open), merged into a non-overlapping sorted set per contig, so a
//! membership query is a single binary search.

use rustc_hash::FxHashMap;
use std::io::{self, BufRead};

pub struct IntervalIndex {
    by_contig: FxHashMap<String, Vec<(i64, i64)>>,
}

impl IntervalIndex {
    /// Load a BED file. Extra columns and blank / `#` / `track` / `browser` lines
    /// are ignored. Malformed coordinate lines are skipped rather than fatal.
    pub fn from_bed(path: &str) -> io::Result<IntervalIndex> {
        let file = std::fs::File::open(path)?;
        let reader = io::BufReader::new(file);
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
        // sort + merge overlapping/adjacent intervals per contig
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
        Ok(IntervalIndex { by_contig })
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
}
