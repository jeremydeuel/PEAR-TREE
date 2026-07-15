//! Local-coverage estimator shared by SPEC-3 (pileup mask) and SPEC-4 (adaptive
//! evidence floor).
//!
//! Coverage is binned by read start (bins of `bin_size` bp). The genome-wide
//! median is taken over populated bins (a deterministic subsample bounds the cost,
//! standing in for the register's "3000 random locations"). A local query returns
//! the bin count at a position; SPEC-3/4 compare that to a multiple of the median.
//!
//! Bin counts are a proxy for depth (reads-started-per-bin), but SPEC-3/4 only use
//! the *ratio* local/median, in which the bin/read-length scaling cancels.

use rustc_hash::FxHashMap;

#[derive(Default)]
pub struct Coverage {
    bin_size: i64,
    by_contig: FxHashMap<String, Vec<u32>>,
    median: f64,
}

impl Coverage {
    pub fn new(bin_size: i64) -> Self {
        Coverage {
            bin_size: bin_size.max(1),
            by_contig: FxHashMap::default(),
            median: 0.0,
        }
    }

    #[inline]
    fn bin(&self, pos: i64) -> usize {
        (pos / self.bin_size) as usize
    }

    /// Install the per-bin counts for one contig (resolved once per contig in the
    /// estimation pass, so no per-read name allocation).
    pub fn set_contig(&mut self, name: String, bins: Vec<u32>) {
        self.by_contig.insert(name, bins);
    }

    /// Compute the median over populated bins; `sample_size` bounds the cost.
    pub fn finalize(&mut self, sample_size: usize) {
        let mut counts: Vec<u32> = self
            .by_contig
            .values()
            .flat_map(|v| v.iter().copied())
            .filter(|&c| c > 0)
            .collect();
        if counts.is_empty() {
            self.median = 0.0;
            return;
        }
        if sample_size > 0 && counts.len() > sample_size {
            let stride = counts.len() / sample_size;
            counts = counts.iter().step_by(stride.max(1)).copied().collect();
        }
        counts.sort_unstable();
        self.median = counts[counts.len() / 2] as f64;
    }

    pub fn median(&self) -> f64 {
        self.median
    }

    /// Local bin count at (contig, pos); 0 if the position was never covered.
    pub fn local(&self, contig: &str, pos: i64) -> u32 {
        if pos < 0 {
            return 0;
        }
        self.by_contig.get(contig).and_then(|v| v.get(self.bin(pos))).copied().unwrap_or(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_and_local() {
        let mut cov = Coverage::new(100);
        // contig "1": bins 0..5 with counts 10,10,10,10,50 (a pileup at bin 4)
        cov.set_contig("1".into(), vec![10, 10, 10, 10, 50]);
        cov.finalize(0);
        assert_eq!(cov.median(), 10.0);
        assert_eq!(cov.local("1", 50), 10); // bin 0
        assert_eq!(cov.local("1", 450), 50); // bin 4 pileup
        assert_eq!(cov.local("2", 50), 0); // unknown contig
        // pileup is 5x median -> at the >5x cutoff it would be kept, at >=5x dropped
        assert!(cov.local("1", 450) as f64 > 4.0 * cov.median());
    }
}
