//! tools/rte/pseudogene.py -- processed-pseudogene proof (exon-exon junction read).
//!
//! WP-TD. Notes:
//! * python caches `cores(gene)` / `mrna(gene)` in the index; the struct here is constructed by
//!   literal (fixed fields, see tests/golden_modules.rs), so there are no cache fields: both are
//!   recomputed per call from the genome (identical results, a few `fetch`es per candidate gene).
//!   The integration stage may add a `Sync` cache if profiling asks for one.
//! * `load_exons_by_gene`: tuple sort of (contig, start, end) then `_merge`.

use crate::config::PseudogeneCfg;
use crate::genome::Genome;
use crate::sequtil::{edlib_best, rc};
use rustc_hash::FxHashMap;
use std::io::BufRead;
use std::path::Path;

/// gene -> merged exons [(contig, start, end)] (0-based half-open).
pub type ExonsByGene = FxHashMap<String, Vec<(String, i64, i64)>>;

fn lines_of(path: &Path) -> Result<Vec<String>, String> {
    let rd = crate::inputs::open_text(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    rd.lines().map(|l| l.map_err(|e| format!("{}: {e}", path.display()))).collect()
}

fn pyint(s: &str) -> Result<i64, String> {
    s.trim().parse::<i64>().map_err(|_| format!("invalid literal for int(): '{s}'"))
}

/// `load_exons_by_gene(path)`.
pub fn load_exons_by_gene(path: &Path) -> Result<ExonsByGene, String> {
    let mut genes: ExonsByGene = FxHashMap::default();
    for line in lines_of(path)? {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        if f.len() < 4 {
            continue;
        }
        genes.entry(f[3].to_string()).or_default().push((f[0].to_string(), pyint(f[1])?, pyint(f[2])?));
    }
    for v in genes.values_mut() {
        v.sort();
        *v = merge(std::mem::take(v));
    }
    Ok(genes)
}

/// `load_gene_strands(path)` -> gene -> '+'/'-' (first seen).
pub fn load_gene_strands(path: &Path) -> Result<FxHashMap<String, char>, String> {
    let mut out: FxHashMap<String, char> = FxHashMap::default();
    for line in lines_of(path)? {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        if f.len() >= 5 && (f[4] == "+" || f[4] == "-") {
            out.entry(f[3].to_string()).or_insert(f[4].chars().next().unwrap());
        }
    }
    Ok(out)
}

fn merge(ivs: Vec<(String, i64, i64)>) -> Vec<(String, i64, i64)> {
    let mut out: Vec<(String, i64, i64)> = Vec::new();
    for (c, s, e) in ivs {
        match out.last_mut() {
            Some(l) if l.0 == c && s <= l.2 => l.2 = l.2.max(e),
            _ => out.push((c, s, e)),
        }
    }
    out
}

/// pseudogene.ExonJunctionIndex
pub struct ExonJunctionIndex<'a> {
    pub cfg: PseudogeneCfg,
    pub exons: ExonsByGene,
    pub genome: Option<&'a dyn Genome>,
    pub strands: FxHashMap<String, char>,
}

impl ExonJunctionIndex<'_> {
    /// `mrna(gene)`: spliced transcript (merged exons, gene sense when the strand is known);
    /// empty without a genome.
    pub fn mrna(&self, gene: &str) -> Vec<u8> {
        let mut seq = Vec::new();
        if let Some(g) = self.genome {
            if let Some(ex) = self.exons.get(gene) {
                for (c, s, e) in ex {
                    seq.extend(g.fetch(c, *s, *e));
                }
            }
        }
        if self.strands.get(gene) == Some(&'-') {
            rc(&seq)
        } else {
            seq
        }
    }

    /// `structure(genes, five_prime_seqs, tol=15, probe=25)`-- `probe` is fixed at 25 here (the
    /// only value the annotator uses).
    pub fn structure(&self, genes: &[String], five_prime_seqs: &[Vec<u8>], tol: i64) -> Option<String> {
        const PROBE: usize = 25;
        for g in genes {
            let m = self.mrna(g);
            if m.is_empty() {
                continue;
            }
            for q in five_prime_seqs {
                if q.len() < PROBE {
                    continue;
                }
                let q = &q[..PROBE];
                let stranded = self.strands.contains_key(g);
                let Some((_, ts, te, strand)) = edlib_best(q, &m, 0.1, !stranded) else { continue };
                let start = if strand > 0 { ts } else { m.len() as i64 - te };
                return Some(if start <= tol { "FULL_LENGTH" } else { "TRUNCATED_5P" }.to_string());
            }
        }
        None
    }

    /// `cores(gene)`: [(label, core)].
    pub fn cores(&self, gene: &str) -> Vec<(String, Vec<u8>)> {
        let ov = self.cfg.exon_junction_overhang;
        let mut out = Vec::new();
        let (Some(g), Some(ex)) = (self.genome, self.exons.get(gene)) else { return out };
        let n = ex.len();
        for i in 0..n {
            let jmax = n.min((i as i64 + 2 + self.cfg.exon_junction_skip).max(0) as usize);
            for j in (i + 1)..jmax {
                let (ci, si, ei) = (&ex[i].0, ex[i].1, ex[i].2);
                let (cj, sj, ej) = (&ex[j].0, ex[j].1, ex[j].2);
                if ci != cj || sj - ei < self.cfg.exon_junction_min_intron {
                    continue;
                }
                if ei - si < ov || ej - sj < ov {
                    continue;
                }
                let mut core = g.fetch(ci, ei - ov, ei);
                core.extend(g.fetch(cj, sj, sj + ov));
                if core.len() as i64 == 2 * ov && !core.contains(&b'N') {
                    out.push((format!("{gene}:e{}-e{}", i + 1, j + 1), core));
                }
            }
        }
        out
    }

    /// `find(genes, seqs)` -> [(junction label, name)].
    pub fn find(&self, genes: &[String], seqs: &[(String, Vec<u8>)]) -> Vec<(String, String)> {
        let k_frac = self.cfg.exon_junction_max_edits as f64 / (2 * self.cfg.exon_junction_overhang) as f64;
        let mut hits = Vec::new();
        for g in genes {
            for (label, core) in self.cores(g) {
                for (name, s) in seqs {
                    if s.len() < core.len() {
                        continue;
                    }
                    if edlib_best(&core, s, k_frac + 1e-9, true).is_some() {
                        hits.push((label.clone(), name.clone()));
                    }
                }
            }
        }
        hits
    }
}
