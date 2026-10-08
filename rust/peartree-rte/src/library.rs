//! tools/rte/library.py -- reference RTE library + minimap2 indices. FOUNDATION (implemented):
//! every other module reads it, so it is not a work package.
//!
//! Layout: see library.py / resources/rte_library/README.md. Every TSV is read by header name
//! (python csv.DictReader over the lines that are non-blank and do not start with "##", the
//! header's leading '#' stripped); FASTA names are the first header word (mappy.fastx_read);
//! iteration order is FILE order everywhere (python dicts), duplicates replace in place.
//!
//! `young_consensus_regex`: the default `^(L1HS|L1PA[23](?![0-9])|ALU_?Y|SVA)` (re.I) uses a
//! lookahead the `regex` crate lacks, so it is implemented by hand ([`default_young`]); a config
//! override is compiled with `regex` (case-insensitive, `search` semantics) and rejected with an
//! error if it needs lookaround.

use crate::genome::read_fasta;
use crate::mm::{Aligner, MapOpts};
use crate::sequtil::trailing_polya_start;
use rustc_hash::{FxHashMap, FxHashSet};
use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

/// library.element_class(name)
pub fn element_class(name: &str) -> &'static str {
    let u = name.to_ascii_uppercase();
    if u.starts_with("L1") || u.starts_with("LINE") {
        "L1"
    } else if u.starts_with("ALU") {
        "ALU"
    } else if u.starts_with("SVA") {
        "SVA"
    } else {
        "OTHER"
    }
}

/// DEFAULT_YOUNG = `^(L1HS|L1PA[23](?![0-9])|ALU_?Y|SVA)`, re.I, `search`.
pub fn default_young(name: &str) -> bool {
    let u = name.to_ascii_uppercase();
    let b = u.as_bytes();
    u.starts_with("L1HS")
        || ((u.starts_with("L1PA2") || u.starts_with("L1PA3")) && !b.get(5).is_some_and(|c| c.is_ascii_digit()))
        || u.starts_with("ALUY")
        || u.starts_with("ALU_Y")
        || u.starts_with("SVA")
}

/// One csv.DictReader row: `get(col)` is None for a column the row is too short for.
#[derive(Clone, Debug)]
pub struct TsvRow {
    cols: Arc<Vec<String>>,
    vals: Vec<String>,
}

impl TsvRow {
    pub fn get(&self, k: &str) -> Option<&str> {
        // DictReader: a duplicated column name -> the LAST one wins
        let i = self.cols.iter().rposition(|c| c == k)?;
        self.vals.get(i).map(|s| s.as_str())
    }
    /// `r.get(k) or ""`
    pub fn s(&self, k: &str) -> &str {
        self.get(k).unwrap_or("")
    }
    /// python `r.get("id") or next(iter(r.values()), None)`
    pub fn id(&self) -> Option<&str> {
        match self.get("id") {
            Some(s) if !s.is_empty() => Some(s),
            _ => self.vals.first().map(|s| s.as_str()).filter(|s| !s.is_empty()),
        }
    }
}

pub fn read_tsv(path: Option<&Path>) -> Result<Vec<TsvRow>, String> {
    let Some(p) = path else { return Ok(Vec::new()) };
    let rdr = crate::inputs::open_text(p).map_err(|e| format!("{}: {e}", p.display()))?;
    let mut lines: Vec<String> = Vec::new();
    for l in rdr.lines() {
        let l = l.map_err(|e| format!("{}: {e}", p.display()))?;
        let l = l.trim_end_matches('\r').to_string();
        if !l.trim().is_empty() && !l.starts_with("##") {
            lines.push(l);
        }
    }
    let Some(first) = lines.first() else { return Ok(Vec::new()) };
    let cols: Arc<Vec<String>> = Arc::new(first.trim_start_matches('#').split('\t').map(String::from).collect());
    Ok(lines[1..]
        .iter()
        .map(|l| TsvRow { cols: cols.clone(), vals: l.split('\t').map(String::from).collect() })
        .collect())
}

/// Insertion-ordered name -> value map (python dict from a FASTA / TSV).
#[derive(Clone, Debug)]
pub struct Ordered<V> {
    pub items: Vec<(String, V)>,
    pos: FxHashMap<String, usize>,
}

impl<V> Default for Ordered<V> {
    fn default() -> Self {
        Ordered { items: Vec::new(), pos: FxHashMap::default() }
    }
}

impl<V> Ordered<V> {
    pub fn insert(&mut self, k: String, v: V) {
        match self.pos.get(&k) {
            Some(&i) => self.items[i].1 = v,
            None => {
                self.pos.insert(k.clone(), self.items.len());
                self.items.push((k, v));
            }
        }
    }
    pub fn get(&self, k: &str) -> Option<&V> {
        self.pos.get(k).map(|&i| &self.items[i].1)
    }
    pub fn contains(&self, k: &str) -> bool {
        self.pos.contains_key(k)
    }
    pub fn len(&self) -> usize {
        self.items.len()
    }
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }
    pub fn iter(&self) -> impl Iterator<Item = (&str, &V)> {
        self.items.iter().map(|(k, v)| (k.as_str(), v))
    }
}

fn read_fa(path: Option<&Path>) -> Result<Ordered<Vec<u8>>, String> {
    let mut o = Ordered::default();
    if let Some(p) = path {
        for (n, s) in read_fasta(p)? {
            o.insert(n, s);
        }
    }
    Ok(o)
}

fn find(root: &Path, name: &str) -> Option<PathBuf> {
    for cand in [name.to_string(), format!("{name}.gz")] {
        let p = root.join(&cand);
        if p.exists() {
            return Some(p);
        }
    }
    None
}

/// library.resolve_library_path: absolute or an existing directory -> as is; else relative to
/// the repository root (compile-time location of this crate, rust/peartree-rte/../..).
pub fn resolve_library_path(root: &str) -> PathBuf {
    let p = Path::new(root);
    if root.is_empty() || p.is_absolute() || p.is_dir() {
        return p.to_path_buf();
    }
    let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let cand = repo.join(root);
    if cand.is_dir() {
        cand
    } else {
        p.to_path_buf()
    }
}

/// Which aligner (`RteLibrary.aligner(kind)`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AlignerKind {
    Consensus,
    Intact,
    Flanks3,
    Flanks5,
}

#[derive(Clone, Debug, Default)]
pub struct LibPaths {
    pub consensus: Option<PathBuf>,
    pub landmarks: Option<PathBuf>,
    pub l1: Option<PathBuf>,
    pub alu: Option<PathBuf>,
    pub sva: Option<PathBuf>,
    pub l1_tsv: Option<PathBuf>,
    pub alu_tsv: Option<PathBuf>,
    pub sva_tsv: Option<PathBuf>,
    pub active: Option<PathBuf>,
    pub sources: Option<PathBuf>,
    pub flanks3: Option<PathBuf>,
    pub flanks5: Option<PathBuf>,
}

enum Young {
    Default,
    Regex(regex::Regex),
}

/// library.RteLibrary
pub struct RteLibrary {
    pub root: PathBuf,
    pub paths: LibPaths,
    /// name -> upper-case sequence, file order
    pub consensus: Ordered<Vec<u8>>,
    pub cons_class: FxHashMap<String, &'static str>,
    /// element 3' end on the consensus (start of its trailing poly-A)
    pub cons_end: FxHashMap<String, usize>,
    pub young: FxHashSet<String>,
    /// consensus -> [(feature, start, end)] 0-based half-open, file order
    pub landmarks: FxHashMap<String, Vec<(String, i64, i64)>>,
    pub intact: Ordered<Vec<u8>>,
    pub intact_class: FxHashMap<String, &'static str>,
    pub intact_meta: Ordered<TsvRow>,
    pub active: Ordered<TsvRow>,
    source_to_intact: FxHashMap<String, String>,
    pub active_classes: FxHashSet<&'static str>,
    pub sources: Ordered<TsvRow>,
    /// Tubio 2014 S7 polymorphic L1 positions on hs1 (contig, pos, id)
    pub polymorphic_l1: Vec<(String, i64, String)>,
    /// case kept: lower = soft-masked
    pub flanks3: Ordered<Vec<u8>>,
    pub flanks5: Ordered<Vec<u8>>,
    aligners: [OnceLock<Option<Aligner>>; 4],
}

impl RteLibrary {
    /// `RteLibrary(root, cfg)` (`young_consensus_regex` from the config).
    pub fn open(root: &str, young_regex: Option<&str>) -> Result<RteLibrary, String> {
        let root = resolve_library_path(root);
        if !root.is_dir() {
            return Err(format!("RTE library directory {:?} not found", root.display()));
        }
        let f = |n: &str| find(&root, n);
        let paths = LibPaths {
            consensus: f("consensus.fa"),
            landmarks: f("consensus_landmarks.tsv"),
            l1: f("l1_intact.fa"),
            alu: f("alu_y_intact.fa"),
            sva: f("sva_intact.fa"),
            l1_tsv: f("l1_intact.tsv"),
            alu_tsv: f("alu_y_intact.tsv"),
            sva_tsv: f("sva_intact.tsv"),
            active: f("active.tsv"),
            sources: f("transduction_sources.tsv"),
            flanks3: f("flanks_3p.fa"),
            flanks5: f("flanks_5p_sva.fa"),
        };
        let Some(cpath) = paths.consensus.clone() else {
            return Err(format!("{}/consensus.fa missing", root.display()));
        };
        let young_rule = match young_regex {
            None => Young::Default,
            Some(r) => Young::Regex(
                regex::RegexBuilder::new(r)
                    .case_insensitive(true)
                    .build()
                    .map_err(|e| format!("young_consensus_regex {r:?}: {e}"))?,
            ),
        };
        let mut consensus = Ordered::default();
        for (n, s) in read_fa(Some(&cpath))?.items {
            consensus.insert(n, s.to_ascii_uppercase());
        }
        let mut cons_class = FxHashMap::default();
        let mut cons_end = FxHashMap::default();
        let mut young = FxHashSet::default();
        for (n, s) in consensus.iter() {
            cons_class.insert(n.to_string(), element_class(n));
            cons_end.insert(n.to_string(), trailing_polya_start(s, 5));
            let y = match &young_rule {
                Young::Default => default_young(n),
                Young::Regex(r) => r.is_match(n),
            };
            if y {
                young.insert(n.to_string());
            }
        }
        let mut landmarks: FxHashMap<String, Vec<(String, i64, i64)>> = FxHashMap::default();
        for r in read_tsv(paths.landmarks.as_deref())? {
            let (Some(c), Some(fe), Some(s), Some(e)) = (r.get("consensus"), r.get("feature"), r.get("start"), r.get("end")) else {
                continue;
            };
            let (Ok(s), Ok(e)) = (s.trim().parse::<i64>(), e.trim().parse::<i64>()) else { continue };
            landmarks.entry(c.to_string()).or_default().push((fe.to_string(), s - 1, e));
        }
        let mut intact = Ordered::default();
        let mut intact_class = FxHashMap::default();
        let mut intact_meta = Ordered::default();
        for (fa, tsv, cls) in [
            (&paths.l1, &paths.l1_tsv, "L1"),
            (&paths.alu, &paths.alu_tsv, "ALU"),
            (&paths.sva, &paths.sva_tsv, "SVA"),
        ] {
            for (n, s) in read_fa(fa.as_deref())?.items {
                intact_class.insert(n.clone(), cls);
                intact.insert(n, s.to_ascii_uppercase());
            }
            for r in read_tsv(tsv.as_deref())? {
                if let Some(id) = r.id().map(String::from) {
                    intact_meta.insert(id, r);
                }
            }
        }
        let mut active = Ordered::default();
        let mut source_to_intact: FxHashMap<String, String> = FxHashMap::default();
        for r in read_tsv(paths.active.as_deref())? {
            let Some(id) = r.id().map(String::from) else { continue };
            if let Some(sid) = r.get("source_id") {
                if !sid.is_empty() && sid != "." && intact.contains(&id) {
                    source_to_intact.insert(sid.to_string(), id.clone());
                }
            }
            active.insert(id, r);
        }
        let mut active_classes = FxHashSet::default();
        for (k, r) in active.iter() {
            let c = intact_class.get(k).copied().unwrap_or_else(|| {
                let ec = r.get("element_class").filter(|s| !s.is_empty());
                let c2 = ec.or(r.get("class").filter(|s| !s.is_empty())).unwrap_or(k);
                element_class(c2)
            });
            active_classes.insert(c);
        }
        let mut sources = Ordered::default();
        for r in read_tsv(paths.sources.as_deref())? {
            let Some(id) = r.id().map(String::from) else { continue };
            for k in ["l1base_id", "intact_id", "element_id"] {
                if let Some(v) = r.get(k) {
                    if !v.is_empty() && v != "." && intact.contains(v) {
                        source_to_intact.entry(id.clone()).or_insert_with(|| v.to_string());
                    }
                }
            }
            sources.insert(id, r);
        }
        let mut polymorphic_l1 = Vec::new();
        for r in read_tsv(find(&root, "polymorphic_l1_candidates.tsv").as_deref())? {
            let (Some(c), Some(p)) = (r.get("hs1_chrom"), r.get("hs1_pos")) else { continue };
            let Ok(pos) = p.trim().parse::<i64>() else { continue };
            if c != "." && !c.is_empty() && pos >= 0 {
                polymorphic_l1.push((c.to_string(), pos, r.get("id").unwrap_or(".").to_string()));
            }
        }
        let flanks3 = read_fa(paths.flanks3.as_deref())?;
        let flanks5 = read_fa(paths.flanks5.as_deref())?;
        Ok(RteLibrary {
            root,
            paths,
            consensus,
            cons_class,
            cons_end,
            young,
            landmarks,
            intact,
            intact_class,
            intact_meta,
            active,
            source_to_intact,
            active_classes,
            sources,
            polymorphic_l1,
            flanks3,
            flanks5,
            aligners: Default::default(),
        })
    }

    /// `lib.cons_class.get(name, "OTHER")`-style lookup ('' when unknown, like `.get(n, "")`).
    pub fn class_of(&self, cons: &str) -> Option<&'static str> {
        self.cons_class.get(cons).copied()
    }

    pub fn is_young(&self, cons: &str) -> bool {
        self.young.contains(cons)
    }

    pub fn is_active_element(&self, intact_id: &str) -> bool {
        if self.active.contains(intact_id) {
            return true;
        }
        matches!(self.intact_class.get(intact_id), Some(c) if !self.active_classes.contains(c))
    }

    fn split_flank_name(flank: &str) -> (String, String) {
        let base = flank.split(['|', ' ', '\t', '\n', '\r', '\x0b', '\x0c']).next().unwrap_or("");
        if let Some(s) = base.strip_suffix("/+") {
            return (s.to_string(), "+".into());
        }
        if let Some(s) = base.strip_suffix("/-") {
            return (s.to_string(), "-".into());
        }
        (base.to_string(), String::new())
    }

    /// `source_for_flank(flank_name)`
    pub fn source_for_flank(&self, flank: &str) -> String {
        if self.sources.contains(flank) {
            return flank.to_string();
        }
        Self::split_flank_name(flank).0
    }

    /// `flank_strand(flank_name)`: "+"/"-" for `<id>/<strand>`, "" otherwise.
    pub fn flank_strand(&self, flank: &str) -> String {
        Self::split_flank_name(flank).1
    }

    /// `source_element(source_id)`
    pub fn source_element(&self, sid: &str) -> Option<String> {
        if let Some(r) = self.sources.get(sid) {
            for k in ["intact_id", "element_id", "element"] {
                if let Some(v) = r.get(k) {
                    if !v.is_empty() && v != "." {
                        return Some(v.to_string());
                    }
                }
            }
        }
        if let Some(v) = self.source_to_intact.get(sid) {
            return Some(v.clone());
        }
        self.intact.contains(sid).then(|| sid.to_string())
    }

    /// `landmark_at(cons, pos)`: feature name or ".".
    pub fn landmark_at(&self, cons: &str, pos: i64) -> &str {
        for (f, s, e) in self.landmarks.get(cons).map(|v| v.as_slice()).unwrap_or(&[]) {
            if *s <= pos && pos < *e {
                return f;
            }
        }
        "."
    }

    /// Lazy mappy index (`RteLibrary.aligner(kind)`, MapOpts::LIBRARY); None if empty/absent.
    pub fn aligner(&self, kind: AlignerKind) -> Option<&Aligner> {
        let slot = &self.aligners[kind as usize];
        slot.get_or_init(|| {
            let al = match kind {
                AlignerKind::Consensus => self.paths.consensus.as_deref().and_then(|p| Aligner::from_path(p, MapOpts::LIBRARY)),
                AlignerKind::Intact => {
                    let ps: Vec<&Path> = [&self.paths.l1, &self.paths.alu, &self.paths.sva].iter().filter_map(|p| p.as_deref()).collect();
                    match ps.len() {
                        0 => None,
                        1 => Aligner::from_path(ps[0], MapOpts::LIBRARY),
                        _ => multi_fasta_aligner(&ps),
                    }
                }
                AlignerKind::Flanks3 => self.paths.flanks3.as_deref().and_then(|p| Aligner::from_path(p, MapOpts::LIBRARY)),
                AlignerKind::Flanks5 => self.paths.flanks5.as_deref().and_then(|p| Aligner::from_path(p, MapOpts::LIBRARY)),
            };
            al.filter(|a| !a.names().is_empty())
        })
        .as_ref()
    }
}

/// `_multi_fasta_aligner`: concatenate the FASTAs (`>{name}\n{seq}\n`) into a temp file and index it.
fn multi_fasta_aligner(paths: &[&Path]) -> Option<Aligner> {
    use std::io::Write;
    let tmp = std::env::temp_dir().join(format!(".peartree-rte-intact.{}.{:?}.fa", std::process::id(), std::thread::current().id()));
    {
        let mut f = std::io::BufWriter::new(std::fs::File::create(&tmp).ok()?);
        for p in paths {
            for (n, s) in read_fasta(p).ok()? {
                f.write_all(b">").ok()?;
                f.write_all(n.as_bytes()).ok()?;
                f.write_all(b"\n").ok()?;
                f.write_all(&s).ok()?;
                f.write_all(b"\n").ok()?;
            }
        }
    }
    let al = Aligner::from_path(&tmp, MapOpts::LIBRARY);
    std::fs::remove_file(&tmp).ok();
    al
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> String {
        format!("{}/../../test/fixtures/rte_library", env!("CARGO_MANIFEST_DIR"))
    }

    #[test]
    fn young_rule() {
        for (n, y) in [("L1HS", true), ("L1PA2", true), ("L1PA3", true), ("L1PA23", false), ("L1PA7", false), ("ALU_Y", true),
                       ("AluYa5", true), ("ALUSX", false), ("SVA_E", true), ("l1hs_5end", true)] {
            assert_eq!(default_young(n), y, "{n}");
        }
    }

    #[test]
    fn fixture_library_loads() {
        let lib = RteLibrary::open(&fixture(), None).unwrap();
        assert!(lib.consensus.contains("L1HS"));
        assert_eq!(lib.class_of("L1HS"), Some("L1"));
        assert!(lib.cons_end["L1HS"] <= lib.consensus.get("L1HS").unwrap().len());
        assert!(lib.is_young("L1HS"));
        assert!(lib.aligner(AlignerKind::Consensus).is_some());
        assert!(lib.aligner(AlignerKind::Intact).is_some());
        let cons = lib.consensus.get("L1HS").unwrap();
        let hits = lib.aligner(AlignerKind::Consensus).unwrap().map(&cons[1000..1150]);
        assert!(hits.iter().any(|h| &*h.ctg == "L1HS" && h.r_st == 1000 && h.r_en == 1150), "{hits:?}");
    }
}
