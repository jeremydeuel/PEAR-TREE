//! tools/rte/record.py -- the structured per-insertion result. FOUNDATION (implemented, tested).
//!
//! [`RteRecord::row`] produces the 17 `COLUMNS` strings byte-identical to python `row()`
//! (`_f` number formatting via pyfmt, `detail_string()` key order). [`Detail`] is the ordered
//! `rec.detail` dict (insertion order kept, `d[k] = v` replaces in place). [`gt_changes`] is
//! annotator.gt_changes (the `gt_changed` column).

use crate::pyfmt::{py_g, py_repr, py_round};
use crate::score::ScoreInput;
use serde_json::Value;

/// A `rec.detail` value. Python stores ints, floats and strings (never bools).
#[derive(Clone, Debug, PartialEq)]
pub enum DetailValue {
    Int(i64),
    Float(f64),
    Str(String),
}

impl DetailValue {
    /// `f"{v}"` (python str())
    pub fn py_str(&self) -> String {
        match self {
            DetailValue::Int(i) => i.to_string(),
            DetailValue::Float(f) => py_repr(*f),
            DetailValue::Str(s) => s.clone(),
        }
    }
    pub fn to_json(&self) -> Value {
        match self {
            DetailValue::Int(i) => Value::from(*i),
            DetailValue::Float(f) => serde_json::Number::from_f64(*f).map(Value::Number).unwrap_or(Value::Null),
            DetailValue::Str(s) => Value::from(s.as_str()),
        }
    }
    pub fn from_json(v: &Value) -> DetailValue {
        match v {
            Value::Number(n) if n.is_i64() => DetailValue::Int(n.as_i64().unwrap()),
            Value::Number(n) => DetailValue::Float(n.as_f64().unwrap_or(0.0)),
            Value::String(s) => DetailValue::Str(s.clone()),
            Value::Bool(b) => DetailValue::Str(if *b { "True" } else { "False" }.into()),
            Value::Null => DetailValue::Str("None".into()),
            other => DetailValue::Str(other.to_string()),
        }
    }
}

impl From<i64> for DetailValue {
    fn from(v: i64) -> Self {
        DetailValue::Int(v)
    }
}
impl From<f64> for DetailValue {
    fn from(v: f64) -> Self {
        DetailValue::Float(v)
    }
}
impl From<&str> for DetailValue {
    fn from(v: &str) -> Self {
        DetailValue::Str(v.to_string())
    }
}
impl From<String> for DetailValue {
    fn from(v: String) -> Self {
        DetailValue::Str(v)
    }
}

/// Insertion-ordered dict (python dict semantics for set / setdefault / pop).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Detail(pub Vec<(String, DetailValue)>);

impl Detail {
    /// `d[k] = v`
    pub fn set(&mut self, k: &str, v: impl Into<DetailValue>) {
        let v = v.into();
        match self.0.iter_mut().find(|(x, _)| x == k) {
            Some(e) => e.1 = v,
            None => self.0.push((k.to_string(), v)),
        }
    }
    /// `d.setdefault(k, v)`
    pub fn setdefault(&mut self, k: &str, v: impl Into<DetailValue>) {
        if !self.contains(k) {
            self.0.push((k.to_string(), v.into()));
        }
    }
    pub fn get(&self, k: &str) -> Option<&DetailValue> {
        self.0.iter().find(|(x, _)| x == k).map(|(_, v)| v)
    }
    pub fn contains(&self, k: &str) -> bool {
        self.0.iter().any(|(x, _)| x == k)
    }
    /// `d.pop(k, None)`
    pub fn pop(&mut self, k: &str) -> Option<DetailValue> {
        let i = self.0.iter().position(|(x, _)| x == k)?;
        Some(self.0.remove(i).1)
    }
    /// `dict.update(other)`
    pub fn update(&mut self, other: &Detail) {
        for (k, v) in &other.0 {
            self.set(k, v.clone());
        }
    }
    pub fn to_json(&self) -> Value {
        Value::Object(self.0.iter().map(|(k, v)| (k.clone(), v.to_json())).collect())
    }
    pub fn from_json(v: &Value) -> Detail {
        let mut d = Detail::default();
        if let Some(o) = v.as_object() {
            for (k, x) in o {
                d.0.push((k.clone(), DetailValue::from_json(x)));
            }
        }
        d
    }
}

/// A number field of the record as python held it (`_f` formats int and float differently).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Num {
    Int(i64),
    Float(f64),
}

/// record._f(v, nd=4): None/'' -> '.'; float -> `f"{v:.4g}"` if |v| < 1 else
/// `f"{round(v, 2):g}"`; anything else `str(v)`.
pub fn fmt_f(v: Option<Num>) -> String {
    match v {
        None => ".".into(),
        Some(Num::Int(i)) => i.to_string(),
        Some(Num::Float(f)) => {
            if f.abs() < 1.0 {
                py_g(f, 4, false)
            } else {
                py_g(py_round(f, 2), 6, false)
            }
        }
    }
}

/// record.RteRecord
#[derive(Clone, Debug, PartialEq)]
pub struct RteRecord {
    pub insertion_id: String,
    pub element: String,
    pub structure: String,
    pub tags: Vec<String>,
    pub covered_5p: Option<i64>,
    pub covered_3p: Option<i64>,
    pub covered_intervals: Vec<(i64, i64)>,
    pub consensus: String,
    pub element_identity: Option<f64>,
    pub nearest_active: String,
    pub tsd_seq: String,
    pub tsd_len: Option<i64>,
    pub en_motif: String,
    pub en_mismatches: Option<i64>,
    pub polya_len: Option<f64>,
    pub beyond_polya: String,
    pub beyond_polya_support: i64,
    pub strand: i32,
    /// python: `round(sum(points), 2)` -- an int 0 when no point fired; `_f` prints both as "0"
    pub tprt_score: f64,
    pub tprt_points: String,
    pub tprt_call: String,
    pub detail: Detail,
    /// score.ScoreInput (kept for the recurrence re-scoring pass)
    pub score_input: Option<ScoreInput>,
    /// (contig, L, R) on the discovery genome
    pub site: (Option<String>, Option<i64>, Option<i64>),
    pub gt_reads: i64,
    pub gt_changed: String,
}

impl RteRecord {
    pub const COLUMNS: [&'static str; 17] = [
        "element", "structure", "tags", "covered_5p", "covered_3p", "element_identity", "nearest_active", "tsd_seq",
        "tsd_len", "en_motif", "en_mismatches", "polya_len", "beyond_polya", "tprt_score", "tprt_points", "tprt_call",
        "rte_detail",
    ];

    /// `RteRecord(insertion_id)` with python's defaults.
    pub fn new(id: &str) -> RteRecord {
        RteRecord {
            insertion_id: id.to_string(),
            element: "UNKNOWN".into(),
            structure: "5P_UNRESOLVED".into(),
            tags: Vec::new(),
            covered_5p: None,
            covered_3p: None,
            covered_intervals: Vec::new(),
            consensus: String::new(),
            element_identity: None,
            nearest_active: ".".into(),
            tsd_seq: String::new(),
            tsd_len: None,
            en_motif: String::new(),
            en_mismatches: None,
            polya_len: None,
            beyond_polya: String::new(),
            beyond_polya_support: 0,
            strand: 0,
            tprt_score: 0.0,
            tprt_points: ".".into(),
            tprt_call: "UNCERTAIN".into(),
            detail: Detail::default(),
            score_input: None,
            site: (None, None, None),
            gt_reads: 0,
            gt_changed: String::new(),
        }
    }

    /// `detail_string()`
    pub fn detail_string(&self) -> String {
        let mut d = self.detail.clone();
        if !self.consensus.is_empty() {
            d.setdefault("consensus", self.consensus.as_str());
        }
        if !self.covered_intervals.is_empty() {
            let c = self.covered_intervals.iter().map(|(a, b)| format!("{a}-{b}")).collect::<Vec<_>>().join(",");
            d.set("covered", c);
        }
        if self.strand != 0 {
            d.set("strand", if self.strand > 0 { "+" } else { "-" });
        }
        if self.beyond_polya_support != 0 {
            d.set("beyond_polya_support", self.beyond_polya_support);
        }
        let s = d.0.iter().map(|(k, v)| format!("{k}={}", v.py_str())).collect::<Vec<_>>().join(";");
        if s.is_empty() {
            ".".into()
        } else {
            s
        }
    }

    /// `row()`: the 17 COLUMNS values (tabs / newlines replaced by spaces).
    pub fn row(&self) -> Vec<String> {
        let or_dot = |s: &str| if s.is_empty() { ".".to_string() } else { s.to_string() };
        let vals = [
            self.element.clone(),
            self.structure.clone(),
            or_dot(&self.tags.join(",")),
            fmt_f(self.covered_5p.map(Num::Int)),
            fmt_f(self.covered_3p.map(Num::Int)),
            fmt_f(self.element_identity.map(Num::Float)),
            or_dot(&self.nearest_active),
            or_dot(&self.tsd_seq),
            fmt_f(self.tsd_len.map(Num::Int)),
            or_dot(&self.en_motif),
            fmt_f(self.en_mismatches.map(Num::Int)),
            fmt_f(self.polya_len.map(Num::Float)),
            or_dot(&self.beyond_polya),
            fmt_f(Some(Num::Float(self.tprt_score))),
            or_dot(&self.tprt_points),
            self.tprt_call.clone(),
            self.detail_string(),
        ];
        vals.into_iter().map(|v| v.replace(['\t', '\n'], " ")).collect()
    }

    /// `empty_row()`
    pub fn empty_row() -> Vec<String> {
        vec![".".to_string(); Self::COLUMNS.len()]
    }
}

/// annotator.gt_changes(base, rec): `field:old>new` items joined by ',' (tags as +added /
/// -removed), then ';' -> ' ' and '=' -> ':'; '' when the calls are the same.
pub fn gt_changes(base: &RteRecord, rec: &RteRecord) -> String {
    let dot = |s: &str| if s.is_empty() { ".".to_string() } else { s.to_string() };
    let opt = |v: Option<i64>| v.map_or(".".to_string(), |x| x.to_string());
    let mut out: Vec<String> = Vec::new();
    for (f, a, b) in [
        ("element", dot(&base.element), dot(&rec.element)),
        ("structure", dot(&base.structure), dot(&rec.structure)),
        ("consensus", dot(&base.consensus), dot(&rec.consensus)),
        ("covered_5p", opt(base.covered_5p), opt(rec.covered_5p)),
        ("covered_3p", opt(base.covered_3p), opt(rec.covered_3p)),
        ("tprt_call", dot(&base.tprt_call), dot(&rec.tprt_call)),
    ] {
        if a != b {
            out.push(format!("{f}:{a}>{b}"));
        }
    }
    let ta: std::collections::BTreeSet<&str> = base.tags.iter().map(|s| s.as_str()).collect();
    let tb: std::collections::BTreeSet<&str> = rec.tags.iter().map(|s| s.as_str()).collect();
    if ta != tb {
        let mut items: Vec<String> = tb.difference(&ta).map(|t| format!("+{t}")).collect();
        items.extend(ta.difference(&tb).map(|t| format!("-{t}")));
        out.push(format!("tags:{}", items.join("/")));
    }
    out.join(",").replace(';', " ").replace('=', ":")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn row_like_python() {
        // python:
        //   r = RteRecord("x"); r.element = "L1"; r.structure = "TRUNCATED_5P"; r.tags = ["TD3P", "TD3P_SOURCE=s1"]
        //   r.covered_5p, r.covered_3p = 4000, 6019; r.covered_intervals = [(4000, 6019)]
        //   r.consensus = "L1HS"; r.element_identity = 0.9876; r.nearest_active = "L1_a"
        //   r.tsd_seq, r.tsd_len = "AAGT", 4; r.en_motif, r.en_mismatches = "TTTTT/AA", 0
        //   r.polya_len = 25.0; r.tprt_score = 11.5; r.tprt_points = "tsd_4_25:+3"; r.tprt_call = "TPRT"
        //   r.strand = -1; r.detail = {"j5": 4000, "source_identity": 0.98761, "inv": "1-2"}
        //   r.beyond_polya_support = 2; "\t".join(r.row())
        let mut r = RteRecord::new("x");
        r.element = "L1".into();
        r.structure = "TRUNCATED_5P".into();
        r.tags = vec!["TD3P".into(), "TD3P_SOURCE=s1".into()];
        r.covered_5p = Some(4000);
        r.covered_3p = Some(6019);
        r.covered_intervals = vec![(4000, 6019)];
        r.consensus = "L1HS".into();
        r.element_identity = Some(0.9876);
        r.nearest_active = "L1_a".into();
        r.tsd_seq = "AAGT".into();
        r.tsd_len = Some(4);
        r.en_motif = "TTTTT/AA".into();
        r.en_mismatches = Some(0);
        r.polya_len = Some(25.0);
        r.tprt_score = 11.5;
        r.tprt_points = "tsd_4_25:+3".into();
        r.tprt_call = "TPRT".into();
        r.strand = -1;
        r.detail.set("j5", 4000i64);
        r.detail.set("source_identity", 0.98761);
        r.detail.set("inv", "1-2");
        r.beyond_polya_support = 2;
        assert_eq!(
            r.row().join("\t"),
            "L1\tTRUNCATED_5P\tTD3P,TD3P_SOURCE=s1\t4000\t6019\t0.9876\tL1_a\tAAGT\t4\tTTTTT/AA\t0\t25\t.\t11.5\t\
             tsd_4_25:+3\tTPRT\tj5=4000;source_identity=0.98761;inv=1-2;consensus=L1HS;covered=4000-6019;strand=-;\
             beyond_polya_support=2"
        );
        let e = RteRecord::new("y");
        // python: "\t".join(RteRecord("y").row())
        assert_eq!(e.row().join("\t"), "UNKNOWN\t5P_UNRESOLVED\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t0\t.\tUNCERTAIN\t.");
    }

    #[test]
    fn number_format() {
        assert_eq!(fmt_f(Some(Num::Float(0.98765))), "0.9877");
        assert_eq!(fmt_f(Some(Num::Float(12.345))), "12.35"); // python round(12.345, 2) = 12.35 (binary above)
        assert_eq!(fmt_f(Some(Num::Float(-4.5))), "-4.5");
        assert_eq!(fmt_f(Some(Num::Float(0.0))), "0");
        assert_eq!(fmt_f(Some(Num::Int(-3))), "-3");
    }

    #[test]
    fn gt_changes_like_python() {
        let mut a = RteRecord::new("x");
        let mut b = RteRecord::new("x");
        a.tags = vec!["TD3P".into(), "X=1".into()];
        b.tags = vec!["X=1".into(), "NOVEL_SOURCE".into()];
        b.element = "L1".into();
        b.covered_5p = Some(5);
        b.consensus = "L1HS".into();
        // python gt_changes(a, b)
        assert_eq!(gt_changes(&a, &b), "element:UNKNOWN>L1,consensus:.>L1HS,covered_5p:.>5,tags:+NOVEL_SOURCE/-TD3P");
        assert_eq!(gt_changes(&a, &a), "");
    }
}
