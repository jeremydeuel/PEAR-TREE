//! Newick parsing (Sanger SNV-tree flavour, see tools/phylo/tree.py). Owner: E.
//!
//! Accepts tip labels (optionally `'quoted'`), internal support labels (`)100:11389`), branch
//! lengths (incl. zero), an optional root length, comments `[...]` (dropped) and the trailing `;`.
//! Unary nodes are collapsed with their lengths summed (as `tree.prune` does before every fit),
//! so each non-root node is a distinct BRANCH whose clade is the set of tips below it.
//!
//! Node ids follow tools/phylo/tree.py `assign_ids`: preorder, root `ROOT`, tips their label,
//! internal nodes `N1, N2, ...` in preorder. `Tree::nodes` is in that preorder (index 0 = root).

#[derive(Clone, Debug, PartialEq)]
pub struct Node {
    /// `ROOT`, the tip label, or `N<k>`
    pub id: String,
    /// tip label or internal support label as written ("" if none)
    pub label: String,
    /// branch length above this node (the root's own length is kept but never used as a branch)
    pub length: f64,
    pub parent: Option<usize>,
    pub children: Vec<usize>,
    /// indices into `Tree::tips` of the tips below this node (its clade), ascending
    pub clade: Vec<usize>,
}

impl Node {
    pub fn is_tip(&self) -> bool {
        self.children.is_empty()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Tree {
    /// preorder; `nodes[0]` is the root
    pub nodes: Vec<Node>,
    /// tip labels in preorder (the colony order used everywhere downstream)
    pub tips: Vec<String>,
    /// `tip_node[t]` = node index of tip t
    pub tip_node: Vec<usize>,
}

/// Intermediate recursive form.
#[derive(Debug, Default)]
struct RawNode {
    label: String,
    length: f64,
    children: Vec<RawNode>,
}

fn tokenize(text: &str) -> Result<Vec<String>, String> {
    // drop comments [...]
    let mut clean = String::with_capacity(text.len());
    let mut depth = 0usize;
    for ch in text.chars() {
        match ch {
            '[' => depth += 1,
            ']' if depth > 0 => depth -= 1,
            _ if depth == 0 => clean.push(ch),
            _ => {}
        }
    }
    let chars: Vec<char> = clean.chars().collect();
    let mut toks = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if c.is_whitespace() {
            i += 1;
        } else if matches!(c, '(' | ')' | ',' | ';' | ':') {
            toks.push(c.to_string());
            i += 1;
        } else if c == '\'' {
            let start = i + 1;
            let mut j = start;
            while j < chars.len() && chars[j] != '\'' {
                j += 1;
            }
            if j >= chars.len() {
                return Err("newick: unterminated quoted label".into());
            }
            // keep the quotes so a quoted label is never mistaken for a delimiter
            toks.push(format!("'{}'", chars[start..j].iter().collect::<String>()));
            i = j + 1;
        } else {
            let start = i;
            while i < chars.len() && !chars[i].is_whitespace() && !matches!(chars[i], '(' | ')' | ',' | ';' | ':') {
                i += 1;
            }
            toks.push(chars[start..i].iter().collect());
        }
    }
    Ok(toks)
}

struct Parser {
    toks: Vec<String>,
    pos: usize,
}

impl Parser {
    fn peek(&self) -> Option<&str> {
        self.toks.get(self.pos).map(|s| s.as_str())
    }

    fn node(&mut self) -> Result<RawNode, String> {
        let mut n = RawNode::default();
        if self.peek() == Some("(") {
            self.pos += 1;
            loop {
                let c = self.node()?;
                n.children.push(c);
                let t = self.peek().ok_or("newick: unexpected end of input inside '(...)'")?.to_string();
                self.pos += 1;
                match t.as_str() {
                    "," => continue,
                    ")" => break,
                    other => return Err(format!("newick: unexpected token {other:?}")),
                }
            }
        }
        if let Some(t) = self.peek() {
            if !matches!(t, "(" | ")" | "," | ";" | ":") {
                n.label = t.trim_matches('\'').to_string();
                self.pos += 1;
            }
        }
        if self.peek() == Some(":") {
            self.pos += 1;
            let t = self.peek().ok_or("newick: missing branch length after ':'")?;
            n.length = t.parse::<f64>().map_err(|_| format!("newick: bad branch length {t:?}"))?;
            if !n.length.is_finite() {
                return Err(format!("newick: non-finite branch length {t:?}"));
            }
            self.pos += 1;
        }
        Ok(n)
    }
}

/// Collapse unary nodes (child length += parent length), as tools/phylo/tree.py `prune`.
fn collapse(mut n: RawNode) -> RawNode {
    n.children = std::mem::take(&mut n.children).into_iter().map(collapse).collect();
    if n.children.len() == 1 {
        let mut kid = n.children.pop().unwrap();
        kid.length += n.length;
        return kid;
    }
    n
}

impl Tree {
    pub fn parse(text: &str) -> Result<Tree, String> {
        let toks = tokenize(text.trim())?;
        if toks.is_empty() {
            return Err("newick: empty input".into());
        }
        let mut p = Parser { toks, pos: 0 };
        let raw = p.node()?;
        match p.peek() {
            None | Some(";") => {}
            Some(t) => return Err(format!("newick: unexpected token {t:?} after the tree")),
        }
        let mut raw = collapse(raw);
        if raw.children.is_empty() {
            // a single colony: wrap so the root still means "all colonies"
            raw = RawNode { label: String::new(), length: 0.0, children: vec![raw] };
        }
        let mut tree = Tree { nodes: Vec::new(), tips: Vec::new(), tip_node: Vec::new() };
        let mut k = 0usize;
        tree.flatten(raw, None, &mut k)?;
        // clades: children have larger preorder indices, so a reverse sweep sees them first
        for i in (0..tree.nodes.len()).rev() {
            if tree.nodes[i].is_tip() {
                continue;
            }
            let mut clade: Vec<usize> = Vec::new();
            for &c in &tree.nodes[i].children {
                clade.extend_from_slice(&tree.nodes[c].clade);
            }
            clade.sort_unstable();
            tree.nodes[i].clade = clade;
        }
        let mut seen = std::collections::HashSet::new();
        for t in &tree.tips {
            if t.is_empty() {
                return Err("newick: a tip has no label".into());
            }
            if !seen.insert(t.as_str()) {
                return Err(format!("newick: duplicate tip label {t:?}"));
            }
        }
        Ok(tree)
    }

    pub fn read(path: &str) -> Result<Tree, String> {
        let text = std::fs::read_to_string(path).map_err(|e| format!("cannot read tree {path}: {e}"))?;
        Tree::parse(&text).map_err(|e| format!("{path}: {e}"))
    }

    fn flatten(&mut self, raw: RawNode, parent: Option<usize>, k: &mut usize) -> Result<usize, String> {
        let idx = self.nodes.len();
        let is_tip = raw.children.is_empty();
        let id = if parent.is_none() {
            "ROOT".to_string()
        } else if is_tip {
            raw.label.clone()
        } else {
            *k += 1;
            format!("N{k}")
        };
        let mut clade = Vec::new();
        if is_tip {
            clade.push(self.tips.len());
            self.tips.push(raw.label.clone());
            self.tip_node.push(idx);
        }
        self.nodes.push(Node { id, label: raw.label, length: raw.length, parent, children: Vec::new(), clade });
        for c in raw.children {
            let ci = self.flatten(c, Some(idx), k)?;
            self.nodes[idx].children.push(ci);
        }
        Ok(idx)
    }

    /// Tip labels of node `i`'s clade, in tip (preorder) order.
    pub fn clade_labels(&self, i: usize) -> Vec<&str> {
        self.nodes[i].clade.iter().map(|&t| self.tips[t].as_str()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const E2E: &str = "((S4:198.668,S8:198.668):801.332,((((S2:3.98929,S6:3.98929):94.5099,(S9:84.7001,S10:84.7001):13.7991):200.992,(S3:175.637,S5:175.637):123.854):39.1833,(S1:310.496,S7:310.496):28.1784):661.325);";

    fn clade_of<'a>(t: &'a Tree, id: &str) -> Vec<&'a str> {
        let i = t.nodes.iter().position(|n| n.id == id).unwrap();
        let mut v = t.clade_labels(i);
        v.sort();
        v
    }

    #[test]
    fn parses_e2e_tree() {
        let t = Tree::parse(E2E).unwrap();
        assert_eq!(t.tips, vec!["S4", "S8", "S2", "S6", "S9", "S10", "S3", "S5", "S1", "S7"]);
        assert_eq!(t.nodes.len(), 19);
        assert_eq!(t.nodes[0].id, "ROOT");
        assert_eq!(t.nodes[0].clade.len(), 10);
        let ids: Vec<&str> = t.nodes.iter().map(|n| n.id.as_str()).collect();
        assert_eq!(
            ids,
            vec!["ROOT", "N1", "S4", "S8", "N2", "N3", "N4", "N5", "S2", "S6", "N6", "S9", "S10", "N7", "S3", "S5", "N8", "S1", "S7"]
        );
        // same ids/clades as the Python fit (tools/phylo/tree.py), e.g. N4 = {S2,S6,S9,S10}, N8 = {S1,S7}
        assert_eq!(clade_of(&t, "N1"), vec!["S4", "S8"]);
        assert_eq!(clade_of(&t, "N2"), vec!["S1", "S10", "S2", "S3", "S5", "S6", "S7", "S9"]);
        assert_eq!(clade_of(&t, "N3"), vec!["S10", "S2", "S3", "S5", "S6", "S9"]);
        assert_eq!(clade_of(&t, "N4"), vec!["S10", "S2", "S6", "S9"]);
        assert_eq!(clade_of(&t, "N5"), vec!["S2", "S6"]);
        assert_eq!(clade_of(&t, "N6"), vec!["S10", "S9"]);
        assert_eq!(clade_of(&t, "N7"), vec!["S3", "S5"]);
        assert_eq!(clade_of(&t, "N8"), vec!["S1", "S7"]);
        assert_eq!(clade_of(&t, "S5"), vec!["S5"]);
        let n5 = t.nodes.iter().position(|n| n.id == "N5").unwrap();
        assert!((t.nodes[n5].length - 94.5099).abs() < 1e-12);
        assert_eq!(t.nodes[n5].parent.map(|p| t.nodes[p].id.as_str()), Some("N4"));
    }

    #[test]
    fn parses_sanger_style_support_zero_lengths_root_length() {
        let s = " ((PD1a:0,PD1b:12.5)100:11389,('PD 1c':3,(PD1d:0,PD1e:0)87:0)95:7):0.0 ; ";
        let t = Tree::parse(s).unwrap();
        assert_eq!(t.tips, vec!["PD1a", "PD1b", "PD 1c", "PD1d", "PD1e"]);
        let ids: Vec<&str> = t.nodes.iter().map(|n| n.id.as_str()).collect();
        assert_eq!(ids, vec!["ROOT", "N1", "PD1a", "PD1b", "N2", "PD 1c", "N3", "PD1d", "PD1e"]);
        assert_eq!(t.nodes[1].label, "100");
        assert_eq!(t.nodes[1].length, 11389.0);
        assert_eq!(t.nodes[6].label, "87");
        assert_eq!(t.nodes[6].length, 0.0);
        assert_eq!(t.nodes[2].length, 0.0);
        assert_eq!(t.clade_labels(4), vec!["PD 1c", "PD1d", "PD1e"]);
    }

    #[test]
    fn collapses_unary_and_drops_comments() {
        let t = Tree::parse("(((A:1,B:2)[&c]:3):4,C:5);").unwrap();
        // the unary node is gone; (A,B) carries 3 + 4
        assert_eq!(t.nodes.len(), 5);
        assert_eq!(t.nodes[1].id, "N1");
        assert_eq!(t.nodes[1].length, 7.0);
        assert_eq!(t.clade_labels(1), vec!["A", "B"]);
        // single tip: wrapped
        let t1 = Tree::parse("A:3;").unwrap();
        assert_eq!(t1.nodes.len(), 2);
        assert_eq!(t1.nodes[1].id, "A");
    }

    #[test]
    fn rejects_malformed() {
        assert!(Tree::parse("((A,B);").is_err());
        assert!(Tree::parse("(A:x,B);").is_err());
        assert!(Tree::parse("(A,A);").is_err());
        assert!(Tree::parse("(A,B)C D;").is_err());
        assert!(Tree::parse("").is_err());
    }
}
