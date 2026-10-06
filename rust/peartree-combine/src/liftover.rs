//! UCSC chain-file liftover with pyliftover 0.4.1 semantics. OWNER: P6.
//!
//! SPEC.md §5.3. pyliftover: `LiftOver(path)` opens gzip iff the path ends with `.gz`
//! (case-insensitive); each `chain` header (score, tName, tSize, tStrand(+ required), tStart,
//! tEnd, qName, qSize, qStrand, qStart, qEnd[, id]) is followed by `size dt dq` lines and a final
//! `size` line; blocks `(sfrom, sfrom+size, tfrom)`. Lines starting with '#', '\n', '\r' between
//! chains are skipped. `convert_coordinate(chrom, pos)` (strand '+'):
//!   * None if `chrom` has no chain;
//!   * else for every block with `start <= pos < end`: `tpos = tfrom + (pos - start)`, if the
//!     chain's target strand is '-': `tpos = target_size - 1 - tpos`; result
//!     `(target_name, tpos, target_strand, score)`, sorted by score descending (stable).
//! combine only asks "does ANY result on the breakpoint's contig lie within 1000 bp", so result
//! order is irrelevant there -- but it is reproduced exactly anyway.
//!
//! Index. pyliftover keeps one `IntervalTree` per source contig and its query order (before the
//! stable score sort) is the tree traversal order. The tree is reproduced faithfully (same
//! float centers, same single-interval leaves, same stable `mid` sorts) in an arena: nodes
//! 40 B, blocks 32 B, no per-node heap allocation, so a millions-of-blocks chain stays at
//! ~200 MB and loads in a couple of seconds.

use flate2::read::MultiGzDecoder;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

const NIL: u32 = u32::MAX;

#[derive(Clone, Copy)]
struct Block {
    start: i64,
    end: i64,
    tfrom: i64,
    chain: u32,
}

struct Chain {
    target_name: String,
    target_size: i64,
    target_minus: bool,
    score: i64,
}

#[derive(Clone, Copy)]
struct Node {
    center: f64,
    left: u32,
    right: u32,
    /// `single_interval` block index while the node holds exactly one interval
    single: u32,
    /// build time: linked list (insertion order) of the blocks stored in this node
    head: u32,
    tail: u32,
    /// after build: (off, len) of the node's mid list in the pools
    off: u32,
    len: u32,
    /// 0 = empty, 1 = single leaf, 2 = normal node
    state: u8,
}

pub struct LiftOver {
    chains: Vec<Chain>,
    blocks: Vec<Block>,
    nodes: Vec<Node>,
    /// mid lists: sorted by start (stable), and by end descending (stable)
    pool_start: Vec<u32>,
    pool_end: Vec<u32>,
    /// source contig -> (root node, source size = root range)
    trees: FxHashMap<String, (u32, i64)>,
}

/// One `convert_coordinate` result.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Lifted<'a> {
    pub chrom: &'a str,
    pub pos: i64,
    pub strand: u8,
    pub score: i64,
}

struct Builder {
    blocks: Vec<Block>,
    nodes: Vec<Node>,
    next: Vec<u32>,
}

impl Builder {
    fn new_node(&mut self, min_p: f64, max_p: f64) -> u32 {
        self.nodes.push(Node {
            center: (min_p + max_p) / 2.0,
            left: NIL,
            right: NIL,
            single: NIL,
            head: NIL,
            tail: NIL,
            off: 0,
            len: 0,
            state: 0,
        });
        (self.nodes.len() - 1) as u32
    }

    /// `IntervalTree.add_interval` for block `b` (length > 0 already checked by the caller).
    fn add_interval(&mut self, n: u32, b: u32, min_p: f64, max_p: f64) {
        match self.nodes[n as usize].state {
            0 => {
                let nd = &mut self.nodes[n as usize];
                nd.state = 1;
                nd.single = b;
            }
            1 => {
                let b0 = self.nodes[n as usize].single;
                self.nodes[n as usize].state = 2;
                self.nodes[n as usize].single = NIL;
                self.add_inner(n, b0, min_p, max_p);
                self.add_inner(n, b, min_p, max_p);
            }
            _ => self.add_inner(n, b, min_p, max_p),
        }
    }

    /// `IntervalTree._add_interval`. Child params follow python: left = (int(min), center),
    /// right = (center, int(max)), with the node's own (float) params.
    fn add_inner(&mut self, n: u32, b: u32, min_p: f64, max_p: f64) {
        let c = self.nodes[n as usize].center;
        let blk = self.blocks[b as usize];
        if (blk.end as f64) <= c {
            let (cmin, cmax) = (min_p.trunc(), c);
            let mut child = self.nodes[n as usize].left;
            if child == NIL {
                child = self.new_node(cmin, cmax);
                self.nodes[n as usize].left = child;
            }
            self.add_interval(child, b, cmin, cmax);
        } else if (blk.start as f64) > c {
            let (cmin, cmax) = (c, max_p.trunc());
            let mut child = self.nodes[n as usize].right;
            if child == NIL {
                child = self.new_node(cmin, cmax);
                self.nodes[n as usize].right = child;
            }
            self.add_interval(child, b, cmin, cmax);
        } else {
            let tail = self.nodes[n as usize].tail;
            if tail == NIL {
                self.nodes[n as usize].head = b;
            } else {
                self.next[tail as usize] = b;
            }
            self.nodes[n as usize].tail = b;
            self.next[b as usize] = NIL;
        }
    }
}

fn parse_i64(s: &[u8], what: &str, header: &str) -> Result<i64, String> {
    std::str::from_utf8(s)
        .ok()
        .and_then(|t| t.parse::<i64>().ok())
        .ok_or_else(|| format!("liftover: bad {what} in chain ({header})"))
}

fn fields(line: &[u8]) -> Vec<&[u8]> {
    line.split(|c| c.is_ascii_whitespace()).filter(|s| !s.is_empty()).collect()
}

impl LiftOver {
    pub fn open(path: &Path) -> Result<LiftOver, String> {
        let f = File::open(path).map_err(|e| format!("liftover: cannot open {}: {e}", path.display()))?;
        let gz = path.to_string_lossy().to_lowercase().ends_with(".gz");
        let r: Box<dyn Read> = if gz { Box::new(MultiGzDecoder::new(BufReader::with_capacity(1 << 20, f))) } else { Box::new(f) };
        LiftOver::from_reader(BufReader::with_capacity(1 << 20, r))
    }

    /// Parse an (already decompressed) chain stream.
    pub fn from_reader<R: BufRead>(mut r: R) -> Result<LiftOver, String> {
        let mut chains: Vec<Chain> = Vec::new();
        let mut blocks: Vec<Block> = Vec::new();
        let mut order: Vec<String> = Vec::new(); // source contigs in first-seen order
        let mut src_size: FxHashMap<String, i64> = FxHashMap::default();
        let mut tgt_size: FxHashMap<String, i64> = FxHashMap::default();
        let mut by_src: FxHashMap<String, Vec<u32>> = FxHashMap::default();
        let mut line: Vec<u8> = Vec::new();
        loop {
            line.clear();
            if r.read_until(b'\n', &mut line).map_err(|e| format!("liftover: read error: {e}"))? == 0 {
                break;
            }
            if !line.starts_with(b"chain") {
                continue; // '#', blank and any other line between chains
            }
            let header = String::from_utf8_lossy(&line).trim_end().to_string();
            let f = fields(&line);
            if f.len() < 12 {
                return Err(format!("liftover: invalid chain format ({header})"));
            }
            let score = parse_i64(f[1], "score", &header)?;
            let sname = String::from_utf8_lossy(f[2]).into_owned();
            let ssize = parse_i64(f[3], "source size", &header)?;
            if f[4] != b"+" {
                return Err(format!("liftover: source strand in an .over.chain file must be + ({header})"));
            }
            let sstart = parse_i64(f[5], "source start", &header)?;
            let send = parse_i64(f[6], "source end", &header)?;
            let tname = String::from_utf8_lossy(f[7]).into_owned();
            let tsize = parse_i64(f[8], "target size", &header)?;
            let tminus = match f[9] {
                b"+" => false,
                b"-" => true,
                _ => return Err(format!("liftover: target strand must be - or + ({header})")),
            };
            let tstart = parse_i64(f[10], "target start", &header)?;
            let tend = parse_i64(f[11], "target end", &header)?;

            let (mut sfrom, mut tfrom) = (sstart, tstart);
            let chain_id = chains.len() as u32;
            let mut mine: Vec<Block> = Vec::new();
            loop {
                line.clear();
                r.read_until(b'\n', &mut line).map_err(|e| format!("liftover: read error: {e}"))?;
                let f = fields(&line);
                if f.len() == 3 {
                    let size = parse_i64(f[0], "block size", &header)?;
                    let sgap = parse_i64(f[1], "block gap", &header)?;
                    let tgap = parse_i64(f[2], "block gap", &header)?;
                    mine.push(Block { start: sfrom, end: sfrom + size, tfrom, chain: chain_id });
                    sfrom += size + sgap;
                    tfrom += size + tgap;
                } else {
                    if f.len() != 1 {
                        return Err(format!("liftover: expecting one number on the last line of alignments block ({header})"));
                    }
                    let size = parse_i64(f[0], "block size", &header)?;
                    mine.push(Block { start: sfrom, end: sfrom + size, tfrom, chain: chain_id });
                    if sfrom + size != send || tfrom + size != tend {
                        return Err(format!("liftover: alignment blocks do not match specified block sizes ({header})"));
                    }
                    break;
                }
            }
            // size consistency (pyliftover `_index_chains`)
            let e = *src_size.entry(sname.clone()).or_insert(ssize);
            if e != ssize {
                return Err(format!("liftover: inconsistent source chromosome size for {sname} ({e} vs {ssize})"));
            }
            let e = *tgt_size.entry(tname.clone()).or_insert(tsize);
            if e != tsize {
                return Err(format!("liftover: inconsistent target chromosome size for {tname} ({e} vs {tsize})"));
            }
            chains.push(Chain { target_name: tname, target_size: tsize, target_minus: tminus, score });
            let list = by_src.entry(sname.clone()).or_insert_with(|| {
                order.push(sname.clone());
                Vec::new()
            });
            for b in mine {
                if b.end - b.start <= 0 {
                    continue; // zero-length blocks are not registered (the contig still gets a tree)
                }
                list.push(blocks.len() as u32);
                blocks.push(b);
            }
        }

        let mut bld = Builder { next: vec![NIL; blocks.len()], blocks, nodes: Vec::new() };
        let mut trees: FxHashMap<String, (u32, i64)> = FxHashMap::default();
        for name in &order {
            let size = src_size[name];
            let root = bld.new_node(0.0, size as f64);
            for &b in &by_src[name] {
                bld.add_interval(root, b, 0.0, size as f64);
            }
            trees.insert(name.clone(), (root, size));
        }

        // flatten the linked lists into the two sorted pools (stable sorts, python semantics:
        // both sorts start from the insertion-ordered list)
        let Builder { blocks, mut nodes, next } = bld;
        let mut pool_start: Vec<u32> = Vec::new();
        let mut pool_end: Vec<u32> = Vec::new();
        let mut tmp: Vec<u32> = Vec::new();
        for nd in nodes.iter_mut() {
            if nd.state != 2 {
                continue;
            }
            tmp.clear();
            let mut b = nd.head;
            while b != NIL {
                tmp.push(b);
                b = next[b as usize];
            }
            nd.off = pool_start.len() as u32;
            nd.len = tmp.len() as u32;
            tmp.sort_by_key(|&b| blocks[b as usize].start);
            pool_start.extend_from_slice(&tmp);
            tmp.clear();
            let mut b = nd.head;
            while b != NIL {
                tmp.push(b);
                b = next[b as usize];
            }
            tmp.sort_by_key(|&b| std::cmp::Reverse(blocks[b as usize].end));
            pool_end.extend_from_slice(&tmp);
        }
        Ok(LiftOver { chains, blocks, nodes, pool_start, pool_end, trees })
    }

    /// pyliftover `convert_coordinate(chrom, pos)` (strand '+'). `None` when `chrom` is not a
    /// source contig of any chain.
    pub fn convert_coordinate(&self, chrom: &str, pos: i64) -> Option<Vec<Lifted<'_>>> {
        let &(root, size) = self.trees.get(chrom)?;
        let mut hits: Vec<u32> = Vec::new();
        self.query(root, pos, 0.0, size as f64, &mut hits);
        let mut out: Vec<Lifted<'_>> = hits
            .iter()
            .map(|&b| {
                let blk = &self.blocks[b as usize];
                let ch = &self.chains[blk.chain as usize];
                let mut p = blk.tfrom + (pos - blk.start);
                if ch.target_minus {
                    p = ch.target_size - 1 - p;
                }
                Lifted { chrom: &ch.target_name, pos: p, strand: if ch.target_minus { b'-' } else { b'+' }, score: ch.score }
            })
            .collect();
        out.sort_by(|a, b| b.score.cmp(&a.score)); // stable, descending
        Some(out)
    }

    /// `IntervalTree._query`: result order = python's traversal order.
    fn query(&self, n: u32, x: i64, min_p: f64, max_p: f64, out: &mut Vec<u32>) {
        let nd = &self.nodes[n as usize];
        match nd.state {
            0 => {}
            1 => {
                let b = &self.blocks[nd.single as usize];
                if b.start <= x && x < b.end {
                    out.push(nd.single);
                }
            }
            _ => {
                let c = nd.center;
                if (x as f64) < c {
                    if nd.left != NIL {
                        self.query(nd.left, x, min_p.trunc(), c, out);
                    }
                    for &b in &self.pool_start[nd.off as usize..(nd.off + nd.len) as usize] {
                        if self.blocks[b as usize].start <= x {
                            out.push(b);
                        } else {
                            break;
                        }
                    }
                } else {
                    for &b in &self.pool_end[nd.off as usize..(nd.off + nd.len) as usize] {
                        if self.blocks[b as usize].end > x {
                            out.push(b);
                        } else {
                            break;
                        }
                    }
                    if nd.right != NIL {
                        self.query(nd.right, x, c, max_p.trunc(), out);
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const CHAIN: &str = "\
# comment
chain 1000 chrA 1000 + 10 200 chrB 500 - 0 190 1
50 10 0
60 20 30
50

chain 900 chrA 1000 + 0 100 chrC 1000 + 500 600 2
100

chain 900 chrA 1000 + 40 70 chrD 1000 + 0 30 3
30

chain 1 chrZ 100 + 0 0 chrZ 100 + 0 0 4
0
";

    #[test]
    fn small_chain() {
        let lo = LiftOver::from_reader(CHAIN.as_bytes()).unwrap();
        assert!(lo.convert_coordinate("chrQ", 5).is_none());
        assert_eq!(lo.convert_coordinate("chrZ", 5).unwrap().len(), 0);
        let r = lo.convert_coordinate("chrA", 45).unwrap();
        assert_eq!(r.len(), 3);
        assert_eq!(r[0], Lifted { chrom: "chrB", pos: 500 - 1 - 35, strand: b'-', score: 1000 });
        assert_eq!((r[1].chrom, r[1].pos), ("chrC", 545));
        assert_eq!((r[2].chrom, r[2].pos), ("chrD", 5));
        assert!(lo.convert_coordinate("chrA", 950).unwrap().is_empty());
    }

    /// Literals generated by pyliftover 0.4.1 (reference venv) from the same chain text.
    #[test]
    fn matches_pyliftover_literals() {
        let lo = LiftOver::from_reader(CHAIN.as_bytes()).unwrap();
        let cases: &[(i64, &[(&str, i64, u8, i64)])] = &[
            (0, &[("chrC", 500, b'+', 900)]),
            (9, &[("chrC", 509, b'+', 900)]),
            (10, &[("chrB", 499, b'-', 1000), ("chrC", 510, b'+', 900)]),
            (39, &[("chrB", 470, b'-', 1000), ("chrC", 539, b'+', 900)]),
            (40, &[("chrB", 469, b'-', 1000), ("chrC", 540, b'+', 900), ("chrD", 0, b'+', 900)]),
            (45, &[("chrB", 464, b'-', 1000), ("chrC", 545, b'+', 900), ("chrD", 5, b'+', 900)]),
            (59, &[("chrB", 450, b'-', 1000), ("chrC", 559, b'+', 900), ("chrD", 19, b'+', 900)]),
            (60, &[("chrC", 560, b'+', 900), ("chrD", 20, b'+', 900)]),
            (69, &[("chrC", 569, b'+', 900), ("chrD", 29, b'+', 900)]),
            (70, &[("chrB", 449, b'-', 1000), ("chrC", 570, b'+', 900)]),
            (99, &[("chrB", 420, b'-', 1000), ("chrC", 599, b'+', 900)]),
            (100, &[("chrB", 419, b'-', 1000)]),
            (129, &[("chrB", 390, b'-', 1000)]),
            (130, &[]),
            (149, &[]),
            (150, &[("chrB", 359, b'-', 1000)]),
            (179, &[("chrB", 330, b'-', 1000)]),
            (199, &[("chrB", 310, b'-', 1000)]),
            (200, &[]),
            (201, &[]),
        ];
        for (pos, exp) in cases {
            let got: Vec<(&str, i64, u8, i64)> =
                lo.convert_coordinate("chrA", *pos).unwrap().iter().map(|x| (x.chrom, x.pos, x.strand, x.score)).collect();
            assert_eq!(&got[..], *exp, "pos {pos}");
        }
    }

    /// Randomised parity vs pyliftover: `LIFTOVER_PARITY_DIR` holds `*.chain.gz` plus
    /// `*.chain.gz.expected` (JSON lines `[chrom, pos, null | [[chrom,pos,strand,score]..]]`
    /// dumped by pyliftover; generator: scratchpad `lo/gen.py`). Skipped when absent.
    #[test]
    fn matches_pyliftover_random() {
        let dir = std::env::var("LIFTOVER_PARITY_DIR").unwrap_or_else(|_| {
            "/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/262d6400-c197-4f82-b61f-06444781e46f/scratchpad/lo".into()
        });
        let Ok(rd) = std::fs::read_dir(&dir) else { return };
        let mut n_files = 0;
        for e in rd.flatten() {
            let p = e.path();
            if !p.to_string_lossy().ends_with(".chain.gz") {
                continue;
            }
            let exp = std::fs::read_to_string(format!("{}.expected", p.display())).unwrap();
            let lo = LiftOver::open(&p).unwrap();
            let mut n = 0;
            for l in exp.lines() {
                let v: serde_json::Value = serde_json::from_str(l).unwrap();
                let (c, pos) = (v[0].as_str().unwrap(), v[1].as_i64().unwrap());
                let got = lo.convert_coordinate(c, pos);
                match (&v[2], got) {
                    (serde_json::Value::Null, None) => {}
                    (serde_json::Value::Array(a), Some(g)) => {
                        assert_eq!(a.len(), g.len(), "{c}:{pos}");
                        for (x, y) in a.iter().zip(g.iter()) {
                            assert_eq!(x[0].as_str().unwrap(), y.chrom, "{c}:{pos}");
                            assert_eq!(x[1].as_i64().unwrap(), y.pos, "{c}:{pos}");
                            assert_eq!(x[2].as_str().unwrap().as_bytes()[0], y.strand, "{c}:{pos}");
                            assert_eq!(x[3].as_i64().unwrap(), y.score, "{c}:{pos}");
                        }
                    }
                    (e, g) => panic!("{c}:{pos}: expected {e:?}, got {g:?}"),
                }
                n += 1;
            }
            assert!(n > 1000);
            n_files += 1;
        }
        eprintln!("liftover parity: {n_files} chain files checked");
    }
}
