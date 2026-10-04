"""Minimal rooted-tree support for the phylogenetic evaluation (no external dependency).

Parses Newick as written by the Sanger SNV-tree pipelines (`patients/<organ>/<patient>/*.tree`):
tip labels, internal support labels (`)100:11389`), branch lengths, zero-length branches, an
optional root length. Every node except the root defines a BRANCH: an insertion that happened on
that branch is carried by exactly the tips below the node (its clade). The root stands for
"present in every colony" (germline, or somatic before the MRCA of all sampled colonies).
"""
from __future__ import annotations

import random
import re
from dataclasses import dataclass, field
from typing import Iterable, List, Optional, Sequence

import numpy as np


@dataclass
class Node:
    name: str = ""                 # tip label; internal: support label or "" (ids are assigned)
    length: float = 0.0
    children: List["Node"] = field(default_factory=list)
    parent: Optional["Node"] = None
    id: str = ""                   # tips: label; internal: N<k> (preorder); root: ROOT

    @property
    def is_tip(self) -> bool:
        return not self.children

    def tips(self) -> List["Node"]:
        if self.is_tip:
            return [self]
        out = []
        for c in self.children:
            out.extend(c.tips())
        return out

    def preorder(self) -> Iterable["Node"]:
        stack = [self]
        while stack:
            n = stack.pop()
            yield n
            stack.extend(reversed(n.children))


_TOKEN = re.compile(r"\s*(\(|\)|,|;|:|'[^']*'|[^(),;:\s]+)")


def parse_newick(text: str) -> Node:
    """Parse one Newick tree. Labels may be quoted; comments `[...]` are dropped."""
    text = re.sub(r"\[[^\]]*\]", "", text.strip())
    toks = [t for t in _TOKEN.findall(text) if t.strip()]
    pos = 0

    def peek():
        return toks[pos] if pos < len(toks) else ";"

    def node() -> Node:
        nonlocal pos
        n = Node()
        if peek() == "(":
            pos += 1
            while True:
                c = node()
                c.parent = n
                n.children.append(c)
                t = peek()
                pos += 1
                if t == ",":
                    continue
                if t == ")":
                    break
                raise ValueError(f"newick: unexpected token {t!r}")
        t = peek()
        if t not in ("(", ")", ",", ";", ":"):
            n.name = t.strip("'")
            pos += 1
        if peek() == ":":
            pos += 1
            n.length = float(toks[pos])
            pos += 1
        return n

    root = node()
    assign_ids(root)
    return root


def read_newick(path: str) -> Node:
    with open(path) as fh:
        return parse_newick(fh.read())


def assign_ids(root: Node) -> None:
    k = 0
    for n in root.preorder():
        if n is root:
            n.id = "ROOT"
        elif n.is_tip:
            n.id = n.name
        else:
            k += 1
            n.id = f"N{k}"


def to_newick(root: Node, digits: int = 6) -> str:
    def rec(n: Node) -> str:
        s = ""
        if n.children:
            s = "(" + ",".join(rec(c) for c in n.children) + ")"
        if n.is_tip:
            s += n.name
        if n.parent is not None:
            s += f":{round(n.length, digits):g}"
        return s
    return rec(root) + ";"


def prune(root: Node, keep: Sequence[str]) -> Node:
    """Copy of the tree restricted to tips in `keep`; unary nodes are collapsed (lengths summed),
    so branch lengths keep their mutation-time meaning. Raises if fewer than 1 tip remains."""
    keep = set(keep)

    def rec(n: Node) -> Optional[Node]:
        if n.is_tip:
            return Node(name=n.name, length=n.length) if n.name in keep else None
        kids = [k for k in (rec(c) for c in n.children) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            kids[0].length += n.length
            return kids[0]
        m = Node(name=n.name, length=n.length, children=kids)
        for k in kids:
            k.parent = m
        return m

    r = rec(root)
    if r is None:
        raise ValueError("prune: no tips left")
    r.parent = None
    if r.is_tip:                       # a single colony: wrap so the root is still "all"
        r = Node(children=[r])
        r.children[0].parent = r
    assign_ids(r)
    return r


@dataclass
class Branches:
    """Branch enumeration of a rooted tree over a fixed colony order.

    ids[b], lengths[b], members[b] (tips below), mask (B x C bool): mask[b, c] = colony c is
    in branch b's clade. Branch 0 is the ROOT (all colonies); tips and internal nodes follow
    in preorder."""
    colonies: List[str]
    ids: List[str]
    lengths: np.ndarray
    is_tip: np.ndarray
    mask: np.ndarray
    members: List[List[str]]

    @property
    def root_index(self) -> int:
        return 0

    def index(self, node_id: str) -> int:
        return self.ids.index(node_id)

    def clade_index(self, colonies: Iterable[str]) -> Optional[int]:
        s = set(colonies)
        for b, m in enumerate(self.members):
            if set(m) == s:
                return b
        return None


def branches(root: Node, colonies: Optional[Sequence[str]] = None) -> Branches:
    tips = [t.name for t in root.tips()]
    if colonies is None:
        colonies = tips
    colonies = list(colonies)
    if set(colonies) != set(tips):
        raise ValueError("branches: colony list must equal the tree's tip set (prune first)")
    col_ix = {c: i for i, c in enumerate(colonies)}
    ids, lens, tipf, rows, members = [], [], [], [], []
    for n in root.preorder():
        mem = [t.name for t in n.tips()]
        row = np.zeros(len(colonies), dtype=bool)
        row[[col_ix[t] for t in mem]] = True
        ids.append(n.id)
        lens.append(n.length if n is not root else 0.0)
        tipf.append(n.is_tip and n is not root)
        rows.append(row)
        members.append(sorted(mem))
    return Branches(colonies, ids, np.array(lens, float), np.array(tipf, bool),
                    np.vstack(rows), members)


def random_coalescent(n: int, rng: random.Random, names: Optional[Sequence[str]] = None,
                      height: float = 1000.0) -> Node:
    """Kingman coalescent with n tips, scaled so the root height is `height` (SNV-like units).
    Deterministic for a given rng state."""
    names = list(names) if names else [f"S{i + 1}" for i in range(n)]
    nodes = [Node(name=nm) for nm in names]
    heights = {id(x): 0.0 for x in nodes}
    t = 0.0
    while len(nodes) > 1:
        k = len(nodes)
        t += rng.expovariate(k * (k - 1) / 2.0)
        i, j = sorted(rng.sample(range(k), 2))
        a, b = nodes[i], nodes[j]
        p = Node(children=[a, b])
        a.parent = b.parent = p
        a.length = t - heights[id(a)]
        b.length = t - heights[id(b)]
        heights[id(p)] = t
        nodes = [x for q, x in enumerate(nodes) if q not in (i, j)] + [p]
    root = nodes[0]
    scale = height / t if t > 0 else 1.0
    for x in root.preorder():
        x.length *= scale
    assign_ids(root)
    return root
