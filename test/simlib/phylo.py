"""Tree-structured presence for the multi-colony simulators (`build_donor.py --types ... --tree`).

A real somatic insertion is inherited by exactly the colonies below the branch it happened on;
a library / mapping artefact (or a mis-merged call) scatters over the tree. This module draws:

  * the tree: a Newick file (tips are mapped to samples S1..Sn in leaf order; the written
    `tree.nwk` uses the sample names) or `random:N` (Kingman coalescent, tools/phylo/tree.py);
  * per-sample colony purity (U(purity range): fraction of the colony's cells descending from
    its founder; contaminating cells carry none of its somatic events) and a depth factor
    (a fraction of colonies is low-depth, to create dropouts);
  * per event, its carriers: a branch (root = every colony with probability root_frac, else a
    non-root branch, uniform or proportional to length), or for a `nonclade` decoy a random
    subset that is NOT a clade (>= 2 carriers, >= 1 non-carrier).

All draws use the caller's `random.Random`, so a run is reproducible from its seed; nothing here
is touched when --tree is not given (the default simulator path is unchanged).
"""
from __future__ import annotations

import os
import random
import sys
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Optional

_REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)
from tools.phylo import tree as T  # noqa: E402


@dataclass
class TreeDesign:
    root: T.Node
    samples: List[str]                      # S1..Sn, index i -> sample i+1
    purity: List[float]
    depth_factor: List[float]
    clades: List[FrozenSet[int]] = field(default_factory=list)   # non-root branches (0-based samples)
    lengths: List[float] = field(default_factory=list)
    ids: List[str] = field(default_factory=list)

    @property
    def n(self) -> int:
        return len(self.samples)

    def clade_set(self) -> set:
        return set(self.clades)


def build_design(spec: str, rng: random.Random, purity_range=(0.7, 1.0), low_depth_frac=0.25,
                 low_depth_range=(0.25, 0.5)) -> TreeDesign:
    if spec.startswith("random:"):
        n = int(spec.split(":", 1)[1])
        root = T.random_coalescent(n, rng)
    else:
        root = T.read_newick(spec)
        for i, t in enumerate(root.tips()):        # map tips to S1..Sn in leaf order
            t.name = f"S{i + 1}"
        T.assign_ids(root)
    tips = [t.name for t in root.tips()]
    samples = sorted(tips, key=lambda s: int(s[1:]))
    purity = [rng.uniform(*purity_range) for _ in samples]
    dfac = [rng.uniform(*low_depth_range) if rng.random() < low_depth_frac else rng.uniform(0.8, 1.2)
            for _ in samples]
    br = T.branches(root, samples)
    clades, lens, ids = [], [], []
    for b in range(1, len(br.ids)):
        clades.append(frozenset(int(m[1:]) - 1 for m in br.members[b]))
        lens.append(float(br.lengths[b]))
        ids.append(br.ids[b])
    return TreeDesign(root, samples, purity, dfac, clades, lens, ids)


def draw_branch(rng: random.Random, d: TreeDesign, root_frac: float = 0.1,
                weight: str = "uniform") -> (str, FrozenSet[int]):
    """(branch id, carriers) for a TP event."""
    if rng.random() < root_frac:
        return "ROOT", frozenset(range(d.n))
    if weight == "length":
        w = [max(x, 1e-9) for x in d.lengths]
        k = rng.choices(range(len(d.clades)), weights=w)[0]
    else:
        k = rng.randrange(len(d.clades))
    return d.ids[k], d.clades[k]


def draw_nonclade(rng: random.Random, d: TreeDesign, tries: int = 2000) -> Optional[FrozenSet[int]]:
    if d.n < 3:
        return None
    cs = d.clade_set()
    for _ in range(tries):
        k = rng.randint(2, d.n - 1)
        s = frozenset(rng.sample(range(d.n), k))
        if s not in cs:
            return s
    return None


def write_design(d: TreeDesign, out_dir: str, placements: Dict[int, tuple]) -> None:
    """tree.nwk, samples.tsv (sample, purity, depth_factor), phylo_truth.tsv (event id,
    placement = ROOT / branch id / NONCLADE / ARTEFACT, carriers as 1-based sample list)."""
    with open(os.path.join(out_dir, "tree.nwk"), "w") as fh:
        fh.write(T.to_newick(d.root) + "\n")
    with open(os.path.join(out_dir, "samples.tsv"), "w") as fh:
        fh.write("sample\tname\tpurity\tdepth_factor\n")
        for i, s in enumerate(d.samples):
            fh.write(f"{i + 1}\t{s}\t{d.purity[i]:.4f}\t{d.depth_factor[i]:.3f}\n")
    with open(os.path.join(out_dir, "phylo_truth.tsv"), "w") as fh:
        fh.write("id\tplacement\tcarriers\n")
        for eid in sorted(placements):
            pl, car = placements[eid]
            fh.write(f"{eid}\t{pl}\t{','.join(str(i + 1) for i in sorted(car))}\n")
