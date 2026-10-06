#!/usr/bin/env python3
"""Score a joint-step table (peartree-genotype2 --step joint) or a tree_fit.py phylo_fit.tsv
against the simulator truth (test/e2e/phylo_truth.py output).

Per truth class: how often the inferred carrier set equals the truth carriers (clade/private;
germline = every tip), and how often a nonclade/none locus is NOT called a clean tree event.

usage: score_joint.py --truth phylo/truth.tsv --joint P.joint.tsv [--python-fit phylo/fit/phylo_fit.tsv]
                      --tips S1,...,S10 [--bf-threshold 1.0]

A Rust-joint private event (one carrier) counts as a clean tree event when the carrier's
P(carrier) >= 0.9 (its BF saturates near log10 1.2-2.2, see the crate README); clades and the
Python fit keep the BF threshold.
"""
import argparse
import collections
import csv

P_CARRIER = 0.9   # tools/genotype2_io.py carrier threshold


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def parse(row, tips, bf_thr):
    """(best id, is_clean_tree_event, carrier set) from a joint.tsv or phylo_fit.tsv row."""
    if "best" in row:                       # Rust joint.tsv
        best = row["best"]
        bf = float(row.get("log10_bf_tree", "nan"))
        carr = set(filter(None, row.get("carriers", "").split(",")))
    else:                                   # Python phylo_fit.tsv
        best = row["best_branch"]
        bf = float(row["log10_bf_tree"])
        if row["class"] == "noise":
            best = "NOISE"
        carr = set(filter(None, row["best_clade"].split(","))) if best != "ROOT" else set(tips)
    if best == "ROOT":
        carr = set(tips)
    tree = best not in ("NOISE", "INDEP") and best != "ROOT" and bf >= bf_thr
    if "best" in row and len(carr) == 1 and best not in ("NOISE", "INDEP", "ROOT"):
        # a private event's BF saturates (INDEP explains one carrier up to a combinatorial factor:
        # log10 BF <= ~1.2 at 10 colonies, ~2.2 at 44) and so does post_best (~0.85): judge it on
        # the carrier's P(carrier), the number the matrix / annotate_v2 threshold at 0.9
        (c,) = tuple(carr)
        try:
            tree = float(row.get("p_" + c, "") or "nan") >= P_CARRIER
        except ValueError:
            tree = False
    return best, tree, carr


def score(rows, truth, tips, bf_thr, label):
    by = {r["locus"]: r for r in rows}
    tips = set(tips)
    classes = ("germline", "clade", "private", "nonclade", "none")
    stats = collections.OrderedDict((k, collections.Counter()) for k in classes)
    for t in truth:
        cls = t["truth"]
        if cls not in stats:
            continue
        r = by.get(t["locus"])
        if r is None:
            stats[cls]["missing"] += 1
            continue
        best, tree, carr = parse(r, tips, bf_thr)
        want = tips if cls == "germline" else set(filter(None, t["carriers"].split(",")))
        stats[cls]["n"] += 1
        if cls == "germline":
            stats[cls]["ok"] += int(best == "ROOT")
        elif cls in ("clade", "private"):
            stats[cls]["carriers_ok"] += int(carr == want and best not in ("NOISE", "INDEP"))
            if tree and carr == want:
                stats[cls]["ok"] += 1
            elif best == "ROOT":
                stats[cls]["->root"] += 1
            elif best in ("NOISE", "INDEP"):
                stats[cls]["->" + best.lower()] += 1
            elif not tree:
                stats[cls]["low_bf"] += 1
            elif carr < want:
                stats[cls]["dropout"] += 1
            elif carr > want:
                stats[cls]["extra"] += 1
            else:
                stats[cls]["wrong_branch"] += 1
        elif cls == "nonclade":
            stats[cls]["ok"] += int(not tree)
        else:  # none: no simulated event -- report what the fit says
            stats[cls]["root" if best == "ROOT" else best.lower() if best in ("NOISE", "INDEP") else ("branch" if tree else "branch_lowbf")] += 1
    print(f"### {label}")
    print("| truth class | n | correct | fraction | breakdown |")
    print("|---|---|---|---|---|")
    for k, c in stats.items():
        n = c.pop("n", 0); ok = c.pop("ok", 0)
        rest = ", ".join(f"{kk} {vv}" for kk, vv in sorted(c.items()))
        frac = f"{ok / n:.3f}" if n and k != "none" else "."
        print(f"| {k} | {n} | {ok if k != 'none' else '.'} | {frac} | {rest} |")
    print("(germline = best is ROOT; clade/private = clean tree event with exactly the truth carriers;")
    print(" nonclade = not a clean tree event; none = no simulated event, shown as the fit's verdict)")
    print()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--joint")
    ap.add_argument("--python-fit")
    ap.add_argument("--tips", required=True)
    ap.add_argument("--bf-threshold", type=float, default=1.0)
    ap.add_argument("--label", default="")
    a = ap.parse_args()
    tips = a.tips.split(",")
    truth = read_tsv(a.truth)
    if a.joint:
        score(read_tsv(a.joint), truth, tips, a.bf_threshold, f"Rust joint {a.label}")
    if a.python_fit:
        score(read_tsv(a.python_fit), truth, tips, a.bf_threshold, f"Python tree_fit {a.label}")


if __name__ == "__main__":
    main()
