#!/usr/bin/env python3
"""Apply the copy-number rule of cluster/somatic_table.py to an EXISTING <P>.somatic.xlsx, in place.

Rows of the `somatic` sheet whose carriers are exactly the tips outside a clade with a loss / LOH of
the locus (patients/<organ>/<P>/copy_number.tsv) move to the `excluded` sheet with the same reason
somatic_table.py gives; nothing else in the workbook changes. For tables built before the rule
existed, where a full rebuild is not possible (genotype / fit / insertions inputs left on the farm)
or would mix in other rule changes. Rerunnable: the README copy-number lines are replaced.

  python cluster/apply_copy_number.py --xlsx ~/Documents/hsc_results/PD49229/PD49229.somatic.xlsx \
      --joint ~/Documents/hsc_results/PD49229/PD49229.joint.tsv --copy-number patients/.../copy_number.tsv \
      --tree patients/.../PD49229.tree --assembly GRCh37

Needs openpyxl (reading).
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from somatic_table import (CN_TOLERANCE, cn_masked, read_centromeres, read_copy_number,  # noqa: E402
                           read_tsv, write_xlsx)

CN_REASON = "germline masked by copy-number loss"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--joint", required=True, help="<P>.joint.tsv (the colony set = its p_ columns)")
    ap.add_argument("--copy-number", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--assembly", required=True, choices=("GRCh38", "GRCh37"))
    a = ap.parse_args()

    import openpyxl
    first = next(iter(read_tsv(a.joint)), None)
    if first is None:
        sys.exit(f"no rows in {a.joint}")
    colonies = [k[2:] for k in first if k.startswith("p_")]
    cen = read_centromeres(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "resources",
                                        "centromeres.tsv"), a.assembly)
    events = read_copy_number(a.copy_number, a.tree, colonies, cen)

    wb = openpyxl.load_workbook(a.xlsx, read_only=False, data_only=True)
    sheets = []
    for ws in wb.worksheets:
        rows = [list(r) for r in ws.iter_rows(values_only=True)]
        widths = [ws.column_dimensions[openpyxl.utils.get_column_letter(i + 1)].width or 10
                  for i in range(max((len(r) for r in rows), default=0))]
        sheets.append([ws.title, rows, widths, ws.title != "README"])
    by = {s[0]: s for s in sheets}
    if "excluded" not in by:          # tables from before the hard rules: add the sheet after `somatic`
        som_at = next(i for i, s in enumerate(sheets) if s[0] == "somatic")
        sheets.insert(som_at + 1, ["excluded", [["excluded_by"] + list(by["somatic"][1][0])],
                                   [60] + by["somatic"][2], True])
        by = {s[0]: s for s in sheets}
    som, exc = by["somatic"][1], by["excluded"][1]
    H = {h: i for i, h in enumerate(som[0])}

    def carriers(r):
        for col in ("carriers", "treefit_legacy_carriers", "treefit_v2_carriers"):
            v = r[H[col]] if col in H else None
            if v:
                return [c for c in str(v).split(",") if c]
        return []

    keep, moved = [som[0]], 0
    for r in som[1:]:
        cn = cn_masked(str(r[H["locus"]]), carriers(r), colonies, events) if events and r[H["tier"]] != "A" else ""
        if cn:
            exc.append([f"{CN_REASON}: carriers = every tip outside {cn}"] + r)
            moved += 1
        else:
            keep.append(r)
    by["somatic"][1] = keep
    rank = {"A": 0, "B": 1, "C": 2, "D": 3}
    loc_i = 1 + H["locus"]
    exc[1:] = sorted(exc[1:], key=lambda r: (rank.get(r[1], 9), str(r[loc_i])))

    n_cn = sum(1 for r in exc[1:] if str(r[0] or "").startswith(CN_REASON))
    readme = [r for r in by["README"][1]
              if not (r and (r[0] == "copy number" or (not r[0] and len(r) > 1 and " = " in str(r[1] or "")
                                                           and str(r[1]).endswith(f"({a.assembly})"))))]
    at = next((i for i, r in enumerate(readme) if r and r[0] == "known insertions"), len(readme) - 1) + 1
    readme[at:at] = [["copy number", f"{os.path.abspath(a.copy_number)}  (applied afterwards by "
                                     f"cluster/apply_copy_number.py; {n_cn} rows moved to excluded, +-"
                                     f"{int(CN_TOLERANCE * 100)} % of tips)"]] + \
                    [["", f"{ev[4]} = {ev[0]}:{ev[1]}-{ev[2]} ({a.assembly})"] for ev in events]
    by["README"][1] = readme
    write_xlsx(a.xlsx + ".tmp", [tuple(s) for s in sheets])
    os.replace(a.xlsx + ".tmp", a.xlsx)
    print(f"{a.xlsx}: {len(events)} events, {moved} rows moved to excluded, {len(keep) - 1} somatic rows left")


if __name__ == "__main__":
    main()
