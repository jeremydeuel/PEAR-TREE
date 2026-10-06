"""Readers for genotyping output, legacy (call strings) and peartree-genotype2 (numeric).

One place for the file formats and the carrier / wild-type thresholds every downstream Python
consumer uses (tools/annotate_v2.py, tools/phylo/tree_fit.py, cluster/tprt/*.py,
test/e2e/score_genotypes.py). The format is always auto-detected from the header / content.

Formats
-------
Per-colony file `<colony>.txt.gz` (TAB-separated):
  * v2      (rust/peartree-genotype2, `types::OUTPUT_HEADER`):
            locus kind status depth n_alt n_ref n_uninf n_art n_disc n_alt_l n_alt_r vaf
            p_absent p_het p_hom pl_absent pl_het pl_hom gq score_alt score_ref
            status in ok / no_reads / high_coverage / error; non-ok rows have empty numeric
            cells beyond `depth`.
  * legacy  (peartree-genotype): insertion genotype score_genotype score_alternative coverage
            n_alt n_ref n_art   (older files lack the count columns)
Patient matrix `<P>.genotypes.csv.gz` (`;`-separated, rows = loci, columns = colonies):
  * numeric (genotype2 joint step): header `;<colony1>;...` (EMPTY first cell), cells
            P(carrier) with 4 decimals, EMPTY = colony without data
  * calls   (src/combine_genotypes.py): header `insertion;<colony1>;...`, cells are call strings
            (heterozygous / homozygous / insertion / wild-type / wild-type? / artefact / ...)
Joint table `<P>.joint.tsv` (genotype2 joint step, next to the numeric matrix): locus locus_kind
  n_colonies_data best carriers n_carriers post_best best_tree best_tree_clade post_tree
  log10_bf_tree log10_L_best_tree log10_L_root log10_L_indep log10_L_noise p_<colony>...

Buckets (`present` / `absent` / `ambiguous`, plus a no-data bucket) are what the consumers
compare; no consumer writes call strings for v2 input.

The stdlib part (detection, row iterators, buckets, the joint table) has no third-party
imports, so the farm's system python can use it; the DataFrame readers import pandas lazily.
"""
from __future__ import annotations

import gzip
import os
import re

# ------------------------------------------------------------------ thresholds (the ONE place)
P_CARRIER = 0.9          # carrier: per-colony p_het + p_hom >= P_CARRIER, matrix P(carrier) >= P_CARRIER
P_ABSENT_ROW = 0.8       # confidently absent: per-colony p_absent >= P_ABSENT_ROW
P_ABSENT_MATRIX = 0.1    # wild-type in the matrix: P(carrier) <= P_ABSENT_MATRIX (= 1 - P_CARRIER)

# ------------------------------------------------------------------ format vocabulary
V2_HEADER = ("locus", "kind", "status", "depth", "n_alt", "n_ref", "n_uninf", "n_art", "n_disc",
             "n_alt_l", "n_alt_r", "vaf", "p_absent", "p_het", "p_hom", "pl_absent", "pl_het",
             "pl_hom", "gq", "score_alt", "score_ref")
V2_NUMERIC = V2_HEADER[3:]
V2_STATUS_OK = "ok"
V2_STATUSES = ("ok", "no_reads", "high_coverage", "error")
LEGACY_HEADER = ("insertion", "genotype", "score_genotype", "score_alternative", "coverage",
                 "n_alt", "n_ref", "n_art")

FMT_V2, FMT_LEGACY = "v2", "legacy"                 # per-colony files
FMT_NUMERIC, FMT_CALLS = "numeric", "calls"         # patient matrix

PRESENT, ABSENT, AMBIGUOUS, NO_DATA = "present", "absent", "ambiguous", "no_data"

LEGACY_CARRIER = frozenset(("heterozygous", "homozygous", "insertion"))
LEGACY_WILDTYPE = frozenset(("wild-type",))
LEGACY_NO_DATA = frozenset(("", "NA", "nan", "None", "no-coverage", "high-coverage", "error"))

_FILE_RE = re.compile(r"\.(txt|tsv)(\.gz)?$")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


# ------------------------------------------------------------------ detection


def colony_header_format(header):
    """Format of a per-colony file from its header fields: FMT_V2, FMT_LEGACY or None."""
    h = [x.strip() for x in header]
    if h[:3] == list(V2_HEADER[:3]):
        return FMT_V2
    if h[:2] == list(LEGACY_HEADER[:2]):
        return FMT_LEGACY
    return None


def colony_format(path):
    """FMT_V2 / FMT_LEGACY of a per-colony genotype file (raises ValueError if neither)."""
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    fmt = colony_header_format(header)
    if fmt is None:
        raise ValueError(f"{path}: not a genotype file (header {header[:4]})")
    return fmt


def is_v2_colony_file(path):
    try:
        return colony_format(path) == FMT_V2
    except (OSError, ValueError, EOFError):
        return False


def cells_numeric(cells):
    """True when every non-empty cell parses as a probability in [0, 1] and at least one does
    (a row of the genotype2 numeric matrix). A row of call strings is False."""
    seen = False
    for c in cells:
        c = c.strip()
        if c == "":
            continue
        try:
            x = float(c)
        except ValueError:
            return False
        if not 0.0 <= x <= 1.0:
            return False
        seen = True
    return seen


def matrix_format(path, max_rows=200):
    """FMT_NUMERIC / FMT_CALLS of a `;` patient matrix. Decided from the first rows with data;
    a matrix without any data cell falls back to the header (numeric has an EMPTY first cell,
    combine_genotypes writes `insertion`)."""
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split(";")
        for i, line in enumerate(fh):
            if i >= max_rows:
                break
            cells = line.rstrip("\n").split(";")[1:]
            if any(c.strip() for c in cells):
                return FMT_NUMERIC if cells_numeric(cells) else FMT_CALLS
    return FMT_NUMERIC if header[0].strip() == "" else FMT_CALLS


def joint_tsv_for(matrix_path):
    """`<P>.joint.tsv` next to `<P>.genotypes.csv.gz` (pipeline phase 4) or
    `<P>.joint_matrix.csv.gz` (the joint step's own default name); None when absent."""
    if not matrix_path:
        return None
    d, base = os.path.split(str(matrix_path))
    for suf in (".genotypes.csv.gz", ".genotypes.csv", ".joint_matrix.csv.gz", ".joint_matrix.csv",
                ".csv.gz", ".csv"):
        if base.endswith(suf):
            cand = os.path.join(d, base[:-len(suf)] + ".joint.tsv")
            return cand if os.path.exists(cand) else None
    return None


# ------------------------------------------------------------------ buckets


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def v2_p_present(row):
    """p_het + p_hom of a v2 per-colony row (dict); NaN when the row has no model result."""
    return _f(row.get("p_het")) + _f(row.get("p_hom"))


def v2_row_bucket(row):
    """Bucket of one v2 per-colony row (dict of strings or numbers): the status word when
    status != ok, else present (p_het + p_hom >= P_CARRIER) / absent (p_absent >= P_ABSENT_ROW)
    / ambiguous."""
    st = str(row.get("status", "") or "")
    if st != V2_STATUS_OK:
        return st or NO_DATA
    if v2_p_present(row) >= P_CARRIER:
        return PRESENT
    if _f(row.get("p_absent")) >= P_ABSENT_ROW:
        return ABSENT
    return AMBIGUOUS


def p_carrier_bucket(p):
    """Bucket of one numeric matrix cell (str / float / None): present at P >= P_CARRIER,
    absent at P <= P_ABSENT_MATRIX, ambiguous between, NO_DATA for an empty cell / NaN."""
    x = _f(p)
    if x != x:                       # NaN / empty
        return NO_DATA
    if x >= P_CARRIER:
        return PRESENT
    if x <= P_ABSENT_MATRIX:
        return ABSENT
    return AMBIGUOUS


def legacy_bucket(call, carrier=LEGACY_CARRIER, wildtype=LEGACY_WILDTYPE):
    """Bucket of a legacy call string: present (carrier calls), absent (wild-type), NO_DATA
    (no-/high-coverage, empty), else ambiguous (wild-type?, insertion?, artefact). `carrier` /
    `wildtype` override the vocabularies (cluster/tprt keeps `wild-type?` as wild-type)."""
    c = "" if call is None else str(call).strip()
    if c in carrier:
        return PRESENT
    if c in wildtype:
        return ABSENT
    if c in LEGACY_NO_DATA:
        return NO_DATA
    return AMBIGUOUS


def matrix_cell_bucket(cell, fmt, carrier=LEGACY_CARRIER, wildtype=LEGACY_WILDTYPE):
    """Bucket of one patient-matrix cell of either format."""
    return p_carrier_bucket(cell) if fmt == FMT_NUMERIC else legacy_bucket(cell, carrier, wildtype)


def colony_row_bucket(row, fmt):
    """Bucket of one per-colony row (dict) of either format (legacy: its `genotype` call)."""
    return v2_row_bucket(row) if fmt == FMT_V2 else legacy_bucket(row.get("genotype", ""))


# ------------------------------------------------------------------ stdlib readers


def iter_colony_rows(path):
    """(format, iterator of row dicts keyed by the header) of a per-colony file. The locus name
    is under `locus` for BOTH formats (legacy `insertion` is copied there)."""
    fh = open_text(path)
    header = fh.readline().rstrip("\n").split("\t")
    fmt = colony_header_format(header)
    if fmt is None:
        fh.close()
        raise ValueError(f"{path}: not a genotype file (header {header[:4]})")

    def gen():
        with fh:
            for line in fh:
                if not line.strip():
                    continue
                r = dict(zip(header, line.rstrip("\n").split("\t")))
                if fmt == FMT_LEGACY:
                    r["locus"] = r.get("insertion", "")
                yield r
    return fmt, gen()


def read_colony_rows(path):
    """(format, {locus: row dict}) of a per-colony file (first row per locus wins)."""
    fmt, rows = iter_colony_rows(path)
    out = {}
    for r in rows:
        out.setdefault(r["locus"], r)
    return fmt, out


def read_matrix(path):
    """(format, colonies, {locus: [cells as str]}) of a `;` patient matrix of either format;
    (None, [], None) when the file does not exist."""
    if not path or not os.path.exists(path):
        return None, [], None
    fmt = matrix_format(path)
    rows = {}
    with open_text(path) as fh:
        colonies = fh.readline().rstrip("\n").split(";")[1:]
        for line in fh:
            f = line.rstrip("\n").split(";")
            if not f or not f[0]:
                continue
            rows[f[0]] = f[1:]
    return fmt, colonies, rows


def read_joint(path):
    """{locus: row dict} of a genotype2 `<P>.joint.tsv` ({} when path is None / missing)."""
    if not path or not os.path.exists(path):
        return {}
    out = {}
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if line.strip():
                r = dict(zip(header, line.rstrip("\n").split("\t")))
                out.setdefault(r.get("locus", ""), r)
    return out


def joint_class(row):
    """Coarse class of a joint-table row: ROOT / clade (a branch with >= 2 carriers) / private
    (a tip) / INDEP / NOISE ('NA' for an empty row)."""
    best = (row or {}).get("best", "")
    if best in ("ROOT", "INDEP", "NOISE", ""):
        return best or "NA"
    try:
        n = int(float(row.get("n_carriers", "0") or 0))
    except ValueError:
        n = 0
    return "clade" if n >= 2 else "private"


def list_colony_files(directory):
    """{colony: path} of the per-colony genotype files (`<colony>.txt[.gz]` / `.tsv[.gz]`) in a
    directory; logs, temporaries and sub-directories are skipped."""
    out = {}
    for f in sorted(os.listdir(directory)):
        p = os.path.join(directory, f)
        if ".tmp" in f or f.endswith((".log", ".md", ".json")) or not _FILE_RE.search(f) \
                or not os.path.isfile(p):
            continue
        out[colony_stem(f)] = p
    return out


def colony_stem(path):
    """Colony name of a per-colony genotype file (`genotypes/<SAMPLE>.txt.gz`)."""
    base = os.path.basename(str(path))
    return re.sub(r"(\.genotypes?)?(\.(txt|csv|tsv))?(\.gz)?$", "", base)


# ------------------------------------------------------------------ pandas readers (lazy import)


def read_colony_df(path):
    """DataFrame of one per-colony file, indexed by locus (index name `locus`), duplicates
    dropped (first wins), `df.attrs['format']` = FMT_V2 / FMT_LEGACY. v2: numeric columns are
    float (NaN on non-ok rows), `status` / `kind` strings. Legacy: as written, `insertion`
    renamed to `locus`."""
    import pandas as pd
    fmt = colony_format(path)
    if fmt == FMT_V2:
        df = pd.read_csv(path, sep="\t", dtype={"locus": str, "kind": str, "status": str},
                         compression="infer", keep_default_na=False, na_values=[""])
        for c in V2_NUMERIC:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")
    else:
        df = pd.read_csv(path, sep="\t", dtype={"insertion": str}, compression="infer")
        df = df.rename(columns={"insertion": "locus"})
    df = df.drop_duplicates("locus").set_index("locus")
    df.attrs["format"] = fmt
    return df


def read_matrix_df(path):
    """DataFrame of a `;` patient matrix (rows loci, columns colonies), `df.attrs['format']` =
    FMT_NUMERIC (float, NaN = no data) or FMT_CALLS (str)."""
    import pandas as pd
    fmt = matrix_format(path)
    df = pd.read_csv(path, sep=";", index_col=0, compression="infer", dtype=str)
    if fmt == FMT_NUMERIC:
        df = df.apply(pd.to_numeric, errors="coerce").astype(float)
    df.index = df.index.astype(str)
    df.index.name = "locus"
    df.attrs["format"] = fmt
    return df


def read_joint_df(path):
    """DataFrame of `<P>.joint.tsv` indexed by locus (None when path is None / missing)."""
    if not path or not os.path.exists(path):
        return None
    import pandas as pd
    df = pd.read_csv(path, sep="\t", dtype={"locus": str, "best": str, "carriers": str,
                                            "best_tree": str, "best_tree_clade": str},
                     keep_default_na=False, na_values=[""])
    return df.drop_duplicates("locus").set_index("locus")
