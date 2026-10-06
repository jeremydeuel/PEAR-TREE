#!/usr/bin/env python3
"""Compare peartree-genotype2 (numeric per-colony rows + joint step) with the legacy genotyper's
calls on the SAME contract, colonies and BAMs — one patient. Standard library only (runs on the
farm's system python3, no venv).

    genotype2_compare.py --patient PD37590 \
        --legacy-dir  <rundir>/genotypes            # legacy <colony>.txt.gz (string calls)
        --v2-dir      <v2>/genotypes                # genotype2 <colony>.txt.gz (numeric rows)
        --joint       <v2>/joint/<P>.joint.tsv      # genotype2 --step joint (length branch prior)
        [--joint-uniform <v2>/joint/<P>.joint.uniform.tsv]
        --tree        patients/<organ>/<P>/<P>.tree
        [--legacy-calls <rundir>/<P>.genotypes.csv.gz]   # combine_genotypes output (after its gates)
        [--legacy-fit  <eval>/<arm>/fit/phylo_fit.tsv]   # tools/phylo/tree_fit.py on the legacy calls
        [--known patients/<organ>/<P>/known_insertions.tsv]
        [--v2-logs <v2>/logs] [--legacy-logs <rundir>/logs]
        --out-dir <v2>/report

Per (locus, colony) the legacy call string is set against the genotype2 posterior, bucketed the
way annotate / score_genotypes read the numeric rows: present when p_het + p_hom >= P_PRESENT,
absent when p_absent >= P_ABSENT, else ambiguous (status != ok is its own bucket).

Writes into --out-dir: report.md, cells.tsv.gz (every locus x colony), per_locus.tsv (carrier
counts per method + the joint step's hypothesis), discordant.tsv (legacy carrier vs v2 absent and
legacy wild-type vs v2 present, with the read counts), known.tsv (when --known is given).
"""
import argparse
import glob
import gzip
import math
import os
import re
import statistics
import sys
from collections import Counter, defaultdict

P_PRESENT = 0.9      # test/e2e/score_genotypes.py V2_P_PRESENT, tools/annotate_v2.py P_CARRIER
P_ABSENT = 0.8       # test/e2e/score_genotypes.py V2_P_ABSENT
LEGACY_CARRIER = ('heterozygous', 'homozygous', 'insertion')
LEGACY_WT = ('wild-type',)
LOCUS_RE = re.compile(r'^(?P<ctg>[^:\s]+):(?P<l>(?:oneside_|polyA_)?\d+)-(?P<r>(?:oneside_|polyA_)?\d+)$')


# ----------------------------------------------------------------------------- small helpers
def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def read_table(path, sep='\t'):
    """-> (header, rows as dicts). '#' comment lines before the header are skipped."""
    with opener(path) as fh:
        header = None
        rows = []
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if header is None:
                if line.startswith('#'):
                    continue
                header = line.split(sep)
                continue
            f = line.split(sep)
            if len(f) < len(header):
                f += [''] * (len(header) - len(f))
            rows.append(dict(zip(header, f)))
    return header or [], rows


def colony_of(path):
    b = os.path.basename(path)
    for ext in ('.txt.gz', '.tsv.gz', '.txt', '.tsv'):
        if b.endswith(ext):
            return b[:-len(ext)]
    return b


def fnum(x, default=float('nan')):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def fint(x, default=0):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return default


def q(vals, p):
    v = sorted(x for x in vals if x is not None and not (isinstance(x, float) and math.isnan(x)))
    if not v:
        return float('nan')
    k = (len(v) - 1) * p
    lo, hi = int(math.floor(k)), int(math.ceil(k))
    return v[lo] if lo == hi else v[lo] + (v[hi] - v[lo]) * (k - lo)


def fmt(x, nd=1):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return 'NA'
    if isinstance(x, int):
        return str(x)
    return f'{x:.{nd}f}'


def md_table(header, rows):
    out = ['| ' + ' | '.join(str(h) for h in header) + ' |', '|' + '---|' * len(header)]
    for r in rows:
        out.append('| ' + ' | '.join(str(c) for c in r) + ' |')
    return '\n'.join(out)


def write_tsv(path, header, rows):
    op = gzip.open(path, 'wt') if path.endswith('.gz') else open(path, 'w')
    with op as fh:
        fh.write('\t'.join(header) + '\n')
        for r in rows:
            fh.write('\t'.join('' if c is None else str(c) for c in r) + '\n')


def tree_tips(path):
    txt = open(path).read()
    txt = re.sub(r'\s', '', txt)
    tips = re.findall(r'[(,]([^(),:;]+)', txt)
    return sorted(set(t.strip("'\"") for t in tips))


class Locus:
    def __init__(self, name):
        self.name = name
        m = LOCUS_RE.match(name)
        self.ctg = self.L = self.R = None
        if not m:
            return
        self.ctg = m.group('ctg')
        lt, rt = m.group('l'), m.group('r')
        if lt.isdigit():
            self.L = int(lt)
        if rt.isdigit():
            self.R = int(rt)

    def near(self, other, tol):
        """'full' when every real side of self lies within tol of the same side of other,
        'partial' when one does, else None."""
        if self.ctg is None or other.ctg is None or self.ctg != other.ctg:
            return None
        hits, sides = 0, 0
        for a, b in ((self.L, other.L), (self.R, other.R)):
            if a is None:
                continue
            sides += 1
            if b is not None and abs(a - b) <= tol:
                hits += 1
        if sides == 0 or hits == 0:
            return None
        return 'full' if hits == sides else 'partial'


# ----------------------------------------------------------------------------- readers
def read_legacy_dir(d):
    """{colony: {locus: row}} from the legacy per-colony files
    (insertion genotype score_genotype score_alternative coverage n_alt n_ref n_art)."""
    out = {}
    for f in sorted(glob.glob(os.path.join(d, '*.txt.gz'))):
        header, rows = read_table(f)
        if 'genotype' not in header:
            print(f'WARNING: {f}: no `genotype` column — not a legacy genotype file, skipped', file=sys.stderr)
            continue
        key = 'insertion' if 'insertion' in header else header[0]
        out[colony_of(f)] = {r[key]: r for r in rows if r.get(key)}
    return out


def read_v2_dir(d):
    """{colony: {locus: row}} from genotype2 per-colony files (header types::OUTPUT_HEADER)."""
    out = {}
    for f in sorted(glob.glob(os.path.join(d, '*.txt.gz'))):
        header, rows = read_table(f)
        if 'p_absent' not in header or 'status' not in header:
            print(f'WARNING: {f}: no numeric genotype2 header — skipped', file=sys.stderr)
            continue
        out[colony_of(f)] = {r['locus']: r for r in rows if r.get('locus')}
    return out


def read_joint(path):
    header, rows = read_table(path)
    cols = [c[2:] for c in header if c.startswith('p_')]
    return {r['locus']: r for r in rows}, cols, header


def read_legacy_calls(path):
    """<P>.genotypes.csv.gz (';'-separated) -> ({locus: [call per colony]}, colonies)."""
    with gzip.open(path, 'rt') as fh:
        header = fh.readline().rstrip('\n').split(';')
        colonies = header[1:]
        calls = {}
        for line in fh:
            f = line.rstrip('\n').split(';')
            if f and f[0]:
                calls[f[0]] = f[1:]
    return calls, colonies


def read_known(path):
    header, rows = read_table(path)
    for r in rows:
        r['carriers'] = [c for c in r.get('carriers', '').split(',') if c]
    return rows


def lsf_reports(logdir, pattern='gt.*.log'):
    """Run time / Max Memory of the LAST LSF report in each -o log (LSF appends on resubmission)."""
    res = []
    for f in sorted(glob.glob(os.path.join(logdir, pattern))):
        try:
            txt = open(f, errors='replace').read()
        except OSError:
            continue
        blocks = [b for b in re.split(r'(?=^Sender: LSF System)', txt, flags=re.M) if 'Resource usage summary' in b]
        if not blocks:
            continue
        b = blocks[-1]
        run = re.search(r'Run time\s*:\s*([\d.]+)\s*sec', b)
        cpu = re.search(r'CPU time\s*:\s*([\d.]+)\s*sec', b)
        mem = re.search(r'Max Memory\s*:\s*([\d.]+)\s*(MB|GB|KB)?', b)
        m = None
        if mem:
            m = float(mem.group(1)) * {'GB': 1024.0, 'KB': 1 / 1024.0}.get(mem.group(2) or 'MB', 1.0)
        res.append({'file': f, 'run_s': fnum(run.group(1)) if run else None,
                    'cpu_s': fnum(cpu.group(1)) if cpu else None, 'mem_mb': m,
                    'ok': 'Successfully completed' in b})
    return res


def v2_own_timing(logdir):
    """genotype2's own summary line: `genotyped N loci in Xs (Y ms/locus, ...` from gt.*.err."""
    out = []
    for f in sorted(glob.glob(os.path.join(logdir, 'gt.*.err'))):
        try:
            txt = open(f, errors='replace').read()
        except OSError:
            continue
        m = re.findall(r'genotyped (\d+) loci in ([\d.]+)s \(([\d.]+) ms/locus', txt)
        if m:
            n, s, ms = m[-1]
            out.append({'file': f, 'n_loci': int(n), 'wall_s': float(s), 'ms_per_locus': float(ms)})
    return out


# ----------------------------------------------------------------------------- bucketing
def legacy_bucket(call):
    if call in LEGACY_CARRIER:
        return 'carrier'
    if call in LEGACY_WT:
        return 'wild-type'
    if call == 'wild-type?':
        return 'wild-type?'
    if call == 'artefact':
        return 'artefact'
    if call in ('', 'NA', 'None'):
        return 'no-call'
    return call


def v2_bucket(row):
    st = row.get('status', '')
    if st != 'ok':
        return st or 'no-status'
    pp = fnum(row.get('p_het')) + fnum(row.get('p_hom'))
    pa = fnum(row.get('p_absent'))
    if pp >= P_PRESENT:
        return 'present'
    if pa >= P_ABSENT:
        return 'absent'
    return 'ambiguous'


V2_ORDER = ('present', 'ambiguous', 'absent', 'no_reads', 'high_coverage', 'error')
LEG_ORDER = ('carrier', 'wild-type', 'wild-type?', 'artefact', 'high-coverage', 'no-call')


def joint_class(jr):
    """ROOT / clade (>= 2 carriers below one branch) / private (1 carrier) / INDEP / NOISE."""
    best = jr.get('best', '')
    if best in ('ROOT', 'INDEP', 'NOISE', ''):
        return best or 'NA'
    n = fint(jr.get('n_carriers'))
    return 'clade' if n >= 2 else 'private'


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--patient', required=True)
    ap.add_argument('--legacy-dir', required=True)
    ap.add_argument('--v2-dir', required=True)
    ap.add_argument('--joint', required=True)
    ap.add_argument('--joint-uniform')
    ap.add_argument('--tree', required=True)
    ap.add_argument('--legacy-calls')
    ap.add_argument('--legacy-fit')
    ap.add_argument('--known')
    ap.add_argument('--known-tol', type=int, default=20)
    ap.add_argument('--v2-logs')
    ap.add_argument('--legacy-logs')
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    P = a.patient
    R = []     # report lines

    legacy = read_legacy_dir(a.legacy_dir)
    v2 = read_v2_dir(a.v2_dir)
    joint, joint_cols, joint_header = read_joint(a.joint)
    joint_u = read_joint(a.joint_uniform)[0] if a.joint_uniform and os.path.exists(a.joint_uniform) else None
    tips = tree_tips(a.tree)

    both = sorted(set(legacy) & set(v2))
    only_leg = sorted(set(legacy) - set(v2))
    only_v2 = sorted(set(v2) - set(legacy))
    loci_leg = set().union(*(set(d) for d in legacy.values())) if legacy else set()
    loci_v2 = set().union(*(set(d) for d in v2.values())) if v2 else set()
    loci = sorted(loci_leg & loci_v2)
    if not both or not loci:
        sys.exit(f'nothing to compare: {len(both)} shared colonies, {len(loci)} shared loci '
                 f'(legacy {len(legacy)} files / {len(loci_leg)} loci; v2 {len(v2)} files / {len(loci_v2)} loci)')

    R.append(f'# {P}: peartree-genotype2 vs legacy genotyper\n')
    R.append('Same contract, same colonies, same BAMs. Legacy = `peartree-genotype` string calls; v2 = '
             '`peartree-genotype2` numeric rows bucketed as present (p_het + p_hom >= '
             f'{P_PRESENT}), absent (p_absent >= {P_ABSENT}), else ambiguous; status != ok is its own bucket.\n')
    R.append('## Inputs\n')
    R.append(md_table(['', 'value'], [
        ['legacy per-colony files', f'{len(legacy)} in `{a.legacy_dir}`'],
        ['v2 per-colony files', f'{len(v2)} in `{a.v2_dir}`'],
        ['colonies compared', f'{len(both)}' + (f' (legacy only: {", ".join(only_leg)})' if only_leg else '')
         + (f' (v2 only: {", ".join(only_v2)})' if only_v2 else '')],
        ['tree tips', f'{len(tips)}; without a v2 file: {len(set(tips) - set(v2))}; v2 files not on the tree: {len(set(v2) - set(tips))}'],
        ['loci compared', f'{len(loci)} (legacy {len(loci_leg)}, v2 {len(loci_v2)})'],
        ['joint step', f'`{a.joint}` ({len(joint)} loci, {len(joint_cols)} colonies)'],
    ]))
    R.append('')

    # ------------------------------------------------------------------ timing / resources
    R.append('## Timing and memory per colony\n')
    trow = []
    if a.v2_logs:
        own = v2_own_timing(a.v2_logs)
        if own:
            trow.append(['v2 (own summary line)', len(own), fmt(q([o["wall_s"] for o in own], .5), 0),
                         fmt(q([o["wall_s"] for o in own], .95), 0), fmt(q([o["ms_per_locus"] for o in own], .5), 2) + ' ms/locus', 'NA'])
        rep = lsf_reports(a.v2_logs)
        if rep:
            trow.append(['v2 (LSF report)', len(rep), fmt(q([r["run_s"] for r in rep], .5), 0),
                         fmt(q([r["run_s"] for r in rep], .95), 0), fmt(q([r["cpu_s"] for r in rep], .5), 0) + ' s CPU',
                         fmt(q([r["mem_mb"] for r in rep], .95), 0) + ' MB'])
    if a.legacy_logs:
        rep = lsf_reports(a.legacy_logs)
        if rep:
            trow.append(['legacy (LSF report)', len(rep), fmt(q([r["run_s"] for r in rep], .5), 0),
                         fmt(q([r["run_s"] for r in rep], .95), 0), fmt(q([r["cpu_s"] for r in rep], .5), 0) + ' s CPU',
                         fmt(q([r["mem_mb"] for r in rep], .95), 0) + ' MB'])
    if trow:
        R.append(md_table(['', 'tasks', 'median wall s', 'p95 wall s', 'per locus / CPU', 'p95 max memory'], trow))
    else:
        R.append('_no LSF logs given (--v2-logs / --legacy-logs)_')
    R.append('')

    # ------------------------------------------------------------------ per-cell cross-tab
    cells = []
    xt = Counter()
    v2_status = Counter()
    leg_calls = Counter()
    depth_diff = []
    alt_leg_only = alt_v2_only = alt_both = 0
    disc_rows = []
    per_locus = {}
    for loc in loci:
        pl = per_locus.setdefault(loc, {'leg_carrier': 0, 'leg_wt': 0, 'leg_art': 0, 'leg_other': 0,
                                        'v2_present': 0, 'v2_absent': 0, 'v2_ambig': 0, 'v2_nodata': 0,
                                        'kind': ''})
        for c in both:
            lr = legacy[c].get(loc)
            vr = v2[c].get(loc)
            if lr is None or vr is None:
                continue
            lb = legacy_bucket(lr.get('genotype', ''))
            vb = v2_bucket(vr)
            xt[(lb, vb)] += 1
            v2_status[vr.get('status', '')] += 1
            leg_calls[lr.get('genotype', '')] += 1
            if not pl['kind']:
                pl['kind'] = vr.get('kind', '')
            pl['leg_carrier' if lb == 'carrier' else 'leg_wt' if lb in ('wild-type', 'wild-type?') else 'leg_art' if lb == 'artefact' else 'leg_other'] += 1
            pl['v2_present' if vb == 'present' else 'v2_absent' if vb == 'absent' else 'v2_ambig' if vb == 'ambiguous' else 'v2_nodata'] += 1
            la, va = fint(lr.get('n_alt')), fint(vr.get('n_alt'))
            if vr.get('status') == 'ok':
                lc, vd = fint(lr.get('coverage')), fint(vr.get('depth'))
                if lc > 0 and vd > 0:
                    depth_diff.append(vd - lc)
                if la > 0 and va > 0:
                    alt_both += 1
                elif la > 0:
                    alt_leg_only += 1
                elif va > 0:
                    alt_v2_only += 1
            pp = fnum(vr.get('p_het')) + fnum(vr.get('p_hom'))
            cells.append((loc, c, lr.get('genotype', ''), lr.get('coverage', ''), lr.get('n_alt', ''), lr.get('n_ref', ''),
                          lr.get('score_genotype', ''), vr.get('status', ''), vr.get('depth', ''), vr.get('n_alt', ''),
                          vr.get('n_ref', ''), vr.get('n_uninf', ''), vr.get('n_art', ''), vr.get('vaf', ''),
                          f'{pp:.4f}' if not math.isnan(pp) else '', vr.get('p_absent', ''), vr.get('gq', ''), lb, vb))
            if (lb == 'carrier' and vb == 'absent') or (lb in ('wild-type',) and vb == 'present'):
                disc_rows.append((loc, c, lr.get('genotype', ''), lr.get('coverage', ''), lr.get('n_alt', ''), lr.get('n_ref', ''),
                                  vr.get('status', ''), vr.get('depth', ''), vr.get('n_alt', ''), vr.get('n_ref', ''),
                                  vr.get('n_uninf', ''), vr.get('vaf', ''), f'{pp:.4f}' if not math.isnan(pp) else '',
                                  'legacy carrier / v2 absent' if lb == 'carrier' else 'legacy wild-type / v2 present'))

    write_tsv(os.path.join(a.out_dir, 'cells.tsv.gz'),
              ['locus', 'colony', 'legacy_genotype', 'legacy_coverage', 'legacy_n_alt', 'legacy_n_ref', 'legacy_score',
               'v2_status', 'v2_depth', 'v2_n_alt', 'v2_n_ref', 'v2_n_uninf', 'v2_n_art', 'v2_vaf', 'v2_p_present',
               'v2_p_absent', 'v2_gq', 'legacy_bucket', 'v2_bucket'], cells)
    n_cells = sum(xt.values())

    R.append(f'## Per (locus, colony): legacy call vs v2 bucket ({n_cells} cells)\n')
    leg_rows = [b for b in LEG_ORDER if any(k[0] == b for k in xt)] + sorted(set(k[0] for k in xt) - set(LEG_ORDER))
    v2_cols = [b for b in V2_ORDER if any(k[1] == b for k in xt)] + sorted(set(k[1] for k in xt) - set(V2_ORDER))
    rows = []
    for lb in leg_rows:
        tot = sum(v for k, v in xt.items() if k[0] == lb)
        rows.append([f'legacy **{lb}**'] + [xt[(lb, vb)] for vb in v2_cols] + [tot])
    rows.append(['total'] + [sum(v for k, v in xt.items() if k[1] == vb) for vb in v2_cols] + [n_cells])
    R.append(md_table(['', *[f'v2 {c}' for c in v2_cols], 'total'], rows))
    R.append('')
    R.append('Legacy call strings: ' + ', '.join(f'{k or "(empty)"} {v}' for k, v in leg_calls.most_common()))
    R.append('')
    R.append('v2 status: ' + ', '.join(f'{k} {v}' for k, v in v2_status.most_common()))
    R.append('')
    R.append('Read counting (cells with v2 status ok):')
    R.append(md_table(['', 'value'], [
        ['median v2 depth − legacy coverage', fmt(q(depth_diff, .5), 0) + f' (p5 {fmt(q(depth_diff, .05), 0)}, p95 {fmt(q(depth_diff, .95), 0)})'],
        ['cells with alt reads in both', alt_both],
        ['alt reads in legacy only', alt_leg_only],
        ['alt reads in v2 only', alt_v2_only],
    ]))
    R.append('')
    R.append('The two denominators differ: v2 `depth` = primary, mapped, non-duplicate reads overlapping the '
             'breakpoint window(s) (PD37590: ~13 reads fewer per cell than legacy `coverage`, which also counted '
             'reads near the locus that never reach the junction register and are neither alt nor ref).')
    R.append('')

    # ------------------------------------------------------------------ discordant cells
    disc_rows.sort(key=lambda r: (r[13], -fint(r[4]) - fint(r[8])))
    write_tsv(os.path.join(a.out_dir, 'discordant.tsv'),
              ['locus', 'colony', 'legacy_genotype', 'legacy_coverage', 'legacy_n_alt', 'legacy_n_ref', 'v2_status',
               'v2_depth', 'v2_n_alt', 'v2_n_ref', 'v2_n_uninf', 'v2_vaf', 'v2_p_present', 'kind'], disc_rows)
    n_ca = sum(1 for r in disc_rows if r[13].startswith('legacy carrier'))
    n_wp = len(disc_rows) - n_ca
    R.append(f'## Hard discordances: legacy carrier → v2 absent ({n_ca}), legacy wild-type → v2 present ({n_wp})\n')
    R.append(f'Full list: `discordant.tsv`. First rows of each kind (sorted by alt reads):\n')
    hdr = ['locus', 'colony', 'legacy', 'leg cov', 'leg alt', 'leg ref', 'v2 depth', 'v2 alt', 'v2 ref', 'v2 uninf', 'v2 vaf', 'v2 P(present)']
    for kind in ('legacy carrier / v2 absent', 'legacy wild-type / v2 present'):
        sub = [r for r in disc_rows if r[13] == kind][:15]
        R.append(f'**{kind}**\n')
        R.append(md_table(hdr, [[r[0], r[1], r[2], r[3], r[4], r[5], r[7], r[8], r[9], r[10], r[11], r[12]] for r in sub]) if sub else '_none_')
        R.append('')

    # ------------------------------------------------------------------ per locus
    jc_all = Counter(joint_class(joint[l]) for l in loci if l in joint)
    pl_rows = []
    cross = Counter()     # (legacy carrier-count band, joint class)
    for loc in loci:
        pl = per_locus[loc]
        jr = joint.get(loc, {})
        jcl = joint_class(jr) if jr else 'NA'
        band = ('0' if pl['leg_carrier'] == 0 else '1' if pl['leg_carrier'] == 1 else '2-5' if pl['leg_carrier'] <= 5
                else f'6-{len(both) - 1}' if pl['leg_carrier'] < len(both) else 'all')
        cross[(band, jcl)] += 1
        pl_rows.append((loc, pl['kind'], pl['leg_carrier'], pl['leg_wt'], pl['leg_art'], pl['leg_other'],
                        pl['v2_present'], pl['v2_ambig'], pl['v2_absent'], pl['v2_nodata'],
                        jr.get('best', ''), jcl, jr.get('n_carriers', ''), jr.get('post_best', ''), jr.get('log10_bf_tree', ''),
                        jr.get('carriers', ''), (joint_u.get(loc, {}).get('best', '') if joint_u else '')))
    write_tsv(os.path.join(a.out_dir, 'per_locus.tsv'),
              ['locus', 'kind', 'legacy_carriers', 'legacy_wt', 'legacy_artefact', 'legacy_other', 'v2_present', 'v2_ambiguous',
               'v2_absent', 'v2_nodata', 'joint_best', 'joint_class', 'joint_n_carriers', 'joint_post_best', 'joint_log10_bf_tree',
               'joint_carriers', 'joint_best_uniform_prior'], pl_rows)

    R.append('## Per locus: legacy carrier count vs the joint step\n')
    R.append('Joint classes: ROOT = every colony (germline or pre-MRCA); clade = one branch with >= 2 tips below it; '
             'private = one tip; INDEP = carriers scatter over the tree; NOISE = absent everywhere.\n')
    jcls = [c for c in ('ROOT', 'clade', 'private', 'INDEP', 'NOISE', 'NA') if any(k[1] == c for k in cross)]
    bands = [b for b in ('0', '1', '2-5', f'6-{len(both) - 1}', 'all') if any(k[0] == b for k in cross)]
    rows = [[f'legacy carriers {b}'] + [cross[(b, c)] for c in jcls] + [sum(v for k, v in cross.items() if k[0] == b)] for b in bands]
    rows.append(['total'] + [jc_all[c] for c in jcls] + [sum(jc_all.values())])
    R.append(md_table(['', *[f'joint {c}' for c in jcls], 'total'], rows))
    R.append('')
    # loci the legacy run never called in any colony but v2 places on a branch with >= 2 carriers
    new_clade = [r for r in pl_rows if r[2] == 0 and r[11] == 'clade']
    new_clade.sort(key=lambda r: -fnum(r[14], -1e9))
    R.append(f'Loci with **no legacy carrier** that the joint step places as a clade event: {len(new_clade)}'
             + (' (top by log10 BF):\n' if new_clade else ''))
    if new_clade:
        R.append(md_table(['locus', 'kind', 'v2 present', 'v2 ambiguous', 'joint carriers', 'post', 'log10 BF'],
                          [[r[0], r[1], r[6], r[7], r[12], r[13], r[14]] for r in new_clade[:20]]))
    R.append('')
    lost = [r for r in pl_rows if r[2] >= 2 and r[11] in ('INDEP', 'NOISE')]
    lost.sort(key=lambda r: -r[2])
    R.append(f'Loci with **>= 2 legacy carriers** that the joint step rejects (INDEP / NOISE): {len(lost)}'
             + (' (most legacy carriers first):\n' if lost else ''))
    if lost:
        R.append(md_table(['locus', 'kind', 'legacy carriers', 'legacy wt', 'v2 present', 'v2 ambiguous', 'joint best', 'post', 'log10 BF'],
                          [[r[0], r[1], r[2], r[3], r[6], r[7], r[10], r[13], r[14]] for r in lost[:20]]))
    R.append('')
    if joint_u:
        chg = Counter((joint_class(joint[l]), joint_class(joint_u[l])) for l in loci if l in joint and l in joint_u)
        moved = [(k, v) for k, v in chg.items() if k[0] != k[1]]
        R.append('Uniform branch prior (`--branch-prior uniform`) vs the length prior: '
                 + (', '.join(f'{k[0]}→{k[1]} {v}' for k, v in sorted(moved, key=lambda kv: -kv[1])) if moved else 'no class changes'))
        R.append('')

    # ------------------------------------------------------------------ legacy combine_genotypes + tree_fit
    if a.legacy_calls and os.path.exists(a.legacy_calls):
        calls, ccol = read_legacy_calls(a.legacy_calls)
        kept = [l for l in calls if l in per_locus]
        jk = Counter(joint_class(joint[l]) if l in joint else 'NA' for l in kept)
        R.append(f'## Legacy combine_genotypes output (`{os.path.basename(a.legacy_calls)}`)\n')
        R.append(f'{len(calls)} loci survived the legacy gates ({len(kept)} in the compared set). The joint step puts them at: '
                 + ', '.join(f'{k} {v}' for k, v in jk.most_common()) + '.')
        R.append('')
    if a.legacy_fit and os.path.exists(a.legacy_fit):
        fh, frows = read_table(a.legacy_fit)
        lcol = 'locus' if 'locus' in fh else fh[0]
        fit = {r[lcol]: r for r in frows}
        ct = Counter()
        for l in loci:
            if l in fit and l in joint:
                ct[(f'{fit[l].get("class", "?")}/{fit[l].get("label", "?")}', joint_class(joint[l]))] += 1
        if ct:
            R.append(f'## Legacy tree_fit (`{a.legacy_fit}`) vs the joint step\n')
            lrows = sorted(set(k[0] for k in ct), key=lambda x: -sum(v for k, v in ct.items() if k[0] == x))
            jcls2 = [c for c in ('ROOT', 'clade', 'private', 'INDEP', 'NOISE') if any(k[1] == c for k in ct)]
            R.append(md_table(['legacy class/label', *[f'joint {c}' for c in jcls2], 'total'],
                              [[lr] + [ct[(lr, c)] for c in jcls2] + [sum(v for k, v in ct.items() if k[0] == lr)] for lr in lrows]))
            R.append('')

    # ------------------------------------------------------------------ known insertions
    if a.known and os.path.exists(a.known):
        known = read_known(a.known)
        R.append(f'## Known insertions (`{os.path.basename(a.known)}`)\n')
        krows = []
        for k in known:
            kl = Locus(k['locus'])
            matches = []
            for loc in loci:
                m = kl.near(Locus(loc), a.known_tol)
                if m:
                    matches.append((loc, m))
            carriers = [c for c in k['carriers'] if c in both]
            missing_c = [c for c in k['carriers'] if c not in both]
            R.append(f'### {k["locus"]} — {k.get("class", "")} {k.get("subclass", "")}, known carriers: '
                     f'{", ".join(k["carriers"]) or "none"}' + (f' ({len(missing_c)} without files)' if missing_c else ''))
            R.append('')
            if not matches:
                R.append('_not in the compared contract (not discovered, or filtered before genotyping)_\n')
                krows.append((k['locus'], k.get('class', ''), ','.join(k['carriers']), '', '', '', '', '', '', ''))
                continue
            for loc, m in matches:
                jr = joint.get(loc, {})
                jcar = set(c for c in jr.get('carriers', '').split(',') if c)
                exact = jcar == set(k['carriers'])
                R.append(f'**contract locus {loc}** ({m} match, kind {per_locus[loc]["kind"]}); joint: best **{jr.get("best", "NA")}** '
                         f'({joint_class(jr) if jr else "NA"}), post {jr.get("post_best", "NA")}, log10 BF {jr.get("log10_bf_tree", "NA")}, '
                         f'carriers {jr.get("n_carriers", "NA")} = ' + (', '.join(sorted(jcar)) or 'none')
                         + (' **(exactly the known carriers)**' if exact else
                            f' (known ∩ joint {len(jcar & set(k["carriers"]))}/{len(k["carriers"])}; extra {len(jcar - set(k["carriers"]))})'))
                R.append('')
                rows = []
                for c in carriers:
                    lr, vr = legacy[c].get(loc, {}), v2[c].get(loc, {})
                    pp = fnum(vr.get('p_het')) + fnum(vr.get('p_hom'))
                    rows.append([c, lr.get('genotype', 'NA'), lr.get('coverage', ''), lr.get('n_alt', ''), lr.get('n_ref', ''),
                                 vr.get('status', 'NA'), vr.get('depth', ''), vr.get('n_alt', ''), vr.get('n_ref', ''), vr.get('n_uninf', ''),
                                 vr.get('vaf', ''), fmt(pp, 3), vr.get('gq', ''), jr.get(f'p_{c}', '')])
                R.append('Known carriers:\n')
                R.append(md_table(['colony', 'legacy', 'cov', 'alt', 'ref', 'v2 status', 'depth', 'alt', 'ref', 'uninf', 'vaf', 'P(present)', 'GQ', 'joint P(carrier)'], rows))
                R.append('')
                others = [c for c in both if c not in k['carriers']]
                o_leg = Counter(legacy_bucket(legacy[c].get(loc, {}).get('genotype', '')) for c in others)
                o_v2 = Counter(v2_bucket(v2[c].get(loc, {'status': 'no-row'})) for c in others)
                o_alt = [(c, fint(v2[c].get(loc, {}).get('n_alt'))) for c in others]
                o_alt = [x for x in o_alt if x[1] > 0]
                R.append(f'Other {len(others)} colonies — legacy: ' + ', '.join(f'{k2} {v}' for k2, v in o_leg.most_common())
                         + '; v2: ' + ', '.join(f'{k2} {v}' for k2, v in o_v2.most_common())
                         + ('; v2 alt reads in: ' + ', '.join(f'{c} ({n})' for c, n in sorted(o_alt, key=lambda x: -x[1])) if o_alt else '; no v2 alt read in any other colony') + '.')
                R.append('')
                krows.append((k['locus'], k.get('class', ''), ','.join(k['carriers']), loc, m,
                              sum(1 for c in carriers if legacy_bucket(legacy[c].get(loc, {}).get('genotype', '')) == 'carrier'),
                              sum(1 for c in carriers if v2_bucket(v2[c].get(loc, {'status': 'no-row'})) == 'present'),
                              len(carriers), jr.get('best', ''), ','.join(sorted(jcar))))
        write_tsv(os.path.join(a.out_dir, 'known.tsv'),
                  ['known_locus', 'class', 'known_carriers', 'contract_locus', 'match', 'legacy_carriers_called',
                   'v2_carriers_present', 'n_known_carriers_with_files', 'joint_best', 'joint_carriers'], krows)

    # ------------------------------------------------------------------ v2 internals
    R.append('## v2 distributions (status ok cells)\n')
    dep, nalt_present, gq_present, vaf_present = [], [], [], []
    for c in both:
        for loc in loci:
            vr = v2[c].get(loc)
            if not vr or vr.get('status') != 'ok':
                continue
            dep.append(fint(vr.get('depth')))
            if v2_bucket(vr) == 'present':
                nalt_present.append(fint(vr.get('n_alt')))
                gq_present.append(fint(vr.get('gq')))
                vaf_present.append(fnum(vr.get('vaf')))
    R.append(md_table(['', 'p5', 'median', 'p95'], [
        ['depth', fmt(q(dep, .05), 0), fmt(q(dep, .5), 0), fmt(q(dep, .95), 0)],
        ['n_alt of present cells', fmt(q(nalt_present, .05), 0), fmt(q(nalt_present, .5), 0), fmt(q(nalt_present, .95), 0)],
        ['VAF of present cells', fmt(q(vaf_present, .05), 2), fmt(q(vaf_present, .5), 2), fmt(q(vaf_present, .95), 2)],
        ['GQ of present cells', fmt(q(gq_present, .05), 0), fmt(q(gq_present, .5), 0), fmt(q(gq_present, .95), 0)],
    ]))
    R.append('')
    R.append('Files: `cells.tsv.gz` (every locus × colony), `per_locus.tsv`, `discordant.tsv`'
             + (', `known.tsv`' if a.known else '') + '.')

    with open(os.path.join(a.out_dir, 'report.md'), 'w') as fh:
        fh.write('\n'.join(R) + '\n')
    print('\n'.join(R))


if __name__ == '__main__':
    main()
