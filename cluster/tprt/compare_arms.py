#!/usr/bin/env python3
"""A vs B report for one patient of the TPRT A/B kit (cluster/tprt/evaluate.sh calls this).

Inputs are the two arms' pipeline run dirs (cluster/pipeline.sh layout) and their evaluation dirs
(tools/phylo/tree_fit.py + discrimination.py output). Either arm may be missing (one-arm report).

    compare_arms.py --patient PD37449 \
        --rundir-a /lustre/.../tprt_ab/PD37449/A/PD37449 --rundir-b /lustre/.../tprt_ab/PD37449/B/PD37449 \
        --eval-a /lustre/.../tprt_ab/PD37449/eval/A --eval-b /lustre/.../tprt_ab/PD37449/eval/B \
        --tree patients/colorectum/PD37449/PD37449_snp_tree_with_branch_length.tree \
        --out-dir /lustre/.../tprt_ab/PD37449/eval

Writes into --out-dir:
  ab_report.md        the report
  stage_counts.tsv    loci per stage (discovery -> contract -> calls -> annotated), per arm
  kind_counts.tsv     locus kinds (TSD / TSD deletion / blunt / L1-mediated del/dup / one-sided)
  class_counts.tsv    annotation class / element / tprt_call of the final calls
  matched_contract.tsv, matched_calls.tsv   fuzzy (+-tol bp) one-to-one A<->B locus matches
  a_only.tsv, b_only.tsv                    unmatched calls with phylo label, tprt_score, carriers
  lsf_resources.tsv   per stage: peak memory / run time / CPU / TERM_* from the LSF job reports

Locus names (SPEC "Pairing modes and locus names"): `contig:L-R` (two-sided; R<L for target-site
and L1-mediated deletions) and `contig:L-oneside_L` / `contig:oneside_R-R` (one-sided; only the
named side is real). Two two-sided loci match when BOTH breakpoints are within --tol; a one-sided
locus matches another locus when its real side is within --tol of that side ("partial").

The phylo-fit table format belongs to tools/phylo/tree_fit.py (plans/tprt_hallmarks/PHYLO_EVAL.md):
its real columns are `locus` and `phylo_label` (= phylo_consistent / ambiguous / phylo_violating for
informative shared loci, noise_violating for the constant-allele-fraction class, else the class:
private / germline / uninformative_depth), which the defaults below pick. Generically it assumes a TSV with
a locus column (locus / insertion / insertion_id / id / name, else the first column) and a label
column (--label-col, else the first of phylo_label / label / phylo_class / class_phylo / verdict /
status / category / fit / consistent / phylo_consistent). Labels are reported verbatim and also
bucketed: violating (viol|incons|conflict|homoplas|discord|reject), consistent
(consist|compat|concord|monophyl|clade|accept|^ok|^fit), other.
"""
import argparse
import bisect
import glob
import gzip
import math
import os
import re
import statistics
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

CARRIER = {'heterozygous', 'homozygous', 'insertion'}
WILDTYPE = {'wild-type', 'wild-type?', 'wildtype', 'wt'}
LABEL_COLS = ('phylo_label', 'label', 'phylo_class', 'class_phylo', 'verdict', 'status', 'category',
              'fit', 'consistent', 'phylo_consistent', 'is_consistent')
LOCUS_COLS = ('locus', 'insertion', 'insertion_id', 'id', 'name')
STAGES = ('sd', 'sd_r2', 'sd_retry', 'combine', 'gt', 'gt_r2', 'gt_retry', 'combine_gt', 'annotate')
STAGE_NAME = {'sd': 'stage+discover (tier1)', 'sd_r2': 'stage+discover (tier2 retry)',
              'sd_retry': 'retry controller (disc)', 'combine': 'combine_insertions',
              'gt': 'genotype (tier1)', 'gt_r2': 'genotype (tier2 retry)',
              'gt_retry': 'retry controller (geno)', 'combine_gt': 'combine_genotypes',
              'annotate': 'annotate'}


# ------------------------------------------------------------------ loci
LOCUS_RE = re.compile(r'^(?P<ctg>[^:\s]+):(?P<l>(?:oneside_|polyA_)?\d+)-(?P<r>(?:oneside_|polyA_)?\d+)$')


class Locus:
    __slots__ = ('name', 'ctg', 'L', 'R', 'kind')

    def __init__(self, name):
        self.name = name
        m = LOCUS_RE.match(name)
        if not m:
            self.ctg, self.L, self.R, self.kind = None, None, None, 'unparsed'
            return
        self.ctg = m.group('ctg')
        lt, rt = m.group('l'), m.group('r')
        if rt.startswith('oneside_'):          # contig:L-oneside_L : LEFT real
            self.L, self.R, self.kind = int(lt), None, 'one-sided (LEFT real)'
        elif lt.startswith('oneside_'):        # contig:oneside_R-R : RIGHT real
            self.L, self.R, self.kind = None, int(rt), 'one-sided (RIGHT real)'
        elif rt.startswith('polyA_'):          # legacy breakpoint + poly-A read (contig:L-polyA_P)
            self.L, self.R, self.kind = int(lt), None, 'breakpoint + poly-A read (legacy)'
        elif lt.startswith('polyA_'):
            self.L, self.R, self.kind = None, int(rt), 'breakpoint + poly-A read (legacy)'
        else:
            self.L, self.R = int(lt), int(rt)
            self.kind = kind_of_gap(self.R - self.L)

    def sides(self):
        out = []
        if self.L is not None:
            out.append(('L', self.L))
        if self.R is not None:
            out.append(('R', self.R))
        return out

    @property
    def two_sided(self):
        return self.L is not None and self.R is not None


def kind_of_gap(gap):
    if gap < -30:
        return 'L1-mediated deletion (gap<-30)'
    if gap < 0:
        return 'target-site deletion (-30..-1)'
    if gap < 2:
        return 'blunt (gap 0..1)'
    if gap <= 40:
        return 'TSD (2..40)'
    return 'L1-mediated dup / long TSD (gap>40)'


def match_loci(names_a, names_b, tol):
    """Greedy one-to-one fuzzy matching. Returns (pairs, a_only, b_only); pairs =
    [(a, b, 'full'|'partial', dL, dR)]."""
    la = [Locus(n) for n in names_a]
    lb = [Locus(n) for n in names_b]
    index = defaultdict(list)                     # (ctg, side) -> sorted [(pos, j)]
    for j, b in enumerate(lb):
        for side, pos in b.sides():
            index[(b.ctg, side)].append((pos, j))
    for v in index.values():
        v.sort()
    cands = []
    for i, a in enumerate(la):
        seen = set()
        for side, pos in a.sides():
            lst = index.get((a.ctg, side), [])
            k = bisect.bisect_left(lst, (pos - tol, -1))
            while k < len(lst) and lst[k][0] <= pos + tol:
                j = lst[k][1]
                k += 1
                if j in seen:
                    continue
                seen.add(j)
                b = lb[j]
                dL = abs(a.L - b.L) if a.L is not None and b.L is not None else None
                dR = abs(a.R - b.R) if a.R is not None and b.R is not None else None
                if a.two_sided and b.two_sided:
                    if dL is not None and dR is not None and dL <= tol and dR <= tol:
                        cands.append((0, dL + dR, i, j, 'full', dL, dR))
                else:
                    d = [x for x in (dL, dR) if x is not None]
                    if d and min(d) <= tol:
                        cands.append((1, min(d), i, j, 'partial', dL, dR))
    cands.sort()
    used_a, used_b, pairs = set(), set(), []
    for _, _, i, j, how, dL, dR in cands:
        if i in used_a or j in used_b:
            continue
        used_a.add(i)
        used_b.add(j)
        pairs.append((la[i].name, lb[j].name, how, dL, dR))
    a_only = [la[i].name for i in range(len(la)) if i not in used_a]
    b_only = [lb[j].name for j in range(len(lb)) if j not in used_b]
    return pairs, a_only, b_only


# ------------------------------------------------------------------ readers
def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def read_tsv(path, sep='\t'):
    """-> (header list, list of row dicts). Tolerates a missing file (None, [])."""
    if not path or not os.path.exists(path):
        return None, []
    with opener(path) as fh:
        header = fh.readline().rstrip('\n').split(sep)
        rows = []
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\n').split(sep)
            rows.append(dict(zip(header, f)))
    return header, rows


def fastq_loci(path):
    """Unique `contig:L-R` locus ids of a 4-line FASTQ (discovery .txt.gz / combined.txt.gz)."""
    loci = set()
    try:
        with gzip.open(path, 'rt') as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0 and line.startswith('@'):
                    parts = line[1:].split(':', 2)
                    if len(parts) >= 2:
                        loci.add(parts[0] + ':' + parts[1].split()[0])
    except (OSError, EOFError) as e:
        print(f"WARNING: cannot read {path}: {e}", file=sys.stderr)
    return path, loci


def contract_loci(path):
    if not path or not os.path.exists(path):
        return None
    out = []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                out.append(line[1:].strip().split()[0])
    return out


def read_calls(path):
    """<P>.genotypes.csv.gz (';'-separated: insertion;<colony>...) -> {locus: Counter of calls}."""
    if not path or not os.path.exists(path):
        return None, []
    calls = {}
    with gzip.open(path, 'rt') as fh:
        header = fh.readline().rstrip('\n').split(';')
        colonies = header[1:]
        for line in fh:
            f = line.rstrip('\n').split(';')
            if not f or not f[0]:
                continue
            calls[f[0]] = f[1:]
    return calls, colonies


def carrier_class(gts, germline_frac):
    n_car = sum(1 for g in gts if g in CARRIER)
    n_wt = sum(1 for g in gts if g in WILDTYPE)
    inf = n_car + n_wt
    if n_car == 0:
        return n_car, inf, 'no carrier'
    if n_car == 1:
        return n_car, inf, 'private'
    if inf and n_car >= max(2, math.ceil(germline_frac * inf)):
        return n_car, inf, 'germline-like'
    return n_car, inf, 'shared'


def bucket_label(v):
    s = str(v).strip().lower()
    if s in ('true', '1', 'yes'):
        return 'consistent'
    if s in ('false', '0', 'no'):
        return 'violating'
    if re.search(r'viol|incons|conflict|homoplas|discord|reject', s):
        return 'violating'
    if re.search(r'consist|compat|concord|monophyl|clade|accept|^ok|^fit', s):
        return 'consistent'
    return 'other'


def read_fit(eval_dir, label_col=None):
    """-> ({locus: raw label}, label column, locus column, path) from <eval>/fit/phylo_fit.tsv."""
    if not eval_dir:
        return None
    path = os.path.join(eval_dir, 'fit', 'phylo_fit.tsv')
    header, rows = read_tsv(path)
    if header is None:
        return None
    lcol = next((c for c in LOCUS_COLS if c in header), header[0])
    if label_col:
        if label_col not in header:
            print(f"WARNING: --label-col {label_col} not in {path} ({header})", file=sys.stderr)
            return {'labels': {}, 'label_col': None, 'locus_col': lcol, 'path': path, 'header': header}
        tcol = label_col
    else:
        tcol = next((c for c in LABEL_COLS if c in header), None)
    labels = {r[lcol]: (r.get(tcol, '') if tcol else '') for r in rows}
    return {'labels': labels, 'label_col': tcol, 'locus_col': lcol, 'path': path, 'header': header}


def lsf_reports(logdir):
    """Parse the LSF job report appended to every -o log. Last report of a file wins (LSF -o
    appends, so a resubmitted element carries several)."""
    res = defaultdict(list)
    for f in sorted(glob.glob(os.path.join(logdir, '*.log'))):
        stage = os.path.basename(f).split('.', 1)[0]
        try:
            txt = open(f, errors='replace').read()
        except OSError:
            continue
        blocks = re.split(r'(?=^Sender: LSF System)', txt, flags=re.M)
        blocks = [b for b in blocks if 'Resource usage summary' in b or 'Sender: LSF System' in b]
        if not blocks:
            continue
        b = blocks[-1]
        mem = re.search(r'Max Memory\s*:\s*([\d.]+)\s*(MB|GB|KB)?', b)
        run = re.search(r'Run time\s*:\s*([\d.]+)\s*sec', b)
        cpu = re.search(r'CPU time\s*:\s*([\d.]+)\s*sec', b)
        m = None
        if mem:
            m = float(mem.group(1)) * {'GB': 1024.0, 'KB': 1 / 1024.0}.get(mem.group(2) or 'MB', 1.0)
        res[stage].append({
            'mem_mb': m,
            'run_s': float(run.group(1)) if run else None,
            'cpu_s': float(cpu.group(1)) if cpu else None,
            'memlimit': 'TERM_MEMLIMIT' in b,
            'runlimit': 'TERM_RUNLIMIT' in b,
            'ok': 'Successfully completed' in b,
        })
    return res


def pct(vals, q):
    v = sorted(x for x in vals if x is not None)
    if not v:
        return None
    k = min(len(v) - 1, max(0, int(math.ceil(q * len(v))) - 1))
    return v[k]


# ------------------------------------------------------------------ one arm
class Arm:
    def __init__(self, name, rundir, evaldir, patient, args):
        self.name, self.rundir, self.evaldir, self.P = name, rundir, evaldir, patient
        self.present = bool(rundir) and os.path.isdir(rundir)
        self.stage = {}
        self.disc_union = set()
        self.contract = None
        self.calls, self.colonies = None, []
        self.ann = {}
        self.ann_header = None
        self.fit = None
        self.lsf = {}
        self.carriers = {}
        self.fail_reasons = Counter()
        if not self.present:
            return
        rd, P = rundir, patient
        st = self.stage
        samples = os.path.join(rd, 'samples.tsv')
        st['colonies listed'] = sum(1 for _ in open(samples)) if os.path.exists(samples) else None
        st['colonies without data (missing/)'] = len(os.listdir(os.path.join(rd, 'missing'))) if os.path.isdir(os.path.join(rd, 'missing')) else None
        for ph in ('discover', 'genotype'):
            _, ex = read_tsv(os.path.join(rd, f'{ph}_excluded.tsv'))
            st[f'colonies excluded at {ph}'] = len(ex)
        disc = sorted(glob.glob(os.path.join(rd, 'discovery', '*.txt.gz')))
        st['discovery files'] = len(disc)
        side = [f for f in glob.glob(os.path.join(rd, 'discovery', '*.txt.gz.evidence.tsv.gz'))]
        st['evidence sidecars (MB on disk)'] = round(sum(os.path.getsize(os.path.realpath(f)) for f in side) / 1e6, 1) if side else 0
        if disc and not args.skip_discovery_scan:
            per = []
            with Pool(max(1, args.threads)) as pool:
                for _, loci in pool.imap_unordered(fastq_loci, disc):
                    per.append(len(loci))
                    self.disc_union |= loci
            st['discovery loci / colony (median)'] = int(statistics.median(per)) if per else 0
            st['discovery loci / colony (max)'] = max(per) if per else 0
            st['discovery loci pooled (unique names)'] = len(self.disc_union)
        comb = os.path.join(rd, 'insertions', f'{P}.combined.txt.gz')
        if os.path.exists(comb):
            st['combined insertions (combined.txt.gz)'] = len(fastq_loci(comb)[1])
        self.contract = contract_loci(os.path.join(rd, 'insertions', f'{P}.genotyping.txt.gz'))
        st['genotyping contract loci'] = len(self.contract) if self.contract is not None else None
        ext = contract_loci(os.path.join(rd, 'insertions', f'{P}.genotyping.tprt.txt.gz'))
        if ext is not None:
            st['contract + one-sided loci (genotyping.tprt)'] = len(ext)
            self.contract = ext
        _, ev = read_tsv(os.path.join(rd, 'insertions', f'{P}.insertions.evidence.tsv.gz'))
        if ev:
            ins_sup = defaultdict(set)
            reasons = Counter()
            for r in ev:
                ins_sup[r.get('insertion_id')].add(r.get('supported'))
                fr = r.get('fail_reason', '')
                if fr and fr not in ('.', ''):
                    reasons[fr.split(':', 1)[0]] += 1
            st['evidence: insertions with an unsupported junction'] = sum(1 for v in ins_sup.values() if '0' in v)
            self.fail_reasons = reasons
        else:
            self.fail_reasons = Counter()
        self.calls, self.colonies = read_calls(os.path.join(rd, f'{P}.genotypes.csv.gz'))
        st['final calls (genotypes.csv.gz)'] = len(self.calls) if self.calls is not None else None
        self.ann_header, rows = read_tsv(os.path.join(rd, f'{P}.annotated.csv.gz'))
        self.ann = {r['locus']: r for r in rows if 'locus' in r}
        st['annotated loci'] = len(self.ann) if self.ann_header else None
        self.carriers = {}
        if self.calls:
            for loc, gts in self.calls.items():
                self.carriers[loc] = carrier_class(gts, args.germline_frac)
        self.fit = read_fit(evaldir, args.label_col)
        self.lsf = lsf_reports(os.path.join(rd, 'logs'))

    # convenience
    def label(self, loc):
        if not self.fit:
            return ''
        return self.fit['labels'].get(loc, '')

    def row_for(self, loc):
        a = self.ann.get(loc, {})
        n_car, inf, cc = self.carriers.get(loc, ('', '', ''))
        lab = self.label(loc)
        return {'locus': loc, 'kind': Locus(loc).kind, 'class': a.get('class', ''),
                'element': a.get('element', ''), 'structure': a.get('structure', ''),
                'tags': a.get('tags', ''), 'tprt_score': a.get('tprt_score', ''),
                'tprt_call': a.get('tprt_call', ''), 'carriers': n_car, 'informative': inf,
                'carrier_class': cc, 'phylo_label': lab,
                'phylo_bucket': bucket_label(lab) if lab != '' else ''}


# ------------------------------------------------------------------ report helpers
def md_table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join('---' for _ in header) + '|']
    for r in rows:
        out.append('| ' + ' | '.join('' if v is None else str(v) for v in r) + ' |')
    return '\n'.join(out)


def write_tsv(path, header, rows):
    with open(path, 'w') as fh:
        fh.write('\t'.join(header) + '\n')
        for r in rows:
            fh.write('\t'.join('' if v is None else str(v) for v in r) + '\n')


def score_key(r):
    try:
        return -float(r['tprt_score'])
    except (TypeError, ValueError):
        return float('inf')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--patient', required=True)
    ap.add_argument('--rundir-a')
    ap.add_argument('--rundir-b')
    ap.add_argument('--eval-a')
    ap.add_argument('--eval-b')
    ap.add_argument('--tree', help='newick (only to count tips for the report)')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--tol', type=int, default=10, help='fuzzy locus match tolerance in bp (default 10)')
    ap.add_argument('--label-col', help='phylo_fit.tsv label column (default: auto-detect)')
    ap.add_argument('--germline-frac', type=float, default=0.8,
                    help='carriers >= this fraction of informative colonies = germline-like (default 0.8)')
    ap.add_argument('--threads', type=int, default=int(os.environ.get('LSB_DJOB_NUMPROC', '2')))
    ap.add_argument('--skip-discovery-scan', action='store_true', help='do not parse discovery FASTQs')
    ap.add_argument('--top', type=int, default=30, help='rows of A-only / B-only shown in the report')
    ap.add_argument('--label-b', default='B', help="name of the second arm in the report (e.g. C); report file a<label>_report.md")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    A = Arm('A', args.rundir_a, args.eval_a, args.patient, args)
    B = Arm(args.label_b, args.rundir_b, args.eval_b, args.patient, args)
    LB = args.label_b
    arms = [x for x in (A, B) if x.present]
    if not arms:
        sys.exit('neither run dir exists')
    P = args.patient
    md = [f'# TPRT A/B report — {P}', '',
          'A = current pipeline (config.discovery.grch38 / config.genotype.grch38 / config.py.grch38); '
          'B = TPRT-hallmark pipeline (the `.tprt` configs). Same commit, same binaries, same BAMs; '
          'both arms annotated with the same annotate config (rte_library + hs1 gene model), so class / '
          'tprt_score differences come from the calls, not the annotator. Arm A has no evidence '
          'sidecars, so its tprt_score is computed from the junction strings alone (weaker by design).', '']
    if args.tree and os.path.exists(args.tree):
        tips = set(re.findall(r'[(,]([^(),:;\s]+)', open(args.tree).read()))
        md.append(f'Tree: `{args.tree}` ({len(tips)} tips).')
    for x in (A, B):
        md.append(f'- arm {x.name}: `{x.rundir}`' + ('' if x.present else ' — **absent**'))
    md.append('')

    # ---- stage counts
    keys = []                                   # union, each new key placed after its predecessor
    for x in arms:
        prev = None
        for k in x.stage:
            if k not in keys:
                keys.insert(keys.index(prev) + 1 if prev in keys else len(keys), k)
            prev = k
    rows = [[k, A.stage.get(k, ''), B.stage.get(k, '')] for k in keys]
    write_tsv(os.path.join(args.out_dir, 'stage_counts.tsv'), ['stage', 'A', LB], rows)
    md += ['## Loci per stage', '', md_table(['stage', 'A', LB], rows), '']
    if B.present and B.fail_reasons:
        md += ['B combine gate (`insertions.evidence.tsv.gz` fail_reason, junction rows):', '',
               md_table(['fail_reason', 'rows'], B.fail_reasons.most_common()), '']

    # ---- kinds
    kind_rows = []
    for label, getter in (('discovery (pooled)', lambda x: x.disc_union),
                          ('contract', lambda x: x.contract or []),
                          ('final calls', lambda x: list(x.calls or {}))):
        ca = Counter(Locus(n).kind for n in getter(A)) if A.present else Counter()
        cb = Counter(Locus(n).kind for n in getter(B)) if B.present else Counter()
        for k in sorted(set(ca) | set(cb)):
            kind_rows.append([label, k, ca.get(k, 0), cb.get(k, 0)])
    write_tsv(os.path.join(args.out_dir, 'kind_counts.tsv'), ['level', 'kind', 'A', LB], kind_rows)
    md += ['## Locus kinds (from the locus name geometry)', '', md_table(['level', 'kind', 'A', LB], kind_rows), '']

    # ---- classes
    class_rows = []
    for col in ('class', 'element', 'tprt_call'):
        ca = Counter(r.get(col, '') for r in A.ann.values() if col in r)
        cb = Counter(r.get(col, '') for r in B.ann.values() if col in r)
        for k in sorted(set(ca) | set(cb), key=lambda k: -(ca.get(k, 0) + cb.get(k, 0))):
            class_rows.append([col, k, ca.get(k, 0), cb.get(k, 0)])
    write_tsv(os.path.join(args.out_dir, 'class_counts.tsv'), ['column', 'value', 'A', LB], class_rows)
    md += ['## Annotation of the final calls', '', md_table(['column', 'value', 'A', LB], class_rows), '']

    # ---- carriers x phylo
    md += [f'## Carriers and phylogeny', '',
           f'Carrier class from `<P>.genotypes.csv.gz` (carrier = heterozygous/homozygous/insertion; '
           f'informative = carrier + wild-type; germline-like = carriers >= {args.germline_frac:g} x informative).', '']
    for x in arms:
        cc = Counter(v[2] for v in x.carriers.values())
        md.append(f'- arm {x.name}: ' + ', '.join(f'{k} {v}' for k, v in sorted(cc.items())))
    md.append('')
    for x in arms:
        if not x.fit:
            md += [f'Arm {x.name}: no phylo fit (`{os.path.join(x.evaldir or "", "fit", "phylo_fit.tsv")}` absent — '
                   'tree_fit not run / failed, see the eval logs).', '']
            continue
        f = x.fit
        md += [f'Arm {x.name} phylo fit: `{f["path"]}` (locus column `{f["locus_col"]}`, label column '
               f'`{f["label_col"]}`; {len(f["labels"])} loci).', '']
        raw = Counter(f['labels'].values())
        md += [md_table(['phylo label (raw)', 'loci', 'bucket'], [[k, v, bucket_label(k)] for k, v in raw.most_common()]), '']
        xt = defaultdict(Counter)
        for loc in (x.calls or {}):
            lab = f['labels'].get(loc)
            xt[x.carriers[loc][2]][bucket_label(lab) if lab is not None else 'not in fit'] += 1
        cols = ['consistent', 'violating', 'other', 'not in fit']
        md += [md_table([f'arm {x.name} carrier class'] + cols,
                        [[k] + [xt[k].get(c, 0) for c in cols] for k in sorted(xt)]), '']
        if x.ann:
            xc = defaultdict(Counter)
            for loc, r in x.ann.items():
                lab = f['labels'].get(loc)
                xc[r.get('tprt_call') or r.get('class', '')][bucket_label(lab) if lab is not None else 'not in fit'] += 1
            md += [md_table([f'arm {x.name} tprt_call (or class)'] + cols,
                            [[k] + [xc[k].get(c, 0) for c in cols] for k in sorted(xc)]), '']
    for x in arms:
        if not x.evaldir:
            continue
        ddir = os.path.join(x.evaldir, 'discrimination')
        files = sorted(glob.glob(os.path.join(ddir, '*')))
        if files:
            md += [f'Arm {x.name} discrimination output (`{ddir}`): ' + ', '.join(os.path.basename(p) for p in files), '']
            for p in files:
                if p.endswith('.md') and os.path.getsize(p) < 20000:
                    md += [f'<details><summary>{os.path.basename(p)}</summary>', '', open(p).read(), '', '</details>', '']
                elif p.endswith('.tsv') and os.path.getsize(p) < 20000:
                    h, rr = read_tsv(p)
                    if h and len(rr) <= 40:
                        md += [f'`{os.path.basename(p)}`', '', md_table(h, [[r.get(c, '') for c in h] for r in rr]), '']

    # ---- overlap
    if A.present and B.present:
        md += [f'## Overlap A vs {LB} (fuzzy, +-{args.tol} bp, one-to-one)', '']
        ov = []
        for level, la, lb, fn in (('contract', A.contract or [], B.contract or [], 'matched_contract.tsv'),
                                  ('final calls', list(A.calls or {}), list(B.calls or {}), 'matched_calls.tsv')):
            pairs, ao, bo = match_loci(la, lb, args.tol)
            full = sum(1 for p in pairs if p[2] == 'full')
            exact = sum(1 for p in pairs if p[0] == p[1])
            ov.append([level, len(la), len(lb), len(pairs), exact, full, len(pairs) - full, len(ao), len(bo)])
            rows = []
            for a, b, how, dL, dR in pairs:
                ra, rb = A.row_for(a), B.row_for(b)
                rows.append([a, b, how, dL, dR, ra['class'], rb['class'], ra['carriers'], rb['carriers'],
                             ra['phylo_label'], rb['phylo_label'], ra['tprt_score'], rb['tprt_score']])
            write_tsv(os.path.join(args.out_dir, fn),
                      ['locus_a', 'locus_b', 'match', 'dL', 'dR', 'class_a', 'class_b', 'carriers_a', 'carriers_b',
                       'phylo_a', 'phylo_b', 'tprt_score_a', 'tprt_score_b'], rows)
            if level == 'final calls':
                a_only, b_only, call_pairs = ao, bo, rows
        md += [md_table(['level', 'A', LB, 'matched', 'identical name', 'full', 'partial (one-sided)', 'A-only', f'{LB}-only'], ov), '']
        cols = ['locus', 'kind', 'class', 'element', 'structure', 'tags', 'tprt_score', 'tprt_call',
                'carriers', 'informative', 'carrier_class', 'phylo_label', 'phylo_bucket']
        for nm, arm, lst in ((f'{LB}-only', B, b_only), ('A-only', A, a_only)):
            rows = sorted((arm.row_for(l) for l in lst), key=score_key)
            write_tsv(os.path.join(args.out_dir, f'{nm[0].lower()}_only.tsv'), cols, [[r[c] for c in cols] for r in rows])
            pb = Counter(r['phylo_bucket'] or 'no label' for r in rows)
            cc = Counter(r['carrier_class'] for r in rows)
            tc = Counter(r['tprt_call'] or '.' for r in rows)
            md += [f'### {nm} calls ({len(rows)}; `{nm[0].lower()}_only.tsv`)', '',
                   'phylo: ' + ', '.join(f'{k} {v}' for k, v in pb.most_common()) + '  ',
                   'carriers: ' + ', '.join(f'{k} {v}' for k, v in cc.most_common()) + '  ',
                   'tprt_call: ' + ', '.join(f'{k} {v}' for k, v in tc.most_common()), '']
            show = ['locus', 'kind', 'class', 'tprt_score', 'tprt_call', 'carriers', 'carrier_class', 'phylo_label']
            md += [md_table(show, [[r[c] for c in show] for r in rows[:args.top]]), '']
            if len(rows) > args.top:
                md += [f'({len(rows) - args.top} more in `{nm[0].lower()}_only.tsv`)', '']
        disagree = [r for r in call_pairs if r[9] and r[10] and bucket_label(r[9]) != bucket_label(r[10])]
        md += [f'Matched calls whose phylo bucket differs between arms: {len(disagree)} (see `matched_calls.tsv`).', '']

    # ---- LSF resources
    lrows = []
    for x in arms:
        for st in STAGES:
            recs = x.lsf.get(st)
            if not recs:
                continue
            mem = [r['mem_mb'] for r in recs]
            run = [r['run_s'] for r in recs]
            cpu = sum(r['cpu_s'] or 0 for r in recs) / 3600.0
            lrows.append([x.name, STAGE_NAME.get(st, st), len(recs),
                          None if pct(mem, .5) is None else round(pct(mem, .5)),
                          None if pct(mem, .95) is None else round(pct(mem, .95)),
                          None if pct(mem, 1) is None else round(pct(mem, 1)),
                          None if pct(run, .5) is None else round(pct(run, .5) / 60, 1),
                          None if pct(run, 1) is None else round(pct(run, 1) / 60, 1),
                          round(cpu, 2), sum(r['memlimit'] for r in recs), sum(r['runlimit'] for r in recs),
                          sum(1 for r in recs if not r['ok'])])
    hdr = ['arm', 'stage', 'jobs', 'mem median MB', 'mem p95 MB', 'mem max MB', 'run median min',
           'run max min', 'CPU h', 'TERM_MEMLIMIT', 'TERM_RUNLIMIT', 'not ok']
    write_tsv(os.path.join(args.out_dir, 'lsf_resources.tsv'), hdr, lrows)
    md += ['## Runtime / memory per stage (LSF job reports in `<rundir>/logs/*.log`)', '',
           'stage+discover includes the stageBam.pl copy when the task staged its BAM (arm A stages, '
           'arm B reuses). A stage absent here has no finished LSF report yet.', '',
           md_table(hdr, lrows) if lrows else '(no LSF job reports found)', '']

    out = os.path.join(args.out_dir, f'a{LB.lower()}_report.md')
    with open(out, 'w') as fh:
        fh.write('\n'.join(md) + '\n')
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
