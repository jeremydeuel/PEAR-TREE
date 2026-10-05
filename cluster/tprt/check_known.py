#!/usr/bin/env python3
"""Trace a patient's known insertions through every stage of one arm (cluster/tprt/evaluate.sh calls this
when patients/<organ>/<P>/known_insertions.tsv exists).

    check_known.py --known patients/colorectum/PD37590/known_insertions.tsv \
        --rundir /lustre/.../tprt_ab/PD37590/C/PD37590 --eval /lustre/.../tprt_ab/PD37590/eval/C \
        --out-dir /lustre/.../tprt_ab/PD37590/eval/C/known

Per known locus (a stage locus matches when one of its real breakpoints lies within --tol bp of the
same breakpoint of the known locus; 'full' when both sides of a two-sided pair do):
  discovery    which carrier colonies' discovery files hold the locus (arm C needs >= 2 fragments in
               each of them, so a single-carrier insertion lives or dies in one colony), and how many
               non-carrier colonies also hold it
  combined     <P>.combined.txt.gz;  contract  <P>.genotyping(.tprt).txt.gz
  genotyped    <P>.genotypes.unfiltered.csv.gz (before the combine_genotypes gates) and
               <P>.genotypes.csv.gz (after): carriers called vs the known carriers
  annotated    class / tprt_call of <P>.annotated.csv.gz;  phylo  label of eval/fit/phylo_fit.tsv
The first stage a known locus is missing from is its loss point.

Also lists the final calls that match no known locus ("potentially more"): every non-germline call
with its class and phylo label, phylo-consistent shared and RTE-classed calls first.

Writes known_report.md, known_trace.tsv, new_calls.tsv into --out-dir.
"""
import argparse
import glob
import os
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_arms import (CARRIER, WILDTYPE, Locus, bucket_label, carrier_class, contract_loci,  # noqa: E402
                          fastq_loci, md_table, read_calls, read_fit, read_tsv, write_tsv)

RTE_CLASSES = ('LINE', 'L1', 'SINE', 'Alu', 'SVA', 'ALU', 'RTE_other', 'pseudogene', 'Pseudogene')


def read_known(path):
    rows = []
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if header is None:
                header = f
                continue
            r = dict(zip(header, f + [''] * (len(header) - len(f))))
            r['carriers'] = [c for c in r.get('carriers', '').split(',') if c]
            rows.append(r)
    return rows


def near(known, names, tol):
    """names within tol of the known locus -> [(name, 'full'|'partial', max side distance)]."""
    k = Locus(known)
    hits = []
    for n in names:
        b = Locus(n)
        if b.ctg != k.ctg:
            continue
        d = {}
        for side, pos in b.sides():
            kp = k.L if side == 'L' else k.R
            if kp is not None:
                d[side] = abs(pos - kp)
        close = {s: v for s, v in d.items() if v <= tol}
        if not close:
            continue
        how = 'full' if k.two_sided and b.two_sided and len(close) == 2 else 'partial'
        hits.append((n, how, max(close.values())))
    hits.sort(key=lambda h: (h[1] != 'full', h[2]))
    return hits


def best(known, names, tol):
    h = near(known, names, tol)
    return h[0] if h else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--known', required=True)
    ap.add_argument('--rundir', required=True)
    ap.add_argument('--eval', help='eval/<arm> dir (fit/phylo_fit.tsv)')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--patient', help='default: basename of --rundir')
    ap.add_argument('--label', default='', help='arm label for the report title')
    ap.add_argument('--tol', type=int, default=30, help='breakpoint tolerance in bp (default 30)')
    ap.add_argument('--near', type=int, default=1000, help='search radius for a missed carrier\'s nearest locus')
    ap.add_argument('--germline-frac', type=float, default=0.8)
    ap.add_argument('--threads', type=int, default=int(os.environ.get('LSB_DJOB_NUMPROC', '4')))
    args = ap.parse_args()

    rd = args.rundir.rstrip('/')
    P = args.patient or os.path.basename(rd)
    os.makedirs(args.out_dir, exist_ok=True)
    known = read_known(args.known)
    if not known:
        sys.exit(f'no loci in {args.known}')

    # ---- discovery, per colony
    disc_files = sorted(f for f in glob.glob(os.path.join(rd, 'discovery', '*.txt.gz'))
                        if not f.endswith('.evidence.tsv.gz'))
    per_colony = {}
    with Pool(max(1, args.threads)) as pool:
        for path, loci in pool.imap_unordered(fastq_loci, disc_files):
            per_colony[os.path.basename(path)[:-len('.txt.gz')]] = loci
    contract_path = os.path.join(rd, 'insertions', f'{P}.genotyping.tprt.txt.gz')
    if not os.path.exists(contract_path):
        contract_path = os.path.join(rd, 'insertions', f'{P}.genotyping.txt.gz')
    comb = os.path.join(rd, 'insertions', f'{P}.combined.txt.gz')
    combined = fastq_loci(comb)[1] if os.path.exists(comb) else None
    contract = contract_loci(contract_path)
    unf, unf_cols = read_calls(os.path.join(rd, f'{P}.genotypes.unfiltered.csv.gz'))
    fin, fin_cols = read_calls(os.path.join(rd, f'{P}.genotypes.csv.gz'))
    ann_header, ann_rows = read_tsv(os.path.join(rd, f'{P}.annotated.csv.gz'))
    ann = {r['locus']: r for r in ann_rows if 'locus' in r}
    fit = read_fit(args.eval) if args.eval else None
    labels = fit['labels'] if fit else {}

    def gt_summary(calls, cols, locus, truth):
        if calls is None:
            return 'n/a', ''
        m = best(locus, calls.keys(), args.tol)
        if not m:
            return 'absent', ''
        gts = dict(zip(cols, calls[m[0]]))
        called = {c for c, g in gts.items() if g in CARRIER}
        hit = sorted(called & set(truth))
        miss = sorted(set(truth) - called)
        extra = sorted(called - set(truth))
        miss_s = ','.join(f'{c}={gts.get(c, "not genotyped")}' for c in miss)
        s = f'{m[0]} ({m[1]}): {len(hit)}/{len(truth)} carriers'
        if miss:
            s += f'; missed {miss_s}'
        if extra:
            s += f'; +{len(extra)} extra ({",".join(extra[:6])}{"..." if len(extra) > 6 else ""})'
        return s, m[0]

    trace, md_rows, detail = [], [], []
    matched_final = set()
    for k in known:
        loc, car = k['locus'], k['carriers']
        disc_car = [c for c in car if c in per_colony and best(loc, per_colony[c], args.tol)]
        no_file = [c for c in car if c not in per_colony]
        disc_other = sorted(c for c, loci in per_colony.items() if c not in car and best(loc, loci, args.tol))
        disc_names = sorted({h[0] for c in disc_car for h in near(loc, per_colony[c], args.tol)})
        other_names = sorted({h[0] for c in disc_other for h in near(loc, per_colony[c], args.tol)})
        # carriers whose (finished) discovery file lacks the locus: nearest locus within --near bp
        missed = []
        for c in car:
            if c in per_colony and c not in disc_car:
                h = near(loc, per_colony[c], args.near)
                missed.append(f'{c}: ' + (f'nearest {h[0][0]} ({h[0][2]} bp)' if h else f'nothing within {args.near} bp'))
        cm = best(loc, combined, args.tol) if combined is not None else None
        ct = best(loc, contract, args.tol) if contract is not None else None
        su, _ = gt_summary(unf, unf_cols, loc, car)
        sf, fname = gt_summary(fin, fin_cols, loc, car)
        if fname:
            matched_final.add(fname)
        am = best(loc, ann.keys(), args.tol) if ann_header else None
        a = ann.get(am[0], {}) if am else {}
        acls = (f"{a.get('class', '')} / {a.get('tprt_call', '') or '-'} (score {a.get('tprt_score', '') or '-'})"
                if am else ('absent' if ann_header else 'n/a'))
        lab = labels.get(am[0], '') if am else ''
        stages = [('discovery', bool(disc_car)), ('combined', cm is not None if combined is not None else None),
                  ('contract', ct is not None if contract is not None else None),
                  ('genotyped', not su.startswith('absent') if unf is not None else None),
                  ('final calls', bool(fname) if fin is not None else None),
                  ('annotated', am is not None if ann_header else None)]
        lost = next((s for s, ok in stages if ok is False), None)
        if lost is None:
            reached = [s for s, ok in stages if ok]
            lost = 'found' if len(reached) == len(stages) else f'pending (reached {reached[-1] if reached else "-"})'
        if lost == 'discovery':
            pending = [c for c in car if c not in per_colony]
            if pending and len(pending) == len(car):
                lost = 'pending (no carrier discovery file yet)'
            elif pending:
                lost = f'discovery so far ({len(pending)} carrier file(s) pending)'

        def mark(ok):
            return 'n/a' if ok is None else ('✓' if ok else '✗')
        trace.append({'locus': loc, 'tier': k.get('tier', ''), 'class_known': f"{k.get('class', '')}/{k.get('subclass', '')}",
                      'n_carriers': len(car), 'discovery_carriers': f'{len(disc_car)}/{len(car)}',
                      'discovery_carrier_ids': ','.join(disc_car), 'carriers_no_discovery_file': ','.join(no_file),
                      'carriers_missed': '; '.join(missed),
                      'discovery_noncarriers': len(disc_other), 'discovery_noncarrier_ids': ','.join(disc_other),
                      'discovery_names': ','.join(disc_names[:5]), 'discovery_noncarrier_names': ','.join(other_names[:5]),
                      'combined': cm[0] if cm else '', 'contract': ct[0] if ct else '',
                      'genotyped_unfiltered': su, 'final_calls': sf, 'annotated': am[0] if am else '',
                      'annotation': acls, 'phylo_label': lab, 'lost_at': lost})
        md_rows.append([loc, k.get('tier', ''), f"{k.get('subclass', '')}", len(car),
                        f'{len(disc_car)}/{len(car)}' + (f' ({len(no_file)} pending)' if no_file else ''),
                        ','.join(disc_other) or '0', mark(stages[1][1]), mark(stages[2][1]), su, sf, acls, lab or '-', lost])
        detail.append((loc, disc_names, missed, other_names))

    # ---- potentially more: final calls matching no known locus
    new = []
    for loc, gts in (fin or {}).items():
        if loc in matched_final:
            continue
        n_car, inf, cc = carrier_class(gts, args.germline_frac)
        if cc in ('germline-like', 'no carrier'):
            continue
        a = ann.get(loc, {})
        lab = labels.get(loc, '')
        carriers = [c for c, g in zip(fin_cols, gts) if g in CARRIER]
        new.append({'locus': loc, 'kind': Locus(loc).kind, 'carriers': n_car, 'informative': inf, 'carrier_class': cc,
                    'class': a.get('class', ''), 'tprt_call': a.get('tprt_call', ''), 'tprt_score': a.get('tprt_score', ''),
                    'phylo_label': lab, 'phylo_bucket': bucket_label(lab) if lab else '',
                    'carrier_ids': ','.join(carriers)})
    new.sort(key=lambda r: (r['phylo_bucket'] != 'consistent', r['class'] not in RTE_CLASSES,
                            r['carrier_class'] != 'shared', -int(r['carriers'] or 0)))

    write_tsv(os.path.join(args.out_dir, 'known_trace.tsv'), list(trace[0].keys()), [list(t.values()) for t in trace])
    if new:
        write_tsv(os.path.join(args.out_dir, 'new_calls.tsv'), list(new[0].keys()), [list(r.values()) for r in new])
    n_known = [t for t in trace if t['tier'] == 'known']
    n_found = sum(1 for t in n_known if t['lost_at'] == 'found')
    shown = [r for r in new if r['phylo_bucket'] == 'consistent' or r['class'] in RTE_CLASSES][:40]
    out = [f'# {P}{" arm " + args.label if args.label else ""}: known insertions',
           '', f'known set: {args.known}  (tolerance ±{args.tol} bp)', '',
           f'**{n_found}/{len(n_known)} known insertions reach the annotated calls.** '
           f'Discovery files: {len(per_colony)} colonies.', '',
           md_table(['known locus', 'tier', 'element', 'carriers', 'disc carriers', 'disc non-carriers', 'combined',
                     'contract', 'genotyped (unfiltered)', 'final calls', 'annotation', 'phylo', 'lost at / status'], md_rows),
           '', '## Discovery detail', '',
           *[f'- **{loc}**: carrier names {", ".join(dn) or "-"}'
             + (f'; missed in {"; ".join(ms)}' if ms else '')
             + (f'; non-carrier names {", ".join(on)}' if on else '') for loc, dn, ms, on in detail],
           '', f'## Potentially more: {len(new)} non-germline final calls match no known locus', '',
           f'{sum(1 for r in new if r["phylo_bucket"] == "consistent")} phylo-consistent, '
           f'{sum(1 for r in new if r["class"] in RTE_CLASSES)} RTE-classed, '
           f'{sum(1 for r in new if r["carrier_class"] == "private")} private. Shown: phylo-consistent or RTE-classed '
           f'(first 40); all in new_calls.tsv.', '']
    if shown:
        out.append(md_table(['locus', 'kind', 'carriers', 'class', 'tprt_call', 'tprt_score', 'phylo', 'carrier ids'],
                            [[r['locus'], r['kind'], r['carriers'], r['class'], r['tprt_call'], r['tprt_score'],
                              r['phylo_label'], r['carrier_ids'][:120]] for r in shown]))
    path = os.path.join(args.out_dir, 'known_report.md')
    with open(path, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    print('\n'.join(out))
    print(f'wrote {path}')


if __name__ == '__main__':
    main()
