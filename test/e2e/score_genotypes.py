#!/usr/bin/env python3
"""E2E: per-colony genotype concordance vs simulator truth, per locus kind.

Joins the E2E event table (`e2e_events.tsv`: event id -> combined locus name `comb_name`) with
the simulator truth (`results_by_event.tsv`: `vaf_by_sample`, one VAF per colony) and the
per-colony genotyper outputs (`<geno-dir>/S<k>.txt.gz`). A colony's call is

  present  heterozygous / homozygous / insertion   (what combine_genotypes counts as carrier)
  absent   wild-type
  nc       anything else (insertion?, wild-type?, artefact, no-/high-coverage, error) or the
           locus is not in the contract at all (`missing`)

and truth is present iff its VAF > 0. The locus kind is the geometry the genotyper sees in the
name (gap = R - L): one-sided, TSD 2-40, target-site deletion [-30,-1], blunt {0,1},
L1DEL (< -30), L1DUP (> 40).

With `--wt-sample S4` the genotype file `<geno-dir>/S4.txt.gz` of a wild-type control colony
(reads simulated from the event-free haplotype only) is scored too: every TP locus is truth
ABSENT there, so a `present` call is a false positive of the genotyper itself.

usage: score_genotypes.py --e2e-dir $SP/work/e2e --geno-dir DIR [--samples 3] [--label X]
                          [--wt-sample S4]
"""
import argparse
import collections
import gzip
import os

PRESENT = {'heterozygous', 'homozygous', 'insertion', 'present'}
ABSENT = {'wild-type', 'absent'}
# numeric (peartree-genotype2) files: posterior thresholds that define present / absent
V2_P_PRESENT = 0.9
V2_P_ABSENT = 0.8
KINDS = ['TSD 2-40', 'target-site deletion', 'blunt 0-1', 'L1DEL (< -30)', 'L1DUP (> 40)',
         'one-sided', 'other']


def locus_kind(name):
    pos = name.rpartition(':')[2]
    left, _, right = pos.rpartition('-')
    if left.startswith('oneside_') or right.startswith('oneside_'):
        return 'one-sided'
    try:
        gap = int(right) - int(left)
    except ValueError:
        return 'other'
    if gap < -30:
        return 'L1DEL (< -30)'
    if gap < 0:
        return 'target-site deletion'
    if gap <= 1:
        return 'blunt 0-1'
    if gap <= 40:
        return 'TSD 2-40'
    return 'L1DUP (> 40)'


def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        return [dict(zip(head, line.rstrip('\n').split('\t'))) for line in fh]


def read_calls(path):
    """Legacy files: the genotype string. Numeric peartree-genotype2 files (header starts with
    `locus\tkind\tstatus`): 'present' if P(het)+P(hom) >= V2_P_PRESENT, 'absent' if
    P(absent) >= V2_P_ABSENT, else the status word or 'p_present=<x>' (a no-call reason)."""
    calls = {}
    with gzip.open(path, 'rt') as fh:
        head = fh.readline().rstrip('\n').split('\t')
        v2 = head[:3] == ['locus', 'kind', 'status']
        col = {h: i for i, h in enumerate(head)}
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if not v2:
                calls[p[0]] = p[1]
                continue
            if p[col['status']] != 'ok':
                calls[p[0]] = p[col['status']]
                continue
            p_abs = float(p[col['p_absent']])
            p_pres = float(p[col['p_het']]) + float(p[col['p_hom']])
            if p_pres >= V2_P_PRESENT:
                calls[p[0]] = 'present'
            elif p_abs >= V2_P_ABSENT:
                calls[p[0]] = 'absent'
            else:
                calls[p[0]] = f'p_present={p_pres:.1f}'
    return calls


def classify(gt):
    if gt in PRESENT:
        return 'present'
    if gt in ABSENT:
        return 'absent'
    return 'nc'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--e2e-dir', required=True)
    ap.add_argument('--geno-dir', required=True)
    ap.add_argument('--samples', type=int, default=3)
    ap.add_argument('--label', default='')
    ap.add_argument('--wt-sample', default=None)
    a = ap.parse_args()
    events = {r['id']: r for r in read_tsv(os.path.join(a.e2e_dir, 'e2e_events.tsv'))}
    truth = {r['id']: r for r in read_tsv(os.path.join(a.e2e_dir, 'results_by_event.tsv'))}
    calls = [read_calls(os.path.join(a.geno_dir, f'S{k + 1}.txt.gz')) for k in range(a.samples)]

    # kind -> counters
    n_loci = collections.Counter()
    c = collections.defaultdict(collections.Counter)
    nc_why = collections.defaultdict(collections.Counter)
    seen = set()
    for eid, ev in events.items():
        name = ev.get('comb_name', 'None')
        if ev.get('role') != 'TP' or name in ('None', '') or name in seen:
            continue
        seen.add(name)
        vafs = [float(v) for v in truth[eid]['vaf_by_sample'].split(',')]
        kind = locus_kind(name)
        n_loci[kind] += 1
        for k in range(a.samples):
            gt = calls[k].get(name, 'missing')
            got = classify(gt)
            want = 'present' if vafs[k] > 0 else 'absent'
            if got == 'nc':
                c[kind][f'nc_{want}'] += 1
                nc_why[kind][gt] += 1
            elif got == want:
                c[kind][f'ok_{want}'] += 1
            else:
                c[kind][f'bad_{want}'] += 1

    wt = read_calls(os.path.join(a.geno_dir, f'{a.wt_sample}.txt.gz')) if a.wt_sample else None
    wt_c = collections.defaultdict(collections.Counter)
    if wt is not None:
        for name in seen:
            gt = wt.get(name, 'missing')
            wt_c[locus_kind(name)][classify(gt)] += 1
            if classify(gt) == 'nc':
                wt_c[locus_kind(name)]['nc:' + gt] += 1

    # unexplained combined loci (not matched to a TP event): carrier-colony histogram per kind
    unexpl = collections.defaultdict(collections.Counter)
    all_names = set().union(*[set(x) for x in calls])
    for name in all_names - seen:
        kind = locus_kind(name)
        k_present = sum(classify(calls[k].get(name, 'missing')) == 'present' for k in range(a.samples))
        unexpl[kind][k_present] += 1

    print(f'### Genotype concordance vs truth {a.label}'.rstrip())
    print()
    print('| locus kind | TP loci | colony calls | present ok | absent ok | present->absent | '
          'absent->present | no call (truth present / absent) | concordance of calls | '
          'no-call reasons |')
    print('|---|---|---|---|---|---|---|---|---|---|')
    tot = collections.Counter()
    for kind in KINDS:
        if not n_loci[kind]:
            continue
        x = c[kind]
        tot.update(x)
        tot['loci'] += n_loci[kind]
        called = x['ok_present'] + x['ok_absent'] + x['bad_present'] + x['bad_absent']
        conc = f"{(x['ok_present'] + x['ok_absent']) / called:.3f}" if called else '.'
        why = ', '.join(f'{g} {n}' for g, n in nc_why[kind].most_common())
        print(f"| {kind} | {n_loci[kind]} | {n_loci[kind] * a.samples} | {x['ok_present']} | "
              f"{x['ok_absent']} | {x['bad_present']} | {x['bad_absent']} | "
              f"{x['nc_present']} / {x['nc_absent']} | {conc} | {why or '.'} |")
    called = tot['ok_present'] + tot['ok_absent'] + tot['bad_present'] + tot['bad_absent']
    conc = f"{(tot['ok_present'] + tot['ok_absent']) / called:.3f}" if called else '.'
    print(f"| **all** | {tot['loci']} | {tot['loci'] * a.samples} | {tot['ok_present']} | "
          f"{tot['ok_absent']} | {tot['bad_present']} | {tot['bad_absent']} | "
          f"{tot['nc_present']} / {tot['nc_absent']} | {conc} | |")
    if wt is not None:
        print()
        print(f'Wild-type control colony {a.wt_sample} (every TP locus truth-absent):')
        print()
        print('| locus kind | TP loci | absent (wild-type) | present (false positive) | no call | '
              'no-call reasons |')
        print('|---|---|---|---|---|---|')
        for kind in KINDS:
            if not n_loci[kind]:
                continue
            x = wt_c[kind]
            why = ', '.join(f'{g[3:]} {n}' for g, n in x.most_common() if g.startswith('nc:'))
            print(f"| {kind} | {n_loci[kind]} | {x['absent']} | {x['present']} | {x['nc']} | "
                  f"{why or '.'} |")
    print()
    print(f'Unexplained contract loci (no TP event), colonies called present (0/1/.../{a.samples}):')
    print()
    print('| locus kind | ' + ' | '.join(str(k) for k in range(a.samples + 1)) + ' |')
    print('|---|' + '---|' * (a.samples + 1))
    for kind in KINDS:
        if unexpl[kind]:
            print(f'| {kind} | ' + ' | '.join(str(unexpl[kind][k]) for k in range(a.samples + 1)) + ' |')


if __name__ == '__main__':
    main()
