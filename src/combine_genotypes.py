# PEAR-TREE - paired ends of aberrant retrotransposons in phylogenetic trees
#
# Copyright (C) 2025 Jeremy Deuel <jeremy.deuel@usz.ch>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.


# sirt v2.0, step 4 collect genotypes
# usage: collect_genotype.py input_folder output_file.csv

import gzip
import os
import re
import numpy as np
import pandas as pd
from config import CONFIG
from multiprocessing import Pool


def _sample_stem(path: str) -> str:
    """Sample name from a genotype file path, stripping the full compound suffix.

    ``file[:-6]`` assumed a 6-char extension and left a trailing '.' on every
    ``.txt.gz`` / ``.csv.gz`` name (a 7-char suffix), corrupting the matrix
    column names. Strip the real suffix instead.
    """
    base = os.path.basename(path)
    return re.sub(r'(\.genotypes)?(\.(txt|csv))?(\.gz)?$', '', base)


ONESIDE_TOKEN = 'oneside_'
LOCUS_KINDS = ('TSD', 'TSD_DELETION', 'BLUNT', 'L1_MED_DELETION', 'L1_MED_DUPLICATION',
               'ONE_SIDED', 'OTHER')


def locus_kind(name: str) -> str:
    """Geometry class of a locus name (plans/tprt_hallmarks/SPEC.md "Pairing modes and locus
    names"), from gap = R - L: `contig:L-oneside_L` / `contig:oneside_R-R` -> ONE_SIDED;
    gap < -30 -> L1_MED_DELETION; -30..-1 -> TSD_DELETION; 0..1 -> BLUNT; 2..40 -> TSD;
    > 40 -> L1_MED_DUPLICATION (or a long-TSD artefact: same geometry). Anything else (e.g. a
    legacy `polyA_` / Feature-A `disc_` token) -> OTHER. Split from the right: contig names
    may contain ':' or '-'.

    Only used for the per-kind summary; the matrix rows are keyed by the name verbatim, which
    is what annotate (`<patient>.genotypes.csv.gz`) joins on."""
    contig, sep, pos = str(name).rpartition(':')
    left, sep2, right = pos.rpartition('-')
    if not sep or not sep2:
        return 'OTHER'
    if left.startswith(ONESIDE_TOKEN) != right.startswith(ONESIDE_TOKEN):
        return 'ONE_SIDED'
    try:
        gap = int(right) - int(left)
    except ValueError:
        return 'OTHER'
    if gap < -30:
        return 'L1_MED_DELETION'
    if gap < 0:
        return 'TSD_DELETION'
    if gap <= 1:
        return 'BLUNT'
    if gap <= 40:
        return 'TSD'
    return 'L1_MED_DUPLICATION'


def kind_summary(index, failed) -> pd.DataFrame:
    """Per locus kind: loci in, removed by any gate, kept (printed by collect_genotype)."""
    kinds = pd.Series([locus_kind(n) for n in index], index=index)
    failed = pd.Series(failed, index=index).astype(bool)
    out = pd.DataFrame({'loci': kinds.value_counts(),
                        'removed': kinds[failed].value_counts()}).fillna(0).astype(int)
    out['kept'] = out['loci'] - out['removed']
    return out.reindex([k for k in LOCUS_KINDS if k in out.index])


def _overdispersion(n_alt: pd.DataFrame, n_ref: pd.DataFrame) -> pd.Series:
    """Pearson chi-square dispersion of a locus's per-colony alt counts about ONE
    constant allele fraction, averaged over its degrees of freedom.

    Colonies are related by descent, so a real clonal insertion -- somatic or
    embryonic -- is present in some colonies and absent from the rest. Its allele
    fraction therefore genuinely DIFFERS between colonies (~0.5 where carried, 0.0
    elsewhere) and the alt counts scatter far wider than any single binomial: the
    statistic runs ~4-5 at this study's depths. An artefact whose alt reads come from
    a fixed, sequence-determined error rate has the SAME fraction in every colony, so
    its counts ARE one binomial and the statistic sits at ~1.

    This is a property of clonality rather than a tuned threshold: dispersion ~1 says
    "no colony differs from any other", which no real clonal event can satisfy. It is
    also what a wild-type-count gate cannot see -- `min_wild-types` removes a locus
    that is het everywhere (p~0.5) and keeps one that is ambiguous everywhere
    (p~0.05-0.30), i.e. it band-passes on the allele fraction and so selects exactly
    the loci that are too weak to call.

    Returns dispersion per locus, NaN where undefined: fewer than two informative
    colonies, or a degenerate pooled fraction of exactly 0 or 1 (no alt reads at all,
    or no ref reads at all -- both already removed by min_insertions/min_wild-types).
    """
    info = n_alt + n_ref
    usable = info > 0
    n_used = usable.sum(axis=1)
    a = n_alt.where(usable)
    i = info.where(usable)
    p = a.sum(axis=1) / i.sum(axis=1)
    # p in {0,1} => binomial variance 0 => the statistic is 0/0, not "no dispersion".
    degenerate = (p <= 0) | (p >= 1) | (n_used < 2)
    p_safe = p.where(~degenerate)
    expected = i.mul(p_safe, axis=0)
    variance = i.mul(p_safe * (1 - p_safe), axis=0)
    chi2 = (((a - expected) ** 2) / variance).sum(axis=1)
    return (chi2 / (n_used - 1)).mask(degenerate)


def mc_import(input_path, stem):
    corrected_artefacts = 0
    insertions = {}
    support_score = {}
    alt_score = {}
    n_alt = {}
    n_ref = {}
    with gzip.open(input_path, 'rt') as input_file:
        # The genotyper emits `insertion genotype score_genotype score_alternative
        # coverage n_alt n_ref n_art`. Locate the count columns by name rather than
        # position: older genotype files predate them, and the dispersion gate is
        # skipped rather than crashing when they are absent.
        header = input_file.readline().strip().split("\t")
        i_alt = header.index('n_alt') if 'n_alt' in header else None
        i_ref = header.index('n_ref') if 'n_ref' in header else None
        for line in input_file:
            line = line.strip()
            if not len(line): continue
            line = line.split("\t")
            gt = line[1]
            gt_support = int(line[2])
            alt_support = int(line[3])
            insertions[line[0]] = line[1]
            support_score[line[0]] = gt_support
            alt_score[line[0]] = alt_support
            if i_alt is not None and i_ref is not None:
                n_alt[line[0]] = int(line[i_alt])
                n_ref[line[0]] = int(line[i_ref])

    d1 = pd.DataFrame({stem: insertions.values()}, index=insertions.keys())
    d2 = pd.DataFrame({stem: support_score.values()}, index=support_score.keys(), dtype=int)
    d3 = pd.DataFrame({stem: alt_score.values()}, index=alt_score.keys(), dtype=int)
    d4 = (pd.DataFrame({stem: n_alt.values()}, index=n_alt.keys(), dtype=int)
          if n_alt else None)
    d5 = (pd.DataFrame({stem: n_ref.values()}, index=n_ref.keys(), dtype=int)
          if n_ref else None)
    print(f"completed {stem} with {d1.shape[0]} insertions.")
    return(d1,d2,d3,d4,d5)

def unfiltered_path(output_file):
    """`<P>.genotypes.csv.gz` -> `<P>.genotypes.unfiltered.csv.gz`."""
    for ext in ('.csv.gz', '.csv'):
        if output_file.endswith(ext):
            return output_file[:-len(ext)] + '.unfiltered' + ext
    return output_file + '.unfiltered'


def collect_genotype(input_files, output_file, threads):
    d = None
    pool = Pool(threads)
    results = []
    for file in input_files:
        stem = _sample_stem(file)
        results.append(pool.apply_async(mc_import, args=(file, stem)))
    print(f"collecting results...")
    insertions = []
    support_score = []
    alt_score = []
    n_alt_l = []
    n_ref_l = []
    for i in results:
        i, s, a, na, nr = i.get()
        insertions.append(i)
        support_score.append(s)
        alt_score.append(a)
        if na is not None and nr is not None:
            n_alt_l.append(na)
            n_ref_l.append(nr)
    print(f"concatenating")
    d = pd.concat(insertions, copy=False, axis=1)
    support_score = pd.concat(support_score, copy=False, axis=1)
    alt_score = pd.concat(alt_score, copy=False, axis=1)
    have_counts = len(n_alt_l) == len(insertions) and len(n_alt_l) > 0
    n_alt = pd.concat(n_alt_l, copy=False, axis=1) if have_counts else None
    n_ref = pd.concat(n_ref_l, copy=False, axis=1) if have_counts else None
    print(f"imported {d.shape[1]} samples covering {d.shape[0]} distinct insertions.")
    print(d)

    cg = CONFIG['combine_genotypes']
    # genotypes that mean "not assessable in this sample" and are counted as NA.
    # These strings must match the GT_* vocabulary emitted by genotyping_insertion.py.
    NA_LIKE = ('high-coverage', 'no-coverage', 'error')
    min_best_score = cg.get('min_best_score', 800)

    n_wt = (d == "wild-type").sum(axis=1)
    # a colony carries the insertion if it is a confident het, hom, OR a zygosity-unclear
    # 'insertion' call (presence certain, hom/het indeterminate) -- the last is essential
    # for clade detection across many colonies, where presence, not zygosity, defines a clade.
    n_insertions = ((d == 'homozygous').sum(axis=1) + (d == 'heterozygous').sum(axis=1)
                    + (d == 'insertion').sum(axis=1))
    # low-confidence calls now actually emitted by the genotyper (GT_*_UNCERTAIN).
    n_uncertain_insertion = (d == 'insertion?').sum(axis=1)
    n_uncertain = (d == 'wild-type?').sum(axis=1) + n_uncertain_insertion
    too_many_artefacts = (d == 'artefact').sum(axis=1) > cg['max_artefact']
    n_na = d.isna().sum(axis=1)
    for s in NA_LIKE:
        n_na = n_na + (d == s).sum(axis=1)
    too_many_nas = n_na > cg['max_na']
    # max() over a row with no het/hom sample is NaN; treat that as failing the
    # score gate rather than silently passing it (NaN < x is False).
    best_ins_score = support_score[(d == "heterozygous") | (d == "homozygous") | (d == "insertion")].max(axis=1).fillna(-np.inf)

    # Clonality gate. Every gate above counts CALLS, and a call is a thresholded
    # allele fraction, so none of them can tell a locus that is genuinely present in
    # some colonies from one whose alt reads are a constant per-locus error rate that
    # the VAF bands slice at random. This one tests the read counts directly: a real
    # clonal event must be over-dispersed relative to a single binomial. Set
    # min_dispersion to 0/None to disable. Requires n_alt/n_ref in the genotype files.
    min_dispersion = cg.get('min_dispersion', 3.0)
    if have_counts and min_dispersion:
        dispersion = _overdispersion(n_alt.reindex(index=d.index, columns=d.columns),
                                     n_ref.reindex(index=d.index, columns=d.columns))
        # NaN => undefined (degenerate p, or too few informative colonies) => fails.
        underdispersed = ~(dispersion >= min_dispersion)
    else:
        if not have_counts:
            print("NOTE: genotype files carry no n_alt/n_ref columns -- the "
                  "dispersion gate is skipped.")
        dispersion = pd.Series(np.nan, index=d.index)
        underdispersed = pd.Series(False, index=d.index)

    print(f"filtering strategy, starting with {d.shape[0]} insertions")
    print(f"- removing {int((n_wt < cg['min_wild-types']).sum())} insertions without at least {cg['min_wild-types']} wild-type colonies")
    print(f"- removing {int((n_insertions < cg['min_insertions']).sum())} insertions without at least {cg['min_insertions']} certain het or hom colony")
    print(f"- removing {int(too_many_artefacts.sum())} insertions with more than {cg['max_artefact']} artefact colonies")
    print(f"- removing {int(too_many_nas.sum())} insertions with more than {cg['max_na']} NA colonies")
    print(f"- removing {int((n_uncertain > n_wt + n_insertions).sum())} insertions with more than half uncertain calls.")
    print(f"- removing {int((n_uncertain_insertion > n_insertions + 1).sum())} insertions with more uncertain than certain insertion calls.")
    print(f"- removing {int((best_ins_score < min_best_score).sum())} with a het/hom score below {min_best_score}")
    if have_counts and min_dispersion:
        print(f"- removing {int(underdispersed.sum())} insertions whose per-colony alt "
              f"counts are one binomial (dispersion < {min_dispersion}): the same "
              f"allele fraction in every colony, so not a clonal event")

    summary_filtering = pd.DataFrame([n_wt < cg['min_wild-types'],
                                      n_insertions < cg['min_insertions'],
                                      too_many_artefacts,
                                      too_many_nas,
                                      n_uncertain > n_wt + n_insertions,
                                      n_uncertain_insertion > n_insertions + 1,
                                      best_ins_score < min_best_score,
                                      underdispersed,
                                      ])
    summary_filtering = summary_filtering.any(axis=0)
    print(f"= removing {sum(summary_filtering)} insertions failing any of these tests.")
    print("per locus kind (TPRT pairing modes; name geometry, gap = R - L):")
    print(kind_summary(d.index, summary_filtering.reindex(d.index).values).to_string())
    if cg.get('write_unfiltered'):
        # every locus, before the gates: what annotate scores when the gates are being evaluated
        # rather than applied (TPRT A/B kit, cluster/tprt/arm_config.py)
        unf = unfiltered_path(output_file)
        d.to_csv(unf, sep=";")
        print(f"wrote all {d.shape[0]} insertions (before the gates) to {unf}")
    d = d.loc[~summary_filtering]
    print(f" applying score filtering")
    print(f"writing a final of {d.shape[0]} filtered insertions")
    d.to_csv(output_file, sep=";")

