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
def mc_import(input_path, stem):
    corrected_artefacts = 0
    insertions = {}
    support_score = {}
    alt_score = {}
    with gzip.open(input_path, 'rt') as input_file:
        input_file.readline()  # ignore
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

    d1 = pd.DataFrame({stem: insertions.values()}, index=insertions.keys())
    d2 = pd.DataFrame({stem: support_score.values()}, index=support_score.keys(), dtype=int)
    d3 = pd.DataFrame({stem: alt_score.values()}, index=alt_score.keys(), dtype=int)
    print(f"completed {stem} with {d1.shape[0]} insertions.")
    return(d1,d2,d3)

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
    for i in results:
        i, s, a = i.get()
        insertions.append(i)
        support_score.append(s)
        alt_score.append(a)
    print(f"concatenating")
    d = pd.concat(insertions, copy=False, axis=1)
    support_score = pd.concat(support_score, copy=False, axis=1)
    alt_score = pd.concat(alt_score, copy=False, axis=1)
    print(f"imported {d.shape[1]} samples covering {d.shape[0]} distinct insertions.")
    print(d)

    cg = CONFIG['combine_genotypes']
    # genotypes that mean "not assessable in this sample" and are counted as NA.
    # These strings must match the GT_* vocabulary emitted by genotyping_insertion.py.
    NA_LIKE = ('high-coverage', 'no-coverage', 'error')
    min_best_score = cg.get('min_best_score', 800)

    n_wt = (d == "wild-type").sum(axis=1)
    n_insertions = (d == 'homozygous').sum(axis=1) + (d == 'heterozygous').sum(axis=1)
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
    best_ins_score = support_score[(d == "heterozygous") | (d == "homozygous")].max(axis=1).fillna(-np.inf)
    print(f"filtering strategy, starting with {d.shape[0]} insertions")
    print(f"- removing {int((n_wt < cg['min_wild-types']).sum())} insertions without at least {cg['min_wild-types']} wild-type colonies")
    print(f"- removing {int((n_insertions < cg['min_insertions']).sum())} insertions without at least {cg['min_insertions']} certain het or hom colony")
    print(f"- removing {int(too_many_artefacts.sum())} insertions with more than {cg['max_artefact']} artefact colonies")
    print(f"- removing {int(too_many_nas.sum())} insertions with more than {cg['max_na']} NA colonies")
    print(f"- removing {int((n_uncertain > n_wt + n_insertions).sum())} insertions with more than half uncertain calls.")
    print(f"- removing {int((n_uncertain_insertion > n_insertions + 1).sum())} insertions with more uncertain than certain insertion calls.")
    print(f"- removing {int((best_ins_score < min_best_score).sum())} with a het/hom score below {min_best_score}")

    summary_filtering = pd.DataFrame([n_wt < cg['min_wild-types'],
                                      n_insertions < cg['min_insertions'],
                                      too_many_artefacts,
                                      too_many_nas,
                                      n_uncertain > n_wt + n_insertions,
                                      n_uncertain_insertion > n_insertions + 1,
                                      best_ins_score < min_best_score
                                      ])
    summary_filtering = summary_filtering.any(axis=0)
    print(f"= removing {sum(summary_filtering)} insertions failing any of these tests.")
    d = d.loc[~summary_filtering]
    print(f" applying score filtering")
    print(f"writing a final of {d.shape[0]} filtered insertions")
    d.to_csv(output_file, sep=";")

