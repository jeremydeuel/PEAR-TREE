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
import pandas as pd
from config import CONFIG
from multiprocessing import Pool
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
        stem = os.path.basename(file[:-7])
        insertions = {}
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
    print(f"{stem}: imported {len(insertions)} insertions.")
    print(d)

    # remove all insertions without at least 20 wild-types
    n_wt = (d == "wild-type").sum(axis=1)
    n_insertions = (d == 'homozygous').sum(axis=1) + (d == 'heterozygous').sum(axis=1)
    n_uncertain_insertion = (d == 'insertion?').sum(axis=1)
    n_uncertain = (d == 'wild-type?').sum(axis=1) + n_uncertain_insertion
    too_many_artefacts = (d == 'artefact').sum(axis=1) > CONFIG['combine_genotypes']['max_artefact']
    too_many_nas = d.isna().sum(axis=1) > CONFIG['combine_genotypes']['max_na']
    best_wt_score = support_score[d=="wild-type"].max(axis=1)
    best_ins_score = support_score[ ( d == "heterozygous" ) | ( d == "homozygous" )].max(axis=1)
    print(f"filtering strategy, starting with {d.shape[0]} insertions")
    print(f"- removing {sum(n_wt < CONFIG['combine_genotypes']['min_wild-types'])} insertions without at least {CONFIG['combine_genotypes']['min_wild-types']} wild-type colonies")
    print(f"- removing {sum(n_insertions < CONFIG['combine_genotypes']['min_insertions'])} insertions without at least {CONFIG['combine_genotypes']['min_insertions']} certain het or hom colony")
    print(f"- removing {sum(too_many_artefacts)} insertions with more than {CONFIG['combine_genotypes']['max_artefact']} artefact colonies")
    print(f"- removing {sum(too_many_nas)} insertions with more than {CONFIG['combine_genotypes']['max_na']} NA colonies")
    print(f"- removing {sum(n_uncertain > n_wt + n_insertions)} insertions with more than half uncertain calls.")
    print(f"- removing {sum(n_uncertain_insertion > n_insertions + 1)} insertions with more uncertain than certain insertion calls.")
    print(f"- removing {sum(best_wt_score<5000)} with a wt score of 10000 or less")
    print(f"- removing {sum(best_ins_score<5000)} with a het/hom score of 10000 or less")

    summary_filtering = pd.DataFrame([n_wt < CONFIG['combine_genotypes']['min_wild-types'],
                                      n_insertions < CONFIG['combine_genotypes']['min_insertions'],
                                      too_many_artefacts,
                                      too_many_nas,
                                      n_uncertain > n_wt + n_insertions,
                                      n_uncertain_insertion > n_insertions + 1,
                                      best_wt_score<1000,
                                      best_ins_score<1000
                                      ])
    summary_filtering = summary_filtering.any(axis=0)
    print(f"= removing {sum(summary_filtering)} insertions failing any of these tests.")
    d = d.loc[~summary_filtering]
    print(f" applying score filtering")
    print(f"writing a final of {d.shape[0]} filtered insertions")
    d.to_csv(output_file, sep=";")

