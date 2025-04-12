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
    with gzip.open(input_path, 'rt') as input_file:
        input_file.readline()  # ignore
        for line in input_file:
            line = line.strip()
            if not len(line): continue
            line = line.split("\t")
            gt = line[1]
            wt = int(line[4])
            if gt == "artefact":  # switched off this part
                right_ins = int(line[2])
                left_ins = int(line[3])
                good = right_ins + left_ins + wt
                art_left = int(line[5])
                art_right = int(line[6])
                art_both = int(line[7])
                bad = art_left + art_right + art_both
                ins_both = int(line[8])
                if not ins_both and not bad:
                    if right_ins or left_ins:
                        if right_ins + left_ins > wt / 15:
                            gt = 'insertion?'
                            corrected_artefacts += 1
            if gt == 'wild-type' and wt < 5:
                gt = 'wild-type?'
            insertions[line[0]] = gt
    if not len(insertions):
        print(f"ignoring {stem}: no insertions found.")
        return None
    d1 = pd.DataFrame({stem: insertions.values()}, index=insertions.keys())
    print(f"completed {stem} with {d1.shape[0]} insertions.")
    return(d1)

def collect_genotype(input_files, output_file, threads):
    d = None
    pool = Pool(threads)
    results = []
    for file in input_files:
        stem = os.path.basename(file[:-7])
        insertions = {}
        results.append(pool.apply_async(mc_import, args=(file, stem)))
    print(f"collecting results...")
    results = [i.get() for i in results]
    print(f"concatenating")
    d = pd.concat(results, copy=False, axis=1)
    print(f"{stem}: imported {len(insertions)} insertions.")
    print(d)

    # remove all insertions without at least 20 wild-types
    n_wt = (d == "wild-type").sum(axis=1)
    n_insertions = (d == 'homozygous').sum(axis=1) + (d == 'heterozygous').sum(axis=1)
    n_uncertain_insertion = (d == 'insertion?').sum(axis=1)
    n_uncertain = (d == 'wild-type?').sum(axis=1) + n_uncertain_insertion
    too_many_artefacts = (d == 'artefact').sum(axis=1) > CONFIG['combine_genotypes']['max_artefact']
    too_many_nas = d.isna().sum(axis=1) > CONFIG['combine_genotypes']['max_na']

    print(f"filtering strategy, starting with {d.shape[0]} insertions")
    print(f"- removing {sum(n_wt < CONFIG['combine_genotypes']['min_wild-types'])} insertions without at least {CONFIG['combine_genotypes']['min_wild-types']} wild-type colonies")
    print(f"- removing {sum(n_insertions < CONFIG['combine_genotypes']['min_insertions'])} insertions without at least {CONFIG['combine_genotypes']['min_insertions']} certain het or hom colony")
    print(f"- removing {sum(too_many_artefacts)} insertions with more than {CONFIG['combine_genotypes']['max_artefact']} artefact colonies")
    print(f"- removing {sum(too_many_nas)} insertions with more than {CONFIG['combine_genotypes']['max_na']} NA colonies")
    print(f"- removing {sum(n_uncertain > n_wt + n_insertions)} insertions with more than half uncertain calls.")
    print(f"- removing {sum(n_uncertain_insertion > n_insertions + 1)} insertions with more uncertain than certain insertion calls.")
    summary_filtering = pd.DataFrame([n_wt < CONFIG['combine_genotypes']['min_wild-types'],
                                      n_insertions < CONFIG['combine_genotypes']['min_insertions'],
                                      too_many_artefacts,
                                      too_many_nas,
                                      n_uncertain > n_wt + n_insertions,
                                      n_uncertain_insertion > n_insertions + 1])
    summary_filtering = summary_filtering.any(axis=0)
    print(f"= removing {sum(summary_filtering)} insertions failing any of these tests.")
    d = d.loc[~summary_filtering]
    print(f"writing a final of {d.shape[0]} filtered insertions")
    d.to_csv(output_file, sep=";")

