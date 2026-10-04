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

"""Extend a genotyping contract with the TPRT one-sided loci.

combine_insertions keeps one-sided loci (`contig:L-oneside_L` / `contig:oneside_R-R`) in
`<stem>.combined.txt.gz` but leaves them out of `<stem>.genotyping.txt.gz`. This writes
`<out>` = the contract unchanged, followed by one entry per surviving one-sided locus that
carries ONLY its real side (`@LEFT_*` for `L-oneside_L`, `@RIGHT_*` for `oneside_R-R`),
built exactly like combine builds that side for a two-sided locus:

  LEFT_INSERTION  = the max_bases clip bases next to the junction (reference-forward)
  LEFT_REFERENCE  = genome[L - max_bases, L)
  RIGHT_INSERTION = the first max_bases clip bases after the junction
  RIGHT_REFERENCE = genome[R, R + max_bases)

and with combine's per-side exclusions (long contig names, missing/N reference, a clip that
is indistinguishable from the reference). The Rust genotyper scores such a locus on its real
side only when `one_sided_loci = true` (cluster/config.genotype.grch38.tprt).

The clip consensus is read back from combined.txt.gz, whose `<name>:L` record is
`lower(clip) + UPPER(aligned)` and `<name>:R` record `UPPER(aligned) + lower(clip)`.

usage: python src/genotyping_contract_oneside.py --contract <stem>.genotyping.txt.gz \
           --combined <stem>.combined.txt.gz --out <stem>.genotyping.tprt.txt.gz [--config-dir DIR]
"""

import argparse
import gzip
import os
import sys

ONESIDE_TOKEN = 'oneside_'


def read_combined_sides(path):
    """`{title: {'L': seq, 'R': seq}}` (insertion order) from a combined.txt.gz fastq."""
    sides = {}
    with gzip.open(path, 'rt') as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            header = header.strip()
            if not header:
                continue
            seq = fh.readline().strip()
            fh.readline()
            fh.readline()
            title, side = header[1:-2], header[-1]
            sides.setdefault(title, {})[side] = seq
    return sides


def one_sided_side(title):
    """('LEFT'|'RIGHT', contig, pos) of a one-sided locus name's REAL side, else None."""
    contig, _, pos = title.rpartition(':')
    left, _, right = pos.rpartition('-')
    if right.startswith(ONESIDE_TOKEN) and not left.startswith(ONESIDE_TOKEN):
        return 'LEFT', contig, int(left)
    if left.startswith(ONESIDE_TOKEN) and not right.startswith(ONESIDE_TOKEN):
        return 'RIGHT', contig, int(right)
    return None


def _clip(seq, side):
    """Lower-case clip part of a combined consensus (prefix for L, suffix for R)."""
    if side == 'L':
        n = len(seq) - len(seq.lstrip('acgtn'))
        return seq[:n]
    n = len(seq) - len(seq.rstrip('acgtn'))
    return seq[len(seq) - n:]


def build_entry(title, sides, get_sequence, max_bases, matching_score):
    """(contract text, None) for one one-sided locus, or (None, exclusion reason)."""
    real = one_sided_side(title)
    if real is None:
        return None, 'not one-sided'
    side, contig, pos = real
    if len(contig) > 5:
        return None, 'contig'
    key = 'L' if side == 'LEFT' else 'R'
    seq = sides.get(key)
    if not seq:
        return None, 'no consensus'
    clip = _clip(seq, key).upper()
    if not clip:
        return None, 'no clip'
    if side == 'LEFT':
        alt = clip[-max_bases:]
        ref = get_sequence(contig, pos - max_bases, pos).upper()
    else:
        alt = clip[:max_bases]
        ref = get_sequence(contig, pos, pos + max_bases).upper()
    if not len(ref):
        return None, 'missing reference'
    if 'N' in ref:
        return None, 'N in reference'
    if matching_score([ref, alt]) > 0:
        return None, 'clip similar to reference'
    return f'>{title}\n@{side}_INSERTION\n{alt}\n@{side}_REFERENCE\n{ref}\n', None


def extend_contract(contract, combined, out, get_sequence, max_bases, matching_score):
    """Write `out` = `contract` + one-sided entries. Returns (n_added, {reason: n_excluded})."""
    sides = read_combined_sides(combined)
    with gzip.open(contract, 'rt') as fh:
        base = fh.read()
    present = {line[1:].strip() for line in base.splitlines() if line.startswith('>')}
    added, excluded = [], {}
    for title, s in sides.items():
        if title in present or one_sided_side(title) is None:
            continue
        text, why = build_entry(title, s, get_sequence, max_bases, matching_score)
        if text is None:
            excluded[why] = excluded.get(why, 0) + 1
        else:
            added.append(text)
    tmp = f'{out}.tmp.{os.getpid()}'
    with gzip.open(tmp, 'wt') as fh:
        fh.write(base)
        if base and not base.endswith('\n'):
            fh.write('\n')
        fh.writelines(added)
    os.replace(tmp, out)
    return len(added), excluded


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--contract', required=True)
    ap.add_argument('--combined', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--config-dir', help='directory holding the config.py to use (default: src/)')
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    if a.config_dir:
        sys.path.insert(0, os.path.abspath(a.config_dir))
    sys.path.insert(1, here)
    from config import CONFIG
    from combine_insertions_get_sequence import get_sequence
    from sequence_checks import sequence_matching_score
    n, excluded = extend_contract(a.contract, a.combined, a.out, get_sequence,
                                  CONFIG['genotyping']['max_bases'], sequence_matching_score)
    print(f'one-sided loci: added {n} to {a.out}; excluded {sum(excluded.values())} {excluded}')


if __name__ == '__main__':
    main()
