# Geminate - Genomic events of mutation through insertional alterations by transposable elements
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


import gzip
import pandas as pd
import os
from collections import Counter
if __name__ == '__main__':
    d = pd.read_excel("all.insertions.xlsx")
    d = d.loc[d['wt'] < 50]
    files = list(os.listdir('PD45517'))+list(os.listdir('PD45534'))+list(os.listdir('PD48402'))
    files = [file[:-8] for file in files if len(file)>6]
    files_per_mouse = Counter([file[:7] for file in files])
    print(files_per_mouse)
    print(files_per_mouse)
    print(d)
    insertions = list(d['insertion'])
    output = {
        'insertion': [],
        'mouse': [],
        'colonies': [],
        'in_first_step': []
    }
    for patient in ('PD45517','PD45534','PD48402'):
        with gzip.open(f'{patient}.insertions.combined.txt.gz', 'rt') as fh:
            insertion = None
            state = None
            files = []
            for line in fh:
                line = line.strip()
                if not line: continue
                if line[0] == '>':
                    files = Counter([file[:7] for file in files])
                    for k, v in files.items():
                        print(f"{insertion}: found {v} of {files_per_mouse[k]}")
                        output['insertion'].append(insertion)
                        output['mouse'].append(k)
                        output['colonies'].append(files_per_mouse[k])
                        output['in_first_step'].append(v)
                    files = []
                    insertion = line[1:]
                    if insertion in insertions:
                        #print(f">{insertion}")
                        state = None
                    else:
                        insertion = None
                if insertion is None:
                    continue
                if line[0] == '@':
                    state = line[1:]
                    continue
                if state == 'FILES':
                    files.append(os.path.basename(line[:-13]))

    pd.DataFrame(output).to_excel('first.step.stats.xlsx')
