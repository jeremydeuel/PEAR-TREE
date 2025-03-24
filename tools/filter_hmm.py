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

# filters hmm files for species
import argparse, sys
clm = argparse.ArgumentParser(
        prog="Filter HMM File for Species",
        description="This program filteres a HMM file for species. (c) by Jeremy Deuel <jeremy.deuel@usz.ch> 2024",
        usage="python filter_hmm.py [--in unfiltered file, e.g. from Dfam] [--out outfile.hmm] [--species \"Mus musculus\"]"
    )
clm.add_argument('-i','--infile', help="input HMM file", nargs=1)
clm.add_argument('-o','--outfile',help="filtered output HMM file",nargs=1)
clm.add_argument('-s','--species', help="species to filter for", nargs=1)
args = clm.parse_args(sys.argv[1:])
if not len(args.infile) or not len(args.outfile) or not len(args.species):
    clm.print_help()
    exit(0)
print(f"This is {clm.prog} with arguments {' '.join(sys.argv[1:])}")
species = f"TaxName:{args.species[0]}"
with open(args.infile[0], 'r') as in_file:
    with open(args.outfile[0],'w') as out_file:
        total_hmms = 0
        filtered_hmms = 0
        current_hmm = ""
        include_flag = False
        for l in in_file:
            if l == "//\n": #sign for next HMM
                if include_flag:
                    filtered_hmms += 1
                    current_hmm += l
                    out_file.write(current_hmm)
                current_hmm = ""
                include_flag = False
                total_hmms += 1
                if total_hmms%100 == 0:
                    print(f"Processed HMM {total_hmms}, included {filtered_hmms} ({(100*filtered_hmms/total_hmms):.2f}%), continuing...")
            else:
                if l[0:2] == "TH" and species in l:
                    include_flag = True
                current_hmm += l
        print(f"Parsed a total of {total_hmms}, thereof {filtered_hmms} passed filter.")
