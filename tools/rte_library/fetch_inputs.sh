#!/usr/bin/env bash
# Download every input of tools/rte_library/build.py into one directory (~2.2 GB, mostly the
# two 2bit genomes and the RepeatMasker tracks). Usage:
#   bash tools/rte_library/fetch_inputs.sh /absolute/path/to/rte_inputs
#   python tools/rte_library/build.py --inputs /absolute/path/to/rte_inputs \
#       --work /absolute/path/to/rte_work --out resources/rte_library
#
# L1Base 2 has no stable file URL: its "export all" page for dataset hsflil1_8438 is fetched
# below; if the server layout changed, export BED + FASTA by hand from
# https://l1base.charite.de/ (dataset "hsflil1_8438", human FLI-L1, GRCh38) into libs/l1base/.
set -euo pipefail
D="${1:?usage: fetch_inputs.sh OUTDIR}"
mkdir -p "$D"/genomes "$D"/libs/l1base "$D"/supp "$D"/ncbi
UCSC=https://hgdownload.soe.ucsc.edu/goldenPath
get() { [ -s "$2" ] || curl -fsSL --retry 3 -o "$2" "$1"; }

# --- genomes, repeat tracks, chains, cytobands (UCSC)
get $UCSC/hg38/bigZips/hg38.2bit                       "$D"/genomes/hg38.2bit
get $UCSC/hs1/bigZips/hs1.2bit                         "$D"/genomes/hs1.2bit
get $UCSC/hg38/database/rmsk.txt.gz                    "$D"/genomes/hg38.rmsk.txt.gz
get $UCSC/hs1/bigZips/hs1.repeatMasker.out.gz          "$D"/genomes/hs1.repeatMasker.out.gz
get $UCSC/hg38/liftOver/hg38ToHs1.over.chain.gz        "$D"/genomes/hg38ToHs1.over.chain.gz
get $UCSC/hg19/liftOver/hg19ToHg38.over.chain.gz       "$D"/genomes/hg19ToHg38.over.chain.gz
get $UCSC/hg38/database/cytoBand.txt.gz                "$D"/genomes/hg38.cytoBand.txt.gz

# --- L1Base 2 (Penzkofer et al. 2017 NAR 45:D68): FLI-L1 = full-length, both ORFs intact
get "https://l1base.charite.de/exportall.php?DBN=hsflil1_8438&TYPE=fasta" "$D"/libs/l1base/hsflil1_8438.fa
get "https://l1base.charite.de/exportall.php?DBN=hsflil1_8438&TYPE=bed"   "$D"/libs/l1base/hsflil1_8438.bed

# --- Dfam consensus sequences (CC0), young human families used as cross-check
: > "$D"/libs/dfam_consensus.fa.tmp
for acc in DF000000002:AluY DF000000053:AluYa5 DF000000055:AluYb8 DF000000056:AluYb9 \
           DF000000057:AluYc DF000001240:AluYe5 DF000001317:AluYf1 DF000000634:AluYg6 \
           DF000001318:AluYh3 DF000000063:AluYh9 DF000000066:AluYk4 DF000001154:AluYm1 \
           DF000000047:AluSx DF000000007:AluJb \
           DF000000225:L1HS_3end DF000000226:L1HS_5end DF000000339:L1PA2_3end \
           DF000000340:L1PA3_3end DF000000341:L1PA4_3end DF000000342:L1PA5_3end \
           DF000000343:L1PA6_3end DF000000344:L1PA7_3end \
           DF000001067:SVA_A DF000001068:SVA_B DF000001069:SVA_C DF000001070:SVA_D \
           DF000001071:SVA_E DF000001072:SVA_F DF003894004:SVA2; do
  id=${acc%%:*}; name=${acc##*:}
  seq=$(curl -fsSL "https://dfam.org/api/families/$id/sequence?format=fasta" | grep -v '^>' | tr -d '\n')
  printf ">%s|%s|len=%d\n%s\n" "$name" "$id" "${#seq}" "$seq" >> "$D"/libs/dfam_consensus.fa.tmp
done
mv "$D"/libs/dfam_consensus.fa.tmp "$D"/libs/dfam_consensus.fa

# --- reference L1 sequences (NCBI E-utilities): L1.3, L1.2, L1RP
for a in L19088 M80343 AF148856; do
  get "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$a&rettype=gb&retmode=text" "$D"/ncbi/$a.gb
  sleep 1
done

# --- published source-element tables
# Nam et al. 2023 Nature 617:540 (PMC10191854): Supp. Table 2 (MOESM5) and 4 (MOESM7)
for i in 5 7; do
  get "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06046-z/MediaObjects/41586_2023_6046_MOESM${i}_ESM.xlsx" \
      "$D"/supp/nam_MOESM${i}.xlsx
done
# Rodriguez-Martin et al. 2020 Nat Genet 52:306 (PMC7058536): Supp. Tables 1-8 (MOESM3)
get "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7058536/supplementaryFiles" "$D"/supp/rm_supp.zip
unzip -o -q "$D"/supp/rm_supp.zip 41588_2019_562_MOESM3_ESM.xlsx -d "$D"/supp
# Gardner et al. 2017 Genome Res 27:1916 (MELT; PMC5668948): Supplemental Table S9
get "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5668948/supplementaryFiles" "$D"/supp/melt_supp.zip
unzip -o -q "$D"/supp/melt_supp.zip supp_gr.218032.116_Supplemental_Table_S9.xlsx -d "$D"/supp
# Tubio et al. 2014 Science 345:1251343 Table S5 is NOT open access (Europe PMC: "not open
# access"; PMC serves the NIHMS supplement only behind a proof-of-work browser check). Its
# per-source transduction counts enter via Gardner 2017 S9 B/C ("Tubio et al. Activity").
echo "inputs ready in $D"
