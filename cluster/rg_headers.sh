#!/bin/bash
# =============================================================================
# cluster/rg_headers.sh — one row per colony of a patient: the @RG / @PG facts of its BAM,
# read from the iRODS replica HEADER ONLY (no staging, no record reads — Sanger policy).
#
#   bash cluster/rg_headers.sh <PATIENT_ID>            > <P>.rg.tsv
#
# Columns: sample proj n_rg rg_id lb pu dt cn pl ds pm aligner_pg  (multi-RG BAMs: fields
# joined by ',' over the RGs). Use it to tell sequencing/library batches apart, e.g. PD45886's
# 15 colonies with a ~5x higher passing-clip rate (2026-10-09).
# =============================================================================
set -uo pipefail
P="${1:?usage: rg_headers.sh <PATIENT_ID>}"
# shellcheck source=cluster/tprt/common.sh
. "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/tprt/common.sh"
load_irods    || { echo "rg_headers: iquest not available (module load IRODS)" >&2; exit 1; }
load_samtools || { echo "rg_headers: samtools not available (module load $SAMTOOLS_MODULE)" >&2; exit 1; }

printf 'sample\tproj\tn_rg\trg_id\tlb\tpu\tdt\tcn\tpl\tds\tpm\taligner_pg\n'
# every colonies.tsv row, whatever its assembly (fleet.sh samples filters to one assembly)
d="$(ls -d "$PATIENTS_DIR"/*/"$P" 2>/dev/null | head -1)"
[ -s "$d/colonies.tsv" ] || { echo "rg_headers: no $PATIENTS_DIR/*/$P/colonies.tsv" >&2; exit 1; }
awk -F'\t' '!/^#/ && $1 != "donor" && NF >= 7 { print $7 "\t" $2 }' "$d/colonies.tsv" | sort -u |
while IFS=$'\t' read -r s proj; do
    bam="$(irods_readable_bam "$s" "$proj")" \
        || { printf '%s\t%s\tNO_READABLE_REPLICA\n' "$s" "$proj"; continue; }
    samtools view -H "$bam" | awk -F'\t' -v s="$s" -v proj="$proj" '
        function add(k, v) { f[k] = (k in f) ? f[k] "," v : v }
        $1 == "@RG" { n++; for (i = 2; i <= NF; i++) { k = substr($i, 1, 2); if (k ~ /^(ID|LB|PU|DT|CN|PL|DS|PM)$/) add(k, substr($i, 4)) } }
        $1 == "@PG" && pg == "" { for (i = 2; i <= NF; i++) if ($i ~ /^(PN|VN):/) pg = pg (pg ? " " : "") substr($i, 4) }
        END { printf "%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", s, proj, n,
                     f["ID"], f["LB"], f["PU"], f["DT"], f["CN"], f["PL"], f["DS"], f["PM"], pg }'
done
