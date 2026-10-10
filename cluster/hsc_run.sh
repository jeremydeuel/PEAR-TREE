#!/bin/bash
# =============================================================================
# cluster/hsc_run.sh — the deployed TPRT pipeline, one patient at a time (farm22).
#
# The configuration validated on PD37590 (arm C + genotype2 with reference-bias correction),
# from ONE checkout of branch hsc-deploy:
#   discovery   config.discovery.grch38.tprt2frag  (>= 2 fragments per junction end IN EACH colony)
#               (GRCh37 patients: config.discovery.grch37.tprt2frag, same keys, hs37d5 contig names)
#   combine     Rust peartree-combine, src/config.py = config.py.grch38.tprt (+ tprt/arm_config.py)
#   genotype    peartree-genotype2, config.genotype2.grch38.refbias, two-sided loci only
#   joint       --ref-bias auto (global / per kind / per colony), colony zygosity, no NOISE cap
#   annotate    annotate_v2 + RTE library (Rust peartree-rte, tools/rte port; spill in $HSC_ROOT/tmp/<P>), locus_class
#   report      tools/phylo/tree_fit.py + cluster/somatic_table.py -> <P>.somatic.xlsx
#
# Assembly: per patient from colonies.tsv (all WGS rows GRCh38, or all hs37d5/GRCh37), override
# HSC_ASSEMBLY=GRCh38|GRCh37. GRCh37 runs NATIVELY (no remap): discovery allowlist 1..Y,
# PT_ASSEMBLY=GRCh37 makes src/config.py use hg19.2bit + the hs1->hg19 chain (combine, annotate),
# genotype2 reads hg19.2bit (chr prefix toggled), and the per-BAM header gate checks hs37d5.
#
# Usage (head node, from anywhere):
#   bash <checkout>/cluster/hsc_run.sh setup              # build binaries, venv link, src/config.py
#   bash <checkout>/cluster/hsc_run.sh populate <P>       # fill patients/*/<P>/colonies.tsv (nst_links headers; iRODS fallback)
#        REPOPULATE=1 redoes a populated one; PREFER_ASSEMBLY (nst_links) / PREFER_PROJECTS (iRODS) pick
#        among a colony's several releases
#   bash <checkout>/cluster/hsc_run.sh submit <P>         # samples.tsv + the whole pipeline
#   bash <checkout>/cluster/hsc_run.sh status <P>
#
# Layout:  $HSC_ROOT/<P>/<P>/      pipeline run dir (discovery/ insertions/ genotypes/ fit/ logs/ ...)
#          $HSC_ROOT/tmp/<P>/      annotate scratch
#          $HOME/results/hsc/<P>/  NFS copies: calls, joint table, annotation, <P>.somatic.xlsx
# Env overrides: HSC_ROOT, VENV, HSC_ASSEMBLY, KEEP_BAMS=1 (skip the staged-BAM cleanup),
# PT_RESTAGE=1 (drop each BAM after discovery, stage it again for genotyping: large patients), plus every
# pipeline.sh knob (GT_THROTTLE, CI_MEM, ...).
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PT_ROOT="$(cd "$HERE/.." && pwd)"
JD="${JD:-/lustre/scratch126/casm/teams/team273/users/jd43}"
TPRT_ROOT="${TPRT_ROOT:-$JD/tprt_ab}"                 # shared resources (hs1.2bit, gene model, minimap2)
TPRT_RES="${TPRT_RES:-$TPRT_ROOT/resources}"
HSC_ROOT="${HSC_ROOT:-$JD/hsc}"
VENV="${VENV:-$TPRT_ROOT/PEAR-TREE/venv}"            # the tprt kit venv (edlib, mappy, pandas, scipy)
RESULTS_ROOT="${RESULTS_ROOT:-$HOME/results/hsc}"
RUST_MODULE="${RUST_MODULE:-rust/1.87.0}"

die()  { echo "hsc_run: $*" >&2; exit 1; }
note() { echo "[$(date +%H:%M:%S)] $*" >&2; }

patient_dir() {
    local d; d=$(ls -d "$PT_ROOT"/patients/*/"$1" 2>/dev/null || true)
    if [ -z "$d" ]; then
        # a PD id that is one half of a shared-tree folder (transplant Pair<N>_<donor>_<recipient>,
        # SDS5_PD42190_PD45888): that folder is the unit -- never analyse one person of a pair alone
        d=$(ls -d "$PT_ROOT"/patients/*/*_"$1" "$PT_ROOT"/patients/*/*_"$1"_* 2>/dev/null || true)
        [ "$(printf '%s\n' "$d" | grep -c .)" -eq 1 ] \
            && die "$1 shares one tree with the other PD id(s) in $(basename "$d"); run that unit: hsc_run.sh ${CMD:-submit} $(basename "$d")"
    fi
    [ "$(printf '%s\n' "$d" | grep -c .)" -eq 1 ] || die "expected exactly one patients/*/$1 in $PT_ROOT, found: ${d:-none}"
    echo "$d"
}

# GRCh38 | GRCh37 for patient dir $1: HSC_ASSEMBLY, else the colonies.tsv assembly of its WGS rows
patient_assembly() {
    if [ -n "${HSC_ASSEMBLY:-}" ]; then echo "$HSC_ASSEMBLY"; return; fi
    local a
    a="$(awk -F'\t' 'NR>2 && $3 ~ /WGS/ && $7!="" {
            if ($6=="GRCh38") g38++; else if ($6=="hs37d5_GRCh37" || $6=="GRCh37(hs37d5/hg19)") g37++; else oth++ }
        END { if (g38 && !g37) print "GRCh38"; else if (g37 && !g38) print "GRCh37";
              else printf "mixed/unknown (GRCh38=%d GRCh37=%d other=%d)\n", g38, g37, oth }' "$1/colonies.tsv")"
    case "$a" in GRCh38|GRCh37) echo "$a" ;; *) die "$1/colonies.tsv: assembly $a -- set HSC_ASSEMBLY=GRCh38 or GRCh37" ;; esac
}

cmd_setup() {
    [ -x "$VENV/bin/python" ] || die "no venv at $VENV (set VENV)"
    [ -e "$PT_ROOT/venv" ] || ln -s "$VENV" "$PT_ROOT/venv"   # install.sh check-config defaults to $PT_ROOT/venv
    note "building Rust binaries in $PT_ROOT"
    module load "$RUST_MODULE" >/dev/null 2>&1 || true
    bash "$PT_ROOT/cluster/build.sh"
    for b in peartree-discovery peartree-combine peartree-genotype2 peartree-rte; do
        [ -x "$PT_ROOT/rust/$b/target/release/$b" ] || die "build did not produce $b"
    done
    cat > "$PT_ROOT/src/config.py" <<EOF
# GENERATED by cluster/hsc_run.sh setup — the TPRT configuration (arm B/C python side:
# cluster/config.py.grch38.tprt + cluster/tprt/arm_config.py). Do not edit.
import os
import sys
os.environ.setdefault('TPRT_ROOT', '$TPRT_ROOT')
os.environ.setdefault('TPRT_RES', '$TPRT_RES')
sys.path.insert(0, '$PT_ROOT/cluster/tprt')
from arm_config import build  # noqa: E402
TPRT_AB_ARM = 'B'
CONFIG = build(TPRT_AB_ARM)
EOF
    (cd "$PT_ROOT" && "$VENV/bin/python" -c "import sys; sys.path.insert(0, 'src'); import config; \
sys.path.insert(0, 'cluster/tprt'); import arm_config; print(arm_config.describe(config.TPRT_AB_ARM))")
    note "setup done at $(git -C "$PT_ROOT" log --oneline -1)"
}

cmd_populate() {
    local P="${1:?usage: hsc_run.sh populate <PATIENT_ID>}" d org
    d="$(patient_dir "$P")"; org="$(basename "$(dirname "$d")")"
    local nst="${NST:-/nfs/cancer_ref01/nst_links/live}"
    # REPOPULATE=1: redo an already-populated colonies.tsv (kept as colonies.tsv.bak.<date>)
    if [ "${REPOPULATE:-0}" = 1 ] && awk '!/^#/ && !/^donor\t/ && NF {f=1; exit} END {exit !f}' "$d/colonies.tsv"; then
        cp "$d/colonies.tsv" "$d/colonies.tsv.bak.$(date +%Y%m%d%H%M%S)"
        printf 'donor\tproj\tds\treadlen\tmapped\tassembly\tsample\n# PENDING -- REPOPULATE=1\n' > "$d/colonies.tsv"
        note "REPOPULATE=1: old colonies.tsv kept as $(ls -t "$d"/colonies.tsv.bak.* | head -1)"
    fi
    if [ -d "$nst" ]; then
        # header reads on nst_links, run HERE (compute nodes may not mount it: PD51635 node-13-14)
        cd "$PT_ROOT"
        ONLY="$org" PAR="${PAR:-4}" PREFER_ASSEMBLY="${PREFER_ASSEMBLY:-}" bash cluster/populate_colonies_tsv.sh
    else
        note "nst_links not visible on $(hostname) ($nst): tip -> project from iRODS (iquest)"
        populate_irods "$P" "$d" "${HSC_ASSEMBLY:-GRCh38}"
    fi
    grep -v '^#' "$d/colonies.tsv" | awk -F'\t' 'NR>1 {n++; a[$6]++} END {printf "%d BAMs:", n; for (k in a) printf " %s=%d", k, a[k]; print ""}'
}

# iRODS fallback (farm22-head2 lost /nfs/cancer_ref01, 2026-10-06): one iquest for the donor's
# *.sample.dupmarked.bam objects, matched to the tree tips. Assay and assembly are NOT read here
# (no header access without nst_links): rows say ds=WGS_unverified / assembly=GRCh38 (or
# hs37d5_GRCh37 with HSC_ASSEMBLY=GRCh37), and `submit` always runs pipeline.sh's per-BAM header
# gate (PT_HEADER_GATE=<assembly>) after staging, which also catches a wrong guess here.
# A tip found in more than one project (WGS + targeted twins, Chapman 2024) is left out, unless
# PREFER_PROJECTS="p1 p2 ..." names one of them (first listed wins): PD45534's GRCh38 releases are
# 2450/3819/3838, its hs37d5 + targeted ones 2515/2566/2902/3191.
populate_irods() {
    local P="$1" d="$2" label tips hits
    case "$3" in GRCh38) label=GRCh38 ;; GRCh37) label=hs37d5_GRCh37 ;; *) die "HSC_ASSEMBLY must be GRCh38 or GRCh37, got '$3'" ;; esac
    module load IRODS >/dev/null 2>&1 || true
    command -v iquest >/dev/null || die "iquest not on PATH (module load IRODS)"
    tips="$(mktemp)"; hits="$(mktemp)"
    grep -oE '[(,][A-Za-z][A-Za-z0-9._-]*' "$d"/*.tree | cut -c2- | sort -u > "$tips"
    # one iquest per PD prefix of the tree tips, not of the patient id: a patient dir can hold one
    # individual sampled under two PD ids (SDS5_PD42190_PD45888: PD42190* + PD45888* tips)
    local pfx pfxs
    pfxs="$(grep -oE '^PD[0-9]+' "$tips" | sort -u)"
    [ -n "$pfxs" ] || pfxs="$P"
    for pfx in $pfxs; do
        iquest --no-page "%s/%s" "SELECT COLL_NAME, DATA_NAME WHERE DATA_NAME like '${pfx}%.sample.dupmarked.bam'"
    done | awk -F/ '$2=="cgp" && $3=="intproj" && $5=="sample" {print $6 "\t" $4}' | sort -u > "$hits"
    [ -s "$hits" ] || { rm -f "$tips" "$hits"; die "iquest found no $(echo $pfxs | tr ' ' '/')*.sample.dupmarked.bam in iRODS"; }
    local tmp="$d/colonies.tsv.tmp.$$"
    {
        printf 'donor\tproj\tds\treadlen\tmapped\tassembly\tsample\n'
        printf '# populated %s by cluster/hsc_run.sh populate from iRODS (iquest; nst_links not mounted). ds/assembly NOT read from headers: verified per BAM after staging (PT_HEADER_GATE=%s).\n' "$(date +%F)" "$3"
        awk -F'\t' -v donor="$P" -v asm="$label" -v pref="${PREFER_PROJECTS:-}" '
            BEGIN {np = split(pref, pp, " "); for (i = 1; i <= np; i++) rank[pp[i]] = i}
            NR==FNR {tip[$1]=1; next}
            ($1 in tip) {n[$1]++; proj[$1]=$2
                         if (($2 in rank) && (!($1 in best) || rank[$2] < rank[best[$1]])) best[$1]=$2}
            END {for (s in n) {
                     if (s in best) print donor "\t" best[s] "\tWGS_unverified\tNA\tNA\t" asm "\t" s
                     else if (n[s]==1 && np == 0) print donor "\t" proj[s] "\tWGS_unverified\tNA\tNA\t" asm "\t" s}}' "$tips" "$hits" | sort -t$'\t' -k7,7
    } > "$tmp" && mv -f "$tmp" "$d/colonies.tsv"
    awk -F'\t' 'NR==FNR {tip[$1]=1; next} ($1 in tip) {n[$1]++; p[$1]=p[$1] " " $2}
        END {for (s in n) if (n[s]>1) print "  in several projects:" , s, p[s]}' "$tips" "$hits" >&2
    [ -z "${PREFER_PROJECTS:-}" ] || awk -F'\t' -v pref="$PREFER_PROJECTS" '
        BEGIN {np = split(pref, pp, " "); for (i = 1; i <= np; i++) ok[pp[i]] = 1}
        NR==FNR {tip[$1]=1; next} ($1 in tip) {seen[$1]=1; if ($2 in ok) good[$1]=1; p[$1]=p[$1] " " $2}
        END {for (s in seen) if (!(s in good)) print "  no release in PREFER_PROJECTS, left out:", s, p[s]}' "$tips" "$hits" >&2
    awk -F'\t' 'NR==FNR {seen[$1]=1; next} !($1 in seen) {print "  tree tip with no BAM in iRODS:", $1}' "$hits" "$tips" >&2
    note "projects: $(awk -F'\t' 'NR==FNR {tip[$1]=1; next} ($1 in tip) {c[$2]++} END {for (k in c) printf "%s=%d ", k, c[k]}' "$tips" "$hits")"
    rm -f "$tips" "$hits"
}

cmd_submit() {
    local P="${1:?usage: hsc_run.sh submit <PATIENT_ID>}" d
    d="$(patient_dir "$P")"
    ls "$d"/*.tree >/dev/null 2>&1 || die "no tree in $d"
    # one awk, no pipe: `grep -q` quitting early SIGPIPEd the writer, and under pipefail a 722-row
    # colonies.tsv (PD49229) read as "no rows"
    awk '!/^#/ && NF && ++n > 1 { found = 1; exit } END { exit !found }' "$d/colonies.tsv" \
        || die "$d/colonies.tsv has no rows — run: bash $PT_ROOT/cluster/hsc_run.sh populate $P"
    [ -s "$PT_ROOT/src/config.py" ] && grep -q TPRT_AB_ARM "$PT_ROOT/src/config.py" \
        || die "no generated src/config.py — run: bash $PT_ROOT/cluster/hsc_run.sh setup"
    for b in peartree-discovery peartree-combine peartree-genotype2 peartree-rte; do
        [ -x "$PT_ROOT/rust/$b/target/release/$b" ] || die "missing $b — run setup"
    done

    local asm disc twobit chain exons
    asm="$(patient_assembly "$d")"
    case "$asm" in
        GRCh38) disc="$PT_ROOT/cluster/config.discovery.grch38.tprt2frag"; twobit="$JD/hg38.2bit"; chain="$JD/hs1.hg38.all.chain.gz" ;;
        GRCh37) disc="$PT_ROOT/cluster/config.discovery.grch37.tprt2frag"; twobit="$JD/hg19.2bit"; chain="$JD/hs1.hg19.all.chain.gz" ;;
        *) die "HSC_ASSEMBLY must be GRCh38 or GRCh37, got '$asm'" ;;
    esac
    # the assembly-specific files, checked now rather than hours in (discovery exits on a missing
    # exon track; a missing 2bit/chain would only surface at combine / genotype)
    exons="$(sed -n 's/^exon_annotation *= *//p' "$disc" | tail -1)"
    for f in "$disc" "$twobit" "$chain" ${exons:+"$exons"}; do
        [ -s "$f" ] || die "$asm run needs $f (missing or empty)"
    done

    local W="$HSC_ROOT/$P"
    mkdir -p "$W" "$HSC_ROOT/tmp/$P"
    PATIENTS_DIR="$PT_ROOT/patients" SAMPLES_ASSEMBLY="$asm" bash "$PT_ROOT/cluster/fleet.sh" samples "$P" > "$W/samples.tsv"
    local n; n=$(grep -c . "$W/samples.tsv" || true)
    [ "$n" -gt 0 ] || die "fleet.sh samples $P listed no $asm WGS colonies"
    note "$P: $n $asm WGS samples -> $W/samples.tsv"
    # capacity: without PT_RESTAGE every BAM stays staged until genotyping ends (~25 GB per colony
    # at ~15-30x). Refuse a run that cannot fit the team's free Lustre quota (PD49229: 722 colonies
    # ~17 TB vs 8.2 TB free, 2026-10-07); PT_RESTAGE=1 bounds the peak to the running jobs.
    if [ "${PT_RESTAGE:-0}" != 1 ]; then
        local need_tb free_tb
        need_tb=$(awk -v n="$n" -v g="${HSC_GB_PER_BAM:-25}" 'BEGIN { printf "%.1f", n * g / 1000 }')
        free_tb=$(lfs quota -g team273 /lustre/scratch126 2>/dev/null \
            | awk 'NR==3 && NF>=4 { printf "%.1f", ($4 - $2) / 1e9 } NR==4 && NF>=3 { printf "%.1f", ($3 - $1) / 1e9 }')
        if [ -n "$free_tb" ] && awk -v a="$need_tb" -v b="$free_tb" 'BEGIN { exit !(a > b * 0.9) }'; then
            die "$P: ~$need_tb TB staged at peak but team273 has $free_tb TB free -- rerun with PT_RESTAGE=1 (re-stage per genotype job)"
        fi
        note "capacity: ~$need_tb TB staged at peak, ${free_tb:-?} TB free"
    else
        note "PT_RESTAGE=1: each BAM is dropped after discovery and staged again for genotyping"
    fi
    (export PT_ASSEMBLY="$asm"; cd "$PT_ROOT" && "$VENV/bin/python" -c "import sys; sys.path.insert(0, 'cluster/tprt'); \
import arm_config; print(arm_config.describe('B'))") >&2

    env WORKROOT="$W" RESULTS_DIR="$RESULTS_ROOT" VENV="$VENV" \
        PT_ASSEMBLY="$asm" GENOME_2BIT="$twobit" DISC_CFG="$disc" \
        GENO_CFG="$PT_ROOT/cluster/config.genotype.grch38.tprt" GENO_ONE_SIDED=0 \
        GENOTYPE_IMPL=v2 GENO2_CFG="$PT_ROOT/cluster/config.genotype2.grch38.refbias" \
        JOINT_ARGS="--ref-bias auto" COMBINE_IMPL=rust \
        TPRT_ANNOT_TMP="$HSC_ROOT/tmp/$P" \
        GT_THROTTLE="${GT_THROTTLE:-20}" GT_MEM_T1="${GT_MEM_T1:-4000}" GT_MEM_T2="${GT_MEM_T2:-16000}" \
        SD_MEM_T1="${SD_MEM_T1:-4000}" CI_MEM="${CI_MEM:-16000}" CI_CORES="${CI_CORES:-8}" \
        AN_MEM="${AN_MEM:-32000}" AN_CORES="${AN_CORES:-8}" \
        PT_HEADER_GATE="$asm" \
        PT_JOB_PREFIX="${P}_hsc" PT_NO_CLEANUP="${KEEP_BAMS:-0}" PT_JOBIDS_FILE="$W/jobids.tsv" \
        bash "$PT_ROOT/cluster/pipeline.sh" submit-list "$P" "$W/samples.tsv"
    note "results will land in $RESULTS_ROOT/$P (final: $P.somatic.xlsx)"
}

cmd_status() {
    local P="${1:?usage: hsc_run.sh status <PATIENT_ID>}"
    PT_JOB_PREFIX="${P}_hsc" WORKROOT="$HSC_ROOT/$P" bash "$PT_ROOT/cluster/pipeline.sh" status "$P"
}

CMD="${1:-}"
CMD="${1:-}"
case "${1:-}" in
    setup)    shift; cmd_setup "$@" ;;
    populate) shift; cmd_populate "$@" ;;
    submit)   shift; cmd_submit "$@" ;;
    status)   shift; cmd_status "$@" ;;
    *) sed -n '2,30p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
