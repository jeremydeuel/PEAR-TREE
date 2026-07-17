#!/bin/bash
# Where are a donor's BAMs? Search every place they could be, in cost order.
#
#   bash cluster/locate_donors.sh PD51632 PD51633 PD41276 ...
#   bash cluster/locate_donors.sh --missing        # every donor in donors.txt with no BAM yet
#
# WHY THIS EXISTS. The catalogue only reads nst_links, and **nst_links is an INCOMPLETE view
# of iRODS** — measured 2026-07: PD57333 has 156 BAMs in iRODS but only 79 linked (49%
# missing); PD51634 1,178 vs 982. So "not in nst_links" NEVER means "not on the farm", and a
# tool that reads only nst_links cannot see its own blind spot. Data can also sit in staging
# or a team scratch dir.
#
# RESULT of the first run (2026-07): every CML donor (PD51632-35, PD56961, PD57332-35,
# PD60243) and every Fabre donor (PD41276, PD34493, PD41305) IS on the farm, in both iRODS
# and nst_links. None needed an EGA download. CML lives in projects 3044/3045 (+3224/3168
# for PD60243); Fabre in 2001/2011/2015/2305/2395/3156.
#
# CAUTION on project spread: PD51634 has 1,178 BAMs across SEVEN projects against a
# published manifest of 151 samples (~8x re-release duplication). Select from the paper's
# manifest, never by globbing a donor across projects.
#
# Search order (cheapest first):
#   1. nst_links   NFS, one glob per project      -- free
#   2. staging     our own staged pulls           -- free
#   3. team lustre a find, bounded depth          -- slowish
#   4. iRODS       the ARCHIVE OF RECORD (nst_links is an incomplete view of it)
#                  -- needs `module load IRODS` (CAPITALISED); SLOW
#
# Env: NST, STAGING, TEAMROOT, DO_IRODS=1 to include the iRODS step
set -uo pipefail

NST=${NST:-/nfs/cancer_ref01/nst_links/live}
STAGING=${STAGING:-/lustre/scratch126/casm/staging/team273/jd43}
TEAMROOT=${TEAMROOT:-/lustre/scratch126/casm/teams/team273/users/jd43}
DO_IRODS=${DO_IRODS:-0}

if [ "${1:-}" = "--missing" ]; then
    CAT=${CAT:-$HOME/catalogue/catalogue.tsv}
    [ -s "$CAT" ] || { echo "no catalogue at $CAT — pass donors explicitly" >&2; exit 1; }
    mapfile -t DONORS < <(sed 's/#.*//' cluster/donors.txt | awk 'NF{print $1}' | while read -r d; do
        awk -F'\t' -v d="$d" 'NR>1 && $3==d && $7!="NO_BAM"{f=1} END{exit f?0:1}' "$CAT" || echo "$d"; done)
    echo "donors in donors.txt with no BAM in the catalogue: ${#DONORS[@]}"
else
    [ $# -gt 0 ] || { echo "usage: $0 <donor> [...]  |  $0 --missing" >&2; exit 1; }
    DONORS=("$@")
fi

for d in "${DONORS[@]}"; do
    echo "=============== $d"
    hit=0

    # Donor ids are NOT prefix-free: PD5163 is a real donor AND a prefix of PD51632/3/4/5.
    # Sample names are <DONOR><lowercase letter>_lo####, so drop any hit where a DIGIT follows
    # the donor id -- that is a different, longer donor. (Same rule as catalogue_scan.sh.)
    not_prefix_of_other() { awk -F/ -v d="$d" '{s=$NF; if (s !~ "^"d"[0-9]") print}'; }

    n=$(ls -d "$NST"/*/"$d"* 2>/dev/null | not_prefix_of_other | wc -l | tr -d ' ')
    if [ "$n" -gt 0 ]; then
        hit=1
        echo "  nst_links: $n sample dirs"
        ls -d "$NST"/*/"$d"* 2>/dev/null | not_prefix_of_other | awk -F/ '{print $(NF-1)}' | sort | uniq -c \
            | awk '{printf "    proj %-6s %s samples\n", $2, $1}'
    else
        echo "  nst_links: none"
    fi

    if [ -d "$STAGING" ]; then
        s=$(ls -d "$STAGING"/*/"$d"* 2>/dev/null | wc -l | tr -d ' ')
        [ "$s" -gt 0 ] && { hit=1; echo "  staging: $s dirs"; ls -d "$STAGING"/*/"$d"* 2>/dev/null | head -3 | sed 's/^/    /'; } \
                       || echo "  staging: none"
    fi

    if [ -d "$TEAMROOT" ]; then
        t=$(find "$TEAMROOT" -maxdepth 4 -name "$d*" -print -quit 2>/dev/null)
        [ -n "$t" ] && { hit=1; echo "  team lustre: $t"; } || echo "  team lustre: none"
    fi

    if [ "$DO_IRODS" = 1 ]; then
        if command -v iquest >/dev/null 2>&1; then
            # Query by DATA_NAME, never by attribute name: `imeta qu -d sample = X` returns
            # "No rows found" in the cgp zone (NPG attribute names do not apply here).
            # iRODS is the ARCHIVE OF RECORD and nst_links is an incomplete view of it —
            # measured: PD57333 has 156 BAMs in iRODS but only 79 linked. So a donor absent
            # from nst_links may still be here.
            # iRODS `like` cannot express "not followed by a digit", so filter after the fact
            # on the collection's sample component -- same prefix rule as above.
            iq=$(iquest --no-page "SELECT COLL_NAME WHERE COLL_NAME like '/cgp/intproj/%/sample/${d}%' AND DATA_NAME like '%.sample.dupmarked.bam'" 2>/dev/null \
                 | grep '^COLL_NAME' | sed 's/^COLL_NAME = //' | awk -F/ -v d="$d" '$NF !~ "^"d"[0-9]"')
            ib=$(printf '%s\n' "$iq" | grep -c '^/cgp/')
            ip=$(printf '%s\n' "$iq" | grep '^/cgp/' | awk -F/ '{print $4}' | sort -u | paste -sd, -)
            [ "${ib:-0}" -gt 0 ] && hit=1
            echo "  iRODS/cgp: $ib BAMs   projects: ${ip:-none}"
            [ "${ib:-0}" -gt "$n" ] && echo "    NB iRODS has MORE BAMs than nst_links ($ib vs $n) — nst_links is incomplete here"
        else
            echo "  iRODS: iquest not on PATH — run: module load IRODS   (CAPITALISED)"
        fi
    fi

    [ "$hit" -eq 1 ] || echo "  ==> NOT FOUND on any local filesystem"
done

echo
echo "If a donor is nowhere: it is EGA-only and must be requested/downloaded, OR it lives"
echo "under a sample naming scheme we did not guess. Check the paper's manifest for the"
echo "actual sample names before concluding it is absent."
