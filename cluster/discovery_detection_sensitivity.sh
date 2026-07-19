#!/usr/bin/env bash
# ============================================================================
# Per-tip DISCOVERY detection sensitivity of young-MEI phylogeny-TP loci,
# compared across the ablation-sweep discovery arms.
#
# Question answered: of the tips we EXPECT a given MEI in (confident het/hom
# carriers in the baseline genotype table), in how many did that arm's per-tip
# discovery actually emit the breakpoint? That ratio, aggregated per arm, is the
# discovery detection sensitivity — the ONLY thing that varies across sweep arms
# (genotyping config is identical; only discovery config differs).
#
# Inputs it needs on the farm:
#   TARGETS  the young_mei_phyloTP_targets.tsv produced locally (scp it up).
#            cols: locus  ann_class  ann_verdict  n_exp_carrier_tips  expected_carrier_tips(csv)
#   REPO     the sweep run dir; arm discovery outputs are $REPO/discovery_grch38*<arm>/<tip>.txt.gz
#            (same FASTQ-style @chr:start-end:{L,R} headers combine_mei.sh consumes).
#
# Self-contained: absolute defaults, no placeholders. Run on a farm22 head node.
#   TARGETS=/path/to/young_mei_phyloTP_targets.tsv bash cluster/discovery_detection_sensitivity.sh
# ============================================================================
set -euo pipefail

REPO="${REPO:-/lustre/scratch126/casm/teams/team273/users/jd43/PEAR-TREE-pd44579}"
TARGETS="${TARGETS:-$REPO/young_mei_phyloTP_targets.tsv}"
WINDOW="${WINDOW:-50}"          # bp: breakpoint-midpoint match tolerance across arms
OUT="${OUT:-$REPO/discovery_detection_sensitivity.tsv}"

[ -s "$TARGETS" ] || { echo "no TARGETS file: $TARGETS (scp it up)" >&2; exit 1; }
[ -d "$REPO" ]    || { echo "no REPO dir: $REPO" >&2; exit 1; }

# every arm discovery dir: discovery_grch38, discovery_grch38_alloff, _loo_*, _aoi_*, _noSPEC8 ...
mapfile -t ARMDIRS < <(find "$REPO" -maxdepth 1 -type d -name 'discovery_grch38*' | sort)
[ "${#ARMDIRS[@]}" -gt 0 ] || { echo "no discovery_grch38* dirs under $REPO" >&2; exit 1; }
echo "arms found (${#ARMDIRS[@]}):" >&2
for d in "${ARMDIRS[@]}"; do echo "  $(basename "$d")  ($(ls "$d"/*.txt.gz 2>/dev/null | wc -l | tr -d ' ') tip files)" >&2; done

export TARGETS WINDOW OUT
python3 - "${ARMDIRS[@]}" <<'PY'
import os, sys, gzip, re
from collections import defaultdict

TARGETS=os.environ['TARGETS']; W=int(os.environ['WINDOW']); OUT=os.environ['OUT']
armdirs=[a for a in sys.argv[1:] if a.strip()]

# --- load targets: locus, class, expected carrier tips ---
loc_re=re.compile(r'^([^:]+):(\d+)-(\d+)$')
targets=[]   # (chrom, mid, cls, [tips])
tips_needed=set()
with open(TARGETS) as f:
    hdr=f.readline()
    for line in f:
        p=line.rstrip('\n').split('\t')
        if len(p)<5: continue
        m=loc_re.match(p[0])
        if not m: continue
        chrom,s,e=m.group(1),int(m.group(2)),int(m.group(3))
        tips=[t for t in p[4].split(',') if t]
        targets.append((chrom,(s+e)//2,p[1],tips))
        tips_needed.update(tips)

def load_disc(path):
    """chrom -> sorted list of breakpoint midpoints for one tip's discovery file."""
    d=defaultdict(list)
    try:
        with gzip.open(path,'rt') as fh:
            for line in fh:
                if line and line[0]=='@':
                    h=line[1:].rstrip('\n')
                    if h.endswith(':L') or h.endswith(':R'): h=h[:-2]
                    m=loc_re.match(h)
                    if m: d[m.group(1)].append((int(m.group(2))+int(m.group(3)))//2)
    except FileNotFoundError:
        return None
    for c in d: d[c].sort()
    return d

import bisect
def hit(mids, mid):
    if not mids: return False
    i=bisect.bisect_left(mids, mid)
    for j in (i-1,i):
        if 0<=j<len(mids) and abs(mids[j]-mid)<=W: return True
    return False

rows=[]
for ad in armdirs:
    arm=os.path.basename(ad)
    cache={}
    def disc_for(tip):
        if tip not in cache: cache[tip]=load_disc(os.path.join(ad,tip+'.txt.gz'))
        return cache[tip]
    # aggregate overall + per class
    agg=defaultdict(lambda:[0,0,0,0])  # cls -> [exp_tips, det_tips, n_loci, loci_fully_detected]
    for chrom,mid,cls,tips in targets:
        exp=0; det=0
        for t in tips:
            d=disc_for(t)
            if d is None: continue          # missing tip file for this arm -> not counted in denom
            exp+=1
            if hit(d.get(chrom,[]), mid): det+=1
        for key in (cls,'ALL'):
            a=agg[key]; a[0]+=exp; a[1]+=det; a[2]+=1; a[3]+= (1 if exp>0 and det==exp else 0)
    for cls in ('ALL','ALU','L1','SVA'):
        if cls in agg:
            exp,det,nl,full=agg[cls]
            sens=det/exp if exp else 0.0
            rows.append((arm,cls,nl,exp,det,f"{sens:.4f}",full))

hdr=["arm","class","n_loci","exp_carrier_tips","detected_tips","tip_sensitivity","loci_all_tips_detected"]
with open(OUT,'w') as o:
    o.write("\t".join(hdr)+"\n")
    for r in rows: o.write("\t".join(map(str,r))+"\n")

wcol=max(len(r[0]) for r in rows) if rows else 4
print("\t".join(hdr))
for r in rows:
    if r[1]=='ALL':
        print(f"{r[0]:<{wcol}}  {r[1]:<4} loci={r[2]:<4} exp={r[3]:<5} det={r[4]:<5} sens={r[5]}  full-detect-loci={r[6]}")
print(f"\nwrote {OUT}", file=sys.stderr)
PY
echo "=== done -> $OUT ===" >&2
