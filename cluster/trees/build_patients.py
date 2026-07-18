#!/usr/bin/env python3
# Build the patients/ tree: patients/<organ_group>/<patient>/{<tree>.tree, colonies.tsv}
# plus a per-organ <organ>.md summary and a top-level PLAN.md.
#
# Driven by donor_table.tsv (organ authority) + donor_trees.tsv (canonical tips/alias).
# Authoritative tree file per patient = the one whose leaf count matches donor_trees tips.
#   - Chapman HSCT: twin donor+recipient share ONE patient folder (one pair tree).
#   - Liver explants PD48367/PD48372: all 8 per-region trees copied into the one patient.
#   - 13 donors have no tree on disk yet -> folder + placeholder + NO_TREE_ON_DISK marker.
# The per-patient TSV is a PLACEHOLDER now (header only); it will be filled from the per-tip
# header probe currently running on the farm. See PLAN.md.
import re,os,glob,json,shutil,collections

ROOT=os.path.dirname(os.path.abspath(__file__))
OUT=os.path.join(os.path.dirname(os.path.dirname(ROOT)),'patients')  # repo_root/patients
TREE_DIRS=['coorens2025stomach','mitchell2025','mitchell2022','chapman2021',
  'chapman2024hsct/trees_main','chapman2025/EM','chapman2025/KY','chapman2025/MF',
  'chapman2025/MSC_BMT','chapman2025/MSC_fetal','chapman2025/NW','chapman2025/SN',
  'robinson2022mutyh','bcl_unpublished','leesix2019colorectal','williams2022mpn']

def leaves(fp):
    s=open(fp).read()
    return len(re.findall(r'[(,]([A-Za-z0-9][A-Za-z0-9._\-]*)\s*:', s))
idx=[]
for r in TREE_DIRS:
    for fp in glob.glob(os.path.join(ROOT,r,'*.tree')):
        idx.append((os.path.basename(fp), r, leaves(fp)))

ORGAN_LABEL={'bronchial_epithelium':'Bronchial epithelium (lung)','tonsil':'Tonsil',
 'cord_blood':'Cord blood','foetal_liver':'Foetal liver','liver':'Liver',
 'bone_marrow':'Bone marrow','blood_bm_spleen':'Blood + BM + spleen','colorectum':'Colorectum / ileum',
 'colon':'Colon','gastric_epithelium':'Gastric epithelium (stomach)','peripheral_blood':'Peripheral blood / BM',
 'unknown_tissue_bcl':'Unknown tissue (BCL, unplaced)'}
def organ(tissue,unit):
    t=tissue.lower()
    if tissue=='?': return 'unknown_tissue_bcl'
    if 'bronchial' in t: return 'bronchial_epithelium'
    if 'tonsil' in t: return 'tonsil'
    if 'cord blood' in t: return 'cord_blood'
    if 'foetal liver' in t: return 'foetal_liver'
    if 'liver' in t: return 'liver'
    if t.startswith('bone marrow'): return 'bone_marrow'
    if 'spleen' in t: return 'blood_bm_spleen'
    if 'colorectal' in t or 'ileum' in t: return 'colorectum'
    if 'colon' in t: return 'colon'
    if 'gastric' in t or 'stomach' in t: return 'gastric_epithelium'
    if 'peripheral blood' in t or 'blood/bone' in t: return 'peripheral_blood'
    return 'unknown_tissue_bcl'

# donor_table
rows={}
with open(os.path.join(ROOT,'donor_table.tsv')) as f:
    f.readline()
    for line in f:
        c=line.rstrip('\n').split('\t'); d=c[0]
        if d in rows: continue
        m=re.match(r'\d+',c[11]); tips=int(m.group()) if m else 0
        rows[d]=dict(tissue=c[2],unit=c[3],tips=tips,alias=c[1],paper=c[7],diag=c[6])
# donor_trees: canonical tips + alias (also pair totals)
dt_tips={}; dt_alias={}
with open(os.path.join(ROOT,'donor_trees.tsv')) as f:
    f.readline()
    for line in f:
        c=line.rstrip('\n').split('\t')
        dt_tips[c[0]]=int(c[1]); dt_alias[c[0]]=c[3] if len(c)>3 else ''

pairs={'Pair11':('PD45792','PD45793'),'Pair13':('PD45794','PD45795'),'Pair21':('PD45798','PD45799'),
'Pair24':('PD45800','PD45801'),'Pair25':('PD45802','PD45803'),'Pair28':('PD45804','PD45805'),
'Pair31':('PD45806','PD45807'),'Pair38':('PD45808','PD45809'),'Pair40':('PD45810','PD45811'),
'Pair41':('PD45812','PD45813')}
inpair={pd:p for p,ds in pairs.items() for pd in ds}

def find_file(keys, want, prefer=None):
    cands=[(b,r,n) for (b,r,n) in idx if any(k and k in b for k in keys)]
    if prefer: cands=[c for c in cands if prefer in c[1]] or cands
    if not cands: return None
    cands.sort(key=lambda x: abs(x[2]-want)); return cands[0]

TSV_HEADER="donor\tproj\tds\treadlen\tmapped\tassembly\tsample\n"
TSV_PENDING="# PENDING — populated from the farm per-tip header probe (see patients/PLAN.md); placeholder only\n"

plan_txt=None; kept_tsv={}
if os.path.isdir(OUT):
    pf=os.path.join(OUT,'PLAN.md')
    if os.path.exists(pf): plan_txt=open(pf).read()   # preserve hand-written PLAN across rebuilds
    for tp in glob.glob(os.path.join(OUT,'*','*','colonies.tsv')):  # preserve POPULATED tsvs
        body=open(tp).read()
        if len([l for l in body.splitlines() if l and not l.startswith('#') and not l.startswith('donor\t')])>0:
            kept_tsv[os.path.relpath(tp,OUT)]=body
    shutil.rmtree(OUT)
os.makedirs(OUT,exist_ok=True)
if plan_txt is not None: open(os.path.join(OUT,'PLAN.md'),'w').write(plan_txt)
patients=[]; done=set(); missing=[]
for d,info in rows.items():
    if d in done: continue
    if info['tips']==0 and d not in inpair: continue
    org=organ(info['tissue'],info['unit'])
    if d in inpair:                                   # Chapman pair = one patient
        p=inpair[d]; d1,d2=pairs[p]; done.update([d1,d2])
        tips=dt_tips.get(p,0)
        ff=find_file([p], tips, prefer='chapman2024hsct')
        pid=f"{p}_{d1}_{d2}"; donors=[d1,d2]; files=[ff] if ff else []
    elif d in ('PD48367','PD48372'):                  # liver explant = all region trees
        files=[(os.path.basename(x),'chapman2025/SN',leaves(x)) for x in sorted(glob.glob(os.path.join(ROOT,'chapman2025/SN',f'tree_{d}?.tree')))]
        pid=d; donors=[d]; tips=info['tips']; done.add(d)
    else:
        keys=[d, info['alias'] if info['alias'] not in('','-') else None, dt_alias.get(d,'') or None]
        ff=find_file(keys, dt_tips.get(d,info['tips']))
        pid=d; donors=[d]; tips=dt_tips.get(d,info['tips']); files=[ff] if ff else []; done.add(d)
    # write folder
    pdir=os.path.join(OUT,org,pid); os.makedirs(pdir,exist_ok=True)
    for fb in files:
        shutil.copy(os.path.join(ROOT,fb[1],fb[0]), os.path.join(pdir,fb[0]))
    with open(os.path.join(pdir,'colonies.tsv'),'w') as t: t.write(TSV_HEADER); t.write(TSV_PENDING)
    if not files:
        open(os.path.join(pdir,'NO_TREE_ON_DISK.txt'),'w').write(
            f"{pid}: {tips} tips per donor_trees.tsv, but no .tree file downloaded yet.\nFetch — see patients/PLAN.md.\n")
        missing.append((org,pid,tips))
    patients.append(dict(org=org,pid=pid,donors=donors,tips=tips,nfiles=len(files),
        paper=info['paper'],diag=info['diag']))

# Stomach (Coorens 2025) — NOT in donor_table.tsv (tracked as a new-organ candidate); trees on disk.
for fp in sorted(glob.glob(os.path.join(ROOT,'coorens2025stomach','PD*.tree'))):
    b=os.path.basename(fp); pid=b[:-5]; tips=leaves(fp)
    pdir=os.path.join(OUT,'gastric_epithelium',pid); os.makedirs(pdir,exist_ok=True)
    shutil.copy(fp,os.path.join(pdir,b))
    with open(os.path.join(pdir,'colonies.tsv'),'w') as t: t.write(TSV_HEADER); t.write(TSV_PENDING)
    patients.append(dict(org='gastric_epithelium',pid=pid,donors=[pid],tips=tips,nfiles=1,
        paper='Coorens 2025',diag=''))

# restore any populated colonies.tsv snapshotted before the wipe
for rel,body in kept_tsv.items():
    dst=os.path.join(OUT,rel)
    if os.path.isdir(os.path.dirname(dst)): open(dst,'w').write(body)

# per-organ summary md
byorg=collections.defaultdict(list)
for p in patients: byorg[p['org']].append(p)
for org,recs in byorg.items():
    recs.sort(key=lambda r:-r['tips'])
    tot=sum(r['tips'] for r in recs)
    with open(os.path.join(OUT,org,f'{org}.md'),'w') as m:
        m.write(f"# {ORGAN_LABEL.get(org,org)} — {len(recs)} patients, {tot} tips total\n\n")
        m.write("Single-cell / single-stem-cell-derived WGS. TSVs are placeholders pending the farm probe (see [PLAN](../PLAN.md)).\n\n")
        m.write("| patient | PD id(s) | tips | tree |\n|---|---|---|---|\n")
        for r in recs:
            tstat='✓' if r['nfiles']>0 else '**missing**'
            pdids=', '.join(r['donors'])
            m.write(f"| {r['pid']} | {pdids} | {r['tips']} | {tstat} |\n")

json.dump(dict(n=len(patients),norgan=len(byorg),missing=missing),open('/tmp/pt/built.json','w'))
print(f"built {len(patients)} patients across {len(byorg)} organ groups at {OUT}")
print(f"trees missing on disk: {len(missing)}")
for org in sorted(byorg): print(f"  {org:24s} {len(byorg[org])} patients, {sum(r['tips'] for r in byorg[org])} tips")
