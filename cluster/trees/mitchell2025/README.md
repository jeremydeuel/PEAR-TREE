# Mitchell et al. 2025 (Nat Genet) trees — the chemotherapy cohort

44 Newick trees, fetched 2026-07-17 from **Mendeley Data `10.17632/2fczcd49yj.1`** (v1,
published 2025-03-21). CC-BY. 216 KB.

## Get the right link — the paper names three, and two are decoys

The paper's Data availability section distinguishes them, and the distinction matters:

| Source | What is actually there |
|---|---|
| Zenodo `10.5281/zenodo.15235476` / github `emily-mitchell/chemotherapy` | **Code + "smaller derived datasets".** Contains exactly ONE donor's tree (PD37580, as `PX007_1`). Not the cohort. |
| **Mendeley `10.17632/2fczcd49yj.1`** | **"The main data needed to reanalyze and reproduce the results"** — all 44 trees + the per-donor annotated mutation matrices. **This is the one.** |
| EGA `EGAD00001015339` (WGS) / `EGAD00001015340` (NanoSeq) | Raw sequencing, restricted access. |

Our memory pointed at the Zenodo DOI. It is the code repo. Use Mendeley for trees.

## Fetching it

The Mendeley public API embeds the file listing in the dataset JSON — there is no working
`/files` endpoint (`/files?version=1` returns `{"error":400}`):

```bash
curl -s "https://data.mendeley.com/public-api/datasets/2fczcd49yj" -H 'Accept: application/json' -o mend.json
python3 -c "
import json
for f in json.load(open('mend.json'))['files']:
    if f['filename'].endswith('.tree'):
        print('curl -sL \"%s\" -o \"%s\"' % (f['content_details']['download_url'], f['filename']))
" | bash
```

The dataset is 691 MB, but the trees are 216 KB of it. The bulk is
`annotated_mut_set_*` (up to 230 MB each) — per-donor mutation matrices we do not need.

## Two file naming schemes, same trees

- `tree_<CODENAME>_<N>_01_standard_rho01.tree` — full phylogeny, `N` = timepoint.
- `<PDID>.tree` — usually the SAME tree; but for donors where the 2025 paper only used a small
  burden subset, this is the SUBSET, not the full tree. **PD50308.tree has 8 tips while
  tree_PX004_1 has 41.** Prefer `tree_*_rho01.tree` where both exist.

## Codename -> PD id, established here by IDENTICAL TIP SETS (not inference)

    KX001=PD40521  KX004=PD45534  KX009=PD49236  KX010=PD49237
    SX001=PD41048  PX001=PD47703  PX002=PD44579  PX003=PD50307
    PX004=PD50308  PX005=PD47537  PX007=PD37580

**AX001 = PD43976 — PROVEN.** `PD43976.tree`'s tips are named `BMH1_TG001_P31_A11` etc. The
tree file named for the PD id contains the codename's colonies. This closes a question we had
flip-flopped on twice.

**PD47703 (PX001) is SERIAL** — two timepoints, `PX001_1` (200 tips, age 48) and `PX001_4`
(59 tips, age 49). Our own 259-tip PD47703 tree is these two POOLED (200+59=259).

## Where these disagree with our own trees — and why it is not a conflict

For the Mitchell 2022 normals, the 2025 paper re-used only a subset, so the Mendeley tree is
SMALLER than the Mitchell 2022 tree we already hold. **Ours are the better ones — keep them:**

    donor     ours (2022)   Mendeley (2025 re-use)
    PD45534       922            200
    PD40521       407             41
    PD48402       367             10
    PD47738       315             10
    PD41048       362            200
    PD40315       195             10

Supp Table 1 explains it exactly: "Number colonies mutation burden analysis" = 10 for
PD40315/PD47738/PD48402, plus "Additional colonies for phylogeny" = 190 for PD41048/PD45534.

## THE ANOMALY — our PD49236/PD49237 trees are BIGGER than the published ones

    donor     ours    Mendeley/published (Supp T1)
    PD49236    422     173  (= 10 burden + 163 phylogeny)
    PD49237    231      44  (= 10 burden +  34 phylogeny)

These two donors were FIRST published in the 2025 paper ("Normal - unpublished" before it), so
there is no earlier, larger tree to have come from. Our trees have ~2.4x and ~5.2x more tips
than anything published. Either they are an unpublished extension Jeremy has access to, or our
tree ingest mixed something in. Worth resolving before either tree is used as truth.

Suggestively: PD49226/27/28/29 (345/424/487/723 tips) are ALSO unpublished, sit 7-10 ids away
from PD49236/37, and are our largest trees. An unpublished organ-donor set would explain both.

## Donors this adds that we had no tree for at all

PD37580(84) PD43976(10) PD47537(44) PD47540(5) PD47541(10) PD47699(10) PD47701(10)
PD50306(9) PD50308(41)

## Donors entirely new to our catalogue (12)

PD47536(5) PD47538(10) PD47539(5) PD47695(5) PD47696(5) PD47697(10) PD47698(10) PD47700(4)
PD47702(10) PD60009(9) PD60010(10) PD60011(9)

All are Mitchell 2025 peripheral-blood HSPC colonies. Note the tip counts: these are 4-10
colony trees — too small to be useful as phylogenetic truth for MEI calling, but they are
real chemo-exposed donors. See ../donor_clinical.tsv for each one's regimen.
