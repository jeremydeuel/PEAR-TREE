# 12. Tool implementations compared at code level (TraFiC, v‑TraFiC, MEIGA‑SR/LR, xTEA, ArtifactsFinder, MEIsimulator)

*Sources: the **actual source code** of seven repositories, read in full (not the papers):
[TraFiC](https://gitlab.com/mobilegenomesgroup/TraFiC),
[v‑TraFiC](https://gitlab.com/mobilegenomesgroup/v-trafic),
[MEIGA‑LR](https://gitlab.com/mobilegenomesgroup/MEIGA_LR),
[MEIGA‑SR](https://gitlab.com/mobilegenomesgroup/MEIGA-SR),
[xTEA](https://github.com/parklab/xTEA),
[ArtifactsFinder](https://github.com/lilicai/ArtifactsFinder),
[MEIsimulator](https://gitlab.com/mobilegenomesgroup/MEIsimulator).
Sections [§3](03_detecting_true_events.md)/[§6](06_delly_review.md) described these tools from the
literature; this section reports what the code actually does, corrects a few things the papers
imply, and extracts concrete, adoptable mechanisms for PEAR‑TREE. File:line citations are to the
cloned repositories. Complements the empirical work in [§10](10_empirical_discovery_patterns.md) and
the cross‑dataset/mouse‑ERV findings in [§11](11_cross_dataset_and_mouse_erv.md).*

---

## 12.1 Why this section exists

The rest of the review reasoned about MELT/xTea/MEIGA/TraFiC from their publications. Reading the
code changed several conclusions and surfaced **specific, portable techniques** that are cheaper and
more concrete than the literature-level recommendations. Three findings recur across *independent*
codebases and are the most important takeaways:

1. **The "both‑ends‑soft‑clipped read" (SMS) filter is the universal, reference‑free palindrome/
   chimera guard.** MEIGA, v‑TraFiC and xTEA all reject a read whose CIGAR is soft/hard‑clipped at
   *both* ends. This is a one‑line drop‑in that catches most of what a reference inverted‑repeat/
   palindrome blacklist ([§4.2.1](04_artefacts_in_discovery.md), R7) is for — PEAR‑TREE does not do it.
   This becomes new recommendation **R15**.
2. **Artefact‑region masking is done from the *local pileup*, not a static blacklist.** v‑TraFiC,
   MEIGA and xTEA each recompute, per candidate locus, the fraction of low‑MAPQ reads and/or the
   fraction of SMS reads (and, in xTEA, local depth vs genome median) and drop the locus if it looks
   like a multimapping/palindromic/pile‑up region. This is exactly the masking PEAR‑TREE has
   **disabled** (R1) — and it needs no blacklist file.
3. **Both mature short‑read callers end in a trained classifier** — MEIGA a logistic regression over
   ~29 features, xTEA a random forest over 15 features. PEAR‑TREE has none. The classifiers' feature
   sets are a ready‑made checklist of "what evidence should decide a call," and their fitted weights
   (extracted in §12.4) confirm the biology: **poly‑A and clip support are the strongest positive
   signals; both‑end‑clip‑rich regions and normal‑sample support are the strongest negatives.** This
   becomes new recommendation **R18** — with a mouse‑ERV caveat from [§11](11_cross_dataset_and_mouse_erv.md).

The rest of the section gives the per‑tool detail, then a master comparison matrix (§12.7), the
ArtifactsFinder blacklist recipe (§12.8), and the MEIsimulator recall benchmark (§12.9). Everything
here is folded into the recommendations in [§8](08_synthesis_and_recommendations.md).

---

## 12.2 TraFiC‑mem and v‑TraFiC (discordant‑pair callers, the opposite design)

TraFiC is a Perl/bash + Snakemake pipeline; v‑TraFiC is its viral‑integration fork. Both are
**discordant‑read‑pair** callers — the philosophical inverse of PEAR‑TREE — and reach base‑pair
resolution only in an *optional* local‑assembly step (TEIBA).

**Discovery (TraFiC).** `samtools.sh` collects discordant FR pairs (`samtools view -F1792` then keeps
flag values `81,161,97,145,65,129,113,177`); the **anchor** is the mate with MAPQ≠0 and the
**candidate element read is explicitly the MAPQ==0 (multi‑mapping) mate**, which is then `bwa mem`‑mapped
to a **56‑sequence TE consensus library** (`sequences.fa`; L1HS split into `L1HS_5`/`L1HS_3`, 22 Alu
subfamilies, SVA, ERVK). This is the exact opposite of PEAR‑TREE, which never trusts a MAPQ0 mate.
Reads are clustered by single‑linkage merge within **200 bp** of the same TE family
(`parser3_TE.pl:46`); a call requires a **positive cluster** met within **≤200 bp** by a **negative
cluster** on the same family (`script005_v3TEs.pl:38`), with **≥4 reads/cluster**
(`tumour_cluster_MGE.sh`). The ≤200 bp reciprocal gap absorbs both TSD and insert‑size slop; TraFiC
does **not** resolve TSD geometry at discovery.

**TraFiC's filter cascade** (all in `tx_pancancer.sh`) is a chain of database subtractions, and is the
mature version of what PEAR‑TREE does with its cohort genotyping:

| Filter (script) | What it removes | Threshold |
|---|---|---|
| `cleaner_v5.pl` | matched‑normal cluster of same family near the call | normal ≥3 reads, ≤200 bp |
| `cleaner_v5_v2.pl` | satellite/simple‑repeat regions (`Satellites.txt`, 3,684 regions) | interval overlap |
| `TEs_poly_cleaner.pl` | **breakpoint inside a *young* reference repeat of the same family** (RepeatMasker `RM.tab`, 5.3 M rows) | **divergence ≤20%**, family match, ±100 bp |
| `parseMGE4C.py` | known germline MEI polymorphism (tabix `pancan.large.bed`, **515,446** family‑labelled entries) | exact overlap **+ family‑name match** |
| PoN | panel‑of‑normals recurrent MEIs (`db_04102015.tab`, 182,269) | ≤200 bp / ≥3 |

**Transductions (TraFiC).** Both mates MAPQ ≥37, same‑chr separation >10 kb → cluster → trace the
donor coordinate back and link it to an **active‑source‑L1 catalog** (`master_db.txt` = 124 validated
active L1 sources; `master.txt` = 311 reference full‑length L1HS) within **20 kb**, with a **500 kb
self‑match guard**. This is a complete, portable 3′‑transduction source‑tracing recipe.

**TEIBA breakpoint analyzer.** Per insertion: pull cluster reads → **velvet local assembly (k=21,
`-min_contig_lgth`/auto‑cov)** → **BLAT** contigs vs reference and family consensus → derive poly‑A/TSD
at bp resolution → filter. Notably it uses **family‑specific assembly‑score floors**: L1 = 2, **Alu = 4,
SVA = 4, ERVK = 5** — i.e. non‑L1 families must clear a higher bar. (The Python that extracts poly‑A/TSD
base‑by‑base, `insertionBkpAnalysis.py`, is not committed to the repo, so those specifics are inferred
from the wrapper contract.)

**v‑TraFiC adds exactly the pieces closest to PEAR‑TREE's gaps.** The viral fork keeps the reciprocal‑
cluster architecture but bolts on:
- **A clip/unmapped‑mate discovery arm** (`samtools view -F1792 -f4` collects unmapped reads with a
  mapped mate — the viral read is the *unmapped* mate) alongside a relaxed‑MAPQ (`<20`) discordant arm.
- **Read‑level quality gates inside the discovery awk**: drop reads containing `N`; **drop both‑ends‑
  clipped reads** via CIGAR `$6 !~ /^[0-9]*S[0-9]*M[0-9]*S/` (the SMS guard — see §12.6); exclude
  `hs/GL/NC`/Y contigs.
- **Low‑complexity masking wired directly into discovery**: NCBI **`dustmasker`** on the candidate‑read
  FASTA, then drop reads that are wholly low‑complexity — *before* clustering. This is precisely
  PEAR‑TREE's "`is_low_complexity` exists but isn't wired into discovery" gap (R2).
- **A local‑pileup artefact mask** (`filterSMSclusters.sh`): re‑pull the BAM over each candidate cluster
  and drop it unless **<30% of reads have MAPQ<10** *and* **<15% of reads are SMS (both‑ends‑clipped)**.
  Reference‑free; catches multimapping/palindromic loci with no blacklist.

## 12.3 MEIGA‑SR / MEIGA‑LR (clip + discordant, ending in a logistic‑regression classifier)

MEIGA‑SR and MEIGA‑LR are Python and share the `GAPI/` library. **Important framing correction:** the
shipped, runnable MEIGA‑SR (`demo/`, README) is a **targeted VAF genotyper** (`scripts/MEIs_vaf.py`,
`RG_vaf.py`) of *previously identified* events; the de‑novo discovery caller (`modules/caller.py`) plus
its ML classifier is present but most of its numeric thresholds come from an **unshipped INI config**
(the concrete numbers below are hardcoded literals or the MEIGA‑LR argparse defaults, which are the best
proxy).

**Discovery signals.** Per‑read gate (`bamtools.py:390`): skip unmapped/dup/low‑MAPQ, then
`filter_alignments(['duplicate','supplementary','mateUnmap','SMS'])`. Clip recruitment
(`collectCLIPPING`, `bamtools.py:458`): soft‑clip = CIGAR op 4, hard‑clip = op 5; LEFT clip at
`reference_start`, RIGHT clip at `reference_end`, each requiring `minCLIPPINGlen`. Discordant recruitment
(`collectDISCORDANT`, `bamtools.py:511`): `not is_proper_pair`, with an insert‑size floor (default
**5000 bp**). **A single split read can substitute for a discordant pair**: `events.py::SA_as_DISCORDANTS`
converts an SA‑tag supplementary clip into a `pseudo=True` DISCORDANT object — a lightweight bridge that
fits a clip‑first design without a full discordant engine.

**Clustering.** Reciprocal‑overlap clustering (`modules/clustering.py`) for the MEI path;
orientation‑aware distance thresholds `equalOrientThreshold` vs the wider `oppOrientThreshold`
(`clustering.py:416`); metaclusters at **200 bp**; right‑clip→PLUS, left‑clip→MINUS, paired into a
metacluster kept if `supportingReads ≥ minReads`.

**Consensus.** No de‑novo assembler — a custom **minimap2‑overlap layout → racon polish → MUSCLE MSA →
EMBOSS `cons`** pipeline (`assembly.py`); the reconstructed insert is re‑aligned to a
retrotransposon consensus DB; `IS_FULL` if `percConsensus ≥ 95`.

**Poly‑A / TSD / EN‑motif — where PEAR‑TREE's plan *leads*.** Poly‑A is detected with tuned monomer
parameters (Illumina clip path: `window 8, maxWindowDist 2, minMonomer 8, purity 95%, ≤1 bp from end`;
insert path: `minMonomer 10, purity 80%, ≤10 bp from end`). But **TSD detection is a commented‑out stub
with no length window** — TSD is only inferred at breakpoint level (`plus_bkp > minus_bkp`), and there is
**no endonuclease‑motif (TTAAAA) code anywhere** in either repo. So PEAR‑TREE's planned explicit
TSD‑length and L1‑EN‑motif scoring (R8) would put it *ahead* of MEIGA, not behind.

**Filter set (`modules/filters.py`, SR).** ~20 filters; the ones PEAR‑TREE lacks or under‑uses:

| Filter | Threshold |
|---|---|
| `filter_discordant_mate_unspecific` | drop if `nbDiscordant/nbProperPair > 0.95` (region captures everything) |
| `area` / AREAMAPQ / AREASMS | drop if `percMAPQ ≥ maxRegionlowMQ` **or** `percSMS ≥ maxRegionSMS` (local pileup, bkp ±100) |
| `filter_clusterRange_reciprocalMeta` | **KS test** of PLUS vs MINUS breakpoint positions; drop if p > 0.05 (positions don't co‑localise) |
| `filter_germline_MEI` | drop if overlaps a known germline MEI of the **same identity** within 150 bp |
| `filter_normal_VAFs` | drop if normal `VAF > 0.05` or `cREF < 5` |
| (MEIGA‑LR) `filter_low_complexity` | poly‑A>0.4 auto‑pass; else drop if `Σ_K|freq−0.25^K| > 0.83` over K∈{1,2,3,4} |
| (MEIGA‑LR) `filter_expanded_repeat` | microsatellite: drop if `|obs−0.25^len| > 0.2` |
| (MEIGA‑LR) `filter_len` | length caps: L1 > 7000, Alu > 500, SVA > 3000 → drop |

`bins_lowMAPQ_SMS.py` is the tool that characterises artefact regions: per window (±25 bp) it emits
`percLowMAPQ` (MAPQ≤1) and `percSMS` (both‑ends‑clipped), which the `area` filter then thresholds.

## 12.4 MEIGA's classifier and xTEA's genotyper (the ML layer PEAR‑TREE lacks)

**MEIGA‑SR** ends in an sklearn `Pipeline` (ColumnTransformer scaling/one‑hot → **LogisticRegression**,
L2, `C=1.0`, class‑weighted `{False:1.16, True:0.88}`), applied to L1/poly‑A rows
(`modules/classifier.py`; model in `modules/meiga.joblib`). Its **29 features** and the **fitted
coefficients** (unpickled) tell you what MEIGA learned distinguishes a true MEI from an artefact
(+ pushes toward *true*):

```
-1.29  svsNormal          (nearby SVs in the normal → artefact)
+0.90  nbCLIPPING         (clip support → real)
-0.86  nbNormal           (normal support → germline/artefact)
+0.79  minus_pA           (poly-A on the 3' side → real)
+0.79  tumourPerc / -0.79 germPerc   (tumour-specific → real)
+0.77  pA                 (poly-A present → real)
-0.72  areaSMS            (both-end-clip-rich region → artefact)
+0.71  src_end=3          (3' transduction signal → real)
+0.55  ins_type=TD2 ...   +0.32 minus_id=L1, +0.32 plus_pA, +0.27 nbDISCORDANT ...
```

Read plainly: **poly‑A (three features) and clip support are the strongest positive evidence; SVs/
support in the normal and SMS‑rich regions are the strongest negatives.** This is quantitative
confirmation of the whole review's thesis — and a template PEAR‑TREE can adapt by swapping MEIGA's
matched‑normal features (`svsNormal`, `nbNormal`, `germPerc`) for **cohort/tree‑derived features**
(cross‑donor sharing fraction, tree‑consistency, per‑branch recurrence). **Mouse‑ERV caveat
([§11](11_cross_dataset_and_mouse_erv.md)):** MEIGA's model is L1/poly‑A‑centric and weights poly‑A
heavily; a PEAR‑TREE classifier trained on a cohort that includes **mouse ERV/LTR insertions (which make
no poly‑A)** must gate poly‑A features on element class, or it will systematically down‑weight every true
ERV insertion. Train per element class, or include an "expects‑poly‑A" indicator.

**xTEA** ends in a `RandomForestClassifier(n_estimators=20)` over **15 features**
(`x_genotype_classify_sklearn.py`), almost all **coverage‑normalised**: left/right clip‑on‑consensus /
depth, left/right disc‑on‑consensus / depth, `polyA/cov`, `clipratio = clip/(clip+fullmap)`,
`discratio`, raw clip / depth, disc/concordant per depth. PEAR‑TREE's quality‑weighted match score is a
plausible drop‑in for these raw normalised counts; the feature *set* is a good sanity checklist.

**xTEA's other distinctive mechanisms** (all things PEAR‑TREE lacks):
- **Coverage‑adaptive count thresholds** (`x_parameter.py`): the minimum clip/disc/clip+disc counts
  **scale with local coverage** via a lookup table (Illumina germline: cov 5→(1,3,0), 30→(3,4,1),
  100→(8,12,3), 300→(25,30,8)), with a lower case‑control table for somatic. PEAR‑TREE uses fixed counts.
  This becomes new recommendation **R16**.
- **Clip↔disc geometric consistency** (`_is_distance_consistency`, `x_clip_disc_filter.py:1231`): a
  left‑clip must be corroborated by **right‑side discordant** reads whose consensus mapping is within one
  insert‑size; "both‑side clips but no discordant support" is an explicit reject. This is xTEA's strongest
  FP filter and depends on having a discordant leg.
- **High‑coverage‑island mask**: local depth computed in two windows (200 bp focal, 900 bp island); any
  locus with depth **> 3× the sample median** (`MAX_COV_TIMES=3`, sample median from 3,000 random sites)
  is moved to `HIGH_COV_ISD` and skipped. This is precisely PEAR‑TREE's disabled masking (R1).
- **Centromere/low‑mappability BED blacklist** (`XBlackList`, IntervalTree) + panel‑of‑normals/gnomAD‑SV
  subtraction in the somatic path.
- **RepeatMasker low‑divergence same‑family reject**: drop a call inside a reference copy of the *same*
  family with divergence < **15%** (`REP_DIVERGENT_CUTOFF`) — same idea as TraFiC's `TEs_poly_cleaner.pl`
  (new recommendation **R17**).
- **Short‑clip poly‑A rescue**: keep a 7–9 bp clip *only if* it is pure poly‑A/T (`MINIMUM_POLYA_CLIP=7`).
- **Orientation‑aware poly‑A** (`is_consecutive_polyA_T_with_oritation`: left‑clip→poly‑A, right‑clip→
  poly‑T) — coded but, like PEAR‑TREE's, only partly wired; and a **TSD confidence tier** (`TWO_SIDE_TPRT_BOTH`
  needs poly‑A **and** `0 < TSD < 100`).
- xTEA keeps clips at **MAPQ ≥ 12** (not PEAR‑TREE's ≥40) and controls FP via the low‑MAPQ‑clip ratio
  (`MAX_LOWQ_CLIP_RATIO=0.65`) plus the disc cross‑check — evidence that a hard MAPQ≥40 floor is stricter
  than necessary if the other guards are present.

## 12.5 ArtifactsFinder (the IVR/palindrome blacklist generator) — with a correction

ArtifactsFinder implements R7's reference inverted‑repeat/palindrome blacklist, but **not** the way the
literature summary in [§4.2.1/§5.1](04_artefacts_in_discovery.md) implied:

- It does **not** emit a merged interval BED. It takes a BED of target regions, and for each writes out
  the **individual base positions** where a near‑palindrome/near‑IVR carries a mismatch, as pseudo‑variant
  rows `chr start end ref alt` — i.e. artefact‑prone *positions*, meant to be subtracted from a variant
  caller's VCF. Collapsing these into padded intervals is left to the user (bedtools).
- **Real parameters differ from the paper's headline numbers.** `ArtifactsFinderIVS.py` (inverted repeats):
  arm‑pair total **≥8 bp (`D_LEN`)**, spacer/loop **≥5 bp (`STEM_LEN`)**, sub‑arm **≥2 bp (`S_LEN`)**,
  **±50 bp** flank, tolerating ≤2‑bp indels/SNPs between arms. `ArtifactsFinderPS.py` (palindromes):
  extends arms allowing **≤1 mismatch (`GAP=1`)**, ±50 bp flank — but its **minimum‑length gate is
  commented out** (`LEN=17` is dead code), so PS emits *every* short ≤1‑mismatch palindrome center unless
  you re‑enable a length filter. `bedReform.py` only **chunks** long regions into 200 bp windows (the
  finder is ~O(N²) per window); it does no merging or padding.

Given §12.6, the read‑level SMS filter (R15) is a much cheaper first line of defence than genome‑wide
ArtifactsFinder blacklisting; the blacklist is worth building only if residual palindrome FPs survive.
The full recipe is in §12.8.

## 12.6 The convergent finding: the SMS "both‑ends‑clipped read" filter (→ R15)

Three independent codebases reject the same read‑level artefact:

| Tool | Rule | Location |
|---|---|---|
| MEIGA‑SR/LR | CIGAR op 4/5 at **first and last** position → `SMS`, dropped | `GAPI/bamtools.py:654` |
| v‑TraFiC | CIGAR `^[0-9]*S[0-9]*M[0-9]*S$` in discovery awk → dropped | `src/bash/*` |
| xTEA | both‑side clips longer than `MAX_CLIP_CLIP_LEN=8` → dropped | `clip_read.py:701,805,920` |

A read soft‑clipped on **both** flanks is the read‑level fingerprint of a reference palindrome/inverted
repeat, an adapter‑dimer, or a spurious multimap — exactly the structure‑specific chimera of
[§5.1](05_sequencing_artefacts.md). PEAR‑TREE currently classifies such a read by its *longer* clip
(`discovery.py`, "if both ends clipped, the longer clip wins") and keeps it. **Adding an SMS reject is a
one‑line, reference‑free filter that three mature callers consider essential** — the cheapest item in this
whole review (recommendation R15). It is complementary to (and cheaper than) the ArtifactsFinder blacklist
(R7). **Note the mouse caveat ([§11](11_cross_dataset_and_mouse_erv.md)): mouse WGS carries ~10× more
palindrome artefacts than human, so R15 is higher‑priority for mouse cohorts — but the SMS test keys on
CIGAR shape, not on "is a repeat," so it does not touch genuine repeat‑consensus clips and is safe for the
mouse ERV events.**

Its natural partner is the **local‑pileup artefact mask** — recompute, per candidate locus, the fraction of
low‑MAPQ reads and the fraction of SMS reads, and drop the locus above a threshold. v‑TraFiC uses
**MAPQ<10 ≤30% and SMS ≤15%**; MEIGA uses configurable `maxRegionlowMQ`/`maxRegionSMS` over `bkp±100`;
xTEA uses **depth > 3× median**. None needs a blacklist file. This is the concrete, adoptable form of
re‑enabling PEAR‑TREE's disabled masking (R1).

## 12.7 Master comparison matrix

| Axis | PEAR‑TREE | TraFiC‑mem | v‑TraFiC | MEIGA‑SR | xTEA | Delly ([§6](06_delly_review.md)) |
|---|---|---|---|---|---|---|
| Primary signal | **clip/split** | discordant (multimap mate) | discordant + unmapped‑mate | clip + discordant | clip + discordant (parallel) | discordant→split |
| Discordant leg | **✗** | ✓ | ✓ | ✓ (+SA→pseudo‑disc) | ✓ (clip↔disc consistency) | ✓ |
| Clip/split leg | ✓ (primary) | only in TEIBA assembly | ✓ (added) | ✓ | ✓ | ✓ (refine) |
| TSD geometry | **✓ (L<R, core)** | ≤200 bp gap only | ≤350 bp gap | bkp‑level (stub, no length) | length‑bounded presence (<100) | ✗ |
| Poly‑A positive score | presence + rescue | merge/subtract track | as TraFiC | detected, top ML weight | detected, ML feature | ✗ |
| TSD‑length / EN‑motif score | **planned (R8) — would lead** | ✗ | ✗ | **✗** | TSD tier, **no EN motif** | ✗ |
| Both‑end‑clip (SMS) reject | **✗ (→ R15)** | ✗ | **✓** | **✓** | **✓** | — |
| Local‑pileup artefact mask | **✗ (disabled)** | ✗ | **✓ (MAPQ/SMS)** | **✓ (area filter)** | **✓ (3× depth)** | read cap (‑x) |
| Low‑complexity in discovery | **✗ (unwired)** | ✗ | **✓ (dustmasker)** | ✓ (LR k‑mer) | ✓ (short‑clip gate) | ✗ |
| Coverage‑adaptive thresholds | ✗ (fixed → R16) | ✗ | ✗ | ✗ | **✓ (table)** | per‑library insert model |
| RepeatMasker self‑mask (young copy) | ✗ (→ R17) | **✓ (≤20% div)** | ✗ | annot only | **✓ (<15% div)** | ✗ |
| Reference blacklist | `len(name)>5` proxy | satellites + PoN + germline DB | + SMS/MAPQ mask | germline‑MEI BED | centromere BED + PoN | telomere/centromere ‑x |
| Transduction source tracing | ✗ | **✓ (124 sources, 20 kb)** | viral DB (BLAST) | ✓ (15 kb src DB) | **✓ (1 kb, dominant‑source)** | ✗ |
| Germline/somatic basis | **cohort/tree genotyping** | matched normal | matched normal | matched normal + germline DB | matched normal / case‑ctrl | matched normal |
| Final classifier | ✗ (quality score → R18) | ✗ | ✗ (BLAST tiers) | **✓ logistic reg (29 feat)** | **✓ random forest (15 feat)** | ✗ |
| MAPQ floor (clip) | **≥40** | ≠0 (multimap mate) | <20 disc | config (≈15) | **≥12** | ≥20 |

**Reading of the matrix.** PEAR‑TREE is the only tool with TSD‑first geometry as a *core* discovery signal
and the only one using cohort/tree genotyping instead of a matched normal — genuine strengths. Its
consistent *omissions* relative to the field are the right‑hand rows marked ✗: **no discordant leg, no
SMS reject, disabled local masking, unwired low‑complexity, no coverage‑adaptive thresholds, no
self‑mask, no final classifier.** Encouragingly, its *planned* work (TSD‑length + EN‑motif positive
scoring, R8) is an area where every other tool is weak or absent, so it is a differentiator rather than a
catch‑up.

## 12.8 R7 recipe — a genome‑wide palindrome/IVR blacklist from ArtifactsFinder

ArtifactsFinder is panel‑scoped; to make a PEAR‑TREE reference blacklist for T2T/hg38:

```bash
# deps: python3, pyfasta, pysam, bedtools; samtools faidx ref.fa
# 1. tile the genome; chunk into 200 bp windows the O(N^2) finder can handle
awk 'BEGIN{OFS="\t"}{print $1,0,$2}' ref.fa.fai > genome.bed
python3 bedReform.py genome.bed                       # -> genome.reform.bed

# 2. run both finders (shard by chromosome / array job)
python3 ArtifactsFinderIVS.py -g ref.fa -b genome.reform.bed -l 50   # -> *.out (arm pairs)
python3 ArtifactsFinderPS.py  -g ref.fa -b genome.reform.bed -l 50 -p PS  # -> PS.backlist.txt

# 3. convert to intervals
cut -f1-3 PS.backlist.txt > ps_sites.bed
awk 'BEGIN{OFS="\t"}{split($3,a,"_");split($4,b,"_");
     print a[1],a[2],a[3]; print b[1],b[2],b[3]}' genome.reform.bed.out > ivs_arms.bed

# 4. merge + pad (small pad: the finder already used a 50 bp structural flank)
cat ps_sites.bed ivs_arms.bed | sort -k1,1 -k2,2n \
  | bedtools slop -b 15 -g ref.fa.fai | bedtools merge > palindrome_IVR_blacklist.bed
```

Apply by dropping/flagging any candidate whose clip‑anchor breakpoint intersects the blacklist. **Before
shipping, re‑enable a PS length gate** (`LEN ≥ 17`) or PS floods short 2–4 bp palindromes; keep the final
`slop` small (10–20 bp) so real neighbouring insertions are not masked. Cost: ~15 M windows genome‑wide,
so parallelise. Given §12.6, treat this as a *second line* behind the read‑level SMS filter (R15). Note
([§11](11_cross_dataset_and_mouse_erv.md)) the mouse genome's higher palindrome load makes a **per‑species**
blacklist (mouse GRCm39 as well as human T2T) worthwhile if mouse cohorts are analysed.

## 12.9 §8.4 recall benchmark — MEIsimulator as the spike‑in engine

The validation plan in [§8.4](08_synthesis_and_recommendations.md) had no spike‑in generator.
MEIsimulator is exactly that: it spikes labelled MEIs into a reference, simulates reads, and emits a
**truth report** (`#CHR,BEG,END,STRAND,MEI,MEI_LENGTH,MODE,rt_size,td_size,td_coord,polyA_size,tsd_size,…`).

**What it models.** Solo L1/Alu/SVA insertions, 3′ partnered/orphan and 5′ transductions (with a database
of 113 named source‑L1 loci laid out as 5 kb up + element + 50 bp polyA + 5 kb down), L1‑mediated
deletions/duplications/inversions/translocations, and viral/generic integrations. TPRT hallmarks:
**5′ truncation drawn from an empirical L1 length distribution** (`rtLen_prob_PCAWG.txt`, 30–6023 bp);
**poly‑A from an empirical distribution** (`polyAdist_prob_PCAWG.txt`, 0–105 bp, mode ~60) *for
transductions* (solo insertions use a fixed length — a realism gap); sequence divergence via `mutate_dna`.
**VAF model:** it builds a MEI‑bearing genome at `cov·clonality` and an unmodified genome at
`cov·(1−clonality)`, simulates reads from each with `art_illumina` (HS25, 150 bp PE, insert 350±10), and
**`samtools merge`s** them — the mixing ratio *is* the VAF. RNG is seeded (`random.seed(20)`), so truth
sets are reproducible across tool versions.

**Benchmark for PEAR‑TREE.** Run MEIsimulator to make a merged BAM + truth report; run PEAR‑TREE; intersect
calls vs truth (`BEG_mockcoord` after SV shifts, ±(TSD + clip tolerance)); then sweep:
- **VAF** (`--clonality` 5/10/25/50/100%) → the key recall‑vs‑VAF curve for clip‑first sensitivity;
- **coverage** (15/30/60×);
- **element size / 5′ truncation** (stratify recall by `rt_size`/`MEI_LENGTH` — tests short vs full‑length L1);
- **event class** (solo vs partnered/orphan transduction vs SV) — tests poly‑A/TSD/transduction handling.

**Caveats** (document these): TSD is hardcoded to a **fixed 15 bp** (no length distribution), solo‑insertion
poly‑A is fixed, and there is **no L1‑EN target motif and no 5′ twin‑priming inversion** — so MEIsimulator
will *not* exercise PEAR‑TREE's EN‑motif or inverted‑truncation logic (validate those on real WGS or the S1
LINE‑call reanalysis), and TSD‑length metrics should be read with caution unless the simulator's TSD is
patched to vary. MEIsimulator is human‑L1‑centric and ships no ERV/LTR model, so **mouse ERV recall
([§11](11_cross_dataset_and_mouse_erv.md)) must be benchmarked from the real mouse cohort, not this
simulator.**

---

## 12.10 What changed in the review because of the code

- **R6/R7 reframed.** The single highest‑value, lowest‑effort palindrome/chimera defence is the
  **read‑level SMS filter (R15)** (three tools agree), not the ArtifactsFinder blacklist. The blacklist
  recipe (§12.8) is retained as a second line.
- **R1 made concrete.** Re‑enable masking as a **local‑pileup MAPQ/SMS/depth mask** (v‑TraFiC/MEIGA/xTEA),
  which needs no blacklist — cheaper than the review originally implied.
- **R2 made concrete.** v‑TraFiC (`dustmasker`) and MEIGA‑LR (k‑mer entropy `Σ|freq−0.25^K|>0.83` with a
  poly‑A escape hatch) give two ready algorithms for wiring low‑complexity rejection into discovery.
- **Three new recommendations** surfaced — **R15** (SMS both‑ends‑clip reject), **R16** (coverage‑adaptive
  count thresholds), **R17** (RepeatMasker self‑mask by divergence) — plus **R18** (a cohort‑adapted final
  classifier, using MEIGA/xTEA's feature schemas). See [§8](08_synthesis_and_recommendations.md).
- **R8 confirmed as a differentiator, not a catch‑up:** neither MEIGA nor xTEA scores the L1‑EN motif, and
  MEIGA's TSD detection is a stub — PEAR‑TREE's planned TSD‑length + EN‑motif scoring would lead the field.
- **A concrete §8.4 recall benchmark** now exists (MEIsimulator, §12.9).
- **ArtifactsFinder description corrected** in [§4](04_artefacts_in_discovery.md)/[§9](09_references.md):
  per‑base positions not intervals; real params D_LEN=8/STEM_LEN=5/±50 bp; PS length gate disabled.
- **Consistency with [§11](11_cross_dataset_and_mouse_erv.md):** the code review's poly‑A findings (MEIGA/
  xTEA weight poly‑A heavily as positive) are reconciled with the mouse‑ERV finding that poly‑A must never be
  *required* — R8 and R18 are both gated on element class so LTR/ERV insertions are scored by consensus
  identity + short TSD instead.
