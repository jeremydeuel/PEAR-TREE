# Inventory — new-organ donors: what exists, where, and what is still unknown

**Status: INVENTORY ONLY. Nothing is staged. No cohort has been selected.**
Surveyed 2026-07-17 (`~/hunt_new_organs.txt` on farm22). This file records what we HAVE.
The decision about which samples to progress is deliberately NOT made here.

## Headline

**All 40 candidate donors are on the farm. None needs an EGA request.**
Manifest coverage in `nst_links`: **stomach 238/238, NF1 909/909 — zero missing.**

| cohort | donors | WGS samples | assembly | project to use | trees | usable? |
|---|---|---|---|---|---|---|
| **NF1 multi-organ** | 3 | 909 (838 LCM + 71 bulk) | GRCh38 | **2571** | **none exist** | **assay/read-len UNKNOWN** |
| **Stomach** | 30 | 238 | dual: GRCh37 + GRCh38 | see below | 30 fetched, validated | **YES — WGS / 151bp (probed)** |
| **Wilms kidney** | 7 | ~40 bulk + ~80 LCM | hs37d5_GRCh37 | 2939/2628/2541/1141 | newick in supplement | mostly TUMOUR |

## NF1 — PASSES THE ASSAY/READ-LENGTH GATE

Probed 2026-07-17 (`~/probe_nf1.txt`), 8 BAMs across all 3 donors, unanimous:
**`WGS_ILLUMINA_short` / 151bp / GRCh38 / 804M-4,098M mapped.** This cohort is usable.
151bp is the answer we needed — clip-based MEI discovery works.

### *** DO NOT "pin to project 2571". THAT WAS MY ERROR — here is the correction. ***

The glob-based step-1 survey reported PD51122 = **413** samples in 2571, which equals the
manifest's 413, and I called that proof the release was complete. **It was a coincidence of
two different quantities.** The MANIFEST-restricted count is **410**: the glob's 413 counted
~3 non-manifest samples while MISSING 3 manifest ones.

Manifest-restricted truth (`probe_new_organs.sh`):

| donor | manifest WGS | in 2571 | in 2789 |
|---|---|---|---|
| PD50297 | 305 | **305** (all) | — |
| PD51122 | 413 | **410** | 24 (21 dual-listed + **3 found NOWHERE ELSE**) |
| PD51123 | 191 | **191** (all) | — |

Reconciles exactly: rows 305+410+24+191 = 930; distinct resolved = 909; overlap = 21;
PD51122 distinct = 410+24-21 = **413** = manifest.

**So 2571 holds 906 of the 909. Three PD51122 WGS samples exist ONLY in project 2789** and
pinning to 2571 would silently drop them. The two 2789 samples probed are block `PD51122e`
= *Frontal cortex, left / BRAIN / NORMAL* — so the missing 3 may be part of the 305 normal
brain set. **They are NOT yet individually identified** (the probe only spot-checks 2 BAMs
per group); `probe_new_organs.sh` now writes the full per-sample map to resolve this.

**LESSON, and it is the same one as PD37590:** a count matching a count is not proof.
Two aggregates can agree while disagreeing member-by-member. Select on the sample->project
map, never on "the project whose total looks right".

**Still ignore these projects for this cohort:**
- `2796` (50+66+65 = **181**) = the **EXOME** arm (paper's WES table has 182 rows). Not WGS.
- `2882` (`_WGMS`), `3334`/`3329`/`3507` (`_tds`), `2745`/`3439`/`3161` (`_ds`) = NO_BAM
  **analysis releases**. Expected — see [[irods-cgp-zone-limits]].
- `2585` = ONE stray hs37d5 BAM for PD50297. Ignore.

Content of the 838 WGS_LCM (`cluster/manifests/oliver2025nf1.wgs.tsv`):

| organ_group | n | of which NORMAL |
|---|---|---|
| BRAIN | 498 | **305** |
| OTHER_VISCERA | 250 | |
| SPINAL_CORD | 36 | |
| PERIPHERAL_NERVE | 20 | |
| unknown site (`?`) | 34 | — written `?`, never defaulted |

193 of the 838 are glioma **TUMOUR**. The 305 normal-brain figure is
PD51122=155, PD51123=83, PD50297=67, across frontal/parietal/occipital cortex, cerebellum,
hippocampus, SVZ, pons, medulla, pineal, pituitary.

## Stomach — PASSES THE ASSAY/READ-LENGTH GATE (probed 2026-07-17)

`bash cluster/probe_new_organs.sh cluster/manifests/coorens2025stomach.wgs.tsv`
(`~/probe_stomach.txt`) — probed all 30 donors across both release projects, unanimous:
**`WGS_ILLUMINA_short` / 151bp / 776M-1,104M mapped.** DS never reads a panel, confirming
the manifest correctly excluded the 829-sample TGS arm. This cohort is usable.

Dual-released, and the panel trap. Same samples released twice, identical counts per pair:
`1911→2405`, `1933→2475`, `2134→2476`, `2398→2477` (left = hs37d5_GRCh37, right = GRCh38).
The paper is GRCh38. The probe confirms 1911=hs37d5, 2405=GRCh38, same sample basenames on
both — pick the GRCh38 (2405-series) release, or remap to hs1 as usual.

*** **PROJECT IS NOT A PROXY FOR ASSAY, AND NAME CANNOT SEPARATE WGS FROM PANEL.** ***
A donor glob returns **~87 samples for PD40293 but the paper ran only 12 as WGS**. The other
~75 are the **829-sample TARGETED PANEL**, named identically (`PD40293c_lo0003`). Per-donor
project counts do NOT split cleanly — PD41759 is 23 WGS but proj 1933 holds 17; PD42790 is
11 WGS but 1911 holds 18. **The only exact discriminator is the manifest: the WGS (238) and
TGS (829) name sets are DISJOINT, overlap 0.** Select by manifest name, never by glob or
project. `survey_assay.sh`'s header records a 10x10 run burned on three targeted panels that
returned 11/7/92 loci against 8,833 for real WGS.

Trees: 30 fetched to `coorens2025stomach/`, 239 real tips vs 238 manifest rows. But they are
**tiny — median ~6 tips/donor (2-23)**, below where our monophyly/carrier-count logic means
anything. `PD41759`(25,600 SNV) and `PD41762`(64,203) are hypermutator CANCER clones, not
normal epithelium. `PD41762`'s tree carries a non-sample tip named `Ancestral` — filter `^PD`.

## THE OPEN QUESTION — RESOLVED for both usable cohorts (2026-07-17)

Both NF1 and stomach have now been header-probed and both are **WGS / 151bp**:
- NF1: `~/probe_nf1.txt` — 8 BAMs, WGS_ILLUMINA_short / 151bp / GRCh38.
- Stomach: `~/probe_stomach.txt` — 30 donors × 2 releases, WGS_ILLUMINA_short / 151bp / dual.

Presence in `nst_links` still says **nothing** about usability, so the gate remains: for any
NEW cohort, probe assay (`@RG DS:` must be WGS) and read length (151bp fine; **75bp badly
compromises clip-based MEI discovery, remapping cannot fix it** — see
[[hg19-remap-to-grch38-plan]]) before progressing.

## Liver LCM (Brunner 2019, 34 donors) — PROBED 2026-07-17, USABLE

`bash cluster/probe_donors.sh -f cluster/liver.donors.txt` (`~/probe_liver.txt`). The colony
WGS is **WGS_ILLUMINA_short / 151bp / GRCh38** for all 34 donors — my sub-151bp worry was
wrong for the WGS. GRCh38 project per donor: 36713-36718,37105,37111,37113,37114,37237 = **2306**;
37107,37110,37115,37116,37118 = **2325**; 37230-37245 = **2324**; 37904-37918 = **2365**;
48367/48372 (multi-region explants) = **2839/2840/2842**. Older hs37d5/hg19 dup releases also
exist (3415/3426/2462/2589) — assembly is per-BAM, re-check @SQ of any BAM you use.

**The donor glob swept in three OTHER assay arms named like the WGS — EXCLUDE them:**
- **proj 1680 = `TARGETED_ILLUMINA_short` / 75bp** (the matched-normal `PDxxxxxa` panel). The
  ONLY 75bp risk here, and it is not the colony WGS.
- **proj 2454 / 3335 = `RNA-Seq_ILLUMINA_short` / 75-90bp** (transcriptome).
- **proj 2847 / 2927 = `WXS_ILLUMINA_short`** (exome, 151bp but exonic-only).
This is the report-not-gate design working as intended: `probe_donors.sh` shows DS per sample
so you read off which release is the WGS. Never fofn straight off the glob.

The only outstanding item now is **Wilms** (mostly tumour) — not yet probed, low priority.

## Excluded from the inventory (verified, not guessed)

- **THYROID** (PD63118/63121/63126, PD66707-66884, PD43850/51, PD56386): **69 targeted + 17
  exome, ZERO WGS**. An uncovered organ with nothing to call.
- **NF1 `PD61044-61053`**: targeted duplex on muscle. No WGS.
- **Lung PMC7617789**: its "PD ids" are **base64 artefacts** from embedded images, not donors.
- **`PD37590` is ALREADY OURS** (55F colorectal, Lee-Six 2019) — the Wilms paper uses it as an
  adult comparator and mistypes it `PD37950`. Do not re-ingest it as a kidney donor.
- **`PD45518` (stomach) is NOT `PD45517`** (=CB002, cord blood).

Related: [[new-organ-candidates-2026-07]], [[never-infer-tissue-from-pd-number]],
[[irods-cgp-zone-limits]].
