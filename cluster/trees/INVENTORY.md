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
| **Stomach** | 30 | 238 | dual: GRCh37 + GRCh38 | see below | 30 fetched, validated | **assay/read-len UNKNOWN** |
| **Wilms kidney** | 7 | ~40 bulk + ~80 LCM | hs37d5_GRCh37 | 2939/2628/2541/1141 | newick in supplement | mostly TUMOUR |

## NF1 — project 2571 is the WGS release (PROVEN, not inferred)

Exact manifest match against the farm, which is why we trust it:

| donor | manifest WGS | proj 2571 |
|---|---|---|
| PD51122 | 413 | **413** |
| PD51123 | 191 | **191** |
| PD50297 | 305 | 306 (+1) |

**Pin to 2571; ignore every other project for this cohort:**
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

## Stomach — dual-released, and the panel trap

Same samples released twice, identical counts per pair:
`1911→2405`, `1933→2475`, `2134→2476`, `2398→2477` (left = hs37d5_GRCh37, right = GRCh38).
The paper is GRCh38.

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

## THE OPEN QUESTION — nothing should be progressed until this is answered

`bash cluster/probe_new_organs.sh cluster/manifests/oliver2025nf1.wgs.tsv`

Presence in `nst_links` says **nothing** about usability. Still unknown for every donor above:
- **assay type** (`@RG DS:`) — must be WGS.
- **read length** — 151bp is fine; **75bp badly compromises clip-based MEI discovery and
  remapping cannot fix it**. See [[hg19-remap-to-grch38-plan]]: a donor can be on the farm,
  correctly assembled, and still worthless because the reads are short.

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
