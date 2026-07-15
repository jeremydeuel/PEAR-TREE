# 7. PEAR-TREE: architecture and code review

*Source reviewed: `~/Documents/PEAR-TREE` (branch `PEAR-TREE2`), commit state of 2026‑07‑15. This section reviews the code as a detection method and maps each of its filters onto the artefact classes discussed in [§4](04_artefacts_in_discovery.md) and [§5](05_sequencing_artefacts.md).*

---

## 7.1 What PEAR-TREE is

PEAR-TREE ("paired ends of aberrant retrotransposons in phylogenetic trees") detects
**insertional mutations** — chiefly retrotransposon (RTE) insertions — from short-read
WGS BAM files, and genotypes them across many samples so that insertions can be placed
on a phylogenetic tree of clonal samples (e.g. colonies/organoids/crypts). It is the
whole-genome sibling of the targeted-sequencing tool PARTRIDGE.

The design philosophy is distinctive and worth stating up front because it shapes every
filter: **PEAR-TREE is a split-read / clipped-read caller, not a discordant-read-pair
caller.** The primary evidence for an insertion is a pair of soft-clipped reads facing
each other, whose clipped portions represent the two ends of inserted (non-reference)
sequence. Paired-end mates are used only to *extend* the clipped consensus, not as the
primary discovery signal. This is the opposite emphasis to TraFiC/MELT/Delly, which lead
with discordant pairs and use split reads for refinement (see [§3](03_detecting_true_events.md),
[§6](06_delly_review.md)).

## 7.2 The four-step pipeline

| Step | Entry point | Genome needed? | Purpose |
|------|-------------|----------------|---------|
| 1. **Discovery** | `discovery.py` | No | Per-BAM scan for clipped-read breakpoint pairs; emits a custom `.txt.gz` |
| 2. **Combine insertions** | `combine_insertions.py` | Yes (2bit + bowtie2 + chain) | Merge candidates across BAMs; genome-aware artefact removal; build genotyping panel |
| 3. **Genotype** | `genotype.py` | BAM only | Re-genotype every sample at every candidate locus (ref / alt / artefact) |
| 4. **Combine genotypes** | `combine_genotypes.py` | No | Assemble genotype matrix; population-level filtering |

The on-disk contract between steps is a custom FASTQ-like text format (`>locus`,
`@RIGHT_CONSENSUS`, sequence, `+`, quality). The PEAR-TREE2 plan (`PEAR-TREE2_PLAN.md`)
explicitly preserves this contract while it speeds discovery up and moves more artefact
removal earlier.

## 7.3 Discovery — the core algorithm (`discovery.py`, `breakpoint.py`)

### 7.3.1 Breakpoint model

A breakpoint is a soft-clip position on one read. Directions are relative to the `+`
strand (BAM convention). Using the module's own diagram:

```
                R                                  L
      ==========*****************....**************============
        RIGHT_CLIPPED                    LEFT_CLIPPED
```

- **R** = reference position of the first clipped base of a *right*-clipped read (the 5′
  end of the insertion on the `+` strand).
- **L** = reference position of the first *mapped* base of a *left*-clipped read (the 3′
  end).
- An insertion site is named `seqname:R-L`. A **target-site duplication (TSD)** makes
  `L < R` (their difference is `−TSD_length`); a target-site *deletion* makes `L > R`.

This TSD-aware pairing is the key biological signal: it is exactly the hallmark that
distinguishes a true TPRT-mediated (target-primed reverse transcription) insertion from
most artefacts ([§2](02_biology_of_retrotransposition.md)).

### 7.3.2 Discovery walk (`extract_chimeric`)

Iterating a **coordinate-sorted** BAM (`_assert_coordinate_sorted` fails loudly
otherwise):

1. **MAPQ gate.** `read.mapping_quality < min_mapq (40)` → the read is not usable as an
   anchor. Such low-MAPQ reads are *not* simply discarded; they are checked for a poly‑A
   tail (`PolyABreakpoint.findPolyA`) and retained as poly‑A evidence if found. This is a
   good design: the anchoring mate of a poly‑A read is often in a repeat and low-MAPQ.
2. **Flag gates.** Secondary, QC-fail and duplicate reads skipped.
3. **Contig gates.** `len(reference_name) > 5` skips alt/random/decoy contigs; `MT`/`chrM`
   skipped. (This is a crude but effective way to drop the artefact-dense contigs — see
   §7.7 caveats.)
4. **Clip classification.** From `cigartuples`, decide `CLIP_LEFT` vs `CLIP_RIGHT`; if
   both ends clipped, the longer clip wins.
5. **Cruciform / micro-indel exclusion.** If the read is *supplementary* and its `SA` tag
   shows a split alignment to the **same contig within `exclude_same_contig_supplementary`
   (1000 bp)**, an `exclude` flag is set. In `Breakpoint.join` a single excluded read
   "poisons" the whole cluster (returns `None`). This is PEAR-TREE's dedicated defence
   against the palindrome/inverted-repeat/foldback artefacts of [§5](05_sequencing_artefacts.md):
   a clipped part that re-maps locally is not an insertion.
6. **Adapter clip.** `is_adapter()` on the clipped end drops read-through adapter
   (NebNext defaults) and simple repeat dimers.
7. Surviving clips are added via `add_breakpoint` as `Breakpoint` objects holding the
   clipped and unclipped `QualitySeq`.

### 7.3.3 Clustering and consensus (`cleanup` + `Breakpoint.join`)

Per contig, left and right breakpoints are sorted and greedily clustered when successive
positions are within **6 bp** (note: this hard-coded 6 is *not* the configurable
`max_bp_window` = 40; see §7.7). Each cluster goes through `Breakpoint.join`, which is the
real artefact gauntlet:

- **Single-read rescue.** A lone read is discarded *unless* its clipped tail is a clean
  poly‑A (`AAAAAAAA`) / poly‑T (`TTTTTTTT`), in which case it is rescued as poly‑A
  evidence. Otherwise `too_few`.
- **Exclusion poisoning.** Any `exclude`-flagged read in the cluster kills it (cruciform).
- **Evidence threshold.** The modal breakpoint position must be supported by
  `min_evidence_reads_per_breakpoint` (2) precise reads with clipped length ≥
  `min_good_bases` (10); otherwise poly‑A rescue is attempted, else `too_few_after_filter`.
- **Consensus.** Clipped and unclipped consensuses are built (`find_consensus`, extends
  only to the first ambiguous base). Clipped consensus must exceed `min_clip_len` (12);
  unclipped must exceed 40 bp, otherwise `clipped_failed` / `unclipped_failed`.
- **Homopolymer / low-period filter.** The unclipped consensus is rejected if its first
  24 bp are an exact 1–4-bp repeat (`polymer`). This removes breakpoints anchored in
  microsatellite/low-complexity DNA — a major artefact source.

The per-side `Breakpoint.stats` counter (too_few / rescued_pA / excluded /
too_few_after_filter / clipped_failed / unclipped_failed / polymer / passed) is an
excellent, under-exploited QC instrument: it is effectively a discovery-time artefact
audit and should be written to the output for every run.

### 7.3.4 Poly‑A handling (`polyABreakpoint.py`)

Poly‑A is the single most specific sequence signature of L1/Alu/SVA retrotransposition.
PEAR-TREE treats it as a first-class breakpoint type:

- A read whose sequence contains ≥12 A (or T) and that is **not** in a proper pair, with a
  mapped mate, becomes a `PolyABreakpoint`. Read1/2 and strand logic decide whether it is
  a `CLIP_RIGHT` (5′) or `CLIP_LEFT` (3′) event, and the poly‑A run is stripped to recover
  the adjacent unique clipped sequence.
- In `discovery.output`, poly‑A breakpoints "rescue" one-sided clipped breakpoints: a
  right-clip with no matching left-clip is paired with a nearby poly‑A within a 12–120 bp
  window (and vice versa). This directly implements the biological expectation that the 3′
  end of an L1 insertion is a poly‑A tail — the same logic TraFiC/MELT use to confirm
  MEIs.

### 7.3.5 Mate extension (`find_mates`, `extend_mates`)

A second full BAM pass collects mate sequences for reads at breakpoints, orients them by
strand, and `extend_mates` lengthens the clipped consensus by walking mate overlaps. This
turns ~150 bp of clipped evidence into a longer contig for better downstream repeat
annotation. **Note:** the plan flags `extend_mates()` as currently iterating an
already-emptied list (a silent no-op) — see §7.7.

## 7.4 Combine insertions — genome-aware artefact removal (`combine_insertions.py`)

This is where PEAR-TREE separates true insertions from local mapping artefacts, and it is
the conceptual heart of "artefact detection that needs the reference." Steps:

1. **Cross-file intersection** (`intersect_insertions`). Identical loci from many BAMs are
   merged; the merge requires the clipped/aligned consensuses to agree
   (`sequence_matching_score ≥ 0.6`), so a locus only survives if independent samples
   report a *consistent* inserted sequence. Random artefacts that differ base-to-base are
   dropped here.
2. **High-insertion-rate region masking.** Loci are binned into 100 bp windows; any window
   with ≥4 candidate insertions is wiped (± the window). This is a blacklist-by-density —
   the pragmatic substitute for the ENCODE/mappability blacklists used elsewhere ([§4](04_artefacts_in_discovery.md)).
3. **"End maps entirely" filter (bowtie2 `--end-to-end --sensitive`).** The full
   aligned+clipped consensus of each end is realigned to a T2T reference. If *either* end
   maps end-to-end, the locus **cannot be chimeric** (there is no true novel junction) and
   is removed. This catches reference-genome structure artefacts and mis-assembly-driven
   soft clips.
4. **"Clipped part maps locally" filter (bowtie2 `--local --very-fast -k 1000` + liftover).**
   The *clipped* portion alone is realigned; if it maps within **1000 bp** of its own
   breakpoint (after chain liftover between the T2T index and the BAM's assembly), the
   locus is a local rearrangement/artefact, not an insertion, and is removed. This is the
   genome-level version of the cruciform check from discovery.
5. **Reference-similarity filter for genotyping panel.** For each surviving locus the tool
   fetches the reference sequence flanking the breakpoint and requires the *alt* (clipped)
   sequence to be **dissimilar** to the *ref* sequence (`sequence_matching_score > 0` →
   excluded). If ref and alt look the same, there is nothing to genotype and it is likely
   an alignment wobble. It also drops loci whose two flanks are identical or
   reverse-complements (palindromic/duplication artefacts).

## 7.5 Genotyping and population filters

**Per-sample genotyping** (`genotype.py`, `genotyping_insertion.py`, `genotype_qscore.py`).
At each locus, every spanning read (MAPQ ≥ 40, not secondary/dup/QC-fail) is scored
against the **ref** and **alt** junction sequences using a quality-weighted match score,
with a ±1 bp offset search (`qscore` tries three registers) to tolerate breakpoint
imprecision. Each read contributes to `q_ref`, `q_alt` or `q_art`:

- A read matching **alt on both sides** (`double alt`) is scored as **artefact**, not
  insertion — a neat internal consistency check: a genuine heterozygous insertion shows
  ref on one allele and alt on the other, whereas an artefact read matches the inserted
  junction on both flanks.
- The read-level "artefact" evidence (`q_art > 60` and dominant) aggregates to an
  **`artefact` genotype** for the sample. This per-sample artefact call is what step 4 then
  counts.

**Population filtering** (`combine_genotypes.py`). The genotype matrix (one row per
insertion, one column per sample) is filtered on: `min_wild-types`, `min_insertions`
(≥1 confident het/hom), `max_artefact` samples, `max_na` samples, "more than half
uncertain," "more uncertain-insertion than certain-insertion," and a best het/hom score ≥
800. Requiring both **wild-type** and **insertion** samples for the same locus is a
powerful, population-level artefact filter: a recurrent library/reference artefact tends to
appear as *artefact* or *NA* across many samples rather than as a clean het/hom-vs-wt
split, so it fails `min_insertions`/`max_artefact`.

## 7.6 How PEAR-TREE's filters map to artefact classes

| Artefact class ([§5](05_sequencing_artefacts.md)) | PEAR-TREE defence | Stage |
|---|---|---|
| Read-through **adapter** | `is_adapter`, `clean_clipped_seq` | Discovery |
| **Poly‑G** (two-colour dark cycles) | `clean_clipped_seq` trims 3′ poly‑G | Discovery |
| **Homopolymer / microsatellite** slippage | `max_homopolymer_len`; n‑polymer unclipped filter; low-complexity | Discovery |
| **Cruciform / inverted-repeat / foldback** chimeras | `SA` same-contig-within-1000 bp exclusion → cluster poisoning | Discovery |
| **Local rearrangement** mimicking insertion | clipped-maps-locally bowtie2 `--local` + liftover | Combine |
| **Reference structure-specific** sequences / segdups | end-to-end remap removal; high-density region mask | Combine |
| **Random / non-reproducible** chimeras | cross-file consensus agreement (≥0.6); ≥2 evidence reads | Discovery + Combine |
| **High-coverage pileup** artefact regions | 100 bp density mask (combine); `max_read_count`/`reads_for_high_coverage` (partly disabled — §7.7) | Combine/Genotype |
| **Recurrent** library/mapping artefacts across samples | `max_artefact`, `max_na`, require wt+ins split | Combine genotypes |

The correspondence is strong: PEAR-TREE already implements, in some form, most of the
artefact-removal strategies the literature recommends. The main gaps are about *when* and
*how completely* they run.

## 7.7 Findings, gaps and risks

These are observations from reading the code, ordered roughly by impact on
sensitivity/specificity. They are consistent with, and extend, the authors' own
`PEAR-TREE2_PLAN.md`.

1. **High-coverage exclusion is effectively off.** `max_read_count` (config/README) is not
   referenced on the discovery path, and genotyping's high-coverage check is guarded by
   `if False and ...` (`genotype.py:54`). High-coverage pileups are, per the artefact
   literature, one of the richest false-positive sources. Re-enabling discovery-time
   coverage masking (Plan §2.3) is the single highest-ROI specificity improvement.
2. **Two clustering distances disagree.** Discovery clusters clipped reads with a
   **hard-coded 6 bp** window (`discovery.py:136,147`), but `max_bp_window = 40` (config)
   is what the README documents as the pairing window and what the TSD logic in
   `output()` uses (`tsd > 40`). A TSD between 6 and 40 bp can therefore have its two
   sides clustered inconsistently. Worth reconciling and making the 6 configurable.
3. **`extend_mates()` is a silent no-op** (Plan §1.3): it iterates
   `self.temporary_breakpoints`, already emptied by the final `cleanup()`. Clipped
   consensuses are therefore *not* being extended by mates in the current code, despite the
   pipeline name emphasising paired ends. Decide whether to point it at the final
   breakpoint lists (longer contigs → better repeat annotation, possible new artefacts) or
   remove it. Needs a real-WGS benchmark.
4. **Poly‑A-only insertions are dropped in combine.** `intersect_insertions` has an
   unconditional `continue` (`combine_insertions_intersect_insertions.py:84`) that skips
   *all* poly‑A-only loci ("this needs a heavy filter"). So the elaborate poly‑A rescue in
   discovery is partly wasted downstream: an insertion seen only via its poly‑A end in a
   given sample will not be promoted to a full insertion unless the other end is also seen.
   This is a deliberate conservative choice but costs sensitivity for 5′-truncated L1s and
   Alu/SVA, whose non-poly‑A end is short.
5. **No discordant-read-pair discovery.** The discordant-mate branch was removed as dead
   code and deferred. PEAR-TREE therefore relies entirely on reads that *span* the junction
   with a clip. Events where no single read crosses the breakpoint (large TSD, low
   coverage, junction in unmappable sequence) are invisible, whereas a TraFiC/Delly-style
   discordant-pair cluster would still flag them. This is the biggest *sensitivity*
   architectural difference from the reference tools ([§3](03_detecting_true_events.md),
   [§6](06_delly_review.md)).
6. **Contig filter by name length is fragile.** `len(reference_name) > 5` is used as a
   proxy for "weird contig." It correctly drops `chr14_...` alt contigs but also silently
   handles only assemblies whose main chromosomes have ≤5-char names; it would misbehave on
   contig-naming schemes like `NC_000014.9`. An explicit main-chromosome allowlist (or a
   proper `.bed` include-region) is safer and doubles as the artefact blacklist mechanism.
7. **Filters exist but are not wired in.** `has_well_defined_breakpoint` and
   `is_low_complexity` (`sequence_checks.py`) are implemented but not called on the
   discovery path (Plan §2.3). Wiring them in would move low-complexity rejection earlier
   and cheaply.
8. **No explicit endonuclease-motif or TSD-length prior.** PEAR-TREE uses TSD *geometry*
   (L vs R ordering, window ≤40) but does not score the L1 EN nick motif (5′‑TTAAAA) or
   prefer the canonical ~2–20 bp TSD length range, both of which are strong positive
   evidence used in the literature to raise precision without losing many true events.
9. **Hard-coded absolute cluster paths in `config.py`** (Sanger `/lustre/...` and
   `/nfs/...`). Fine for the authors, but it means the shipped default config is not
   runnable elsewhere; the plan's move to `--config file.yaml` (§2.4) is the right fix and
   also makes the artefact/adapter tables reproducible per platform.

## 7.8 Strengths worth preserving

- **TSD-first, clip-first discovery** gives base-resolution breakpoints natively — no
  separate split-read refinement pass needed.
- **Poly‑A as a first-class, rescuing signal** matches the biology and recovers one-sided
  events.
- **Layered, escalating artefact removal** (cheap sequence filters → cross-sample
  consensus → genome remapping → population genotype consistency) is exactly the
  defence-in-depth the artefact literature argues for, and it puts the genome-free filters
  first, which is efficient.
- **The per-side `stats` audit** is a ready-made discovery-time artefact dashboard.
- **Population/phylogenetic genotyping** provides an orthogonal, powerful artefact filter
  that single-sample callers lack: true somatic insertions partition samples cleanly;
  artefacts do not.

See [§8](08_synthesis_and_recommendations.md) for how these findings combine into a
concrete recommendation list for PEAR-TREE2.
