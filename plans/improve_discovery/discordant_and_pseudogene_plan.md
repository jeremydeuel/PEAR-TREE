# Plan: discordant-mate anchoring + processed-pseudogene annotation

Two new discovery capabilities, both built the same way every prior toggle was:
**off by default → default output byte-identical** (proven by `differential_all.sh`),
ON-state validated against the VAL-1 synthetic harness, with real-WGS / orthogonal
confirmation flagged as the pre-enable gate.

- **Feature A — `discordant_anchor`**: use discordant read pairs as an evidence
  source, so a one-sided junction (esp. a lone poly-A clip) can be emitted when a
  cluster of discordant mates supplies the missing reciprocal side; optionally label
  whether those mates' destinations fall in an RTE element ("RTE-origin").
- **Feature B — `splice_hallmark`**: a non-gating annotation (like SENS-5) that flags
  candidates whose mate reads map across an exon–exon junction of a single reference
  gene — the processed-pseudogene signature.

Neither exists today. This plan is the build order; per-phase commit+push as before.

---

## Background: what the code does now (grounding)

- Clipped breakpoints are harvested from **high-MAPQ** reads in `handle_record`
  (`discovery.rs:399-489`); **poly-A** breakpoints only from **low-MAPQ discordant**
  reads (`discovery.rs:363-375`, gated on `PolyABreakpoint::find_polya` →
  `!is_proper_pair && mate_is_mapped`, `polya.rs:132-136`).
- `set_mate` (`polya.rs:112-127`) anchors a poly-A breakpoint at the **mate's**
  mapping coordinate — but is **blind to what the mate maps to**.
- `output()` pairs left+right within the TSD window; a poly-A can fill **one** slot
  but the other slot **must** be a real clipped breakpoint (`discovery.rs:875`,
  `while il < l.len() && ir < r.len()`). **A one-sided junction is dropped.**
- Mate **sequences** are attached to breakpoints (`bp.mate_seqs`) and emitted as
  `:MATE` reads; mate **coordinates** are recorded only under `mate_fetch`
  (`mate_coords`, `discovery.rs:451-458`).

The three gaps this plan closes: (1) no discordant-pair evidence class; (2) no
one-sided emission; (3) mate destination never used to establish RTE origin or to
detect a spliced (processed) source.

---

## Feature A — `discordant_anchor`

### A0. What counts as a "discordant anchor read"
A read collected as discordant-anchor support when `config.discordant_anchor` is on:
- primary, mapped, `mapq >= min_mapq`, `contig_ok_extract` passes (a **confidently
  placed** reference anchor — not the element read);
- `!is_proper_pair`; and its mate is either on a **different contig**, or on the same
  contig with `|template_length| > discordant_max_tlen` (default **1000**), or in a
  non-FR orientation.

The **side + approximate breakpoint** from an anchor read: forward read ⇒ insertion is
downstream ⇒ `side=RIGHT`, `pos=reference_end`; reverse read ⇒ `side=LEFT`,
`pos=reference_start` (classic MEI innermost-coordinate logic). Record
`DiscordantObs { contig_id, side, pos, mate_ref, mate_pos }`.

This is a **new read class** entering discovery (these reads are not necessarily
soft-clipped), collected in a `Vec<DiscordantObs>` on `Discovery`, only when enabled.

### A1. Discordant-partner rescue (the core behaviour change)
Cluster `DiscordantObs` by `(contig, side, pos within discordant_window)` (reuse
`cluster_window`), keep clusters with `>= discordant_min_reads` distinct pairs
(default **3**). In `output()`, when a real split/poly-A breakpoint on one side has
**no** reciprocal partner within `[tsd_min, tsd_max]`, allow a discordant cluster on
the opposite side within the window to serve as the partner. This directly answers
"we only have a poly-A clip at one end": the discordant mates supply the other end.

**Pipeline-contract flag (must decide):** the emitted call name needs a token for the
discordant end. Mirror the existing `polyA_<pos>` with **`disc_<pos>`**
(`print_output`, `discovery.rs:981`). Downstream `combine_insertions` (Python) must
tolerate this token — a small matching change there, or confirm it already ignores
unknown coordinate prefixes. **Feature A touches the 4-step pipeline contract; Feature
B does not.**

Discordant reads are emitted as `:LEFT/RIGHT:DISCORDANT<i>` FASTQ entries (new tag)
so the downstream realignment can use them.

### A2. RTE-origin label (and optional gate)
Reuse `IntervalIndex::from_repeatmasker` (`intervals.rs`, family-agnostic +
divergence-gated — the same mechanism as SPEC-7, but applied to the **mate
destination** instead of the breakpoint). For each discordant cluster, compute the
fraction of mates whose `(mate_ref, mate_pos)` lands in the RTE track.
- **Default: annotation only** — add an `rte_origin` fraction column to the sidecar;
  never gates.
- **Opt-in gate** `discordant_rte_only`: require `rte_origin >= discordant_rte_min`
  (default **0.5**) before a discordant cluster may act as a partner — kills random
  structural-variant discordant clusters, keeps genuine element-origin ones.

### A-config keys (all OFF/neutral by default)
| key | default | effect |
|---|---|---|
| `discordant_anchor` | false | master switch for A |
| `discordant_max_tlen` | 1000 | same-contig isize beyond which a pair is discordant |
| `discordant_min_reads` | 3 | min distinct pairs per discordant cluster |
| `discordant_window` | (=`cluster_window`, 6) | cluster width |
| `discordant_rte_track` | none | RTE `.out[.gz]` for mate-origin check (may equal `rm_track`) |
| `discordant_rte_divergence_max` | 20.0 | keep RTE copies ≤ this % divergence for the origin test |
| `discordant_rte_only` | false | gate: require RTE-origin ≥ `discordant_rte_min` |
| `discordant_rte_min` | 0.5 | origin fraction for the gate |

### A-stats
Extend OBS-1 with discordant counters: `disc_obs`, `disc_clusters`,
`disc_paired`, `disc_rejected_rte`.

---

## Feature B — `splice_hallmark` (processed-pseudogene annotation)

### B1. Mate-destination capture
Splice detection needs mate **coordinates** at output time, which today exist only
under `mate_fetch`. Add a per-breakpoint `mate_dests: Vec<(usize /*ref_id*/, i64 /*pos*/)>`
populated in **both** `find_mates_scan` (`discovery.rs:547`) and `find_mates_fetch`
(we already visit the mate record there). Cheap; gated on `splice_hallmark ||
discordant_anchor` so the default path is untouched. Validate default byte-identical.

### B2. Exon-annotation index
New gene-aware index (extend `intervals.rs` or a small `exons.rs`): load a BED/GTF of
exons carrying `gene_id`, build `gene_id -> sorted exon intervals` + a position→(gene,
exon-rank) lookup. The user supplies the annotation for the BAM's assembly (external
input, like the RM track) — flagged as a required file, not shipped.

### B3. Splice/pseudogene signal
For each emitted candidate, gather its mates' destinations (`bp.mate_dests` for both
sides, plus poly-A mate). Map each to `(gene, exon-rank)`. Flag
`processed_pseudogene_candidate = true` when the cluster's mates hit **≥2 distinct
exons of the same gene** whose genomic span **exceeds the summed exon lengths** (an
intron was skipped — the defining processed-copy signature). Emit non-gating:
`<out>.splice.tsv` with `contig, left, right, gene, n_exons, intron_bp_skipped, score`
(or extra columns on `<out>.hallmarks.tsv`). **Main breakpoint output is unchanged.**

Fallback note (no annotation available): a weaker annotation-free heuristic —
mates clustering into ≥2 well-separated same-contig loci — can flag "multi-locus
source" but cannot distinguish a processed gene from an SV; documented, not default.

### B-config keys
| key | default | effect |
|---|---|---|
| `splice_hallmark` | false | write `<out>.splice.tsv`; non-gating |
| `exon_annotation` | none | BED/GTF of exons w/ gene_id (required when on) |
| `splice_min_exons` | 2 | distinct same-gene exons to flag a candidate |

---

## Validation (VAL-1 extensions)

Extend `test/val1/simulate.py` with two labelled classes + matched artefacts:
- `--n-discordant N`: one-sided junctions (real clip/poly-A on one flank only) with a
  cluster of discordant mates pointing to the element on the other flank → tests A1
  recall. Plus `--n-disc-artefact`: one-sided junctions whose discordant mates point to
  **random genome** → must stay rejected when `discordant_rte_only` is on (A2 precision).
- `--n-pseudogene N`: insertions whose mate reads are drawn from **≥2 exons of a
  synthetic multi-exon gene** with the intron removed → tests B3 flags them; and a
  single-exon insertion that must **not** be flagged.

`score.py`/`run_matrix.py` already report recall/precision by class & VAF; add the new
classes. Per-feature exit criterion: recall gain on the target class at preserved
precision, default matrix cell still byte-identical.

**Not-yet-done gates (same as the rest of the register):** real-WGS differential on the
cluster before enabling A on production; orthogonal (long-read/PCR/IGV) confirmation of
new discordant-anchored calls via `new_calls.py` before defaulting any A gate **on**;
RTE track + exon annotation must match the BAM assembly.

---

## Phasing (commit + push per phase)

| Phase | Content | Default-identical? |
|---|---|---|
| **D0** | config keys (all off), VAL-1 simulator + scorer classes | yes (no code path touched) |
| **D1** | discordant-observation collection in scan + stats counters | yes (gated, no output change) |
| **D2** | A1 discordant-partner rescue in `output()` + `disc_`/`:DISCORDANT` tokens | ON-state only; OFF identical |
| **D3** | A2 mate-origin RTE label + `discordant_rte_only` gate | ON-state only; OFF identical |
| **D4** | B1 mate-destination capture (both mate paths) | yes (gated) |
| **D5** | B2 exon index + B3 splice sidecar | ON-state only; OFF identical (sidecar-only) |
| **D6** | README config table, val1 README, decision-register entries | docs |

---

## Decisions to confirm before D2 (flagged inconsistencies / external inputs)

1. **Pipeline contract (A1).** OK to add a `disc_<pos>` coordinate token + a
   `:DISCORDANT` FASTQ tag, and make the matching `combine_insertions` change so the
   Python steps accept them? (Feature A is the first discovery change that isn't purely
   internal to step 1.)
2. **Discordance definition (A0).** Accept the default "different contig **or**
   `|isize| > 1000` **or** non-FR"? Adjust the 1000 bp threshold?
3. **RTE-origin (A2).** Reuse the SPEC-7 RepeatMasker loader (divergence-gated,
   family-agnostic) pointed at an RTE track, annotation-by-default with an opt-in gate?
4. **Exon annotation (B2).** You supply a GENCODE/RefSeq exon GTF/BED (with gene_id)
   for the target assembly — confirm, and name the file/format you want to feed.
5. **Poly-A harvest scope.** Today poly-A breakpoints come only from **low-MAPQ** reads
   (`discovery.rs:363`). Discordant-anchor reads are **high-MAPQ** anchors, so no
   conflict — but confirm you don't also want high-MAPQ poly-A reads harvested (a
   separate, larger change I'd scope only if you want it).
