# TPRT-hallmark overhaul — shared spec (branch `tprt-hallmarks`, base `PEAR-TREE2`)

This file is the **contract between parallel work packages**. If you need to deviate
from a format below, change this file in the same commit and say so in your report.

## Decisions already made by Jeremy (do not re-litigate)

- Base is `PEAR-TREE2`. Human only (mouse later). No UMIs, no hybcap.
- Discovery BAMs may be **GRCh38, hg19 or hs1**; we never lift BAMs. All *remapping*
  (combine/annotate bowtie2/minimap2) is against **hs1**. Libraries therefore need to be
  sequence-based, with coordinates recorded in hs1 (primary) and hg38 (documentation).
- New tools allowed: `minimap2`, `edlib`, `parasail` (and their Python bindings
  `mappy`, `edlib`, `parasail`). Keep bowtie2/HMMER where they already work.
- ">= 2 independent fragments per breakpoint" is enforced **after combine, pooled across all
  colonies of a patient**, and applies to **every junction including the poly-A/3' end**.
  Discovery must therefore pass single-fragment per-sample evidence through (configurable).
- Independence is decided from the reads, **never from the 0x400 dup flag**.
- Always fetch the mate of every evidence read.
- Consensus building must be **indel-aware**, especially in homopolymers / poly-A
  (Illumina SBS gives poly-A length jitter; the 3'-beyond-poly-A sequence must survive).
- Translocation bridges / chimeric bridges (Zumalave 2026) are **out of scope** this round.
- Old analysis/ artefact data are gone; calibrate on the simulators + test/fullstack.

## Insertion-type vocabulary (used by docs, simulator truth, and annotate output)

`element`: `L1`, `ALU`, `SVA`, `PSEUDOGENE`, `POLYA_ONLY` (solitary poly(A/T)),
`ORPHAN_TD` (3' transduction with no element sequence), `NON_TPRT`, `UNKNOWN`.

`structure` (5' end of the element, needs 5' junction coverage, else `5P_UNRESOLVED`):
`FULL_LENGTH`, `TRUNCATED_5P`, `INVERTED_5P` (twin priming), `INVERTED_5P_SWITCH`
(twin priming + additional internal inversion/switch), `5P_UNRESOLVED`.

`tags` (comma list, any combination): `TD3P` (partnered 3' transduction), `TD5P` (SVA 5'
transduction), `TD3P_SOURCE=<source_id>`, `NOVEL_SOURCE`, `TEMPLATED_LOCAL` (short non-element
template from <=15 bp of the site embedded), `PREMRNA_COINSERT`, `TSD_DELETION` (target-site
deletion instead of duplication), `EN_INDEPENDENT` (no TSD, no poly-A, both ends truncated),
`L1_MED_DELETION`, `L1_MED_DUPLICATION`, `FOLDBACK_INVDUP_5P`, `CHIMERIC_ENDS` (5'/3' ends
from incompatible elements — artefact indicator), `EXON_JUNCTION` (pseudogene proof).

Literature catalogue (16 types) lives in `docs/insertion_types.html`; out-of-scope types are
shown there but marked.

## Data flow and new files

```
discovery (Rust, per sample)
   <sample>.txt.gz                (UNCHANGED format)
   <sample>.txt.gz.evidence.tsv.gz  NEW sidecar, one row per evidence read (incl. mates)
combine_insertions (Python, per patient)
   <patient>.combined.txt.gz      (UNCHANGED format; consensus now indel-aware)
   <patient>.genotyping.txt.gz    (UNCHANGED)
   <patient>.insertions.evidence.tsv.gz  NEW: per insertion junction pooled stats + consensus
   <patient>.insertions.reads.fa.gz      NEW: pooled reads+mates per insertion (for annotate)
annotate_v2 (+ new tools/rte/ package)
   existing table + new columns (see below)
```

### `<sample>.txt.gz.evidence.tsv.gz` (discovery → combine)

Header line, tab-separated, one row per read record:

| col | meaning |
|---|---|
| `locus` | the discovery locus id exactly as in the `.txt.gz` record name (`contig:L-R`) |
| `side` | `LEFT` / `RIGHT` |
| `role` | `CLIP` (junction-clipped read), `POLYA` (poly-A read placed by mate), `MATE` (mate of a CLIP/POLYA/DISC read), `DISC` (discordant anchor), `SPAN` |
| `frag` | 64-bit hash of qname (hex) — links a read to its mate; same fragment ⇒ same value |
| `r12` | 1 or 2 |
| `flag` | raw SAM flag (dup bit kept for diagnostics only) |
| `ref` / `pos` / `strand` | alignment of this record (0-based pos; `*`/-1 if unmapped) |
| `outer` | unclipped 5' end of this read on the reference (soft clips included) |
| `mref` / `mpos` / `mstrand` | mate alignment from the record (`*`/-1 if unmapped) |
| `tlen` | SAM TLEN |
| `mapq` | MAPQ |
| `cigar` | CIGAR |
| `clip_at` | read offset of the junction (CLIP only, else -1) |
| `seq` / `qual` | full read sequence in **reference-forward orientation** + phred+33 |

Caps: `max_mates_per_breakpoint` (default 50) and `max_evidence_reads_per_breakpoint`
(default 200) per side; when capped, keep a deterministic subset.

**As implemented (discovery worker, Rust) — read this before parsing the sidecar:**

- **Path**: `<out>.evidence.tsv.gz` next to the other sidecars, i.e.
  `<sample>.txt.gz.evidence.tsv.gz` (the wrappers rename `$TMP.<ext>` → `$OUT.<ext>`; a
  `<sample>.evidence.tsv.gz` name could not survive that rename). gzip, one header line.
- `frag`: FNV-1a 64-bit hash of the qname, 16 lowercase hex digits.
- `r12`: 1/2 from flag 0x40/0x80 (0 for unpaired input).
- `pos`/`outer`/`mpos`: 0-based. `outer` = 0-based *inclusive* coordinate of the read's 5'
  base incl. soft clips (hard clips skipped): `pos - leading_softclip` on `+`,
  `aln_end - 1 + trailing_softclip` on `-`. Unmapped record: `ref *`, `pos -1`, `strand *`,
  `outer -1` (its `seq`/`qual` are still written). Mate unmapped: `mref *`, `mpos -1`, `mstrand *`.
- `seq`/`qual`: exactly as stored in the BAM (= reference-forward for mapped records); `*` if absent.
- `clip_at`: offset of the junction in the stored `seq` (LEFT clip: leading soft-clip length;
  RIGHT clip: `len - trailing soft-clip length`); `-1` for non-CLIP rows.
- Rows per locus, LEFT then RIGHT; per breakpoint side `CLIP` rows, then `DISC`, then `MATE`.
  A poly-A end writes `POLYA` + its anchoring `MATE` for the paired poly-A read **and** every
  other poly-A read of the same clip side within `cluster_window` (pooled poly-A end).
- `CLIP` = every read of the clustered junction (incl. supplementary records and reads off the
  modal position), not only those that passed the consensus floor.
- `DISC` = primary, MAPQ ≥ `min_mapq`, non-proper pair whose mate is unmapped / on another
  contig / > `discordant_max_tlen` away / same strand, pointing into the insertion: a REVERSE
  anchor starting ≥ B−5 and ending ≤ B+`sidecar_disc_span` for a LEFT junction at B, a FORWARD
  anchor ending ≤ B+5 and starting ≥ B−`sidecar_disc_span` for a RIGHT junction. Independent of
  `discordant_anchor`. Fragments already present as CLIP are not repeated as DISC.
- `MATE` = the primary mate record (not supplementary). Without `fetch_all_mates`: mates of the
  legacy `has_mate` orientation only (LEFT reverse / RIGHT forward CLIP reads). With it: mates of
  every CLIP and DISC read. A mate already present as a primary CLIP record is not repeated.
- `SPAN` is reserved but **not emitted yet** (needs the D2 fragment-spanning-pair work, not on
  this base).
- Caps per side: CLIP+DISC ≤ `max_evidence_reads_per_breakpoint` (CLIP first), MATE ≤
  `max_mates_per_breakpoint`; subset = lowest (frag, flag) / (frag, r12) — deterministic and
  independent of `--threads`.
- Extra key: `sidecar_disc_span` (500).
- Locus ids are exactly the `.txt.gz` record-name prefixes; the set of loci in the sidecar
  equals the set in the `.txt.gz` (sub-floor partners of the dormant Feature-A rescue carry
  CLIP rows but no MATE/DISC).
- `ignore_dup_flag = false` keeps the legacy asymmetry: the clip path and mate pass drop 0x400,
  the low-MAPQ poly-A path never checked it (changing that would break byte-identity). With
  `true`, dups are kept everywhere (clip path, poly-A, mate pass, sidecar).
- `min_evidence_fragments_per_sample`: distinct qname hashes among the reads supporting the
  modal position (± `evidence_window`); `adaptive_evidence` scales it like the read floor
  (`max(base, round(base·local/median))`). With floor ≤ 1 a single-read cluster takes the normal
  consensus path (not only the poly-A rescue). Also fixes the LEFT `n_reads` double count (only
  in this mode; the legacy value only feeds the off-by-default Feature-A rescue).

New discovery config keys (default values keep the current FASTQ output byte-identical):
`evidence_sidecar` (false), `fetch_all_mates` (false), `ignore_dup_flag` (false),
`min_evidence_fragments_per_sample` (unset ⇒ current read floor). `cluster/config.discovery.*`
for the new pipeline mode sets: sidecar=true, fetch_all_mates=true, ignore_dup_flag=true,
min_evidence_fragments_per_sample=1. The poly-A single-read rescue stays in discovery (it is
the pooled rule in combine that enforces >=2).

### Independence rule (combine)

Per insertion junction, pooled over colonies:
1. Collapse records by (`sample`, `frag`) → one fragment (mates and supplementary of one
   template count once).
2. Within one sample, two fragments are PCR/optical duplicates if their outer coordinates
   match on both ends (read `outer` ±2 bp and mate `outer`/`mpos` ±2 bp, same strands) **or**
   (same read `outer` ±2 and junction-clip sequences with Hamming/edit distance ≤ 2 over the
   overlap **and** mate placement within ±2). Duplicates collapse to one.
3. Fragments in different samples are independent libraries ⇒ independent, except exact
   identity of both outer coordinates AND sequence (flag as `cross_sample_identical`).
4. Junction supported iff `n_independent >= min_independent_fragments` (default 2).
   Applies to LEFT, RIGHT and the poly-A end. Config key under
   `CONFIG['combine_insertions']`.

### `<patient>.insertions.evidence.tsv.gz` (combine → annotate)

One row per insertion junction: `insertion_id` (as in combined.txt.gz), `side`,
`n_reads`, `n_fragments`, `n_independent`, `n_samples`, `n_mates`, `supported` (0/1),
`clip_consensus` (indel-aware, lowercase clip / uppercase reference, same convention as
combined.txt.gz), `consensus_depth` (comma list per clip base), `polya_len_median`,
`polya_len_range`, `beyond_polya` (sequence 3' of the poly-A recovered by consensus, may be
empty), `beyond_polya_support` (independent fragments covering it).

Implementation notes (combine worker; `src/combine_insertions_evidence.py`, `src/indel_consensus.py`):
- Orientation: `clip_consensus`, `consensus_depth` (one value per lowercase clip base, same
  order) and `beyond_polya` are all **reference-forward**. "Beyond" = 3' of the poly-A in
  element sense: for a poly-A seen as an A-run going outward from the junction it is the
  sequence further out; for a T-run it is the sequence between junction and run (empty when
  the T-run starts at the junction). `polya_len_range` = `min-max` of trusted per-read runs.
- `supported` is `NA` when an input file contributing to the insertion has no sidecar
  (such insertions are never gated).
- `clip_consensus` here uses junction reads **and** overlapping mates (may extend through a
  short insertion into the far flank); the clip written to `combined.txt.gz`
  (`indel_aware_consensus=true`) uses junction reads only, so the clipped-remap filter is not
  fooled by flank sequence.
- Additive trailing columns (after the SPEC ones, order fixed): `polya_end` (0/1),
  `n_duplicates` (within-sample merges), `n_cross_sample_identical`, `fail_reason`,
  `consensus_stop` (`end`/`depth`/`disagreement`/`empty`).
- Rows: surviving insertions first (combined.txt.gz order), then gated-out ones (supported=0).
- Config (`CONFIG['combine_insertions']`): `require_independent_fragments` (False; True in
  `cluster/config.py.grch38.tprt`), `min_independent_fragments` (2), `indel_aware_consensus`
  (False; True in .tprt), `dup_coord_tolerance` (2), `dup_max_edit` (2), `polya_min_len` (8).

### `<patient>.insertions.reads.fa.gz` (combine → annotate)

FASTA, record name `insertion_id|side|role|sample|frag|r12`, sequence in reference-forward
orientation. Mates included. annotate builds the covered-element consensus from these.

### Reference libraries — `resources/rte_library/` (committed; small) built by `tools/rte_library/build.py`

| file | content |
|---|---|
| `l1_intact.fa` / `.tsv` | L1Base hsflil1_8438 (146, GRCh38), oriented sense, with subfamily (rmsk), Ta/pre-Ta, hs1+hg38 coords |
| `alu_y_intact.fa` / `.tsv` | young near-full-length AluY* elements (rmsk filtered), sense |
| `sva_intact.fa` / `.tsv` | near-full-length SVA_A..F, sense |
| `consensus.fa` | per-class consensus built from the intact sets (`L1HS`, `L1PA2`, `ALU_Y`, `ALU_YA5`, `ALU_YB8`, `SVA_E`, `SVA_F`, …), Dfam consensus as cross-check |
| `active.tsv` | subset regarded as active/hot (L1HS-Ta, known hot sources) + per-element identity to class consensus |
| `transduction_sources.tsv` | source elements: id, class, hs1+hg38 coords, strand, reference/non-reference, evidence (paper), hotness |
| `flanks_3p.fa` | 0–15 kb downstream of each source (sense of the element), repeats soft-masked |
| `flanks_5p_sva.fa` | upstream flanks of SVA sources (SVA 5' transductions) |
| `README.md` | provenance, licences, rebuild command |

Documentation: `docs/transduction_sources.html` — includes how **novel sources** are accepted.

### Novel source rule (annotate)

A unique (non-repeat, MAPQ≥20 on hs1) inserted segment not in `flanks_3p.fa` is a *credible
novel 3' transduction source* if, on hs1, it lies within 15 kb **downstream** (strand-aware)
of a reference L1 that is ≥5.5 kb and ≥ 95 % identical to the L1HS consensus (or of an L1
insertion called elsewhere in the same cohort), and the insertion carries TPRT hallmarks
(poly-A after the tag, TSD/EN). Report `TD3P_SOURCE=novel:<hs1 coords>` + `NOVEL_SOURCE` and the
source element's identity to consensus.

### annotate output — new columns

`element`, `structure`, `tags`, `covered_5p`, `covered_3p` (consensus coords of covered part),
`element_identity` (to nearest active element), `nearest_active`, `tsd_seq`, `tsd_len`,
`en_motif` (the 7-mer at the nick, strand-corrected), `en_mismatches`, `polya_len`,
`beyond_polya`, `tprt_score`, `tprt_points` (semicolon list `feature:+n`), `tprt_call`
(`TPRT` / `LIKELY_TPRT` / `UNCERTAIN` / `ARTEFACT_LIKE`).

### TPRT point system (annotate, `tools/rte/score.py`) — design targets

Additive, every feature reported. Positive: TSD 4–25 bp with sequence identity (strong) or
TSD_DELETION ≤ 20 bp (weak); poly-A ≥ 10 at the 3' junction on the strand-consistent side;
**covered 3'-beyond-poly-A sequence (consensus extends past the poly-A into flank/transduction
— extra points, more if ≥2 independent fragments)**; EN motif at the nick (TTTT/AA-like, by
mismatch bin 0–1 / 2 / 3); 5' and 3' junction element classes concordant (or the 3' tag is
the transduced flank of exactly the element at the 5' end); high identity to an active element;
twin-priming junction ≥ 590 bp into L1; exon–exon junction for pseudogenes; independent
fragment count and number of colonies. Negative (artefact): poly-A at both junctions / slippage
context; chimeric ends; TSD > 50 bp; old/inactive subfamily only; clip identical to the
adjacent reference (fold-back); recurrence across many unrelated loci; single-sample with
cross-sample-identical fragments. Calibrate thresholds on simulator truth (TP vs artefact).
