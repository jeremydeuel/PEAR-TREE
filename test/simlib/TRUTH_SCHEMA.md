# Simulator truth schema (insertion-type catalogue)

Shared by `test/val1/simulate.py --types ...` and `test/fullstack/build_donor.py --types ...`.
Labels use the SPEC vocabulary (`plans/tprt_hallmarks/SPEC.md`).

## Coordinate columns

**val1** (`--out-truth`): the 8 legacy columns are kept first, unchanged in meaning, so
`score.py`/`run_matrix.py` and older readers keep working:

| col | meaning |
|---|---|
| `contig` | synthetic contig (`1`/`2`/`3`) |
| `left` | LEFT-clip breakpoint (0-based; reads *start* here, clipped on their left). For a `+` insertion this is the 3'/poly-A junction |
| `right` | RIGHT-clip breakpoint (reads *end* here). For a `+` insertion the 5' junction |
| `class` | `<element>_<variant>` (legacy rows: the old class string) |
| `tsd` | `right - left`: TSD length (>0), **negative** for target-site / L1-mediated deletions, 0 for none |
| `alt_reads` | reads crossing a junction (>=10 bp each side), all samples, incl. PCR duplicates |
| `ref_reads` | 0 for catalogue rows (legacy rows: reference reads) |
| `vaf` | mean VAF over samples |

**fullstack** (`truth_types_hs1.tsv`): `id`, `hs1_contig`, `hs1_left`, `hs1_right`, `tsd`,
`donor_record`, then the label columns. `score_types.py` lifts each event to the discovery
reference by mapping its 300 bp hs1 flanks (`flanks.fa`, `<id>_L` ends at min(left,right),
`<id>_R` starts at max(left,right)) and writes `results_by_event.tsv`.

Matching is order-free on the two breakpoints (deletions have `left > right`).

## Label columns (appended, same order in both)

| col | values |
|---|---|
| `role` | `TP` (real insertion) or `ARTEFACT` (library/mapping artefact; a call here is a false positive, reported separately) |
| `type_id` | literature catalogue id (see below); artefacts `0`; legacy val1 rows `0` |
| `variant` | catalogue key, e.g. `L1_INV`, `ORPHAN_TD3P`, `ART_LIGATION_PCR` (legacy rows `legacy:<class>`) |
| `element` | `L1` `ALU` `SVA` `PSEUDOGENE` `POLYA_ONLY` `ORPHAN_TD` `NON_TPRT` `UNKNOWN` |
| `structure` | `FULL_LENGTH` `TRUNCATED_5P` `INVERTED_5P` `INVERTED_5P_SWITCH` `5P_UNRESOLVED` (no element 5' end: poly-A-only, orphan TD, decoy, artefacts without one) |
| `tags` | comma list: `TD3P`, `TD3P_SOURCE=<id>`, `TD5P`, `TD5P_SOURCE=<id>`, `TEMPLATED_LOCAL`, `PREMRNA_COINSERT`, `TSD_DELETION`, `EN_INDEPENDENT`, `L1_MED_DELETION`, `L1_MED_DUPLICATION`, `FOLDBACK_INVDUP_5P`, `CHIMERIC_ENDS`, `EXON_JUNCTION`; `.` if none |
| `strand` | insertion (element-sense) strand; `-` = poly-T at the RIGHT breakpoint |
| `ins_len` | length of the inserted sequence between the breakpoints (incl. poly-A, templated/fold-back/pre-mRNA extras) |
| `polya_len` | simulated poly-A length (reads see per-read jittered lengths) |
| `polya_side` | `LEFT` / `RIGHT` / `-` (no poly-A) |
| `tsd_seq` | TSD sequence on the reference + strand |
| `en_motif` | 7-mer around the nick, element orientation (3 bp 5' + 4 bp 3' of the nick; consensus `.TT|AAAA`) |
| `en_mismatches` | mismatches of the 6-mer vs `TT|AAAA` (0-3 sampled: 25/35/25/15 %) |
| `mh_seq` | 5' microhomology (L1-mediated deletion) |
| `element_id` / `source_id` / `subfamily` | library element (rmsk/L1Base id), transduction source, subfamily |
| `samples` | 1-based samples carrying the event (artefacts: exactly one) |
| `vaf_by_sample` | per-sample VAF (0 = absent; artefacts 1.0 where present) |
| `frags_R` / `frags_L` | per-sample count of **independent fragments** crossing the RIGHT / LEFT junction (>=10 bp each side). val1 only (fullstack: `<sample>.support.tsv` from `simulate_reads.py`, merged by `score_types.py`) |
| `reads_R` / `reads_L` | same, counting reads incl. unflagged PCR duplicates |
| `parts` | `label:start-end[info]` segments of the inserted sequence (element sense): `L1`, `L1_INV`, `L1_SWITCH`, `ALU`, `SVA`, `TD3P[source:start-end]`, `TD5P`, `SRC_TAIL`, `EXON<n>[gene]`, `INTRON`, `POLYA`, `TEMPLATED[ref:a-b±]`, `PREMRNA[...]`, `FOLDBACK[...]` |
| `info` | `k=v;...` model parameters (twin-priming breakpoint/junction, transduction end, deletion size, artefact kind, ...) |

`<truth>.ins.fa` (val1) / `truth_types.ins.fa` (fullstack) hold the inserted sequences.

## type_id catalogue

| id | type | keys |
|---|---|---|
| 1 | solo L1 | `L1_FULL`, `L1_TRUNC`, `L1_INV` |
| 2 | partnered 3' transduction | `L1_TD3P` |
| 3 | orphan 3' transduction | `ORPHAN_TD3P` |
| 4 | Alu / SVA (incl. SVA TDs) | `ALU_YA5`, `ALU_YB8`, `SVA_E`, `SVA_F`, `SVA_TD5P`, `SVA_TD3P` |
| 5 | processed pseudogene (+ decoy) | `PSEUDOGENE`, `PSEUDOGENE_DECOY` (element `UNKNOWN`, must NOT be called PSEUDOGENE) |
| 6 | solitary poly(A/T) | `POLYA_ONLY` |
| 7 | L1-mediated deletion | `L1_MED_DELETION` |
| 8 | L1-mediated tandem duplication | `L1_MED_DUPLICATION` |
| 9 | RT-mediated rearrangement | only its deletion/duplication-like forms, as 7/8 |
| 10-12 | reciprocal translocation bridge / chimeric bridge / reciprocal inversion-complex | out of scope (not simulated) |
| 13 | twin priming + 5' switching | `L1_INV_SWITCH` |
| 14 | templated local insertion | `TEMPLATED_LOCAL` |
| 15 | co-inserted local pre-mRNA | `PREMRNA_COINSERT` |
| 16 | fold-back inverted duplication 5' of the site | `FOLDBACK_INVDUP_5P` |
| A | target-site deletion instead of a TSD | `L1_TSD_DELETION` |
| B | EN-independent insertion | `EN_INDEPENDENT` |
| 0 | artefacts | `ART_LIGATION_CHIMERA`, `ART_LIGATION_PCR`, `ART_LONG_TSD`, `ART_CHIMERIC_ENDS`, `ART_POLYA_SLIPPAGE`, `ART_SUBFAMILY_MISMAP` (val1 only), `ART_FOLDBACK_PALINDROME` |

Ids match the section ids of `docs/insertion_types.html` (1-16, A, B); `TYPE_IDS`/`_KEY_TYPE_ID`
in `test/simlib/models.py` are the single place to change them.
