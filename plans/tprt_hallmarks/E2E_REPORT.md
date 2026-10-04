# TPRT-hallmark pipeline — end-to-end report (integration worker, 2026-10-04)

Branch `tprt-hallmarks` + integration commits (worktree `agent-abc17d3ebb5372ddc`), including the
merged discovery pairing modes (`5b70145`) and the Tubio-2014 library update (`b36237b`).

## What was run

One simulated patient, **3 colonies**, fullstack simulator (`test/fullstack/run_multisample.sh`):
all 22 in-scope catalogue TP types x 8 + 6 artefact classes x 8 (224 events) planted in hs1
`chr22:19,000,000-33,000,000`, real elements / transduction flanks from `resources/rte_library/`,
15x per colony, Illumina noise incl. poly-A jitter, **10 % unflagged PCR duplicates with +-3 bp
start/end jitter** (new `--pcr-dup-jitter`), bwa-mem to a reduced GRCh38 (chr22 + hs1 source
decoys) -> Rust discovery with `cluster/config.discovery.grch38.tprt` (only `splice_hallmark` off:
the GRCh38 exon track lives on the farm) -> `combine_insertions` with `cluster/config.py.grch38.tprt`
(gate >= 2 independent fragments pooled over colonies, lenient 2nd dedup, indel-aware consensus,
SHORT overhang reads, far-flank trim) -> `annotate_v2` + `tools/rte` with `resources/rte_library`
(relative path, resolved at runtime) -> per-type scoring + `tools/rte/calibrate.py`.
Wall time ~8 min (simulation + mapping + discovery 6 min, combine 11 s, annotate 75 s).

Local substitutions (see `test/e2e/make_config.py`): combine's bowtie2 index and the clip re-map
genome are the reduced reference itself (no hs1 bowtie2 index locally) with an identity liftover
chain; annotate's rmsk is hg38 RepeatMasker converted to `.out`; Dfam via `nhmmscan` on the
5-family test HMM; **no exon track** (PSEUDOGENE cannot be proven); genotyping is stubbed (every
combined locus het in every colony). One-sided loci are not genotyped yet (combine excludes
`TYPE_*_DISC` from genotyping.txt.gz) — not implemented here.

### Exact commands

```bash
cd /Users/jeremy/Documents/PEAR_TREE/.claude/worktrees/agent-abc17d3ebb5372ddc
rm -rf /private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/fde0700f-e325-4651-8daf-0cdd52bd072b/scratchpad/work/e2e
bash /Users/jeremy/Documents/PEAR_TREE/.claude/worktrees/agent-abc17d3ebb5372ddc/test/e2e/run_e2e.sh
```

`run_e2e.sh` defaults: `SP=/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/fde0700f-e325-4651-8daf-0cdd52bd072b/scratchpad`,
`PY=$SP/venv/bin/python`, genomes in `$SP/genomes` (hs1.2bit, hg38.2bit, hs1.repeatMasker.out.gz,
hg38.rmsk.txt.gz), `OUT=$SP/work/e2e`, `SAMPLES=3 N_PER_TYPE=8 DEPTH=15 PCR_DUP=0.10
PCR_DUP_JITTER=3 SEED=7 REGION=chr22:19000000-33000000`. A re-run skips simulation/discovery when
`$OUT/score.txt` exists — that is how the combine/annotate iterations below were measured on
identical discovery output. Outputs: `$OUT/e2e_tables.md`, `$OUT/e2e_events.tsv` (per event),
`$OUT/calibrate.{fit,eval}.txt`, `$OUT/combine/combine.log`, `$OUT/annot/P1.annotated.tsv`.

## Per-type results (final run)

`disc` discovery recall (some colony has both breakpoints within 30 bp; brackets: found only as a
one-sided `oneside_` locus) · `ev` reached combine's evidence evaluation · `pooled>=2` both
junctions >= 2 independent fragments pooled over colonies · `comb` survives combine (brackets:
one-sided) · elem/struct/strand annotate correct, of `comb` · tag recall (names before `=`) ·
`src id` exact `TD3P_SOURCE` · beyond = annotate recovered >= 10 bp across the poly-A (of events
with a poly-A), and with >= 2 independent fragments · T/L/U/A = TPRT / LIKELY_TPRT / UNCERTAIN /
ARTEFACT_LIKE.

| variant | role | n | lifted | disc (1-sided) | ev | pooled>=2 | comb (1-sided) | elem | struct | strand | tag recall | extra tags | src id | beyond (>=10bp) | beyond >=2 frag | T/L/U/A |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ALU_YA5 | TP | 8 | 8 | 8 (0) | 8 | 8 | 8 (1) | 8/8 | 7/8 | 8/8 | . | 8 | . | 7/8 | 6/8 | 7/0/1/0 |
| ALU_YB8 | TP | 8 | 7 | 7 (0) | 7 | 7 | 7 (0) | 7/7 | 7/7 | 7/7 | . | 0 | . | 7/7 | 6/7 | 7/0/0/0 |
| EN_INDEPENDENT | TP | 8 | 8 | 8 (0) | 8 | 8 | 8 (0) | 8/8 | 8/8 | 8/8 | 4/8 | 1 | . | . | . | 0/4/4/0 |
| FOLDBACK_INVDUP_5P | TP | 8 | 7 | 5 (0) | 3 | 2 | 3 (1) | 3/3 | 1/3 | 3/3 | 0/3 | 2 | . | 3/3 | 2/3 | 3/0/0/0 |
| L1_FULL | TP | 8 | 6 | 6 (0) | 6 | 6 | 6 (0) | 6/6 | 6/6 | 6/6 | . | 1 | . | 2/6 | 2/6 | 6/0/0/0 |
| L1_INV | TP | 8 | 8 | 7 (0) | 7 | 7 | 7 (0) | 7/7 | 5/7 | 5/7 | . | 0 | . | 4/7 | 3/7 | 5/2/0/0 |
| L1_INV_SWITCH | TP | 8 | 8 | 8 (0) | 8 | 8 | 8 (0) | 8/8 | 0/8 | 8/8 | . | 4 | . | 7/8 | 7/8 | 8/0/0/0 |
| L1_MED_DELETION | TP | 8 | 8 | 8 (1) | 8 | 6 | 7 (1) | 7/7 | 6/7 | 7/7 | 6/7 | 0 | . | 4/7 | 3/7 | 5/1/1/0 |
| L1_MED_DUPLICATION | TP | 8 | 7 | 5 (1) | 5 | 4 | 5 (2) | 4/5 | 3/5 | 5/5 | 3/5 | 0 | . | 0/5 | 0/5 | 1/1/3/0 |
| L1_TD3P | TP | 8 | 8 | 8 (0) | 8 | 8 | 8 (1) | 7/8 | 7/8 | 8/8 | 15/16 | 0 | 7/8 | 6/8 | 4/8 | 7/0/1/0 |
| L1_TRUNC | TP | 8 | 8 | 7 (0) | 7 | 7 | 7 (0) | 7/7 | 7/7 | 7/7 | . | 2 | . | 6/7 | 5/7 | 7/0/0/0 |
| L1_TSD_DELETION | TP | 8 | 8 | 7 (0) | 7 | 6 | 6 (0) | 6/6 | 5/6 | 6/6 | 5/6 | 4 | . | 4/6 | 4/6 | 4/1/1/0 |
| ORPHAN_TD3P | TP | 8 | 7 | 5 (0) | 5 | 5 | 5 (0) | 5/5 | 5/5 | 5/5 | 5/10 | 2 | 5/5 | 4/5 | 4/5 | 5/0/0/0 |
| POLYA_ONLY | TP | 8 | 8 | 5 (4) | 4 | 1 | 4 (3) | 1/4 | 3/4 | 4/4 | . | 0 | . | 2/4 | 1/4 | 2/1/1/0 |
| PREMRNA_COINSERT | TP | 8 | 7 | 6 (0) | 5 | 5 | 5 (0) | 5/5 | 0/5 | 5/5 | 0/5 | 6 | . | 5/5 | 4/5 | 4/0/1/0 |
| PSEUDOGENE | TP | 8 | 8 | 6 (0) | 6 | 6 | 6 (1) | 0/6 | 0/6 | 6/6 | 0/6 | 6 | . | 4/6 | 3/6 | 5/1/0/0 |
| PSEUDOGENE_DECOY | TP | 8 | 8 | 8 (0) | 8 | 8 | 8 (0) | 3/8 | 6/8 | 8/8 | . | 6 | . | 5/8 | 4/8 | 8/0/0/0 |
| SVA_E | TP | 8 | 7 | 6 (0) | 6 | 6 | 6 (0) | 6/6 | 6/6 | 6/6 | . | 0 | . | 5/6 | 4/6 | 5/0/1/0 |
| SVA_F | TP | 8 | 7 | 7 (0) | 7 | 7 | 7 (0) | 7/7 | 7/7 | 7/7 | . | 3 | . | 7/7 | 6/7 | 7/0/0/0 |
| SVA_TD3P | TP | 8 | 7 | 6 (0) | 6 | 5 | 5 (0) | 5/5 | 5/5 | 5/5 | 10/10 | 3 | 5/5 | 5/5 | 3/5 | 3/0/2/0 |
| SVA_TD5P | TP | 8 | 8 | 7 (0) | 7 | 7 | 7 (0) | 6/7 | 0/7 | 7/7 | 6/14 | 5 | . | 5/7 | 5/7 | 5/1/1/0 |
| TEMPLATED_LOCAL | TP | 8 | 8 | 5 (0) | 2 | 2 | 2 (0) | 2/2 | 2/2 | 2/2 | 1/2 | 1 | . | 1/2 | 1/2 | 1/0/0/1 |
| **all TP** | TP | 176 | 166 | 145 (6) | 138 | 129 | 135 (10) | 118/135 | 96/135 | 133/135 | 55/92 | 54 | 17/18 | 93/127 | 77/127 | 105/12/17/1 |

Discovery recall went from 110/166 (66 %, before merging the pairing modes) to **145/166 (87 %)**;
the new modes recover L1_TSD_DELETION 1→7, L1_MED_DELETION 0→8, L1_MED_DUPLICATION 0→5,
EN_INDEPENDENT 0→8, POLYA_ONLY 0→5 (4 of them one-sided). Post-combine 135/166 (81 %).
Strand is right for 133/135; element 118/135 (misses: PSEUDOGENE 0/6 — no exon track in this run,
so PSEUDOGENE_CANDIDATE at best; PSEUDOGENE_DECOY truth is UNKNOWN, called UNKNOWN 3/8).

## Are artefacts still being called? (per ART_* class)

| class | events | called at discovery | survive >=2 gate | survive combine | tprt_call | why it died / survived |
|---|---|---|---|---|---|---|
| ART_LIGATION_CHIMERA | 8 | 0 | 0 | 0 | — | never paired at discovery |
| ART_LIGATION_PCR | 8 | 2 (both one-sided) | 0 | 0 | — | 5 resp. 4 fragments = 1 molecule: lenient dedup collapsed 4+3 copies -> 1 independent |
| ART_CHIMERIC_ENDS | 8 | 2 | 0 | 0 | — | 2nd locus survived in the first run (2+2 independent, LIKELY_TPRT): two PCR copies whose multi-mapped mates landed on different paralogs on opposite strands. Fixed (mate sequences compared allele-forward, sequence-direction-aware poly-A cut) |
| ART_LONG_TSD | 8 | 2 (paired as L1DUP) | 0 | 0 | — | 2nd locus survived the first runs (2+2, ARTEFACT_LIKE via `tsd_gt50`): its PCR copy carried the evidence on the OTHER read (CLIP in one copy, DISC + MAPQ-0 mate in the other) and the dedup loop pruned by the primary read's outer. Fixed (`Fragment.swapped`, no outer pruning) |
| ART_POLYA_SLIPPAGE | 8 | 0 | 0 | 0 | — | simulator gives 0 junction fragments (under-represented; see organic slippage below) |
| ART_FOLDBACK_PALINDROME | 8 | 0 | 0 | 0 | — | never paired, also not as one-sided locus -> the proposed combine fold-back reject was **not** added (nothing to remove) |
| ART_SUBFAMILY_MISMAP | — | — | — | — | — | only in the val1 simulator, **not evaluated here** |

**Answer: none of the 48 simulated artefacts survives combine** (0/6 discovered artefact loci pass
the pooled >= 2 gate). No SHORT read rescued an artefact junction. Systematic artefacts are
the remaining problem — see "organic" calls below.

### Unexplained ("organic") calls

Discovery emits 1,758 loci over 3 colonies, 1,448 not matching any planted event; combine keeps
**242** unexplained loci (of 387; the breakdown below counts every combined name not matched to an event, 252):

| geometry (gap = R - L) | 1 colony | 2 colonies | 3 colonies |
|---|---|---|---|
| TSD 2-40 | 47 | 35 | 20 |
| target-site deletion | 19 | 13 | 7 |
| blunt 0-1 | 0 | 2 | 4 |
| L1DEL (< -30) | 45 | 7 | 0 |
| L1DUP (> 40) | 20 | 0 | 0 |
| one-sided | 31 | 1 | 1 |

Loci present in all 3 colonies (32) are very likely real CHM13-vs-GRCh38 germline differences
(reads are simulated from hs1 haplotypes and mapped to GRCh38; e.g. `chr22:21913684-21913685`, an
Alu with 24+16 fragments in all 3 colonies) — not artefacts. Single-colony loci are mostly
**multi-fragment poly-A slippage at reference Alu tails** (simulator slippage tracts: clip =
poly-T + shifted reference) and, new with the pairing modes, **far L1DEL/L1DUP pairs of unrelated
breakpoints** (45 + 20). Annotate calls of the unexplained loci: TPRT 25, LIKELY 43, UNCERTAIN
135, ARTEFACT_LIKE 39; the L1DEL class is the worst (15 TPRT). The >= 2 gate cannot remove them
(they have >= 2 independent molecules by construction).

## Evidence pooling, dedup, SHORT reads

* Pooling follows combine's own grouping (`Insertion.member_loci`); in this run only 4 events were
  split into colony-specific breakpoint variants (`chr22:29090478-491` / `-497` etc.), because
  `intersect_insertions` still merges exact names only (I2-owned); 5 loci were dropped by its
  clip-agreement check (score < 0.6, poly-A clips of different lengths).
* Lenient 2nd dedup vs qname truth (simulator names PCR copies `_d<k>`), same discovery output:

| combine version | true PCR-copy frags | dups missed | false merges | artefacts surviving combine |
|---|---|---|---|---|
| lenient dedup as first specified | 1,577 | 164 | 64 | 2 |
| + multi-mapped mate -> sequence; sequence-direction-aware poly-A cut; raw OR RLE edit | 1,552 | 44 | 137 | 1 |
| + `Fragment.swapped`, no outer pruning | 1,554 | 14 | 298 | 0 |
| + outer-in-clip tolerance 2x tol instead of skipped (**final**) | 1,554 | **14** | **217** | **0** |

  Final: n_independent equals the true molecule count at 2,567/2,723 junctions; 217 false merges
  out of 11,210 evidence fragments (1.9 %) cost no TP (135 combined in every row). The residual
  false merges are inherent: two molecules whose reads start within +-5 bp and whose mates are
  unplaceable inside the same element are indistinguishable from PCR copies.
* SHORT overhang reads: first implementation compared the overhang to the reference
  position-by-position; an unclipped SHORT read carrying a slipped poly-T as a 4-5 bp deletion
  next to the junction then looked "non-reference", and **157 junctions (53 unexplained combined
  loci) reached >= 2 only thanks to SHORT reads**. With an indel-aware (edlib infix) reference
  test: 182 SHORT fragments used, 5,747 rejected (matches_reference 2,338, consensus_mismatch
  2,022, overhang_too_short 1,357, ref_homopolymer 11), **18 junctions** reach >= 2 only thanks to
  them (1 TP locus, 7 unexplained, 0 artefacts).

## TPRT score separation (held-out half)

No simulated artefact reaches annotate any more, so separation is measured TP vs the
unexplained loci (an impure negative set: it contains real germline assembly differences).
Split: events by id parity, unexplained loci by locus hash.

| half | AUC TP vs ARTEFACT+unexplained | TPR / FPR at TPRT (>= 7) | TPR / FPR at LIKELY (>= 4) |
|---|---|---|---|
| fit | 0.931 | 0.83 / 0.09 | 0.91 / 0.25 |
| eval (held out) | 0.894 | 0.76 / 0.12 | 0.88 / 0.32 |

TP calls (combined): TPRT 105, LIKELY 12, UNCERTAIN 17 (incl. EN_INDEPENDENT 4 — by design not
TPRT, and L1_MED_DUPLICATION 3), ARTEFACT_LIKE 1.

**No weights were recalibrated.** The fit half shows `tsd_deletion_le20` (+1, log-odds -1.05) and
`tsd_other_1_50` (+1, -0.45) pointing the wrong way, but the negatives are contaminated by real
germline insertions and by pairing-mode FPs that have no artefact model in the simulator; changing
weights on that would be overfitting. Two structural (not weight) fixes were made instead:
SVA 5' junction pieces in the SVA Alu-like domain no longer raise CHIMERIC_ENDS (4/12 combined SVA_TD3P/TD5P
TPs were capped by `chimeric_ends:-6`), and the L1DUP-vs-long-TSD rule (below).

## Strand / LEFT-RIGHT convention (task 3)

Verified consistent end to end, nothing to fix: discovery names loci `contig:L-R` (Rust
`locus_name`: LEFT-clip breakpoint first; the `src/discovery.py` header comment saying `R-L` is
stale), LEFT = left-clipped reads (`[clip][aligned]`, clip = insertion end), so a `+` element has
its poly-A at the LEFT junction and a `-` element poly-T at the start of the RIGHT clip; the
simulator truth (`models.apply_event`: `polya_side=LEFT` for `+`, `RIGHT` for `-`) and annotate
(`hallmarks.polya_info`) agree; strand correct 133/135 in the E2E. The reads.fa orientation was
**not** as SPEC promised: combine wrote MATE/POLYA rows as stored in the BAM (a mate multi-mapped
onto a paralog is stored relative to that paralog). Fixed in combine (`allele_forward_seq`, pair
geometry from the 0x20 / 0x10 flags).

## Remaining problems, ranked by impact

1. **Organic multi-fragment FPs** (242 unexplained combined loci, 68 called TPRT/LIKELY): poly-A
   slippage at reference A-tracts with >= 2 slipping molecules, and far L1DEL/L1DUP pairs of
   unrelated breakpoints from the new pairing modes. Need a combine-level reference-context
   slippage reject (clip = homopolymer + shifted local reference) and a stricter polarity /
   element-class check for far pairs; the score's `slippage_context` fires on only 30 %.
2. **Structure calls**: L1_INV_SWITCH 0/8 (called INVERTED_5P/TRUNCATED), SVA_TD5P 0/7,
   PREMRNA_COINSERT 0/5, FOLDBACK_INVDUP_5P 1/3; tag precision is low for TD3P (38 extra),
   TEMPLATED_LOCAL (24 extra), L1_MED_DELETION/TSD_DELETION on unexplained loci.
3. **intersect_insertions exact-name merging** (I2): colony-specific breakpoint variants stay
   separate insertions (pooling code is ready to follow any fuzzy merge); 5 loci lost to the
   clip-agreement check on poly-A clips.
4. Dedup false merges 1.9 % (inherent at +-5 bp with unplaceable mates); consider tol 3.
5. ART_SUBFAMILY_MISMAP not evaluated (val1-only artefact; the fullstack simulator lacks it).
   Annotate features that should catch it: `inactive_only` does not apply (young subfamily);
   proposal: flag a clip whose flank part maps uniquely to a reference L1 3' flank elsewhere while
   the element part sits on a reference L1 copy here (mismap signature), rather than a score tweak.
6. ~~One-sided loci are not genotyped~~ — now genotyped, see "Genotyping new locus kinds" below;
   `keep_polya_one_sided` can now be enabled (the gate skips the open side) but was left off
   pending a measurement.

## Genotyping new locus kinds (genotyping worker, 2026-10-04)

The E2E above stubbed genotyping. This section genotypes the same combine output with the Rust
genotyper, per colony, and scores present/absent per colony against the simulator truth
(`vaf_by_sample`; present = VAF > 0).

### How the genotyper builds alleles, and what each locus kind needs

Per read, per junction: `qleft` compares the <= 12 read bases 5' of L with `LEFT_REFERENCE`
(= genome[L-12, L)) and `LEFT_INSERTION` (the left clip's 12 junction-adjacent bases); `qright`
the bases 3' of R with `RIGHT_REFERENCE` (genome[R, R+12)) and `RIGHT_INSERTION`. A read with
an alt junction is one alt vote, a read with only ref junctions one ref vote; VAF bands decide.
No allele is built "between" L and R, so:

* **target-site deletion** (`R < L`, <= 30 bp) and **blunt** (`R - L` in {0, 1}): correct
  unchanged — an alt read covers only its own junction (no false double-alt), a reference read
  spans both and votes ref. Verified on synthetic reads (Rust unit tests) and in the E2E.
* **L1-mediated deletion / duplication** (`|R - L|` up to 50 kb): the legacy single window
  `[min, max]` counts every read of the span as depth -> `high-coverage` (6/18 TP colony calls,
  5/6 in the wild-type control), and a mate inside the span can claim a spanning read's qname.
  Fixed by `split_breakpoint_span = 40`: one 1-bp window per breakpoint, shared qname dedup.
  Two further biases, both from the VAF bands being calibrated on a TSD locus (alt reads from
  two junctions per reference span, het VAF = 2f/(2f+1), f = junction-read yield):
  far pairs have two *disjoint* reference spans and one-sided loci one junction, so a het reads
  f/(1+f) (E2E pooled VAF 0.30-0.38) -> `halve_single_junction_ref` (ref votes count half);
  and a far **duplication**'s alt haplotype `ref[..R) + element + ref[L..)` still carries both
  reference junctions, so alt-haplotype molecules also yield reference reads -> 
  `dup_ref_discount_min_span = 150` (`n_ref -= min(n_ref, n_alt)`; wild-type colonies unchanged).
* **one-sided** (`contig:L-oneside_L` / `contig:oneside_R-R`): combine leaves them out of
  `<patient>.genotyping.txt.gz` (I don't own combine), so `src/genotyping_contract_oneside.py`
  appends them -> `<patient>.genotyping.tprt.txt.gz`: real side only, built exactly as combine
  builds that side (clip from `combined.txt.gz`, reference from the 2bit, combine's per-side
  exclusions; 42 added, 1 excluded as clip == reference). The genotyper (`one_sided_loci`)
  parses the token and never scores the open end. The missing end's junction reads (soft clip
  on the open side, <= 50 bp from the real breakpoint) cross the real breakpoint in reference
  configuration and would be false ref votes -> skipped (`one_sided_open_window = 50`).

`cluster/pipeline.sh` builds the extended contract in the (single) combine task when
`GENO_CFG` sets `one_sided_loci = true`, and the genotype tasks read it. Use
`GENO_CFG=cluster/config.genotype.grch38.tprt`. `combine_genotypes.py` needed no change to
accept the names (rows are keyed verbatim — what annotate joins on); it now prints a per-kind
(`locus_kind`) removed/kept summary; `<patient>.genotypes.csv.gz` is unchanged.

Byte-identity: default config and `config.genotype.grch38` produce md5-identical output to
the pre-change binary on S1-S3, a wild-type control colony and `test_data/test.bam`; the `.tprt`
config on the legacy contract changes only the 80 far-pair rows (TSD / TSD-deletion / blunt
rows identical).

### Commands

```bash
cd /Users/jeremy/Documents/PEAR_TREE/.claude/worktrees/agent-a76527686f1a0a72a
bash /Users/jeremy/Documents/PEAR_TREE/.claude/worktrees/agent-a76527686f1a0a72a/test/e2e/run_genotype_e2e.sh
```

(needs the `run_e2e.sh` output in `$SP/work/e2e`; simulates a 4th, **wild-type control colony
S4** from the event-free haplotype with the same simulator/mapping, extends the contract,
genotypes S1-S4 with the legacy and the `.tprt` config, writes
`$SP/work/e2e/genotype/genotype_score.md`.)

### Results (135 combined TP loci x 3 colonies; S4 = wild-type control)

Colony calls: present = het/hom/`insertion` (combine_genotypes' carrier set), absent =
`wild-type`, no call = anything else (all no-calls here are truth-present colonies; most are
subclonal VAF 0.125-0.375 reading `wild-type?`).

| locus kind | TP loci | legacy: present ok / ->absent / no call | **tprt: present ok / ->absent / no call** | legacy S4 FP / no call | **tprt S4 FP / no call** | het (VAF 0.5) pooled VAF legacy -> tprt |
|---|---|---|---|---|---|---|
| TSD 2-40 | 101 | 188 / 18 / 55 (+42 absent ok) | 188 / 18 / 55 (+42 absent ok) | 0 / 1 | 0 / 1 | 0.449 -> 0.449 |
| target-site deletion | 7 | 13 / 2 / 2 (+2 ok, 2 absent->present) | identical | 1 / 0 | 1 / 0 | 0.458 -> 0.458 |
| blunt 0-1 | 8 | 20 / 1 / 3 | identical | 0 / 0 | 0 / 0 | 0.535 -> 0.535 |
| L1DEL (< -30) | 6 | 10 / 0 / 8 (6 `high-coverage`) | **15 / 0 / 3** | 0 / 5 (`high-coverage`) | **0 / 0** | (high-cov) -> 0.544 |
| L1DUP (> 40) | 3 | 2 / 1 / 6 | **8 / 0 / 1** | 0 / 0 | 0 / 0 | 0.273 -> 0.509 |
| one-sided | 10 | 0 / 0 / 30 (not in contract) | **20 / 6 / 4** | — | **0 / 0** | — -> 0.446 |
| **all** | 135 | 233 / 22 / 104 | **264 / 27 / 68** | 1 / 16 | **1 / 1** | |

Key-by-key on the same data (one-sided present ok of 30 / L1DEL of 18 / L1DUP of 9):
`one_sided_loci` + `split_breakpoint_span` only 10 / 12 / 2 (one-sided het pooled VAF 0.27)
-> + `one_sided_open_window` and `dup_ref_discount_min_span` 16 / 12 / 3 (0.30)
-> + `halve_single_junction_ref` **20 / 15 / 8** (0.45). No new-kind TP locus is called present
in the wild-type control (0/19), so none of the corrections buys sensitivity with false carriers.
The 2 TSD-deletion absent->present calls and the S4 FP are one locus, `chr22:26440278-26440276`
(matched to an SVA_TD5P event present only in S1), homozygous in all four colonies including
the control: a germline CHM13-vs-GRCh38 difference at the event site, identical under both
configs.

Remaining one-sided misses (6 present->absent, 4 no-call) are alt-read mappability, not
genotyper logic: the real junction is the poly-A side, its junction reads start in the poly-A /
element and lose MAPQ 60 in Alu-rich flanks (e.g. `chr22:32246706-oneside_32246706`: alt reads
at MAPQ 34 / 4, 0 of them pass `min_mapq = 60` in S3), plus subclonal colonies.

Unexplained (no TP event) far pairs and one-sided loci are now genotyped too: of 52 L1DEL /
19 L1DUP unexplained loci, 12 / 7 are present in all three colonies **and** in the wild-type
control (germline-like CHM13-vs-GRCh38 differences — a real cohort's `min_wild-types` removes
them), single-colony ones are mostly absent in S4; the 3-colony E2E cannot run the clade gates
meaningfully (`min_wild-types 20`), but `combine_genotypes` was run on the 4 genotype files and
accepts every new name.

Caveats: (1) with `halve_single_junction_ref` the `n_ref` column of one-sided / far-pair rows
is a half-weight count — the dispersion gate sees it as such (fewer effective trials,
conservative). (2) annotate_v2's `read_genotyping` counts only het/hom as carriers, not
`insertion` (presence certain, zygosity unclear), unlike combine_genotypes — low-coverage
one-sided loci reading `insertion` everywhere would be dropped there (annotate-owned).
(3) The Python oracle genotyper (`src/genotype.py`) was not extended and rejects `oneside_`
names; the cluster uses the Rust binary.
