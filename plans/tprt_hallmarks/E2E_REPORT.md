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

1. ~~**Organic multi-fragment FPs**~~ (242 unexplained combined loci, 68 called TPRT/LIKELY) —
   addressed, see "Combine filters (problems #1 / #3)" below: 242 -> 97 (seed 7), 229 -> 113
   (seed 11); what is left is ~90 % genuine hs1-vs-GRCh38 germline differences.
2. **Structure calls**: L1_INV_SWITCH 0/8 (called INVERTED_5P/TRUNCATED), SVA_TD5P 0/7,
   PREMRNA_COINSERT 0/5, FOLDBACK_INVDUP_5P 1/3; tag precision is low for TD3P (38 extra),
   TEMPLATED_LOCAL (24 extra), L1_MED_DELETION/TSD_DELETION on unexplained loci.
3. ~~**intersect_insertions exact-name merging**~~ — fuzzy merge (`merge_tolerance_bp`) and
   poly-A-aware clip agreement added, see below.
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

## Annotate round 2 (annotate worker, 2026-10-04)

Same discovery + combine output as the final run above (`$SP/work/e2e/combine`, copied); only
`annotate_v2` + `tools/rte` re-run (`test/e2e/run_annotate_e2e.py` + `score_e2e.py`). Harness fix
first: `score_e2e.py` picked the combined name of an event from a Python `set` (hash order ->
element/structure counts moved by +-3 between identical runs); it now prefers two-sided names, then
sorts by name. "base" below = `tprt-hallmarks` 4681f72 re-scored with that deterministic scorer
(120/99/56 instead of the 118/96/55 quoted above).

### Pseudogene parents: the same exons for truth and annotate

The simulator built its pseudogene parents as synthetic GT..AG genes on real hs1 chr22 sequence
(no `--gene-model`), so no public track contains them. `test/fullstack/donor_types.py` now writes
them (`donor/genes_hs1.tsv`, build_gene_model format); `test/e2e/make_gene_track.py` reads that, or
for older runs recovers the exons from the truth `EXON<n>[gene]` parts (34 exons / 8 genes here),
and writes (a) the hs1 track for the exon-exon junction cores (new key `rte_exon_annotation`,
because `remap_2bit` is hs1) and (b) the same exons placed on the E2E clip-remap reference (reduced
GRCh38) for annotate_v2's clip-exon candidates (`exon_annotation`). `run_e2e.sh` does both.
**Real runs** use ONE hs1 track for both (the clip-remap genome is hs1; `rte_exon_annotation`
defaults to `exon_annotation`):

```bash
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz
python tools/build_gene_model.py --curated hs1.ncbiRefSeq.gtf.gz hs1.gene_model.tsv.gz
# CONFIG['annotate']['exon_annotation'] = <staged>/hs1.gene_model.tsv.gz ; ['remap_2bit'] = <staged>/hs1.2bit
```

(28,776 genes / 252,903 merged exon rows; checked locally, e.g. the GAPDH transcript 5' end; the
file lives in `$SP/genes`, not the repo, and still has to be staged on the farm.) PSEUDOGENE
requires an exon-exon junction read; the no-junction decoy gets `PSEUDOGENE_CANDIDATE` only (0/8
decoys called PSEUDOGENE).

### Per-type accuracy (135 combined TPs; base -> round 2)

| variant | comb | element | structure | tag recall | extra tags | T/L/U/A (r2) |
|---|---|---|---|---|---|---|
| ALU_YA5 | 8 | 8/8 -> 8/8 | 8/8 -> 8/8 | . | 6 -> 0 | 8/0/0/0 |
| ALU_YB8 | 7 | 7/7 -> 7/7 | 7/7 -> 7/7 | . | 0 -> 0 | 7/0/0/0 |
| EN_INDEPENDENT | 8 | 8/8 -> 8/8 | 8/8 -> 8/8 | 4/8 -> 8/8 | 1 -> 0 | 0/0/8/0 |
| FOLDBACK_INVDUP_5P | 3 | 3/3 -> 3/3 | 1/3 -> 1/3 | 0/3 -> 2/3 | 2 -> 0 | 3/0/0/0 |
| L1_FULL | 6 | 6/6 -> 6/6 | 6/6 -> 6/6 | . | 1 -> 0 | 6/0/0/0 |
| L1_INV | 7 | 7/7 -> 7/7 | 5/7 -> 5/7 | . | 0 -> 0 | 5/2/0/0 |
| L1_INV_SWITCH | 8 | 8/8 -> 8/8 | 0/8 -> 4/8 | . | 4 -> 0 | 8/0/0/0 |
| L1_MED_DELETION | 7 | 7/7 -> 7/7 | 6/7 -> 6/7 | 6/7 -> 6/7 | 0 -> 0 | 5/1/1/0 |
| L1_MED_DUPLICATION | 5 | 5/5 -> 5/5 | 4/5 -> 4/5 | 4/5 -> 4/5 | 0 -> 0 | 2/0/3/0 |
| L1_TD3P | 8 | 8/8 -> 8/8 | 8/8 -> 8/8 | 15/16 -> 14/16 | 0 -> 0 | 8/0/0/0 |
| L1_TRUNC | 7 | 7/7 -> 7/7 | 7/7 -> 7/7 | . | 2 -> 0 | 7/0/0/0 |
| L1_TSD_DELETION | 6 | 6/6 -> 6/6 | 5/6 -> 5/6 | 5/6 -> 5/6 | 4 -> 0 | 4/1/1/0 |
| ORPHAN_TD3P | 5 | 5/5 -> 5/5 | 5/5 -> 5/5 | 5/10 -> 10/10 | 2 -> 0 | 5/0/0/0 |
| POLYA_ONLY | 4 | 1/4 -> 1/4 | 3/4 -> 4/4 | . | 0 -> 0 | 2/1/1/0 |
| PREMRNA_COINSERT | 5 | 5/5 -> 5/5 | 0/5 -> 0/5 | 0/5 -> 5/5 | 6 -> 1 | 3/1/1/0 |
| PSEUDOGENE | 6 | 0/6 -> 5/6 | 0/6 -> 5/6 | 0/6 -> 5/6 | 8 -> 2 | 6/0/0/0 |
| PSEUDOGENE_DECOY | 8 | 3/8 -> 5/8 | 6/8 -> 6/8 | . | 6 -> 7 (4 = PSEUDOGENE_CANDIDATE) | 8/0/0/0 |
| SVA_E | 6 | 6/6 -> 6/6 | 6/6 -> 6/6 | . | 0 -> 0 | 5/0/1/0 |
| SVA_F | 7 | 7/7 -> 7/7 | 7/7 -> 7/7 | . | 3 -> 0 | 7/0/0/0 |
| SVA_TD3P | 5 | 5/5 -> 5/5 | 5/5 -> 5/5 | 10/10 -> 10/10 | 3 -> 3 | 3/0/2/0 |
| SVA_TD5P | 7 | 6/7 -> 6/7 | 0/7 -> 6/7 | 6/14 -> 12/14 | 5 -> 2 | 5/1/1/0 |
| TEMPLATED_LOCAL | 2 | 2/2 -> 2/2 | 2/2 -> 2/2 | 1/2 -> 1/2 | 1 -> 1 | 1/0/0/1 |
| **all TP** | 135 | **120 -> 127** | **99 -> 115** | **56/92 -> 82/92** | **54 -> 16** | 108/7/19/1 |

Strand unchanged (133/135). TPRT-score AUC (TP vs unexplained) fit 0.940 -> 0.928, eval 0.900 ->
0.900 (no weight changed; the 8 EN_INDEPENDENT TPs are UNCERTAIN by design).

### Per-tag precision / recall (combined TPs) and tags on the 252 unexplained loci

| tag | base TP/FP/FN | r2 TP/FP/FN | base P / R | r2 P / R | on unexplained (base -> r2) |
|---|---|---|---|---|---|
| TD3P | 13/18/5 | 17/3/1 | 0.42 / 0.72 | 0.85 / 0.94 | 57 -> 20 |
| TD3P_SOURCE | 17/8/1 | 17/0/1 | 0.68 / 0.94 | 1.00 / 0.94 | 48 -> 17 |
| TD5P | 6/1/1 | 6/1/1 | 0.86 / 0.86 | 0.86 / 0.86 | 1 -> 1 |
| TEMPLATED_LOCAL | 1/16/1 | 1/0/1 | 0.06 / 0.50 | 1.00 / 0.50 | 39 -> 17 |
| PREMRNA_COINSERT | 0/0/5 | 5/2/0 | - / 0.00 | 0.71 / 1.00 | 0 -> 7 |
| EN_INDEPENDENT | 4/0/4 | 8/0/0 | 1.00 / 0.50 | 1.00 / 1.00 | 2 -> 2 |
| L1_MED_DUPLICATION | 4/7/1 | 4/1/1 | 0.36 / 0.80 | 0.80 / 0.80 | 9 -> 7 |
| FOLDBACK_INVDUP_5P | 0/0/3 | 2/0/1 | - / 0.00 | 1.00 / 0.67 | 0 -> 1 |
| EXON_JUNCTION | 0/0/6 | 5/0/1 | - / 0.00 | 1.00 / 0.83 | 0 -> 7 (see "Remaining" 1) |
| TSD_DELETION / L1_MED_DELETION / CHIMERIC_ENDS | 5/2/1, 6/0/1, 0/2/0 | unchanged | | | |

The 3 remaining TD3P FPs: 2 PSEUDOGENE_DECOY events called ALU (an exon starting with an Alu 3'
end + tail, then gene sequence before the insertion's poly-A) and 1 PSEUDOGENE without a
candidate gene.

### What was wrong, per systematic miss (and the fix)

* **L1_INV_SWITCH 0/8 -> 4/8.** The simulated shape is `sense switch piece | inverted piece
  (further 3') | sense body | polyA`; the 5' junction reads REF|sense, so it was TRUNCATED_5P. New:
  a read joining a sense piece (at/after the 5' junction piece) to an anti-sense piece lying further
  3' on the consensus -> INVERTED_5P_SWITCH (`detail switch=`). The other 4 have a 300-480 bp switch
  piece that no read crosses (fragments ~350 bp): undecidable from these reads.
* **SVA_TD5P structure 0/7 -> 6/7.** The 5' junction reads REF | source 5' flank, longer than a
  fragment, so no read joins it to the hexamer -> 5P_UNRESOLVED. A 5' transduction means
  transcription started upstream, so the SVA is complete: FULL_LENGTH when the 5' junction piece is a
  sense FLANK5P. `TD5P_SOURCE=<id>` is now emitted (truth had it, annotate never did: 6 tag misses).
* **PREMRNA_COINSERT tags 0/5 -> 5/5** (structure still 0/5). The simulator copies local sequence
  400-2400 bp from the site to the 5' end (no gene involved); annotate only looked for a host-gene hit
  between element and poly-A, and only with a gene model. New: unexplained pieces are looked up in a
  +-10 kb window (`rte_wide_window`); a local template whose near end is > 250 bp away ->
  PREMRNA_COINSERT (`premrna=local:...`), <= 250 bp -> TEMPLATED_LOCAL. That template also stopped
  producing `L1_MED_DUPLICATION` (3 FPs) via annotate_v2's intrachromosomal-SV clip partner.
  Structure: no read joins the 225-620 bp template to the L1 -> 5P_UNRESOLVED.
* **PSEUDOGENE 0/6 -> 5/6.** Besides the missing track, the synthetic parents were cut from
  repeat-rich chr22, so exons carry Alu pieces and 3/6 were called ALU. An exon-exon junction read
  now wins over an element class (`rte_in_mrna=` keeps the class; TD/templated/pre-mRNA tags read
  from mRNA pieces are dropped); structure from the 5' insert vs the spliced transcript (5/6).
  Candidates = every gene with an exon under a clip remap, not only `_pseudogene()`'s (that one
  needs a clip poly-A). Remaining miss: event 105 (no candidate gene).
* **FOLDBACK_INVDUP_5P 0/3 -> 2/3**: a 5' junction read REF | LOCAL on the opposite strand whose
  template ends within 12 bp of the junction on the flank side (>= 2 fragments) is a fold-back of the
  5' flank (was TEMPLATED_LOCAL). **EN_INDEPENDENT 4/8 -> 8/8**: a blunt pair has gap 0 **or 1**
  (SPEC) and the 3' end had to be >= 60 bp short; now >= 20 bp. **ORPHAN_TD3P** now carries TD3P.

### Over-tagging: causes and the tightened rules

* **Poly-A tail noise** (most TD3P and TEMPLATED_LOCAL FPs): SBS jitter + low-quality bases after a
  long poly-A split one tail into `POLYA | junk | POLYA`; the junk became the "unexplained segment
  before the poly-A" (TD3P), and read poly-A aligned to a reference A-run 100-500 bp away became a
  LOCAL "template". `assembly._smooth_polya` merges a tail (pieces <= 20 bp between same-base runs,
  or non-flank pieces >= 60 % the tail base); templates/tags must not be low complexity.
* **TD3P** (unexplained tag): >= 30 bp, complex, directly before the TAIL poly-A (followed by REF /
  read end, not an A-run inside genomic sequence), in >= 2 fragments (`td_frags`).
  **TD3P_SOURCE**: a flank hit >= 30 bp, or >= 20 bp at >= 95 % identity in the tail position; not
  an element end (short hits that also match a consensus: an L1 3' end inside some 15 kb flank);
  seen from the 3' junction side; no class element after it. (Tried and rejected: requiring unmasked
  flank sequence -- `flanks_3p.fa` is ~53 % soft-masked in its first kb; it killed 8 true sources.)
* **TEMPLATED_LOCAL**: >= 20 bp, >= 90 % identity, complex, not inside the TSD +- 5 bp, near end
  <= 250 bp from a breakpoint, >= 2 fragments; a REF split by a read indel / homopolymer is merged
  back into one flank (was a "template"); a LOCAL piece >= 30 bp that also matches an element
  consensus (the insert next to a reference copy of its family) is ELEMENT.
* 5' junction reads whose REF piece does not end at a breakpoint (a reference element elsewhere in
  the window) no longer pick the 5' element class (a CHIMERIC_ENDS source).

### Threshold changes on the held-out split (events by id parity; P = TP/called, R = TP/truth)

Each row re-runs annotate with ONE threshold reverted / loosened (rest = round 2 final):

| threshold | value | fit | eval | unexplained loci |
|---|---|---|---|---|
| (base, all old logic) | - | TD3P 7/16 P, 7/11 R | TD3P 6/15 P, 6/7 R | TD3P 57 |
| **round 2 final** | td_min_fragments 2, td_min_bp 30, flank 30 or 20@0.95 | TD3P 11/11, 11/11 | TD3P 6/9, 6/7 | TD3P 20, TD3P_SOURCE 17 |
| td_min_fragments | 1 | TD3P 10/12, 10/11 | TD3P 6/12, 6/7 | TD3P 30 |
| td_min_bp | 20 | TD3P 10/12, 10/11 | TD3P 6/10, 6/7 | TD3P 27 |
| td_min_flank_bp | 20 (any identity) | TD3P 11/11, 11/11 | TD3P 6/9, 6/7 | TD3P_SOURCE 29 |
| templated_min_fragments | 1 (final 2) | TEMPLATED 0/4 P (final 0/0) | 2/3 P, 2/2 R (final 1/1, 1/2) | TEMPLATED 29 (final 17) |
| en_independent_3p_tolerance | 60 (old; final 20) | EN_IND 3/4 R (final 4/4) | 4/4 R (final 4/4) | 2 |

(The td_* ablations were run before the short-flank rule was added; with it, final fit TD3P is
11/11.) All final values were set a priori from the task spec (N = 30 bp, >= 2 fragments, 20 bp,
90 %, 250 bp) and confirmed by these ablations; the one rule found by looking at an event is the
short tail-position flank hit (fit-half event 42): +1 TP on fit, nothing on eval, and no extra
unexplained TD3P_SOURCE (the plain 20 bp floor adds 12). Unchanged: TSD_DELETION (2 FP: templated /
SVA TD5P breakpoint geometry), the decoy element (2 truth-UNKNOWN decoys called ALU), POLYA_ONLY
(3/4 one-sided loci called UNKNOWN/ALU).

### Remaining

1. **Parent-gene splice ghosts**: 6 unexplained far L1DEL/L1DUP loci now read PSEUDOGENE +
   EXON_JUNCTION (3 TPRT/LIKELY). Their breakpoints sit on exon boundaries of the parent gene: reads
   of the processed pseudogene map back to the parent and are clipped at the splice sites. Flag them
   (both breakpoints = exon boundaries of the candidate gene on the DISCOVERY genome; needs the
   discovery-genome gene model) -- not done.
2. L1_INV 2/7 5P_UNRESOLVED, L1_INV_SWITCH 4/8, PREMRNA structure 0/5: the informative junction is
   farther than a fragment from every evidence read.
3. annotate_v2 now counts the genotyper's `insertion` call as a carrier in TPRT mode (`rte_library`
   set, or explicit `count_insertion_call`); the legacy default is unchanged (point (2) above).
   One-sided (`contig:L-oneside_L`) and far-pair names join `genotypes.csv.gz` by exact title
   (tested). The E2E numbers above use stub genotypes (every locus het), so they are unaffected.

## Combine filters (problems #1 / #3) (combine worker, 2026-10-04)

### What the 242 unexplained calls are (read-level truth)

`test/e2e/classify_unexplained.py` (now step 7 of `run_e2e.sh`; writes `$OUT/unexplained_classes.tsv`
and `.md`) traces every junction read of an unexplained combined locus to its source. The simulator
qname names the haplotype, and the read is re-aligned (mappy sr) to that hs1-derived haplotype.
`contiguous` means the read is a faithful copy of hs1, so its GRCh38 clip is an hs1-vs-GRCh38
difference. `slipped` means an indel, end clip or partial alignment on its own source (SBS
homopolymer slippage or post-homopolymer phasing junk). `planted_junction` means the read spans a
planted event's junction. Per locus the table also gives the hs1-vs-GRCh38 alignment of an 800 bp
window, the reference repeat at each junction, and library hits of each clip.

Seed 7, before (= `.tprt` with the new keys off, byte-identical to the integration run):

| cause | 1 colony | 2 | 3 | total | TPRT / LIKELY / UNC / ART |
|---|---|---|---|---|---|
| far pair: planted TP junction + unrelated breakpoint | 28 | 3 | 0 | 31 | 14 / 6 / 4 / 7 |
| far pair: two unrelated breakpoints | 30 | 4 | 0 | 34 | 2 / 1 / 8 / 23 |
| planted-event reads (displaced junction) | 1 | 2 | 0 | 3 | 0 / 1 / 2 / 0 |
| sequencing slippage at a 13-25 bp reference A/T tract (25 one-sided) | 54 | 15 | 2 | 71 | 1 / 16 / 50 / 4 |
| germline: hs1 carries an Alu/L1 that GRCh38 lacks (TRUE insertion) | 0 | 3 | 8 | 11 | 5 / 3 / 2 / 1 |
| germline: STR / VNTR length difference between the assemblies | 21 | 9 | 7 | 37 | 2 / 8 / 24 / 3 |
| germline: other hs1-vs-GRCh38 indel / divergent block | 21 | 18 | 14 | 53 | 1 / 8 / 43 / 1 |
| other | 1 | 0 | 1 | 2 | 0 / 0 / 2 / 0 |

Corrections to the earlier diagnosis:

* The 3-colony loci are indeed assembly differences, but so are ~40 % of the single-colony ones
  (101/242 germline in total). They show up in a single colony because discovery paired them
  differently per colony, and intersect merged exact names only.
* Slippage is mostly not "clip = poly-T + shifted reference". The reads cross a long (13-25 bp)
  reference A/T tract and soft-clip the post-homopolymer phasing junk, which is still dominated
  by the tract base. Both junctions of a "TSD" call often sit at the two ends of the tract
  (TSD = tract length).
* Far pairs: 31/65 pair a real TP poly-A junction with an unrelated breakpoint (an assembly
  difference or a slipped tract within 50 kb). The rest pair two unrelated breakpoints. Almost
  all are in one colony.

### New keys (default off in src, on in `cluster/config.py.grch38.tprt`; details in SPEC.md)

* `slippage_reject`: applies when a reference tract (homopolymer >= 8 bp or STR >= 12 bp) touches
  the junction. After stripping the tract continuation, the clip is slippage if it is
  `repeat_only` (< 10 structured bases), `repeat_junk` (>= 50 % tract base) or
  `repeat_shifted_reference`. Such an insertion is dropped unless another, non-slipped junction
  carries element or transduction-flank sequence (mappy k=11 vs `rte_library`, or >= 2
  inside-insertion mates when the clip is < 20 bp).
* `far_pair_strict` (+ `far_pair_split`): a pair with gap < -30 or > 40 needs all of:
  * a sense element hit on the complex clip (or >= 2 inside mates);
  * a poly-T-led clip on the other side that is not itself slippage;
  * no conflicting element class beyond the tail;
  * >= 2 independent fragments per junction;
  * consistent colonies (discovery breakpoints within tol: at least one shared colony, and the
    sets differ by <= max(1, 20 %)).

  A failing pair keeps its poly-A junction as a one-sided locus. This is done in combine, where
  the library, pooled fragments and colony set live; there is no Rust change.
* `merge_tolerance_bp` (8): fuzzy cross-sample merge in `intersect_insertions`, with a
  shift-tolerant clip check; member evidence is re-anchored by the junction offset. One-sided
  loci of the same side merge with each other. A surviving one-sided locus is folded into a
  surviving two-sided call of the same junction only AFTER the remap filters
  (`EvidencePool.absorb_one_sided`). Folding earlier lost FOLDBACK #169, whose two-sided record
  died in the clipped-remap filter.
* `polya_aware_clip_agreement`: clips are compared homopolymer-compressed up to their poly-A
  tail. This recovers the poly-A loci that the legacy column check dropped.

### Before / after per type (comb, brackets = one-sided), merged deterministic scorer

| variant | seed 7 before | seed 7 after | seed 11 before | seed 11 after |
|---|---|---|---|---|
| ALU_YA5 | 8 | 8 | 8 | 7 |
| ALU_YB8 | 7 | 7 | 7 | 8 (1) |
| EN_INDEPENDENT | 8 | 8 | 7 | 7 |
| FOLDBACK_INVDUP_5P | 3 (1) | 3 (1) | 2 (1) | 3 (2) |
| L1_FULL | 6 | 6 | 8 | 8 |
| L1_INV | 7 | 7 | 5 | 5 |
| L1_INV_SWITCH | 8 | 8 | 6 | 7 (1) |
| L1_MED_DELETION | 7 (1) | 8 (2) | 7 | 7 |
| L1_MED_DUPLICATION | 5 (1) | 7 (3) | 2 (2) | 7 (7) |
| L1_TD3P | 8 | 8 | 7 | 7 |
| L1_TRUNC | 7 | 7 | 7 (1) | 7 (1) |
| L1_TSD_DELETION | 6 | 5 | 5 | 6 (1) |
| ORPHAN_TD3P | 5 | 6 (1) | 7 | 8 (1) |
| POLYA_ONLY | 4 (3) | 5 (4) | 4 (3) | 5 (4) |
| PREMRNA_COINSERT | 5 | 5 | 7 | 7 |
| PSEUDOGENE | 6 | 6 | 5 | 5 |
| PSEUDOGENE_DECOY | 8 | 8 | 7 | 7 |
| SVA_E | 6 | 7 (1) | 6 | 6 |
| SVA_F | 7 | 7 | 8 | 8 |
| SVA_TD3P | 5 | 6 (1) | 7 (1) | 7 (1) |
| SVA_TD5P | 7 | 7 | 6 | 6 |
| TEMPLATED_LOCAL | 2 | 4 (3) | 2 | 2 |
| **all TP** (of 166 / 165 lifted) | **135 (6)** | **143 (16)** | **130 (8)** | **140 (19)** |
| TP calls T / L / U / A | 108/7/19/1 | 118/8/16/1 | 99/15/16/0 | 108/13/18/1 |
| simulated artefacts combined | 0 | 0 | 1 (ART_POLYA_SLIPPAGE) | 0 |
| combined insertions | 387 | 241 | 371 | 253 |
| **unexplained combined** | **242** | **97** | **229** | **113** |
| - far pair (TP junction / unrelated) | 31 / 34 | 0 / 1 | 31 / 30 | 5 / 1 |
| - planted-event reads | 3 | 0 | 2 | 2 |
| - slippage | 71 | 4 | 53 | 6 |
| - germline RTE / STR / other | 11 / 37 / 53 | 12 / 28 / 52 | 15 / 41 / 56 | 15 / 31 / 53 |

Seed-7 ablation (comb TP / unexplained):

| configuration | comb TP | unexplained |
|---|---|---|
| keys off | 135 | 242 |
| merge only | 135 | 244 |
| merge + slippage reject | 133 | 122 |
| merge + far-pair strictness | 145 | 294 |
| all four keys | 143 | 97 |

Far-pair splitting without the slippage reject turns slipped poly-A junctions into one-sided
calls (slippage 71 -> 159), so the two keys belong together.

TP calls lost to the slippage reject, 3 over both seeds (seed 11 net per-type changes include
other moves). Each is a "TP" call that pairs a real junction with a slipped tract end, or an
insertion INTO a 20-22 bp tract where every clip is poly-A + junk:

* seed 7 L1_TSD_DELETION #32 (`30845878-30845898`): TP LEFT + slipped RIGHT at a 20T tract.
* seed 7 TEMPLATED_LOCAL #159: the exact call dies in the clipped-remap filter, as before.
* seed 11 ALU_YA5 #59: TSD = a 22A tract, with no element on either clip.

No POLYA_ONLY or ORPHAN_TD3P was lost, because none sits at a reference tract. The gains are
far-pair poly-A junctions kept as one-sided loci, plus loci that the poly-A-aware agreement no
longer drops. The extra one-sided TP calls are one-sided because the complex junction had < 2
fragments; before, the whole locus was dropped.

What is left:

* ~90 % germline assembly differences. This is acceptable: annotate classes most of them as
  `non_RTE_SV` / `microsatellite` / `unknown`, and the 11-15 hs1-only Alus are real insertions.
* 4-6 slippage loci whose other junction carries reference Alu sequence.
* 1-6 far pairs whose TP poly-A junction pairs with a germline Alu junction in the same colonies.

Byte identity:

* Legacy `config.py.grch38` without sidecars: combined/genotyping md5-identical to
  `tprt-hallmarks` HEAD.
* `.tprt` with the four new keys off: all four combine outputs identical to the integration run.

Tests: `test/test_tprt_combine_filters.py` (14).

Commands. Seed 11 used `SEED=11 OUT=$SP/work/e2e_i3_s11`. Ablations use `CI_OVERRIDES`, new in
`run_e2e.sh` / `make_config.py --override`:

```bash
cd /Users/jeremy/Documents/PEAR_TREE/.claude/worktrees/agent-a1f8bdc6b79981c09
OUT=/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/fde0700f-e325-4651-8daf-0cdd52bd072b/scratchpad/work/e2e_i3 bash test/e2e/run_e2e.sh
OUT=/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/fde0700f-e325-4651-8daf-0cdd52bd072b/scratchpad/work/e2e_i3 CI_OVERRIDES="merge_tolerance_bp=0 polya_aware_clip_agreement=False slippage_reject=False far_pair_strict=False" bash test/e2e/run_e2e.sh
```
