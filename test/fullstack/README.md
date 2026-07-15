# Full-stack test case — hs1 donor → reads → bwa-mem → GRCh38 BAM

A realistic end-to-end test: simulated short reads carrying implanted retrotransposon
insertions (true positives) **plus** deliberate assembly-discordance false positives,
mapped by the real aligner. Unlike the synthetic clip-signal harness in `../val1`, this
exercises the whole stack — real bwa-mem MAPQ, real soft-clipping, real repeat mismapping —
so it tests PEAR-TREE's *combine*-step genome-aware filtering, not just discovery.

## Design

Reads are simulated from a **donor built from hs1 (T2T-CHM13v2.0)** and mapped to
**GRCh38 (hg38)**. This cross-assembly mapping is deliberate: every hs1↔hg38 difference
(T2T centromeres/telomeres, resolved segdups, novel sequence) becomes a mapping-driven
false positive the pipeline must reject — exactly the specificity challenge real data poses.

- **True positives:** real young elements extracted from hs1 — L1HS (chrX:11.29 Mb, div 0.3),
  AluYa5 (chrX:2.76 Mb, div 0.0), HERVK (chr10:5.04 Mb, div 0.4), and **SVA_F** (chr18, div 1.1)
  / **SVA_E** (chr3, div 2.0), the youngest active SVA subfamilies (SVA_A–D also registered in
  `ELEMENT_LOCI` as positive controls) — implanted into two euchromatin windows
  (chr21:20–26 Mb, chr22:25–31 Mb) at **L1 endonuclease TTAAAA motif sites** (preferential
  retrotransposition), flanked by a target-site duplication (TSD), full-length and as variants
  (5′ truncation, 5′ inversion / twin priming, 3′ transduction; HERVK as a no-poly-A LTR
  provirus; SVA full + 5′-truncated with poly-A/TSD). Clips therefore carry genuine
  Alu/L1/HERV-K/SVA sequence. SVA loci come from the hs1 (= chm13v2.0) RepeatMasker track
  (`hs1.repeatMasker.out.gz`, local equivalent of UCSC `chm13v2.0_rmsk.bb`).
- **False-positive sources:** three hs1 FP windows included with no implants — acrocentric
  p-arm satellite/rDNA (chr21:3–6 Mb), pericentromeric alpha-satellite (chr21:10–13 Mb),
  and a telomeric end (chr22:50.8 Mb) — plus the assembly discordance everywhere.
- **Library artefacts** injected into the reads (read-intrinsic, so mapping can't create
  them): adapter read-through, poly-G/dark-cycle, fold-back palindrome chimeras, PCR
  duplicates (`inject_artefacts.py`; marked reads keep an `:art_<kind>` name suffix).

## Files

| script | does |
|---|---|
| `build_donor.py` | extracts hs1 regions + real elements, finds TTAAAA sites, implants MEIs+TSD+poly-A, writes `donor.fa`, `truth_hs1.tsv`, `flanks.fa` |
| `inject_artefacts.py` | adds adapter / poly-G / chimera / duplicate reads to the wgsim FASTQ |
| `lift_truth.py` | maps the hs1 flanks to hg38 to express the truth in **hg38 coordinates** (`truth_hg38.tsv`); implants whose hs1 flanks have no clean hg38 homolog are marked `unmappable` |
| `run_fullstack.sh` | orchestrates: donor → wgsim (~30×) → artefacts → `bwa-mem` hg38 → fixmate/sort/markdup → index → lift truth |

Outputs (written to `$OUT`, default `~/Downloads/fullstack/`): `full_stack.bam` (+`.bai`),
`truth_hg38.tsv` (score against `status==scoreable` rows), `truth_hs1.tsv`, `donor.fa`.

## Run

```bash
# needs ~/Downloads/hs1.fa (donor source) and a bwa-indexed ~/Downloads/hg38.fa
bash test/fullstack/run_fullstack.sh          # env: HS1= HG38= OUT= THREADS= DEPTH=
```

## Measured (seed 1, 30×, 34 implants)

- BAM: 3.9 M reads, 105 k duplicates marked.
- Truth: **34 scoreable** in hg38, all placed by chain-lifting the insertion point
  (`lift_truth.py --chain`). The earlier flank-mapping lift wrongly marked 4 implants
  "unmappable" — their 300 bp flanks fall in LINE/segdup-dense hs1 regions and don't bwa-map,
  though the breakpoint lifts cleanly via the chain.
- Discovery (`peartree-discovery --step discover`): **30/34 true positives recovered**,
  **65 assembly-discordance false positives**. FPs cluster in the pericentromere (chr21
  ~10 Mb), telomere (chr22 ~50 Mb) and the coordinate-shifted hg38 homolog of the implanted
  windows — exactly where hs1 and hg38 disagree. The 4 discovery misses are honest
  false-negatives: id 32 (an L1 3′-transduction recoverable with `mate_anchor_rescue`) and
  ids 6/7/8 (full L1/Alu in the LINE-dense chr21:22.1–22.2 Mb homolog).
- **combine (stage 2)** — end-to-end remap to T2T hs1, with the clean-alignment Filter A
  (`clean_remap_*`): removes **65/65** assembly-discordance FPs and keeps every real
  insertion → **0 genuine false positives**, recall unchanged at 30/34. (id 31, a real L1
  3′-transduction, is correctly retained; only its consensus can't map end-to-end cleanly.)

> `truth_hg38.tsv` junctions are ±TSD-precise; score with a small window. The `lift_method`
> column records how each row was placed (chain vs flank). SVA_E/F were added after this run
> (build_donor now implants 42); rerun `run_fullstack.sh` to fold them in.
