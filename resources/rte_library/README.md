# RTE reference library (`resources/rte_library/`)

Sequence-based references for the TPRT-hallmark annotator (see
`plans/tprt_hallmarks/SPEC.md`, section "Reference libraries", and
`docs/transduction_sources.html`). All remapping in PEAR-TREE is against **hs1 (T2T-CHM13
v2.0)**, so every element carries hs1 coordinates (primary) and hg38 coordinates
(documentation); published GRCh37 coordinates are kept verbatim in `hg19_published`.

Everything here is produced by one command (inputs: `tools/rte_library/fetch_inputs.sh`):

```
bash tools/rte_library/fetch_inputs.sh /absolute/path/rte_inputs
python tools/rte_library/build.py --inputs /absolute/path/rte_inputs \
    --work /absolute/path/rte_work --out resources/rte_library
```

Needs Python with `py2bit pyliftover edlib pysam openpyxl numpy` and `mafft` on
PATH (the MSA step). Runtime ~1 min after the RepeatMasker caches exist. Deterministic for
fixed inputs (`manifest.tsv` holds md5s).

## Files

| file | records | content |
|---|---|---|
| `l1_intact.fa` / `.tsv` | 146 | L1Base 2 `hsflil1_8438` (full-length, both ORFs intact, GRCh38), trimmed to the element, **sense** |
| `alu_y_intact.fa` / `.tsv` | 1,802 | near-full-length young AluY* copies (rmsk hg38), youngest ≤150 per subfamily, sense |
| `sva_intact.fa` / `.tsv` | 360 | near-full-length SVA_A–F (rmsk hg38, fragments chained), youngest ≤60 per subfamily, sense |
| `consensus.fa` | 9 | majority-rule consensus built from the intact sets: `L1HS L1PA2 L1PA3 ALU_Y ALU_YA5 ALU_YB8 SVA_D SVA_E SVA_F` (sense, poly-A stripped) |
| `consensus_landmarks.tsv` | 52 | consensus coordinates (1-based, inclusive) of 5'UTR/ORF1/ORF2/EN/RT/3'UTR/poly-A signal/Ta sites (L1), A box/B box/monomers (Alu), hexamer/SINE-R/pA (SVA) |
| `consensus_crosscheck.tsv` | 14 | identity of each consensus to Dfam and to L1.3 / L1.2 / L1RP |
| `dfam_young.fa` | 29 | Dfam consensus sequences of young human Alu/L1/SVA families (CC0) — cross-check / fallback only |
| `active.tsv` | 178 | L1s regarded as active: every L1HS-class FLI-L1 (Ta / pre-Ta) + every published source with ≥5 daughters, with identity to the L1HS consensus |
| `transduction_sources.tsv` | 859 | 3' transduction sources (L1 and SVA), one row per locus |
| `flanks_3p.fa.gz` (+`.fai`, `.gzi`) | 881 | downstream flank of every source in **element sense**, soft-masked (hs1 RepeatMasker); bgzip, `pysam.FastaFile`/`mappy` read it directly |
| `flanks_5p_sva.fa.gz` (+`.fai`, `.gzi`) | 257 | 5 kb upstream of every SVA source, element sense (ends at the SVA 5' end) |
| `transduction_stats.tsv` | 6 | transduction length / distal-end distance distribution from the published catalogues |
| `manifest.tsv` | — | file sizes, record counts, md5 |

Total ≈ 6.8 MB. Sources: 371 reference L1 (hg38), 198 non-reference L1, 33 hs1-only L1HS,
257 SVA_E/F; 22 non-reference sources have no published strand, so both candidate flanks
are shipped (`<id>/+`, `<id>/-`). (SPEC lists `flanks_3p.fa` / `flanks_5p_sva.fa`; they are shipped
bgzip-compressed as `.fa.gz` with faidx indices to stay small — SPEC updated.)

### Coordinate conventions
* `*_start` / `*_end` are **1-based inclusive** for elements.
* `hs1_status=insertion_point` (non-reference source, or a reference-hg38 element absent from
  CHM13): `hs1_start = hs1_end =` the **0-based junction offset** — the 3' flank begins at that
  offset (`+`) or ends there (`-`).
* `strand` is the element's own orientation; a flank is always reported in element sense, so
  read position 1 of a `flanks_3p` record is the first base downstream of the source's 3' end.
* FASTA description of every flank: `hs1:chrN:start-end(strand)` (1-based) — map a hit at
  flank offset *k* back to the genome with it.

## Provenance and verification

**L1Base.** The exported FASTA covers the BED interval *with ~1 kb of flank on each side* and
L1Base BED starts are 1-based. Orientation was **verified**: 146/146 exported sequences equal
hg38[start−1, end) reverse-complemented for `-` elements, i.e. the export is already in element
sense. Each sequence is trimmed to the L1 by an infix alignment of L1.3 (poly-A stripped); the
genomic A-rich tail is recorded as `polya_tail`. FLI-L1 means *both ORFs intact*, not
full-length 5'UTR: 5 copies start >50 bp into the consensus (`l1hs_cons_start`, up to 776). Subfamily = the RepeatMasker (hg38) L1 record
with the largest overlap. "Intact ORFs" ≠ young: 35 of the 146 are L1PA2 by RepeatMasker; the
`young`/`subfamily_call` columns say which are L1HS.

**Ta / pre-Ta** is typed from two 3'UTR sites (L1.3 numbering): 5931 (`ACA` Ta vs `ACG`
pre-Ta — the Boissinot et al. 2000 "ACA at 5930–5932" diagnostic) and 5712 (A vs G). Both
separate the long-read-typed reference sources of Nam et al. 2023 Supp. Table 4 (50/50 Ta-0/
Ta-1 carry A/A; 28/29 pre-Ta carry G/G). Result for the 146: 78 L1HS-Ta, 32 L1HS-pre-Ta,
34 L1PA2, 2 ambiguous. Ta-0 vs Ta-1 (and the Ta-1d 5'UTR deletion) is **not** typed here; the
long-read subfamily from Nam 2023 is shown in brackets in `transduction_sources.subfamily`.

**AluY thresholds.** RepeatMasker record with consensus coverage ≥280 bp of the ~311 bp
consensus (≥90 %), ≤15 consensus bases missing at the 5' end (keeps the A box at 5–15 intact,
which Pol III transcription needs), genomic span 270–360 bp and milliDiv ≤30 ‰ (≤3 %
divergence from the subfamily consensus — a copy that old has lost most CpGs/promoter
fidelity, while the youngest 150 copies of the active AluY / AluYa5 / AluYb8 sit at 0–6 ‰, so
the cap only bites for the older/rarer Y subfamilies). 25 AluY* subfamilies; 1,812 AluY, 3,134
AluYa5, 2,138 AluYb8 copies pass; the youngest (lowest milliDiv) 150 per subfamily are shipped. Consensus: youngest 120 of AluY,
AluYa5, AluYb8 → identity to Dfam 99.6 %, 100 %, 100 %.

**SVA thresholds.** RepeatMasker fragments of one subfamily chained across ≤150 bp gaps;
≥1,300 bp of the ~1,380 bp (VNTR-collapsed) Dfam consensus covered and ≤60 consensus bases
missing at the 5' end (the Dfam 5' end is the (CCCTCT)n hexamer, often annotated as simple
repeat). Youngest 60 per subfamily shipped. Consensus identity to Dfam is 65–74 % overall but
98.8–99.9 % ungapped (the difference is the VNTR, collapsed in Dfam, ~400–700 bp expanded in
our majority consensus).

**L1 consensus.** L1HS = 40 L1Base L1HS-Ta copies (mafft --auto, gap-majority columns dropped,
so minority insertions/homopolymer slippage do not enter). L1HS vs L1.3 99.68 %, L1RP 99.80 %,
L1.2 99.65 %; vs Dfam L1HS_5end 98.7 %, L1HS_3end 99.7 %. L1PA2/L1PA3 = youngest 26/40
5'-complete full-length RepeatMasker copies (Dfam *_3end 100 % / 99.9 %). L1PA2's consensus
has no intact ORF1, so its ORF1 span is transferred from L1HS (noted in the landmark row).
ORF2_EN / ORF2_RT spans are approximate (literature amino-acid ranges).

**Transduction sources** (`transduction_sources.tsv`), merged per locus:

| list | obtained? | used as |
|---|---|---|
| Nam et al. 2023 *Nature* 617:540, Supp. Table 4 (276 sources, GRCh37) | **yes** (Springer static content) | sources + daughter counts (sum over clones) + long-read subfamily/truncations |
| Rodriguez-Martin et al. 2020 *Nat Genet* 52:306, Supp. Table 5 (124, hg19) | **yes** (Europe PMC) | sources + daughter counts; strand from Supp. Table 2 transduction geometry |
| Gardner et al. 2017 *Genome Res* 27:1916 (MELT), Table S9 B (38) and C (literature compendium) | **yes** (Europe PMC) | sources, offspring counts, Tubio 2014 counts, cell-culture activity (Brouha/Beck), literature membership |
| Tubio et al. 2014 *Science* 345:1251343, Table S5 | **no** — not open access; PMC serves the supplement only behind a browser proof-of-work check | Tubio counts enter via Gardner S9 B/C |
| Brouha et al. 2003 *PNAS* 100:5280 hot L1s | **no** coordinates in the open text | their activities enter via Gardner S9 B (`cell_culture_activity`) |
| Damert et al. 2009 *Genome Res* (SVA 5' transduction groups) | no supplement with coordinates | SVA sources are seeded from RepeatMasker instead |
| seeds: every hg38 L1HS ≥5.9 kb; L1Base-intact L1PA2/3; hs1-only L1HS ≥5.9 kb; near-full-length SVA_E/F | computed | `seed` column |

All named hot loci are present with their published counts, e.g. 22q12.1 (TTC28; 922 daughters
across studies), 14q23.1, 6p24.1, 6p22.1, Xp22.2 (two sources), 2q24.1 (LRE3), 9q32, 3q21.1,
1p12, 7p12.3, 12p13.32, 1p22.3. Reference status is decided against hg38 RepeatMasker
(young L1 ≥4 kb within 300 bp); published strands agree with RepeatMasker for 97/97 Nam and
43/44 Rodriguez-Martin sources. One Nam source (X:140515257, hg19) does not lift to hg38 and is
absent. `hotness`: hot ≥20 daughters, strong 5–19, active 1–4, none_reported (published, 0),
candidate (seed only). Counts are summed across independent cohorts.

Liftover: hg19→hg38→hs1 with the UCSC chains (pyliftover); only same-chromosome hits are
accepted and an element's lift must be bracketed by its ±1 kb anchors (the over.chain files
contain small paralog chains — e.g. a chr6 SVA lifting onto a chr1 SVA in hs1). Elements absent
from CHM13 get `hs1_status=insertion_point` at the lifted 3' junction (their downstream flank
is still reference sequence, which is what transduces).

**Flank length.** 15 kb for L1 (Tubio 2014: transductions reach up to ~12 kb downstream). From
the published catalogues (`transduction_stats.tsv`, n = 3,676): median transduced length 296 bp,
95 % ≤ 991 bp; distal end median 666 bp downstream, 99.2 % ≤ 10 kb, **99.8 % ≤ 15 kb**. SVA:
5 kb 3' and 5 kb 5' flanks. Validation: 383/400 randomly drawn Rodriguez-Martin transduced
segments map to `flanks_3p.fa.gz` in sense orientation (1 antisense, 16 no hit; median offset
316 bp, 95th percentile 3.7 kb into the flank).

`pas_hexamers_3p`: first five AATAAA / ATTAAA positions in each flank (1-based, element sense)
— transductions from a given source terminate at its downstream pA sites (Zumalave et al.
2024, bioRxiv: all 83 transductions of the 2q24.1 source end 234 bp downstream, at an alternative pA
site, the source lacking the canonical L1 pA signal); `canonical_pas_3p` says whether the
source keeps the canonical AATAAA at its 3' end.

## Licences and citation
* L1Base 2 — Penzkofer T. et al. (2017) *Nucleic Acids Res* 45:D68–D73. Cite when using
  `l1_intact.*`.
* Dfam — consensus sequences are CC0 (Storer J. et al. 2021 *Mob DNA* 12:2).
* RepeatMasker tracks, genomes, chains — UCSC Genome Browser downloads (hg38, hs1 /
  T2T-CHM13v2.0, Nurk S. et al. 2022 *Science* 376:44).
* Source tables — Nam C.H. et al. 2023 *Nature*; Rodriguez-Martin B. et al. 2020 *Nat Genet*;
  Gardner E.J. et al. 2017 *Genome Res*; Tubio J.M.C. et al. 2014 *Science* (via Gardner).
  These are facts (coordinates/counts) re-derived and lifted, not redistributed tables.
* GenBank L19088 (L1.3; Dombroski et al. 1993), M80343 (L1.2), AF148856 (L1RP).

## Extending
Accepted novel sources are appended with `tools/rte_library/add_source.py` (tier A: ≥5.5 kb and
identity ≥0.98 to the L1HS consensus; tier B: 0.95–0.98 with ≥2 independent daughters; criteria
and procedure: `docs/transduction_sources.html`, rendered from this directory by
`tools/rte_library/make_doc.py`). Re-running `build.py` regenerates the library
from scratch and therefore drops appended sources — re-apply them or add them to a seed list.
