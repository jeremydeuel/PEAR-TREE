# 5. Common sequencing artefacts (focus question 3)

*Sources: Chen et al. BMC Genomics 2024;25:227 ("Characterization and mitigation of
artifacts derived from NGS library preparation due to structure-specific sequences");
Tanaka et al. PLOS ONE 2020;15(1):e0227427 ("Sequencing artifacts derived from a library
preparation method using enzymatic fragmentation"); plus artefact observations from the
retrotransposition literature ([§3](03_detecting_true_events.md)) and the Illumina
two-colour chemistry. This section catalogues the artefacts; [§4](04_artefacts_in_discovery.md)
turns each into a discovery-time filter.*

---

The artefacts below are ordered by how strongly they mimic a retrotransposition/SV
breakpoint, because a clipped-read caller like PEAR-TREE is most endangered by anything
that produces **soft-clipped reads, split reads, or discordant pairs**.

> For **measured prevalence of each of these classes across all 560 real PEAR-TREE discovery
> files** (5.8M clipped sequences, with concrete example records), see
> [§10](10_empirical_discovery_patterns.md). In that data the heaviest contributors are
> homopolymer/low-complexity clips (§5.6; ~12–14%) and recurrent repeat-element/mapping clips
> (§5.9; ~45% of clips are non-unique, Alu-dominated), with structure-specific chimeras
> (§5.1; 0.65%) and poly‑G (§5.4; 0.66%) forming the tail.

## 5.1 Structure-specific chimeras (the dominant SV mimic)

This is the central finding of both artefact papers and the artefact most dangerous to a
clipped-read RTE caller.

**Mechanism (Chen's PDSM model — "Pairing of partial single strands derived from a similar
molecule").** The human genome is dense with **inverted repeats (IVRs)** and
**palindromes/hairpins**. During library prep:

1. the duplex template is cleaved (randomly by sonication, or **at specific palindromic
   sites** by the fragmentation endonuclease);
2. a partial single strand folds back on itself, or reverse-complement-pairs with the
   matching half of the same IVR/palindrome carried on a *different* partial molecule;
3. the 3′ overhang is trimmed and the gap is **filled in by polymerase during
   end-repair / A-tailing** — this fill-in is where the erroneous bases are actually
   incorporated (the endonuclease itself cannot mutate DNA);
4. PCR amplifies the chimeric molecule.

**Read-level signature.** A **soft-clipped read whose clipped segment, re-mapped locally,
aligns back nearby in reverse-complement (inverted) orientation.** At first pass this is
*indistinguishable* from a foldback-inversion junction, a small inversion, or the clipped
signature of an insertion breakpoint. Key discriminating features:

- The clipped sequence is itself **palindromic**; the variant/junction sits at the
  **centre of a palindrome** ("SNV-centred palindrome," SCP; usually odd-length ≥ 5 bp).
- **Heavy soft-clipping**: artefact reads had a mean soft-clip ratio of **50.8%** vs
  **5.0%** for genuine reads (Tanaka Fig. 3D).
- **Positional bias**: the artefactual base/junction clusters **10–15 bp from the read
  end**; in **90.4%** of SCPs the whole palindrome lies within **30 bp of the read edge**.
- **Recurrence** at the *same* coordinate across many samples (because enzymatic cutting is
  site-specific).
- Chen notes that when the reverse-complement pairing is between **more distant** loci, the
  chimeric read instead mimics a **fusion / translocation junction** — directly relevant to
  transduction/translocation calling ([§3](03_detecting_true_events.md)).

**Prevalence.** 371–655 (median 568) distinct artefact palindromes per sample (Tanaka);
> 1,000 per sample (Chen). Category-[a] chimeric artefacts: 11,731 vs 2,984 genuine
variants in Tanaka's data. They concentrate at natural inverted-repeat/palindromic loci in
the reference — the *same* structure-rich regions where SV/MEI callers are most prone to
false positives.

## 5.2 Enzymatic vs mechanical fragmentation bias

Enzymatic fragmentation produces **substantially more** structure-specific chimeras than
mechanical shearing:

- Chen (54 tumours): median **115 (26–278)** artefactual variants enzymatic vs **61
  (6–187)** sonication; 5,544 enzymatic-only vs 2,599 sonication-only vs 682 shared.
- Tanaka (HyperPlus vs SureSelect): enzymatic gave **2.3–9.9× more** calls (median 2,308
  SNVs, 89 indels).
- **Mechanistic reason:** sonication cuts *randomly* → longer single strands, and any
  fill-in error is spread across positions so it never accumulates at one coordinate. The
  endonuclease cuts *specific palindromic sites* → shorter single strands and **recurrent
  centre-of-palindrome errors**. Controls confirm causality: a same-vendor *ultrasonication*
  kit produced no excess; using the enzymatic method for both tumour and normal did **not**
  remove the recurrent palindrome artefacts.

**Implication for PEAR-TREE:** know your fragmentation chemistry. Enzymatically fragmented
WGS (increasingly common for low-input/automated prep) will hand the discovery step many
more inverted-repeat chimeras, making the same-contig-supplementary/`--local` map-back
checks ([§7](07_peartree_code_review.md)) essential rather than optional.

## 5.3 Poly-A dropout (the RTE-specific library artefact)

The most consequential artefact for *retrotransposition* calling specifically. Low-input /
enzymatic library preps (notably LCM-based WGS) **deplete poly‑A-carrying reads**:

- soL1Rs from LCM had **0.11 vs 0.33** poly‑A reads per adjusted depth and a poly‑A/L1 read
  ratio of **0.26 vs 1.07** (P = 2.2×10⁻⁷²; Nature 2023, Suppl. Figs 2–3).
- Because the poly‑A tail is the imprint of reverse transcription and the primary
  supporting evidence for an L1 insertion, its loss **cripples MELT, TraFiC-mem, DELLY and
  xTea** on such data.

This is why the Nature 2023 authors argue LCM/low-input WGS is unsuitable for soL1R
calling, and why PEAR-TREE's poly‑A-rescue logic is only as good as the poly‑A content the
library retains.

## 5.4 Two-colour chemistry poly-G / dark-cycle artefacts

On Illumina two-channel systems (NextSeq, NovaSeq, iSeq), **"no signal" is read as G**.
When a cluster's template is exhausted or synthesis fails, the tail of the read becomes a
run of **high-quality-looking G's**. These produce spurious G-homopolymer clipped tails.
PEAR-TREE explicitly trims 3′ poly‑G in `clean_clipped_seq` — a necessary defence, since a
poly‑G clip could otherwise masquerade as inserted sequence.

## 5.5 Adapter read-through

When the insert is shorter than the read length, sequencing runs into the 3′ adapter. The
adapter appears as a clipped tail. Uncorrected, it is a chimeric-looking junction. PEAR-TREE
carries explicit NebNext adapter sequences (and common 1-bp-error variants) and clips them
(`is_adapter`, `clean_clipped_seq`). Adapter tables **must** match the actual library kit
and platform.

## 5.6 Homopolymer / microsatellite slippage

Polymerase slippage in homopolymer runs and short tandem repeats creates length-variable
indels and mismatched tails, especially near read ends. These generate low-complexity
clipped sequences and imprecise breakpoints. Defences: minimum-complexity checks, rejecting
breakpoints whose flank is a short-period repeat, and homopolymer-length caps (PEAR-TREE's
`max_homopolymer_len`, n-polymer unclipped filter, `is_low_complexity`).

## 5.7 PCR duplicates and PCR chimeras

- **Duplicates** inflate apparent read support and can make a single artefactual molecule
  look like independent evidence. Mark/remove with Picard/samblaster **before** discovery;
  PEAR-TREE skips `is_duplicate` reads, so upstream duplicate marking is a prerequisite.
- **PCR chimeras / template switching** join unrelated molecules, producing false junctions.
  Mitigated by fewer PCR cycles (PCR-free preps), duplicate removal, and requiring
  reproducibility across independent samples.

## 5.8 Index hopping (patterned flow cells)

On ExAmp/patterned flow cells (NovaSeq), free adapters cause **reads to be assigned to the
wrong sample**. In a multi-sample clonal study this can make a real insertion from sample A
appear at low level in sample B — directly corrupting the population/phylogenetic filters
that PEAR-TREE relies on. Mitigations: **unique dual indexes (UDI)**, and treating very
low-fraction cross-sample calls with suspicion (PEAR-TREE's `max_artefact`/VAF-style score
thresholds partly absorb this).

## 5.9 Mapping artefacts (not library, but indistinguishable downstream)

- **Reference gaps / mis-assembly / segmental duplications:** reads spanning a
  mis-assembled or collapsed-duplication boundary soft-clip and look like a novel junction.
- **Mismapping near existing repeats:** a read from one repeat copy mapped to another
  clips at the divergence point. This is why the retrotransposition papers require **TSD +
  poly‑A** and matched-normal/panel support before trusting a clipped junction near a
  repeat, and why PEAR-TREE re-maps the clipped part to check it does not land locally.
- **Oxidative damage (8-oxoG)** from mechanical shearing gives the classic C:G>A:T /
  C:G>G:C low-VAF SNV artefacts — distinct from the palindrome chimeras, mitigated by
  antioxidants and by low-VAF filtering.

## 5.10 Summary table

| # | Artefact | Molecular cause | Read-level signature | Mimics |
|---|---|---|---|---|
| 5.1 | Structure-specific chimera | IVR/palindrome fold-back + fill-in error (PDSM) | clip re-maps locally in **inverted** orientation; palindromic clip; clip ratio ~50%; ≤30 bp from read end; recurrent | inversion / insertion / (distal) fusion |
| 5.2 | Enzymatic-fragmentation bias | site-specific endonuclease cutting | as 5.1 but **recurrent at fixed coordinate** | as 5.1, amplified |
| 5.3 | Poly‑A dropout | low-input/enzymatic prep loses poly‑A reads | **missing** poly‑A evidence | *false negatives* for L1 |
| 5.4 | Poly‑G / dark cycle | 2-colour "no signal"=G | high-Q 3′ poly‑G clip | inserted sequence |
| 5.5 | Adapter read-through | insert < read length | adapter seq as clipped tail | insertion junction |
| 5.6 | Homopolymer/STR slippage | polymerase slippage | low-complexity clip, imprecise bp | small indel / insertion |
| 5.7 | PCR duplicate / chimera | amplification | inflated support / false junction | any (spurious support) |
| 5.8 | Index hopping | free adapter on patterned cell | low-level cross-sample calls | corrupts population filter |
| 5.9 | Mapping / assembly | segdup, mis-assembly, repeat mismap | clip at reference discontinuity | insertion / SV |
