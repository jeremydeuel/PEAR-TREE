# Spencer Chapman et al. 2025 — the aggregated tree collection

**"Prolonged persistence of mutagenic DNA lesions in somatic cells", Nature, 2025-01-15,
doi 10.1038/s41586-024-08423-8.** 104 Newick trees, 11,496 tips, 89 donors, 568 KB.

Source: `github.com/mspencerchapman/Prolonged_persistence_of_DNA_lesions`,
`Data/input_data/<COHORT>/`. Fetched 2026-07-17.

## THIS IS WHERE OUR TREES CAME FROM

Verified 2026-07-17 by exact tip-count match against `donor_table.tsv` on 16 donors:

    PD45534=922  PD49236=422  PD40521=407  PD40667=380  PD45517=379  PD48402=367
    PD41048=362  PD43974=328  PD47738=315  PD43947=277  PD47703=259  PD41768=234
    PD49237=231  PD40315=195  PD44579=174  ...

That settles what looked like an anomaly. This paper **re-analysed and re-built** the source
cohorts' trees, and its trees are BIGGER than the originals:

    donor     original paper          Chapman 2025 (= ours)
    PD45534   451  (Mitchell 2022)          922
    PD49236   173  (Mitchell 2025)          422
    PD49237    44  (Mitchell 2025)          231
    PD45517   390  (Mitchell 2022)          379   (smaller)
    PD40315   216  (Mitchell 2022)          195   (smaller)

So "our PD49236 tree has 422 tips but the paper says 173" was never a bad ingest — we had the
2025 re-analysis and were comparing it to the 2025 *chemo* paper's subset. **When a tip count
disagrees with a source paper, check this repo before assuming a bug.**

## Layout — one directory per source cohort

    EM/         14   Mitchell 2022 ageing (KX/CB/SX/AX codenames)
    KY/         16   Yoshida 2020 bronchial epithelium
    MF/          3   Fabre 2022 clonal haematopoiesis
    MSC_BMT/    10   transplant donor-recipient PAIRS -- each tree is TWO people
    MSC_fetal/   2   Spencer Chapman 2021 foetal (8 pcw, 18 pcw)
    NW/         11   Williams 2022 MPN
    SN/         48   liver LCM (Brunner 2019 + Ng 2021)

## Parsing gotchas — both of these will silently give you zero

1. **These newicks end in an INTERNAL NODE LABEL**, e.g. `)Cl.PD;` — not `);`. A parser that
   requires a trailing `);` returns 0 tips for MF and SN without erroring.
2. **SN (liver) tip labels are `Cl.N` clone ids, NOT sample names.** The PD id exists only in
   the FILENAME (`tree_PD36713b.tree`). These 48 trees are the `Cl.N-NEEDS-KEY` rows in
   `donor_trees.tsv`. Attributing them means trusting the filename — the exact thing that bit
   us in commit 3f21262 and again on PD48367b..h. There is no clone->sample key in the repo.
3. `MSC_BMT/` Pair trees each hold **two individuals** (donor + recipient), e.g. Pair11 =
   PD45792 (179 tips) + PD45793 (146). Do not treat a Pair as one donor.
4. `SN/` includes `tree_PD48367b..i` and `tree_PD48372b..i` — 8 per-region liver trees for
   TWO people, not 16 donors. See donor_clinical.tsv.

## Donor manifest

`41586_2024_8423_MOESM4_ESM.xlsx`, sheet `Individual_metadata` — 103 phylogeny rows / 89
unique individuals. Cohort sizes: Liver 48, Lung 16, Ageing 12, MPN 10, Transplant 10, CH 3,
Chemo 2, Foetal 2. Every donor is RE-USED from an earlier paper except the two chemo donors
(PX001=PD47703, PX002=PD44579), whose citation cell literally reads
`TO ADD WHEN HAVE PREPRINT REFERENCE`.

Validation: tip counts match the manifest's independent "Number of samples" column on
**103/103 rows**.

## An undocumented 104th tree

`NW/tree_PD5147.tree` (67 tips) is in the repo but **PD5147 appears nowhere in the manifest or
the full text**. It is a Williams 2022 MPN donor (81F, PV). Present in the data, absent from
the paper's own donor list.

## The six unplaceable donors are NOT here

PD49226, PD49227, PD49228, PD49229, PD53373, PD53374 — **zero hits** across the whole 311 MB
repo and all supplements. Checked explicitly. They remain unpublished. Beware the near-miss:
this repo DOES contain PD49236 (KX009) and PD49237 (KX010). Adjacent numbers, different people.

## Why the paper matters beyond the trees

818 DNA lesions that persist ACROSS CELL CYCLES (~2.2 years mean, 15-25% >=3 years, ~8 per
cell). A persistent lesion is replicated repeatedly, so daughter lineages independently
misincorporate opposite it, yielding **PVVs -- phylogeny-violating variants**: the same
mutation in lineages the consensus tree says are not sisters.

**Relevance to us: low for MEI calling, but it punctures an axiom.** It is strictly a
single-base-substitution mechanism -- no soft-clips, no discordant pairs, no insertion
signature -- so it cannot manufacture an MEI artefact. BUT our filtering treats
non-monophyletic support as an FP signature (see [[carrier-count-is-binomial-tail]], the
0/1395 monophyletic result). This paper is published proof that **real biology can violate the
consensus phylogeny**. It is far too rare (<1 in 1e5 lesions yields a detectable PVV) to
explain our FP bands, and it does not apply to insertions, so the monophyly gate stands. But
"phylogeny-violating => artefact" is not an axiom, and this is the paper that says so.
