# Memory-budgeted storage of PENDING disc_agree junctions (`mem_budget_mb`)

Branch `discovery-mem-budget` (off `discovery-mem-profile` b5a333c; see REPORT.md for the
diagnosis this answers). Crate `rust/peartree-discovery`. Commits: 1c192a2 (the feature),
fbb4b32 (default-off `cluster/pipeline.sh` wiring), plus the BUDGET.md / ord-u32 commit.

## 1. What it does

The 18.5 M single-molecule junctions that `disc_agree_second_fragment` holds between a
contig's `cleanup` and the mate pass cost ~1.75 KB each as fully built `Breakpoint`s
(consensus, mates, sidecar lite lists, `AgreePending`), i.e. ~32 GB on PD45808b_lo0003,
and 85 % of them are then rejected. With `mem_budget_mb = N` (config key; env
`PEARTREE_MEM_BUDGET_MB` wins; unset = unlimited = exactly today's behaviour):

* pending junctions are kept **FULL** while the *accounted* bytes of the FULL pending pool
  (`crate::mem::breakpoint_breakdown`, the same arithmetic as the `[mem]` lines — never RSS,
  so the switch point is a pure function of the input) stay `<= N MiB`;
* the first time a newly admitted pending junction pushes the total over `N`, **every FULL
  pending junction is converted in place to the COMPACT form** and every later pending
  junction is stored compact ("trip"); `N = 0` therefore means everything compact;
* after the mate pass and `resolve_disc_agree`, every **promoted** compact junction is
  **rebuilt from the BAM** (one forward pass over the sorted record addresses) through the
  same `clip_breakpoint` -> `join()` / lite constructors / `mate_seq_of` that the FULL path
  ran, and merged back into its slot; rejected compact junctions are dropped unread.

The representation is never observable: the acceptance test below shows byte-identical
`.txt.gz`, `.evidence.tsv.gz` and `stats.json` for every budget and thread count. `stats.json`
counters are incremented once (at the first build; the rebuild joins into a sink `Stats`).

One stderr line summarises the run, always printed when a budget is set:

```
[mem-budget] pending-junction budget 20 MiB: overflowed at 1:3457550 LEFT; 98092 junction(s) compact (10585 converted in place, 15.9 MiB arena, 169 B/junction), 14152 rebuilt from the BAM (105881 record requests, 60353 seeks, 103932 records read, 2.7 s)
```

CRAM input has no BGZF record addresses: a budget is ignored with a note.

## 2. The compact form (`src/compact.rs`)

One `Vec<Hdr>` (56 B per junction) plus one byte arena — no per-junction heap allocation.
Per junction the arena holds:

| part | bytes | why it is needed |
|------|------:|------------------|
| clip record addresses (BGZF virtual positions, u64 each), **every** record of the cluster in scan order, duplicates included | 8 x reads (typ. 1–2) | rebuild: duplicate collapsing happens inside `join()`, so the whole cluster is re-fed to it |
| DISC sidecar read addresses, in `attach_disc` list order | 8 x (typ. 3–4) | rebuild of `EvExtra::disc_lite` by re-reading exactly the selected records |
| SHORT sidecar read addresses, in `attach_short` list order | 8 x (typ. 4–5) | rebuild of `EvExtra::short_lite` |
| anchors: qname hash (u64) + read1/strand bits (u8) | 9 x (<= 16) | the mate pass: this is precisely what the FULL path's `AgreeWait` carries (it also keys on the hash alone) |
| clip inserts, 4-bit packed (`=ACMGRSVTWYHKDBN`, lossless; raw-byte fallback) | 3 + n/2 each (<= 4) | `agrees()` in the mate pass |
| `has_mate` clip qnames (raw bytes) + mate-read1 bit | 3 + len each | `get_mates`: the compact junction's qnames enter the mate-capture maps in slot order |

Measured on the synthetic stress BAMs: **169 B/junction** (stress1) and **200 B/junction**
(stress2, more anchors) = header + arena + 1 B promoted flag. The mate records the mate pass
captures for a compact junction are stored as `(junction, address)` pairs (12–16 B) instead
of a cleaned 150-bp `QualitySeq` (~350 B) and re-read at rebuild. Accounted pool size on
stress1: 98,092 pending junctions = 170.9 MiB FULL (1,797 B/bp) vs 23.1 MiB compact (incl.
Vec capacity).

Slot identity: every final breakpoint now carries `ref_id` and a per-run creation counter
`ord` (assigned in `cleanup`); `final_*_breakpoints` are in `(ref_id, ord)` order per side
(contigs in scan order, junctions in cleanup order). A compact junction leaves the vector and
is merged back by that key, so `Emitted`, `sidecar_prepare`, `output`, the rescue etc. see the
usual order. During the mate pass the index-based structures carry a `Compact(i)` variant
(`BpRef::Compact`, `WaitAt::Compact`), and `get_mates` walks FULL slots and compact entries
in merged slot order so `qmap` overwrite semantics are unchanged.

### Conversion order (deterministic)

At the trip: LEFT finals (slot order), LEFT staging of the current contig, RIGHT finals,
RIGHT staging — i.e. the order of `final_left ++ final_right` as they would have been. Then
every later pending junction in admission order. The trip is tested per junction inside
`cleanup` (each junction is now finished — anchors selected, no-anchor drop, DISC/SHORT
attached — and admitted one at a time; these steps depend only on the junction and the
contig's complete collectors, so this is the same computation as the former slice passes),
so the overflow point can and does fall mid-contig (e.g. `1:670700 RIGHT` for 50 MiB on
stress1, `2:2452100 RIGHT` for 150 MiB).

SPD-3 (`--threads > 1`): a worker admits against `budget / threads` (its own share; it holds
at most that much FULL before returning), and the parent re-admits every worker's junctions
in contig order against the full budget as it merges, converting as needed. The trip point
with `--threads > 1` is therefore a different (still deterministic) junction than with
`--threads 1`; the output is identical either way (tested). The farm runs `--threads 1`.

### Why the lite lists and the mate decision are exactly reproduced (approach (a))

* CLIP lite list: `merge_ev` inside `join()` on the re-read cluster — same records, same
  order (records come back in file order = scan order; the group is stable-sorted by
  position exactly as `cleanup` sorts the contig's reads, and restricting that sort to the
  group's members preserves their relative order).
* DISC / SHORT lite lists: the *selected* entries' addresses are stored in list order and
  each record is re-read through the same `disc_lite_of` / `short_lite_of` constructor the
  extract used, so the rebuilt `Vec` has the same elements in the same order. No window,
  cap, dedup or ordering has to be re-derived.
* Promotion: the FULL path matches an anchor's inside mate by `frag_hash(qname)` alone
  (`agree_mate_map` is keyed by `DiscLite::frag`; the anchor's qname string was never kept).
  The compact form stores that hash and the same two bits, so a 64-bit collision is
  resolved identically in both forms — tested with a simulated colliding anchor
  (`compact_anchor_map_equals_full_including_colliding_hashes`). Verifying the name on
  re-read would *change* output relative to the FULL path on a collision, so it is not done.
* Mates: `mate_seq_of` / landing site from the re-read mate record, pushed in capture order.
* `n_frags = max(n_frags, 2)` and `agree = None` as `resolve_disc_agree` leaves a promoted
  FULL junction.

The end-to-end unit test `budget_zero_rebuilds_identical_breakpoints_from_the_bam` writes a
small BAM, runs discovery with the budget unset and 0, and asserts equality of every field
(consensus seq+qual, mates, mate_seqs, mate_dests, CLIP/DISC/SHORT lite lists, n_reads,
n_frags, n_frags_span_polya, ref_id, ord) plus identical `stats.json`.

## 3. Test matrix (local, macOS, jemalloc; `/usr/bin/time -l` peak RSS)

Baseline = binary built from `discovery-mem-profile` b5a333c, budget unset, `--threads 1`.
Every row is `cmp` of the decompressed `.txt.gz`, decompressed `.evidence.tsv.gz`,
`stats.json` (and hallmarks/splice when present) against that baseline. `main` 08d51e5 was
also run on stress1 and test.bam: identical to the baseline.

stress1 = `tests/gen_pending_stress.py --contigs 2 --len 8000000 --cov 25 --junction-every 150 --seed 1`
(105,589 pending / 7,497 no-anchor / 14,152 promoted / 83,940 rejected; 4,748 txt records,
28,189 evidence lines), config `tests/config.pending_stress`.

| run | peak RSS | wall | compact / converted in place | identical |
|-----|---------:|-----:|------------------------------|-----------|
| baseline b5a333c | 244.9 MB | 11.9 s | – | – |
| main 08d51e5 | 272.1 MB | 11.9 s | – | yes |
| new, unset, t1 | 268.1 MB | 11.7 s | 0 | yes |
| new, 0 MiB, t1 | 132.5 MB | 14.2 s | 98,092 / 1 (trip 1:2150 LEFT) | yes |
| new, 5 MiB, t1 | 147.4 MB | 14.2 s | 98,092 / 2,656 (1:878150 LEFT) | yes |
| new, 20 MiB, t1 | 158.6 MB | 14.2 s | 98,092 / 10,585 (1:3457550 LEFT) | yes |
| new, 50 MiB, t1 | 183.5 MB | 14.3 s | 98,092 / 26,577 (1:670700 RIGHT, mid-contig) | yes |
| new, 100 MiB, t1 | 227.4 MB | 14.7 s | 98,092 / 54,230 (2:1684850 LEFT) | yes |
| new, 150 MiB, t1 | 259.6 MB | 14.7 s | 98,092 / 81,077 (2:2452100 RIGHT) | yes |
| new, unset, t2 | 353.1 MB | 8.0 s | 0 | yes |
| new, 0 MiB, t2 | 236.4 MB | 12.7 s | 98,092 / 2 | yes |
| new, 5 MiB, t2 | 239.8 MB | 12.0 s | 98,092 / 2,650 | yes |
| new, 20 MiB, t2 | 251.2 MB | 11.1 s | 98,092 / 10,594 | yes |
| new, 50 MiB, t2 | 255.2 MB | 14.0 s | 98,092 / 26,475 | yes |
| new, 100 MiB, t2 | 304.0 MB | 12.1 s | 98,092 / 53,135 | yes |
| new, 150 MiB, t2 | 343.5 MB | 10.5 s | 98,092 / 81,078 | yes |

Rebuild pass on stress1: 105,881 record requests, 60,353 seeks, 103,932 records read, 2.7–3.3 s
(= the whole extra wall time, +2.4 s on 11.9 s).

stress2 = `--contigs 3 --len 6000000 --cov 30 --junction-every 100 --anchor-frac 0.95 --agree-frac 0.3 --seed 7`
(denser: 118,498 pending / 912 no-anchor / 33,605 promoted / 83,981 rejected with the
`1,2` allowlist; 3,588 txt records, 26,631 evidence lines). Same config.

| run | peak RSS | wall | compact / converted | identical |
|-----|---------:|-----:|---------------------|-----------|
| baseline | 303.0 MB | 17.8 s | – | – |
| new, unset, t1 | 331.5 MB | 17.4 s | 0 | yes |
| new, 0 MiB, t1 | 191.9 MB | 23.2 s | 117,586 / 1 | yes |
| new, 5 MiB, t1 | 192.4 MB | 23.3 s | 117,586 / 2,396 | yes |
| new, 20 MiB, t1 | 198.9 MB | 23.0 s | 117,586 / 9,572 | yes |
| new, 50 MiB, t1 | 216.0 MB | 23.3 s | 117,586 / 23,898 | yes |
| new, 100 MiB, t1 | 247.5 MB | 23.4 s | 117,586 / 48,813 (1:3959200 RIGHT) | yes |
| new, 150 MiB, t1 | 298.0 MB | 23.8 s | 117,586 / 73,263 (2:2949700 LEFT) | yes |
| new, unset, t2 | 406.7 MB | 11.7 s | 0 | yes |
| new, 0 MiB, t2 | 253.9 MB | 19.3 s | 117,586 / 2 | yes |
| new, 5 MiB, t2 | 258.1 MB | 18.9 s | 117,586 / 2,399 | yes |
| new, 20 MiB, t2 | 272.6 MB | 19.1 s | 117,586 / 9,572 | yes |
| new, 50 MiB, t2 | 300.1 MB | 22.8 s | 117,586 / 23,940 | yes |
| new, 100 MiB, t2 | 338.6 MB | 19.7 s | 117,586 / 47,849 | yes |
| new, 150 MiB, t2 | 376.7 MB | 23.9 s | 117,586 / 72,445 (1:1389000 RIGHT) | yes |

Rebuild pass on stress2: 333,049 requests, 169,363 seeks, 310,243 records read, 7.2–7.7 s
(33.6 k promoted junctions; the extra wall time is +5.5 s on 17.8 s).

stress2 with all three contigs allowed (`contig_allowlist = 1,2,3`; 176,371 pending,
50,407 promoted; 5,382 txt records, 40,298 evidence lines) at 1, 2 and 3 threads
(exercises the parent merge over three worker results): see §3b.

`test_data/test.bam` (210 reads, mouse, no pending junctions) with the stress config minus
the allowlist, budget unset / 0, 1 / 2 threads, vs baseline and main: see §3b.

`cargo test --release`: 83 passed (79 before + `pack4_roundtrips_the_bam_alphabet_and_rejects_others`,
`roundtrip_every_field_and_no_per_junction_allocation`,
`compact_anchor_map_equals_full_including_colliding_hashes`,
`budget_zero_rebuilds_identical_breakpoints_from_the_bam`), 0 warnings.

### 3b. stress2 with all three contigs, and test.bam

| run | peak RSS | wall | compact / converted | identical |
|-----|---------:|-----:|---------------------|-----------|
| baseline | 427.1 MB | 21.9 s | – | – |
| new, unset, t1 | 476.7 MB | 21.8 s | 0 | yes |
| new, 0 MiB, t1 | 262.5 MB | 32.9 s | 176,371 / 1 | yes |
| new, 30 MiB, t1 | 261.6 MB | 33.6 s | 176,371 / 14,348 (1:2931900 LEFT) | yes |
| new, 120 MiB, t1 | 316.9 MB | 34.9 s | 176,371 / 58,895 (2:21700 LEFT) | yes |
| new, unset, t2 | 513.0 MB | 19.2 s | 0 | yes |
| new, 0 MiB, t2 | 257.8 MB | 28.4 s | 176,371 / 3 | yes |
| new, 30 MiB, t2 | 283.8 MB | 30.4 s | 176,371 / 21,539 | yes |
| new, 120 MiB, t2 | 355.5 MB | 28.2 s | 176,371 / 86,054 | yes |
| new, unset, t3 | 599.0 MB | 13.9 s | 0 | yes |
| new, 0 MiB, t3 | 386.6 MB | 24.6 s | 176,371 / 3 | yes |
| new, 30 MiB, t3 | 400.2 MB | 24.8 s | 176,371 / 14,359 | yes |
| new, 120 MiB, t3 | 482.4 MB | 24.7 s | 176,371 / 57,387 | yes |

Rebuild pass here: 500,118 requests, 254,271 seeks, 465,942 records read, 11.1–12.6 s for
50,407 rebuilt junctions (+11 s on 21.9 s).

`test_data/test.bam` (5 txt records, 23 evidence lines; 3 pending junctions, none promoted):
baseline, main, new unset / 0 MiB at 1 and 2 threads — all identical.

All 46 rows above were produced by the binary before the final `ord: u64 -> u32` shrink; the
confirmation matrix after it (stress1 unset/0/50 MiB x t1/t2, stress2-all-contigs 0/120 MiB
x t1/t3) is in §3c.

### 3c. Confirmation after the `ord` shrink (final binary, `Breakpoint` = 328 B)

| run | peak RSS | wall | compact / converted | identical |
|-----|---------:|-----:|---------------------|-----------|
| stress1 baseline | 244.9 MB | 11.7 s | – | – |
| stress1 unset t1 | 267.4 MB | 12.0 s | 0 | yes |
| stress1 0 MiB t1 | 145.8 MB | 14.3 s | 98,092 / 1 | yes |
| stress1 50 MiB t1 | 183.5 MB | 14.5 s | 98,092 / 26,691 (1:707900 RIGHT) | yes |
| stress1 unset t2 | 339.3 MB | 7.9 s | 0 | yes |
| stress1 0 MiB t2 | 243.5 MB | 10.5 s | 98,092 / 2 | yes |
| stress1 50 MiB t2 | 266.8 MB | 10.5 s | 98,092 / 26,584 | yes |
| stress2-all baseline | 427.1 MB | 21.5 s | – | – |
| stress2-all 0 MiB t1 | 261.5 MB | 31.8 s | 176,371 / 1 | yes |
| stress2-all 120 MiB t1 | 323.8 MB | 32.8 s | 176,371 / 59,108 (2:65300 LEFT) | yes |
| stress2-all 0 MiB t3 | 386.7 MB | 45.7 s (disk-contended rebuild, 18.6 s) | 176,371 / 3 | yes |
| stress2-all 120 MiB t3 | 486.6 MB | 24.4 s | 176,371 / 57,596 | yes |

(The 50 MiB trip point moved from 1:670700 to 1:707900 because the accounted size per FULL
junction shrank with the struct — expected; output identical.)

RSS noise check, stress1, three sequential runs each: baseline 245.0 / 245.0 / 245.0 MB,
new (unset) 267.3 / 267.4 / 262.8 MB, wall 11.8–11.9 vs 11.8–12.0 s. So the unset-mode
overhead on this input is a real +18–22 MB (+8 %), ~9 MB of which is the accounted struct
growth; see §4.

## 4. Memory: what the budget covers and what it does not

Per pending junction, FULL vs compact (accounted): 1,797 B vs ~170–200 B (9–10x). Peak RSS
on stress1 goes 268 -> 132 MB at budget 0; the residual is everything the budget does not
cover:

* `temporary_breakpoints` — every soft-clipped read of the current contig as a per-read
  `Breakpoint` (now 328 B + its sequences, ~600 B each): ~3.4 GB on chr1 of PD45808b.
  Transient per contig; chr1 comes first, when the pending pool is still small, so it does
  not add to the end-of-extract peak (which is dominated by the pool).
* the per-contig collectors: the short-overhang `kept` list (`(ShortLite, u64)` = 64 B, was
  56: +8 B per entry for the record address, ~ every read within 20 bp of a clipped read)
  and `sc_disc_tmp` (`(DiscLite, u64)` = 56 B, was 48). Transient per contig.
* ordinary (non-pending) finals: ~1 M x ~1.8 KB ≈ 1.8 GB on PD45808b — unchanged.
* `qmap` + the two qname sets in the mate pass (3 copies of ~8.9 M qnames ≈ 1.8 GB on
  PD45808b): unchanged; the compact junctions' qnames still enter them (they must, for
  identical output).
* the `agree_map` (`AgreeWait` per anchor, ~16 B + Vec overhead per entry) — unchanged.
* coverage bins, BAM decode buffers, the poly-A list: < 1 GB.
* the rebuild's transient: the FULL form of the promoted compact junctions (what the FULL
  path keeps anyway) plus the sorted request list (24 B per request).

FULL-mode cost of the plumbing (budget unset): `Breakpoint` 304 -> 328 B (`vp` u64, `ref_id`
i32 + `ord` u32, `src` 8 B) = +24 B per final and per temporary breakpoint (~0.6 GB on a
19.5 M-final run), and +8 B per collector entry. Measured locally: peak 245 -> 268 MB on
stress1, 303 -> 332 MB on stress2 (+9 %; the accounted structures grow by far less, so part
of this is allocator layout). With any budget set this is more than paid back; if the
unset-mode growth matters, the next step is to keep `vp` out of `Breakpoint` (parallel
vector in `temporary_breakpoints`) and the collector addresses in parallel vectors only
when a budget is set.

### Recommended budget fraction

On PD45808b_lo0003 the end-of-extract RSS was 38.3 GB of which the pending pool was ~32 GB
(REPORT §2); everything else at that moment ~6 GB, and the find_mates maps add ~1.8 GB.
Budget `N` bounds the FULL pending pool at `N` + one junction; the compact pool itself is
~200 B x pending (3.7 GB for 18.5 M) and is NOT inside `N`. So the peak is roughly
`N + 0.2 KB x pending + 6 GB (+1.8 GB in the mate pass)`. For a 24 GB tier
(`SD_MEM_T2`-style) `SD_MEM_BUDGET_FRAC = 0.4` (budget 9.6 GB) leaves ~8 GB for the rest
plus the compact pool; **0.4–0.5 of the job's -M** is the recommendation for dirty libraries,
0.6 for the usual 1–3 M-pending colonies (where the pool never trips and the fraction is
irrelevant). The farm validation below uses an explicit `PEARTREE_MEM_BUDGET_MB=8000` at
`-M 24000`, which the model above puts at ~18–20 GB peak.

Wall-time cost: the rebuild pass re-reads only the records behind promoted junctions, with
sorted addresses (sequential within a block; a seek per non-adjacent record). Locally 2.7 s
for 104 k records / 60 k seeks. On PD45808b (~3 M promoted x ~8 records ≈ 25 M records
spread over the genome) expect it to approach one more linear decode of the BAM, i.e. on the
order of the mate pass (~10–20 min single-threaded), against a 2–3 h discovery.

## 5. Farm validation (Jeremy runs; no placeholders)

Mac — push the branch (I never push):

```
git -C /Users/jeremy/Documents/PEAR_TREE push origin discovery-mem-budget
```

farm22 — separate worktree, writable CARGO_HOME inside it, build the discovery crate only:

```
git -C /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-hsc fetch git@github.com:jeremydeuel/PEAR-TREE.git discovery-mem-budget
git -C /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-hsc worktree add /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-membudget FETCH_HEAD
module load rust/1.87.0
export CARGO_HOME=/lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-membudget/.cargo-home
cargo build --release --manifest-path /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-membudget/rust/peartree-discovery/Cargo.toml
```

Is the test BAM still staged? (the staging area is cleaned by the pipeline's phase 5)

```
ls -la /lustre/scratch126/casm/staging/team273/jd43/2569/PD45808b_lo0003/mapped_sample/PD45808b_lo0003.sample.dupmarked.bam /lustre/scratch126/casm/staging/team273/jd43/2569/PD45808b_lo0003/mapped_sample/PD45808b_lo0003.sample.dupmarked.bam.bai
```

If it is gone, re-stage it the way `cluster/pipeline.sh` does (`stage-discover` uses
`stageBam.pl`-based staging into `$STAGING_ROOT/<project>/<sample>/mapped_sample/`); the
`.bai` is not needed for `--threads 1` (the rebuild seeks by virtual position, not by index).

Two runs on the same BAM and config (`--threads 1` as `discover_one.sh`), both with the memory
trace on: budget unset at 64 GB (the reference) and `PEARTREE_MEM_BUDGET_MB=8000` at 24 GB:

```
mkdir -p /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget
bsub -G team273 -q long -n 1 -M 64000 -R "select[mem>64000] rusage[mem=64000]" -J membudget_unset -o /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/unset.%J.log -e /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/unset.%J.err "PEARTREE_MEM_DEBUG=1 /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-membudget/rust/peartree-discovery/target/release/peartree-discovery --step discover --bam /lustre/scratch126/casm/staging/team273/jd43/2569/PD45808b_lo0003/mapped_sample/PD45808b_lo0003.sample.dupmarked.bam --out /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/unset.txt.gz --threads 1 --config /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-membudget/cluster/config.discovery.grch37.tprt2frag"
bsub -G team273 -q long -n 1 -M 24000 -R "select[mem>24000] rusage[mem=24000]" -J membudget_8000 -o /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/b8000.%J.log -e /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/b8000.%J.err "PEARTREE_MEM_DEBUG=1 PEARTREE_MEM_BUDGET_MB=8000 /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-membudget/rust/peartree-discovery/target/release/peartree-discovery --step discover --bam /lustre/scratch126/casm/staging/team273/jd43/2569/PD45808b_lo0003/mapped_sample/PD45808b_lo0003.sample.dupmarked.bam --out /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/b8000.txt.gz --threads 1 --config /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-membudget/cluster/config.discovery.grch37.tprt2frag"
```

(Optional third point, the baseline binary of the previous branch on the same BAM, is the
`memprof_base` job in REPORT.md §6 — its `base.txt.gz` can be `cmp`'d the same way.)

When both are done — byte-identity, the budget line, and the peaks:

```
cmp <(gzip -dc /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/unset.txt.gz) <(gzip -dc /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/b8000.txt.gz) && echo TXT_IDENTICAL
cmp <(gzip -dc /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/unset.txt.gz.evidence.tsv.gz) <(gzip -dc /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/b8000.txt.gz.evidence.tsv.gz) && echo EVIDENCE_IDENTICAL
cmp /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/unset.txt.gz.stats.json /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/b8000.txt.gz.stats.json && echo STATS_IDENTICAL
cmp /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/unset.txt.gz.splice.tsv /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/b8000.txt.gz.splice.tsv && echo SPLICE_IDENTICAL
grep -h -E '^\[mem\]|^\[mem-budget\]|disc-agree second|MAX MEM|Max Memory|Run time|CPU time' /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/unset.*.err /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/unset.*.log /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/b8000.*.err /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/membudget/b8000.*.log
```

(`[mem]` / `[mem-budget]` lines are on stderr = the `.err` files; LSF's `Max Memory` and run
time are in the `.log`. `splice.tsv` exists only if the config turns `splice_hallmark` on;
`cmp` of two missing files prints an error and no `SPLICE_IDENTICAL`, which is fine.)

Then the pipeline route, default off: `SD_MEM_BUDGET_FRAC=0.4` in the environment of
`cluster/pipeline.sh submit` (it is frozen into `run.env`) makes every stage-discover job
read its own `-M` back from LSF and run with `PEARTREE_MEM_BUDGET_MB = floor(0.4 x -M)`; the
per-sample log line `discovery memory budget <mb> MB (0.4 x job memlimit <lim> MB)` confirms
it. Unset = the previous invocation byte for byte.

## 6. Files

* `rust/peartree-discovery/src/compact.rs` — the pool (header vector + arena), 4-bit packing, tests.
* `rust/peartree-discovery/src/discovery.rs` — `BamIn` (readers with `virtual_position`),
  `query_intersects` + the hand-rolled chunk walk of the SPD-3 worker, `clip_breakpoint` /
  `disc_lite_of` / `short_lite_of` / `mate_seq_of` (shared by extract and rebuild), per-junction
  `cleanup` with `admit` / `trip_budget`, `get_mates` / `agree_mate_map` / `apply_mate_record`
  with the compact variants, `resolve_disc_agree` + `rebuild_compact` + `merge_by_slot`, the
  SPD-3 share / parent re-admission.
* `rust/peartree-discovery/src/model.rs` — `Breakpoint::{vp, ref_id, ord, src}`, `SrcOffsets`;
  `join()` records the cluster's addresses for a pending junction when a budget is set.
* `rust/peartree-discovery/src/evidence.rs` — `ShortCollector` entries carry the address.
* `rust/peartree-discovery/src/config.rs` — `mem_budget_mb`, `PEARTREE_MEM_BUDGET_MB`.
* `rust/peartree-discovery/src/mem.rs` — `src` bytes, compact pool in the `[mem]` line.
* `cluster/pipeline.sh` — `SD_MEM_BUDGET_FRAC` (default off), `job_memlimit_mb`.
