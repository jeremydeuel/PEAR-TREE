# Discovery memory blow-up on Pair38 (PD45808b_lo0003: 38–40 GB RSS) — diagnosis and fixes

Branch `discovery-mem-profile` (off local `main` 08d51e5). Crate `rust/peartree-discovery`.
Config under investigation: `cluster/config.discovery.grch37.tprt2frag` (`disc_agree_second_fragment = true`,
`evidence_sidecar = true`, `short_overhang_evidence = true`, `min_evidence_fragments_per_sample = 2`).

## 1. Diagnosis in one paragraph

The 19.5 M breakpoints that survive per-contig `cleanup` on PD45808b_lo0003 are not "retained clipped
reads" — they are **`disc_agree_second_fragment` PENDING junctions**: the trace itself says
`21,026,361 single-molecule junction(s) pending, 2,530,631 without a candidate anchor` (dropped at
cleanup), leaving 18.5 M pending + ~1 M ordinary finals = the 19.5 M in `final_L/final_R`. Every
single-molecule clip cluster whose clip is complex and >= 25 bp and that has **any** discordant anchor
within 500 bp is turned by `join()` into a *complete* consensus `Breakpoint` (clipped + unclipped
consensus, mates, sidecar `EvExtra` with CLIP/DISC/SHORT lite lists, `AgreePending` with inserts, clip
identities and up to 16 anchors) and held until the mate pass decides promotion; 85 % (15.5 M) are then
rejected. The old `PEARTREE_MEM_DEBUG` accounting counted only the struct + consensus + mate-vector
header (~680 B) and ignored `ev`, `agree` and the mate qname strings; measured with full accounting a
pending junction costs **~1.75–1.8 KB** (2.7x), which — together with ~4.3 GB of `Vec<Breakpoint>`
doubling slack — reproduces the 38.3 GB end-of-extract RSS from 19.5 M breakpoints to within a few
percent. The 2.5x "RSS vs accounted" gap is therefore live data, not allocator retention (jemalloc *is*
the default allocator, `dirty_decay_ms:0`). The ~10x memory for ~2x input is a product of three
factors that all grow with library "dirtiness" (chimeric / discordant pair rate): more single-molecule
complex clips -> more pending candidates; a denser discordant-anchor background -> a larger fraction
of them find an anchor within 500 bp and are kept (plus more anchors and DISC-lite entries each); and
a denser read background near each junction -> more SHORT-lite entries each. Pending count x payload
per pending both rise, so memory is super-linear in "junctions judged".

## 2. Where the bytes are: per retained breakpoint

Struct sizes (release build, `cargo test print_struct_sizes`): `Breakpoint` 304, `QualitySeq` 48,
`EvExtra` 168 (now 80), `AgreePending` 80, `ClipMol` 40, `DiscLite` 48, `ClipLite` 32, `ShortLite` 56,
`EvRec` 104, `PolyABreakpoint` 168.

Measured on the synthetic stress BAM (`tests/gen_pending_stress.py`, 98 k pending junctions, 25x
background, 1–3 anchors each), baseline code, `[mem]   final breakdown` line:

| field (per pending breakpoint)                    | bytes | counted by old `breakpoint_bytes`? |
|---------------------------------------------------|------:|------------------------------------|
| `Breakpoint` struct (Vec slot)                    |   304 | yes |
| clipped + unclipped consensus (seq+qual, 4 Vecs)  |  ~355 | yes |
| `mates` qname strings + `mate_seqs`               |   ~70 | header only |
| `ev` = `Box<EvExtra>` 168 + clip_lite 32 + disc_lite ~195 + short_lite ~286 | ~680 | **no** |
| `agree` = `Box<AgreePending>` 80 + inserts ~80 + clips ~60 + anchors ~195   | ~466 | **no** |
| **total as stored**                               | **~1,806** | 659 counted |

Farm cross-check: at `cleanup Y` RSS 38,350 MiB with 19.53 M finals and `tmp_bp` empty. Subtracting the
Vec doubling slack (len 9.77 M / 9.76 M -> capacity 2^24 = 16.78 M per side: 2 x 7.0 M x 304 B = 4.3 GB)
leaves 34.0 GB / 19.53 M = **1,743 B per breakpoint — the same number as the synthetic measurement
(1,747–1,806 B)**. Everything else (polyA 170 MiB, coverage bins, BAM buffers) is < 1 GB.

Transient per-contig extras (also now reported): `temporary_breakpoints` (every soft-clipped read, even
1-bp clips, as a full 600-B `Breakpoint`: 3.4 GB on contig 1), the short-overhang collector (`kept` is
essentially every read within 20 bp of any clipped read: ~56 B x reads of the contig), and
`sc_disc_tmp` (48 B per discordant primary). The find_mates peak (+1.8 GB) is the qname maps (qmap +
two qname sets = 3 copies of 8.9 M qnames) plus the `agree_map` (one `AgreeWait` per pending anchor).

## 3. Why Pair38 is ~10x Pair31 for ~2x input

`passed` (median 1.37 M vs 3.15 M) only counts promoted + ordinary junctions; the number that matters
for memory is `disc_agree.pending` minus `disc_agree.no_anchor` in `stats.json`, which is not in the
per-sample medians quoted. The mechanism: pending candidates ∝ single-molecule complex clips (a dirtier
library has more), the kept fraction ∝ P(a discordant anchor within 500 bp) which rises steeply with the
discordant-pair density (PD45808b: only 12 % had *no* anchor, i.e. anchors are near-ubiquitous), and the
per-pending payload (disc_lite, anchors, short_lite) also grows with that density. Three multiplicative
~1.5–2x factors give ~5–10x. To confirm on the farm, compare the `disc_agree` block of the two cohorts'
`stats.json` (command in section 6).

## 4. Fixes

### Implemented on this branch (output byte-identical, see section 5)

| fix | what | saving on PD45808b_lo0003 (19.5 M finals) | risk |
|-----|------|------------------------------------------:|------|
| A | `cleanup` stages each contig's new breakpoints and appends with `reserve_exact` instead of pushing one by one into `final_*_breakpoints` (which doubled capacity) | **~4.3 GB** at the end-of-extract peak (2 x (16.78 M − 9.77 M) x 304 B) | none to output: order and indices unchanged; one realloc per contig |
| B | `EvExtra` no longer carries four empty `Vec<EvRec>` per breakpoint; the full records live in a lazily boxed `EvRecs` (allocated only for the ~0.5 % emitted loci) — 168 -> 80 B | **~1.7 GB** (88 B x 19.5 M) | mechanical; sidecar rows unchanged |
| C | `AgreePending::clips` (needed only for the same-molecule check while selecting anchors in `cleanup`) is freed right after | **~0.8–1.0 GB** (Vec + 1–2 `ClipMol` of 40 B, jemalloc 48-B class) | none |
| D | `AgreePending::anchors` move into the mate-pass `agree_map` instead of being copied, so 48 B x anchors per pending junction are freed before the mate scan | **~3–4 GB at the find_mates peak** (40.1 GB -> below the extract peak) | none; `anchors` is never read after the map is built |
| E | `PEARTREE_MEM_DEBUG` accounting now includes `ev`, `agree`, mate strings, Vec slack and the per-contig collectors, and prints a per-field breakdown line | diagnosis only | none |

Expected end-of-extract peak: 38.3 GB -> **~31–32 GB** (−18 %); find_mates peak 40.1 GB -> ~28–30 GB.
Measured on the small synthetic case (where the slack share is tiny): peak RSS 278 -> 256 MB (−8 %),
per-pending breakpoint 1,806 -> 1,747 B before the mate pass. These do NOT bring a 65-GB sample under
48 GB on their own; they are the safe part. The structural fix is below.

Note: the SPD-3 contig-parallel path (`--threads > 1`) merges worker vectors with `append` and still
grows by doubling; the farm runs `--threads 1` (`cluster/discover_one.sh`), where fix A applies.

### Proposed (not implemented — each needs its own byte-identity validation)

1. **Spill pending junctions (the real fix, ~15–20 GB).** 85 % of pending payload is thrown away by
   `resolve_disc_agree`. Serialise each pending `Breakpoint` at the end of `cleanup` into one
   contiguous byte blob (or an append-only temp file on the job's scratch) and keep in the final vector
   only a stub with the fields the mate pass touches: side, position, `agree.inserts` (for `agrees()`),
   `mates` qnames (for the qmap) and a blob offset; restore the 15 % promoted ones after the scan.
   Byte-identical by construction if the (de)serialiser round-trips every field (a unit test can assert
   that). Live cost per pending junction drops from ~1.75 KB to ~100–150 B (in-memory blob: ~500 B,
   because one allocation replaces ~20 and the padding/Vec headers go). Complexity: a serialiser for
   `Breakpoint` + `EvExtra` lite lists + `AgreePending` (~250 lines) and index bookkeeping; `mate_seqs`
   captured for a spilled stub during the scan must be kept on the stub and re-attached.
2. **Re-derive SHORT-lite in the sidecar pass instead of attaching it to every final breakpoint
   (~5 GB).** `short_lite` (~286 B/bp) is consumed only for emitted loci (97 k of 19.5 M here), and the
   sidecar pass already does indexed fetches around exactly those loci. The collector keeps a superset
   of `short_candidate` matches and `attach_short` re-applies the exact geometry, then `dedup` +
   `select_lowest`, so recomputing from a `[b − reach, b + window]` fetch (same start bound, sort by
   `(start, frag, flag)` before dedup) yields the identical set. Needs a two-stage sidecar fetch
   (SHORT reads first, then their mates) because `build_short_requests` currently runs before the fetch.
3. **Share DISC-lite and agree anchors (~2–3 GB).** For a pending junction `agree.anchors` is (almost
   always) a subset of `ev.disc_lite` selected by the same geometry from the same per-contig list; store
   anchors as indices into `disc_lite` with an owned fallback when the subset property fails (different
   caps / `same_molecule` tolerance).
4. **Compact `ShortLite`/`DiscLite`** (56/48 B -> ~32 B with i32 positions relative to the junction and
   u16 clip lengths): ~2 GB; edge cases (positions > 2^31, clips > 65535) need guards.
5. **Do not materialise 1–11 bp junk clips as full `Breakpoint`s in `temporary_breakpoints`** (3.4 GB
   per-contig transient): they can only ever contribute to a cluster's consensus, so a compact raw form
   suffices; output-neutral but `stats.json` counters must be kept identical.
6. **Default-off operational knobs** (change output, so for a dirty-library tier only): a cap on pending
   junctions per contig/sample, or a stricter anchor requirement (e.g. anchor mate MAPQ or
   `disc_agree_span` 500 -> 300). With `disc_agree_second_fragment = false` the retained set collapses to
   the ordinary finals (~1 M) and memory returns to the Pair31 regime — at the cost of the
   single-molecule+agreeing-discordant junctions (2.99 M promoted here, before downstream gates).

## 5. Test evidence

* `cargo test --release` in `rust/peartree-discovery`: 79 passed, 0 failed (includes the disc_agree
  end-to-end tests through `handle_record` / `cleanup` / `apply_mate_record` / `resolve_disc_agree`).
* Synthetic stress BAM (`tests/gen_pending_stress.py --contigs 2 --len 8000000 --cov 25
  --junction-every 150 --seed 1`, config `tests/config.pending_stress` = the tprt2frag keys without the
  farm-only 2bit/exon tracks): 105,589 pending / 7,497 no-anchor / 14,152 promoted / 83,940 rejected,
  1,596 loci emitted, 28,188 sidecar rows.
  * baseline (main 08d51e5) vs this branch, `--threads 1`: `.txt.gz` (18,992 lines), `.evidence.tsv.gz`
    (28,189 lines) and `stats.json` **byte-identical**; peak RSS 278.2 MB -> 256.0 MB.
  * `--threads 2` (SPD-3 parallel extract): baseline vs branch byte-identical, and serial vs parallel
    byte-identical.
* Not measured here: RSS on Linux/jemalloc at farm scale (macOS has no `/proc`, `rss_bytes()` is None;
  `/usr/bin/time -l` peaks were used). The farm run below is the real before/after.

## 6. Farm validation (Jeremy runs; no placeholders)

Mac — push the branch (I never push):

```
git -C /Users/jeremy/Documents/PEAR_TREE push origin discovery-mem-profile
```

farm22 — a separate worktree so the deployed checkout is untouched, build the discovery crate only:

```
git -C /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-hsc fetch git@github.com:jeremydeuel/PEAR-TREE.git discovery-mem-profile
git -C /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-hsc worktree add /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-memprof FETCH_HEAD
module load rust/1.87.0
export CARGO_HOME=/lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-memprof/.cargo-home
cargo build --release --manifest-path /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-memprof/rust/peartree-discovery/Cargo.toml
```

Baseline (deployed binary) and branch binary on the same staged BAM, same config, `--threads 1` as
`discover_one.sh`, with the memory trace on (64 GB so the baseline cannot be killed):

```
mkdir -p /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof
bsub -G team273 -q long -n 1 -M 64000 -R "select[mem>64000] rusage[mem=64000]" -J memprof_base -o /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/base.%J.log -e /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/base.%J.err "PEARTREE_MEM_DEBUG=1 /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-hsc/rust/peartree-discovery/target/release/peartree-discovery --step discover --bam /lustre/scratch126/casm/staging/team273/jd43/2569/PD45808b_lo0003/mapped_sample/PD45808b_lo0003.sample.dupmarked.bam --out /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/base.txt.gz --threads 1 --config /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-hsc/cluster/config.discovery.grch37.tprt2frag"
bsub -G team273 -q long -n 1 -M 64000 -R "select[mem>64000] rusage[mem=64000]" -J memprof_new -o /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/new.%J.log -e /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/new.%J.err "PEARTREE_MEM_DEBUG=1 /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-memprof/rust/peartree-discovery/target/release/peartree-discovery --step discover --bam /lustre/scratch126/casm/staging/team273/jd43/2569/PD45808b_lo0003/mapped_sample/PD45808b_lo0003.sample.dupmarked.bam --out /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/new.txt.gz --threads 1 --config /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-memprof/cluster/config.discovery.grch37.tprt2frag"
```

When both are done — byte-identity and peaks:

```
cmp <(gzip -dc /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/base.txt.gz) <(gzip -dc /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/new.txt.gz) && echo TXT_IDENTICAL
cmp <(gzip -dc /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/base.txt.gz.evidence.tsv.gz) <(gzip -dc /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/new.txt.gz.evidence.tsv.gz) && echo EVIDENCE_IDENTICAL
cmp /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/base.txt.gz.stats.json /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/new.txt.gz.stats.json && echo STATS_IDENTICAL
cmp /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/base.txt.gz.splice.tsv /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/new.txt.gz.splice.tsv && echo SPLICE_IDENTICAL
grep -h -E '^\[mem\]|MAX MEM|Max Memory' /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/base.*.err /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/base.*.log /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/new.*.err /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/memprof/new.*.log
```

(`[mem]` lines go to stderr = the `.err` files; LSF's `Max Memory` is in the `.log`.) The new binary's
`[mem]   final breakdown:` lines give the real per-field split on PD45808b, which decides between
proposals 1–3 above.

Pair31 vs Pair38 pending counts (the super-linearity check), from the existing run outputs under the
hsc run root (`$HSC_ROOT/<P>/<P>/discovery/`, `HSC_ROOT` defaults to `$JD/hsc` in `cluster/hsc_run.sh`):

```
for f in $(find /lustre/scratch126/casm/teams/team273/users/jd43 -path '*Pair3[18]*' -name '*.txt.gz.stats.json' 2>/dev/null); do printf '%s\t' "$f"; grep -o '"disc_agree": {[^}]*}' "$f" || echo no_disc_agree_block; done
```

## 7. Files

* `rust/peartree-discovery/src/discovery.rs` — fixes A, C, D; collector sizes in the mem report.
* `rust/peartree-discovery/src/evidence.rs` — fix B (`EvExtra::recs`, `EvRecs`), `ShortCollector::mem_stats`.
* `rust/peartree-discovery/src/mem.rs` — full accounting + breakdown line + struct-size probe test.
* `rust/peartree-discovery/tests/gen_pending_stress.py`, `tests/config.pending_stress` — the stress harness.
