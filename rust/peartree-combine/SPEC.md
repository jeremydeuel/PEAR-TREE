# peartree-combine — behavioural specification

Rust port of PEAR-TREE step 2, `python src/main.py --step combine_insertions` (driver
`src/combine_insertions.py`). Normative reference: the Python at commit d74a897. Where this
document and the Python disagree, **the Python wins** — fix this document.

Python file abbreviations: `CI` = src/combine_insertions.py, `INS` = …_insertion.py,
`ISX` = …_intersect_insertions.py, `RF` = …_region_filter.py, `GS` = …_get_sequence.py,
`EV` = …_evidence.py, `TF` = …_tprt_filters.py, `IC` = src/indel_consensus.py,
`SC` = src/sequence_checks.py, `QS` = src/quality_seq.py. `file:N` = line N.

---------------------------------------------------------------------------------------------

## 0. Scope and decisions

**Goal.** Same output files, byte-identical after decompression (gzip headers/levels are not
specified), for the TPRT config (`cluster/config.py.grch38.tprt`, production) and the legacy
configs (`cluster/config.py.grch38` / `.grch37`, with and without evidence sidecars).
`tests/equiv.sh` is the acceptance test.

**Dropped: the pooled ≥2-independent-fragment gate.** Not ported: `require_independent_fragments`
(EV:789, EV:988-997 — fail reason `n_independent<k`, reasons `no_evidence_reads` / `LEFT(polyA)+…`),
the `n_ind` argument / `"few_fragments"` reason of `far_pair_verdict` (TF:293), and the gate
clause of the far-pair split condition (EV:962, `not gate or …` is always true). The rule is now
enforced by discovery (`min_evidence_fragments_per_sample`). `Config::load` returns an error when
`require_independent_fragments` is truthy (the Rust output would silently differ from Python).

**`n_independent` decision — the lenient dedup IS ported, in full.** `independent_clusters`
(EV:472) is not just the gate's input:
1. its clusters set the consensus vote weights (`1/len(cluster)`) and the depth unit (`group` =
   cluster index; consensus stops below `min_independent_fragments` agreeing clusters) → they shape
   `clip_consensus`, `consensus_depth`, `polya_*`, `beyond_polya*` in evidence.tsv AND, with
   `indel_aware_consensus`, the clips written to `combined.txt.gz`, the remap FASTQs and the
   genotyping contract;
2. `n_independent` → `supported` (0/1) is read by `tools/rte/annotator.py:245` →
   `score.py` `junction_supported` (+0.5 per supported junction); `n_independent`,
   `n_duplicates`, `n_cross_sample_identical` (as `cross_sample_identical`, a score penalty),
   `n_short_*` are read by `tools/rte/inputs.py`; `cluster/tprt/compare_arms.py` reads
   `supported` and `fail_reason`.
So every column is computed exactly as in Python. Only the *dropping* is gone.

**Dead data dropped.** Discovery `MATE` FASTQ records (`Insertion.left_mates/right_mates`) are
never written or inspected by any combine output (only `Insertion.__str__`, unused). The parser
validates and discards them. Sidecar MATE *rows* are kept (they are evidence-file content).

**Duplicate input basenames** are rejected (Python keys files by basename and would silently
conflate two colonies' sidecars).

---------------------------------------------------------------------------------------------

## 1. Interface

### 1.1 CLI (src/main.py:148-188)

```
peartree-combine [--step combine_insertions] --config <config.py|config.json>
                 --discovery_files F1 [F2 ...] --out <STEM> --threads <N>
```
* `--discovery-files` accepted as an alias; values run until the next `--option`.
* `--out`: a trailing `.gz` is stripped (main.py:153). `--config` default: `$PEARTREE_CONFIG`,
  else `src/config.py` relative to the cwd.
* Checks before any work (main.py:156-170): `threads <= cpu count`; every input exists;
  `samtools_executable`, `bowtie2_executable` exist; `{bowtie2_index}.1.bt2` exists. Then the
  genome is opened (GS:22 opens `genome_2bit` at import → a missing 2bit is fatal at startup).
* Cluster switch: `cluster/combine_mei.sh:94` and `cluster/pipeline.sh:523` call
  `"$VENV/bin/python" [-u] "$PT_ROOT/src/main.py" --step combine_insertions --discovery_files … --out … --threads …`;
  the drop-in is `"$PT_ROOT/rust/peartree-combine/target/release/peartree-combine" --config "$PT_ROOT/src/config.py" <same args>`
  with `PEARTREE_PYTHON="$VENV/bin/python"` (the config loader runs python once).

### 1.2 Config

`--config X.json`: JSON of the whole `CONFIG`. Otherwise X is executed by
`$PEARTREE_PYTHON` (else `python3`) running
`json.dumps(runpy.run_path(X)["CONFIG"], default=lambda o: None)` (generated configs `exec()`
another file — `runpy` handles that). Keys consumed (defaults = Python `.get` defaults;
"required" = Python indexes without default):

| section.key | default | used by |
|---|---|---|
| combine_insertions.genome_2bit | required | genome (§2.3) |
| …exclude_files_with_many_insertions | required | import (§3.1) |
| …samtools_executable, bowtie2_executable, bowtie2_index, bowtie2_index2, bowtie2_index2_lo | required | remap (§5) |
| …clean_remap_max_insertion | `discovery.min_clip_len` (then required) | §5.1 |
| …clean_remap_min_as | -15 | §5.1 |
| …trim_far_flank_before_remap | False | §6.2 |
| …keep_polya_one_sided | False | §3.3 |
| …merge_tolerance_bp | `int(v or 0)` → 0 | §3.3, §4.5, §4.7 |
| …polya_aware_clip_agreement | False | §3.3 |
| …require_independent_fragments | False | **must be falsy** (§0) |
| …min_independent_fragments | 2 | §4 |
| …indel_aware_consensus | False | §4 |
| …dup_coord_tolerance / dup_max_edit / dup_max_edit_frac / dup_mate_min_mapq | 5 / 3 / 0.02 / 20 | §4.2 |
| …polya_min_len | 8 | §4.2-4.6 |
| …count_short_overhang | False | §4.3 |
| …short_overhang_min_bases / short_overhang_min_ref_mismatch | 5 / 2 | §4.3 |
| …short_mate_max_dist / short_mate_min_mapq | 1000 / 20 | §4.3 |
| …slippage_reject / far_pair_strict / far_pair_split | False / False / True | §4 |
| …rte_library | `v or "resources/rte_library"` | §4.6 |
| …slippage_min_ref_run_combine / slippage_min_str_len / slippage_min_structured / slippage_max_period / slippage_junk_frac | 8 / 12 / 10 / 6 / 0.5 | §4.6 |
| …far_pair_max_tsd_deletion / far_pair_tsd_max / far_pair_min_polya / far_pair_allow_antisense / far_pair_colony_frac / far_pair_colony_tol | 30 / 40 / 10 / False / 0.2 / 5 | §4.5-4.6 |
| genotyping.max_bases | required | §6.3 |
| discovery.min_clip_len | (only as fallback above) | §5.1 |

Flags use Python truthiness. `rte_library` relative path resolution (TF:173 `_resolve`): as given
if it exists relative to the cwd, else `<repo>/<path>` (Rust: `<config dir>/..`, `$PT_ROOT`,
binary-relative `../../../..`; first existing). The harness runs from the repo root.

### 1.3 Output files (STEM = `--out` minus `.gz`)

| file | when | §  |
|---|---|---|
| `STEM.combined.txt.gz` | always | 6.1 |
| `STEM.genotyping.txt.gz` | always | 6.3 |
| `STEM.fq.gz` | always; final content = the 2nd remap input (clip FASTQ) when the 2nd remap ran, else the consensus FASTQ | 6.1/6.2 |
| `STEM.bam` | bowtie2 end-to-end; **not regenerated if it exists** | 5.1 |
| `STEM.insertionsonly.bam` | bowtie2 local; **not regenerated (and `STEM.fq.gz` not rewritten) if it exists** | 5.2 |
| `STEM.insertions.evidence.tsv.gz`, `STEM.insertions.reads.fa.gz` | only when ≥1 accepted file has a sidecar | 6.4, 6.5 |
| `STEM.combined.splice.tsv` (uncompressed) | only when ≥1 input `<f>.splice.tsv` exists | 6.6 |
| `STEM.evidence_shards/` | scratch, created fresh, removed after the evidence outputs | 4.0 |

---------------------------------------------------------------------------------------------

## 2. Core types and primitives

### 2.1 Insertion (INS:40-128) → `model::Insertion`, parser `insertion.rs`

Discovery file: gzip FASTQ, title `@contig:start-end:SIDE:FIELD`. Parser (INS:166): read
4-line records; title line `strip()`ped, skipped if empty; `l[1:].split(":")` must give exactly
4 parts and `positions.split("-")` exactly 2; seq/plus/qual lines `strip()`ped; plus must be `+`;
qual = `ord(c)-33`. A record group = consecutive titles with the same (contig, start, end);
repeated FIELD overwrites; FIELD containing `MATE` → mate (discarded here). Groups are yielded in
file order; a non-consecutive repeat of an id is a separate group.

Decoding one group (INS:41-98), RIGHT first then LEFT (LEFT overwrites `type`):
* RIGHT: if `RIGHT:ALIGNED` and `RIGHT:CLIPPED` present → `right_aligned = ALIGNED.revcomp()`,
  `right_clipped = CLIPPED`, `right_pos = int(end)`; elif end starts `oneside_` → no seq,
  `right_pos = int(end[8:])`, type 4, `open_side = RIGHT`; elif `disc_` → no seq,
  `right_pos = int(end[5:])`, type 4; else (must start `polyA_`) `right_clipped = RIGHT:CLIPPED_POLYA`,
  `right_pos = None`, type 1.
* LEFT: `LEFT:ALIGNED`+`LEFT:CLIPPED` → `left_clipped = CLIPPED.revcomp()`, `left_aligned =
  ALIGNED`, `left_pos = int(start)`; `oneside_` → type 5, `open_side = LEFT`; `disc_` → type 5;
  else `left_clipped = LEFT:CLIPPED_POLYA.revcomp().lower()`, `left_pos = None`, type 2.
* no type set → FULL_INFO (3). `files = [basename]`, `member_loci = [(basename, name)]`,
  `name = f"{contig}:{start}-{end}"` with the RAW tokens (`model::Tok` round-trips them).

`left_consensus = left_clipped.revcomp().lower() + left_aligned`;
`right_consensus = right_aligned.revcomp() + right_clipped.lower()` (asserts: not on the
open/poly-A side). `QualitySeq.lower()/upper()` change bases only; `revcomp` reverses
qualities. FASTQ (QS:78): `@{title}\n{seq}\n+\n{chr(min(q+33,126)) per q}\n`.

### 2.2 sequence_matching_score (SC:89) → `seq::sequence_matching_score`

`min_len = min(24, min(len(s)))`; for i < min_len: tally `bases[base.upper()] += w` (w = quality
for a QualitySeq, 1 for a str; keys A,T,G,C,N — another char raises KeyError); `+1` if
`sum>0 and max/sum > 0.5` else `-2`; return `score/min_len` (ZeroDivisionError if min_len = 0).
Rust panics in both error cases (no real input triggers them).

### 2.3 Genome (GS) → `genome.rs`

py2bit, default `storeMasked=False` → uppercase; N blocks → `N`. `get_sequence(name, s, e)`:
name mapping (GS:37-52; index-order substring fallback), `s >= e` → `""`; py2bit clamps `e` to
the sequence length; `s < 0` or `s >= clamped e` → exception → `""` (measured: (L-5, L+5) → 5
bases; (L+1, L+5) → ""; (-3, 5) → ""). All callers uppercase.

### 2.4 edlib → `align.rs` (implemented, parity-tested 3000/3000 vs python-edlib 1.3.9)

Vendored edlib 1.2.7 (the C library python-edlib 1.3.9 wraps). Call sites:

| site | mode | task | k | equalities |
|---|---|---|---|---|
| ISX:72 `_shift_tolerant_agree` | HW | distance | `max(1, len(q)//4)` | – |
| EV:323 `_semi_close` | SHW | distance | `budget(m)` | – |
| EV:716 SHORT `ed_ref` / EV:717 `ed_cons` | HW / SHW | distance | -1 | – |
| TF:165 `slippage_junction` | SHW | distance | `budget` | – |
| IC:114 `_align` | SHW (anchored) / HW | path | -1 | WILDCARD_N |
| IC:394 `_pick_seed` | SHW | distance | -1 | WILDCARD_N |

### 2.5 Numerics (normative)

* `sum()` of floats = CPython ≥3.12 Neumaier sum (`pyfmt::py_sum`); only IC `_decide`'s `others`.
  The farm runs python-3.12.0. All other float accumulations are explicit `+=` loops (naive).
* `round()` = half-to-even (`pyfmt::py_round`); `int()` truncates.
* `budget(n) = max(max_edit, ceil(frac*n))` on the IEEE product (`pyfmt::dedup_budget`).
* `_median` (`pyfmt::py_median`). Expression evaluation order must match Python
  (`r.weight * wsum / len(ins_b)` = `(w*s)/n`).
* No float is ever formatted into an output file.

---------------------------------------------------------------------------------------------

## 3. Stages before the evidence

### 3.1 Import (CI:116-129)

For each input in command-line order: parse; keep records with `len(contig) < 6` and contig
not `MT`/`chrM`; if the kept count `> exclude_files_with_many_insertions` the file is excluded
(not "accepted"), else its records are appended and the file is accepted. Parsing may run in
parallel; the concatenation order is fixed.

### 3.2 discovery_breakpoints (EV:1306; only when `far_pair_strict`)

Over all accepted records in order (before intersect): `(contig, LEFT)` gets
`(left_pos, sample(files[0]))` when type ∉ {2,5} and left_pos set; `(contig, RIGHT)` gets
`(right_pos, …)` when type ∉ {1,4} and right_pos set.

### 3.3 intersect_insertions (ISX:160) — fixes the order of every output

1. Bucket every record in input order into three insertion-ordered maps keyed by
   `(contig, left_pos, right_pos)` (None allowed): `full` (type 3), `polyA` (1, 2), `disc` (4, 5).
2. For each `full` key in order: one record → itself. Several → agreement checks, first failure
   → value None: if `polya_aware_clip_agreement`: `clips_agree(right_clipped…)` and
   `clips_agree(left_clipped…)` (defaults 0.6/8/6, plain bases); else
   `sms(right_clipped QualSeqs) < 0.6`, then `sms(left_clipped) < 0.6`; then (both modes)
   `sms(left_aligned) < 0.6`, `sms(right_aligned) < 0.6` (QualSeq-weighted). Pass → fold
   `combined += h` over the list in order (INS:130; for FULL+FULL: longer clipped replaces,
   longer aligned replaces, strictly longer; files/member_loci concatenated).
3. `merge_tolerance_bp > 0` — fuzzy merge of live `full` values (ISX:254):
   `n_rec[k] = len(member_loci)` (before any absorb); `_fuzzy_clusters` (see `intersect.rs`
   doc; sort key (-n_rec, contig NAME, L, R)); per cluster rep = first; each other member `o`
   is absorbed (`_absorb(rep, o, both)`, value → None) iff `_shift_tolerant_agree` holds for
   LEFT and RIGHT clips.
4. `keep_polya_one_sided`: every polyA record (in polyA key order, records in order) is turned
   into a one-sided locus (`_polya_to_one_sided`, ISX:27: poly-A end gets its coordinate from the
   NAME token `polyA_<P>`, clipped/aligned None, type 4/5, open_side) and appended to `disc` under
   its ORIGINAL polyA key (new keys at the end); polyA emptied. Otherwise polyA records are
   dropped (ISX:285 `continue` — the rest of that loop is dead code).
5. For each `disc` key in order: rep = first record, replaced by a later one whose real-side clip
   (`right_clipped` for type 5 else `left_clipped`) is non-None and strictly longer than rep's (or
   rep's is None). If rep has `open_side`: for every other record (in order): append its files
   not yet present (**files only — member_loci are NOT extended**; mates not stored); keep the
   strictly longer real-side aligned flank. Feature-A `disc_` reps (no open_side) pool nothing.
   Insert `full[rep.name] = rep` (string key → appended after all tuple keys, in disc order).
6. `merge_tolerance_bp > 0` — `_fuzzy_one_sided` (ISX:133): live values with `open_side`,
   grouped by (contig, real side) [groups are independent]; each group sorted by
   (-len(files), real-side pos, NAME string); greedy reps: an item joins the FIRST rep with
   |Δpos| ≤ tol and `_shift_tolerant_agree(clip(r), clip(v))` → `_absorb(r, v, (side,))`,
   value → None.
7. Output = non-None values in map order. `uid` = output index.

`_absorb(target, other, sides)` (ISX:116): `mem = other.member_loci + [(f, other.name) for f in
other.files]`, deduplicated in order; each `m` not in target.member_loci is appended and
`member_sides[m] = sides` (member_sides starts as a copy of target's or empty); files of other not
in target.files appended.

### 3.4 filter_dense_regions(bin_range=100, ins_cutoff=4) (RF:34) — see `region_filter.rs`.

---------------------------------------------------------------------------------------------

## 4. Evidence (only when some accepted file has a sidecar)

Sidecar of `F` = `F.evidence.tsv.gz` (Rust discovery), else `<F minus .txt.gz>.evidence.tsv.gz`
if only that exists (EV:61). Presence test at CI:149 and `have` (EV:800) agree.

### 4.0 apply_evidence / _judge_chunk (EV:766-948)

1. `members[k] = _member_loci(ins_k)` = member_loci + [(f, ins.name) for f in files], dedup in
   order. (Python maps through a name-keyed dict; names are unique here.)
2. Route sidecar rows: for each accepted file with a sidecar, rows whose (file, locus) is a
   member of some insertion. Row parsing §4.1; a row is kept for every chunk that needs it;
   per (member, side) the order is the sidecar line order. Chunking is free (outputs never depend
   on it); chunks may be evaluated in parallel; results merged in insertion order.
3. Matcher (`LibraryMatcher(rte_library)`) built iff `slippage_reject or far_pair_strict`.
   `ref_fetch` = the genome (always available — GS imports cleanly when the 2bit exists).
4. Per insertion `ins` (in order): `open = _open_side(ins)`; for side in (LEFT, RIGHT) except
   `open`: `recs.append(evaluate(ins, side))` where evaluate (EV:894) =
   * pooled rows = for m in members (order) if `side in member_sides.get(m, SIDES)`: rows(m, side);
   * `_reanchor(pooled, side, ins junction)` (EV:1257);
   * `rec = evaluate_junction(ins.name, side, pooled)` (§4.3);
   * `rec.aligned = _aligned_part(ins, side)`; `loci = sorted({locus NAME for _, locus in members})`;
     if `loci != [ins.name]`: `rec.member_loci = ",".join(loci)`.
5. `ok = True`. If `far_pair_strict` and open is None and 2 recs: `v = _far_pair_check(…)`
   (§4.5). If v = (reason, pside): if pside is not None and `far_pair_split`: `old = ins.name`;
   `_to_one_sided(ins, pside)`; `ins.member_sides = {m: (pside,) for m in members}`; the pside
   record gets `insertion_id = ins.name`, `member_loci = ",".join(sorted(loci ∪ {old}))`,
   `fail_reason = "split_from_far_pair:" + reason`; `recs = [that record]` (ok stays True).
   Else every rec `fail_reason = "far_pair:" + reason`, ok = False.
6. If ok and `slippage_reject`: `why = _slippage_check(…)`; non-empty → every rec
   `fail_reason = why`, ok = False.
7. If ok: if every `ins.files` entry has a sidecar: `supported = 0 if n_independent <
   min_independent_fragments else 1` per rec (else stays `NA`). **No drop** (gate removed).
8. If ok and `indel_aware_consensus`: `_replace_clips(ins, recs)` (EV:1283).
9. Records' reads are rendered (`_render_reads(rec.insertion_id, rec)`) and may be moved to disk.
10. Collect: kept (ok) and failed lists in insertion order. `records[name]`: setdefault over
    kept, then over failed whose name is not a kept name (those names go to `failed`).

### 4.1 Rows (EV:80-148) → `evidence/row.rs`

Header line `rstrip("\n").split("\t")`; data lines likewise; a row whose field count ≠ header's
is skipped. Columns by name with the defaults listed in `row.rs`. `outward_clip`, `quals`,
`mapped`, `mate_mapped`, `allele_forward_seq` as documented there.

### 4.2 Fragments and dedup (EV:178-519) → `evidence/dedup.rs`

Port literally; doc comments in `dedup.rs` give the order rules. Notes:
`_is_dup` tries `(a, b.swapped())` then `(a.swapped(), b)` only when `a.strand != b.strand` and
the oriented test failed, skipping pairs where either primary is unmapped. `_identical` compares
`(mate_outer() or mate_coord())` tuples, `primary.seq.upper()`, `mate_seq()`.
`_seq_close` reversal uses Python string reversal of the raw strings.

### 4.3 evaluate_junction (EV:569) → `evidence/junction.rs`

1. `frags = collapse_fragments(rows)`; `short` = frags with SHORT primary; `main` = the rest.
2. `clusters = independent_clusters(main)`.
3. Consensus reads: for gi, cluster in order; `w = 1.0/len(cluster)`; for each fragment, for
   each of its rows in row order: SHORT → skip; CLIP → `ClipRead(*outward_clip(), gi, w, True)`;
   other with non-empty seq → `ClipRead(seq, quals(), gi, w, False)`.
4. `consensus = IAC(reads, min_depth=min_ind, polya_min_len)`. If `indel_aware_consensus`:
   `combined_consensus = consensus if all reads anchored else IAC(anchored reads, …)`.
5. `n_independent_no_short = len(clusters)`.
6. If `count_short_overhang` and short: `has_clip` = any main primary CLIP; `vcons = IAC(anchored,
   min_depth=1).seq if anchored else ""`; each short fragment: `_short_overhang_check` (EV:666)
   → reason (counted) or used (+ `n_short_mate_inside += _mate_inside`); `n_short_used/rejected`;
   if used: re-run `independent_clusters(main + used)` (replaces clusters, n_dup, n_cross, stats).
7. Dropped rows: fragments in short − used, and every SHORT-role row not of a used fragment.
   `n_mates` = MATE rows kept; `n_reads` = kept − mates; `n_fragments = len(main + used)`;
   `n_samples` = distinct samples of those; `n_independent = len(clusters)`.
8. `polya_end = consensus.polya_base is not None or any POLYA row or
   _clips_start_with_polya(CLIP rows)` — over the ORIGINAL pooled rows (incl. dropped SHORT).

### 4.4 Indel-aware consensus (IC) → `consensus.rs`

Port line by line. Ordering/numerics that must be reproduced:
* `_Col.w`, `_Col.groups`, `_Col.ins` are insertion-ordered maps (first vote creates the key);
  `max(x.w, key=x.w.get)` / `max(x.ins.items(), key=w)` = first key with the maximal value;
  group sets only need cardinality.
* `others = py_sum(v for k, v in x.w.items() if k != best)` (Neumaier, map order).
* weights: `r.weight * _qw(rq)`; D: `_qw((qn + qp) / 2.0)`; I: `wsum += _qw(rq)` then
  `e[0] += r.weight * wsum / len(ins_b)`.
* `margin = max(1, int(round((wb - others) * 40)))` (half-even); `score = min(93, margin)`.
* `_pick_seed`: candidates = first 8 distinct RLE strings in a STABLE sort by raw length
  descending; probes = all anchored RLE strings if ≤ 40, else `bases[::max(1, n//40)][:40]`;
  key `(total distance, -len(c))`, strict `<`.
* `_vote` mates: try orientation 0 then 1 (`rle(revcomp(seq), qual[::-1])`), keep the strictly
  larger overlap bp among `ok` ones.
* `_decide`, `_expand`, `_expand_bases`, `_polya` exactly as written (note `_expand_bases`
  re-RLEs with quality 30, and the `cols` dict is indexed by integer column; `last = max(cols)`).

### 4.5 Far-pair / slippage checks (EV:1329-1452) → `evidence/filters.rs`

`_far_pair_check`: `gap = right_pos - left_pos`; not `far_geometry(gap)` → None.
`tol = max(int(merge_tolerance_bp or 0), int(far_pair_colony_tol))`. Per side present (LEFT,
RIGHT order): clips = dedup-in-order of [`_outward_clip(rec, ins)`, old discovery clip
uppercased or ""], empties removed; colonies; inside mates. `(reason, pside) =
far_pair_verdict(…)`; if reason == "": `line, j = outward_reference(genome, contig,
junction(pside), pside)`; `slippage_junction(_outward_clip(rec[pside]), line, j)` non-empty →
reason `polya_side_slippage`. Return None if reason == "" else (reason, pside).

`_slippage_check`: per rec (order): `line, j = outward_reference(…, rec.side)`; `slip =
slippage_junction(_outward_clip(rec), line, j)`; if slip and `rec.consensus.seq` is longer than
`_outward_clip(rec)` and `slippage_junction(consensus.seq.upper(), …)` is empty → slip = "". Then
for each side with slip (rec order): `others` = recs of the other side whose slip is empty; if
none of them `_carries_element` → return `f"slippage:{SIDE}({why})"`. Return "".

### 4.6 TPRT primitives and LibraryMatcher (TF) → `tprt.rs`, `library.rs`

Port literally (`repeat_at`, `strip_repeat`, `structured_len`, `slippage_junction`,
`leading_polyt`, `after_polyt`, `far_geometry`, `colonies_consistent`, `far_pair_verdict`).
Python `str.count`, `max("ACGT", key=over.count)` (first max), `in unit` membership.
`slippage_junction` reads `polya_min_len` (8) too.
LibraryMatcher: mappy `Aligner(path, k=11, w=3, min_cnt=1, min_chain_score=15, min_dp_score=20,
best_n=3)`; `hit()` per `library.rs`. `far_pair_verdict` order: poly-T test per side → exactly one
→ pside; complex side hit = first non-None `hit(c)` over its clips; none → inside-mate hits
(non-FLANK, `+` unless antisense allowed) ≥ 2 → best by mlen (first max); `h.strand != '+'` and
not antisense → `element_antisense`; tail = first non-None `hit(after_polyt(c))` over pside clips,
class ≠ h.class → `element_class_conflict`; colonies → `colony_mismatch`; else ("", pside).

### 4.7 absorb_one_sided (EV:1193; after the remap filters)

`tol = int(merge_tolerance_bp or 0)`; if `far_pair_strict`: `tol = max(tol, 5)`; tol ≤ 0 → no-op.
`one` = insertions with an open side, sorted by (-len(files), NAME). For each x still alive:
side = x's real side; target = `_absorb_target`: among alive insertions ≠ x on the same contig
whose open side ≠ side with `|junction(t, side) − junction(x, side)| ≤ tol`, the minimum of
`(t is one-sided, d, t.NAME)`. Skip if none, or both clips present and not
`clips_agree([t clip, x clip])`. Else: `members[t] = dedup(members[t] + members[x])`; t.member_sides
(copy or {}) gets `(side,)` for every x member not previously in members[t]; t.files += x's new
files; `new = evaluate(t, side)` (rows via the parent-side store lookup; §4.0 step 4); t's record
list = old list with the `side` record replaced by `new`; **every** record of t gets
`supported = 1 if n_independent >= min_ind else 0` (no sidecar-presence test here);
`records[t.name] = that list`; if `indel_aware_consensus`: `_replace_clips(t, [new])`; x dies.
Finally survivors keep order; for each dead x whose name is no survivor's name:
`records.pop(name)`.

---------------------------------------------------------------------------------------------

## 5. Remap (bowtie2 + samtools subprocesses; `sh -c`, return codes ignored like `os.system`)

### 5.1 Clean-remap filter (CI:160-201)

Write `STEM.fq.gz` (§6.1 text). If `STEM.bam` does not exist run
`{bt2} {STEM.fq.gz} -x {bowtie2_index} --end-to-end --sensitive --threads {N} --qc-filter | {samtools} view -F 4 -b -o {STEM.bam}`.
Read every record (`samtools view STEM.bam`): skip unmapped; `max_ins` = longest CIGAR `I`
(0 if none); `AS:i` (missing → -999); if `max_ins < clean_remap_max_insertion and AS >=
clean_remap_min_as` → filter name `qname[:-2]`. Drop insertions whose name is filtered.

### 5.2 Clipped-remap filter (CI:203-261)

If `STEM.insertionsonly.bam` does not exist: write the clip FASTQ (§6.2) to `STEM.fq.gz`, run
`{bt2} {STEM.fq.gz} -k 1000 -x {bowtie2_index2} --local --very-fast --threads {N} --qc-filter | {samtools} view -F 2308 -b -o {STEM.insertionsonly.bam}`.
Load the chain (`bowtie2_index2_lo`, §5.3). For each record: skip secondary (0x100) /
supplementary (0x800); if mapped: `contig, pos, side = qname.split(":")` (3 parts); contig gets a
`chr` prefix unless it has ≥ 3 chars starting `chr`; `L, R = pos.split("-")`; `p = int(R)` if side
`R` else `int(L)`; lift `(rname, reference_start if forward else reference_end)` (0-based start /
exclusive end = start + M/D/N/=/X lengths); any lifted result with chrom == contig and
`|p − lifted| < 1000` → filter `qname[:-2]`. Drop filtered insertions (by name).

### 5.3 Chain liftover — see `liftover.rs` (pyliftover 0.4.1 semantics, verified from its source).

---------------------------------------------------------------------------------------------

## 6. Output formats (decompressed bytes)

### 6.1 `combined.txt.gz` and the first `fq.gz` (CI:163-175, 269-281)

Five passes over the insertion list (final list for combined, post-evidence list for fq.gz):
FULL_INFO → `left_consensus.fastq("{name}:L")` + `right_consensus.fastq("{name}:R")`;
RIGHT_POLYA → `:L`; LEFT_POLYA → `:R`; RIGHT_DISC → `:L` (left consensus); LEFT_DISC → `:R`.

### 6.2 Clip FASTQ (2nd `fq.gz`, CI:206-218)

Per insertion in order: `lc.fastq("{name}:L")` if lc else `""`, then `rc.fastq("{name}:R")` if rc.
`trim_far_flank_before_remap` → lc/rc = `_far_flank_trimmed(ins, "L"/"R")` (CI:78: clip; far =
`str(left_aligned).upper()[:20]` for R, `revcomp(str(right_aligned.revcomp()).upper()[-20:])` for
L, `""` when that aligned is None; `len(far) < 20` → clip; `k = str(clip).upper().find(far)`;
`k >= 10` → `clip[:k]` else clip). Else raw `left_clipped` / `right_clipped`.

### 6.3 `genotyping.txt.gz` (CI:298-365)

For each final insertion (order): skip if name in the clipped-remap filter set; skip types 4/5;
`right = right_clipped.upper()[:max_bases]`, `left = left_clipped[:max_bases].upper().revcomp()`;
contig len > 5 → skip; `right_ref = get_sequence(contig, right_pos, right_pos+max_bases).upper()`
unless type 1; `left_ref = get_sequence(contig, left_pos-max_bases, left_pos).upper()` unless
type 2; exclusions in order: empty right_ref (type ≠ 1), empty left_ref (type ≠ 2), `N` in
right_ref, `N` in left_ref, type 1 and `sms([left_ref, left]) > 0`, type 2 and
`sms([right_ref, right]) > 0`, FULL and both sms > 0, FULL and (`left_ref == revcomp(right_ref)`
or `left_ref == right_ref`). Included:
`>{name}\n` + (type ≠ 1: `@RIGHT_INSERTION\n{right}\n@RIGHT_REFERENCE\n{right_ref}\n`) +
(type ≠ 2: `@LEFT_INSERTION\n{left}\n@LEFT_REFERENCE\n{left_ref}\n`). (`sms` on plain strings.)

### 6.4 `insertions.evidence.tsv.gz` (EV:546, EV:1307)

Header: the 26 names of `EVIDENCE_TSV_COLUMNS` (EV:49) tab-joined + `\n`. Rows for each name in
[final insertion names in order] + sorted(failed) (byte order), one row per record of
`records[name]` (absent → none; a repeated name repeats). Row = `str()` of, tab-joined:
`insertion_id, side, n_reads, n_fragments, n_independent, n_samples, n_mates, supported (NA|0|1),
clip_consensus, consensus_depth, polya_len_median ("" if None), polya_len_range ("" or "a-b"),
beyond_polya.lower(), beyond_polya_support if beyond else 0, polya_end (0|1), n_duplicates,
n_cross_sample_identical, fail_reason, consensus_stop, n_dup_coord, n_dup_seq, member_loci or ".",
n_short_used, n_short_rejected, n_short_mate_inside, n_independent_no_short`.
`clip = consensus.seq.lower()`; RIGHT: `clip_consensus = aligned + clip`, depth as is, beyond as
is; LEFT: `revcomp(clip) + aligned`, depth reversed, beyond revcomp'd. depth comma-joined.

### 6.5 `insertions.reads.fa.gz`

Same record order; per record `>{insertion_id}|{side}|{role}|{sample}|{frag}|{r12}\n{allele_forward_seq}\n`
for each kept row in row order (detached records were rendered with `rec.insertion_id` at
detach time, attached ones with the output name — identical in practice).

### 6.6 `combined.splice.tsv` (CI:37)

Inputs: `<f>.splice.tsv` of every input file (CI uses `input_files`, not accepted). Header line
skipped; rows split on `\t`, < 7 fields skipped, int parse failure skipped; `src[contig]` lists in
file then line order. Output header `insertion\tgene\tside\tn_exons\tintron_bp\tspan_bp\n`; for
each final insertion, for each src row of its contig: `pos = right_pos if side == "RIGHT" else
left_pos`; None → skip; `|bp − pos| <= 25` → `{name}\t{gene}\t{side}\t{n_exons}\t{intron}\t{span}\n`
(ints re-formatted by `str(int)`).

---------------------------------------------------------------------------------------------

## 7. Ordering inventory (Python dict/set iteration that reaches outputs)

| where | Python structure | Rust requirement |
|---|---|---|
| ISX full/polyA/disc maps | dict insertion order | insertion-ordered map (Vec + FxHashMap index) |
| ISX `_fuzzy_clusters`, `_fuzzy_one_sided` | explicit sorts on (…, NAME str) | sort on formatted names |
| ISX `_absorb` | `dict.fromkeys(mem)` | dedup preserving order |
| EV `_member_loci` | `dict.fromkeys` | same |
| EV collapse / independent_clusters | dict groups by first appearance; stable sorts | same |
| EV `_inside_mates` | dict keyed by fragment, first MATE wins | same |
| EV far-pair `clips` | `dict.fromkeys(c)` | dedup preserving order |
| EV `{l for _, l in ms}` | set → `sorted` | sort NAME strings |
| EV records / recmap | dict by name / id | name → uid map (§4.0 step 10, §4.7) |
| EV absorb `one`, `_absorb_target` | sorts / min on (…, NAME) | same keys |
| IC `_Col.w/.groups/.ins` | dict insertion order + first-max | ordered small maps |
| IC `_pick_seed` | stable sort, first 8 distinct | same |
| TF far_pair `pa` | dict order LEFT, RIGHT | same |
| failed names | `sorted(set)` | byte-order sort |
Sets used only for membership/cardinality (`have`, groups, colonies, filter names) are free.

---------------------------------------------------------------------------------------------

## 8. Fidelity risks

1. **mappy 2.31 vs minimap2 2.30** (the `minimap2` crate bundles 2.30; the reference venv has
   mappy 2.31; 2.31 fixed secondary/supplementary flagging and an inversion-alignment OOB). The
   matcher takes the max-mlen hit over ALL returned hits, so flag changes should not matter, but
   hit sets could differ in rare cases → affects slippage `_carries_element` and far-pair verdicts.
   Check the farm venv's `mappy.__version__`; if hits diverge, vendor minimap2 2.31 C sources.
2. edlib — eliminated (same C library, parity-tested).
3. Neumaier `sum` (Python ≥ 3.12 only). A pre-3.12 interpreter would produce naive sums; the
   reference venv and the farm are 3.12.
4. Consensus first-max ties over dict order (§4.4) — the most likely source of residual diffs.
5. bowtie2 output is identical for identical FASTQ input; the harness compares BAMs as sorted SAM.
6. Raw-token names: a discovery id with leading zeros would not round-trip (`Tok`); none exist.

---------------------------------------------------------------------------------------------

## 9. Module map

| Rust | Python | owner |
|---|---|---|
| `align.rs` | edlib.align | done |
| `pyfmt.rs` | builtins sum/round/int, DedupParams.budget, IC._median | done |
| `model.rs`, `context.rs` | Insertion attrs, tuples, ids | done (consensus helpers: P1) |
| `config.rs` | CONFIG | P1 |
| `seq.rs` | QS, revcomp, SC:89, TF:316-349 | P1 |
| `insertion.rs` | INS | P1 |
| `intersect.rs` | ISX | P1 |
| `region_filter.rs` | RF | P1 |
| `genome.rs` | GS / py2bit | P1 |
| `consensus.rs` | IC | P2 |
| `evidence/row.rs` | EV:80-148, 245-263, 1265 | P3 |
| `evidence/dedup.rs` | EV:178-519 | P3 |
| `evidence/junction.rs` | EV:523-755, 1301 | P3 |
| `tprt.rs` | TF:45-311 | P4 |
| `library.rs` | TF:171-229 | P4 |
| `evidence/filters.rs` | EV:358, 1294-1452 | P4 |
| `evidence/mod.rs` | EV:766-948, 1183-1298 | P5 |
| `evidence/store.rs` | EV:838-980 (purpose) | P5 |
| `evidence/output.rs` | EV:1307, CI:103 | P5 |
| `remap.rs` | CI:78-100, 160-261 | P6 |
| `liftover.rs` | pyliftover | P6 |
| `genotyping_out.rs` | CI:298-365 | P6 |
| `splice.rs` | CI:32-75 | P6 |
| `pipeline.rs`, `main.rs` | CI:113-365, main.py:148-188 | P6 |
