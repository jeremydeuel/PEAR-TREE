# peartree-rte port plan

Branch `rte-rust`. Python reference: `tools/rte/*.py` at `eaa2718`. Contract: [SPEC.md](SPEC.md).

## Ground rules (every package)

* **File ownership is exclusive.** Edit only the files your package owns (table below). The
  FOUNDATION files are read-only for packages; if you need something there, implement it
  privately in your own file and say so in your report (the integration stage consolidates).
  `tests/golden_modules.rs` is shared: remove the `#[ignore]` line of YOUR tests only.
* **Signatures are fixed.** The stubs' public signatures and the shared types (in their
  foundation files) are the interface other packages compile against. Private helpers: free.
* **Done =** `cargo build --release` + `cargo test --release` pass, your golden tests are
  un-ignored and pass on the committed pytest golden AND the e2e golden
  (`PEARTREE_RTE_E2E=<manifest>`; default path is the foundation session's scratchpad), no
  `todo!()` left in your files, `cargo clippy` clean for your files. Commit with
  `git commit --no-gpg-sign`, never push.
* Port python semantics literally (SPEC §6: dict order, first-max ties, stable sort, `int()`
  truncation, `pyfmt` for sum/round/formatting). Floats are compared exactly.
* Thread safety: annotator runs loci in parallel; anything cached must be `Sync` (Mutex /
  OnceLock). `RteLibrary::aligner` and `SiteContext::aligner` are already lazily `Sync`.
* Quick reference of foundation APIs: `sequtil::{rc, polya_runs, trailing_polya_start,
  homopolymer_at, hamming, shannon, edlib_best, edlib_path, base_fraction, low_complexity}`,
  `mm::{Aligner::{from_path, from_seq, map}, MapOpts::{LIBRARY, LOCAL, SR}, Hit}`,
  `library::RteLibrary` (all python attributes, `aligner(AlignerKind)`, `source_for_flank`,
  `flank_strand`, `source_element`, `landmark_at`, `is_young`, `is_active_element`,
  `class_of`), `genome::{Genome, open_genome, FastaGenome, TwoBit}`, `config::*Cfg`,
  `record::{RteRecord, Detail}`, `pyfmt::*`, `transduction::known_source`.

## Status

| file | python | owner | status |
|---|---|---|---|
| align.rs, pyfmt.rs, packed.rs, sequtil.rs, mm.rs, genome.rs, config.rs, inputs.rs, stream.rs, record.rs, library.rs, io.rs, golden.rs, main.rs, lib.rs | inputs, record, sequtil, genome, library (+ edlib/mappy/py2bit) | FOUNDATION | done, tested |
| hallmarks.rs, score.rs | hallmarks, score | **WP-HALL** | done (golden) |
| assembly.rs | assembly | **WP-ASM** | done (golden) |
| structure.rs | structure | **WP-STRUCT** | done (golden) |
| transduction.rs (L1Rmsk, cons_identity, NovelSourceFinder), pseudogene.rs, genemodel.rs | transduction, pseudogene, annotate_v2.GeneModel subset | **WP-TD** | done (golden) |
| annotator.rs (+ main.rs wiring, annotate_v2 switch) | annotator | **WP-INT** | done: golden_annotate (445 events), golden_e2e_binary, annotate_v2 table byte-identical rust vs python (e2e, with/without GT reads) |
| — | locus_class.py, calibrate.py | stay python | — |

The four packages WP-HALL, WP-ASM, WP-STRUCT, WP-TD are independent (they only use foundation
code and each other's TYPES, which are already defined) and can run in parallel. WP-INT runs
after them.

---------------------------------------------------------------------------------------------

## WP-HALL — hallmarks + score  (simple → Sonnet)

* Owns: `src/hallmarks.rs`, `src/score.rs`.
* Port: hallmarks.py `split_junction` (27), `parse_locus` (37, regex `^(.+):(-?\d+)-(-?\d+)$`,
  greedy contig), `edge_run` (44), `polya_info` (84), `_find_near` (124: ALL overlapping exact
  occurrences via lookahead `finditer`, nearest to the hint, FIRST on ties; else
  `edlib_best(probe, win, max_frac=2/len(probe), both_strands=False)`; probe < 12 → None),
  `locate_site` (140, probe_len 30, slack 1500), `tsd_from_flanks` (158), `target_site` (172),
  `en_motif` (194), `en_bin` (209), `slippage_context` (221), `foldback` +
  `_longest_prefix_match` (247/262). score.py `WEIGHTS`, `THRESHOLDS`, `score` (incl. the hard
  cap, `l1dup`, `if v:` zero-weight skip, `round(py_sum(..), 2)`, points `f"{n}:{v:+g}"`).
* `edge_run` returns base 0 for python's `""`; `polya_info` ints vs floats: `la/rt` become
  `max(int, float)` → keep python's value for `length` (`float(la)`).
* Golden: `golden_hallmarks` (> 1,000 events incl. e2e hg38 2bit), `golden_score`.
* Effort: ~400 lines.

## WP-ASM — assembly  (complex → Opus)

* Owns: `src/assembly.rs` (types are there; implement `Assembler::layout`, `assemble`,
  `strand_from_polya` and private helpers).
* Port assembly.py: `_cut_hit_at` (159), `Assembler.layout` (178), `_gaps` (258), `_resolve`
  (270; sort by -score is STABLE over insertion order), `_edlib_targets_list` (306; consensus file
  order, cut at cons_end), `_rescue` (312; local first with the -0.02 bonus, then consensus, then
  flanks (file order) when `rescue_flanks`), `_mark_local` (347), `_local_is_element` (373),
  `_element_is_local` (392), `_merge_ref_runs` (423), `_smooth_polya` (450; restart-on-change
  loop), `_wide_local` (520), `assemble` (548; layouts: junction strings first in dict order,
  then reads), `AssemblyResult.build` (590), `_terminal_3p` (635), `_strand_from_polya` (676),
  `_strand_from_segments` (705), `_pileup` (730; votes per position, insertion votes,
  `most_common(1)` / `max(items)` FIRST max; `ins[t-1][""] -= 1`), `_nearest` (807; intact
  aligner, per-ctg first hit only, key (identity, blen), first max).
* mappy usage: `lib.aligner(Consensus|Flanks3|Flanks5|Intact)` (MapOpts::LIBRARY),
  `ctx.aligner()` (LOCAL), `ctx.wide_aligner()` (SR) — all done in foundation.
* Python aliasing to reproduce (see assembly.rs header): strand ≥ 0 sense layouts share Segment
  objects with raw_layouts, so `_pileup`'s updates appear in `raw_layouts`; strand < 0 not.
  Store layouts so that both lists in the result match the golden `out`.
* `Segment.matches` is `int(...)` truncation; `_resolve` recomputes proportional t-coords.
* Golden: `golden_assemble` (109 pytest + e2e events; compares the full AssemblyResult incl.
  raw/sense layouts, covered, covered_seqs, segments_on_cons, identities, class_bp order).
* Effort: ~900 lines. Hardest package; start with `layout` on single reads (the golden
  raw_layouts give per-read expected segments), then build/pileup/nearest.

## WP-STRUCT — structure.classify  (complex → Opus)

* Owns: `src/structure.rs`.
* Port structure.py: `classify` (106) and every helper (`_first_element_after_ref`,
  `_last_element_before_ref`, `_chain`, `_inversion`, `_element_after`, `_unexplained_tail`,
  `_template_abs`, `_foldback_5p`, `_local_templates`, `_sense_switch`, `_inverted_tail_5p`,
  `_flank_masked_frac`, `_source_class_ok`, `_at_rich`, `_td5_at_junction`,
  `_explained_by_consensus`, `_ref_at_breakpoint`). `known_source` is in transduction.rs
  (foundation, done).
* Callbacks are traits/closures (`SourceFinder`, `PremrnaFn`, `PseudogeneArg`) — the golden test
  replays python's recorded answers, so no dependency on WP-TD / the gene model.
* Identity vs equality: `exclude={id(x)}` → exclude by (layout, segment) index; `e in f5` →
  `Segment ==`. `call.detail` is an ordered `Detail` (insertion order = output order!).
  `call.tags.remove(t)` while iterating a copy; `detail.pop("premrna")`.
* Golden: `golden_classify` (109 pytest + e2e events; input AssemblyResult from the golden).
* Effort: ~700 lines.

## WP-TD — transduction + pseudogene + gene model  (moderate → Sonnet)

* Owns: `src/transduction.rs` (only `L1Rmsk`, `cons_identity`, `NovelSourceFinder::available`,
  `SourceFinder for NovelSourceFinder::find`; keep `known_source` / `MappyLocator` as they are),
  `src/pseudogene.rs`, `src/genemodel.rs`.
* Port: transduction.py `L1Rmsk.__init__` (both .out and rmsk.txt layouts, LINE/L1 ≥ min_len,
  strand C→-), `upstream_of`, `NovelSourceFinder.available` / `_l1_identity` (Mutex cache) /
  `find` (unique hit with mapq ≥ min_mapq, tag identity, rmsk → cohort_l1 → polymorphic_l1
  tier B, `round(ident, 4)`, detail string); `cons_identity` = tools/rte_library/common.py
  `cons_identity` → `identity(a, b, mode="HW")` (edlib path, `m / (m+x+i+d)`), max of both
  directions. pseudogene.py `load_exons_by_gene` (+`_merge`, gene → sorted merged exons),
  `load_gene_strands`, `ExonJunctionIndex.mrna/structure/cores/find` (thread-safe caches).
  annotate_v2.py `GeneModel._load/_resolve/_candidates/_splice_class/_genic_feature`
  (`genemodel.rs`; windows from `config::GeneModelCfg`).
* Golden: `golden_novel_find` (locator answers replayed; rmsk = fixture
  `novel_source.rmsk.out`), `golden_cons_identity`, `golden_exon_junctions`. The gene model has
  NO golden coverage (pytest uses a duck-typed model, the e2e has none): write unit tests whose
  expected values you generate by running annotate_v2's `GeneModel` in python on a small track
  (e.g. `test/fullstack` / `tools/build_gene_model.py` output) — put the expectations in the test.
* Effort: ~500 lines.

## WP-INT — annotator + switch-over  (complex → Opus; after the four above)

* Owns: `src/annotator.rs`, `src/main.rs`, `tests/golden_modules.rs` (whole), new
  `tests/golden_annotate.rs`, `tools/annotate_v2.py` (run_rte), `tools/rte/rust_bridge.py`,
  `cluster/*` wiring, SPEC.md / PORT_PLAN.md updates; may consolidate foundation helpers.
* Port annotator.py: `RteAnnotator.__init__` → `Resources::load` (library, genome/remap via
  `open_genome`, `MappyLocator` (remap_index), `L1Rmsk` (remap_rmsk), exon index
  (`exon_track()`, strands), gene model), `annotate` (all of SPEC §1.2), `_annotate_key`,
  `_premrna_fn` (window fetch ±rte_premrna_window, `Aligner::from_seq(.., MapOpts::SR)` built
  lazily per locus, `mlen/blen ≥ 0.9`, ≤ 250 bp → skip, gene model candidates/feature),
  `_pseudogene_structure_fn`, `_cap_reads` (= `stream::cap_reads`), `_score`,
  `_supported_tails`, `_beyond_polya`, `annotate_all` (stream::map_loci with `--chunk`, default
  e.g. 64; per-locus `catch_unwind` → error record), `cohort_source_pass` (needs `&mut`
  cohort_l1 before the re-run), `recurrence_pass`.
* Golden: `golden_e2e_binary` (binary over the e2e sidecars == python `final` rows + gt
  columns); add per-insertion `annotate` checks from the pytest golden where reconstructible
  (events whose test passed duck-typed objects — `_GeneModel`, mock libraries — may be skipped,
  report how many).
* Then: annotate_v2 `run_rte` → write inputs/config with rust_bridge, run the binary
  (path from config, e.g. `CONFIG['annotate']['rte_impl'] = 'rust'`, python default until
  validated), read records with `read_records`; set `exon_junction_proven`; keep `self.rte_lib`
  (RteLibrary, for locus_class); stream annotate_v2's own `_read_gt_fasta` too. Validate on the
  farm (PD49229) — Jeremy runs farm commands.
