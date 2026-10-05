#!/usr/bin/env bash
# shellcheck disable=SC2015  # ok/warn/fail always return 0, so "test && ok || fail" is safe
# =============================================================================
# cluster/tprt/preflight.sh — head-node check before a TPRT A/B run of one patient.
#
#   bash cluster/tprt/preflight.sh PD37449           # full check (header reads of every BAM)
#   bash cluster/tprt/preflight.sh PD37449 --quick   # hard blockers only (run_ab.sh calls this)
#
# Reads only: BAM headers (samtools view -H / quickcheck — head-node safe and allowed on
# nst_links), config files, md5s. Submits nothing, writes nothing.
# Exit 0 = no FAIL (WARNs are listed); exit 1 = at least one FAIL.
#
#  1 checkouts   PT_ROOT_A worktree exists and is at the SAME commit as PT_ROOT_B
#  2 binaries    built, from this commit (build.stamp + mtime vs last rust/ commit), and the
#                binaries ACCEPT every key of both arms' Rust configs (probed by running them on
#                each config: config.rs ignores unknown keys with a warning, so a stale binary
#                would silently run arm B as arm A)
#  3 python      venv imports (pysam pyliftover py2bit numpy pandas edlib mappy scipy matplotlib),
#                src/config.py of each worktree imports and is the right arm, tools/phylo present
#  4 minimap2    binary present (only needed to build the hs1 index)
#  5 library     resources/rte_library md5 vs manifest.tsv
#  6 resources   hs1 gene model, hs1.2bit, hs1 .mmi (optional), discovery exon tracks,
#                combine/annotate references named in src/config.py
#  7 patient     tree tips vs GRCh38 WGS colonies (both directions)
#  8 BAMs        staged & quickcheck-clean & indexed, else header-readable in nst_links
#  9 assembly    install.sh check-config --bam <a real BAM> for BOTH worktrees
# =============================================================================
set -uo pipefail
# shellcheck source=cluster/tprt/common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

P="${1:?usage: preflight.sh <PATIENT_ID> [--quick]}"
QUICK=0; [ "${2:-}" = "--quick" ] && QUICK=1

NFAIL=0; NWARN=0
ok()   { printf '  OK    %s\n' "$*"; }
warn() { printf '  WARN  %s\n' "$*"; NWARN=$((NWARN+1)); }
fail() { printf '  FAIL  %s\n' "$*"; NFAIL=$((NFAIL+1)); }
hdr()  { printf '== %s\n' "$*"; }

PY="$VENV/bin/python"

# ---------------------------------------------------------------- 1 checkouts
hdr "checkouts"
HEAD_B="$(git -C "$PT_ROOT_B" rev-parse HEAD 2>/dev/null)" || HEAD_B=""
[ -n "$HEAD_B" ] && ok "arm B checkout $PT_ROOT_B @ ${HEAD_B:0:12} ($(git -C "$PT_ROOT_B" rev-parse --abbrev-ref HEAD))" \
    || fail "arm B checkout $PT_ROOT_B is not a git checkout"
if [ -d "$PT_ROOT_A" ]; then
    HEAD_A="$(git -C "$PT_ROOT_A" rev-parse HEAD 2>/dev/null || true)"
    if [ "$HEAD_A" = "$HEAD_B" ]; then ok "arm A worktree $PT_ROOT_A @ same commit"
    else fail "arm A worktree $PT_ROOT_A is at ${HEAD_A:0:12}, arm B at ${HEAD_B:0:12} — run: bash $TPRT_KIT_DIR/setup.sh worktree"; fi
else
    fail "arm A worktree missing: $PT_ROOT_A — run: bash $TPRT_KIT_DIR/setup.sh worktree"
fi
for R in "$PT_ROOT_B" "$PT_ROOT_A"; do
    [ -d "$R" ] || continue
    d="$(git -C "$R" status --porcelain --untracked-files=no 2>/dev/null | head -5)"
    [ -z "$d" ] || warn "uncommitted changes to tracked files in $R (the arms must differ by config only): $(printf '%s' "$d" | tr '\n' ' ')"
done

# ---------------------------------------------------------------- 2 binaries
hdr "binaries"
# Ask the binary itself which keys it does not know: it parses --config, prints
# "warning: ignoring unknown config key '<k>'" per unknown key, then exits on the missing --out /
# --bam (no output, no side effects). Grepping the binary for the key strings does NOT work:
# LLVM lowers short `match` string literals to immediate compares, so implemented keys are often
# absent from the binary's bytes.
probe_unknown_keys() {   # probe_unknown_keys <discovery|genotype> <binary> <cfg>
    local kind="$1" bin="$2" cfg="$3" tmp
    if [ "$kind" = discovery ]; then
        "$bin" --step discover --bam "$cfg" --config "$cfg" 2>&1 >/dev/null </dev/null
    else
        tmp="$(mktemp "${TMPDIR:-/tmp}/tprt_probe.XXXXXX")"
        "$bin" --step genotype --insertions "$tmp" --config "$cfg" 2>&1 >/dev/null </dev/null
        rm -f "$tmp"
    fi | sed -n "s/.*ignoring unknown config key '\([^']*\)'.*/\1/p" | sort -u
}
check_bin() {   # check_bin <label> <binary> <rust crate dir> <cfg...>
    local label="$1" bin="$2" crate="$3"; shift 3
    local kind="${label#peartree-}"
    if [ ! -x "$bin" ]; then fail "$label not built: $bin — run: bash $TPRT_KIT_DIR/setup.sh build"; return; fi
    local last_rust bin_mtime
    last_rust="$(git -C "$PT_ROOT_B" log -1 --format=%ct -- "$crate" 2>/dev/null || echo 0)"
    bin_mtime="$(stat -c %Y "$bin" 2>/dev/null || stat -f %m "$bin")"
    if [ "${bin_mtime:-0}" -lt "${last_rust:-0}" ]; then
        fail "$label is OLDER than the last commit touching $crate — rebuild: bash $TPRT_KIT_DIR/setup.sh build"
    else ok "$label newer than the last $crate commit"; fi
    local c unknown
    for c in "$@"; do
        unknown="$(probe_unknown_keys "$kind" "$bin" "$c" | paste -sd' ' -)"
        if [ -z "$unknown" ]; then ok "$label knows every key of $(basename "$c")"
        else fail "$label IGNORES keys of $(basename "$c"): $unknown — stale build; run: bash $TPRT_KIT_DIR/setup.sh build"; fi
    done
}
check_bin peartree-discovery "$DISCOVER_BIN" rust/peartree-discovery "$(arm_disc_cfg A)" "$(arm_disc_cfg B)"
check_bin peartree-genotype "$GENOTYPE_BIN" rust/peartree-genotype "$(arm_geno_cfg A)" "$(arm_geno_cfg B)"
if [ -s "$BUILD_STAMP" ]; then
    sc="$(awk -F'\t' '$1=="commit"{print $2}' "$BUILD_STAMP")"
    if [ "$sc" = "$HEAD_B" ]; then ok "build.stamp commit = HEAD"
    else
        changed="$(git -C "$PT_ROOT_B" diff --name-only "$sc" "$HEAD_B" -- rust 2>/dev/null | head -3 | paste -sd' ' -)"
        if [ -z "$changed" ] && git -C "$PT_ROOT_B" cat-file -e "$sc" 2>/dev/null; then ok "build.stamp commit ${sc:0:12} != HEAD, but rust/ is unchanged since"
        else fail "binaries were built at ${sc:0:12}, HEAD is ${HEAD_B:0:12} and rust/ changed ($changed) — run: bash $TPRT_KIT_DIR/setup.sh build"; fi
    fi
    for b in "$DISCOVER_BIN" "$GENOTYPE_BIN"; do
        [ -x "$b" ] || continue
        want="$(awk -F'\t' -v b="$b" '$1=="md5" && $3==b{print $2}' "$BUILD_STAMP")"
        [ "$want" = "$(md5_of "$b")" ] && ok "$(basename "$b") md5 matches build.stamp" \
            || fail "$(basename "$b") changed since setup.sh build (md5) — rebuild through setup.sh build"
    done
else
    warn "no $BUILD_STAMP (binaries not built through setup.sh build; commit provenance unknown)"
fi

# ---------------------------------------------------------------- 3 python
hdr "python ($PY)"
if [ -x "$PY" ]; then
    miss=""
    for m in pysam pyliftover py2bit numpy pandas edlib mappy scipy matplotlib; do
        "$PY" -c "import $m" >/dev/null 2>&1 || miss="$miss $m"
    done
    [ -z "$miss" ] && ok "imports: pysam pyliftover py2bit numpy pandas edlib mappy scipy matplotlib" \
        || fail "venv lacks:$miss — run: bash $TPRT_KIT_DIR/setup.sh venv"
    for ARM in A B; do
        R="$(arm_root "$ARM")"; [ -d "$R" ] || continue
        if [ ! -s "$R/src/config.py" ]; then fail "no $R/src/config.py — run: bash $TPRT_KIT_DIR/setup.sh configs"; continue; fi
        got="$(cd "$R" && TPRT_ROOT="$TPRT_ROOT" TPRT_RES="$TPRT_RES" "$PY" -c "import sys; sys.path.insert(0,'src'); import config; print(getattr(config,'TPRT_AB_ARM','none'), config.CONFIG['combine_insertions'].get('indel_aware_consensus'))" 2>&1 | tail -1)"
        case "$ARM:$got" in
            "A:A False"|"B:B True") ok "$R/src/config.py = arm $ARM (indel_aware_consensus=${got#* })" ;;
            *) fail "$R/src/config.py is not the arm-$ARM config (got: $got) — run: bash $TPRT_KIT_DIR/setup.sh configs" ;;
        esac
    done
    for t in tree_fit discrimination; do
        f="$PT_ROOT_B/tools/phylo/$t.py"
        if [ ! -s "$f" ]; then warn "tools/phylo/$t.py not in this checkout yet — the arms run, the evaluation job will fail until it lands (re-run evaluate.sh then)"
        elif (cd "$PT_ROOT_B" && "$PY" "$f" --help >/dev/null 2>&1); then ok "tools/phylo/$t.py --help"
        else fail "tools/phylo/$t.py --help fails in the venv (missing dependency?): $(cd "$PT_ROOT_B" && "$PY" "$f" --help 2>&1 | tail -1)"; fi
    done
else
    fail "no venv python $PY — run: bash $TPRT_KIT_DIR/setup.sh venv"
fi

# ---------------------------------------------------------------- 4 minimap2
hdr "minimap2"
if [ -x "$MINIMAP2" ]; then ok "minimap2 $("$MINIMAP2" --version 2>/dev/null) ($MINIMAP2)"
elif command -v minimap2 >/dev/null 2>&1; then ok "minimap2 $(minimap2 --version 2>/dev/null) ($(command -v minimap2))"
else warn "no minimap2 (only needed to build $TPRT_RES/hs1.sr.mmi) — bash $TPRT_KIT_DIR/setup.sh minimap2"; fi

# ---------------------------------------------------------------- 5 library
hdr "resources/rte_library (md5 vs manifest.tsv)"
for R in "$PT_ROOT_B" "$PT_ROOT_A"; do
    L="$R/resources/rte_library"; [ -d "$L" ] || { [ "$R" = "$PT_ROOT_A" ] && continue; fail "no $L"; continue; }
    bad=""; n=0
    while IFS=$'\t' read -r f _rec _bytes md5; do
        [ "$f" = file ] && continue; [ -n "$f" ] || continue
        n=$((n+1))
        if [ ! -s "$L/$f" ]; then bad="$bad $f(missing)"
        elif [ "$(md5_of "$L/$f")" != "$md5" ]; then bad="$bad $f(md5)"; fi
    done < "$L/manifest.tsv"
    [ -z "$bad" ] && ok "$n files match ($L)" || fail "$L:$bad"
done

# ---------------------------------------------------------------- 6 resources
hdr "reference resources"
gm="$TPRT_RES/hs1.gene_model.tsv.gz"
if [ -s "$gm" ] && gzip -t "$gm" 2>/dev/null; then
    gmrows="$(gzip -dc "$gm" | wc -l | tr -d ' ')"
    if [ "$gmrows" -gt 1000 ]; then ok "hs1 gene model $gm ($gmrows rows)"
    else fail "hs1 gene model $gm has only $gmrows rows — run: bash $TPRT_KIT_DIR/setup.sh resources"; fi
elif [ -e "$gm" ]; then fail "hs1 gene model $gm is not valid gzip — run: bash $TPRT_KIT_DIR/setup.sh resources (rebuilds it)"
else fail "no hs1 gene model $gm — run: bash $TPRT_KIT_DIR/setup.sh resources"; fi
[ -s "$TPRT_RES/hs1.2bit" ] && ok "hs1.2bit $TPRT_RES/hs1.2bit" || fail "no $TPRT_RES/hs1.2bit — run: bash $TPRT_KIT_DIR/setup.sh resources"
[ -s "$TPRT_RES/hs1.sr.mmi" ] && ok "hs1 minimap2 index $TPRT_RES/hs1.sr.mmi" \
    || warn "no $TPRT_RES/hs1.sr.mmi — annotate's novel-source locator is off (known sources still work); bash $TPRT_KIT_DIR/setup.sh index"
for ARM in A B; do
    c="$(arm_disc_cfg "$ARM")"
    if grep -qiE '^[[:space:]]*splice_hallmark[[:space:]]*=[[:space:]]*(true|1)' "$c"; then
        ex="$(grep -E '^[[:space:]]*exon_annotation[[:space:]]*=' "$c" | head -1 | sed 's/^[^=]*=[[:space:]]*//; s/[[:space:]]*#.*//')"
        [ -s "$ex" ] && ok "arm $ARM discovery exon track $ex" || fail "arm $ARM: $(basename "$c") enables splice_hallmark but exon_annotation is missing: $ex (every discovery task exits 1)"
    fi
done
if [ -x "$PY" ] && [ -s "$PT_ROOT_B/src/config.py" ]; then
    for ARM in A B; do
        R="$(arm_root "$ARM")"; [ -s "$R/src/config.py" ] || continue
        paths="$(cd "$R" && TPRT_ROOT="$TPRT_ROOT" TPRT_RES="$TPRT_RES" "$PY" - <<'PY' 2>/dev/null
import sys; sys.path.insert(0, 'src'); from config import CONFIG
c, a = CONFIG['combine_insertions'], CONFIG['annotate']
for k in ('genome_2bit', 'bowtie2_index2_lo', 'bowtie2_executable', 'samtools_executable'):
    print(f"combine.{k}\t{c[k]}")
print(f"combine.bowtie2_index2\t{c['bowtie2_index2']}.1.bt2")
for k in ('rmsk', 'hmm', 'dfamscan', 'hmmer', 'genome_2bit', 'remap_rmsk', 'exon_annotation', 'remap_2bit'):
    if a.get(k):
        print(f"annotate.{k}\t{a[k]}")
print(f"annotate.tmp_dir\t{a['tmp']('x')('p').rsplit('/', 1)[0]}")
PY
)"
        while IFS=$'\t' read -r k v; do
            [ -n "$k" ] || continue
            if [ "$k" = annotate.tmp_dir ]; then
                [ -d "$v" ] && ok "arm $ARM $k $v" || fail "arm $ARM $k $v does not exist — run: bash $TPRT_KIT_DIR/setup.sh configs"
            elif [ "$k" = combine.bowtie2_index2 ]; then
                { [ -s "$v" ] || [ -s "${v%.bt2}.bt2l" ]; } && ok "arm $ARM $k ${v%.1.bt2}" || fail "arm $ARM $k not built: ${v%.1.bt2}"
            else
                [ -e "$v" ] && ok "arm $ARM $k $v" || fail "arm $ARM $k missing: $v"
            fi
        done <<<"$paths"
    done
fi

# ---------------------------------------------------------------- 7 patient
hdr "patient $P"
PDIR="$(patient_dir "$P")" || exit 1
TREE="$(patient_tree "$P")" || exit 1
ok "patient dir $PDIR, tree $(basename "$TREE")"
SAMP="$(patient_samples "$P")" || { fail "fleet.sh samples $P"; SAMP=""; }
[ -n "$SAMP" ] || fail "$P has no GRCh38 WGS colony in colonies.tsv"
tips="$(tree_tips "$TREE")"
t_not_s="$(comm -23 <(printf '%s\n' "$tips") <(printf '%s\n' "$SAMP" | cut -f1 | sort -u))"
s_not_t="$(comm -13 <(printf '%s\n' "$tips") <(printf '%s\n' "$SAMP" | cut -f1 | sort -u))"
nt="$(printf '%s\n' "$tips" | grep -c .)"; ns="$(printf '%s\n' "$SAMP" | grep -c .)"
if [ -z "$t_not_s" ]; then ok "all $nt tree tips are GRCh38 WGS colonies ($ns colonies)"
else warn "$(printf '%s\n' "$t_not_s" | grep -c .) of $nt tree tips have no GRCh38 WGS colony (tree_fit sees them as missing): $(printf '%s' "$t_not_s" | tr '\n' ' ')"; fi
[ -z "$s_not_t" ] || warn "$(printf '%s\n' "$s_not_t" | grep -c .) colonies are not tree tips (run, but unplaceable on the tree): $(printf '%s' "$s_not_t" | tr '\n' ' ')"

# ---------------------------------------------------------------- 8 BAMs + 9 assembly
GATE_BAM=""
if [ "$QUICK" = 0 ] && [ -n "$SAMP" ]; then
    hdr "BAMs ($ns colonies; staged -> quickcheck + index, else nst_links header)"
    if load_samtools; then
        n_st=0; n_irr=0; n_iro=0; n_nst=0; n_none=0; n_noidx=0
        have_irods=0; load_irods && have_irods=1
        [ "$have_irods" = 1 ] || warn "iquest not available (module load IRODS) — falling back to nst_links header reads"
        while IFS=$'\t' read -r s pj; do
            b="$(staged_bam "$s" "$pj")"
            if [ -s "$b" ]; then
                if samtools quickcheck "$b" 2>/dev/null; then
                    n_st=$((n_st+1)); GATE_BAM="${GATE_BAM:-$b}"
                    [ -s "$b.bai" ] || [ -s "${b%.bam}.bai" ] || n_noidx=$((n_noidx+1))
                else fail "$s: staged BAM fails samtools quickcheck (truncated?): $b — delete it so the run re-stages"; fi
            elif [ "$have_irods" = 1 ] && rb="$(irods_readable_bam "$s" "$pj")"; then
                n_irr=$((n_irr+1)); GATE_BAM="${GATE_BAM:-$rb}"
            elif [ "$have_irods" = 1 ] && [ -n "$(irods_replicas "$s" "$pj")" ]; then
                n_iro=$((n_iro+1))   # in iRODS but no replica header-readable from this node: stageBam.pl still works
            elif [ -d "$NST" ] && samtools view -H "$(nst_bam "$s" "$pj")" >/dev/null 2>&1; then
                n_nst=$((n_nst+1)); GATE_BAM="${GATE_BAM:-$(nst_bam "$s" "$pj")}"
            else
                n_none=$((n_none+1)); fail "$s: not staged and not found in iRODS as /cgp/intproj/$pj/sample/$s/$s.*sample.dupmarked.bam — check colonies.tsv project id"
            fi
        done <<<"$SAMP"
        ok "$n_st staged (quickcheck OK; $n_noidx without .bai — genotype indexes on demand), $n_irr in iRODS with a readable replica, $n_iro in iRODS (no replica readable here; staging still works), $n_nst via nst_links, $n_none missing"
        n_nst=$((n_nst + n_irr + n_iro))   # everything not yet staged, for the footprint note
        [ "$n_st" -gt 0 ] && [ "$n_nst" -eq 0 ] || note "staging footprint: ~$(( (n_nst + n_none) * 45 )) GB to stage into $STAGING_ROOT (45 GB/BAM estimate)"
    else
        fail "samtools not available (module load $SAMTOOLS_MODULE)"
    fi
    hdr "assembly gate (install.sh check-config, both worktrees)"
    if [ -n "$GATE_BAM" ]; then
        for ARM in A B; do
            R="$(arm_root "$ARM")"; [ -s "$R/src/config.py" ] || continue
            if out="$(TPRT_ROOT="$TPRT_ROOT" TPRT_RES="$TPRT_RES" VENV="$VENV" bash "$PT_ROOT_B/cluster/install.sh" check-config --bam "$GATE_BAM" --pt-root "$R" 2>&1)"; then
                ok "arm $ARM: $(printf '%s\n' "$out" | grep -E 'assembly of' | sed 's/^ *//') — config matches"
            else
                fail "arm $ARM check-config:"; printf '%s\n' "$out" | sed 's/^/        /'
            fi
        done
    elif [ "${n_iro:-0}" -gt 0 ]; then
        warn "assembly gate deferred: colonies are in iRODS but no replica is header-readable from $(hostname) — re-run preflight after the first BAM is staged (run_ab.sh stages them)"
    else
        fail "no readable BAM to run the assembly gate against"
    fi
fi

echo
if [ "$NFAIL" -eq 0 ]; then echo "preflight $P: PASS ($NWARN warning(s))"; exit 0
else echo "preflight $P: $NFAIL FAIL(s), $NWARN warning(s)"; exit 1; fi
