#!/bin/bash
# Equivalence harness: python combine_insertions (reference) vs the Rust peartree-combine binary
# on the same inputs and config, diffing EVERY output decompressed.
#
#   bash rust/peartree-combine/tests/equiv.sh                 # all variants
#   VARIANTS="tprt legacy" bash rust/peartree-combine/tests/equiv.sh
#   NO_BUILD=1 / THREADS=8 / THREAD_CHECK=1 (also run Rust with 1 thread and diff vs THREADS)
#
# Variants (fixture = 3 simulated colonies of the TPRT E2E, test/e2e/run_e2e.sh):
#   tprt         cluster/config.py.grch38.tprt, evidence sidecars       (production target)
#   tprt_more    tprt + count_short_overhang + keep_polya_one_sided      (SHORT reads, poly-A one-sided)
#   tprt_basecfg cluster/config.py.grch38 (legacy keys) WITH sidecars    (evidence defaults, no TPRT filters)
#   legacy       cluster/config.py.grch38, NO sidecars                   (legacy path)
# Every variant also gets a synthetic discovery splice sidecar for S1 (combined.splice.tsv).
#
# The python reference is cached per variant (stamp = config + python sources + inputs); it is
# re-run only when one of them changes. Exit status: 0 = every variant identical.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
CRATE="$(cd "$HERE/.." && pwd)"
REPO="$(cd "$CRATE/../.." && pwd)"
SP="${SP:-/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/262d6400-c197-4f82-b61f-06444781e46f/scratchpad/rustcombine}"
# fixture sources (read-only): the TPRT E2E work dir + genomes + a venv with edlib/mappy/pysam/
# pyliftover/py2bit. Regenerate with `bash test/e2e/run_e2e.sh` (SP=<dir>) if they are gone.
E2E_SP="${E2E_SP:-/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/fde0700f-e325-4651-8daf-0cdd52bd072b/scratchpad}"
SRC_E2E="${SRC_E2E:-$E2E_SP/work/e2e}"
GENOME_2BIT="${GENOME_2BIT:-$E2E_SP/genomes/hg38.2bit}"
BT2_INDEX="${BT2_INDEX:-$SRC_E2E/ref/reduced}"
CHAIN="${CHAIN:-$SRC_E2E/ref/identity.chain.gz}"
PY="${PY:-$E2E_SP/venv/bin/python}"
THREADS="${THREADS:-4}"
VARIANTS="${VARIANTS:-tprt tprt_more tprt_basecfg legacy}"
BIN="$CRATE/target/release/peartree-combine"
FIX="$SP/fixture"
WORK="$SP/equiv"
OUTPUTS="combined.txt.gz genotyping.txt.gz fq.gz insertions.evidence.tsv.gz insertions.reads.fa.gz combined.splice.tsv bam insertionsonly.bam"

die() { echo "equiv.sh: $*" >&2; exit 2; }
for f in "$PY" "$GENOME_2BIT" "$BT2_INDEX.1.bt2" "$CHAIN" "$SRC_E2E/S1.discovery.txt.gz"; do
    [ -e "$f" ] || die "missing fixture resource: $f"
done
command -v samtools >/dev/null || die "samtools not on PATH"
command -v bowtie2 >/dev/null || die "bowtie2 not on PATH"
"$PY" -c "import edlib, mappy, pysam, pyliftover, py2bit" || die "$PY lacks edlib/mappy/pysam/pyliftover/py2bit"

if [ -z "${NO_BUILD:-}" ]; then
    (cd "$CRATE" && cargo build --release 2>&1 | grep -E "^error|Finished" ) || die "cargo build failed"
fi

# ---------------------------------------------------------------- fixture
mkdir -p "$FIX/sidecar" "$FIX/nosidecar"
for i in 1 2 3; do
    s="$SRC_E2E/S$i.discovery.txt.gz"
    cp -p "$s" "$FIX/sidecar/"
    cp -p "$s.evidence.tsv.gz" "$FIX/sidecar/"
    cp -p "$s" "$FIX/nosidecar/"
done
# synthetic splice-hallmark sidecar for S1 (discovery splice_hallmark is off in the E2E): rows at
# the first 30 loci' junctions (+-0..30 bp so the 25 bp window both matches and misses)
"$PY" - "$FIX/sidecar/S1.discovery.txt.gz" > "$FIX/S1.splice.tsv" <<'EOF'
import gzip, sys
seen = []
with gzip.open(sys.argv[1], "rt") as fh:
    for k, line in enumerate(fh):
        if k % 4 == 0 and line.startswith("@"):
            locus = line[1:].rsplit(":", 2)[0]
            if locus not in seen:
                seen.append(locus)
        if len(seen) >= 30:
            break
print("contig\tbreakpoint\tside\tgene\tn_exons\tintron_bp\tspan_bp")
for n, locus in enumerate(seen):
    c, pos = locus.rsplit(":", 1)
    a, b = pos.split("-")
    side, tok = ("LEFT", a) if n % 2 else ("RIGHT", b)
    p = int(tok.split("_")[-1]) + (n % 4) * 10
    print(f"{c}\t{p}\t{side}\tGENE{n}\t{2 + n % 3}\t{100 * n}\t{150 + n}")
print("short\tline")
EOF
cp "$FIX/S1.splice.tsv" "$FIX/sidecar/S1.discovery.txt.gz.splice.tsv"
cp "$FIX/S1.splice.tsv" "$FIX/nosidecar/S1.discovery.txt.gz.splice.tsv"

variant_cfg() {   # name -> "base|fixture|overrides..."
    case "$1" in
        tprt)         echo "config.py.grch38.tprt|sidecar|" ;;
        tprt_more)    echo "config.py.grch38.tprt|sidecar|count_short_overhang=True keep_polya_one_sided=True" ;;
        tprt_basecfg) echo "config.py.grch38|sidecar|" ;;
        legacy)       echo "config.py.grch38|nosidecar|" ;;
        *) die "unknown variant $1" ;;
    esac
}

show_diff() {  # a b label
    echo "    first differences ($3):"
    diff <(cat "$1") <(cat "$2") | head -12 | sed 's/^/      /'
}

decomp() {  # file kind -> text on stdout
    case "$1" in
        *.bam) samtools view "$1" | LC_ALL=C sort ;;   # bowtie2 -p record order is not deterministic
        *.gz) gzip -dc "$1" ;;
        *) cat "$1" ;;
    esac
}

stamp_of() {  # cfgdir fixture -> stamp
    {
        cat "$1/config.py" "$REPO/cluster/$2"
        cat "$REPO"/src/combine_insertions*.py "$REPO"/src/indel_consensus.py "$REPO"/src/quality_seq.py \
            "$REPO"/src/revcomp.py "$REPO"/src/sequence_checks.py
        ls -l "$FIX/$3"
    } | shasum | cut -c1-16
}

FAILED=0
for v in $VARIANTS; do
    IFS='|' read -r base fx ovr <<< "$(variant_cfg "$v")"
    VD="$WORK/$v"; mkdir -p "$VD/cfg" "$VD/py" "$VD/rs"
    OV=(); for o in $ovr; do OV+=(--override "$o"); done
    "$PY" "$HERE/make_config.py" --base "$REPO/cluster/$base" --out "$VD/cfg" --genome-2bit "$GENOME_2BIT" \
        --bt2-index "$BT2_INDEX" --chain "$CHAIN" ${OV[@]+"${OV[@]}"} > /dev/null || die "make_config failed"
    FILES=("$FIX/$fx/S1.discovery.txt.gz" "$FIX/$fx/S2.discovery.txt.gz" "$FIX/$fx/S3.discovery.txt.gz")
    st="$(stamp_of "$VD/cfg" "$base" "$fx")"
    echo "== $v  (base $base, fixture $fx${ovr:+, overrides: $ovr})"
    if [ "$(cat "$VD/py/.stamp" 2>/dev/null)" != "$st" ]; then
        rm -rf "$VD/py"; mkdir -p "$VD/py"
        t0=$(date +%s)
        (cd "$REPO" && "$PY" test/e2e/run_combine.py --config-dir "$VD/cfg" --out-stem "$VD/py/P1" \
            --threads "$THREADS" "${FILES[@]}") > "$VD/py/combine.log" 2>&1 \
            || { echo "   python reference FAILED (see $VD/py/combine.log)"; FAILED=1; continue; }
        echo "$st" > "$VD/py/.stamp"
        echo "   python reference: $(( $(date +%s) - t0 )) s"
    else
        echo "   python reference: cached"
    fi
    run_rust() {  # outdir threads
        rm -rf "$1"; mkdir -p "$1"
        local t0=$(date +%s)
        (cd "$REPO" && PEARTREE_PYTHON="$PY" /usr/bin/time -l -o "$1/time.txt" "$BIN" --step combine_insertions --config "$VD/cfg/config.py" \
            --discovery_files "${FILES[@]}" --out "$1/P1" --threads "$2") > "$1/combine.log" 2>&1
        local rc=$?
        echo "   rust (threads $2): exit $rc, $(( $(date +%s) - t0 )) s, peak RSS $(awk '/maximum resident/ {printf "%.0f MB", $1/1048576}' "$1/time.txt" 2>/dev/null)"
        return $rc
    }
    run_rust "$VD/rs" "$THREADS" || { echo "   rust FAILED (tail of $VD/rs/combine.log):"; tail -5 "$VD/rs/combine.log" | sed 's/^/      /'; FAILED=1; continue; }
    vfail=0
    for o in $OUTPUTS; do
        a="$VD/py/P1.$o"; b="$VD/rs/P1.$o"
        if [ ! -e "$a" ] && [ ! -e "$b" ]; then continue; fi
        if [ ! -e "$a" ] || [ ! -e "$b" ]; then
            echo "   DIFF  P1.$o exists only in $([ -e "$a" ] && echo python || echo rust)"; vfail=1; continue
        fi
        decomp "$a" > "$VD/py.$o.txt"; decomp "$b" > "$VD/rs.$o.txt"
        if cmp -s "$VD/py.$o.txt" "$VD/rs.$o.txt"; then
            echo "   same  P1.$o ($(wc -l < "$VD/py.$o.txt" | tr -d ' ') lines)"
            rm -f "$VD/py.$o.txt" "$VD/rs.$o.txt"
        else
            echo "   DIFF  P1.$o"; show_diff "$VD/py.$o.txt" "$VD/rs.$o.txt" "$o"; vfail=1
        fi
    done
    [ -e "$VD/rs/P1.evidence_shards" ] && { echo "   DIFF  rust left P1.evidence_shards behind"; vfail=1; }
    if [ -n "${THREAD_CHECK:-}" ]; then
        run_rust "$VD/rs1" 1 || vfail=1
        for o in $OUTPUTS; do
            [ -e "$VD/rs/P1.$o" ] || continue
            cmp -s <(decomp "$VD/rs/P1.$o") <(decomp "$VD/rs1/P1.$o") || { echo "   DIFF  P1.$o: threads $THREADS vs 1"; vfail=1; }
        done
    fi
    [ $vfail -eq 0 ] && echo "   => IDENTICAL" || { echo "   => DIFFERENT"; FAILED=1; }
done
exit $FAILED
