#!/bin/bash
# =============================================================================
# PEAR-TREE — install / preflight check for a cluster run.
#
#   install.sh detect-assembly <BAM|CRAM>   print the assembly of a BAM's @SQ header
#   install.sh check   [--bam B | --assembly A] [--annotate]   verify, download nothing
#   install.sh install [--bam B | --assembly A] [--annotate]   verify and fetch what is missing
#   install.sh build-index <hs1|mm39>       download the fasta and bowtie2-build it (slow)
#   install.sh build-exons <asm> [--bam B]  build the Feature B exon model from Ensembl (~2 min)
#
# WHY THE ASSEMBLY MATTERS. Discovery and genotyping are assembly-agnostic — they
# only ever read the BAM. combine_insertions is NOT: it remaps clipped consensuses
# to the most complete assembly of the species (hs1 for human, mm39 for mouse) and
# then lifts those hits back onto the BAM's own coordinates with a pyliftover chain,
# and it pulls reference flanks from a 2bit of the BAM's assembly. Point either of
# those at the wrong assembly and you do not get an error — you get silently wrong
# coordinates. So the assembly is DETECTED FROM THE BAM HEADER and the matching
# 2bit + chain are asserted to exist before combine_insertions is allowed to run.
#
# Supported BAM assemblies: hg19 (== GRCh37/hs37d5), hg38, mm10, mm39.
# =============================================================================
set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PT_ROOT="$(cd "$SELF_DIR/.." && pwd)"

RES="${RES:-/lustre/scratch126/casm/teams/team273/users/jd43}"
VENV="${VENV:-$PT_ROOT/venv}"
BOWTIE2="${BOWTIE2:-/nfs/users/nfs_j/jd43/software/bowtie2-2.5.4-linux-x86_64/bowtie2}"
SAMTOOLS_MODULE="${SAMTOOLS_MODULE:-samtools-1.19}"
PYTHON_MODULE="${PYTHON_MODULE:-python/3.12.3}"
RUST_MODULE="${RUST_MODULE:-rust/1.87.0}"
UCSC="${UCSC:-https://hgdownload.soe.ucsc.edu/goldenPath}"

MISSING=0
FIX=0                       # 0 = check only, 1 = download/build

red()  { printf '\033[31m%s\033[0m\n' "$*"; }
grn()  { printf '\033[32m%s\033[0m\n' "$*"; }
ok()   { grn "  OK    $*"; }
bad()  { red  "  MISS  $*"; MISSING=$((MISSING+1)); }
warn() { printf '\033[33m  WARN  %s\033[0m\n' "$*"; }
die()  { red "$*" >&2; exit 1; }

# -----------------------------------------------------------------------------
# assembly fingerprints — chr1 LENGTH, not contig NAME.
#
# Length is the honest key: GRCh37 BAMs name it `1`, hg19 BAMs name it `chr1`, and
# both are the same assembly. hs37d5 == hg19 for 1..22,X,Y (it only adds decoys), so
# it takes the hg19 resources.
# -----------------------------------------------------------------------------
assembly_from_len() {
    case "$1" in
        249250621) echo hg19 ;;   # GRCh37 / hs37d5
        248956422) echo hg38 ;;   # GRCh38
        248387328) echo hs1  ;;   # T2T-CHM13v2.0 (a remap target, not a BAM assembly we support)
        195471971) echo mm10 ;;   # GRCm38
        195154279) echo mm39 ;;   # GRCm39
        *)         echo "" ;;
    esac
}

# Per-assembly plan. Each line: <2bit> <remap-index-assembly> <chain-basename-at-ucsc>
#   2bit   = genome_2bit        (reference flanks, in the BAM's coordinates)
#   remap  = bowtie2_index(2)   (best assembly of the species; clips are mapped here)
#   chain  = bowtie2_index2_lo  (remap-assembly -> BAM-assembly, reconciles the two)
remap_for()  { case "$1" in hg19|hg38) echo hs1 ;; mm10|mm39) echo mm39 ;; *) echo "" ;; esac; }
chain_ucsc() {
    case "$1" in
        hg19) echo "hs1/liftOver/hs1ToHg19.over.chain.gz" ;;
        hg38) echo "hs1/liftOver/hs1ToHg38.over.chain.gz" ;;
        mm10) echo "mm39/liftOver/mm39ToMm10.over.chain.gz" ;;
        mm39) echo "" ;;   # remap target IS the BAM assembly -> identity chain, synthesised
        *)    echo "" ;;
    esac
}
# Accept the file under either the UCSC name or the name this site already uses.
chain_candidates() {
    case "$1" in
        hg19) echo "$RES/hs1.hg19.all.chain.gz $RES/hs1ToHg19.over.chain.gz" ;;
        hg38) echo "$RES/pt_hu_trees/hs1.hg38.all.chain.gz $RES/hs1.hg38.all.chain.gz $RES/hs1ToHg38.over.chain.gz" ;;
        mm10) echo "$RES/mm39ToMm10.over.chain.gz $RES/mm39.mm10.all.chain.gz" ;;
        mm39) echo "$RES/mm39.identity.chain.gz" ;;
    esac
}
index_candidates() {
    case "$(remap_for "$1")" in
        hs1)  echo "$RES/pt_hu_trees/hs1/hs1 $RES/hs1/hs1" ;;
        mm39) echo "$RES/mm39/mm39" ;;
    esac
}
rmsk_ucsc() { case "$(remap_for "$1")" in hs1) echo "hs1/bigZips/hs1.repeatMasker.out.gz" ;; mm39) echo "mm39/bigZips/mm39.fa.out.gz" ;; esac; }
rmsk_candidates() {
    case "$(remap_for "$1")" in
        hs1)  echo "$RES/pt_hu_trees/hs1.repeatMasker.out.gz $RES/hs1.repeatMasker.out.gz" ;;
        mm39) echo "$RES/mm39.fa.out.gz" ;;
    esac
}

# -----------------------------------------------------------------------------
# Feature B (splice / processed-pseudogene) exon model.
#
# Only needed when the discovery config sets splice_hallmark=true — but then EVERY
# discovery job dies at startup ("cannot read exon_annotation") without it, so it is
# checked here rather than found out 174 array elements later.
#
# The model is queried with the MATE read's contig name, so its contigs must match the
# BAM's naming exactly. Ensembl is the right source: its GTFs use numeric names
# (1,2,..,X,Y), matching hs37d5/GRCh37 and GRCm38/39 BAMs directly; build-exons adds a
# `chr` prefix when the BAM is UCSC-style.
# -----------------------------------------------------------------------------
exon_candidates() {
    case "$1" in
        hg19) echo "$RES/grch37.exons.bed.gz $RES/hg19.exons.bed.gz" ;;
        hg38) echo "$RES/grch38.exons.bed.gz $RES/hg38.exons.bed.gz" ;;
        mm10) echo "$RES/grcm38.exons.bed.gz $RES/mm10.exons.bed.gz" ;;
        mm39) echo "$RES/grcm39.exons.bed.gz $RES/mm39.exons.bed.gz" ;;
    esac
}

# release-75 is the last Ensembl on GRCh37; 102 the last on GRCm38.
ensembl_gtf() {
    case "$1" in
        hg19) echo "https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz" ;;
        hg38) echo "https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz" ;;
        mm10) echo "https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz" ;;
        mm39) echo "https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz" ;;
    esac
}

first_present() { for c in $1; do [ -s "$c" ] && { echo "$c"; return 0; }; done; return 1; }
first_index()   { for c in $1; do [ -s "$c.1.bt2" ] || [ -s "$c.1.bt2l" ] && { echo "$c"; return 0; }; done; return 1; }

fetch() {  # fetch <url> <dest>
    local url="$1" dest="$2"
    mkdir -p "$(dirname "$dest")"
    echo "        downloading $url"
    if command -v curl >/dev/null; then curl -fL --retry 3 -o "$dest.part" "$url"
    else wget -O "$dest.part" "$url"; fi
    mv -f "$dest.part" "$dest"
}

# -----------------------------------------------------------------------------
# detect-assembly — the gate that protects combine_insertions
# -----------------------------------------------------------------------------
cmd_detect_assembly() {
    local BAM="${1:?usage: install.sh detect-assembly <BAM>}"
    [ -s "$BAM" ] || die "no such file: $BAM"
    command -v samtools >/dev/null || module load "$SAMTOOLS_MODULE" >/dev/null 2>&1 || true
    command -v samtools >/dev/null || die "samtools not on PATH (module load $SAMTOOLS_MODULE)"

    local hdr len asm
    hdr="$(samtools view -H "$BAM")" || die "cannot read header of $BAM"
    len="$(awk -F'\t' '/^@SQ/ {
              n=""; l="";
              for (i=2;i<=NF;i++) { if ($i ~ /^SN:/) n=substr($i,4); if ($i ~ /^LN:/) l=substr($i,4) }
              if (n=="chr1" || n=="1") { print l; exit }
          }' <<<"$hdr")"
    [ -n "$len" ] || die "no chr1/1 @SQ line in $BAM — cannot identify the assembly"
    asm="$(assembly_from_len "$len")"
    [ -n "$asm" ] || die "unrecognised assembly: chr1 length $len in $BAM
Supported: hg19/GRCh37 (249250621), hg38 (248956422), mm10 (195471971), mm39 (195154279)."
    echo "$asm"
}

# -----------------------------------------------------------------------------
# toolchain
# -----------------------------------------------------------------------------
check_tools() {
    echo "== toolchain"

    if command -v samtools >/dev/null || module load "$SAMTOOLS_MODULE" >/dev/null 2>&1; then
        ok "samtools ($(command -v samtools 2>/dev/null || echo "module $SAMTOOLS_MODULE"))"
    else bad "samtools — module load $SAMTOOLS_MODULE"; fi

    if [ -x "$BOWTIE2" ]; then ok "bowtie2 $BOWTIE2"
    elif command -v bowtie2 >/dev/null; then warn "bowtie2 on PATH but not at BOWTIE2=$BOWTIE2"
    else bad "bowtie2 — set BOWTIE2= or install bowtie2 >= 2.5"; fi

    for b in "$PT_ROOT/rust/peartree-discovery/target/release/peartree-discovery" \
             "$PT_ROOT/rust/peartree-genotype/target/release/peartree-genotype"; do
        if [ -x "$b" ]; then ok "$(basename "$b")"
        elif [ "$FIX" = 1 ]; then
            echo "        building rust binaries (cluster/build.sh)"
            bash "$PT_ROOT/cluster/build.sh" || die "build.sh failed"
            [ -x "$b" ] && ok "$(basename "$b")" || bad "$(basename "$b")"
        else bad "$(basename "$b") — run cluster/build.sh (module load $RUST_MODULE)"; fi
    done

    if [ -x "$VENV/bin/python" ]; then
        local miss=""
        for m in pysam pyliftover py2bit numpy pandas; do
            "$VENV/bin/python" -c "import $m" >/dev/null 2>&1 || miss="$miss $m"
        done
        if [ -z "$miss" ]; then ok "venv $VENV (pysam pyliftover py2bit numpy pandas)"
        elif [ "$FIX" = 1 ]; then
            "$VENV/bin/pip" install -q -r "$PT_ROOT/requirements.txt" && ok "venv deps installed" || bad "venv deps:$miss"
        else bad "venv missing:$miss — $VENV/bin/pip install -r requirements.txt"; fi
    elif [ "$FIX" = 1 ]; then
        echo "        creating venv at $VENV"
        module load "$PYTHON_MODULE" >/dev/null 2>&1 || true
        python3 -m venv "$VENV" && "$VENV/bin/pip" install -q --upgrade pip \
            && "$VENV/bin/pip" install -q -r "$PT_ROOT/requirements.txt" && ok "venv created"
    else bad "venv $VENV — module load $PYTHON_MODULE && python3 -m venv $VENV"; fi
}

# -----------------------------------------------------------------------------
# an identity chain, for when the remap target IS the BAM assembly (mm39)
#
# pyliftover always wants a chain, and there is no mm39->mm39 chain to download; a
# per-contig identity chain is the correct, exact answer (every base maps to itself).
# -----------------------------------------------------------------------------
make_identity_chain() {
    local asm="$1" dest="$2" sizes="$RES/$asm.chrom.sizes"
    [ -s "$sizes" ] || fetch "$UCSC/$asm/bigZips/$asm.chrom.sizes" "$sizes"
    awk 'BEGIN{OFS=" "} {print "chain", 1000, $1, $2, "+", 0, $2, $1, $2, "+", 0, $2, NR; print $2; print ""}' \
        "$sizes" | gzip -c > "$dest.part"
    mv -f "$dest.part" "$dest"
}

# -----------------------------------------------------------------------------
# per-assembly resources
# -----------------------------------------------------------------------------
check_assembly() {
    local asm="$1" want_annotate="$2"
    local remap; remap="$(remap_for "$asm")"
    [ -n "$remap" ] || die "unsupported assembly: $asm (want hg19|hg38|mm10|mm39)"

    echo "== assembly $asm  (combine_insertions remaps clips to $remap, then lifts back to $asm)"

    # 1. genome_2bit — reference flanks in the BAM's coordinates
    local twobit="$RES/$asm.2bit"
    if [ -s "$twobit" ]; then ok "genome_2bit $twobit"
    elif [ -s "$RES/pt_hu_trees/$asm.2bit" ]; then twobit="$RES/pt_hu_trees/$asm.2bit"; ok "genome_2bit $twobit"
    elif [ "$FIX" = 1 ]; then fetch "$UCSC/$asm/bigZips/$asm.2bit" "$twobit"; ok "genome_2bit $twobit"
    else bad "genome_2bit $twobit — $UCSC/$asm/bigZips/$asm.2bit"; fi

    # 2. bowtie2 index of the remap assembly (big; built, not downloaded)
    local idx
    if idx="$(first_index "$(index_candidates "$asm")")"; then ok "bowtie2_index $idx"
    else bad "bowtie2 index for $remap — none of: $(index_candidates "$asm")
        build it once with:  cluster/install.sh build-index $remap   (~3 h, 64 GB)"; fi

    # 3. the chain — the piece that makes the coordinates mean what you think
    local chain
    if chain="$(first_present "$(chain_candidates "$asm")")"; then ok "bowtie2_index2_lo $chain"
    elif [ "$FIX" = 1 ]; then
        local rel; rel="$(chain_ucsc "$asm")"
        if [ -n "$rel" ]; then
            chain="$RES/$(basename "$rel")"; fetch "$UCSC/$rel" "$chain"
        else
            chain="$RES/$asm.identity.chain.gz"
            echo "        synthesising identity chain ($remap == $asm, nothing to lift)"
            make_identity_chain "$asm" "$chain"
        fi
        ok "bowtie2_index2_lo $chain"
    else bad "chain $remap -> $asm — none of: $(chain_candidates "$asm")"; fi

    # 4. annotate (optional step 3)
    if [ "$want_annotate" = 1 ]; then
        echo "== annotate (optional)"
        local rmsk
        if rmsk="$(first_present "$(rmsk_candidates "$asm")")"; then ok "rmsk $rmsk"
        elif [ "$FIX" = 1 ]; then
            rmsk="$RES/$(basename "$(rmsk_ucsc "$asm")")"; fetch "$UCSC/$(rmsk_ucsc "$asm")" "$rmsk"; ok "rmsk $rmsk"
        else bad "rmsk for $remap — $UCSC/$(rmsk_ucsc "$asm")"; fi
        # rmsk must describe the REMAP assembly, not the BAM assembly: annotate reads
        # the clip mappings, which live in $remap space.
        [ -n "${rmsk:-}" ] && echo "        (note: this is $remap-space rmsk — correct for annotate; a $asm rmsk would be wrong)"

        for f in "$RES/Dfam-curated_only-hs.hmm" "$RES/dfamscan.pl"; do
            [ -s "$f" ] && ok "$(basename "$f")" || bad "$(basename "$f") — Dfam library / dfamscan.pl (see tools/annotate_v2.py)"
        done
        [ -s "$RES/Dfam-curated_only-hs.hmm.h3i" ] || warn "Dfam HMM not hmmpress'd — run: hmmpress $RES/Dfam-curated_only-hs.hmm"
        [ -x "$RES/hmmer-3.3.2/bin/nhmmscan" ] || command -v nhmmscan >/dev/null \
            || bad "nhmmscan — HMMER 3.x (hmmer.org)"
    fi

    if [ "$MISSING" -eq 0 ]; then
        echo
        grn "config['combine_insertions'] for $asm:"
        cat <<EOF
    'genome_2bit':       '$twobit',
    'bowtie2_index':     '$idx',
    'bowtie2_index2':    '$idx',
    'bowtie2_index2_lo': '$chain',
EOF
    fi
}

# -----------------------------------------------------------------------------
# check-config — the assertion pipeline.sh runs before combine_insertions.
#
# Detects the assembly from a real staged BAM and proves that src/config.py's
# genome_2bit and bowtie2_index2_lo actually belong to it. This is the check that
# catches the realistic mistake: copying config.py.grch37 for an hg38 run. Nothing
# downstream would complain — combine_insertions would happily lift hs1 hits through
# an hg19 chain onto hg38 breakpoints and emit plausible, wrong coordinates.
# Compares BASENAMES, so it does not care where you keep your references.
# Escape hatch for a deliberately unusual layout: PT_SKIP_ASSEMBLY_CHECK=1.
# -----------------------------------------------------------------------------
cmd_check_config() {
    local BAM="" CFG_ROOT="$PT_ROOT"
    while [ $# -gt 0 ]; do
        case "$1" in
            --bam) BAM="$2"; shift 2 ;;
            --pt-root) CFG_ROOT="$2"; shift 2 ;;
            *) die "check-config: unknown option $1" ;;
        esac
    done
    [ -n "$BAM" ] || die "usage: install.sh check-config --bam <BAM>"

    if [ "${PT_SKIP_ASSEMBLY_CHECK:-0}" = 1 ]; then
        warn "PT_SKIP_ASSEMBLY_CHECK=1 — skipping the assembly/config assertion"; return 0
    fi

    local asm; asm="$(cmd_detect_assembly "$BAM")"
    echo "  assembly of $(basename "$BAM"): $asm"

    [ -x "$VENV/bin/python" ] || die "no venv python at $VENV/bin/python"
    local vals
    vals="$("$VENV/bin/python" - "$CFG_ROOT" <<'PY'
import sys, os
sys.path.insert(0, os.path.join(sys.argv[1], 'src'))
from config import CONFIG
c = CONFIG['combine_insertions']
print(c['genome_2bit']); print(c['bowtie2_index2_lo']); print(c['bowtie2_index2'])
PY
)" || die "cannot import $CFG_ROOT/src/config.py (did you copy cluster/config.py.<asm> to src/config.py?)"

    local twobit chain idx
    twobit="$(sed -n 1p <<<"$vals")"; chain="$(sed -n 2p <<<"$vals")"; idx="$(sed -n 3p <<<"$vals")"

    # genome_2bit must BE this assembly
    if [ "$(basename "$twobit")" = "$asm.2bit" ]; then ok "genome_2bit is $asm ($twobit)"
    else bad "genome_2bit = $twobit, but the BAM is $asm — expected a file named $asm.2bit.
        Reference flanks would be read from the wrong assembly."; fi
    [ -s "$twobit" ] || bad "genome_2bit does not exist: $twobit"

    # the chain must lift the remap assembly onto THIS assembly
    local want_base="" c
    for c in $(chain_candidates "$asm"); do want_base="$want_base $(basename "$c")"; done
    if grep -qw -- "$(basename "$chain")" <<<"$want_base"; then ok "chain lifts $(remap_for "$asm") -> $asm ($chain)"
    else bad "bowtie2_index2_lo = $(basename "$chain"), but the BAM is $asm — expected one of:$want_base
        A wrong chain yields silently wrong coordinates, not an error."; fi
    [ -s "$chain" ] || bad "chain does not exist: $chain"

    [ -s "$idx.1.bt2" ] || [ -s "$idx.1.bt2l" ] || bad "bowtie2_index2 not built: $idx"

    [ "$MISSING" -eq 0 ] || die "
combine_insertions BLOCKED: config does not match the data ($MISSING problem(s) above).
Fix src/config.py (start from cluster/config.py.$asm if it exists) or run:
    cluster/install.sh install --bam $BAM"
    grn "  config matches the data ($asm) — combine_insertions may run."
}

# -----------------------------------------------------------------------------
# build-exons — one-off per assembly, ~2 min
#
# Builds the Feature B exon model: per-gene MERGED exons (the union of every
# transcript's exons collapsed into non-overlapping blocks), which is the right input
# for the "mates span >= N exons of one gene, intron skipped" signature — raw
# transcript exons overlap each other and blur it.
# Format (matches exons.rs): `contig begin end gene_id`, 0-based half-open, gzipped.
# -----------------------------------------------------------------------------
cmd_build_exons() {
    local asm="" bam=""
    while [ $# -gt 0 ]; do
        case "$1" in
            --bam) bam="${2:?--bam needs a path}"; shift 2 ;;
            *)     asm="$1"; shift ;;
        esac
    done
    [ -n "$asm" ] || [ -n "$bam" ] || die "usage: install.sh build-exons <hg19|hg38|mm10|mm39> [--bam <BAM>]"
    [ -n "$asm" ] || asm="$(cmd_detect_assembly "$bam")"

    local url out
    url="$(ensembl_gtf "$asm")"; [ -n "$url" ] || die "no Ensembl GTF known for assembly '$asm'"
    out="$(set -- $(exon_candidates "$asm"); echo "${1:-}")"
    [ -n "$out" ] || die "no exon-model path known for assembly '$asm'"
    [ -s "$out" ] && { echo "exon model already present: $out ($(zcat "$out" | wc -l) rows)"; return 0; }

    # Contig naming must match the BAM (the model is looked up by the mate's contig).
    # Ensembl is numeric; prefix with chr only if the BAM is UCSC-style.
    local prefix=""
    if [ -n "$bam" ]; then
        command -v samtools >/dev/null || module load "$SAMTOOLS_MODULE" >/dev/null 2>&1 || true
        if samtools view -H "$bam" 2>/dev/null | grep -qE '^@SQ.*SN:chr'; then prefix="chr"; fi
    fi

    local gtf="$RES/$(basename "$url")"
    [ -s "$gtf" ] || fetch "$url" "$gtf"

    echo "building per-gene merged exon model for $asm -> $out (contig prefix: '${prefix:-none}')"
    zcat "$gtf" \
      | awk -F'\t' '$3=="exon"{a=$9; sub(/.*gene_id "/,"",a); sub(/".*/,"",a); print a"\t"$1"\t"($4-1)"\t"$5}' \
      | sort -k1,1 -k2,2 -k3,3n \
      | awk -F'\t' -v p="$prefix" 'BEGIN{OFS="\t"}
            {if($1==g && $2==c && $3<=e){ if($4>e) e=$4 }
             else { if(g!="") print p c,s,e,g; g=$1; c=$2; s=$3; e=$4 }}
            END{ if(g!="") print p c,s,e,g }' \
      | gzip > "$out.part"
    mv -f "$out.part" "$out"
    echo "  $(zcat "$out" | wc -l) merged exon rows -> $out"
    echo "  point config.discovery.* at it:  exon_annotation = $out"
}

# check the exon model, but only if the discovery config actually turns Feature B on
check_exons() {
    local asm="$1"
    local cfg="${DISC_CFG:-$PT_ROOT/cluster/config.discovery.grch37}"
    [ -s "$cfg" ] || return 0
    grep -qiE '^[[:space:]]*splice_hallmark[[:space:]]*=[[:space:]]*(true|1)' "$cfg" || return 0

    echo "== Feature B exon model (splice_hallmark=true in $(basename "$cfg"))"
    local want
    want="$(grep -E '^[[:space:]]*exon_annotation[[:space:]]*=' "$cfg" | head -1 | sed 's/^[^=]*=[[:space:]]*//' | tr -d '[:space:]')"
    if [ -n "$want" ] && [ -s "$want" ]; then
        ok "exon_annotation -> $want ($(zcat "$want" 2>/dev/null | wc -l) rows)"
        return 0
    fi
    if [ "$FIX" = 1 ]; then
        cmd_build_exons "$asm" ${BAM_FOR_EXONS:+--bam "$BAM_FOR_EXONS"} && return 0
    fi
    bad "exon_annotation -> ${want:-<unset in config>}"
    echo "        every discovery job will die at startup without it. Build it with:"
    echo "        install.sh build-exons $asm --bam <a real BAM>"
    return 0
}

# -----------------------------------------------------------------------------
# build-index — one-off, expensive
# -----------------------------------------------------------------------------
cmd_build_index() {
    local asm="${1:?usage: install.sh build-index <hs1|mm39>}"
    case "$asm" in hs1|mm39) ;; *) die "build-index takes hs1 or mm39 (the remap targets)" ;; esac
    local dir="$RES/$asm" fa="$RES/$asm/$asm.fa"
    mkdir -p "$dir"
    [ -s "$fa" ] || { fetch "$UCSC/$asm/bigZips/$asm.fa.gz" "$fa.gz"; gunzip -f "$fa.gz"; }
    [ -x "$BOWTIE2" ] || die "need bowtie2 at $BOWTIE2"
    local cmd="${BOWTIE2}-build --threads 16 '$fa' '$dir/$asm'"
    if command -v bsub >/dev/null; then
        mkdir -p "$RES/logs"
        bsub -J "bt2build_$asm" -n 16 -q long -M 64000 \
             -R "select[mem>64000] rusage[mem=64000] span[hosts=1]" \
             -o "$RES/logs/bt2build_$asm.%J.log" -e "$RES/logs/bt2build_$asm.%J.err" "$cmd"
    else
        echo "no bsub — running locally (hours):"; eval "$cmd"
    fi
}

# -----------------------------------------------------------------------------
main() {
    local action="${1:-check}"; shift || true
    case "$action" in
        detect-assembly) cmd_detect_assembly "$@"; exit 0 ;;
        check-config)    cmd_check_config "$@"; exit 0 ;;
        build-index)     cmd_build_index "$@"; exit 0 ;;
        build-exons)     cmd_build_exons "$@"; exit 0 ;;
        install)         FIX=1 ;;
        check)           FIX=0 ;;
        -h|--help|help)  sed -n '2,20p' "$0"; exit 0 ;;
        *) die "unknown action: $action (check|install|detect-assembly|build-index|build-exons)" ;;
    esac

    local asm="" bam="" annotate=0
    while [ $# -gt 0 ]; do
        case "$1" in
            --assembly) asm="$2"; shift 2 ;;
            --bam)      bam="$2"; shift 2 ;;
            --annotate) annotate=1; shift ;;
            *) die "unknown option: $1" ;;
        esac
    done

    check_tools

    if [ -n "$bam" ]; then
        asm="$(cmd_detect_assembly "$bam")"
        echo "== detected assembly from $(basename "$bam"): $asm"
    fi
    if [ -z "$asm" ]; then
        warn "no --bam/--assembly given: checked the toolchain only, NOT the reference resources."
        warn "combine_insertions needs an assembly-matched 2bit + chain; re-run with --bam <a real BAM>."
    else
        check_assembly "$asm" "$annotate"
        # Feature B: no-op unless the discovery config enables splice_hallmark.
        BAM_FOR_EXONS="$bam" check_exons "$asm"
    fi

    echo
    if [ "$MISSING" -eq 0 ]; then grn "all present — ready to run."; exit 0
    elif [ "$FIX" = 1 ]; then red "$MISSING item(s) still missing after install (see above)."; exit 1
    else red "$MISSING item(s) missing — re-run as: cluster/install.sh install ${bam:+--bam $bam}${asm:+ --assembly $asm}"; exit 1; fi
}
main "$@"
