# PEAR-TREE - paired ends of aberrant retrotransposons in phylogenetic trees
#
# Copyright (C) 2025 Jeremy Deuel <jeremy.deuel@usz.ch>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.


import gzip
import os
from revcomp import revcomp
import pysam
import pyliftover
from config import CONFIG
from combine_insertions_insertion import Insertion, TYPE_LEFT_POLYA, TYPE_RIGHT_POLYA, TYPE_FULL_INFO, TYPE_LEFT_DISC, TYPE_RIGHT_DISC
from combine_insertions_intersect_insertions import intersect_insertions
from combine_insertions_region_filter import filter_dense_regions
from combine_insertions_get_sequence import get_sequence
from sequence_checks import sequence_matching_score
from collections import Counter

def _splice_sidecar(path):
    """`<combined>.txt.gz` -> `<combined>.splice.tsv` (shared by combine + annotate)."""
    return (path[:-7] if path.endswith(".txt.gz") else path) + ".splice.tsv"


def write_combined_splice(input_files, insertions, combined_insertions, window=25):
    """Re-key discovery's per-file `<file>.splice.tsv` sidecars (Feature B / splice_hallmark:
    mate reads spanning >= splice_min_exons exons of one gene, introns skipped) onto the
    surviving COMBINED insertions, by the combined insertion name, so stage-4 annotate can
    attach the processed-pseudogene evidence. Writes `<combined>.splice.tsv`. No-op (writes
    nothing) when no discovery splice sidecar exists (splice_hallmark off) -> fully backward
    compatible."""
    src = {}  # contig -> [(breakpoint, side, gene, n_exons, intron_bp, span_bp)]
    found = False
    for f in input_files:
        sp = f + ".splice.tsv"
        if not os.path.exists(sp):
            continue
        found = True
        with open(sp) as fh:
            fh.readline()  # header: contig breakpoint side gene n_exons intron_bp span_bp
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 7:
                    continue
                try:
                    src.setdefault(p[0], []).append((int(p[1]), p[2], p[3], int(p[4]), int(p[5]), int(p[6])))
                except ValueError:
                    continue
    if not found:
        return
    out_path = _splice_sidecar(combined_insertions)
    n = 0
    with open(out_path, "w") as out:
        out.write("insertion\tgene\tside\tn_exons\tintron_bp\tspan_bp\n")
        for i in insertions:
            for (bp, side, gene, nex, intron, span) in src.get(i.reference_name, []):
                pos = i.right_pos if side == "RIGHT" else i.left_pos
                if pos is None:
                    continue
                if abs(bp - pos) <= window:
                    out.write(f"{i.name}\t{gene}\t{side}\t{nex}\t{intron}\t{span}\n")
                    n += 1
    print(f"aggregated {n} discovery splice-hallmark rows -> {out_path}")


def _far_flank_trimmed(ins, side, probe=20, min_insert=10):
    """Clip of `side` for the clipped-remap filter, cut where it runs into the OTHER
    junction's reference flank. A short insertion is spanned completely by junction reads, so
    the (indel-aware) clip = insert + far flank; that flank remaps right next to the breakpoint
    and the filter would discard a real insertion. RIGHT clip (outward = reference-forward):
    insert + ref[L:]; LEFT clip (outward): rc(insert) + rc(ref[..R]). The far flank is taken
    from the junction records themselves (left_aligned starts at L, right_aligned ends at R)."""
    clip = ins.right_clipped if side == "R" else ins.left_clipped
    if clip is None:
        return None
    try:
        if side == "R":
            far = str(ins.left_aligned).upper()[:probe] if ins.left_aligned is not None else ""
        else:
            far = revcomp(str(ins.right_aligned.revcomp()).upper()[-probe:]) if ins.right_aligned is not None else ""
    except Exception:
        return clip
    if len(far) < probe:
        return clip
    k = str(clip).upper().find(far)
    if k >= min_insert:
        return clip[:k]
    return clip


def _evidence_paths(combined_insertions):
    """`<stem>.combined.txt.gz` -> (`<stem>.insertions.evidence.tsv.gz`, `<stem>.insertions.reads.fa.gz`)."""
    stem = combined_insertions
    for suffix in (".combined.txt.gz", ".txt.gz"):
        if stem.endswith(suffix):
            stem = stem[:-len(suffix)]
            break
    return f"{stem}.insertions.evidence.tsv.gz", f"{stem}.insertions.reads.fa.gz"


def combine_insertions(input_files, insertions_genotyping_file, combined_insertions, insertions_fasta, insertion_bam, threads,
                       evidence_tsv=None, reads_fa=None):

    all_insertions = []
    accepted_files = []
    j_cutoff = 0
    for f in input_files:
        j_cutoff += 1
        #if j_cutoff > 5: break
        fb = os.path.basename(f)
        file_insertions = [i for i in Insertion.parseFile(f) if len(i.reference_name) < 6 and i.reference_name != "MT" and i.reference_name !="chrM"]
        print(f"\033[37mFile {fb}: imported {len(file_insertions)} insertions\033[0m")
        if len(file_insertions) > CONFIG['combine_insertions']['exclude_files_with_many_insertions']:
            print(f"\033[31mremoved file {f} since it contains too many insertions.\033[0m")
        else:
            all_insertions += file_insertions
            accepted_files.append(f)
    print(f"intersecting insertions from {len(input_files)} files...")
    # per-sample discovery breakpoints, taken before intersect merges records (far-pair
    # colony-consistency test, CONFIG['combine_insertions']['far_pair_strict'])
    breakpoints = None
    if CONFIG['combine_insertions'].get('far_pair_strict', False):
        from combine_insertions_evidence import discovery_breakpoints
        breakpoints = discovery_breakpoints(all_insertions)
    insertions = intersect_insertions(all_insertions)
    #remove insertions in regions with far too high count
    bin_range = 100
    ins_cutoff = 4
    insertions, removed, regions = filter_dense_regions(insertions, bin_range, ins_cutoff)
    print(f"filtering regions with very high insertion rate of {ins_cutoff} or higher per {bin_range} bases ,removed {removed} insertions in {regions} regions, {len(insertions)} insertions are remaining")
    # TPRT-hallmarks: pooled per-patient junction evidence from discovery's optional
    # `<sample>.evidence.tsv.gz` sidecars (independent-fragment gate + indel-aware clip
    # consensus). Returns None when no sidecar exists -> legacy behaviour, byte-identical.
    # Imported lazily so the legacy path needs neither the module nor its edlib dependency.
    evidence = None
    # sidecar is `<sample>.txt.gz.evidence.tsv.gz` (Rust) or `<sample>.evidence.tsv.gz`
    if any(os.path.exists(f + ".evidence.tsv.gz")
           or os.path.exists((f[:-7] if f.endswith(".txt.gz") else f) + ".evidence.tsv.gz")
           for f in accepted_files):
        from combine_insertions_evidence import apply_evidence
        evidence = apply_evidence(insertions, accepted_files, CONFIG['combine_insertions'],
                                  breakpoints=breakpoints)
    if evidence is not None:
        insertions, evidence_records, evidence_failed, _ = evidence
    print(f"writing summarised insertions fasta file {insertions_fasta}")
    # compresslevel=1: this is a scratch file bowtie2 reads back immediately, so the
    # default level 9 spends CPU shrinking bytes nothing keeps.
    with gzip.open(insertions_fasta, 'wt', compresslevel=1) as f:
        f.writelines(
            [f'{i.left_consensus.fastq(f"{i.name}:L")}{i.right_consensus.fastq(f"{i.name}:R")}' for i in insertions if i.type is TYPE_FULL_INFO])
        f.writelines(
            [f'{i.left_consensus.fastq(f"{i.name}:L")}' for i in insertions if i.type is TYPE_RIGHT_POLYA])
        f.writelines(
            [f'{i.right_consensus.fastq(f"{i.name}:R")}' for i in insertions if i.type is TYPE_LEFT_POLYA])
        # Feature A: a discordant-anchored call has one real side; remap that consensus
        # so the clean-remap filter can still reject it if it aligns contiguously.
        f.writelines(
            [f'{i.left_consensus.fastq(f"{i.name}:L")}' for i in insertions if i.type is TYPE_RIGHT_DISC])
        f.writelines(
            [f'{i.right_consensus.fastq(f"{i.name}:R")}' for i in insertions if i.type is TYPE_LEFT_DISC])
    if not os.path.exists(insertion_bam):
        print(f"running bowtie2 {CONFIG['combine_insertions']['bowtie2_executable']} with index {CONFIG['combine_insertions']['bowtie2_index']}")
        cmd = f"{CONFIG['combine_insertions']['bowtie2_executable']} {insertions_fasta} -x {CONFIG['combine_insertions']['bowtie2_index']} --end-to-end --sensitive --threads {threads} --qc-filter | {CONFIG['combine_insertions']['samtools_executable']} view -F 4 -b -o {insertion_bam}"
        os.system(cmd)
    # A consensus (aligned + clipped) counts as "maps entirely to the reference" only if
    # it aligns end-to-end CLEANLY. A genuine insertion junction can be forced end-to-end
    # against a reference that lacks the insertion only by opening an insertion (the
    # element) >= the clip length, at a poor alignment score; a real assembly-discordance /
    # reference-contiguous read aligns with no such gap and a near-perfect score. Requiring
    # a clean alignment stops us discarding full-length Alu/SVA and 3'-transduction junctions
    # (their consensus otherwise force-aligns with a big I to a paralogous copy).
    max_clean_ins = CONFIG['combine_insertions'].get('clean_remap_max_insertion',
                                                     CONFIG['discovery']['min_clip_len'])
    min_clean_as = CONFIG['combine_insertions'].get('clean_remap_min_as', -15)
    filter_reads = set()
    with pysam.AlignmentFile(insertion_bam) as f:
        for read in f:
            if read.is_unmapped:
                continue
            max_ins = max((length for op, length in (read.cigartuples or []) if op == pysam.CINS), default=0)
            align_score = read.get_tag('AS') if read.has_tag('AS') else -999
            if max_ins < max_clean_ins and align_score >= min_clean_as:
                filter_reads.add(read.query_name[:-2])
    print(
        f"detected {len(filter_reads)} insertions where at least one end maps cleanly (no inserted block) to the reference genome, removing these (since they can not be chimeric)...")
    insertions = [i for i in insertions if i.name not in filter_reads]

    print(f"now re-mapping in local mode all clipped parts of reads")
    clipped_bam = insertion_bam[:-4] + ".insertionsonly.bam"
    if not os.path.exists(clipped_bam):
        with gzip.open(insertions_fasta, 'wt', compresslevel=1) as f:
            # Feature A: a discordant-anchored call has a None clipped side; emit only the
            # side(s) that carry sequence (full-info calls still emit both).
            if CONFIG['combine_insertions'].get('trim_far_flank_before_remap', False):
                # TPRT mode: a short insertion's clip runs into the far flank; remap only the
                # inserted part (see _far_flank_trimmed)
                clips = [(_far_flank_trimmed(i, "L"), _far_flank_trimmed(i, "R")) for i in insertions]
            else:
                clips = [(i.left_clipped, i.right_clipped) for i in insertions]
            f.writelines(
                [(lc.fastq(f"{i.name}:L") if lc is not None else "")
                 + (rc_.fastq(f"{i.name}:R") if rc_ is not None else "")
                 for i, (lc, rc_) in zip(insertions, clips)])
        print(f"running bowtie2 {CONFIG['combine_insertions']['bowtie2_executable']} with index {CONFIG['combine_insertions']['bowtie2_index2']}")
        # -F 2308 = unmapped (4) + secondary (256) + supplementary (2048). Only the primary
        # alignment is ever read below, but -k 1000 emits up to 1000 records per clip, so
        # dropping the rest here keeps them out of the BAM instead of compressing, storing
        # and re-parsing records the loop discards. The equivalent Python skip is retained
        # below, so a clipped_bam written by an older version still yields the same calls.
        cmd = f"{CONFIG['combine_insertions']['bowtie2_executable']} {insertions_fasta} -k 1000 -x {CONFIG['combine_insertions']['bowtie2_index2']} --local --very-fast --threads {threads} --qc-filter | {CONFIG['combine_insertions']['samtools_executable']} view -F 2308 -b -o {clipped_bam}"
        os.system(cmd)
    filter_reads = set()
    print(f"removing all reads where one of the clipped ends maps within 1000bp of the breakpoint.")
    lo = pyliftover.LiftOver(CONFIG['combine_insertions']['bowtie2_index2_lo'])
    filter_reads = set()
    delta_sampler = []
    with pysam.AlignmentFile(clipped_bam) as f:
        for read in f:
            # Only a PRIMARY (best-scoring) clip alignment near the breakpoint indicates a
            # genuine reference-contiguous junction. A real MEI clip's best hit is a distant
            # element paralog; bowtie2 -k also emits weak SECONDARY multimapper hits, and one
            # of those can land near the breakpoint by chance and wrongly discard the real
            # insertion. Ignoring secondaries here recovers those calls (~33/10k on the
            # full-stack harness) while still catching true local misalignments (whose best
            # hit IS near the breakpoint).
            if read.is_secondary or read.is_supplementary:
                continue
            if read.is_mapped:
                reference_name, pos, side = read.query_name.split(":")
                if len(reference_name) < 3 or reference_name[:3] != 'chr':
                    reference_name = f"chr{reference_name}"
                left_pos, right_pos = pos.split("-")
                if side == 'R':
                    pos = int(right_pos)
                else:
                    pos = int(left_pos)
                map = lo.convert_coordinate(read.reference_name, read.reference_start if read.is_forward else read.reference_end)
                if map is not None:
                    for map_rn, map_coord, map_strand, map_len in map:
                        if reference_name == map_rn:
                            if abs(pos-map_coord) < 1000:
                                delta_sampler.append(abs(pos-map_coord))
                                filter_reads.add(read.query_name[:-2])
    #print(Counter(delta_sampler))
    print(f"detected {len(filter_reads)} insertions where the clipped part maps near the breakpoint. Removing these")
    insertions = [i for i in insertions if i.name not in filter_reads]
    if evidence is not None and evidence[3].get("pool") is not None:
        # TPRT (merge_tolerance_bp / far_pair_split): fold surviving one-sided loci into the
        # surviving call of the same junction (evidence pooled there), only now that both
        # passed every filter
        insertions, n_abs = evidence[3]["pool"].absorb_one_sided(insertions)
        if n_abs:
            print(f"folded {n_abs} one-sided loci into a surviving call of the same junction")
    with gzip.open(combined_insertions, 'wt') as f:
        f.writelines(
            [f'{i.left_consensus.fastq(f"{i.name}:L")}{i.right_consensus.fastq(f"{i.name}:R")}' for i in insertions if
             i.type is TYPE_FULL_INFO])
        f.writelines(
            [f'{i.left_consensus.fastq(f"{i.name}:L")}' for i in insertions if i.type is TYPE_RIGHT_POLYA])
        f.writelines(
            [f'{i.right_consensus.fastq(f"{i.name}:R")}' for i in insertions if i.type is TYPE_LEFT_POLYA])
        # Feature A: emit the one real side of each surviving discordant-anchored call.
        f.writelines(
            [f'{i.left_consensus.fastq(f"{i.name}:L")}' for i in insertions if i.type is TYPE_RIGHT_DISC])
        f.writelines(
            [f'{i.right_consensus.fastq(f"{i.name}:R")}' for i in insertions if i.type is TYPE_LEFT_DISC])
    print(f'wrote {len([i for i in insertions if i.name not in filter_reads])} insertions to insertions.txt.gz')
    # carry discovery's splice-hallmark (Feature B) evidence forward, re-keyed onto the
    # combined insertion names, for stage-4 processed-pseudogene annotation.
    write_combined_splice(input_files, insertions, combined_insertions)
    if evidence is not None:
        from combine_insertions_evidence import write_evidence_outputs
        default_tsv, default_fa = _evidence_paths(combined_insertions)
        # surviving insertions first (combined.txt.gz order), then the gated-out ones
        # (supported=0) for diagnostics.
        write_evidence_outputs(evidence_records,
                               [i.name for i in insertions] + sorted(evidence_failed),
                               evidence_tsv or default_tsv, reads_fa or default_fa)
    n_excluded = 0
    n_included = 0
    with gzip.open(insertions_genotyping_file, 'wt') as f:
        for i in insertions:
            if i.name in filter_reads:
                continue
            # Feature A: discordant-anchored calls carry only one real side and are not
            # genotyped in this prototype (they are still emitted to combined.txt.gz above).
            if i.type is TYPE_LEFT_DISC or i.type is TYPE_RIGHT_DISC:
                n_excluded += 1
                continue
            right = i.right_clipped.upper()[:CONFIG['genotyping']['max_bases']]
            left = i.left_clipped[:CONFIG['genotyping']['max_bases']].upper().revcomp()
            if len(i.reference_name) > 5:
                print(f"Excluding {i.name} since the breakpoint is on {i.reference_name}.")
                n_excluded += 1
                continue
            if i.type is not TYPE_RIGHT_POLYA:
                right_ref = get_sequence(i.reference_name, i.right_pos,
                                         i.right_pos + CONFIG['genotyping']['max_bases']).upper()
            else:
                right_ref = None
            if i.type is not TYPE_LEFT_POLYA:
                left_ref = get_sequence(i.reference_name, i.left_pos - CONFIG['genotyping']['max_bases'],
                                    i.left_pos).upper()
            else:
                left_ref = None
            if i.type is not TYPE_RIGHT_POLYA and not len(right_ref):
                print(f"Excluding {i.name} due to missing reference")
                n_excluded += 1
                continue
            if i.type is not TYPE_LEFT_POLYA and not len(left_ref):
                print(f"Excluding {i.name} due to missing reference")
                n_excluded += 1
                continue
            if right_ref is not None and 'N' in right_ref:
                print(f"Excluding {i.name} due to Ns in right reference")
                n_excluded += 1
                continue
            if left_ref is not None and 'N' in left_ref:
                print(f"Excluding {i.name} due to Ns in left reference")
                n_excluded += 1
                continue
            if i.type is TYPE_RIGHT_POLYA  is not None and sequence_matching_score([left_ref, str(left)]) > 0:
                print(f"Excluding {i.name} due to similar left sequence between ref and alt {left_ref} vs {left}")
                n_excluded += 1
                continue
            if i.type is TYPE_LEFT_POLYA and sequence_matching_score([right_ref, str(right)]) > 0:
                print(f"Excluding {i.name} due to similar right sequence between ref and alt: {right_ref} vs {right}")
                n_excluded += 1
                continue
            if i.type is TYPE_FULL_INFO and sequence_matching_score([right_ref, str(right)]) > 0 and sequence_matching_score([left_ref, str(left)]) > 0:
                print(f"Excluding {i.name} due to similar right and left sequence between ref and alt: {right_ref} vs {right} (right) and {left_ref} vs {left} (left)")
                n_excluded += 1
                continue
            if i.type is TYPE_FULL_INFO and (left_ref == revcomp(right_ref) or left_ref == right_ref):
                print(f"Excluding {i.name} due to identical clipped sequences: {left_ref} vs {right_ref}")
                n_excluded += 1
                continue
            n_included += 1
            f.write(f'>{i.name}\n')
            if i.type is not TYPE_RIGHT_POLYA:
                f.write(f'@RIGHT_INSERTION\n{right}\n')
                f.write(f'@RIGHT_REFERENCE\n{right_ref}\n')
            if i.type is not TYPE_LEFT_POLYA:
                f.write(f'@LEFT_INSERTION\n{left}\n')
                f.write(f'@LEFT_REFERENCE\n{left_ref}\n')
    print(f"wrote {n_included} of {n_included + n_excluded}, excluded {n_excluded}")
