//! Port of src/discovery.py (Discovery).

use std::fs::File;
use std::io::{self, Read, Write};
use std::num::NonZeroUsize;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::config::*;
use crate::filters::{clean_clipped_seq, is_adapter};
use crate::intervals::IntervalIndex;
use crate::model::{join, Breakpoint};
use crate::polya::PolyABreakpoint;
use crate::qseq::{revcomp_bytes, QualitySeq};
use crate::read::BamRead;
use crate::stats::Stats;

#[derive(Clone, Copy)]
enum BpRef {
    Left(usize),
    Right(usize),
    PolyA(usize),
}

enum Emit<'a> {
    Bp(&'a Breakpoint),
    Pa(&'a PolyABreakpoint),
}

/// Open the BAM as a record reader. With `bam_threads > 1` the BGZF blocks are
/// decoded on a worker pool (SPD-2); otherwise the single-threaded decoder is used,
/// preserving the documented single-core footprint. Both yield identical records in
/// identical order, so output stays byte-identical.
fn open_bam(path: &str, bam_threads: usize) -> io::Result<bam::io::Reader<Box<dyn Read>>> {
    let file = File::open(path)?;
    let inner: Box<dyn Read> = if bam_threads > 1 {
        let wc = NonZeroUsize::new(bam_threads).unwrap_or(NonZeroUsize::MIN);
        Box::new(bgzf::io::MultithreadedReader::with_worker_count(wc, file))
    } else {
        Box::new(bgzf::io::Reader::new(file))
    };
    Ok(bam::io::Reader::from(inner))
}

/// True if an alternative-alignment CIGAR (XA or SA tag) covers essentially the
/// whole read end-to-end (clips fewer than MIN_CLIP_LEN bases). Mirrors the
/// Python `_alt_is_full_length`.
fn alt_is_full_length(cigar: &str) -> bool {
    let mut aligned: usize = 0;
    let mut clip: usize = 0;
    let mut num: usize = 0;
    let mut saw_digit = false;
    for ch in cigar.bytes() {
        if ch.is_ascii_digit() {
            num = num * 10 + (ch - b'0') as usize;
            saw_digit = true;
        } else {
            if saw_digit {
                match ch {
                    b'S' | b'H' => clip += num,
                    b'M' | b'=' | b'X' => aligned += num,
                    _ => {}
                }
            }
            num = 0;
            saw_digit = false;
        }
    }
    aligned > 0 && clip < MIN_CLIP_LEN
}

/// True if the read has any XA/SA alternative alignment spanning the whole read.
/// Such reads map contiguously elsewhere and are not genuine junctions.
fn maps_fully_elsewhere(read: &BamRead<'_>) -> bool {
    // (tag string, index of the CIGAR field within a comma-separated entry)
    for (tag, cigar_idx) in [(read.xa(), 2usize), (read.sa(), 3usize)] {
        if let Some(s) = tag {
            for entry in s.split(';') {
                if entry.is_empty() {
                    continue;
                }
                let fields: Vec<&str> = entry.split(',').collect();
                if fields.len() > cigar_idx && alt_is_full_length(fields[cigar_idx]) {
                    return true;
                }
            }
        }
    }
    false
}

pub struct Discovery {
    temporary_breakpoints: Vec<Breakpoint>,
    final_left_breakpoints: Vec<Breakpoint>,
    final_right_breakpoints: Vec<Breakpoint>,
    polya: Vec<PolyABreakpoint>,
    reference_name: Option<String>,
    filepath: String,
    config: DiscoveryConfig,
    exclude: Option<IntervalIndex>,
    stats: Stats,
    bam_threads: usize, // BGZF decode workers (SPD-2); 1 = single-threaded
}

impl Discovery {
    pub fn new(filepath: String, bam_threads: usize, config: DiscoveryConfig, exclude: Option<IntervalIndex>) -> Self {
        Discovery {
            temporary_breakpoints: Vec::new(),
            final_left_breakpoints: Vec::new(),
            final_right_breakpoints: Vec::new(),
            polya: Vec::new(),
            reference_name: None,
            filepath,
            config,
            exclude,
            stats: Stats::default(),
            bam_threads: bam_threads.max(1),
        }
    }

    /// OBS-1 reject-counter sidecar, serialised as JSON.
    pub fn stats_json(&self) -> String {
        self.stats.to_json()
    }

    /// Contig filter for the read-scan (extract_chimeric). With an allowlist set the
    /// name must be in it; otherwise the legacy `len(name) <= 5` + not-MT heuristic.
    fn contig_ok_extract(&self, name: &str) -> bool {
        match &self.config.contig_allowlist {
            Some(set) => set.contains(name),
            None => name.len() <= 5 && name != "MT" && name != "chrM",
        }
    }

    /// Contig filter at output. Legacy path checks only `len(name) <= 5` (MT/chrM
    /// clipped breakpoints were already dropped in the scan, but MT polyA mates were
    /// historically emitted), so keep that exact behaviour when no allowlist is set.
    fn contig_ok_output(&self, name: &str) -> bool {
        match &self.config.contig_allowlist {
            Some(set) => set.contains(name),
            None => name.len() <= 5,
        }
    }

    /// True if a breakpoint at (name, pos) falls in an exclude-BED region (SPEC-5).
    fn excluded(&self, name: &str, pos: i64) -> bool {
        self.exclude.as_ref().is_some_and(|ix| ix.contains(name, pos))
    }

    #[allow(clippy::too_many_arguments)]
    fn add_breakpoint(
        &mut self,
        side: i32,
        breakpoint: i64,
        qname: &str,
        clipped: QualitySeq,
        unclipped: QualitySeq,
        is_read1: bool,
        is_forward: bool,
        exclude: bool,
    ) {
        let bp = Breakpoint::new(
            side,
            self.reference_name.clone().unwrap_or_default(),
            breakpoint,
            Some(qname.to_string()),
            clipped.upper(),
            unclipped.upper(),
            Some(is_read1),
            Some(is_forward),
            exclude,
        );
        self.temporary_breakpoints.push(bp);
    }

    fn cleanup(&mut self) {
        if self.reference_name.is_none() {
            return;
        }
        if self.temporary_breakpoints.is_empty() {
            return;
        }
        let mut left_bps: Vec<Breakpoint> = Vec::new();
        let mut right_bps: Vec<Breakpoint> = Vec::new();
        for bp in self.temporary_breakpoints.drain(..) {
            if bp.side == CLIP_LEFT {
                left_bps.push(bp);
            } else {
                right_bps.push(bp);
            }
        }
        // stable sort by breakpoint position
        left_bps.sort_by_key(|b| b.breakpoint);
        right_bps.sort_by_key(|b| b.breakpoint);

        for (side_bps, out) in [
            (left_bps, true),  // true -> final_left
            (right_bps, false),
        ] {
            let mut groups: Vec<Vec<Breakpoint>> = Vec::new();
            let mut current: Vec<Breakpoint> = Vec::new();
            for bp in side_bps {
                if let Some(last) = current.last() {
                    if (bp.breakpoint - last.breakpoint).abs() < self.config.cluster_window {
                        current.push(bp);
                        continue;
                    } else {
                        groups.push(std::mem::take(&mut current));
                    }
                }
                current = vec![bp];
            }
            if !current.is_empty() {
                groups.push(current);
            }
            for g in groups {
                if let Some(joined) = join(g, &self.config, &mut self.stats) {
                    if out {
                        self.final_left_breakpoints.push(joined);
                    } else {
                        self.final_right_breakpoints.push(joined);
                    }
                }
            }
        }
    }

    pub fn extract_chimeric(&mut self) -> io::Result<()> {
        let mut reader = open_bam(&self.filepath, self.bam_threads)?;
        let header = reader.read_header()?;
        let reject_fullmap = self.config.reject_fully_mapping_reads;
        let min_mapq = self.config.min_mapq;
        self.reference_name = None;
        let mut current_ref_id: Option<usize> = None;
        let mut record = bam::Record::default();
        while reader.read_record(&mut record)? != 0 {
            let read = BamRead::from_record(&record, &header)?;

            if read.mapq < min_mapq {
                if let Some(b) = PolyABreakpoint::find_polya(&read) {
                    self.polya.push(b);
                }
                continue;
            }
            if read.is_secondary || read.is_qcfail || read.is_duplicate {
                continue;
            }
            // SPD-5a: detect contig change by integer reference id; resolve the name
            // string only when it changes, not for every read.
            let ref_id = match read.reference_sequence_id {
                Some(id) => id,
                None => continue, // high-mapq unmapped: skip (would crash pysam)
            };
            if current_ref_id != Some(ref_id) {
                self.cleanup();
                current_ref_id = Some(ref_id);
                self.reference_name = read.reference_name();
            }
            let ref_name = self.reference_name.as_deref().unwrap_or_default();
            if !self.contig_ok_extract(ref_name) {
                continue;
            }
            if !read.has_cigar {
                continue;
            }

            let left_soft = read.left_is_soft;
            let right_soft = read.right_is_soft;
            let left_len = read.left_len;
            let right_len = read.right_len;

            let clip: Option<i32> = if left_soft && !right_soft {
                Some(CLIP_LEFT)
            } else if right_soft && !left_soft {
                Some(CLIP_RIGHT)
            } else if right_soft && left_soft {
                if left_len > right_len {
                    Some(CLIP_LEFT)
                } else if right_len > left_len {
                    Some(CLIP_RIGHT)
                } else {
                    None
                }
            } else {
                None
            };

            let Some(clip) = clip else { continue };

            // drop reads that map contiguously elsewhere in the reference
            // (XA/SA full-length alt): not real chimeric junctions.
            if reject_fullmap && maps_fully_elsewhere(&read) {
                continue;
            }

            // cruciform / short-indel exclusion via SA tag
            let mut exclude_flag = false;
            if read.is_supplementary {
                if let Some(sa) = read.sa() {
                    for part in sa.split(';') {
                        if part.is_empty() {
                            continue;
                        }
                        let fields: Vec<&str> = part.split(',').collect();
                        if fields.len() >= 2 && fields[0] == ref_name {
                            if let Ok(start) = fields[1].parse::<i64>() {
                                if (read.reference_start - start).abs() < self.config.exclude_same_contig_supplementary {
                                    exclude_flag = true;
                                }
                            }
                        }
                    }
                }
            }

            // surviving clip candidate: now decode the heavy fields (SPD-1).
            let seq = read.seq();
            let qname = read.query_name();
            let full = QualitySeq::new(seq.clone(), read.qual());
            let n = seq.len();

            if clip == CLIP_LEFT {
                // adapter check on the last MIN_CLIP_LEN of the clipped part, revcomped
                let clip_part = &seq[..left_len];
                let sub = &clip_part[clip_part.len().saturating_sub(MIN_CLIP_LEN)..];
                if is_adapter(&revcomp_bytes(sub)) {
                    continue;
                }
                let clipped = full.pyslice(None, Some(left_len as isize));
                let unclipped = if right_soft {
                    full.pyslice(Some(left_len as isize), Some(-(right_len as isize)))
                } else {
                    full.pyslice(Some(left_len as isize), None)
                };
                self.add_breakpoint(clip, read.reference_start, &qname, clipped, unclipped, read.is_read1, read.is_forward(), exclude_flag);
            } else {
                // CLIP_RIGHT: adapter check on the first MIN_CLIP_LEN of the clipped part
                let clip_part = &seq[n - right_len..];
                let sub = &clip_part[..clip_part.len().min(MIN_CLIP_LEN)];
                if is_adapter(sub) {
                    continue;
                }
                let clipped = full.pyslice(Some(-(right_len as isize)), None);
                let unclipped = if left_soft {
                    full.pyslice(Some(left_len as isize), Some(-(right_len as isize)))
                } else {
                    full.pyslice(None, Some(-(right_len as isize)))
                };
                self.add_breakpoint(clip, read.reference_end, &qname, clipped, unclipped, read.is_read1, read.is_forward(), exclude_flag);
            }
        }
        self.cleanup();
        Ok(())
    }

    fn get_mates(&self) -> (FxHashSet<String>, FxHashSet<String>, FxHashMap<String, BpRef>) {
        let mut read1_bp: FxHashSet<String> = FxHashSet::default();
        let mut read2_bp: FxHashSet<String> = FxHashSet::default();
        let mut qmap: FxHashMap<String, BpRef> = FxHashMap::default();

        for (idx, bp) in self.final_left_breakpoints.iter().enumerate() {
            for (read_1, qname) in &bp.mates {
                if *read_1 {
                    read1_bp.insert(qname.clone());
                } else {
                    read2_bp.insert(qname.clone());
                }
                qmap.insert(qname.clone(), BpRef::Left(idx));
            }
        }
        for (idx, bp) in self.final_right_breakpoints.iter().enumerate() {
            for (read_1, qname) in &bp.mates {
                if *read_1 {
                    read1_bp.insert(qname.clone());
                } else {
                    read2_bp.insert(qname.clone());
                }
                qmap.insert(qname.clone(), BpRef::Right(idx));
            }
        }
        for (idx, pa) in self.polya.iter().enumerate() {
            if qmap.contains_key(&pa.qname) {
                continue; // real breakpoints have priority
            }
            qmap.insert(pa.qname.clone(), BpRef::PolyA(idx));
            if pa.is_read1 {
                read2_bp.insert(pa.qname.clone());
            } else {
                read1_bp.insert(pa.qname.clone());
            }
        }
        let union: FxHashSet<String> = read1_bp.intersection(&read2_bp).cloned().collect();
        for q in &union {
            read1_bp.remove(q);
            read2_bp.remove(q);
        }
        (read1_bp, read2_bp, qmap)
    }

    pub fn find_mates(&mut self) -> io::Result<()> {
        let (read1_mates, read2_mates, qmap) = self.get_mates();
        let min_mapq = self.config.min_mapq;
        let mut reader = open_bam(&self.filepath, self.bam_threads)?;
        let _header = reader.read_header()?;
        let mut record = bam::Record::default();
        while reader.read_record(&mut record)? != 0 {
            let read = BamRead::from_record(&record, &_header)?;
            if read.is_secondary || read.is_qcfail || read.is_duplicate {
                continue;
            }
            // query_name is needed for every read (membership check), so decode it here.
            let qname = read.query_name();
            if read.is_read1 {
                if !read1_mates.contains(&qname) {
                    continue;
                }
            } else if !read2_mates.contains(&qname) {
                continue;
            }
            let Some(&bpref) = qmap.get(&qname) else { continue };
            match bpref {
                BpRef::PolyA(i) => {
                    self.polya[i].set_mate(&read, min_mapq);
                }
                BpRef::Left(i) | BpRef::Right(i) => {
                    let seq = clean_clipped_seq(&QualitySeq::new(read.seq(), read.qual()));
                    let seq = if read.is_forward() { seq } else { seq.revcomp() };
                    let bp = match bpref {
                        BpRef::Left(_) => &mut self.final_left_breakpoints[i],
                        BpRef::Right(_) => &mut self.final_right_breakpoints[i],
                        _ => unreachable!(),
                    };
                    bp.mate_seqs.push(seq);
                }
            }
        }
        Ok(())
    }

    pub fn discovery(&mut self) -> io::Result<()> {
        self.extract_chimeric()?;
        self.find_mates()?;
        // extend_mates() is a no-op in the Python (operates on the already-emptied
        // temporary_breakpoints); intentionally omitted.
        Ok(())
    }

    pub fn output<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        // group by reference, preserving the Python dict-key ordering of polyA
        let mut polya_order: Vec<String> = Vec::new();
        let mut polya_map: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        for (i, pa) in self.polya.iter().enumerate() {
            let (rn, bp) = (pa.reference_name.as_ref(), pa.breakpoint);
            let (Some(rn), Some(_)) = (rn, bp) else { continue };
            if !polya_map.contains_key(rn) {
                polya_order.push(rn.clone());
                polya_map.insert(rn.clone(), Vec::new());
            }
            polya_map.get_mut(rn).unwrap().push(i);
        }
        let mut left_map: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        for (i, bp) in self.final_left_breakpoints.iter().enumerate() {
            left_map.entry(bp.reference_name.clone()).or_default().push(i);
        }
        let mut right_map: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        for (i, bp) in self.final_right_breakpoints.iter().enumerate() {
            right_map.entry(bp.reference_name.clone()).or_default().push(i);
        }

        // union of all reference names, sorted (lexicographic, matching Python sorted())
        let mut union: Vec<String> = polya_map
            .keys()
            .chain(left_map.keys())
            .chain(right_map.keys())
            .cloned()
            .collect();
        union.sort();
        union.dedup();
        for rn in union {
            if !self.contig_ok_output(&rn) {
                continue;
            }
            if !polya_map.contains_key(&rn) {
                polya_order.push(rn.clone());
                polya_map.insert(rn, Vec::new());
            }
        }

        let empty: Vec<usize> = Vec::new();
        for rn in &polya_order {
            if !self.contig_ok_output(rn) {
                continue;
            }
            let mut l: Vec<&Breakpoint> = left_map.get(rn).unwrap_or(&empty).iter().map(|&i| &self.final_left_breakpoints[i]).collect();
            let mut r: Vec<&Breakpoint> = right_map.get(rn).unwrap_or(&empty).iter().map(|&i| &self.final_right_breakpoints[i]).collect();
            let mut p: Vec<&PolyABreakpoint> = polya_map.get(rn).unwrap_or(&empty).iter().map(|&i| &self.polya[i]).collect();
            l.sort_by_key(|b| b.breakpoint);
            r.sort_by_key(|b| b.breakpoint);
            p.sort_by_key(|b| b.breakpoint.unwrap());

            // SPEC-5: drop breakpoints inside exclude-BED regions (no-op if unset)
            if self.exclude.is_some() {
                l.retain(|b| !self.excluded(rn, b.breakpoint));
                r.retain(|b| !self.excluded(rn, b.breakpoint));
                p.retain(|b| !self.excluded(rn, b.breakpoint.unwrap()));
            }

            let (mut il, mut ir, mut ip) = (0usize, 0usize, 0usize);
            while il < l.len() && ir < r.len() {
                let tsd = r[ir].breakpoint - l[il].breakpoint;
                if tsd < self.config.tsd_min {
                    while ip < p.len() && p[ip].breakpoint.unwrap() - l[il].breakpoint < self.config.polya_near_dist {
                        ip += 1;
                    }
                    if ip != p.len() && p[ip].breakpoint.unwrap() - l[il].breakpoint < self.config.polya_far_dist && p[ip].clip == CLIP_RIGHT {
                        print_output(writer, Emit::Bp(l[il]), Emit::Pa(p[ip]))?;
                        il += 1;
                        continue;
                    }
                    ir += 1;
                    continue;
                } else if tsd > self.config.tsd_max {
                    while ip < p.len() && r[ir].breakpoint - p[ip].breakpoint.unwrap() > self.config.polya_far_dist {
                        ip += 1;
                    }
                    if ip != p.len() && r[ir].breakpoint - p[ip].breakpoint.unwrap() > self.config.polya_near_dist && p[ip].clip == CLIP_LEFT {
                        print_output(writer, Emit::Pa(p[ip]), Emit::Bp(r[ir]))?;
                        ir += 1;
                        continue;
                    }
                    il += 1;
                    continue;
                } else {
                    print_output(writer, Emit::Bp(l[il]), Emit::Bp(r[ir]))?;
                    il += 1;
                }
            }
        }
        Ok(())
    }
}

fn print_output<W: Write>(writer: &mut W, left: Emit, right: Emit) -> io::Result<()> {
    let (reference_name, left_str, right_str) = match (&left, &right) {
        (Emit::Pa(pa), Emit::Bp(rb)) => (
            pa.reference_name.clone().unwrap(),
            format!("polyA_{}", pa.breakpoint.unwrap()),
            format!("{}", rb.breakpoint),
        ),
        (Emit::Bp(lb), Emit::Pa(pa)) => (
            pa.reference_name.clone().unwrap(),
            format!("{}", lb.breakpoint),
            format!("polyA_{}", pa.breakpoint.unwrap()),
        ),
        (Emit::Bp(lb), Emit::Bp(rb)) => (
            lb.reference_name.clone(),
            format!("{}", lb.breakpoint),
            format!("{}", rb.breakpoint),
        ),
        (Emit::Pa(_), Emit::Pa(_)) => unreachable!("two polyA breakpoints are never paired"),
    };
    let bp_name = format!("{}:{}-{}", reference_name, left_str, right_str);

    match &left {
        Emit::Pa(pa) => {
            let c = pa.clipped.as_ref().unwrap().revcomp();
            writer.write_all(c.fastq(&format!("{}:LEFT:CLIPPED_POLYA", bp_name)).as_bytes())?;
        }
        Emit::Bp(lb) => {
            writer.write_all(lb.clipped.revcomp().fastq(&format!("{}:LEFT:CLIPPED", bp_name)).as_bytes())?;
            writer.write_all(lb.unclipped.fastq(&format!("{}:LEFT:ALIGNED", bp_name)).as_bytes())?;
            for (i, mate) in lb.mate_seqs.iter().enumerate() {
                writer.write_all(mate.revcomp().fastq(&format!("{}:LEFT:MATE{}", bp_name, i)).as_bytes())?;
            }
        }
    }
    match &right {
        Emit::Pa(pa) => {
            let c = pa.clipped.as_ref().unwrap();
            writer.write_all(c.fastq(&format!("{}:RIGHT:CLIPPED_POLYA", bp_name)).as_bytes())?;
        }
        Emit::Bp(rb) => {
            writer.write_all(rb.unclipped.revcomp().fastq(&format!("{}:RIGHT:ALIGNED", bp_name)).as_bytes())?;
            writer.write_all(rb.clipped.fastq(&format!("{}:RIGHT:CLIPPED", bp_name)).as_bytes())?;
            for (i, mate) in rb.mate_seqs.iter().enumerate() {
                writer.write_all(mate.revcomp().fastq(&format!("{}:RIGHT:MATE{}", bp_name, i)).as_bytes())?;
            }
        }
    }
    Ok(())
}
