//! COMPACT form of a `disc_agree_second_fragment` PENDING junction (`mem_budget_mb`;
//! analysis/discovery_mem/BUDGET.md).
//!
//! A pending junction is a fully built consensus `Breakpoint` (~1.75 KB with its sidecar
//! lite lists and `AgreePending`) that waits from the end of its contig's `cleanup` until
//! the mate pass has looked for the inside mate of one of its discordant anchors; ~85 % are
//! then rejected. Under a memory budget the pool holds instead, per junction, ONLY
//!
//! * what the mate pass needs to decide promotion: the raw clip inserts `agrees()` compares
//!   (4-bit packed, lossless for the 16 BAM base symbols), each candidate anchor's qname
//!   hash + read1/strand bits (exactly what `AgreeWait` carries in the FULL path, which also
//!   matches anchors by hash only), and the `has_mate` clip qnames that `get_mates` feeds
//!   into the mate-capture maps;
//! * the BGZF virtual position of every record needed to REBUILD the junction after the
//!   decision: every clip record of its cluster (scan order, duplicates included, because
//!   duplicate collapsing happens inside `join()`), then its selected DISC and SHORT sidecar
//!   reads in list order; the mate records captured for it in the mate pass are appended to
//!   `captured_mates` as they are met.
//!
//! Everything lives in one header vector plus one byte arena: no per-junction heap
//! allocation. A promoted junction is rebuilt by `Discovery::rebuild_compact` from exactly
//! these records through the same `join()` / lite constructors as the FULL path; a
//! rejected one is dropped without ever being read again. The form is never observable in
//! any output.

use crate::model::Breakpoint;

/// BAM 4-bit base alphabet (`=ACMGRSVTWYHKDBN`), the symbols noodles decodes a BAM
/// sequence to. A clip insert is `upper()`-cased read sequence, so it is always
/// representable; anything else (CRAM/odd input) falls back to raw bytes (`INS_RAW`).
const BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
const INS_RAW: u8 = 0;
const INS_PACKED: u8 = 1;

#[inline]
fn code_of(b: u8) -> Option<u8> {
    BASES.iter().position(|&x| x == b).map(|p| p as u8)
}

/// Pack a sequence to 4 bits/base. None when a byte is outside the BAM alphabet.
pub fn pack4(seq: &[u8]) -> Option<Vec<u8>> {
    let mut out = Vec::with_capacity(seq.len().div_ceil(2));
    for pair in seq.chunks(2) {
        let hi = code_of(pair[0])?;
        let lo = if pair.len() == 2 { code_of(pair[1])? } else { 0 };
        out.push((hi << 4) | lo);
    }
    Some(out)
}

/// Inverse of `pack4` for a sequence of `n` bases.
pub fn unpack4(packed: &[u8], n: usize, out: &mut Vec<u8>) {
    out.clear();
    out.reserve(n);
    for i in 0..n {
        let byte = packed[i / 2];
        let code = if i % 2 == 0 { byte >> 4 } else { byte & 0x0f };
        out.push(BASES[code as usize]);
    }
}

/// Fixed-size part of one compact junction; the variable part is `blob[off..off+len]`:
/// clip vps (n_clip x u64 LE), disc vps, short vps, anchors (n_anchor x (u64 frag, u8 bits:
/// 1 = read1, 2 = reverse)), inserts (n_ins x (u8 kind, u16 n_bases, bytes)), mates
/// (n_mate x (u8 read1, u16 len, qname bytes)).
#[derive(Clone, Copy, Debug)]
pub struct Hdr {
    off: u64,
    pub ord: u32,
    pub pos: i64,
    pub ref_id: i32,
    len: u32,
    n_clip: u32,
    n_disc: u32,
    n_short: u32,
    n_mate: u32,
    n_anchor: u16,
    n_ins: u8,
    pub left: bool,
    pub floor: u8,
}

/// The compact pending-junction pool (see the module doc).
#[derive(Default)]
pub struct CompactPool {
    hdr: Vec<Hdr>,
    blob: Vec<u8>,
    promoted: Vec<bool>,
    /// (junction index, virtual position of a mate record captured for it by the mate pass),
    /// in scan order — the compact counterpart of `Breakpoint::mate_seqs` / `mate_dests`.
    pub captured_mates: Vec<(u32, u64)>,
}

/// One candidate discordant anchor of a compact junction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Anchor {
    pub frag: u64,
    pub read1: bool,
    pub reverse: bool,
}

struct Cursor<'a> {
    b: &'a [u8],
    i: usize,
}

impl<'a> Cursor<'a> {
    fn u8(&mut self) -> u8 {
        let v = self.b[self.i];
        self.i += 1;
        v
    }
    fn u16(&mut self) -> u16 {
        let v = u16::from_le_bytes([self.b[self.i], self.b[self.i + 1]]);
        self.i += 2;
        v
    }
    fn u64(&mut self) -> u64 {
        let v = u64::from_le_bytes(self.b[self.i..self.i + 8].try_into().unwrap());
        self.i += 8;
        v
    }
    fn bytes(&mut self, n: usize) -> &'a [u8] {
        let s = &self.b[self.i..self.i + n];
        self.i += n;
        s
    }
}

impl CompactPool {
    pub fn len(&self) -> usize {
        self.hdr.len()
    }

    /// Accounted heap bytes of the pool (headers, arena, flags, captured mates; capacity).
    pub fn bytes(&self) -> usize {
        self.hdr.capacity() * std::mem::size_of::<Hdr>()
            + self.blob.capacity()
            + self.promoted.capacity()
            + self.captured_mates.capacity() * std::mem::size_of::<(u32, u64)>()
    }

    /// Arena bytes actually used (for the per-junction average in the log line).
    pub fn used_bytes(&self) -> usize {
        self.hdr.len() * std::mem::size_of::<Hdr>() + self.blob.len() + self.promoted.len()
    }

    /// Encode a FULL pending junction (`agree` and `src` present) and return its index.
    /// `left` = LEFT side. Panics if the breakpoint is not a budgeted pending junction (an
    /// internal invariant: `join()` records `src` for every pending junction when a budget
    /// is set).
    pub fn push_from(&mut self, bp: &Breakpoint, left: bool, floor: usize) -> u32 {
        let ag = bp.agree.as_deref().expect("compact: not a pending junction");
        let src = bp.src.as_deref().expect("compact: pending junction without record addresses");
        let off = self.blob.len() as u64;
        let b = &mut self.blob;
        for v in src.clip.iter().chain(&src.disc).chain(&src.short) {
            b.extend_from_slice(&v.to_le_bytes());
        }
        for d in &ag.anchors {
            b.extend_from_slice(&d.frag.to_le_bytes());
            b.push(u8::from(d.flag & 0x40 != 0) | (u8::from(d.is_reverse()) << 1));
        }
        for ins in &ag.inserts {
            assert!(ins.len() <= u16::MAX as usize, "compact: insert longer than u16");
            match pack4(ins) {
                Some(p) => {
                    b.push(INS_PACKED);
                    b.extend_from_slice(&(ins.len() as u16).to_le_bytes());
                    b.extend_from_slice(&p);
                }
                None => {
                    b.push(INS_RAW);
                    b.extend_from_slice(&(ins.len() as u16).to_le_bytes());
                    b.extend_from_slice(ins);
                }
            }
        }
        for (r1, q) in &bp.mates {
            assert!(q.len() <= u16::MAX as usize, "compact: qname longer than u16");
            b.push(u8::from(*r1));
            b.extend_from_slice(&(q.len() as u16).to_le_bytes());
            b.extend_from_slice(q.as_bytes());
        }
        let len = (self.blob.len() as u64 - off) as u32;
        self.hdr.push(Hdr {
            off,
            ord: bp.ord,
            pos: bp.breakpoint,
            ref_id: bp.ref_id,
            len,
            n_clip: src.clip.len() as u32,
            n_disc: src.disc.len() as u32,
            n_short: src.short.len() as u32,
            n_mate: bp.mates.len() as u32,
            n_anchor: ag.anchors.len() as u16,
            n_ins: ag.inserts.len() as u8,
            left,
            floor: floor as u8,
        });
        self.promoted.push(false);
        (self.hdr.len() - 1) as u32
    }

    pub fn hdr(&self, i: u32) -> &Hdr {
        &self.hdr[i as usize]
    }

    /// (ref_id, ord): the slot key shared with FULL final breakpoints.
    pub fn key(&self, i: u32) -> (i32, u32) {
        let h = &self.hdr[i as usize];
        (h.ref_id, h.ord)
    }

    pub fn promoted(&self, i: u32) -> bool {
        self.promoted[i as usize]
    }

    pub fn set_promoted(&mut self, i: u32) {
        self.promoted[i as usize] = true;
    }

    fn cur(&self, i: u32) -> Cursor<'_> {
        let h = &self.hdr[i as usize];
        Cursor { b: &self.blob[h.off as usize..(h.off + h.len as u64) as usize], i: 0 }
    }

    /// Record addresses: (clip records in scan order, DISC reads, SHORT reads).
    pub fn vps(&self, i: u32) -> (Vec<u64>, Vec<u64>, Vec<u64>) {
        let h = self.hdr[i as usize];
        let mut c = self.cur(i);
        let clip = (0..h.n_clip).map(|_| c.u64()).collect();
        let disc = (0..h.n_disc).map(|_| c.u64()).collect();
        let short = (0..h.n_short).map(|_| c.u64()).collect();
        (clip, disc, short)
    }

    fn skip_vps(&self, i: u32, c: &mut Cursor<'_>) {
        let h = &self.hdr[i as usize];
        c.i += 8 * (h.n_clip + h.n_disc + h.n_short) as usize;
    }

    pub fn anchors(&self, i: u32) -> Vec<Anchor> {
        let h = self.hdr[i as usize];
        let mut c = self.cur(i);
        self.skip_vps(i, &mut c);
        (0..h.n_anchor)
            .map(|_| {
                let frag = c.u64();
                let bits = c.u8();
                Anchor { frag, read1: bits & 1 != 0, reverse: bits & 2 != 0 }
            })
            .collect()
    }

    fn skip_anchors(&self, i: u32, c: &mut Cursor<'_>) {
        c.i += 9 * self.hdr[i as usize].n_anchor as usize;
    }

    /// Decode the clip inserts into `out` (reused scratch buffers; `out.len()` = n inserts).
    pub fn inserts_into(&self, i: u32, out: &mut Vec<Vec<u8>>) {
        let h = self.hdr[i as usize];
        let mut c = self.cur(i);
        self.skip_vps(i, &mut c);
        self.skip_anchors(i, &mut c);
        out.resize_with(h.n_ins as usize, Vec::new);
        for slot in out.iter_mut() {
            let kind = c.u8();
            let n = c.u16() as usize;
            if kind == INS_PACKED {
                let packed = c.bytes(n.div_ceil(2));
                unpack4(packed, n, slot);
            } else {
                slot.clear();
                slot.extend_from_slice(c.bytes(n));
            }
        }
    }

    /// (mate_is_read1, qname) of the junction's `has_mate` clip reads, in `join()` order.
    pub fn mates(&self, i: u32) -> Vec<(bool, &[u8])> {
        let h = self.hdr[i as usize];
        let mut c = self.cur(i);
        self.skip_vps(i, &mut c);
        self.skip_anchors(i, &mut c);
        for _ in 0..h.n_ins {
            let kind = c.u8();
            let n = c.u16() as usize;
            c.i += if kind == INS_PACKED { n.div_ceil(2) } else { n };
        }
        (0..h.n_mate)
            .map(|_| {
                let r1 = c.u8() != 0;
                let n = c.u16() as usize;
                (r1, c.bytes(n))
            })
            .collect()
    }

    /// Indices of one side's junctions in (ref_id, ord) order.
    pub fn side_order(&self, left: bool) -> Vec<u32> {
        let mut v: Vec<u32> = (0..self.hdr.len() as u32).filter(|&i| self.hdr[i as usize].left == left).collect();
        v.sort_by_key(|&i| self.key(i));
        v
    }

    /// Move another pool's junctions (and captured mates) into this one (SPD-3 merge).
    pub fn append(&mut self, mut other: CompactPool) {
        let base = self.hdr.len() as u32;
        let boff = self.blob.len() as u64;
        for h in &mut other.hdr {
            h.off += boff;
        }
        self.hdr.append(&mut other.hdr);
        self.blob.append(&mut other.blob);
        self.promoted.append(&mut other.promoted);
        self.captured_mates.extend(other.captured_mates.drain(..).map(|(i, vp)| (i + base, vp)));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::DiscLite;
    use crate::model::{AgreePending, SrcOffsets};
    use crate::qseq::QualitySeq;

    #[test]
    fn pack4_roundtrips_the_bam_alphabet_and_rejects_others() {
        for s in [&b"ACGTN"[..], b"=ACMGRSVTWYHKDBN", b"A", b"", b"NNACGTAC"] {
            let p = pack4(s).unwrap();
            assert_eq!(p.len(), s.len().div_ceil(2));
            let mut out = Vec::new();
            unpack4(&p, s.len(), &mut out);
            assert_eq!(out, s);
        }
        assert!(pack4(b"ACGTa").is_none(), "lower case is not in the alphabet");
        assert!(pack4(b"AC.G").is_none());
    }

    fn pending(ord: u32) -> Breakpoint {
        let mut b = Breakpoint::new(
            crate::config::CLIP_LEFT,
            "1".into(),
            1000 + ord as i64,
            None,
            QualitySeq::new(b"ACGT".to_vec(), vec![30; 4]),
            QualitySeq::new(b"ACGT".to_vec(), vec![30; 4]),
            None,
            None,
            false,
            0,
        );
        b.ref_id = 0;
        b.ord = ord;
        b.mates = vec![(true, format!("q{ord}/x")), (false, "other".into())];
        let d = |frag, flag| DiscLite { frag, flag, ref_id: 0, start: 1, end: 2, mref: -1, mpos: -1 };
        b.agree = Some(Box::new(AgreePending {
            inserts: vec![b"ACGTTTGCANNACGATCGATCGTAGCTAGCTAGC".to_vec(), b"acgtX".to_vec()],
            clips: Vec::new(),
            anchors: vec![d(11, 0x1 | 0x40 | 0x10), d(12, 0x1 | 0x80), d(13, 0x10)],
            promoted: false,
        }));
        b.src = Some(Box::new(SrcOffsets { clip: vec![100 + ord as u64, 200], disc: vec![300, 301, 302], short: vec![400] }));
        b
    }

    #[test]
    fn roundtrip_every_field_and_no_per_junction_allocation() {
        let mut p = CompactPool::default();
        let i = p.push_from(&pending(5), true, 2);
        let j = p.push_from(&pending(6), false, 2);
        assert_eq!((i, j, p.len()), (0, 1, 2));
        assert_eq!(p.vps(0), (vec![105, 200], vec![300, 301, 302], vec![400]));
        assert_eq!(p.vps(1).0, vec![106, 200]);
        assert_eq!(
            p.anchors(0),
            vec![
                Anchor { frag: 11, read1: true, reverse: true },
                Anchor { frag: 12, read1: false, reverse: false },
                Anchor { frag: 13, read1: false, reverse: true }
            ]
        );
        let mut ins = Vec::new();
        p.inserts_into(1, &mut ins);
        assert_eq!(ins, vec![b"ACGTTTGCANNACGATCGATCGTAGCTAGCTAGC".to_vec(), b"acgtX".to_vec()], "packed and raw fallback");
        let m = p.mates(0);
        assert_eq!(m, vec![(true, &b"q5/x"[..]), (false, &b"other"[..])]);
        assert_eq!(p.key(1), (0, 6));
        assert!(p.hdr(0).left && !p.hdr(1).left);
        assert_eq!((p.hdr(0).pos, p.hdr(0).floor), (1005, 2));
        assert!(!p.promoted(0));
        p.set_promoted(0);
        assert!(p.promoted(0) && !p.promoted(1));
        // the arena is one Vec: a second pool appended keeps every field
        let mut q = CompactPool::default();
        q.push_from(&pending(7), true, 2);
        q.captured_mates.push((0, 999));
        p.captured_mates.push((1, 888));
        p.append(q);
        assert_eq!(p.len(), 3);
        assert_eq!(p.vps(2).0, vec![107, 200]);
        assert_eq!(p.mates(2)[0].1, b"q7/x");
        assert_eq!(p.captured_mates, vec![(1, 888), (2, 999)]);
        assert_eq!(p.side_order(true), vec![0, 2]);
        assert_eq!(p.side_order(false), vec![1]);
        // ~ header + arena only
        assert!(p.used_bytes() < 3 * 400, "used {}", p.used_bytes());
    }
}
