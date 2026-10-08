//! Compact read sequence storage. FOUNDATION (implemented).
//!
//! A read is stored 4 bits per base in the BAM nibble alphabet `=ACMGRSVTWYHKDBN`, so any
//! upper-case IUPAC sequence round-trips losslessly at half the bytes; anything else (a byte
//! outside that alphabet, e.g. '*' or '.') falls back to raw bytes. Case is NOT kept: every
//! consumer of `EvidenceRead.seq` in tools/rte upper-cases first (assembly.layout,
//! sequtil.edlib_best in pseudogene.find), so `decode()` returns upper case.

const NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

fn code(b: u8) -> Option<u8> {
    let u = b.to_ascii_uppercase();
    NT16.iter().position(|&c| c == u).map(|i| i as u8)
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub enum PackedSeq {
    /// (number of bases, 2 bases per byte, high nibble first)
    Nibble(u32, Box<[u8]>),
    /// raw bytes as read (upper-cased)
    Raw(Box<[u8]>),
}

impl PackedSeq {
    pub fn encode(seq: &[u8]) -> PackedSeq {
        if seq.len() <= u32::MAX as usize && seq.iter().all(|&b| code(b).is_some()) {
            let mut out = vec![0u8; seq.len().div_ceil(2)];
            for (i, &b) in seq.iter().enumerate() {
                let c = code(b).unwrap();
                out[i / 2] |= if i % 2 == 0 { c << 4 } else { c };
            }
            PackedSeq::Nibble(seq.len() as u32, out.into_boxed_slice())
        } else {
            PackedSeq::Raw(seq.to_ascii_uppercase().into_boxed_slice())
        }
    }

    pub fn len(&self) -> usize {
        match self {
            PackedSeq::Nibble(n, _) => *n as usize,
            PackedSeq::Raw(b) => b.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The sequence, upper case.
    pub fn decode(&self) -> Vec<u8> {
        match self {
            PackedSeq::Nibble(n, b) => (0..*n as usize)
                .map(|i| {
                    let byte = b[i / 2];
                    NT16[(if i % 2 == 0 { byte >> 4 } else { byte & 15 }) as usize]
                })
                .collect(),
            PackedSeq::Raw(b) => b.to_vec(),
        }
    }

    /// Serialized form used by the spill file: tag byte, u32 length, payload.
    pub(crate) fn write_to(&self, out: &mut Vec<u8>) {
        match self {
            PackedSeq::Nibble(n, b) => {
                out.push(0);
                out.extend_from_slice(&n.to_le_bytes());
                out.extend_from_slice(b);
            }
            PackedSeq::Raw(b) => {
                out.push(1);
                out.extend_from_slice(&(b.len() as u32).to_le_bytes());
                out.extend_from_slice(b);
            }
        }
    }

    /// Inverse of `write_to`; returns the value and the bytes consumed.
    pub(crate) fn read_from(buf: &[u8]) -> Option<(PackedSeq, usize)> {
        let tag = *buf.first()?;
        let n = u32::from_le_bytes(buf.get(1..5)?.try_into().ok()?) as usize;
        match tag {
            0 => {
                let nb = n.div_ceil(2);
                let b = buf.get(5..5 + nb)?;
                Some((PackedSeq::Nibble(n as u32, b.into()), 5 + nb))
            }
            1 => {
                let b = buf.get(5..5 + n)?;
                Some((PackedSeq::Raw(b.into()), 5 + n))
            }
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip() {
        for s in [&b""[..], b"A", b"ACGTN", b"acgtnRYKM", b"ACGT*", b"NNNNNNN"] {
            let p = PackedSeq::encode(s);
            assert_eq!(p.decode(), s.to_ascii_uppercase(), "{:?}", String::from_utf8_lossy(s));
            assert_eq!(p.len(), s.len());
            let mut buf = Vec::new();
            p.write_to(&mut buf);
            let (q, used) = PackedSeq::read_from(&buf).unwrap();
            assert_eq!(used, buf.len());
            assert_eq!(q, p);
        }
        assert!(matches!(PackedSeq::encode(b"ACGT"), PackedSeq::Nibble(4, _)));
        assert!(matches!(PackedSeq::encode(b"AC*T"), PackedSeq::Raw(_)));
    }
}
