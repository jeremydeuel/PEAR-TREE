//! peartree-rte: Rust port of tools/rte (TPRT-hallmark annotation of PEAR-TREE insertions).
//!
//! Module map (python source -> status) is in PORT_PLAN.md; the annotate_v2 contract in SPEC.md.
//! FOUNDATION modules (implemented): align, pyfmt, packed, sequtil, mm, genome, config, inputs,
//! stream, record, library, io, golden. Stubs (`todo!`) for the work packages: assembly,
//! structure, transduction (partly), hallmarks, score, pseudogene, genemodel, annotator.

pub mod align;
pub mod annotator;
pub mod assembly;
pub mod config;
pub mod genemodel;
pub mod genome;
pub mod golden;
pub mod hallmarks;
pub mod inputs;
pub mod io;
pub mod library;
pub mod mm;
pub mod packed;
pub mod pseudogene;
pub mod pyfmt;
pub mod record;
pub mod score;
pub mod sequtil;
pub mod stream;
pub mod structure;
pub mod transduction;
