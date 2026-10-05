//! PEAR-TREE combine_insertions (Rust port). See SPEC.md (behaviour, formats) and PLAN.md
//! (work packages, memory/parallelism design).

// skeleton phase: todo!() bodies leave parameters unused. Remove at integration (PLAN.md WP7).
#![allow(unused_variables, unused_imports, dead_code)]

pub mod align;
pub mod config;
pub mod consensus;
pub mod context;
pub mod evidence;
pub mod genome;
pub mod genotyping_out;
pub mod insertion;
pub mod intersect;
pub mod library;
pub mod liftover;
pub mod model;
pub mod pipeline;
pub mod pyfmt;
pub mod region_filter;
pub mod remap;
pub mod seq;
pub mod splice;
pub mod tprt;
