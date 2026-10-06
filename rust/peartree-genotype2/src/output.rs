//! Output row formatting (SPEC "Output"). Owner: C.
//!
//! Column order is `types::OUTPUT_HEADER`; every row has all 15 columns:
//! `insertion genotype score_genotype score_alternative coverage n_alt n_ref n_art vaf gq
//! pl_ref pl_het pl_hom n_uninf n_disc`.

use crate::types::{Call, GT_HIGH_COVERAGE};

/// Format one output row (ends with '\n'). `coverage` and `n_disc` come from the driver.
// live once driver.rs (owner D) calls it; remove at integration
#[allow(dead_code)]
pub fn format_row(name: &str, call: &Call, coverage: i64, n_disc: i64) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        name,
        call.genotype,
        call.score_genotype,
        call.score_alternative,
        coverage,
        call.n_alt,
        call.n_ref,
        call.n_art,
        call.vaf,
        call.gq,
        call.pl[0],
        call.pl[1],
        call.pl[2],
        call.n_uninf,
        n_disc
    )
}

/// Rows for the driver-decided states: `high-coverage` (coverage known), `error`, `no-coverage`.
/// As in the legacy genotyper, a `high-coverage` row carries the (capped) coverage as its
/// `score_genotype`; every other count / score is 0, vaf `0.000`, gq 0, PL `0 0 0`.
// live once driver.rs (owner D) calls it; remove at integration
#[allow(dead_code)]
pub fn format_simple_row(name: &str, genotype: &'static str, coverage: i64) -> String {
    let score_genotype = if genotype == GT_HIGH_COVERAGE { coverage } else { 0 };
    format!("{name}\t{genotype}\t{score_genotype}\t0\t{coverage}\t0\t0\t0\t0.000\t0\t0\t0\t0\t0\t0\n")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Config;
    use crate::model::call_locus;
    use crate::types::{ReadClass, ReadObs, GT_ERROR, GT_NO_COVERAGE, OUTPUT_HEADER};

    fn obs(class: ReadClass, llr: f64) -> ReadObs {
        ReadObs { ll_ref: -40.0, ll_alt: -40.0 + llr, class, explained_frac: 1.0, crosses_junction: true }
    }

    #[test]
    fn hand_computed_row() {
        // see model::tests::scores_and_hand_computed_case for the numbers
        let reads = [
            obs(ReadClass::Alt, 10.0),
            obs(ReadClass::Alt, 15.0),
            obs(ReadClass::Ref, -20.0),
            obs(ReadClass::Uninformative, 1.0),
            obs(ReadClass::Unexplained, 0.0),
        ];
        let c = call_locus(&reads, 2, &Config::default());
        let row = format_row("chr1:100-112", &c, 9, 2);
        assert_eq!(row, "chr1:100-112\theterozygous\t108\t87\t9\t2\t1\t1\t0.670\t15\t36\t0\t14\t1\t2\n");
        assert_eq!(row.trim_end_matches('\n').split('\t').count(), OUTPUT_HEADER.trim_end().split('\t').count());
    }

    #[test]
    fn simple_rows() {
        let n_cols = OUTPUT_HEADER.trim_end().split('\t').count();
        assert_eq!(n_cols, 15);
        let hc = format_simple_row("chr2:5-9", GT_HIGH_COVERAGE, 181);
        assert_eq!(hc, "chr2:5-9\thigh-coverage\t181\t0\t181\t0\t0\t0\t0.000\t0\t0\t0\t0\t0\t0\n");
        let er = format_simple_row("chr2:5-9", GT_ERROR, 0);
        assert_eq!(er, "chr2:5-9\terror\t0\t0\t0\t0\t0\t0\t0.000\t0\t0\t0\t0\t0\t0\n");
        let nc = format_simple_row("chr2:5-9", GT_NO_COVERAGE, 3);
        assert_eq!(nc, "chr2:5-9\tno-coverage\t0\t0\t3\t0\t0\t0\t0.000\t0\t0\t0\t0\t0\t0\n");
        for r in [hc, er, nc] {
            assert!(r.ends_with('\n'));
            assert_eq!(r.trim_end_matches('\n').split('\t').count(), n_cols);
        }
    }
}
