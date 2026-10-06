//! Output row formatting (`types::OUTPUT_HEADER`): numeric columns + a status word.

use crate::types::{Call, LocusKind, Status};

/// Format one `ok` row (ends with '\n'). `depth` and `n_disc` come from the driver.
pub fn format_row(name: &str, kind: LocusKind, call: &Call, depth: i64, n_disc: i64) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        name,
        kind.as_str(),
        Status::Ok.as_str(),
        depth,
        call.n_alt,
        call.n_ref,
        call.n_uninf,
        call.n_art,
        n_disc,
        call.n_alt_l,
        call.n_alt_r,
        call.vaf,
        call.post[0],
        call.post[1],
        call.post[2],
        call.pl[0],
        call.pl[1],
        call.pl[2],
        call.gq,
        call.score_alt,
        call.score_ref
    )
}

/// Row for a locus without a model result (`no_reads`, `high_coverage`, `error`): the depth
/// (capped for `high_coverage`), zeros everywhere else; the posterior columns are empty so a
/// consumer cannot mistake them for a measurement.
pub fn format_simple_row(name: &str, kind: LocusKind, status: Status, depth: i64) -> String {
    debug_assert!(status != Status::Ok);
    format!("{name}\t{}\t{}\t{depth}\t0\t0\t0\t0\t0\t0\t0\t0.000\t\t\t\t0\t0\t0\t0\t0\t0\n", kind.as_str(), status.as_str())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Config;
    use crate::model::call_locus;
    use crate::types::{AltSide, ReadClass, ReadObs, OUTPUT_HEADER};

    fn obs(class: ReadClass, llr: f64) -> ReadObs {
        ReadObs { ll_ref: -40.0, ll_alt: -40.0 + llr, class, explained_frac: 1.0, crosses_junction: true, alt_side: AltSide::None }
    }

    fn n_cols() -> usize {
        OUTPUT_HEADER.trim_end().split('\t').count()
    }

    #[test]
    fn ok_row_has_every_column() {
        let reads = [obs(ReadClass::Alt, 10.0), obs(ReadClass::Alt, 15.0), obs(ReadClass::Ref, -20.0), obs(ReadClass::Uninformative, 1.0), obs(ReadClass::Unexplained, 0.0)];
        let c = call_locus(&reads, 2, &Config::default());
        let row = format_row("chr1:100-112", LocusKind::Tsd, &c, 9, 2);
        let f: Vec<&str> = row.trim_end_matches('\n').split('\t').collect();
        assert_eq!(f.len(), n_cols(), "{row}");
        assert_eq!(&f[..11], &["chr1:100-112", "TSD", "ok", "9", "2", "1", "1", "1", "2", "0", "0"]);
        assert_eq!(f[11], "0.720");
        let p: Vec<f64> = f[12..15].iter().map(|x| x.parse().unwrap()).collect();
        assert!((p.iter().sum::<f64>() - 1.0).abs() < 2e-4);
        assert_eq!(f[16], "0"); // pl_het is the best
        assert_eq!((f[19], f[20]), ("108", "87"));
    }

    #[test]
    fn simple_rows() {
        for st in [Status::NoReads, Status::HighCoverage, Status::Error] {
            let row = format_simple_row("x:1-2", LocusKind::Blunt, st, 181);
            let f: Vec<&str> = row.trim_end_matches('\n').split('\t').collect();
            assert_eq!(f.len(), n_cols(), "{row}");
            assert_eq!((f[1], f[2], f[3]), ("BLUNT", st.as_str(), "181"));
            assert_eq!((f[12], f[13], f[14]), ("", "", ""));
        }
    }
}
