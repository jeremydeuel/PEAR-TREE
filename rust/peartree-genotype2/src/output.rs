//! Output row formatting (SPEC "Output"). Owner: C.

use crate::types::Call;

/// Format one output row (ends with '\n'). `coverage` and `n_disc` come from the driver.
pub fn format_row(name: &str, call: &Call, coverage: i64, n_disc: i64) -> String {
    let _ = (name, call, coverage, n_disc);
    todo!("owner C: output::format_row")
}

/// Rows for the driver-decided states: `high-coverage` (coverage known), `error`, `no-coverage`.
pub fn format_simple_row(name: &str, genotype: &'static str, coverage: i64) -> String {
    let _ = (name, genotype, coverage);
    todo!("owner C: output::format_simple_row")
}
