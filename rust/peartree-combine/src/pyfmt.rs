//! Python-numerics compatibility helpers.
//!
//! IMPLEMENTED (architect, phase 1). Use these instead of ad-hoc arithmetic wherever the
//! Python source uses the corresponding builtin -- the outputs depend on them bit-for-bit
//! (SPEC.md "Numerics").

/// CPython >= 3.12 builtin `sum()` over floats (Neumaier-compensated; the farm runs
/// python-3.12.0). `sum(v for k, v in x.w.items() if k != best)` in
/// indel_consensus._decide is the only float `sum()` on an output path. Items must be fed in
/// the Python iteration order (dict insertion order).
///
/// CPython: an int start (0) plus the first float is exact, then the compensated loop runs;
/// that is identical to running the loop from 0.0.
pub fn py_sum(items: impl IntoIterator<Item = f64>) -> f64 {
    let mut f = 0.0f64;
    let mut c = 0.0f64;
    for x in items {
        let t = f + x;
        if f.abs() >= x.abs() {
            c += (f - t) + x;
        } else {
            c += (x - t) + f;
        }
        f = t;
    }
    if c != 0.0 && c.is_finite() {
        f += c;
    }
    f
}

/// Python `round(x)` for a float -> int: round-half-to-even on the exact double value.
pub fn py_round(x: f64) -> i64 {
    x.round_ties_even() as i64
}

/// Python `int(x)` for a float: truncation toward zero.
pub fn py_int(x: f64) -> i64 {
    x.trunc() as i64
}

/// `DedupParams.budget(n) = max(max_edit, int(-(-max_edit_frac * n // 1)))`
/// (float floor division by 1 == floor, so this is `ceil(frac * n)` on the IEEE product).
pub fn dedup_budget(max_edit: i64, max_edit_frac: f64, n: usize) -> i64 {
    let v = -((-max_edit_frac * n as f64).floor());
    max_edit.max(v as i64)
}

/// indel_consensus._median: middle element; for an even count
/// `int((a + b) / 2.0 + 0.5)` of the two middle elements.
pub fn py_median(xs: &[i64]) -> i64 {
    let mut s = xs.to_vec();
    s.sort_unstable();
    let m = s.len();
    assert!(m > 0, "median of empty list");
    if m % 2 == 1 {
        s[m / 2]
    } else {
        py_int((s[m / 2 - 1] + s[m / 2]) as f64 / 2.0 + 0.5)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_half_even() {
        assert_eq!(py_round(0.5), 0);
        assert_eq!(py_round(1.5), 2);
        assert_eq!(py_round(2.5), 2);
        assert_eq!(py_round(-0.5), 0);
        assert_eq!(py_round(2.4999), 2);
    }

    #[test]
    fn neumaier() {
        // naive left-to-right summation gives 0.0 here; CPython 3.12 sum() gives 1.0
        assert_eq!(py_sum([1e100, 1.0, -1e100]), 1.0);
        assert_eq!(py_sum([0.1, 0.2, 0.3]), 0.6);
        assert_eq!(py_sum(std::iter::empty()), 0.0);
    }

    #[test]
    fn budget() {
        assert_eq!(dedup_budget(3, 0.02, 0), 3);
        assert_eq!(dedup_budget(3, 0.02, 150), 3);
        assert_eq!(dedup_budget(3, 0.02, 151), 4); // 3.02 -> ceil 4
        assert_eq!(dedup_budget(3, 0.02, 200), 4);
    }

    #[test]
    fn median() {
        assert_eq!(py_median(&[3, 1, 2]), 2);
        assert_eq!(py_median(&[1, 2]), 2);
        assert_eq!(py_median(&[1, 4]), 3);
    }
}
