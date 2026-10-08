//! Python numerics / formatting compatibility. FOUNDATION (implemented).
//!
//! Use these wherever the Python source uses the corresponding builtin -- output strings
//! (record.row(), detail_string(), score points) depend on them byte-for-byte.
//!
//! * [`py_sum`] -- builtin `sum()` over floats (CPython >= 3.12: Neumaier-compensated; the farm
//!   and the golden venvs run 3.12)
//! * [`py_round`] -- `round(x, ndigits)` (correctly rounded, half-to-even on the exact value)
//! * [`py_repr`] -- `str(float)` / `repr(float)`
//! * [`py_g`] -- `format(x, ".{p}g")`, `f"{x:+g}"`
//! * [`py_int`] -- `int(float)` (truncation)

/// CPython >= 3.12 builtin `sum()` over floats, fed in Python iteration order.
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

/// `int(x)` for a float: truncation toward zero.
pub fn py_int(x: f64) -> i64 {
    x.trunc() as i64
}

/// `round(x, nd)` for a float (nd >= 0): the exact binary value rounded half-to-even to `nd`
/// decimals, as CPython does via `_Py_dg_dtoa` mode 3.
pub fn py_round(x: f64, nd: usize) -> f64 {
    if !x.is_finite() {
        return x;
    }
    // Rust's fixed-precision float formatting is exact and rounds half-to-even on exact ties
    // (verified in the tests against CPython values).
    format!("{:.*}", nd, x).parse::<f64>().unwrap_or(x)
}

/// (digits, exponent) of the shortest round-trip representation: x = 0.d1d2d3... * 10^(exp+1),
/// i.e. `d1.d2d3 * 10^exp`.
fn shortest_digits(x: f64) -> (String, i32) {
    let s = format!("{:e}", x.abs()); // e.g. "1.2345e-5", "1e16"
    let (m, e) = s.split_once('e').unwrap();
    let digits: String = m.chars().filter(|c| *c != '.').collect();
    (digits, e.parse().unwrap())
}

/// `repr(x)` / `str(x)` of a Python float.
pub fn py_repr(x: f64) -> String {
    if x.is_nan() {
        return "nan".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".into() } else { "-inf".into() };
    }
    let neg = x.is_sign_negative();
    if x == 0.0 {
        return if neg { "-0.0".into() } else { "0.0".into() };
    }
    let (digits, exp) = shortest_digits(x);
    let mut out = String::new();
    if neg {
        out.push('-');
    }
    if (-4..16).contains(&exp) {
        let nd = digits.len() as i32;
        if exp >= 0 {
            let int_len = (exp + 1) as usize;
            if digits.len() <= int_len {
                out.push_str(&digits);
                out.push_str(&"0".repeat(int_len - digits.len()));
                out.push_str(".0");
            } else {
                out.push_str(&digits[..int_len]);
                out.push('.');
                out.push_str(&digits[int_len..]);
            }
        } else {
            out.push_str("0.");
            out.push_str(&"0".repeat((-exp - 1) as usize));
            out.push_str(&digits);
        }
        let _ = nd;
    } else {
        out.push_str(&digits[..1]);
        if digits.len() > 1 {
            out.push('.');
            out.push_str(&digits[1..]);
        }
        out.push('e');
        out.push(if exp < 0 { '-' } else { '+' });
        out.push_str(&format!("{:02}", exp.abs()));
    }
    out
}

/// `format(x, f".{prec}g")`, with a leading '+' for non-negative values when `plus`
/// (`f"{x:+g}"` is `py_g(x, 6, true)`).
pub fn py_g(x: f64, prec: usize, plus: bool) -> String {
    let p = prec.max(1);
    let sign = if x.is_sign_negative() && !x.is_nan() {
        "-"
    } else if plus {
        "+"
    } else {
        ""
    };
    if x.is_nan() {
        return format!("{}nan", if plus { "+" } else { "" });
    }
    if x.is_infinite() {
        return format!("{sign}inf");
    }
    let a = x.abs();
    if a == 0.0 {
        return format!("{sign}0");
    }
    // round to p significant digits first; the exponent AFTER rounding decides the notation
    let s = format!("{:.*e}", p - 1, a);
    let (m, e) = s.split_once('e').unwrap();
    let exp: i32 = e.parse().unwrap();
    let body = if exp >= -4 && exp < p as i32 {
        let decimals = (p as i32 - 1 - exp).max(0) as usize;
        strip_zeros(&format!("{:.*}", decimals, a))
    } else {
        format!("{}e{}{:02}", strip_zeros(m), if exp < 0 { '-' } else { '+' }, exp.abs())
    };
    format!("{sign}{body}")
}

fn strip_zeros(s: &str) -> String {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repr_like_python() {
        // python: str(x)
        for (x, want) in [
            (0.0, "0.0"),
            (1.0, "1.0"),
            (0.9876, "0.9876"),
            (0.1 + 0.2, "0.30000000000000004"),
            (1e-5, "1e-05"),
            (0.0001, "0.0001"),
            (1e16, "1e+16"),
            (123456789012345.0, "123456789012345.0"),
            (-2.5, "-2.5"),
            (1.5e300, "1.5e+300"),
            (0.98, "0.98"),
        ] {
            assert_eq!(py_repr(x), want, "{x}");
        }
    }

    #[test]
    fn g_like_python() {
        // python: format(x, ".4g") / f"{x:+g}" / f"{x:g}"
        assert_eq!(py_g(0.98765, 4, false), "0.9877");
        assert_eq!(py_g(0.5, 4, false), "0.5");
        assert_eq!(py_g(0.0, 4, false), "0");
        assert_eq!(py_g(0.00012345, 4, false), "0.0001234");
        assert_eq!(py_g(0.000012345, 4, false), "1.234e-05");
        assert_eq!(py_g(12.0, 6, false), "12");
        assert_eq!(py_g(12.5, 6, false), "12.5");
        assert_eq!(py_g(1234567.0, 6, false), "1.23457e+06");
        assert_eq!(py_g(2.0, 6, true), "+2");
        assert_eq!(py_g(-1.5, 6, true), "-1.5");
        assert_eq!(py_g(0.5, 6, true), "+0.5");
        assert_eq!(py_g(-4.0, 6, true), "-4");
        assert_eq!(py_g(0.99995, 4, false), "1");
        assert_eq!(py_g(100000.0, 6, false), "100000");
        assert_eq!(py_g(1000000.0, 6, false), "1e+06");
    }

    #[test]
    fn round_like_python() {
        // python: round(x, 2)
        assert_eq!(py_round(0.125, 2), 0.12);
        assert_eq!(py_round(0.375, 2), 0.38);
        assert_eq!(py_round(2.675, 2), 2.67);
        assert_eq!(py_round(7.004999, 2), 7.0);
        assert_eq!(py_round(-1.005, 2), -1.0);
        assert_eq!(py_round(0.98765, 4), 0.9877);
        assert_eq!(py_round(2.5, 0), 2.0);
        assert_eq!(py_round(3.5, 0), 4.0);
    }

    #[test]
    fn neumaier_sum() {
        assert_eq!(py_sum([0.1, 0.2, 0.3]), 0.6);
        assert_eq!(py_sum([1e100, 1.0, -1e100]), 1.0);
    }
}
