//! Compile the vendored edlib (vendor/edlib, v1.2.7, MIT) -- the C++ library that
//! python-edlib 1.3.9 wraps. Using the identical library (rather than a re-implementation)
//! guarantees identical edit distances, start locations and CIGAR paths, on which the
//! indel-aware consensus and the lenient dedup depend bit-for-bit.
fn main() {
    println!("cargo:rerun-if-changed=vendor/edlib/edlib.cpp");
    println!("cargo:rerun-if-changed=vendor/edlib/edlib.h");
    cc::Build::new()
        .cpp(true)
        .file("vendor/edlib/edlib.cpp")
        .include("vendor/edlib")
        .flag_if_supported("-std=c++11")
        .opt_level(3)
        .warnings(false)
        .compile("edlib");
}
