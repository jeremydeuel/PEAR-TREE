//! `peartree-combine` -- drop-in for `python src/main.py --step combine_insertions ...`.

#[cfg(feature = "jemalloc")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let args = match peartree_combine::pipeline::Args::parse(&argv) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("peartree-combine: {e}");
            std::process::exit(1);
        }
    };
    if let Err(e) = peartree_combine::pipeline::run(&args) {
        eprintln!("peartree-combine: {e}");
        std::process::exit(1);
    }
}
