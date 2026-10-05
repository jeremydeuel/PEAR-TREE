//! Stage driver -- python `combine_insertions()` (combine_insertions.py:113-365) and the
//! `--step combine_insertions` branch of src/main.py:148-188. OWNER: P6.
//! SPEC.md §1 (CLI, outputs) and §3-§6 (stages in order).

use std::path::PathBuf;

/// Parsed command line (SPEC.md §1.1).
#[derive(Clone, Debug)]
pub struct Args {
    pub discovery_files: Vec<String>,
    /// `--out` with a trailing `.gz` stripped (python main.py)
    pub out_stem: String,
    pub threads: usize,
    pub config: PathBuf,
}

impl Args {
    /// Accepts the python spelling (`--step combine_insertions` ignored / validated,
    /// `--discovery_files F...`, `--out STEM`, `--threads N`) plus `--config PATH`
    /// (default `$PEARTREE_CONFIG`, else `src/config.py` relative to the cwd). Multi-value
    /// `--discovery_files` consumes arguments until the next `--option`. Also accepts
    /// `--discovery-files`.
    pub fn parse(argv: &[String]) -> Result<Args, String> {
        todo!("P6: SPEC.md §1.1")
    }
}

/// Run the whole step. Stage order (SPEC.md §3-§6):
///  1. load config, open genome (error if missing), validate inputs (python main.py checks:
///     files exist, samtools/bowtie2 executables exist, `{bowtie2_index}.1.bt2` exists,
///     threads <= cpu count), reject duplicate basenames;
///  2. parse discovery files in parallel, import loop / exclusion in command-line order (§3.1);
///  3. far_pair_strict -> discovery_breakpoints over all accepted records (§3.2);
///  4. intersect_insertions (§3.3); 5. filter_dense_regions(100, 4) (§3.4);
///  6. evidence (if any accepted sidecar) -> apply_evidence (§4);
///  7. write `<stem>.fq.gz` consensus FASTQ, bowtie2 end-to-end -> `<stem>.bam` unless it exists,
///     clean-remap filter (§5.1);
///  8. unless `<stem>.insertionsonly.bam` exists: write clip FASTQ to `<stem>.fq.gz`, bowtie2
///     local; clipped-remap filter with liftover (§5.2);
///  9. evidence -> absorb_one_sided (§4.7);
/// 10. write `<stem>.combined.txt.gz` (§6.1); 11. `<stem>.combined.splice.tsv` (§6.6);
/// 12. evidence -> write evidence TSV + reads FASTA (§6.4-6.5), cleanup shard dir;
/// 13. write `<stem>.genotyping.txt.gz` (§6.3).
pub fn run(args: &Args) -> Result<(), String> {
    todo!("P6: SPEC.md §1-§6")
}
