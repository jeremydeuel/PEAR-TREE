//! `PEARTREE_MEMLOG=1` diagnostics: stage timestamps, current and peak RSS on stderr. Never
//! affects outputs.

use std::sync::OnceLock;
use std::time::Instant;

static T0: OnceLock<Instant> = OnceLock::new();

/// start the clock (first call wins)
pub fn start() {
    let _ = T0.get_or_init(Instant::now);
}

pub fn enabled() -> bool {
    std::env::var_os("PEARTREE_MEMLOG").is_some()
}

/// peak RSS in MB (getrusage), -1 on failure
pub fn peak_rss_mb() -> f64 {
    #[repr(C)]
    struct RUsage {
        times: [i64; 4],
        maxrss: i64,
        rest: [i64; 13],
    }
    extern "C" {
        fn getrusage(who: i32, usage: *mut RUsage) -> i32;
    }
    let mut ru = RUsage { times: [0; 4], maxrss: 0, rest: [0; 13] };
    // SAFETY: plain libc call with a correctly sized out-struct (struct rusage = 2 timevals + 14 longs)
    let rc = unsafe { getrusage(0, &mut ru) };
    // ru_maxrss: bytes on macOS, KiB on Linux
    if rc != 0 {
        -1.0
    } else if cfg!(target_os = "macos") {
        ru.maxrss as f64 / 1048576.0
    } else {
        ru.maxrss as f64 / 1024.0
    }
}

/// current RSS in MB, -1 when unknown
#[cfg(target_os = "linux")]
pub fn cur_rss_mb() -> f64 {
    let Ok(s) = std::fs::read_to_string("/proc/self/statm") else { return -1.0 };
    let pages: f64 = s.split_whitespace().nth(1).and_then(|v| v.parse().ok()).unwrap_or(-1.0);
    if pages < 0.0 {
        return -1.0;
    }
    pages * 4096.0 / 1048576.0
}

#[cfg(target_os = "macos")]
pub fn cur_rss_mb() -> f64 {
    // mach_task_basic_info (MACH_TASK_BASIC_INFO = 20)
    #[repr(C)]
    struct Info {
        virtual_size: u64,
        resident_size: u64,
        resident_size_max: u64,
        user_time: [i32; 2],
        system_time: [i32; 2],
        policy: i32,
        suspend_count: i32,
    }
    extern "C" {
        static mach_task_self_: u32;
        fn task_info(task: u32, flavor: i32, info: *mut Info, count: *mut u32) -> i32;
    }
    let mut info = std::mem::MaybeUninit::<Info>::zeroed();
    let mut count = (std::mem::size_of::<Info>() / 4) as u32;
    // SAFETY: mach call with a correctly sized out-struct
    let rc = unsafe { task_info(mach_task_self_, 20, info.as_mut_ptr(), &mut count) };
    if rc != 0 {
        return -1.0;
    }
    // SAFETY: filled by task_info
    unsafe { info.assume_init() }.resident_size as f64 / 1048576.0
}

#[cfg(not(any(target_os = "linux", target_os = "macos")))]
pub fn cur_rss_mb() -> f64 {
    -1.0
}

/// one `[memlog]` line (no-op unless PEARTREE_MEMLOG is set)
pub fn memlog(stage: &str) {
    if !enabled() {
        return;
    }
    let t = T0.get_or_init(Instant::now).elapsed().as_secs_f64();
    eprintln!("[memlog] {t:8.2}s  peak {:7.0} MB  cur {:7.0} MB  {stage}", peak_rss_mb(), cur_rss_mb());
}
