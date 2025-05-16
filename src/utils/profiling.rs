use anyhow::Result;
use std::time::Duration;

#[cfg(feature = "profiling")]
use pprof::ProfilerGuard;

/// Start CPU profiling with pprof.
/// 
/// Returns a guard that must be kept alive during the profiling period.
/// When the guard is dropped, profiling stops.
#[cfg(feature = "profiling")]
pub fn start_profiling(name: &str) -> Result<ProfilerGuard<'static>> {
    // Create profiles directory if it doesn't exist
    let profile_dir = Path::new("profiles");
    fs::create_dir_all(profile_dir)?;
    
    // Create profiler with 100 Hz sampling rate
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(100)
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build()?;
    
    println!("Profiling started: {}", name);
    Ok(guard)
}

/// Finish profiling and write results to files.
#[cfg(feature = "profiling")]
pub fn finish_profiling(guard: Option<ProfilerGuard<'static>>, name: &str) -> Result<()> {
    if let Some(guard) = guard {
        // Get report
        let report = guard.report().build()?;
        
        // Create profiles directory if it doesn't exist
        let profile_dir = PathBuf::from("profiles");
        fs::create_dir_all(&profile_dir)?;
        
        // Write flamegraph
        let flamegraph_path = profile_dir.join(format!("{}_flamegraph.svg", name));
        let mut file = fs::File::create(&flamegraph_path)?;
        report.flamegraph(&mut file)?;
        println!("Flamegraph written to: {}", flamegraph_path.display());
        
        // Write pprof proto file - use write_protobuf instead of pprof_proto
        let proto_path = profile_dir.join(format!("{}_profile.pb", name));
        let file = fs::File::create(&proto_path)?;
        report.write_protobuf(file)?;
        println!("Profile proto written to: {}", proto_path.display());
    }
    Ok(())
}

// Dummy implementations for when profiling is disabled
#[cfg(not(feature = "profiling"))]
pub fn start_profiling(_name: &str) -> Result<()> {
    println!("Profiling not enabled. Compile with --features=profiling to enable.");
    Ok(())
}

#[cfg(not(feature = "profiling"))]
pub fn finish_profiling(_guard: Option<()>, _name: &str) -> Result<()> {
    Ok(())
}

/// Time a function execution and print the elapsed time.
pub fn time_function<F, T>(name: &str, f: F) -> T 
where 
    F: FnOnce() -> T 
{
    let start = std::time::Instant::now();
    let result = f();
    let elapsed = start.elapsed();
    println!("{} completed in {:.2?}", name, elapsed);
    result
}

/// Format a duration as a human-readable string.
pub fn format_duration(duration: Duration) -> String {
    if duration.as_secs() > 60 {
        let minutes = duration.as_secs() / 60;
        let seconds = duration.as_secs() % 60;
        let millis = duration.subsec_millis();
        format!("{}m {}s {}ms", minutes, seconds, millis)
    } else if duration.as_secs() > 0 {
        let seconds = duration.as_secs();
        let millis = duration.subsec_millis();
        format!("{}s {}ms", seconds, millis)
    } else {
        let millis = duration.as_millis();
        if millis > 0 {
            format!("{}ms", millis)
        } else {
            format!("{}µs", duration.as_micros())
        }
    }
}

/// Measure memory usage using pprof.
#[cfg(feature = "profiling")]
pub fn measure_memory() -> Result<(usize, usize)> {
    let mut allocated: usize = 0;
    let mut active: usize = 0;
    
    unsafe {
        tikv_jemallocator::stats::get_allocator_stats(&mut allocated, &mut active)?;
    }
    
    Ok((allocated, active))
}

/// Print current memory usage.
#[cfg(feature = "profiling")]
pub fn print_memory_usage(label: &str) -> Result<()> {
    let (allocated, active) = measure_memory()?;
    println!(
        "Memory usage at {}: allocated={}MB, active={}MB", 
        label, 
        allocated / 1_000_000, 
        active / 1_000_000
    );
    Ok(())
}

#[cfg(not(feature = "profiling"))]
pub fn print_memory_usage(label: &str) -> Result<()> {
    println!("Memory usage at {} not available (jemalloc not enabled)", label);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread;
    
    #[test]
    fn test_format_duration() {
        assert_eq!(format_duration(Duration::from_secs(125)), "2m 5s 0ms");
        assert_eq!(format_duration(Duration::from_secs(45)), "45s 0ms");
        assert_eq!(format_duration(Duration::from_millis(750)), "750ms");
        assert_eq!(format_duration(Duration::from_micros(500)), "500µs");
    }
    
    #[test]
    fn test_time_function() {
        let result = time_function("sleep test", || {
            thread::sleep(Duration::from_millis(10));
            42
        });
        assert_eq!(result, 42);
    }
} 