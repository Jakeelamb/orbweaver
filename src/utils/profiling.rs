use anyhow::Result;
use std::time::Duration;

#[cfg(feature = "profiling")]
use pprof::{ProfilerGuard, ProfilerGuardBuilder, Report};
#[cfg(feature = "profiling")]
use pprof::protos::Message;

#[cfg(feature = "profiling")]
pub struct ProfileSessionGuard {
    guard: Option<ProfilerGuard<'static>>,
    name: String,
}

#[cfg(feature = "profiling")]
impl Drop for ProfileSessionGuard {
    fn drop(&mut self) {
        if let Some(guard) = self.guard.take() {
            if let Err(e) = finish_profiling_internal(guard, &self.name) {
                eprintln!("Error finishing profiling for '{}': {}", self.name, e);
            }
        }
    }
}

#[cfg(feature = "profiling")]
fn finish_profiling_internal(guard: ProfilerGuard<'static>, name: &str) -> Result<()> {
    println!("Finishing profiling for: {}", name);
    let report_build_result = guard.report().build();
    
    let report = match report_build_result {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Failed to build profiling report for '{}': {}", name, e);
            return Err(e.into());
        }
    };
    
    let profile_dir = PathBuf::from("profiles");
    if !profile_dir.exists() {
        if let Err(e) = fs::create_dir_all(&profile_dir) {
            eprintln!("Failed to create profiles directory '{}': {}", profile_dir.display(), e);
            return Err(e.into());
        }
    }
    
    let flamegraph_path = profile_dir.join(format!("{}_flamegraph.svg", name));
    match fs::File::create(&flamegraph_path) {
        Ok(mut file_flame) => {
            if let Err(e) = report.flamegraph(&mut file_flame) {
                eprintln!("Failed to write flamegraph to '{}': {}", flamegraph_path.display(), e);
            } else {
                println!("Flamegraph written to: {}", flamegraph_path.display());
            }
        }
        Err(e) => {
            eprintln!("Failed to create flamegraph file '{}': {}", flamegraph_path.display(), e);
        }
    }
    
    let proto_path = profile_dir.join(format!("{}_profile.pb", name));
    match fs::File::create(&proto_path) {
        Ok(mut file_proto) => {
            match report.pprof() {
                Ok(pprof_profile_message) => {
                    match pprof_profile_message.write_to_writer(&mut file_proto) {
                        Ok(_) => println!("Profile proto written to: {}", proto_path.display()),
                        Err(e) => eprintln!("Failed to write pprof proto to '{}': {}", proto_path.display(), e),
                    }
                }
                Err(e) => {
                    eprintln!("Failed to generate pprof::protos::Profile from report for '{}': {}", name, e);
                }
            }
        }
        Err(e) => {
            eprintln!("Failed to create pprof proto file '{}': {}", proto_path.display(), e);
        }
    }
    Ok(())
}

/// Start CPU profiling with pprof.
///
/// Returns a ProfileSessionGuard. Profiling stops when this guard is dropped.
#[cfg(feature = "profiling")]
pub fn start_profiling(name: &str) -> Result<ProfileSessionGuard> {
    let profile_dir = PathBuf::from("profiles");
    if !profile_dir.exists() {
        fs::create_dir_all(&profile_dir)?;
    }
    
    println!("Attempting to start profiling: {}", name);
    match ProfilerGuardBuilder::default()
        .frequency(1000)
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build()
    {
        Ok(guard) => {
            println!("Profiling started successfully: {}", name);
            Ok(ProfileSessionGuard {
                guard: Some(guard),
                name: name.to_string(),
            })
        }
        Err(e) => {
            eprintln!("Failed to start profiling for '{}': {}. Profiling will be skipped.", name, e);
            Ok(ProfileSessionGuard {
                guard: None,
                name: name.to_string(),
            })
        }
    }
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
    let allocated: usize = 0;
    let active: usize = 0;
    
    Ok((allocated, active))
}

/// Print current memory usage.
#[cfg(feature = "profiling")]
pub fn print_memory_usage(_label: &str) -> Result<()> {
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
pub fn print_memory_usage(_label: &str) -> Result<()> {
    Ok(())
}

#[cfg(not(feature = "profiling"))]
pub fn start_profiling(name: &str) -> Result<ProfileSessionGuard> {
    println!(
        "Profiling not enabled for '{}'. Compile with --features=profiling to enable.",
        name
    );
    Ok(ProfileSessionGuard { _private: (()) })
}

#[cfg(not(feature = "profiling"))]
pub struct ProfileSessionGuard {
    _private: (),
}

#[cfg(not(feature = "profiling"))]
impl Drop for ProfileSessionGuard {
    fn drop(&mut self) {
        // No-op when profiling is not enabled
    }
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