use std::time::{Duration, Instant};
use std::fmt;

/// A simple utility for tracking and reporting progress of long-running operations
pub struct ProgressTracker {
    /// Total items to process
    total: usize,
    /// Items processed so far
    current: usize,
    /// When the operation started
    start_time: Instant,
    /// Last update time
    last_update: Instant,
    /// Update interval
    update_interval: Duration,
    /// Operation name
    title: String,
}

impl ProgressTracker {
    /// Creates a new progress tracker
    pub fn new(total: usize, title: &str) -> Self {
        let now = Instant::now();
        Self {
            total,
            current: 0,
            start_time: now,
            last_update: now,
            update_interval: Duration::from_secs(1),
            title: title.to_string(),
        }
    }
    
    /// Set how frequently progress should be reported
    pub fn with_update_interval(mut self, interval: Duration) -> Self {
        self.update_interval = interval;
        self
    }
    
    /// Increment the progress counter, returns true if it's time to display an update
    pub fn increment(&mut self) -> bool {
        self.current += 1;
        let now = Instant::now();
        
        if now.duration_since(self.last_update) >= self.update_interval {
            self.last_update = now;
            true
        } else {
            false
        }
    }
    
    /// Marks the operation as complete
    pub fn finish(&mut self) {
        self.current = self.total;
        self.last_update = Instant::now();
    }
    
    /// Get a status message showing the current progress
    pub fn status(&self) -> String {
        let elapsed = self.last_update.duration_since(self.start_time);
        let progress_pct = (self.current as f64 / self.total as f64) * 100.0;
        
        let items_per_sec = if elapsed.as_secs_f64() > 0.0 {
            self.current as f64 / elapsed.as_secs_f64()
        } else {
            0.0
        };
        
        // Estimate time remaining
        let eta = if self.current > 0 && items_per_sec > 0.0 {
            let remaining_items = self.total - self.current;
            let seconds_remaining = remaining_items as f64 / items_per_sec;
            
            if seconds_remaining < 60.0 {
                format!("{:.1}s", seconds_remaining)
            } else if seconds_remaining < 3600.0 {
                format!("{:.1}m", seconds_remaining / 60.0)
            } else {
                format!("{:.1}h", seconds_remaining / 3600.0)
            }
        } else {
            "Unknown".to_string()
        };
        
        format!(
            "{}: {}/{} ({:.1}%) - {:.1} items/sec, ETA: {}",
            self.title,
            self.current,
            self.total,
            progress_pct,
            items_per_sec,
            eta
        )
    }
}

impl fmt::Display for ProgressTracker {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.status())
    }
}

/// Format a duration as HH:MM:SS
fn format_duration(duration: Duration) -> String {
    let total_seconds = duration.as_secs();
    let hours = total_seconds / 3600;
    let minutes = (total_seconds % 3600) / 60;
    let seconds = total_seconds % 60;
    
    format!("{:02}:{:02}:{:02}", hours, minutes, seconds)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_progress_tracker() {
        let mut tracker = ProgressTracker::new(100, "Test");
        assert_eq!(tracker.current, 0);
        
        // Test increment
        let first_update = tracker.increment();
        assert_eq!(tracker.current, 1);
        
        // Since we're not waiting, update_interval won't have elapsed
        assert!(!first_update);
        
        // Force an update by manipulating the last_update time
        tracker.last_update = Instant::now() - Duration::from_secs(2);
        let should_update = tracker.increment();
        assert!(should_update);
        
        // Test finish
        tracker.finish();
        assert_eq!(tracker.current, 100);
        
        // Verify status string contains expected elements
        let status = tracker.status();
        assert!(status.contains("Test:"));
        assert!(status.contains("100/100"));
        assert!(status.contains("100.0%"));
    }
    
    #[test]
    fn test_format_duration() {
        let duration = Duration::from_secs(3661); // 1h 1m 1s
        assert_eq!(format_duration(duration), "01:01:01");
        
        let duration2 = Duration::from_secs(59); // 59s
        assert_eq!(format_duration(duration2), "00:00:59");
    }
} 