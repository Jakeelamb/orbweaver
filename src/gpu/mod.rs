//! GPU acceleration module for digram finding and suffix array construction.
//!
//! This module provides GPU-accelerated implementations of key algorithms used in grammar construction,
//! particularly for processing large genomic sequences efficiently.

use anyhow::{Result, anyhow};
use ocl::{Platform, Device, Context, Queue, Program};
use std::time::Instant;

// Import the kernels from the workspace crate
extern crate orbweaver_kernels;

pub mod digram;
pub mod suffix_array;

/// Represents a GPU context for accelerating specific operations
/// in the Orbweaver genomic grammar construction process.
#[derive(Debug)]
pub struct GpuContext {
    pub platform: Platform,
    pub device: Device,
    pub context: Context,
    pub queue: Queue,
    pub program: Option<Program>,
}

/// Allow cloning the GPU context
impl Clone for GpuContext {
    fn clone(&self) -> Self {
        // Try to create a fresh context instead of attempting to clone OpenCL objects
        match GpuContext::new() {
            Ok(ctx) => ctx,
            Err(e) => {
                eprintln!("Warning: Failed to clone GPU context: {}. Creating placeholder instead.", e);
                
                // Create a placeholder with default values
                // This should only be used for function signatures and never accessed directly
                let default_platform = Platform::default();
                
                // Create a minimal device, context, and queue
                let device = match Device::first(default_platform) {
                    Ok(d) => d,
                    Err(_) => panic!("Failed to create device for placeholder GPU context"),
                };
                
                let context = match Context::builder()
                    .platform(default_platform)
                    .devices(device)
                    .build() {
                    Ok(c) => c,
                    Err(_) => panic!("Failed to create context for placeholder GPU context"),
                };
                
                let queue = match Queue::new(&context, device, None) {
                    Ok(q) => q,
                    Err(_) => panic!("Failed to create queue for placeholder GPU context"),
                };
                
                // Return a minimally valid context
                Self {
                    platform: default_platform,
                    device,
                    context,
                    queue,
                    program: None,
                }
            }
        }
    }
}

impl GpuContext {
    /// Create a new GPU context with OpenCL
    pub fn new() -> Result<Self> {
        let start = Instant::now();

        println!("Initializing GPU context...");

        // Get all available platforms
        let platforms = Platform::list();
        if platforms.is_empty() {
            return Err(anyhow!("No OpenCL platforms found"));
        }

        // Choose the first platform (can be made smarter with preference logic)
        let platform = platforms[0];
        println!("Using platform: {}", platform.name()?);

        // Get devices for this platform
        let devices = Device::list_all(platform)?;
        if devices.is_empty() {
            return Err(anyhow!("No OpenCL devices found for platform"));
        }

        // Prefer GPUs, but fall back to CPU if necessary
        let mut gpu_device = None;
        let mut cpu_device = None;

        for device in &devices {
            let device_type = device.info(ocl::enums::DeviceInfo::Type)?;
            let device_name = device.info(ocl::enums::DeviceInfo::Name)?;
            
            if device_type.to_string().contains("GPU") {
                gpu_device = Some(*device);
                println!("Found GPU device: {}", device_name);
            } else if device_type.to_string().contains("CPU") {
                cpu_device = Some(*device);
                println!("Found CPU device: {}", device_name);
            }
        }

        // Choose GPU first, fallback to CPU
        let device = gpu_device.unwrap_or_else(|| {
            println!("No GPU found, falling back to CPU");
            cpu_device.expect("No CPU device found either")
        });

        let device_name = device.info(ocl::enums::DeviceInfo::Name)?;
        println!("Selected device: {}", device_name);

        // Log device capabilities
        println!("Device version: {}", device.info(ocl::enums::DeviceInfo::Version)?);
        println!("OpenCL C version: {}", device.info(ocl::enums::DeviceInfo::OpenclCVersion)?);
        println!("Max compute units: {}", device.info(ocl::enums::DeviceInfo::MaxComputeUnits)?);
        let max_wg_size = device.info(ocl::enums::DeviceInfo::MaxWorkGroupSize)?;
        println!("Max work group size: {}", max_wg_size);
        
        // Get global memory size as string and convert to number
        let global_mem_size_str = device.info(ocl::enums::DeviceInfo::GlobalMemSize)?.to_string();
        let global_mem_size: u64 = global_mem_size_str.parse().unwrap_or(0);
        println!("Global memory size: {} MB", global_mem_size / (1024 * 1024));
        
        // Check minimum requirements
        let min_global_mem: u64 = 256 * 1024 * 1024; // 256 MB
        if global_mem_size < min_global_mem {
            println!("WARNING: Device has limited memory ({} MB), may encounter CL_OUT_OF_RESOURCES errors", 
                    global_mem_size / (1024 * 1024));
        }
        
        // Create context with just the selected device
        let context = Context::builder()
            .platform(platform)
            .devices(device)
            .build()?;

        // Create command queue with profiling enabled for debugging
        let queue_properties = ocl::CommandQueueProperties::new().profiling();
        let queue = Queue::new(&context, device, Some(queue_properties))?;
        
        // Create a GpuContext without program (will be initialized when needed)
        let mut gpu_context = GpuContext {
            platform,
            device,
            context,
            queue,
            program: None,
        };
        
        // Try to load kernels
        match gpu_context.load_kernels() {
            Ok(_) => println!("OpenCL kernels loaded successfully"),
            Err(e) => println!("Warning: Failed to load kernels at initialization: {}", e),
        }

        println!("GPU context initialized in {:?}", start.elapsed());
        
        Ok(gpu_context)
    }
    
    /// Load and compile the OpenCL kernels
    fn load_kernels(&mut self) -> Result<()> {
        // Get kernel source code for digram kernel
        let digram_kernel_src = orbweaver_kernels::get_digram_kernel();
        
        // Get suffix array kernel from the kernels crate
        let suffix_array_kernel_src = orbweaver_kernels::get_suffix_array_kernel();
        
        // Combine kernels with proper header (OpenCL C version pragma)
        let combined_src = format!(
            "#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable\n\
             #pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable\n\
             \n\
             {}\n\
             \n\
             {}", 
            digram_kernel_src, 
            suffix_array_kernel_src
        );
        
        // Build program with proper error handling
        let program = match Program::builder()
            .devices(self.device)
            .src(combined_src)
            .build(&self.context) {
                Ok(p) => p,
                Err(e) => {
                    eprintln!("\n###################### OPENCL PROGRAM BUILD DEBUG OUTPUT ######################");
                    
                    // Extract build log for better diagnostics
                    let build_log = e.to_string();
                    eprintln!("{}", build_log);
                    
                    eprintln!("###############################################################################\n");
                    
                    // Print warning but don't fail - we'll try CPU fallbacks
                    println!("Warning: Failed to load kernels at initialization: {}", e);
                    return Err(anyhow!("Failed to build OpenCL program: {}", e));
                }
            };
        
        // Store the successfully built program
        self.program = Some(program);
        println!("OpenCL kernels loaded successfully");
        
        Ok(())
    }
    
    /// Check if we have a GPU device available
    pub fn is_gpu() -> bool {
        // Try to initialize an OpenCL context with a GPU
        let platforms = Platform::list();
        if platforms.is_empty() {
            return false;
        }
        
        let platform = platforms[0];
        
        match Device::list_all(platform) {
            Ok(devices) => {
                for device in devices {
                    if let Ok(device_type) = device.info(ocl::enums::DeviceInfo::Type) {
                        let device_type_str = device_type.to_string();
                        if device_type_str.contains("GPU") {
                            return true;
                        }
                    }
                }
                false
            },
            Err(_) => false,
        }
    }
    
    /// Get the maximum work group size for the device
    pub fn get_max_work_group_size(&self) -> Result<usize> {
        let max_wg_size = self.device.info(ocl::enums::DeviceInfo::MaxWorkGroupSize)?;
        // Convert string to usize safely
        let size_str = max_wg_size.to_string();
        let size = size_str.parse::<usize>().unwrap_or(64);
        Ok(size)
    }
    
    /// Get the recommended work group size based on device type
    pub fn get_recommended_work_group_size(&self) -> Result<usize> {
        let max_wg_size = self.device.info(ocl::enums::DeviceInfo::MaxWorkGroupSize)?;
        let max_size = max_wg_size.to_string().parse::<usize>().unwrap_or(64);
        
        // Choose work group size based on device type
        let vendor = self.device.info(ocl::enums::DeviceInfo::Vendor)?;
        let vendor_str = vendor.to_string();
        
        let device_type = self.device.info(ocl::enums::DeviceInfo::Type)?;
        let device_type_str = device_type.to_string();
        
        if vendor_str.contains("NVIDIA") {
            // NVIDIA GPUs typically work well with 256
            if device_type_str.contains("GPU") {
                return Ok(256.min(max_size));
            }
        }
        
        if vendor_str.contains("AMD") {
            // AMD GPUs
            if device_type_str.contains("GPU") {
                return Ok(256.min(max_size));
            }
        }
        
        if vendor_str.contains("Intel") {
            // Intel GPUs
            if device_type_str.contains("GPU") {
                return Ok(128.min(max_size));
            }
            // Intel CPUs 
            else {
                return Ok(8.min(max_size));
            }
        }
        
        // Generic device - use a conservative value
        Ok(64.min(max_size))
    }
} 