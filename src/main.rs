//! # trajin CLI Tool
//!
//! This simple command-line tool demonstrates the capabilities of the trajin library.
//! It opens a DCD file, reads its contents, and displays information about the trajectory.

// Import the necessary components from our library
use std::path::Path;
use std::time::Instant;
use trajin::DcdReaderBuilder; // For performance measurement

/// Main entry point for the CLI application
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Get the first command line argument as the input file path
    // If no argument is provided, display a usage message
    let path = std::env::args().nth(1).expect("Usage: trajin <file.dcd>");

    // Start a timer to measure performance
    let start_time = Instant::now();

    // Create a DcdReader using the builder pattern for more flexibility
    let mut reader = if Path::new(&path).exists() {
        // Use the builder pattern which allows fine-grained control
        DcdReaderBuilder::new(&path)
            // Uncomment to force endianness if needed
            // .with_endianness(Endianness::Little)
            .build()?
    } else {
        // Handle the case where the file doesn't exist
        eprintln!("File not found: {}", path);
        return Ok(());
    };

    // Display basic information about the DCD file
    println!("DCD File Information:");
    println!("  Number of frames: {}", reader.num_frames());
    println!("  Number of atoms: {}", reader.num_atoms());
    println!("  Starting timestep: {}", reader.start_timestep());
    println!("  Timestep frequency: {}", reader.timestep_freq());
    println!("  Last timestep: {}", reader.end_timestep());
    println!("  Timestep size: {}", reader.timestep_size());
    println!("  Endianness: {:?}", reader.endianness());

    // Display title information if present
    if !reader.titles().is_empty() {
        println!("  Titles:");
        for title in reader.titles() {
            println!("    {}", title);
        }
    }

    // Demonstrate how to use the Iterator implementation to process all frames
    let mut frame_count = 0;

    // Iterator-based processing - for every frame in the trajectory
    for frame_result in &mut reader {
        // Handle any potential errors
        let frame = frame_result?;
        frame_count += 1;

        // Uncomment to see detailed information about each frame
        println!(
            "Frame {} (timestep {}): {} atoms, box dimensions: {:?}",
            frame_count,
            frame.timestep,
            frame.positions.len(),
            frame.box_dimensions
        );
    }

    // Demonstrate seeking to a specific frame
    if reader.num_frames() > 1 {
        println!("\nSeeking to first frame:");

        // The seek_frame method allows random access to frames by index
        reader.seek_frame(0)?;

        // Read the frame we sought to
        if let Some(frame) = reader.read_next()? {
            println!("  First frame timestep: {}", frame.timestep);

            // Show the position of the first atom
            if !frame.positions.is_empty() {
                println!("  First atom position: {:?}", frame.positions[0]);
            }
        }
    }

    // Calculate and display performance metrics
    let elapsed = start_time.elapsed();
    println!("\nRead {} frames in {:.2?}", frame_count, elapsed);

    // Calculate frames per second for benchmarking
    let fps = frame_count as f64 / elapsed.as_secs_f64();
    println!("Performance: {:.1} frames/sec", fps);

    Ok(())
}
