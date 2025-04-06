use std::time::Instant;
use trajin::DcdReader;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start = Instant::now();
    let path = std::env::args().nth(1).expect("Usage: resmol <file.dcd>");
    let mut reader = DcdReader::open(&path)?;

    println!("DCD File Information:");
    println!("  Number of frames: {}", reader.num_frames());
    println!("  Number of atoms: {}", reader.num_atoms());
    println!("  Starting timestep: {}", reader.start_timestep());
    println!("  Timestep frequency: {}", reader.timestep_freq());
    println!("  Last timestep: {}", reader.end_timestep());
    println!("  Timestep size: {}", reader.timestep_size());

    if !reader.titles().is_empty() {
        println!("  Titles:");
        for title in reader.titles() {
            println!("    {}", title);
        }
    }

    let mut frame_count = 0;
    while let Some(frame) = reader.read_next()? {
        frame_count += 1;
        // println!(
        //     "Frame {} (timestep {}): {} atoms, box dimensions: {:?}",
        //     frame_count,
        //     frame.timestep,
        //     frame.positions.len(),
        //     frame.box_dimensions
        // );

        // Optional: You can process the frame data here
    }
    let elapsed = start.elapsed();
    println!("Processed {} frames in {:.2?}", frame_count, elapsed);
    println!("Read {} frames total.", frame_count);
    Ok(())
}
