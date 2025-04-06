use resmol::DcdReader;
use std::error::Error;
use std::time::Instant;

fn main() -> Result<(), Box<dyn Error>> {
    let path = std::env::args()
        .nth(1)
        .expect("Usage: cargo run --example performance <file.dcd>");

    let start = Instant::now();
    let mut reader = DcdReader::open(&path)?;

    println!(
        "DCD: {} frames, {} atoms",
        reader.num_frames(),
        reader.num_atoms()
    );

    // Process all frames
    let mut total_frames = 0;
    let mut total_coords = 0;

    while let Some(frame) = reader.read_next()? {
        total_frames += 1;
        total_coords += frame.positions.len();
    }

    let elapsed = start.elapsed();
    println!(
        "Read {} frames ({} coordinates) in {:.2?}",
        total_frames, total_coords, elapsed
    );

    let fps = total_frames as f64 / elapsed.as_secs_f64();
    println!("Performance: {:.1} frames/sec", fps);

    Ok(())
}
