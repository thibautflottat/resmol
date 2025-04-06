use resmol::DcdReader;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let path = std::env::args()
        .nth(1)
        .expect("Usage: cargo run --example basic <file.dcd>");
    let mut reader = DcdReader::open(&path)?;

    println!(
        "DCD contains {} frames with {} atoms each",
        reader.num_frames(),
        reader.num_atoms()
    );

    // Count atoms in specific regions of the simulation box
    let mut frame_count = 0;
    let mut atoms_in_top_half = 0;

    while let Some(frame) = reader.read_next()? {
        frame_count += 1;

        // Count atoms in top half of z-dimension
        let top_half = frame.positions.iter().filter(|pos| pos.z > 0.0).count();

        atoms_in_top_half += top_half;
    }

    let avg_top_half = atoms_in_top_half as f64 / frame_count as f64;
    println!("Average atoms in top half: {:.1}", avg_top_half);

    Ok(())
}
