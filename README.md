# resmol

A Rust library for working with molecular dynamics simulation files.

## Features

- High-performance DCD trajectory file reader
- Memory-mapped file access for efficient frame seeking
- Support for both little and big endian files
- Iterator-based API
- Zero dependencies except for glam Vec3 implementation

## Usage

```rust
use resmol::DcdReader;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Open a DCD file
    let mut reader = DcdReader::open("trajectory.dcd")?;
    
    println!("DCD contains {} frames with {} atoms", 
             reader.num_frames(), reader.num_atoms());
             
    // Iterate through all frames
    for frame_result in reader {
        let frame = frame_result?;
        
        // Process frame data...
        let positions = frame.positions;
        let box_dim = frame.box_dimensions;
        
        // Find atoms in regions of interest
        let atoms_in_region = positions.iter()
            .filter(|pos| pos.z > 10.0 && pos.z < 20.0)
            .count();
            
        println!("Frame {}: {} atoms in region", 
                 frame.timestep, atoms_in_region);
    }
    
    Ok(())
}