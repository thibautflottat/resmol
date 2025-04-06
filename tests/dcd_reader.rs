use std::error::Error;
use std::path::Path;
use trajin::{DcdReader, DcdReaderBuilder, Endianness};

#[test]
fn test_builder_pattern() -> Result<(), Box<dyn Error>> {
    // This test would need a real DCD file
    // Skip if file doesn't exist
    let path = "tests/data/test.dcd";
    if !Path::new(path).exists() {
        return Ok(());
    }

    let reader = DcdReaderBuilder::new(path)
        .with_endianness(Endianness::Little)
        .skip_validation(true)
        .build()?;

    assert_eq!(reader.endianness(), Endianness::Little);
    Ok(())
}

#[test]
fn test_iterator_interface() -> Result<(), Box<dyn Error>> {
    // This test would need a real DCD file
    // Skip if file doesn't exist
    let path = "tests/data/test.dcd";
    if !Path::new(path).exists() {
        return Ok(());
    }

    let reader = DcdReader::open(path)?;
    let expected_frames = reader.num_frames();

    let frame_count = reader.count();
    assert_eq!(frame_count, expected_frames as usize);

    Ok(())
}

#[test]
fn test_seek_frame() -> Result<(), Box<dyn Error>> {
    // This test would need a real DCD file
    // Skip if file doesn't exist
    let path = "tests/data/traj.dcd";
    if !Path::new(path).exists() {
        return Ok(());
    }

    let mut reader = DcdReader::open(path)?;
    if reader.num_frames() < 2 {
        return Ok(());
    }

    // Get first frame timestep
    let first_frame = reader.read_next()?.unwrap();
    let first_timestep = first_frame.timestep;

    // Seek to last frame
    reader.seek_frame(reader.num_frames() - 1)?;
    let last_frame = reader.read_next()?.unwrap();

    // Seek back to first
    reader.seek_frame(0)?;
    let new_first_frame = reader.read_next()?.unwrap();

    assert_eq!(first_timestep, new_first_frame.timestep);
    assert_ne!(first_timestep, last_frame.timestep);

    Ok(())
}
