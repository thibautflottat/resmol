use std::error::Error;
use std::path::Path;
use trajin::DcdReader;

#[test]
fn test_builder_pattern() -> Result<(), Box<dyn Error>> {
    // This test would need a real DCD file
    // Skip if file doesn't exist
    let path = "tests/data/test.dcd";
    if !Path::new(path).exists() {
        return Ok(());
    }

    let reader = DcdReader::open(&path)?;
    assert_eq!(reader.num_frames(), 100);
    Ok(())
}
