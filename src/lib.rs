use glam::Vec3;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::Path;

/// A reader for CHARMM/NAMD DCD trajectory files
pub struct DcdReader {
    // File reader
    reader: std::io::BufReader<std::fs::File>,
    // File size
    filesize: u64,
    // Size of each frame in bytes
    framesize: u64,
    // Titles from the DCD file
    titles: Vec<String>,
    // Number of atoms
    natoms: i32,
    // Number of frames
    nframes: i32,
    // First timestep
    istart: i32,
    // Timestep frequency
    nevery: i32,
    // Last timestep
    iend: i32,
    // Timestep size
    timestep: f32,
    // Current position in file
    current_frame: i32,
}

/// Represents a single frame from a DCD file
pub struct Frame {
    /// Atom positions (x, y, z)
    pub positions: Vec<glam::Vec3>,
    /// Box dimensions [A, B, C, alpha, beta, gamma]
    pub box_dimensions: [f64; 6],
    /// Timestep of this frame
    pub timestep: i32,
}

/// Possible errors when reading DCD files
#[derive(Debug)]
pub enum DcdError {
    Io(std::io::Error),
    Format(String),
}

impl From<std::io::Error> for DcdError {
    fn from(err: std::io::Error) -> Self {
        DcdError::Io(err)
    }
}

impl std::fmt::Display for DcdError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DcdError::Io(err) => write!(f, "I/O error: {}", err),
            DcdError::Format(msg) => write!(f, "DCD format error: {}", msg),
        }
    }
}

impl std::error::Error for DcdError {}

impl DcdReader {
    /// Opens a DCD file for reading
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, DcdError> {
        let file = File::open(path)?;
        let filesize = file.metadata()?.len();
        let mut reader = BufReader::new(file);

        // Check for magic number and string
        let magic_number = read_i32(&mut reader)?;
        let mut magic_string = [0u8; 4];
        reader.read_exact(&mut magic_string)?;

        if magic_number != 84 || &magic_string != b"CORD" {
            // Try swapping endianness (not implemented here - would need custom byte swapping)
            return Err(DcdError::Format(
                "This DCD file format is not supported, or the file header is corrupt.".into(),
            ));
        }

        // Check if the file identifies as CHARMM (LAMMPS pretends to be CHARMM v. 24)
        reader.seek(SeekFrom::Start(84))?;
        let charmm_version = read_i32(&mut reader)?;
        if charmm_version == 0 {
            return Err(DcdError::Format(
                "DCD file indicates it is not CHARMM. Only CHARMM-style DCD files are supported."
                    .into(),
            ));
        }

        // Check for extra unit cell block
        reader.seek(SeekFrom::Start(48))?;
        let has_extra_block = read_i32(&mut reader)?;
        if has_extra_block != 1 {
            return Err(DcdError::Format("DCD file indicates it does not have unit cell information. Only DCD files with unit cell information are supported.".into()));
        }

        // Check for 4D support (we don't support it)
        let four_dimensions = read_i32(&mut reader)?;
        if four_dimensions == 1 {
            return Err(DcdError::Format("DCD file indicates it has four dimensions. Only DCD files with three dimensions are supported.".into()));
        }

        let mut dcd_reader = DcdReader {
            reader,
            filesize,
            framesize: 0,
            titles: Vec::new(),
            natoms: 0,
            nframes: 0,
            istart: 0,
            nevery: 0,
            iend: 0,
            timestep: 0.0,
            current_frame: 0,
        };

        dcd_reader.read_header()?;
        Ok(dcd_reader)
    }

    /// Reads the header of the DCD file
    fn read_header(&mut self) -> Result<(), DcdError> {
        // Get nframes, istart, nevery, iend
        self.reader.seek(SeekFrom::Start(8))?;
        self.nframes = read_i32(&mut self.reader)?;
        self.istart = read_i32(&mut self.reader)?;
        self.nevery = read_i32(&mut self.reader)?;
        self.iend = read_i32(&mut self.reader)?;

        // Get timestep
        self.reader.seek(SeekFrom::Start(44))?;
        self.timestep = read_f32(&mut self.reader)?;

        // Get titles
        self.reader.seek(SeekFrom::Start(96))?;
        let ntitle = read_i32(&mut self.reader)?;

        self.titles = Vec::with_capacity(ntitle as usize);
        for _ in 0..ntitle {
            let mut title_bytes = [0u8; 80];
            self.reader.read_exact(&mut title_bytes)?;

            // Find first null character or use the whole string
            let null_pos = title_bytes.iter().position(|&b| b == 0).unwrap_or(80);
            let title = String::from_utf8_lossy(&title_bytes[0..null_pos]).to_string();
            self.titles.push(title);
        }

        // Get natoms
        // Skip 8 bytes (2 int32s)
        read_i32(&mut self.reader)?;
        read_i32(&mut self.reader)?;
        self.natoms = read_i32(&mut self.reader)?;
        read_i32(&mut self.reader)?; // Skip 4 bytes (1 int32)

        // Calculate current position
        let current_pos = self.reader.stream_position()? as i64;

        // Calculate frame size: natoms*3 coordinates (4 bytes each) = natoms*12
        // plus 6 box dimensions (8 bytes each) = 48
        // Plus 32 bytes of file information in each frame
        self.framesize = (self.natoms as u64) * 12 + 80;

        // Recalculate nframes based on file size
        let nframes2 = (self.filesize - (current_pos as u64)) / self.framesize;
        if (nframes2 as i32) != self.nframes {
            eprintln!(
                "WARNING: Header indicates {} frames, but file size indicates {}.",
                self.nframes, nframes2
            );
            self.nframes = nframes2 as i32;
        }

        Ok(())
    }

    /// Reads the next frame from the DCD file
    pub fn read_next(&mut self) -> Result<Option<Frame>, DcdError> {
        if self.current_frame >= self.nframes {
            return Ok(None);
        }

        // Read box dimensions (6 double precision values)
        let dummy = read_i32(&mut self.reader)?;
        if dummy != 48 {
            return Err(DcdError::Format("Problem reading in DCD snapshot.".into()));
        }

        let mut box_dimensions = [0.0; 6];
        box_dimensions[0] = read_f64(&mut self.reader)?; // A
        box_dimensions[5] = read_f64(&mut self.reader)?; // gamma
        box_dimensions[1] = read_f64(&mut self.reader)?; // B
        box_dimensions[4] = read_f64(&mut self.reader)?; // beta
        box_dimensions[3] = read_f64(&mut self.reader)?; // alpha
        box_dimensions[2] = read_f64(&mut self.reader)?; // C

        if box_dimensions[0] < 0.0 || box_dimensions[1] < 0.0 || box_dimensions[2] < 0.0 {
            return Err(DcdError::Format(
                "Problem reading in DCD snapshot box dimensions.".into(),
            ));
        }

        // Read coordinates
        let nbytes = (self.natoms * 4) as i32;
        let mut positions = Vec::with_capacity(self.natoms as usize);

        // X coordinates
        let dummy1 = read_i32(&mut self.reader)?;
        let dummy2 = read_i32(&mut self.reader)?;

        if dummy1 != 48 {
            return Err(DcdError::Format(
                "Problem reading in DCD snapshot coordinates.".into(),
            ));
        }

        if dummy2 != nbytes {
            return Err(DcdError::Format(
                "Number of bytes in DCD snapshot is incorrect for size of xyz array passed.".into(),
            ));
        }

        let mut x_coords = vec![0.0f32; self.natoms as usize];
        for i in 0..self.natoms as usize {
            x_coords[i] = read_f32(&mut self.reader)?;
        }

        // Y coordinates
        let dummy3 = read_i32(&mut self.reader)?;
        let dummy4 = read_i32(&mut self.reader)?;

        if dummy3 != nbytes || dummy4 != nbytes {
            return Err(DcdError::Format(
                "Number of bytes in DCD snapshot is incorrect for size of xyz array passed.".into(),
            ));
        }

        let mut y_coords = vec![0.0f32; self.natoms as usize];
        for i in 0..self.natoms as usize {
            y_coords[i] = read_f32(&mut self.reader)?;
        }

        // Z coordinates
        let dummy5 = read_i32(&mut self.reader)?;
        let dummy6 = read_i32(&mut self.reader)?;

        if dummy5 != nbytes || dummy6 != nbytes {
            return Err(DcdError::Format(
                "Number of bytes in DCD snapshot is incorrect for size of xyz array passed.".into(),
            ));
        }

        let mut z_coords = vec![0.0f32; self.natoms as usize];
        for i in 0..self.natoms as usize {
            z_coords[i] = read_f32(&mut self.reader)?;
        }

        // Final read to verify
        let final_dummy = read_i32(&mut self.reader)?;
        if final_dummy != nbytes {
            return Err(DcdError::Format("Problem reading in DCD snapshot.".into()));
        }

        // Construct positions vector
        for i in 0..self.natoms as usize {
            positions.push(Vec3::new(x_coords[i], y_coords[i], z_coords[i]));
        }

        self.current_frame += 1;

        Ok(Some(Frame {
            positions,
            box_dimensions,
            timestep: self.istart + (self.current_frame - 1) * self.nevery,
        }))
    }

    /// Skips the next frame(s) instead of reading them
    pub fn skip_next(&mut self, n: Option<i32>) -> Result<(), DcdError> {
        let n = n.unwrap_or(1);
        if self.current_frame + n > self.nframes {
            return Err(DcdError::Format(
                "Attempted to skip beyond the end of the file".into(),
            ));
        }

        let pos = self.reader.stream_position()?;
        let new_pos = pos + (self.framesize * n as u64);
        self.reader.seek(SeekFrom::Start(new_pos))?;

        self.current_frame += n;
        Ok(())
    }

    /// Returns the total number of frames in the file
    pub fn num_frames(&self) -> i32 {
        self.nframes
    }

    /// Returns the number of atoms in each frame
    pub fn num_atoms(&self) -> i32 {
        self.natoms
    }

    /// Returns the titles from the DCD file
    pub fn titles(&self) -> &[String] {
        &self.titles
    }

    /// Returns the starting timestep of the trajectory
    pub fn start_timestep(&self) -> i32 {
        self.istart
    }

    /// Returns how often (in timesteps) frames were saved
    pub fn timestep_freq(&self) -> i32 {
        self.nevery
    }

    /// Returns the last timestep of the trajectory
    pub fn end_timestep(&self) -> i32 {
        self.iend
    }

    /// Returns the timestep size in simulation units
    pub fn timestep_size(&self) -> f32 {
        self.timestep
    }
}

// Helper functions for reading different data types
fn read_i32<R: Read>(reader: &mut R) -> Result<i32, DcdError> {
    let mut buffer = [0u8; 4];
    reader.read_exact(&mut buffer)?;
    Ok(i32::from_le_bytes(buffer))
}

fn read_f32<R: Read>(reader: &mut R) -> Result<f32, DcdError> {
    let mut buffer = [0u8; 4];
    reader.read_exact(&mut buffer)?;
    Ok(f32::from_le_bytes(buffer))
}

fn read_f64<R: Read>(reader: &mut R) -> Result<f64, DcdError> {
    let mut buffer = [0u8; 8];
    reader.read_exact(&mut buffer)?;
    Ok(f64::from_le_bytes(buffer))
}
