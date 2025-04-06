//! # DCD Trajectory File Reader
//!
//! This module provides functionality to read DCD files produced by CHARMM, NAMD,
//! and other molecular dynamics simulation packages. DCD is a binary format that
//! stores atomic coordinates for molecular dynamics trajectories.
//!
//! ## Performance Optimizations
//!
//! This implementation includes several optimizations:
//! - Memory mapping for efficient file access
//! - Pre-allocated buffers to reduce memory allocations
//! - Zero-copy reading for little endian files (most common case)
//! - Chunked processing for big endian files to improve cache efficiency
//! - Efficient seeking directly to frames
//!
//! ## Example
//!
//! ```no_run
//! use resmol::DcdReader;
//!
//! let mut reader = DcdReader::open("trajectory.dcd").unwrap();
//! for frame_result in reader {
//!     let frame = frame_result.unwrap();
//!     // Process frame...
//! }
//! ```

// Import necessary dependencies
use glam::Vec3; // For 3D vector representation of atom positions
use memmap2::{Mmap, MmapOptions}; // For memory-mapped file access
use std::convert::TryInto; // For converting byte slices to arrays
use std::fs::File; // For file operations
use std::path::{Path, PathBuf}; // For file path handling
use thiserror::Error; // For better error handling

/// Endianness of the DCD file
///
/// DCD files can be either little endian (most common on x86 systems)
/// or big endian (historically from older UNIX workstations).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Endianness {
    /// Little endian (most common)
    Little,
    /// Big endian
    Big,
}

/// A reader for CHARMM/NAMD DCD trajectory files
///
/// This struct handles reading and parsing DCD files, providing methods
/// to access header information and iterate through frames.
pub struct DcdReader {
    /// Memory-mapped file for efficient random access
    ///
    /// Using memory mapping instead of standard file I/O provides:
    /// - Better performance for random access (seeking)
    /// - The OS handles caching of the file data
    /// - Allows zero-copy access to file contents
    mmap: Mmap,

    /// File size in bytes
    filesize: u64,

    /// Size of each frame in bytes
    /// Used for seeking and validation
    framesize: u64,

    /// Titles from the DCD file header
    /// These usually contain information about the simulation
    titles: Vec<String>,

    /// Number of atoms in each frame
    natoms: i32,

    /// Number of frames in the file
    nframes: i32,

    /// First timestep of the trajectory
    istart: i32,

    /// Timestep frequency (how often frames were saved)
    nevery: i32,

    /// Last timestep of the trajectory
    iend: i32,

    /// Timestep size in simulation units
    timestep: f32,

    /// Current position in file (in bytes)
    current_pos: usize,

    /// Current frame number (0-based index)
    current_frame: i32,

    /// Endianness of the file
    endianness: Endianness,

    /// Pre-allocated buffer for positions (improved performance)
    /// This avoids reallocating memory for each frame
    positions_buffer: Vec<Vec3>,

    /// Pre-allocated buffers for coordinate processing
    /// These store the raw coordinate data before constructing Vec3 objects
    x_buffer: Vec<f32>,
    y_buffer: Vec<f32>,
    z_buffer: Vec<f32>,
}

/// Represents a single frame from a DCD file
///
/// Contains the atomic positions and unit cell information for a snapshot
/// from a molecular dynamics trajectory.
pub struct Frame {
    /// Atom positions (x, y, z) for each atom
    pub positions: Vec<Vec3>,

    /// Box dimensions [A, B, C, alpha, beta, gamma]
    /// A, B, C are lengths of the cell vectors
    /// alpha, beta, gamma are angles between the vectors (in degrees)
    pub box_dimensions: [f64; 6],

    /// Timestep of this frame
    /// This is calculated as istart + (frame_index * nevery)
    pub timestep: i32,
}

/// Possible errors when reading DCD files
///
/// This enum represents the different types of errors that can occur
/// when reading and parsing DCD files.
#[derive(Error, Debug)]
pub enum DcdError {
    /// I/O errors from the underlying file system
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// Format errors from invalid DCD file structure
    #[error("DCD format error: {0}")]
    Format(String),

    /// Endianness errors when byte order cannot be determined
    #[error("Endianness error: {0}")]
    Endianness(String),

    /// Validation errors when file doesn't meet expected criteria
    #[error("Validation error: {0}")]
    Validation(String),
}

/// Builder for DcdReader to allow flexible construction
///
/// This implements the Builder pattern to provide a more flexible
/// and user-friendly way to construct a DcdReader with various options.
pub struct DcdReaderBuilder {
    /// Path to the DCD file
    path: PathBuf,

    /// Force a specific endianness (optional)
    /// If not specified, it will be auto-detected
    endianness: Option<Endianness>,

    /// Skip validation checks for file format
    /// Useful for non-standard DCD files that don't strictly follow the format
    skip_validation: bool,

    /// Pre-allocate buffers for better performance
    /// Allows customizing the initial buffer size
    buffer_size: Option<usize>,
}

impl DcdReaderBuilder {
    /// Create a new builder
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the DCD file
    ///
    /// # Returns
    ///
    /// A new DcdReaderBuilder instance
    pub fn new<P: AsRef<Path>>(path: P) -> Self {
        Self {
            path: path.as_ref().to_path_buf(),
            endianness: None,
            skip_validation: false,
            buffer_size: None,
        }
    }

    /// Force a specific endianness
    ///
    /// By default, the endianness is auto-detected, but you can
    /// force a specific endianness if needed.
    ///
    /// # Arguments
    ///
    /// * `endianness` - The endianness to use (Little or Big)
    ///
    /// # Returns
    ///
    /// The builder with the endianness set
    pub fn with_endianness(mut self, endianness: Endianness) -> Self {
        self.endianness = Some(endianness);
        self
    }

    /// Skip validation checks
    ///
    /// This can be useful for non-standard DCD files that don't
    /// strictly follow the format specification.
    ///
    /// # Arguments
    ///
    /// * `skip` - Whether to skip validation (true) or not (false)
    ///
    /// # Returns
    ///
    /// The builder with the validation setting updated
    pub fn skip_validation(mut self, skip: bool) -> Self {
        self.skip_validation = skip;
        self
    }

    /// Set the buffer size for pre-allocation (for performance)
    ///
    /// # Arguments
    ///
    /// * `size` - The buffer size to pre-allocate
    ///
    /// # Returns
    ///
    /// The builder with the buffer size set
    pub fn with_buffer_size(mut self, size: usize) -> Self {
        self.buffer_size = Some(size);
        self
    }

    /// Build the reader
    ///
    /// This method builds the DcdReader with the specified options.
    ///
    /// # Returns
    ///
    /// A Result containing either the DcdReader or an error
    pub fn build(self) -> Result<DcdReader, DcdError> {
        // Open the file
        let file = File::open(&self.path)?;
        let filesize = file.metadata()?.len();

        // Memory map the file for efficient access
        // Using unsafe because memory mapping is inherently unsafe
        // (we're directly accessing file bytes as memory)
        let mmap = unsafe { MmapOptions::new().map(&file)? };

        // Determine endianness - either use the forced value or auto-detect
        let mut pos = 0;
        let endianness = if let Some(e) = self.endianness {
            // Use the forced endianness
            e
        } else {
            // Auto-detect by checking the magic number
            // First try little endian (most common)
            let magic_number = read_i32_at(&mmap, pos, Endianness::Little)?;

            if magic_number == 84 {
                // 84 is the expected value in the correct endianness
                Endianness::Little
            } else {
                // If little endian didn't work, try big endian
                let magic_number = read_i32_at(&mmap, pos, Endianness::Big)?;
                if magic_number == 84 {
                    Endianness::Big
                } else {
                    // If neither endianness works, the file is invalid
                    return Err(DcdError::Endianness(
                        "Could not determine file endianness from magic number".to_string(),
                    ));
                }
            }
        };

        pos += 4;

        // Check magic string "CORD" which identifies a DCD file
        let magic_string = &mmap[pos..pos + 4];
        if magic_string != b"CORD" {
            return Err(DcdError::Format(
                "This DCD file format is not supported, or the file header is corrupt.".to_string(),
            ));
        }

        // Perform validation checks if not skipped
        if !self.skip_validation {
            // Check if the file identifies as CHARMM
            // (LAMMPS pretends to be CHARMM v. 24)
            let charmm_version = read_i32_at(&mmap, 84, endianness)?;
            if charmm_version == 0 {
                return Err(DcdError::Validation(
                    "DCD file indicates it is not CHARMM. Only CHARMM-style DCD files are supported.".to_string()
                ));
            }

            // Check for unit cell information
            // We require this for proper trajectory analysis
            let has_extra_block = read_i32_at(&mmap, 48, endianness)?;
            if has_extra_block != 1 {
                return Err(DcdError::Validation(
                    "DCD file indicates it does not have unit cell information. Only DCD files with unit cell information are supported.".to_string()
                ));
            }

            // Check for 4D support (we don't support 4D trajectories)
            let four_dimensions = read_i32_at(&mmap, 52, endianness)?;
            if four_dimensions == 1 {
                return Err(DcdError::Validation(
                    "DCD file indicates it has four dimensions. Only DCD files with three dimensions are supported.".to_string()
                ));
            }
        }

        // Create the reader with default values
        // The actual values will be filled in by read_header()
        let mut dcd_reader = DcdReader {
            mmap,
            filesize,
            framesize: 0,
            titles: Vec::new(),
            natoms: 0,
            nframes: 0,
            istart: 0,
            nevery: 0,
            iend: 0,
            timestep: 0.0,
            current_pos: 0,
            current_frame: 0,
            endianness,
            positions_buffer: Vec::new(),
            x_buffer: Vec::new(),
            y_buffer: Vec::new(),
            z_buffer: Vec::new(),
        };

        // Read the header to fill in the values
        dcd_reader.read_header()?;

        // Pre-allocate buffers based on natoms for better performance
        // This avoids reallocating memory for each frame
        let buffer_size = self.buffer_size.unwrap_or(dcd_reader.natoms as usize);
        dcd_reader.positions_buffer = Vec::with_capacity(buffer_size);
        dcd_reader.x_buffer = Vec::with_capacity(buffer_size);
        dcd_reader.y_buffer = Vec::with_capacity(buffer_size);
        dcd_reader.z_buffer = Vec::with_capacity(buffer_size);

        Ok(dcd_reader)
    }
}

impl DcdReader {
    /// Opens a DCD file for reading
    ///
    /// This is a convenience method that uses the builder pattern internally.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the DCD file
    ///
    /// # Returns
    ///
    /// A Result containing either the DcdReader or an error
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, DcdError> {
        DcdReaderBuilder::new(path).build()
    }

    /// Reads the header of the DCD file
    ///
    /// This parses the DCD header and fills in the reader's fields.
    ///
    /// # Returns
    ///
    /// A Result indicating success or an error
    fn read_header(&mut self) -> Result<(), DcdError> {
        // Get nframes, istart, nevery, iend from the header
        // These values are at fixed offsets in the header
        self.nframes = read_i32_at(&self.mmap, 8, self.endianness)?;
        self.istart = read_i32_at(&self.mmap, 12, self.endianness)?;
        self.nevery = read_i32_at(&self.mmap, 16, self.endianness)?;
        self.iend = read_i32_at(&self.mmap, 20, self.endianness)?;

        // Get timestep
        self.timestep = read_f32_at(&self.mmap, 44, self.endianness)?;

        // Get titles
        // First get the number of title lines
        let ntitle = read_i32_at(&self.mmap, 96, self.endianness)?;
        let mut pos = 100;

        // Read each title line (80 characters each)
        self.titles = Vec::with_capacity(ntitle as usize);
        for _ in 0..ntitle {
            let title_end = pos + 80;
            if title_end > self.mmap.len() {
                return Err(DcdError::Format(
                    "Unexpected end of file reading titles".to_string(),
                ));
            }

            let title_bytes = &self.mmap[pos..title_end];

            // Find first null character or use the whole string
            // Titles are null-terminated C strings
            let null_pos = title_bytes.iter().position(|&b| b == 0).unwrap_or(80);
            let title = String::from_utf8_lossy(&title_bytes[0..null_pos]).to_string();
            self.titles.push(title);

            pos += 80;
        }

        // Get number of atoms
        // Skip 8 bytes (2 int32s) to get to the natoms field
        pos += 8;
        self.natoms = read_i32_at(&self.mmap, pos, self.endianness)?;
        pos += 8; // Skip 8 bytes (2 int32s) more

        // Calculate frame size:
        // - natoms*3 coordinates (4 bytes each) = natoms*12
        // - plus 6 box dimensions (8 bytes each) = 48
        // - plus 32 bytes of file information in each frame
        self.framesize = (self.natoms as u64) * 12 + 80;

        // Recalculate nframes based on file size
        // This is a sanity check - sometimes the header nframes doesn't match the actual file
        let remaining_size = self.filesize - (pos as u64);
        let nframes2 = remaining_size / self.framesize;
        if (nframes2 as i32) != self.nframes {
            eprintln!(
                "WARNING: Header indicates {} frames, but file size indicates {}.",
                self.nframes, nframes2
            );
            // Use the calculated value as it's based on actual file size
            self.nframes = nframes2 as i32;
        }

        // Set current position to the end of the header
        self.current_pos = pos;
        Ok(())
    }

    /// Reads the next frame from the DCD file
    ///
    /// # Returns
    ///
    /// A Result containing either:
    /// - Some(Frame) if a frame was successfully read
    /// - None if we've reached the end of the file
    /// - An error if something went wrong
    pub fn read_next(&mut self) -> Result<Option<Frame>, DcdError> {
        // Check if we've reached the end of the file
        if self.current_frame >= self.nframes {
            return Ok(None);
        }

        // Make sure we have enough data left to read a frame
        if self.current_pos + (self.framesize as usize) > self.mmap.len() {
            return Err(DcdError::Format(
                "Unexpected end of file reading frame".to_string(),
            ));
        }

        // Read box dimensions (6 double precision values)
        let mut pos = self.current_pos;

        // Check for the box dimension marker (should be 48)
        let dummy = read_i32_at(&self.mmap, pos, self.endianness)?;
        if dummy != 48 {
            return Err(DcdError::Format(
                "Problem reading in DCD snapshot.".to_string(),
            ));
        }
        pos += 4;

        // Read the 6 box dimension values
        // A, B, C are lengths, alpha, beta, gamma are angles
        let mut box_dimensions = [0.0; 6];
        box_dimensions[0] = read_f64_at(&self.mmap, pos, self.endianness)?; // A
        pos += 8;
        box_dimensions[5] = read_f64_at(&self.mmap, pos, self.endianness)?; // gamma
        pos += 8;
        box_dimensions[1] = read_f64_at(&self.mmap, pos, self.endianness)?; // B
        pos += 8;
        box_dimensions[4] = read_f64_at(&self.mmap, pos, self.endianness)?; // beta
        pos += 8;
        box_dimensions[3] = read_f64_at(&self.mmap, pos, self.endianness)?; // alpha
        pos += 8;
        box_dimensions[2] = read_f64_at(&self.mmap, pos, self.endianness)?; // C
        pos += 8;

        // Validate box dimensions - cell lengths should be positive
        if box_dimensions[0] < 0.0 || box_dimensions[1] < 0.0 || box_dimensions[2] < 0.0 {
            return Err(DcdError::Format(
                "Problem reading in DCD snapshot box dimensions.".to_string(),
            ));
        }

        // Read coordinate data
        let natoms = self.natoms as usize;

        // Use optimized, endianness-aware coordinate reading
        // Different paths for little vs big endian for better performance
        if self.endianness == Endianness::Little {
            // Little endian is the most common case
            pos = self.read_coordinates_little_endian(pos, natoms)?;
        } else {
            // Big endian requires byte swapping
            pos = self.read_coordinates_big_endian(pos, natoms)?;
        }

        // Generate the timestep and update position
        let timestep = self.istart + (self.current_frame) * self.nevery;
        self.current_frame += 1;
        self.current_pos = pos;

        // Create frame with the processed positions
        // Use std::mem::replace to avoid cloning the positions vector
        let positions = std::mem::replace(&mut self.positions_buffer, Vec::with_capacity(natoms));

        Ok(Some(Frame {
            positions,
            box_dimensions,
            timestep,
        }))
    }

    /// Optimized reading of coordinates for little endian files
    ///
    /// This is a specialized, highly optimized path for little endian files
    /// which are the most common case. It uses zero-copy casting where possible.
    ///
    /// # Arguments
    ///
    /// * `pos` - Current position in the file
    /// * `natoms` - Number of atoms to read
    ///
    /// # Returns
    ///
    /// A Result containing the new position after reading
    fn read_coordinates_little_endian(
        &mut self,
        mut pos: usize,
        natoms: usize,
    ) -> Result<usize, DcdError> {
        // Calculate the size of each coordinate block
        let nbytes = (natoms * 4) as i32;

        // Clear and reserve space in buffers
        // This avoids reallocations during reading
        self.positions_buffer.clear();
        self.positions_buffer.reserve(natoms);

        self.x_buffer.clear();
        self.x_buffer.reserve(natoms);

        self.y_buffer.clear();
        self.y_buffer.reserve(natoms);

        self.z_buffer.clear();
        self.z_buffer.reserve(natoms);

        // Read X coordinates header
        // First int32 should be 48, second should be nbytes
        let dummy1 = read_i32_at(&self.mmap, pos, self.endianness)?;
        pos += 4;
        let dummy2 = read_i32_at(&self.mmap, pos, self.endianness)?;
        pos += 4;

        // Validate the header values
        if dummy1 != 48 || dummy2 != nbytes {
            return Err(DcdError::Format(
                "Problem reading X coordinate block headers".to_string(),
            ));
        }

        // For little endian, we can directly cast the memory to f32 slices
        // This is a zero-copy operation for better performance
        let x_end = pos + natoms * 4;
        if x_end > self.mmap.len() {
            return Err(DcdError::Format(
                "Unexpected end of file reading X coordinates".to_string(),
            ));
        }

        // Use bytemuck for zero-copy memory casting
        // This interprets the raw bytes as f32 values without copying
        let x_slice = bytemuck::cast_slice::<u8, f32>(&self.mmap[pos..x_end]);
        self.x_buffer.extend_from_slice(x_slice);
        pos = x_end;

        // Read Y coordinates header
        let dummy3 = read_i32_at(&self.mmap, pos, self.endianness)?;
        pos += 4;
        let dummy4 = read_i32_at(&self.mmap, pos, self.endianness)?;
        pos += 4;

        // Validate the header values
        if dummy3 != nbytes || dummy4 != nbytes {
            return Err(DcdError::Format(
                "Problem reading Y coordinate block headers".to_string(),
            ));
        }

        // Read Y coordinates using the same zero-copy approach
        let y_end = pos + natoms * 4;
        if y_end > self.mmap.len() {
            return Err(DcdError::Format(
                "Unexpected end of file reading Y coordinates".to_string(),
            ));
        }

        let y_slice = bytemuck::cast_slice::<u8, f32>(&self.mmap[pos..y_end]);
        self.y_buffer.extend_from_slice(y_slice);
        pos = y_end;

        // Read Z coordinates header
        let dummy5 = read_i32_at(&self.mmap, pos, self.endianness)?;
        pos += 4;
        let dummy6 = read_i32_at(&self.mmap, pos, self.endianness)?;
        pos += 4;

        // Validate the header values
        if dummy5 != nbytes || dummy6 != nbytes {
            return Err(DcdError::Format(
                "Problem reading Z coordinate block headers".to_string(),
            ));
        }

        // Read Z coordinates using the same zero-copy approach
        let z_end = pos + natoms * 4;
        if z_end > self.mmap.len() {
            return Err(DcdError::Format(
                "Unexpected end of file reading Z coordinates".to_string(),
            ));
        }

        let z_slice = bytemuck::cast_slice::<u8, f32>(&self.mmap[pos..z_end]);
        self.z_buffer.extend_from_slice(z_slice);
        pos = z_end;

        // Final read to verify the block
        let final_dummy = read_i32_at(&self.mmap, pos, self.endianness)?;
        pos += 4;

        // Validate the final marker
        if final_dummy != nbytes {
            return Err(DcdError::Format(
                "Problem validating coordinate block".to_string(),
            ));
        }

        // Efficiently construct Vec3 positions from the x, y, z buffers
        for i in 0..natoms {
            self.positions_buffer.push(Vec3::new(
                self.x_buffer[i],
                self.y_buffer[i],
                self.z_buffer[i],
            ));
        }

        Ok(pos)
    }

    /// Optimized reading of coordinates for big endian files
    ///
    /// This is a specialized path for big endian files, which uses
    /// chunked processing for better cache utilization.
    ///
    /// # Arguments
    ///
    /// * `pos` - Current position in the file
    /// * `natoms` - Number of atoms to read
    ///
    /// # Returns
    ///
    /// A Result containing the new position after reading
    fn read_coordinates_big_endian(
        &mut self,
        mut pos: usize,
        natoms: usize,
    ) -> Result<usize, DcdError> {
        // Calculate the size of each coordinate block
        let nbytes = (natoms * 4) as i32;

        // Clear and reserve space in buffers
        // This avoids reallocations during reading
        self.positions_buffer.clear();
        self.positions_buffer.reserve(natoms);

        self.x_buffer.clear();
        self.x_buffer.reserve(natoms);

        self.y_buffer.clear();
        self.y_buffer.reserve(natoms);

        self.z_buffer.clear();
        self.z_buffer.reserve(natoms);

        // Read X coordinates header
        // First int32 should be 48, second should be nbytes
        let dummy1 = read_i32_at(&self.mmap, pos, self.endianness)?;
        pos += 4;
        let dummy2 = read_i32_at(&self.mmap, pos, self.endianness)?;
        pos += 4;

        // Validate the header values
        if dummy1 != 48 || dummy2 != nbytes {
            return Err(DcdError::Format(
                "Problem reading X coordinate block headers".to_string(),
            ));
        }

        // Read X coordinates in chunks for better cache utilization
        // Processing data in chunks that fit in CPU cache improves performance
        const CHUNK_SIZE: usize = 1024; // Adjust based on cache size
        let x_end = pos + natoms * 4;

        // Make sure we have enough data
        if x_end > self.mmap.len() {
            return Err(DcdError::Format(
                "Unexpected end of file reading X coordinates".to_string(),
            ));
        }

        // Process each chunk of X coordinates
        for chunk_start in (0..natoms).step_by(CHUNK_SIZE) {
            let chunk_end = std::cmp::min(chunk_start + CHUNK_SIZE, natoms);
            let chunk_size = chunk_end - chunk_start;

            // Read a chunk of coordinates
            for i in 0..chunk_size {
                let coord_pos = pos + (chunk_start + i) * 4;
                let value = read_f32_at(&self.mmap, coord_pos, self.endianness)?;
                self.x_buffer.push(value);
            }
        }
        pos = x_end;

        // Read Y coordinates header
        let dummy3 = read_i32_at(&self.mmap, pos, self.endianness)?;
        pos += 4;
        let dummy4 = read_i32_at(&self.mmap, pos, self.endianness)?;
        pos += 4;

        // Validate the header values
        if dummy3 != nbytes || dummy4 != nbytes {
            return Err(DcdError::Format(
                "Problem reading Y coordinate block headers".to_string(),
            ));
        }

        // Read Y coordinates in chunks
        let y_end = pos + natoms * 4;
        if y_end > self.mmap.len() {
            return Err(DcdError::Format(
                "Unexpected end of file reading Y coordinates".to_string(),
            ));
        }

        // Process each chunk of Y coordinates
        for chunk_start in (0..natoms).step_by(CHUNK_SIZE) {
            let chunk_end = std::cmp::min(chunk_start + CHUNK_SIZE, natoms);
            let chunk_size = chunk_end - chunk_start;

            // Read a chunk of coordinates
            for i in 0..chunk_size {
                let coord_pos = pos + (chunk_start + i) * 4;
                let value = read_f32_at(&self.mmap, coord_pos, self.endianness)?;
                self.y_buffer.push(value);
            }
        }
        pos = y_end;

        // Read Z coordinates header
        let dummy5 = read_i32_at(&self.mmap, pos, self.endianness)?;
        pos += 4;
        let dummy6 = read_i32_at(&self.mmap, pos, self.endianness)?;
        pos += 4;

        // Validate the header values
        if dummy5 != nbytes || dummy6 != nbytes {
            return Err(DcdError::Format(
                "Problem reading Z coordinate block headers".to_string(),
            ));
        }

        // Read Z coordinates in chunks
        let z_end = pos + natoms * 4;
        if z_end > self.mmap.len() {
            return Err(DcdError::Format(
                "Unexpected end of file reading Z coordinates".to_string(),
            ));
        }

        // Process each chunk of Z coordinates
        for chunk_start in (0..natoms).step_by(CHUNK_SIZE) {
            let chunk_end = std::cmp::min(chunk_start + CHUNK_SIZE, natoms);
            let chunk_size = chunk_end - chunk_start;

            // Read a chunk of coordinates
            for i in 0..chunk_size {
                let coord_pos = pos + (chunk_start + i) * 4;
                let value = read_f32_at(&self.mmap, coord_pos, self.endianness)?;
                self.z_buffer.push(value);
            }
        }
        pos = z_end;

        // Final read to verify the block
        let final_dummy = read_i32_at(&self.mmap, pos, self.endianness)?;
        pos += 4;

        // Validate the final marker
        if final_dummy != nbytes {
            return Err(DcdError::Format(
                "Problem validating coordinate block".to_string(),
            ));
        }

        // Efficiently construct Vec3 positions from the x, y, z buffers
        for i in 0..natoms {
            self.positions_buffer.push(Vec3::new(
                self.x_buffer[i],
                self.y_buffer[i],
                self.z_buffer[i],
            ));
        }

        Ok(pos)
    }

    /// Skips the next frame(s) instead of reading them
    ///
    /// This is more efficient than reading frames you don't need,
    /// as it just updates the file position without processing the data.
    ///
    /// # Arguments
    ///
    /// * `n` - Number of frames to skip (defaults to 1)
    ///
    /// # Returns
    ///
    /// A Result indicating success or an error
    pub fn skip_next(&mut self, n: Option<i32>) -> Result<(), DcdError> {
        let n = n.unwrap_or(1);
        if self.current_frame + n > self.nframes {
            return Err(DcdError::Format(
                "Attempted to skip beyond the end of the file".to_string(),
            ));
        }

        // Calculate new position directly - more efficient than reading headers
        self.current_pos += (self.framesize as usize) * (n as usize);
        self.current_frame += n;

        // Ensure we're at a valid position
        if self.current_pos > self.mmap.len() {
            return Err(DcdError::Format(
                "Unexpected end of file when skipping frames".to_string(),
            ));
        }

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

    /// Returns the endianness of the file
    pub fn endianness(&self) -> Endianness {
        self.endianness
    }

    /// Reset to the beginning of the trajectory
    ///
    /// This positions the reader at the first frame after the header.
    ///
    /// # Returns
    ///
    /// A Result indicating success or an error
    pub fn reset(&mut self) -> Result<(), DcdError> {
        self.current_frame = 0;
        // Position right after the header
        let header_size = self.read_header_size()?;
        self.current_pos = header_size;
        Ok(())
    }

    /// Get the size of the header in bytes
    ///
    /// This calculates the position immediately after the header.
    ///
    /// # Returns
    ///
    /// A Result containing the header size or an error
    fn read_header_size(&self) -> Result<usize, DcdError> {
        // Size calculation based on the DCD format
        // Header size depends on the number of title lines
        let ntitle = read_i32_at(&self.mmap, 96, self.endianness)?;
        let pos = 100 + (ntitle as usize) * 80 + 16;
        Ok(pos)
    }

    /// Skip to a specific frame in the trajectory
    ///
    /// This allows random access to any frame in the file.
    ///
    /// # Arguments
    ///
    /// * `frame_idx` - The frame index to seek to (0-based)
    ///
    /// # Returns
    ///
    /// A Result indicating success or an error
    pub fn seek_frame(&mut self, frame_idx: i32) -> Result<(), DcdError> {
        // Validate frame index
        if frame_idx < 0 || frame_idx >= self.nframes {
            return Err(DcdError::Format(format!(
                "Frame index {} out of range (0-{})",
                frame_idx,
                self.nframes - 1
            )));
        }

        // For efficiency, calculate the position directly for fast seeking
        // This is much faster than sequentially skipping frames
        let header_size = self.read_header_size()?;
        self.current_pos = header_size + (frame_idx as usize) * (self.framesize as usize);
        self.current_frame = frame_idx;

        // Validate position
        if self.current_pos > self.mmap.len() {
            return Err(DcdError::Format(
                "Seek position is beyond the end of file".to_string(),
            ));
        }

        Ok(())
    }
}

// Implement Iterator for DcdReader
// This allows using the reader in for loops and with other iterator methods
impl Iterator for DcdReader {
    type Item = Result<Frame, DcdError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.read_next() {
            Ok(Some(frame)) => Some(Ok(frame)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

// Helper functions for reading different data types from the memory-mapped file

/// Read an i32 value at a specific position with the specified endianness
///
/// # Arguments
///
/// * `mmap` - The memory-mapped file
/// * `pos` - Position to read from
/// * `endianness` - Endianness to use
///
/// # Returns
///
/// A Result containing the read value or an error
#[inline(always)]
fn read_i32_at(mmap: &[u8], pos: usize, endianness: Endianness) -> Result<i32, DcdError> {
    // Make sure we have enough data
    if pos + 4 > mmap.len() {
        return Err(DcdError::Format("Unexpected end of file".to_string()));
    }

    // Get the 4 bytes at the position
    let bytes = &mmap[pos..pos + 4];

    // Convert bytes to i32 using the appropriate endianness
    let value = match endianness {
        Endianness::Little => i32::from_le_bytes(bytes.try_into().unwrap()),
        Endianness::Big => i32::from_be_bytes(bytes.try_into().unwrap()),
    };

    Ok(value)
}

/// Read an f32 value at a specific position with the specified endianness
///
/// # Arguments
///
/// * `mmap` - The memory-mapped file
/// * `pos` - Position to read from
/// * `endianness` - Endianness to use
///
/// # Returns
///
/// A Result containing the read value or an error
#[inline(always)]
fn read_f32_at(mmap: &[u8], pos: usize, endianness: Endianness) -> Result<f32, DcdError> {
    // Make sure we have enough data
    if pos + 4 > mmap.len() {
        return Err(DcdError::Format("Unexpected end of file".to_string()));
    }

    // Get the 4 bytes at the position
    let bytes = &mmap[pos..pos + 4];

    // Convert bytes to f32 using the appropriate endianness
    let value = match endianness {
        Endianness::Little => f32::from_le_bytes(bytes.try_into().unwrap()),
        Endianness::Big => f32::from_be_bytes(bytes.try_into().unwrap()),
    };

    Ok(value)
}

/// Read an f64 value at a specific position with the specified endianness
///
/// # Arguments
///
/// * `mmap` - The memory-mapped file
/// * `pos` - Position to read from
/// * `endianness` - Endianness to use
///
/// # Returns
///
/// A Result containing the read value or an error
#[inline(always)]
fn read_f64_at(mmap: &[u8], pos: usize, endianness: Endianness) -> Result<f64, DcdError> {
    // Make sure we have enough data
    if pos + 8 > mmap.len() {
        return Err(DcdError::Format("Unexpected end of file".to_string()));
    }

    // Get the 8 bytes at the position
    let bytes = &mmap[pos..pos + 8];

    // Convert bytes to f64 using the appropriate endianness
    let value = match endianness {
        Endianness::Little => f64::from_le_bytes(bytes.try_into().unwrap()),
        Endianness::Big => f64::from_be_bytes(bytes.try_into().unwrap()),
    };

    Ok(value)
}
