//! # resmol - Rust Molecular Dynamics Library
//!
//! `resmol` provides efficient readers for common molecular dynamics file formats
//! such as DCD trajectories produced by CHARMM, NAMD, and LAMMPS.
//!
//! ## Features
//!
//! - High-performance DCD trajectory file reader
//! - Memory-mapped file access for efficient frame seeking
//! - Support for both little and big endian files
//! - Iterator-based API for idiomatic Rust usage
//!
//! ## Example
//!
//! ```no_run
//! use resmol::DcdReader;
//!
//! // Open a DCD trajectory file
//! let mut reader = DcdReader::open("trajectory.dcd").unwrap();
//!
//! // Iterate through all frames using the Iterator trait
//! for frame_result in reader {
//!     let frame = frame_result.unwrap();
//!     
//!     // Process the frame data
//!     println!("Frame timestep: {}", frame.timestep);
//!     println!("Contains {} atoms", frame.positions.len());
//! }
//! ```

// Feature flags to enable modular compilation
// The 'dcd' feature enables DCD trajectory file support
// This allows users to only include the code they need, reducing
// binary size and compile time for applications that don't need all features

// Re-export the main components for easier access
// This allows users to write 'use resmol::DcdReader' instead of 'use resmol::dcd::DcdReader'
#[cfg(feature = "dcd")]
pub use dcd::{DcdError, DcdReader, DcdReaderBuilder, Endianness, Frame};

// Modules that contain the actual implementation
// Each file format has its own module, conditionally compiled based on features
#[cfg(feature = "dcd")]
pub mod dcd;

// In the future, other modules can be added for different file formats:
// #[cfg(feature = "xtc")]
// pub mod xtc;
//
// #[cfg(feature = "pdb")]
// pub mod pdb;
