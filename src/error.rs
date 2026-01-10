//! Error types for the Hull-White library

use std::fmt;

/// Errors that can occur in the Hull-White model calculations
#[derive(Debug, Clone, PartialEq)]
pub enum HullWhiteError {
    /// Error in root finding algorithms
    RootFindingError(String),
    /// Invalid input parameters
    InvalidInput(String),
    /// Numerical computation error
    NumericalError(String),
}

impl fmt::Display for HullWhiteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            HullWhiteError::RootFindingError(msg) => write!(f, "Root finding error: {}", msg),
            HullWhiteError::InvalidInput(msg) => write!(f, "Invalid input: {}", msg),
            HullWhiteError::NumericalError(msg) => write!(f, "Numerical error: {}", msg),
        }
    }
}

impl std::error::Error for HullWhiteError {}