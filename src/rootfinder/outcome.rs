//! The result types of a solve: the converged [`Solution`] and the [`SolverError`] saying why
//! there isn't one.
//!
//! These live apart from the loop itself because the harness code (`solver_ab`, and every pricer
//! that forwards a failure as `HullWhiteError::RootFindingError`) consumes them without caring how
//! the search was run.  They are re-exported from `rootfinder`, so `rootfinder::SolverError` and
//! the crate-root `hull_white::SolverError` paths both keep working.

use core::fmt;

use super::MAX_EXPANSIONS;

/// Why the solve did not produce a root.  Every variant says enough to tell a bracketing failure
/// apart from a precision failure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SolverError {
    /// The bracket never straddled a sign change, however far it was widened.  Usually means the
    /// objective cannot take both signs (so the instrument has no critical rate).
    NoSignChange {
        lower: f64,
        upper: f64,
        f_lower: f64,
        f_upper: f64,
    },
    /// The objective (or its derivative) produced a NaN, so no ordering is possible.
    NonFiniteEvaluation { point: f64, value: f64 },
    /// Bracketed and converging, but the iteration cap was hit first.
    Exhausted {
        iterations: u32,
        lower: f64,
        upper: f64,
        residual: f64,
    },
}

impl fmt::Display for SolverError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SolverError::NoSignChange {
                lower,
                upper,
                f_lower,
                f_upper,
            } => write!(
                f,
                "no sign change in [{lower}, {upper}] after {MAX_EXPANSIONS} widening steps \
                 (f(lower) = {f_lower}, f(upper) = {f_upper}): the objective never crosses zero, so \
                 there is no critical rate for this instrument"
            ),
            SolverError::NonFiniteEvaluation { point, value } => write!(
                f,
                "the objective produced {value} at {point}; cannot bracket a root around a NaN"
            ),
            SolverError::Exhausted {
                iterations,
                lower,
                upper,
                residual,
            } => write!(
                f,
                "not converged in {iterations} iterations; bracket is [{lower}, {upper}] \
                 (width {:.3e}) with residual {residual:.3e} — raise max_iterations or loosen the \
                 tolerance",
                (upper - lower).abs()
            ),
        }
    }
}

impl std::error::Error for SolverError {}

/// A converged root, with enough metadata that a caller can see how hard it was to get.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Solution {
    pub root: f64,
    pub iterations: u32,
    pub residual: f64,
}
