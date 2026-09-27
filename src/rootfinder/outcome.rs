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
        /// Low end of the final bracket, in root (rate) units.
        lower: f64,
        /// High end of the final bracket, in root (rate) units.
        upper: f64,
        /// Objective evaluated at [`crate::SolverError::NoSignChange::lower`]; same sign as `f_upper`.
        f_lower: f64,
        /// Objective evaluated at [`crate::SolverError::NoSignChange::upper`]; same sign as `f_lower`.
        f_upper: f64,
    },
    /// The objective (or its derivative) produced a NaN, so no ordering is possible.
    NonFiniteEvaluation {
        /// Point the objective was evaluated at, in root (rate) units.
        point: f64,
        /// The non-finite value that came back (`NaN` or `±inf`).
        value: f64,
    },
    /// Bracketed and converging, but the iteration cap was hit first.
    Exhausted {
        /// Iterations spent before the cap (`crate::SolverSettings::max_iterations`) was reached.
        iterations: u32,
        /// Low end of the bracket as it stood when the cap was hit, in root (rate) units.
        lower: f64,
        /// High end of the bracket as it stood when the cap was hit, in root (rate) units.
        upper: f64,
        /// Objective magnitude at the last iterate — how close the miss was.
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
    /// The converged root, in the units of the variable being solved for (for the Jamshidian solve,
    /// a short rate).
    pub root: f64,
    /// Iterations the search took to get there — useful for spotting a solve that is converging on
    /// the last of its budget rather than comfortably inside it.
    pub iterations: u32,
    /// Objective magnitude at [`crate::Solution::root`].  Read it as a price residual, not as accuracy:
    /// the solve exits on the *rate*-space bracket width, so a small price residual on a steep
    /// objective can still hide a rate that is a long way from the root.
    pub residual: f64,
}
