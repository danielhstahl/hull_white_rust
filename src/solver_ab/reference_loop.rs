//! The pre-port solve loop, kept verbatim as an independent reference implementation.
//!
//! This is `rootfinder::solve` as it shipped before the step-rule port: `inside && makes_progress`,
//! with no correction-based exit.  It is not a candidate and it is not fast -- 53.1 iterations /
//! 9,092 evaluations per solve on the production grid, because it rejects the step that has
//! already landed on the root and then bisects a stale far bracket end down to the rate tolerance.
//! It is kept because it is a *separate* loop with the same contract, so agreement between it and
//! the ported solver is evidence about the root rather than an artefact of shared code.

use crate::solver_ab::variants::bracket_for;
use crate::{SolverSettings, rootfinder};

/// The solve loop exactly as it shipped before the step-rule port: `inside && makes_progress`, with
/// no correction-based exit.
///
/// This is not a candidate and it is not fast -- 53.1 iterations / 9,092 evaluations per solve on
/// the production grid, because it rejects the step that has already landed on the root and then
/// bisects a stale far bracket end down to the rate tolerance.  It is kept as the *reference
/// implementation* for the guard below: an independent loop with the same contract, so agreement
/// between it and the ported `rootfinder::solve` is evidence about the root rather than an
/// artefact of shared code.
pub(super) fn solve_v0(
    f: &dyn Fn(f64) -> f64,
    df: &dyn Fn(f64) -> f64,
    bracket: Option<(f64, f64)>,
    seed: f64,
    settings: &SolverSettings,
) -> Result<rootfinder::Solution, rootfinder::SolverError> {
    let (mut a, mut b) = bracket_for(f, bracket, seed)?;
    let mut f_a = f(a);
    let mut x = seed.clamp(a.min(b), a.max(b));
    for iteration in 1..=settings.max_iterations {
        let f_x = f(x);
        if f_x.is_nan() {
            return Err(rootfinder::SolverError::NonFiniteEvaluation {
                point: x,
                value: f_x,
            });
        }
        if f_x == 0.0 {
            return Ok(rootfinder::Solution {
                root: x,
                iterations: iteration,
                residual: 0.0,
            });
        }
        if (f_x < 0.0) == (f_a < 0.0) {
            a = x;
            f_a = f_x;
        } else {
            b = x;
        }
        let (low, high) = (a.min(b), a.max(b));
        if (high - low) <= settings.tolerance * x.abs().max(1.0) {
            return Ok(rootfinder::Solution {
                root: x,
                iterations: iteration,
                residual: f_x.abs(),
            });
        }
        let slope = df(x);
        let width = high - low;
        let newton = if slope.is_finite() && slope != 0.0 {
            x - f_x / slope
        } else {
            f64::NAN
        };
        let makes_progress = newton.is_finite() && (newton - x).abs() >= 0.1 * width;
        x = if newton.is_finite() && newton > low && newton < high && makes_progress {
            newton
        } else {
            0.5 * (low + high)
        };
    }
    Err(rootfinder::SolverError::Exhausted {
        iterations: settings.max_iterations,
        lower: a.min(b),
        upper: a.max(b),
        residual: f_a.abs(),
    })
}
