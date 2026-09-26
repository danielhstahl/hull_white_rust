//! The candidate solve loops the harness measures against the shipped `rootfinder::solve`.
//!
//! Harness-local, never shipped.  `solve_v3` is kept deliberately unsafe as the documented
//! counter-example for why an exit must stay in rate space; `solve_v2` and `solve_v4` are the
//! steps that led to the shipped rule, and `guards` pins the shipped rule against them.

use crate::solver_ab::measurement::straddles;
use crate::{SolverSettings, rootfinder};

/// Bracket us or widen around the seed, shared by the candidate loops.
pub(super) fn bracket_for(
    f: &dyn Fn(f64) -> f64,
    bracket: Option<(f64, f64)>,
    seed: f64,
) -> Result<(f64, f64), rootfinder::SolverError> {
    if let Some((lo, hi)) = bracket
        && lo.is_finite()
        && hi.is_finite()
        && lo < hi
    {
        let (f_lo, f_hi) = (f(lo), f(hi));
        if !f_lo.is_nan() && !f_hi.is_nan() && straddles(f_lo, f_hi) {
            return Ok((lo, hi));
        }
    }
    let mut width = 1e-4 * seed.abs().max(1.0);
    for _ in 0..400 {
        width *= 2.0;
        let (lo, hi) = (seed - width, seed + width);
        let (f_lo, f_hi) = (f(lo), f(hi));
        if f_lo.is_nan() || f_hi.is_nan() {
            return Err(rootfinder::SolverError::NonFiniteEvaluation {
                point: lo,
                value: f_lo,
            });
        }
        if straddles(f_lo, f_hi) {
            return Ok((lo, hi));
        }
    }
    Err(rootfinder::SolverError::NoSignChange {
        lower: seed - width,
        upper: seed + width,
        f_lower: f(seed - width),
        f_upper: f(seed + width),
    })
}

/// A candidate second-round solve: the same bracketed, safeguarded loop as `rootfinder::solve`, plus
/// the two changes the measurements say are missing.
///
/// 1. **Stop on the Newton correction, not only on bracket width.**  `solve` exits when the
///    bracket is `tolerance * max(1, |x|)` wide.  Near the root the objective is very nearly
///    linear, so a Newton correction there is *tiny* -- which the min-progress guard reads as
///    "no progress" and refuses.  The solver then bisects the remaining width one halving at a
///    time: ~40 passes to take a 0.4-wide bracket to 1e-12, to confirm an answer it already had.
///    A correction below tolerance *is* convergence, so return it.
/// 2. **Let a residual-reducing Newton step through once it is inside the bracket and the
///    bracket itself is already tight in relative terms.**  Without this the tail-crawl guard and
///    the final convergence look identical to the solver.
pub(super) fn solve_v2(
    f: &dyn Fn(f64) -> f64,
    df: &dyn Fn(f64) -> f64,
    bracket: Option<(f64, f64)>,
    seed: f64,
    settings: &SolverSettings,
) -> Result<rootfinder::Solution, rootfinder::SolverError> {
    let seed_root = |s: f64| -> Option<(f64, f64)> {
        if !s.is_finite() {
            return None;
        }
        if let Some((lo, hi)) = bracket
            && lo.is_finite()
            && hi.is_finite()
            && lo < hi
        {
            let (f_lo, f_hi) = (f(lo), f(hi));
            if !f_lo.is_nan() && !f_hi.is_nan() && straddles(f_lo, f_hi) {
                return Some((lo, hi));
            }
        }
        //No usable bracket: widen geometrically around the seed, as the shipped solver does.
        let mut width = 1e-4 * s.abs().max(1.0);
        for _ in 0..400 {
            width *= 2.0;
            let (lo, hi) = (s - width, s + width);
            let (f_lo, f_hi) = (f(lo), f(hi));
            if f_lo.is_nan() || f_hi.is_nan() {
                return None;
            }
            if straddles(f_lo, f_hi) {
                return Some((lo, hi));
            }
        }
        None
    };
    let (mut a, mut b) = match seed_root(seed) {
        Some(range) => range,
        None => {
            return Err(rootfinder::SolverError::NoSignChange {
                lower: seed,
                upper: seed,
                f_lower: f(seed),
                f_upper: f(seed),
            });
        }
    };
    let mut f_a = f(a);
    let mut f_b = f(b);
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
            f_b = f_x;
        }
        let (low, high) = (a.min(b), a.max(b));
        let scale = x.abs().max(1.0);
        let width = high - low;
        if width <= settings.tolerance * scale {
            return Ok(rootfinder::Solution {
                root: x,
                iterations: iteration,
                residual: f_x.abs(),
            });
        }
        let slope = df(x);
        let correction = if slope.is_finite() && slope != 0.0 {
            f_x / slope
        } else {
            f64::NAN
        };
        //Change (1): a sub-tolerance correction is a converged root, not a stalled one.
        if correction.is_finite() && correction.abs() <= settings.tolerance * scale {
            return Ok(rootfinder::Solution {
                root: x - correction,
                iterations: iteration,
                residual: f_x.abs(),
            });
        }
        let newton = if correction.is_finite() {
            x - correction
        } else {
            f64::NAN
        };
        let inside = newton.is_finite() && newton > low && newton < high;
        //Change (2): inside the bracket is enough once the bracket is already tight; the
        //min-progress requirement stays while it is wide, which is what kills the tail crawl.
        let tight = width <= 1e-4 * scale;
        let makes_progress = newton.is_finite() && correction.abs() >= 0.1 * width;
        x = if inside && (tight || makes_progress) {
            newton
        } else {
            0.5 * (low + high)
        };
    }
    Err(rootfinder::SolverError::Exhausted {
        iterations: settings.max_iterations,
        lower: a.min(b),
        upper: a.max(b),
        residual: f_a.abs().min(f_b.abs()),
    })
}

/// v3 = shipped loop + a residual exit.  **Kept as a documented counter-example: it is unsafe.**
///
/// The Jamshidian condition *is* `f(r*) = 0`, so a price-space exit looks natural.  It is not
/// adequate: `f` can be steep, so a tiny price residual bounds nothing about the rate.  On the
/// deep-tail fixture a 1e-8 residual exit stopped 7e-3 away from the root, which would put every
/// leg strike in the decomposition wrong even though the objective "zeroed".
pub(super) fn solve_v3(
    f: &dyn Fn(f64) -> f64,
    df: &dyn Fn(f64) -> f64,
    bracket: Option<(f64, f64)>,
    seed: f64,
    settings: &SolverSettings,
    resid_tol: f64,
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
        if f_x.abs() <= resid_tol {
            return Ok(rootfinder::Solution {
                root: x,
                iterations: iteration,
                residual: f_x.abs(),
            });
        }
        if (f_x < 0.0) == (f_a < 0.0) {
            a = x;
            f_a = f_x;
        } else {
            b = x;
        }
        let (low, high) = (a.min(b), a.max(b));
        let scale = x.abs().max(1.0);
        let width = high - low;
        if width <= settings.tolerance * scale {
            return Ok(rootfinder::Solution {
                root: x,
                iterations: iteration,
                residual: f_x.abs(),
            });
        }
        let slope = df(x);
        let correction = if slope.is_finite() && slope != 0.0 {
            f_x / slope
        } else {
            f64::NAN
        };
        if correction.is_finite() && correction.abs() <= settings.tolerance * scale {
            return Ok(rootfinder::Solution {
                root: x - correction,
                iterations: iteration,
                residual: f_x.abs(),
            });
        }
        let newton = if correction.is_finite() {
            x - correction
        } else {
            f64::NAN
        };
        let inside = newton.is_finite() && newton > low && newton < high;
        let makes_progress = newton.is_finite() && correction.abs() >= 0.1 * width;
        x = if inside && makes_progress {
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

/// v4 = the shipped loop with the min-progress guard demoted to what it is actually for.
///
/// The guard exists to stop a far-tail Newton crawling a fixed small step across a huge bracket.
/// It is implemented as a *lower bound on the step relative to the bracket width*, which is the
/// wrong quantity on realistic inputs: the first Newton step from `mu_r` lands within ~2e-4 of
/// the root, gets rejected for being "too small", and the rest of the run is bisection grinding
/// the stale far end down to the tolerance.
///
/// This variant keeps the same exits -- all in rate space -- but accepts a Newton step whenever it
/// lands inside the bracket AND strictly reduces the residual, in addition to the width rule.  The
/// residual decrease is a step-acceptance test, never an exit test; see `solve_v3` for why an exit
/// on residual is unsafe.
///
/// What actually ships is this rule with a *tighter* acceptance threshold -- `|f(x_newton)| <=
/// 0.25 * |f(x)|` rather than any decrease -- because "any decrease" admits the `exp(-1) = 0.368`
/// per-step contraction of a tail and so leaves the crawl open.  See the threshold scan in the
/// module docs and `MAX_MODEL_ERROR_RATIO` in `rootfinder.rs`.
pub(super) fn solve_v4(
    f: &dyn Fn(f64) -> f64,
    df: &dyn Fn(f64) -> f64,
    bracket: Option<(f64, f64)>,
    seed: f64,
    settings: &SolverSettings,
    _unused: f64,
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
        let scale = x.abs().max(1.0);
        let width = high - low;
        if width <= settings.tolerance * scale {
            return Ok(rootfinder::Solution {
                root: x,
                iterations: iteration,
                residual: f_x.abs(),
            });
        }
        let slope = df(x);
        let correction = if slope.is_finite() && slope != 0.0 {
            f_x / slope
        } else {
            f64::NAN
        };
        if correction.is_finite() && correction.abs() <= settings.tolerance * scale {
            return Ok(rootfinder::Solution {
                root: x - correction,
                iterations: iteration,
                residual: f_x.abs(),
            });
        }
        let newton = if correction.is_finite() {
            x - correction
        } else {
            f64::NAN
        };
        let inside = newton.is_finite() && newton > low && newton < high;
        let f_newton = if inside { f(newton) } else { f64::NAN };
        let improves = f_newton.is_finite() && f_newton.abs() < f_x.abs();
        let makes_progress = correction.is_finite() && correction.abs() >= 0.1 * width;
        x = if inside && (improves || makes_progress) {
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
