//! Bracketed, safeguarded one-dimensional root finding.
//!
//! The Jamshidian decomposition of a coupon-bond option hinges on one scalar equation: find the
//! *critical rate* `r*` at which the bond's value at option expiry equals the strike.  That
//! equation is solved here rather than with an open Newton iteration, because an unbracketed solve
//! is what made the decomposition fragile:
//!
//! * a fixed starting rate (`R_INIT = 0.03`) tells the iteration nothing about the curve, so on a
//!   steep curve or a long tenor it walks a long way before it lands, and nothing checks where it
//!   landed;
//! * when the iteration struggled, the whole instrument came back as an error whose text was the
//!   solver's own numeric residue (`RootFindingError("NaN")`), which is indistinguishable from a
//!   model blow-up;
//! * the caller had no way to ask for more accuracy or more iterations.
//!
//! [`solve`] takes a *proposed* bracket plus a seed, verifies (and if necessary widens) the bracket
//! until it genuinely straddles a sign change, and then runs Newton guarded by bisection so the
//! iterate can never leave the bracket.  Failure is therefore rare, and when it happens it carries
//! the bracket, the iteration count and the residual instead of a bare number.

use core::fmt;

/// Convergence target on the root, relative to its magnitude: the bracket is closed to
/// `tolerance * max(1, |r|)`.  1e-12 on a rate is ~1e-8 basis points, which is far below any
/// price a Hull-White model is quoted to and far above f64 round-off (~1e-16).
pub const DEFAULT_TOLERANCE: f64 = 1e-12;

/// Iteration cap.  Bisection alone halves the bracket each pass, so this covers a bracket of
/// `1e15 * tolerance` wide; Newton typically finishes in well under ten.
pub const DEFAULT_MAX_ITERATIONS: u32 = 100;

/// How much of the bracket a Newton step has to cover to be worth taking *on distance alone*.
///
/// Without this, a safeguarded Newton can be *bracketed and still useless*: for an objective that
/// looks like `exp(-B r)` far from the root, the Newton step is `1/B` no matter how far away the
/// root is, so on a bracket hundreds of rate-units wide it crawls.  Requiring each step to cross a
/// tenth of the bracket means bisection does the long-distance work and Newton only finishes what
/// it is actually good at — the region where it converges quadratically.
///
/// This is a necessary rule but not a sufficient one on its own: the same test also rejects the step
/// that has *already landed on the root*, because there the correction is tiny.  See
/// [`MAX_MODEL_ERROR_RATIO`] for the other half of the acceptance test.
const MIN_NEWTON_FRACTION: f64 = 0.1;

/// How badly the tangent line is allowed to be wrong over a Newton step for the step to be trusted:
/// the step is accepted when `|f(x_newton)| <= MAX_MODEL_ERROR_RATIO * |f(x)|`.
///
/// Newton's step is built from the tangent, which predicts `f = 0` at `x_newton`; the *actual*
/// `|f(x_newton)|` is therefore the tangent's error, expressed in the same units as the residual
/// being corrected.  Requiring that error to be a quarter of the residual says "the linear model of
/// the objective is good to 25% over the distance this step covers", which is exactly the property
/// that separates the two cases the width rule cannot tell apart:
///
/// * **near the root** the step error is quadratic in the distance left, so the residual drops by
///   orders of magnitude in one step — far inside 25%.  Accepted, and the solve is over in a few
///   passes instead of bisecting an answer it already had.
/// * **deep in a tail** the objective is locally exponential, and an exponential contracts by
///   `exp(-B * step) = exp(-1) = 0.368` per Newton step *no matter how far away the root is*.  That
///   is outside 25%, so the crawl is rejected and bisection keeps halving the distance.
///
/// The number is a trust threshold, not an accuracy promise — accuracy is still demanded in rate
/// space, by [`SolverSettings::tolerance`].  A step that fails this test is not discarded, it just
/// falls back to bisection, which cannot be talked into a crawl.
///
/// `0.25` is a measured value, not a guess: the scan in the `solver_ab` harness module docs shows
/// every threshold above `1/e = 0.368` reopens the tail crawl, while the root agrees with a
/// far-tighter reference across the whole range.  Going well below the `0.33` knee only costs
/// evaluations (`0.10` costs ~19% more on the production grid), so `0.25` buys distance from the
/// contraction rate almost for free.  The extra `f` evaluation this test needs is paid only when
/// the width rule has already failed, so a step accepted on distance never pays for it.
const MAX_MODEL_ERROR_RATIO: f64 = 0.25;

/// How far the geometric search will widen around the seed looking for a sign change.  Starting
/// from `1e-4 * max(1, |seed|)` and doubling, this reaches ~1e116 before giving up.
const MAX_EXPANSIONS: u32 = 400;

/// Knobs for [`solve`], carried on the model so callers can tune the Jamshidian solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolverSettings {
    /// Convergence tolerance on the root, relative to its magnitude (see [`DEFAULT_TOLERANCE`]).
    pub tolerance: f64,
    /// Hard cap on solver iterations (see [`DEFAULT_MAX_ITERATIONS`]).
    pub max_iterations: u32,
    /// Overrides the starting point.  `None` means "let the caller pass the seed", which for the
    /// Jamshidian solve is the model's own expected short rate at expiry.
    pub initial_guess: Option<f64>,
}

impl Default for SolverSettings {
    fn default() -> Self {
        Self {
            tolerance: DEFAULT_TOLERANCE,
            max_iterations: DEFAULT_MAX_ITERATIONS,
            initial_guess: None,
        }
    }
}

impl SolverSettings {
    /// Checks the knobs before a solve starts, so a bad `tolerance` is reported as bad input
    /// rather than silently returning an unconverged root.
    pub fn validate(&self) -> Result<(), &'static str> {
        if !self.tolerance.is_finite() || self.tolerance <= 0.0 {
            return Err("solver tolerance must be finite and > 0");
        }
        if self.max_iterations < 1 {
            return Err("solver max_iterations must be >= 1");
        }
        if let Some(guess) = self.initial_guess
            && !guess.is_finite()
        {
            return Err("solver initial_guess must be finite");
        }
        Ok(())
    }
}

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

fn straddles(f_lower: f64, f_upper: f64) -> bool {
    //NaN comparisons are false, so a NaN end point simply fails to straddle.  Infinities are fine:
    //they still carry a usable sign.
    (f_lower <= 0.0 && f_upper >= 0.0) || (f_lower >= 0.0 && f_upper <= 0.0)
}

/// Solves `f(x) = 0` for `x` inside a bracket, using `df` for Newton steps with bisection as the
/// guard.
///
/// * `proposed_bracket` is an optional `[lo, hi]` the caller believes contains the root.  It is
///   verified, not trusted; if it does not straddle, it is discarded in favour of widening.
/// * `seed` is the preferred starting point.  The result is not sensitive to it — the guard keeps
///   every iterate inside the bracket, and bisection converges regardless.
///
/// The function must be continuous on the bracket.  Monotonicity is not required, but with more
/// than one root in the bracket the bracket decides *which* root is returned.
pub fn solve(
    f: &dyn Fn(f64) -> f64,
    df: &dyn Fn(f64) -> f64,
    proposed_bracket: Option<(f64, f64)>,
    seed: f64,
    settings: &SolverSettings,
) -> Result<Solution, SolverError> {
    if !seed.is_finite() {
        return Err(SolverError::NonFiniteEvaluation {
            point: seed,
            value: seed,
        });
    }

    let f_seed = f(seed);
    if f_seed.is_nan() {
        return Err(SolverError::NonFiniteEvaluation {
            point: seed,
            value: f_seed,
        });
    }
    if f_seed == 0.0 {
        return Ok(Solution {
            root: seed,
            iterations: 0,
            residual: 0.0,
        });
    }

    //Try the analytic bracket first: it costs one evaluation at each end and, when it works, the
    //root is already trapped inside a tight interval.
    let mut bracket = match proposed_bracket {
        Some((lo, hi)) if lo.is_finite() && hi.is_finite() && lo < hi => {
            let (f_lo, f_hi) = (f(lo), f(hi));
            if f_lo.is_nan() || f_hi.is_nan() {
                return Err(SolverError::NonFiniteEvaluation {
                    point: if f_lo.is_nan() { lo } else { hi },
                    value: if f_lo.is_nan() { f_lo } else { f_hi },
                });
            }
            if f_lo == 0.0 {
                return Ok(Solution {
                    root: lo,
                    iterations: 0,
                    residual: 0.0,
                });
            }
            if f_hi == 0.0 {
                return Ok(Solution {
                    root: hi,
                    iterations: 0,
                    residual: 0.0,
                });
            }
            if straddles(f_lo, f_hi) {
                Some((lo, hi, f_lo, f_hi))
            } else {
                None
            }
        }
        _ => None,
    };

    //No usable bracket: widen geometrically around the seed until the ends disagree in sign.
    let mut width = match bracket {
        Some((lo, hi, _, _)) => (hi - lo) * 0.5,
        None => 1e-4 * seed.abs().max(1.0),
    };
    let (mut last_lo, mut last_hi, mut f_last_lo, mut f_last_hi) =
        (seed - width, seed + width, f_seed, f_seed);
    let mut expansions = 0;
    while bracket.is_none() && expansions < MAX_EXPANSIONS {
        expansions += 1;
        width *= 2.0;
        last_lo = seed - width;
        last_hi = seed + width;
        f_last_lo = f(last_lo);
        f_last_hi = f(last_hi);
        if f_last_lo.is_nan() || f_last_hi.is_nan() {
            return Err(SolverError::NonFiniteEvaluation {
                point: if f_last_lo.is_nan() { last_lo } else { last_hi },
                value: if f_last_lo.is_nan() {
                    f_last_lo
                } else {
                    f_last_hi
                },
            });
        }
        if straddles(f_last_lo, f_last_hi) {
            bracket = Some((last_lo, last_hi, f_last_lo, f_last_hi));
        }
    }
    let (mut a, mut b, mut f_a, mut f_b) = match bracket {
        Some(bracket) => bracket,
        None => {
            return Err(SolverError::NoSignChange {
                lower: last_lo,
                upper: last_hi,
                f_lower: f_last_lo,
                f_upper: f_last_hi,
            });
        }
    };

    let mut x = seed.clamp(a.min(b), a.max(b));
    for iteration in 1..=settings.max_iterations {
        let f_x = f(x);
        if f_x.is_nan() {
            return Err(SolverError::NonFiniteEvaluation {
                point: x,
                value: f_x,
            });
        }
        if f_x == 0.0 {
            return Ok(Solution {
                root: x,
                iterations: iteration,
                residual: 0.0,
            });
        }
        //Re-trap the bracket around the iterate.
        if (f_x < 0.0) == (f_a < 0.0) {
            a = x;
            f_a = f_x;
        } else {
            b = x;
            f_b = f_x;
        }
        let (low, high) = (a.min(b), a.max(b));
        let scale = x.abs().max(1.0);
        if (high - low) <= settings.tolerance * scale {
            return Ok(Solution {
                root: x,
                iterations: iteration,
                residual: f_x.abs(),
            });
        }
        let slope = df(x);
        let width = high - low;
        //Newton's own estimate of how far the root still is.  Everything below is written in terms of
        //this rather than the raw step so that the two acceptance tests stay in comparable units.
        let correction = if slope.is_finite() && slope != 0.0 {
            f_x / slope
        } else {
            f64::NAN
        };
        //Exit: converged on the correction.  With an honest derivative `correction` *is* the distance
        //left, so a sub-tolerance correction is a converged root, not a stalled iteration.  Reading
        //a tiny correction as "no progress" is what made the shipped loop bisect an answer it had
        //already computed.  Note this is a rate-space exit -- it says the *unknown* is pinned down,
        //not that the price residual is small; see the note on `solve_v3` in the A/B harness for why
        //a price-space exit here would be unsafe.
        if correction.is_finite() && correction.abs() <= settings.tolerance * scale {
            let root = x - correction;
            //Never hand back a point outside the bracket that was actually verified to straddle.
            if root >= low && root <= high {
                return Ok(Solution {
                    root,
                    iterations: iteration,
                    residual: f_x.abs(),
                });
            }
        }
        let newton = if correction.is_finite() {
            x - correction
        } else {
            f64::NAN
        };
        //Step acceptance.  The step must stay inside the bracket -- that is the safeguard, and it is
        //not negotiable: it stops a lying derivative or a bad seed from throwing the iterate into
        //another basin.  Given that, either of two things makes the Newton step the right one:
        //
        //  * `makes_progress` -- it crosses a tenth of the bracket, so it is doing real long-range
        //    work by distance (the original guard); or
        //  * `improves` -- it lands where the tangent line predicted the objective to within a
        //    quarter of the residual being corrected, which is the signature of being inside the
        //    quadratic region even though the step is far too small to cross the bracket.
        //
        //Neither holds -> bisect, which is slow but cannot be argued into a crawl.
        let inside = newton.is_finite() && newton > low && newton < high;
        let makes_progress = inside && correction.abs() >= MIN_NEWTON_FRACTION * width;
        let improves = if inside && !makes_progress {
            let f_newton = f(newton);
            f_newton.is_finite() && f_newton.abs() <= MAX_MODEL_ERROR_RATIO * f_x.abs()
        } else {
            false
        };
        x = if makes_progress || improves {
            newton
        } else {
            0.5 * (low + high)
        };
    }

    Err(SolverError::Exhausted {
        iterations: settings.max_iterations,
        lower: a.min(b),
        upper: a.max(b),
        //The tightest end of the bracket is the honest thing to quote as the residual.
        residual: f_a.abs().min(f_b.abs()),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn settings(tolerance: f64, max_iterations: u32) -> SolverSettings {
        SolverSettings {
            tolerance,
            max_iterations,
            initial_guess: None,
        }
    }

    #[test]
    fn solves_a_simple_monotone_equation() {
        //f(x) = 2x - 1, root at 0.5, no bracket offered.
        let solution = solve(
            &|x| 2.0 * x - 1.0,
            &|_| 2.0,
            None,
            0.0,
            &settings(1e-12, 100),
        )
        .unwrap();
        assert!((solution.root - 0.5).abs() < 1e-12, "{solution:?}");
    }

    #[test]
    fn a_seed_far_from_the_root_lands_on_the_same_answer() {
        //f(x) = exp(-x) - 0.3, root at ln(1/0.3) = 1.2039728043259361
        let expected = (1.0 / 0.3_f64).ln();
        let f = |x: f64| (-x).exp() - 0.3;
        let df = |x: f64| -(-x).exp();
        for seed in [-5.0, 0.0, 40.0, 500.0] {
            let solution = solve(&f, &df, None, seed, &settings(1e-12, 100)).unwrap();
            assert!(
                (solution.root - expected).abs() < 1e-10,
                "seed {seed} gave {solution:?}, expected {expected}"
            );
        }
    }

    #[test]
    fn the_bracket_guard_stops_a_lying_derivative() {
        //A derivative of the wrong sign would send a plain Newton iteration off to infinity; the
        //bisection guard has to keep it on the root anyway.
        let solution = solve(
            &|x| x - 0.25,
            &|_| -1.0, //deliberately wrong sign
            Some((0.0, 1.0)),
            0.9,
            &settings(1e-12, 200),
        )
        .unwrap();
        assert!((solution.root - 0.25).abs() < 1e-10, "{solution:?}");
    }

    #[test]
    fn newton_cannot_crawl_an_exponential_tail() {
        //f(x) = exp(-x) - 1e-6 has root x = ln(1e6) = 13.8155.  Deep in the tail the Newton step
        //is exactly +1 per pass regardless of how far the root is, so an iteration that is merely
        //"bracketed" still needs ~140 passes to cross a bracket this wide -- and a million for a
        //strike this extreme in the real model.  The minimum-progress guard hands the long distance
        //to bisection instead, which halves the bracket every pass.
        let root = 1e6_f64.ln();
        let f = |x: f64| (-x).exp() - 1e-6;
        let df = |x: f64| -(-x).exp();
        let solution = solve(&f, &df, None, -100.0, &settings(1e-12, 100)).unwrap();
        assert!(
            (solution.root - root).abs() < 1e-10,
            "{solution:?} vs root {root}"
        );
        assert!(solution.iterations <= 70, "{solution:?}");
    }

    #[test]
    fn a_bracket_that_does_not_straddle_is_widened_not_trusted() {
        //The proposed bracket is entirely on the positive side of the root at -3.
        let solution = solve(
            &|x| x + 3.0,
            &|_| 1.0,
            Some((0.0, 1.0)),
            0.5,
            &settings(1e-12, 100),
        )
        .unwrap();
        assert!((solution.root + 3.0).abs() < 1e-9, "{solution:?}");
    }

    #[test]
    fn running_out_of_iterations_says_so_with_context() {
        let error = solve(
            &|x| (-x).exp() - 1e-6,
            &|x| -(-x).exp(),
            None,
            0.0,
            &settings(1e-14, 2),
        )
        .unwrap_err();
        assert!(
            matches!(error, SolverError::Exhausted { iterations: 2, .. }),
            "{error:?}"
        );
        let text = error.to_string();
        assert!(text.contains("2 iterations"), "{text}");
        assert!(text.contains("residual"), "{text}");
    }

    #[test]
    fn a_never_straddling_objective_reports_no_sign_change() {
        //Strictly positive: no root exists at any rate.
        let error = solve(
            &|x| x * x + 1.0,
            &|x| 2.0 * x,
            None,
            0.0,
            &settings(1e-12, 100),
        )
        .unwrap_err();
        assert!(
            matches!(error, SolverError::NoSignChange { .. }),
            "{error:?}"
        );
        assert!(error.to_string().contains("no sign change"), "{error}");
    }

    #[test]
    fn a_nan_objective_is_reported_as_such() {
        let error = solve(
            &|x| if x > 0.0 { f64::NAN } else { x + 1.0 },
            &|_| 1.0,
            None,
            0.0,
            &settings(1e-12, 100),
        )
        .unwrap_err();
        assert!(
            matches!(error, SolverError::NonFiniteEvaluation { .. }),
            "{error:?}"
        );
    }

    #[test]
    fn an_exact_endpoint_or_seed_root_short_circuits() {
        let settings = settings(1e-12, 100);
        let solution = solve(&|x| x - 1.0, &|_| 1.0, Some((1.0, 2.0)), 5.0, &settings).unwrap();
        assert_eq!(solution.root, 1.0);
        assert_eq!(solution.iterations, 0);
        let solution = solve(&|x| x - 7.0, &|_| 1.0, None, 7.0, &settings).unwrap();
        assert_eq!(solution.root, 7.0);
        assert_eq!(solution.iterations, 0);
    }

    #[test]
    fn settings_are_validated_before_the_solve() {
        assert!(settings(0.0, 10).validate().is_err());
        assert!(settings(f64::NAN, 10).validate().is_err());
        assert!(settings(1e-9, 0).validate().is_err());
        assert!(
            SolverSettings {
                tolerance: 1e-9,
                max_iterations: 10,
                initial_guess: Some(f64::INFINITY),
            }
            .validate()
            .is_err()
        );
        assert!(settings(1e-9, 10).validate().is_ok());
    }

    #[test]
    fn a_looser_tolerance_stops_looser() {
        //A nonlinear objective, so the number of iterations actually depends on the tolerance:
        //exp(x) = 2, root at ln 2.
        let f = |x: f64| x.exp() - 2.0;
        let df = |x: f64| x.exp();
        let root = 2.0_f64.ln();
        let loose = solve(&f, &df, Some((0.0, 2.0)), 0.0, &settings(1e-2, 100)).unwrap();
        let tight = solve(&f, &df, Some((0.0, 2.0)), 0.0, &settings(1e-14, 100)).unwrap();
        assert!(
            (loose.root - root).abs() <= 1e-2 && (tight.root - root).abs() <= 1e-13,
            "loose {loose:?} / tight {tight:?} vs root {root}"
        );
        assert!(loose.iterations < tight.iterations);
    }
}
