//! Acceptance-threshold scan: the measurement behind the `MAX_MODEL_ERROR_RATIO` table in the
//! harness module docs.
//!
//! For each candidate ratio `|f(x_newton)| <= r * |f(x)|`, run the whole production grid plus the
//! deep-tail fixture and print the mean iterations, total evaluations, worst root error against a
//! 1e-15 reference, and the tail iteration count.  Ratios above `1/e = 0.368` reopen the tail
//! crawl; the shipped `0.25` is the tightest threshold that costs no accuracy.
use crate::solver_ab::fixtures::{Counted, grid};
use crate::solver_ab::measurement::{derivative, objective};
use crate::solver_ab::variants::bracket_for;
use crate::{HullWhite, SolverSettings, rootfinder};
use std::sync::atomic::{AtomicUsize, Ordering};

fn solve_ratio(
    max_ratio: f64,
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
        let newton = if correction.is_finite() {
            x - correction
        } else {
            f64::NAN
        };
        let inside = newton.is_finite() && newton > low && newton < high;
        let makes_progress = inside && correction.abs() >= 0.1 * width;
        let improves = if inside && !makes_progress {
            let f_newton = f(newton);
            if max_ratio < 0.0 {
                f_newton.is_finite() && f_newton.abs() < f_x.abs()
            } else {
                f_newton.is_finite() && f_newton.abs() <= max_ratio * f_x.abs()
            }
        } else {
            false
        };
        if correction.is_finite() && correction.abs() <= settings.tolerance * scale {
            let root = x - correction;
            if root >= low && root <= high {
                return Ok(rootfinder::Solution {
                    root,
                    iterations: iteration,
                    residual: f_x.abs(),
                });
            }
        }
        x = if makes_progress || improves {
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

fn grid_run(r: f64) -> (f64, usize, f64) {
    let settings = SolverSettings {
        tolerance: 1e-12,
        max_iterations: 500,
        initial_guess: None,
    };
    let tight = SolverSettings {
        tolerance: 1e-15,
        max_iterations: 500,
        initial_guess: None,
    };
    let mut iters = 0u64;
    let mut evals = 0usize;
    let mut maxerr = 0f64;
    let mut n = 0usize;
    for case in grid() {
        let counted = Counted::new(case.curr, case.a, case.b, case.sig);
        let hw = HullWhite::init(
            case.a,
            case.sig,
            &counted.yield_curve,
            &counted.forward_curve,
        )
        .unwrap();
        let brk = hw.critical_rate_bracket(
            case.option_maturity,
            &case.schedule,
            case.coupon_rate,
            case.strike,
        );
        let seed = hw.mu_r(case.r_t, case.t, case.option_maturity).unwrap();
        let rfc = AtomicUsize::new(0);
        let rdc = AtomicUsize::new(0);
        let robj = objective(&hw, &case, &rfc);
        let rder = derivative(&hw, &case, &rdc);
        let ro1 = |x: f64| robj(x);
        let rd1 = |x: f64| rder(x);
        let _reference_only = rootfinder::solve(&ro1, &rd1, brk, seed, &tight).unwrap();
        println!(
            "   reference-only evals {}",
            rfc.load(Ordering::Relaxed) + rdc.load(Ordering::Relaxed)
        );
        let fc = AtomicUsize::new(0);
        let dc = AtomicUsize::new(0);
        let obj = objective(&hw, &case, &fc);
        let der = derivative(&hw, &case, &dc);
        let o1 = |x: f64| obj(x);
        let d1 = |x: f64| der(x);
        let reference = rootfinder::solve(&o1, &d1, brk, seed, &tight).unwrap();
        let o2 = |x: f64| obj(x);
        let d2 = |x: f64| der(x);
        let out = solve_ratio(r, &o2, &d2, brk, seed, &settings).unwrap();
        iters += out.iterations as u64;
        evals += fc.load(Ordering::Relaxed) + dc.load(Ordering::Relaxed);
        maxerr = maxerr.max((out.root - reference.root).abs());
        n += 1;
    }
    (iters as f64 / n as f64, evals, maxerr)
}

fn tail_run(r: f64) -> u32 {
    let settings = SolverSettings {
        tolerance: 1e-12,
        max_iterations: 500,
        initial_guess: None,
    };
    let tail_f = |x: f64| (-x).exp() - 1e-6;
    let tail_df = |x: f64| -(-x).exp();
    solve_ratio(r, &tail_f, &tail_df, None, -100.0, &settings)
        .map(|s| s.iterations)
        .unwrap_or(u32::MAX)
}

#[test]
fn scan_thresholds() {
    let tail_f = |x: f64| (-x).exp() - 1e-6;
    let tail_df = |x: f64| -(-x).exp();
    let settings = SolverSettings {
        tolerance: 1e-12,
        max_iterations: 500,
        initial_guess: None,
    };
    let shipped = rootfinder::solve(&tail_f, &tail_df, None, -100.0, &settings).unwrap();
    println!("SHIPPED tail iters {}", shipped.iterations);
    let (si, se, serr) = grid_run(-2.0); // sentinel: not used as a ratio, only for shipped below
    println!("SHIPPED grid (recompute) avg {si:.1} evals {se} maxerr {serr:.3e}");
    for r in [
        -1.0f64, 0.05, 0.10, 0.20, 0.25, 0.30, 0.33, 0.34, 0.3678, 0.37, 0.40, 0.5,
    ] {
        let (it, ev, err) = grid_run(r);
        let t = tail_run(r);
        println!("ratio {r:>8} -> grid it avg {it:.2}  evals {ev:<6} maxerr {err:.3e}  tail {t}");
    }
}
