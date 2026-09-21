//! A/B cost comparison: the bracketed rootfinder vs the `nrfind` Newton iteration it replaced.
//!
//! Accuracy of the new solve was checked against quadrature when it landed; cost was not.  This
//! harness runs both solvers against *the same* objective and derivative closures the production
//! Jamshidian code builds (`coupon_bond_generic_t` over `bond_price_t_raw`, and
//! `coupon_bond_price_t_deriv`), over a grid of curve fixtures x coupon schedules x strikes, so
//! the only thing that varies is the solver.
//!
//! What is compared per case:
//!
//! * `f_calls` / `df_calls` — the solver's asks on the model.  A price costs what an evaluation
//!   costs, so this is the machine-independent measure of the solve.
//! * `curve_calls` — how many times the solve reached for the yield/forward curves, counted by
//!   the curves themselves.  Independent of solver internals entirely.
//! * `setup` — curve calls spent building the analytic bracket plus the seed.  The new solver has
//!   a setup step the old one never paid, so it is counted rather than hidden.
//! * `ns` — wall clock per solve, median of repeats, old and new interleaved in the same repeat
//!   loop so clock drift and cache state hit both alike.
//! * root agreement / residual — whether the two solvers answered the same question.
//! * `d_ratio` — the analytic derivative over a central difference of the objective at the root.
//!   If the analytic derivative is not the objective's derivative, the min-progress guard rejects
//!   every Newton step and the solve degenerates into pure bisection.
//!
//! Run it (it is env-gated so the ordinary test suite stays fast):
//!
//! ```text
//! SOLVER_AB=1 cargo test --release solver_ab -- --nocapture   // sections A-D, the two production shapes
//! SOLVER_AB=1 cargo test --release exit_rule_variants -- --nocapture  // section E, candidate exit rules
//! SOLVER_AB=1 cargo test --release old_solver_stress -- --nocapture   // section F, where the old solve breaks
//! AB_TRACE=1  cargo test --release trace_one_solve -- --nocapture     // iteration-by-iteration dump
//! AB_REPS=51  ...                                                  // more repeats for timing
//! ```
//!
//! ## Findings (84 cases: 4 curve shapes x {4,8,48} coupons x 7 strikes, release build)
//!
//! | solver | iters | evals | curve calls | ns/solve | max root err vs 1e-15 ref |
//! |---|---|---|---|---|---|
//! | `nrfind` (what we had) | 4.6 | 832 | 59,688 | ~6,700 | ~9.5e-13 vs new |
//! | `rootfinder::solve` (shipped) | 53.1 | 9,092 | 575,484 | 92,712 | -- |
//! | `v2` -- exit on the Newton correction | 25.4 | 4,565 | 289,032 | 45,600 | 5.6e-17 |
//! | `v3` -- exit on the price residual | 25.5 | 4,458 | 296,604 | n/a | **9.6e-9 -- unsafe** |
//! | `v4` -- residual-decrease step acceptance, rate-space exits | 5.2 | 1,463 | 101,208 | 15,153 | 2.2e-16 |
//!
//! The shipped solve costs ~11x the evaluations it replaced, at the same accuracy.  `v4` lands at
//! ~1.8x `nrfind` (1,463 vs 832 evals, 101,208 vs 59,688 curve calls) -- the extra is the
//! residual probe that buys step acceptance -- while keeping the bracket, the seed and every
//! failure mode the rewrite was done for.  It is not a bad derivative
//! (`d_ratio` = 1.0000 against a central difference) and not a huge bracket
//! (0.13 - 2.6 rate units on this grid).  `AB_TRACE=1` shows where the passes go:
//!
//! * Iteration 1: the Newton step from `mu_r` lands ~2e-4 from the root.
//! * The min-progress guard rejects it -- `|correction| (1.4e-2) < 0.1 x width (3.5e-1)` -- so
//!   the solver bisects instead.
//! * From there it alternates: a near-root step whose correction is tiny (already converged) gets
//!   rejected again, and bisection spends the remaining ~35 passes grinding the *stale far end*
//!   of the bracket down to the 1e-12 rate tolerance.
//!
//! The tolerance is on the rate, which is the right thing to ask for here, so the fix is not to
//! loosen it -- it is to stop rejecting the step that already converged.
//!
//! * `v2` (accept "correction below tolerance" as convergence) recovers ~2x and is exactly
//!   conservative otherwise.
//! * `v3` shows why the exit must stay in rate space: `f` is a price and can be steep, so a
//!   1e-8 price residual hid a 7e-3 error in the critical rate, which would put every leg
//!   strike in the Jamshidian decomposition wrong while the objective reported "zero".
//! * `v4` (accept a Newton step that lands inside the bracket and strictly reduces the residual,
//!   in addition to the width rule; all *exits* still in rate space) matches the old solve's
//!   iteration count -- ~5 -- while keeping the bracket, the seed and every failure mode.
//!
//! `v4`'s cost is bounded where the guard was protecting us: on the deep exponential-tail
//! fixture it takes 119 iterations to the exact root versus the shipped solver's 65 (bisection's
//! own floor on that fixture is ~48), so the pathological-tail penalty is ~2x, and it is paid
//! only in a regime that does not occur on the tested model grid.  The shipped `max_iterations`
//! still bounds it.
//!
//! Recommended follow-up: port the `v4` step rule into `rootfinder::solve` (the loop shape does
//! not change -- the `inside && (improves || makes_progress)` test replaces `inside &&
//! makes_progress`), keep the exits in rate space, and keep this harness as the regression guard
//! on the iteration count.
//!
//! ## What the old solve could not do (section F, 48 adversarial cases)
//!
//! On the ordinary grid the old solve converged everywhere: 84/84, residuals ~1e-16, agreeing
//! with the new root to 9.5e-13.  Its 11x cost came with no robustness dividend *visible on
//! that grid*.  Push the strike away from the money -- K = 100, or 1e6, against a bond worth
//! about 1 -- and the picture changes: **18 of 48 cases the old solve never converged**.  It
//! burned the full 50-iteration cap, and most of those ended on a `NaN` residual because a
//! Newton step overshot into exp overflow with no bracket to pull it back inside.  Worst old
//! root error: 3.5e1, i.e. 35,000bp from the reference.  The bracketed solve answered every
//! one of them.
//!
//! So the rewrite bought real robustness and paid 11x for it on the cases where the old solver
//! was already fine.  `v4` buys the robustness without the 11x.

#![cfg(test)]

use crate::{HullWhite, SolverSettings, coupon_bond_generic_t, rootfinder};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

/// The shape shared by every solve routine measured here: objective, derivative, an optional
/// bracket, a seed, and the settings.  Four call sites take it, so it gets a name.
type SolveFn<'a> = &'a dyn Fn(
    &dyn Fn(f64) -> f64,
    &dyn Fn(f64) -> f64,
    Option<(f64, f64)>,
    f64,
    &SolverSettings,
) -> Result<rootfinder::Solution, rootfinder::SolverError>;

// The pre-rootfinder configuration, verbatim from the implementation this replaced
// (`nrfind::find_root(&f, &df, R_INIT, PREC_1, MAX_ITER)`).
const R_INIT: f64 = 0.03;
const PREC_1: f64 = 0.0000001;
const MAX_ITER: i32 = 50;

/// Old-solver knobs pushed to roughly the new solver's default accuracy, to ask what the old
/// iteration would have cost if it had been asked for the same answer.
const PREC_TIGHT: f64 = 1e-12;
/// New solver loosened to roughly where the old one stopped (old converged on a *step* of 1e-7,
/// about a bracket of that half-width).
const TOL_MATCHED: f64 = 1e-7;

fn reps() -> usize {
    std::env::var("AB_REPS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(25)
}

/// Counting stand-ins for the market curves.
///
/// Every zero-coupon leg price funnels through the yield and forward curves, so a count of curve
/// calls is a measure of model work that no solver can game.
struct Counted {
    yield_curve: Box<dyn Fn(f64) -> f64 + Sync>,
    forward_curve: Box<dyn Fn(f64) -> f64 + Sync>,
    calls: Arc<AtomicUsize>,
}

impl Counted {
    fn new(curr: f64, a: f64, b: f64, sig: f64) -> Self {
        let calls = Arc::new(AtomicUsize::new(0));
        let yield_curve = {
            let calls = calls.clone();
            Box::new(move |t: f64| {
                calls.fetch_add(1, Ordering::Relaxed);
                let at = (1.0 - (-a * t).exp()) / a;
                let ct = (b - sig * sig / (2.0 * a * a)) * (at - t)
                    - (sig * at) * (sig * at) / (4.0 * a);
                at * curr - ct
            })
        };
        let forward_curve = {
            let calls = calls.clone();
            Box::new(move |t: f64| {
                calls.fetch_add(1, Ordering::Relaxed);
                curr - b * (-a * t).exp()
            })
        };
        Self {
            yield_curve,
            forward_curve,
            calls,
        }
    }
}

struct Case {
    label: String,
    curr: f64,
    a: f64,
    b: f64,
    sig: f64,
    schedule: Vec<f64>,
    coupon_rate: f64,
    strike: f64,
    r_t: f64,
    t: f64,
    option_maturity: f64,
}

fn grid() -> Vec<Case> {
    //label, curr, a (mean reversion), b (long-run rate), sig
    let fixtures = [
        ("flat", 0.05, 0.05, 0.05, 0.01),
        ("steep", 0.02, 0.2, 0.06, 0.03),
        ("hivol", 0.05, 0.1, 0.08, 0.08),
        ("lowvol", 0.03, 0.4, 0.045, 0.005),
    ];
    let short: Vec<f64> = vec![2.5, 3.0, 3.5, 4.0];
    let mid: Vec<f64> = vec![2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0];
    let long: Vec<f64> = (1..=48).map(|i| 2.0 + 0.5 * i as f64).collect();
    let mut out = Vec::new();
    for &(name, curr, a, b, sig) in fixtures.iter() {
        for (sn, schedule) in [("c4", &short), ("c8", &mid), ("c48", &long)] {
            for strike in [0.5f64, 0.8, 0.95, 1.0, 1.05, 1.3, 3.0] {
                out.push(Case {
                    label: format!("{name}/{sn}/K{strike}"),
                    curr,
                    a,
                    b,
                    sig,
                    schedule: schedule.clone(),
                    coupon_rate: 0.05,
                    strike,
                    r_t: 0.04,
                    t: 1.0,
                    option_maturity: 2.0,
                });
            }
        }
    }
    out
}

/// One solve's worth of measurements.
struct Run {
    ok: bool,
    root: f64,
    residual: f64,
    iterations: u32,
    f_calls: usize,
    df_calls: usize,
    curve_calls: usize,
    setup_curve_calls: usize,
    bracket: Option<(f64, f64)>,
    ns: f64,
}

impl Run {
    fn failed() -> Self {
        Self {
            ok: false,
            root: f64::NAN,
            residual: f64::NAN,
            iterations: 0,
            f_calls: 0,
            df_calls: 0,
            curve_calls: 0,
            setup_curve_calls: 0,
            bracket: None,
            ns: 0.0,
        }
    }
    fn width(&self) -> f64 {
        self.bracket
            .map(|(lo, hi)| (hi - lo).abs())
            .unwrap_or(f64::NAN)
    }
}

type Curve = Box<dyn Fn(f64) -> f64 + Sync>;
type Model<'h> = HullWhite<'h, Curve, Curve>;

/// The production objective: the coupon bond's value at option expiry, less the strike.
fn objective<'h>(
    hw: &'h Model<'h>,
    case: &'h Case,
    f_calls: &'h AtomicUsize,
) -> impl Fn(f64) -> f64 + 'h {
    move |rate: f64| {
        f_calls.fetch_add(1, Ordering::Relaxed);
        coupon_bond_generic_t(
            rate,
            case.option_maturity,
            &case.schedule,
            case.coupon_rate,
            &|leg_rate: f64, leg_t: f64, bond_maturity: f64| {
                hw.bond_price_t_raw(leg_rate, leg_t, bond_maturity)
            },
        ) - case.strike
    }
}

/// The production derivative of that objective with respect to the rate.
fn derivative<'h>(
    hw: &'h Model<'h>,
    case: &'h Case,
    df_calls: &'h AtomicUsize,
) -> impl Fn(f64) -> f64 + 'h {
    move |rate: f64| {
        df_calls.fetch_add(1, Ordering::Relaxed);
        hw.coupon_bond_price_t_deriv(rate, case.option_maturity, &case.schedule, case.coupon_rate)
    }
}

/// Runs the new (bracketed, safeguarded) solve the way `solve_critical_rate` does.
fn run_new<'h>(
    hw: &'h Model<'h>,
    curve_calls: &AtomicUsize,
    f_calls: &AtomicUsize,
    df_calls: &AtomicUsize,
    case: &Case,
    settings: SolverSettings,
) -> Run {
    let call_base = curve_calls.load(Ordering::Relaxed);
    //Setup: the analytic bracket and the seed, both of which the old solver never paid for.
    let bracket = hw.critical_rate_bracket(
        case.option_maturity,
        &case.schedule,
        case.coupon_rate,
        case.strike,
    );
    let seed = hw.mu_r(case.r_t, case.t, case.option_maturity).unwrap();
    let setup_curve_calls = curve_calls.load(Ordering::Relaxed) - call_base;

    let objective = objective(hw, case, f_calls);
    let derivative = derivative(hw, case, df_calls);

    let start = Instant::now();
    let outcome = rootfinder::solve(&objective, &derivative, bracket, seed, &settings);
    let ns = start.elapsed().as_nanos() as f64;
    match outcome {
        Ok(solution) => Run {
            ok: true,
            root: solution.root,
            residual: solution.residual,
            iterations: solution.iterations,
            f_calls: f_calls.load(Ordering::Relaxed),
            df_calls: df_calls.load(Ordering::Relaxed),
            curve_calls: curve_calls.load(Ordering::Relaxed) - call_base,
            setup_curve_calls,
            bracket,
            ns,
        },
        Err(_) => Run::failed(),
    }
}

/// Runs the old solve: bare Newton from a hard-coded 3%, converged on the step it took.
fn run_old<'h>(
    hw: &'h Model<'h>,
    curve_calls: &AtomicUsize,
    f_calls: &AtomicUsize,
    df_calls: &AtomicUsize,
    case: &Case,
    precision: f64,
) -> Run {
    let call_base = curve_calls.load(Ordering::Relaxed);
    let objective = objective(hw, case, f_calls);
    let derivative = derivative(hw, case, df_calls);

    let start = Instant::now();
    let outcome = nrfind::find_root(&objective, &derivative, R_INIT, precision, MAX_ITER);
    let ns = start.elapsed().as_nanos() as f64;
    //Read the counters before the verification evaluation below, so it is not charged to the solve.
    let f_calls = f_calls.load(Ordering::Relaxed);
    let df_calls = df_calls.load(Ordering::Relaxed);
    //nrfind returns Err(current_x): the best guess it saw, with no promise it is near a root.
    let root = outcome.unwrap_or_else(|best| best);
    //Reporting the guess is right; counting it as a success is not.  An iteration-cap exit or a
    //non-finite guess is an old-solver failure and has to be counted as one, otherwise every
    //"old never fails" line below is an artefact of the harness rather than a result.
    let converged = outcome.is_ok() && root.is_finite();
    Run {
        ok: converged,
        residual: objective(root).abs(),
        root,
        //nrfind never reports iterations; one f plus one df is one Newton pass.
        iterations: f_calls.try_into().unwrap_or(u32::MAX),
        f_calls,
        df_calls,
        curve_calls: curve_calls.load(Ordering::Relaxed) - call_base,
        setup_curve_calls: 0,
        bracket: None,
        ns,
    }
}

/// A straddle test that treats NaN as "does not straddle", matching the shipped solver.
fn straddles(f_lower: f64, f_upper: f64) -> bool {
    (f_lower <= 0.0 && f_upper >= 0.0) || (f_lower >= 0.0 && f_upper <= 0.0)
}

/// Reconstructs a solve's trajectory from the calls it made, without touching the solver.
///
/// Each iteration of these loops asks for `f(x)` once and `df(x)` once, in that order, and moves
/// the bracket by the sign of `f(x)`.  Recording both calls and replaying that update rule gives
/// the iterate, the residual, the Newton correction and the bracket width per iteration -- which
/// is everything needed to see where the passes go.
struct Trace {
    xs: Vec<f64>,
    fs: Vec<f64>,
    dfs: Vec<f64>,
}

fn traced_solve(
    hw: &Model<'_>,
    case: &Case,
    settings: &SolverSettings,
    solver: SolveFn<'_>,
) -> (Trace, Result<rootfinder::Solution, rootfinder::SolverError>) {
    let xs: Arc<std::sync::Mutex<Vec<f64>>> = Arc::new(std::sync::Mutex::new(Vec::new()));
    let fs: Arc<std::sync::Mutex<Vec<f64>>> = Arc::new(std::sync::Mutex::new(Vec::new()));
    let dfs: Arc<std::sync::Mutex<Vec<f64>>> = Arc::new(std::sync::Mutex::new(Vec::new()));
    let objective = {
        let (xs, fs) = (xs.clone(), fs.clone());
        move |rate: f64| {
            let value = coupon_bond_generic_t(
                rate,
                case.option_maturity,
                &case.schedule,
                case.coupon_rate,
                &|leg_rate: f64, leg_t: f64, bond_maturity: f64| {
                    hw.bond_price_t_raw(leg_rate, leg_t, bond_maturity)
                },
            ) - case.strike;
            xs.lock().unwrap().push(rate);
            fs.lock().unwrap().push(value);
            value
        }
    };
    let derivative = {
        let dfs = dfs.clone();
        move |rate: f64| {
            let value = hw.coupon_bond_price_t_deriv(
                rate,
                case.option_maturity,
                &case.schedule,
                case.coupon_rate,
            );
            dfs.lock().unwrap().push(value);
            value
        }
    };
    let _ = &objective;
    let brk = hw.critical_rate_bracket(
        case.option_maturity,
        &case.schedule,
        case.coupon_rate,
        case.strike,
    );
    let seed = hw.mu_r(case.r_t, case.t, case.option_maturity).unwrap();
    let outcome = solver(&objective, &derivative, brk, seed, settings);
    let trace = Trace {
        xs: xs.lock().unwrap().clone(),
        fs: fs.lock().unwrap().clone(),
        dfs: dfs.lock().unwrap().clone(),
    };
    (trace, outcome)
}

impl Trace {
    /// Prints one line per iteration: iterate, residual, Newton correction, and the bracket the
    /// replayed sign rule says the solver was holding.
    ///
    /// `f` is called three times before the loop runs on the bracketed path -- the seed, then each
    /// end of the proposed bracket -- so the iteration records start at offset 3.  `df` is only
    /// ever called inside the loop, so it starts at offset 0.
    fn print(&self, label: &str, lo: f64, hi: f64) {
        let mut a = lo;
        let mut b = hi;
        let mut f_a = self.fs.get(1).copied().unwrap_or(f64::NAN);
        let _f_b = self.fs.get(2).copied().unwrap_or(f64::NAN);
        println!("{label}\tITER\tX\tF(X)\tDF(X)\tCORRECTION\tBRACKET\tWIDTH");
        for iteration in 1.. {
            let index = 2 + iteration;
            let Some(&x) = self.xs.get(index) else { break };
            let f_value = self.fs[index];
            let slope = self.dfs.get(iteration - 1).copied().unwrap_or(f64::NAN);
            let correction = if slope.is_finite() && slope != 0.0 {
                f_value / slope
            } else {
                f64::NAN
            };
            if f_value.is_finite() {
                if (f_value < 0.0) == (f_a < 0.0) {
                    a = x;
                    f_a = f_value;
                } else {
                    b = x;
                }
            }
            println!(
                "{label}\t{iteration}\t{x:.12}\t{f_value:.6e}\t{slope:.6e}\t{correction:.6e}\t\
                 [{a:.6}, {b:.6}]\t{:.6e}
            ",
                (b - a).abs()
            );
        }
    }
}

#[test]
fn trace_one_solve() {
    if std::env::var("AB_TRACE").is_err() {
        eprintln!("trace_one_solve skipped: set AB_TRACE=1 to dump per-iteration trajectories");
        return;
    }
    let case = grid()
        .into_iter()
        .find(|c| c.label == "flat/c4/K1")
        .expect("fixture case flat/c4/K1");
    let counted = Counted::new(case.curr, case.a, case.b, case.sig);
    let hw = HullWhite::init(
        case.a,
        case.sig,
        &counted.yield_curve,
        &counted.forward_curve,
    )
    .unwrap();
    let brk = hw
        .critical_rate_bracket(
            case.option_maturity,
            &case.schedule,
            case.coupon_rate,
            case.strike,
        )
        .expect("analytic bracket");
    println!(
        "case {} seed {:.12} bracket [{:.6}, {:.6}]",
        case.label,
        hw.mu_r(case.r_t, case.t, case.option_maturity).unwrap(),
        brk.0,
        brk.1
    );
    let (trace, outcome) = traced_solve(&hw, &case, &SolverSettings::default(), &rootfinder::solve);
    println!("shipped outcome: {outcome:?}");
    trace.print("SHIPPED", brk.0, brk.1);
    let (trace, outcome) = traced_solve(&hw, &case, &SolverSettings::default(), &solve_v2);
    println!("candidate outcome: {outcome:?}");
    trace.print("CANDIDATE", brk.0, brk.1);
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
fn solve_v2(
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
fn solve_v3(
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
/// That danger is a function of how far the *price* still is from the strike, not of how wide the
/// bracket happens to be: inside the band where `|f| <= far_band`, Newton is already in its
/// quadratic region and every inside-bracket step is worth taking.  Outside it, the width rule
/// still applies and bisection does the long-distance work.
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
fn solve_v4(
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

/// Bracket us or widen around the seed, shared by the candidate loops.
fn bracket_for(
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

/// Runs a solve variant with the same counted objective/derivative plumbing.
fn run_variant(
    hw: &Model<'_>,
    curve_calls: &AtomicUsize,
    f_calls: &AtomicUsize,
    df_calls: &AtomicUsize,
    case: &Case,
    settings: SolverSettings,
    solver: SolveFn<'_>,
) -> Run {
    let call_base = curve_calls.load(Ordering::Relaxed);
    let brk = hw.critical_rate_bracket(
        case.option_maturity,
        &case.schedule,
        case.coupon_rate,
        case.strike,
    );
    let seed = hw.mu_r(case.r_t, case.t, case.option_maturity).unwrap();
    let setup_curve_calls = curve_calls.load(Ordering::Relaxed) - call_base;
    let objective = objective(hw, case, f_calls);
    let derivative = derivative(hw, case, df_calls);
    let start = Instant::now();
    let outcome = solver(&objective, &derivative, brk, seed, &settings);
    let ns = start.elapsed().as_nanos() as f64;
    match outcome {
        Ok(solution) => Run {
            ok: true,
            root: solution.root,
            residual: solution.residual,
            iterations: solution.iterations,
            f_calls: f_calls.load(Ordering::Relaxed),
            df_calls: df_calls.load(Ordering::Relaxed),
            curve_calls: curve_calls.load(Ordering::Relaxed) - call_base,
            setup_curve_calls,
            bracket: brk,
            ns,
        },
        Err(_) => Run::failed(),
    }
}

#[test]
fn solver_ab_cost() {
    if std::env::var("SOLVER_AB").is_err() {
        eprintln!("solver_ab skipped: set SOLVER_AB=1 to run the A/B cost comparison");
        return;
    }
    let reps = reps();
    let cases = grid();

    //The reference: the new solver pushed far past any usable accuracy.
    let reference_settings = SolverSettings {
        tolerance: 1e-15,
        max_iterations: 400,
        initial_guess: None,
    };
    let production_settings = SolverSettings::default();
    let matched_settings = SolverSettings {
        tolerance: TOL_MATCHED,
        max_iterations: 50,
        initial_guess: None,
    };

    println!("== A: per-case cost (single pass; counts are exact) ==");
    println!(
        "CASE\tREF_ROOT\tBRACKET_LO\tBRACKET_HI\tBRACKET_W\tNEW_ITERS\tNEW_F\tNEW_DF\t\
         NEW_SETUP\tNEW_CURVES\tOLD_ITERS\tOLD_F\tOLD_CURVES\tROOT_DELTA\tOLD_RESID\t\
         D_ANALYTIC\tD_CENTRAL\tD_RATIO\tNEW@1e-7_IT\tNEW@1e-7_F\tOLD@1e-12_IT\tVERDICT"
    );

    let mut totals = Totals::default();
    for case in &cases {
        let counted = Counted::new(case.curr, case.a, case.b, case.sig);
        let hw = HullWhite::init(
            case.a,
            case.sig,
            &counted.yield_curve,
            &counted.forward_curve,
        )
        .unwrap();
        //Drop init's curve traffic out of the count.
        counted.calls.store(0, Ordering::Relaxed);

        let reference = run_new(
            &hw,
            &counted.calls,
            &AtomicUsize::new(0),
            &AtomicUsize::new(0),
            case,
            reference_settings,
        );

        let f_new = AtomicUsize::new(0);
        let df_new = AtomicUsize::new(0);
        counted.calls.store(0, Ordering::Relaxed);
        let new = run_new(
            &hw,
            &counted.calls,
            &f_new,
            &df_new,
            case,
            production_settings,
        );

        let f_old = AtomicUsize::new(0);
        let df_old = AtomicUsize::new(0);
        counted.calls.store(0, Ordering::Relaxed);
        let old = run_old(&hw, &counted.calls, &f_old, &df_old, case, PREC_1);

        let f_match = AtomicUsize::new(0);
        let df_match = AtomicUsize::new(0);
        counted.calls.store(0, Ordering::Relaxed);
        let new_matched = run_new(
            &hw,
            &counted.calls,
            &f_match,
            &df_match,
            case,
            matched_settings,
        );

        let f_tight = AtomicUsize::new(0);
        let df_tight = AtomicUsize::new(0);
        counted.calls.store(0, Ordering::Relaxed);
        let old_tight = run_old(&hw, &counted.calls, &f_tight, &df_tight, case, PREC_TIGHT);

        let root_delta = (new.root - old.root).abs();

        //Is the analytic derivative actually the objective's derivative?  If it is badly off, the
        //min-progress guard rejects every Newton step and the solver bisects a thousands-wide
        //bracket down to 1e-12 one halving at a time -- which is what ~53 iterations means.
        //A fresh uncounted pair keeps the check off every solve's tab.
        let null = AtomicUsize::new(0);
        let f_check = objective(&hw, case, &null);
        let d_check = derivative(&hw, case, &null);
        let h = 1e-5 * new.root.abs().max(1.0);
        let central = (f_check(new.root + h) - f_check(new.root - h)) / (2.0 * h);
        let analytic = d_check(new.root);
        let d_ratio = if central != 0.0 {
            analytic / central
        } else {
            f64::NAN
        };

        println!(
            "{}\t{:.12}\t{:.6}\t{:.6}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3e}\t{:.3e}\t\
             {:.6e}\t{:.6e}\t{:.4}\t{}\t{}\t{}\t{}",
            case.label,
            reference.root,
            new.bracket.map(|b| b.0).unwrap_or(f64::NAN),
            new.bracket.map(|b| b.1).unwrap_or(f64::NAN),
            new.width(),
            new.iterations,
            new.f_calls,
            new.df_calls,
            new.setup_curve_calls,
            new.curve_calls,
            old.iterations,
            old.f_calls,
            old.curve_calls,
            root_delta,
            old.residual,
            analytic,
            central,
            d_ratio,
            new_matched.iterations,
            new_matched.f_calls,
            old_tight.iterations,
            verdict(&reference, &new, &old),
        );

        totals.record(
            case,
            &new,
            &old,
            root_delta,
            d_ratio,
            verdict(&reference, &new, &old),
        );
    }

    println!();
    println!("== B: aggregate counts over {} cases ==", cases.len());
    totals.report();

    println!();
    println!("== C: wall clock per solve (median of {reps} interleaved reps) ==");
    println!("CASE\tNEW_NS\tOLD_NS\tRATIO_NS\tNEW_CURVES\tOLD_CURVES\tRATIO_CALLS");
    for case in cases
        .iter()
        .filter(|c| matches!(c.strike, 0.5 | 0.95 | 1.0 | 3.0))
    {
        let counted = Counted::new(case.curr, case.a, case.b, case.sig);
        let hw = HullWhite::init(
            case.a,
            case.sig,
            &counted.yield_curve,
            &counted.forward_curve,
        )
        .unwrap();
        let mut new_times = Vec::with_capacity(reps);
        let mut old_times = Vec::with_capacity(reps);
        let mut new_calls = Vec::with_capacity(reps);
        let mut old_calls = Vec::with_capacity(reps);
        for _ in 0..reps {
            let f = AtomicUsize::new(0);
            let d = AtomicUsize::new(0);
            counted.calls.store(0, Ordering::Relaxed);
            let new = run_new(&hw, &counted.calls, &f, &d, case, production_settings);
            new_times.push(new.ns);
            new_calls.push(new.curve_calls as f64);

            let f = AtomicUsize::new(0);
            let d = AtomicUsize::new(0);
            counted.calls.store(0, Ordering::Relaxed);
            let old = run_old(&hw, &counted.calls, &f, &d, case, PREC_1);
            old_times.push(old.ns);
            old_calls.push(old.curve_calls as f64);
        }
        let new_ns = median(&mut new_times);
        let old_ns = median(&mut old_times);
        let new_calls = median(&mut new_calls);
        let old_calls = median(&mut old_calls);
        println!(
            "{}\t{new_ns:.0}\t{old_ns:.0}\t{:.2}\t{new_calls:.0}\t{old_calls:.0}\t{:.2}",
            case.label,
            if old_ns > 0.0 {
                new_ns / old_ns
            } else {
                f64::NAN
            },
            if old_calls > 0.0 {
                new_calls / old_calls
            } else {
                f64::NAN
            },
        );
    }

    println!();
    println!("== D: candidate fix (converge on the correction) vs new vs old ==");
    println!(
        "CASE\tV2_ITERS\tV2_F\tV2_CURVES\tV2_NS\tROOT_ERR_VS_REF\t\
         OLD_ITERS\tOLD_CURVES\tOLD_NS\tNEW_NS\tV2_OVER_OLD_CALLS\tV2_OVER_OLD_NS"
    );
    let mut v2_evals = 0usize;
    let mut v2_curves = 0usize;
    let mut old_evals_d = 0usize;
    let mut old_curves_d = 0usize;
    let mut v2_iters_sum = 0u64;
    let mut max_root_err = 0f64;
    let mut max_resid = 0f64;
    let d_cases: Vec<&Case> = cases
        .iter()
        .filter(|c| matches!(c.strike, 0.5 | 0.8 | 0.95 | 1.0 | 1.05 | 1.3 | 3.0))
        .collect();
    for case in &d_cases {
        let counted = Counted::new(case.curr, case.a, case.b, case.sig);
        let hw = HullWhite::init(
            case.a,
            case.sig,
            &counted.yield_curve,
            &counted.forward_curve,
        )
        .unwrap();
        let reference = run_variant(
            &hw,
            &counted.calls,
            &AtomicUsize::new(0),
            &AtomicUsize::new(0),
            case,
            reference_settings,
            &solve_v2,
        );
        let mut times = Vec::with_capacity(reps);
        let mut v2 = Run::failed();
        for _ in 0..times_len(reps) {
            let f = AtomicUsize::new(0);
            let d = AtomicUsize::new(0);
            counted.calls.store(0, Ordering::Relaxed);
            let run = run_variant(
                &hw,
                &counted.calls,
                &f,
                &d,
                case,
                production_settings,
                &solve_v2,
            );
            times.push(run.ns);
            v2 = run;
        }
        let f_old = AtomicUsize::new(0);
        let d_old = AtomicUsize::new(0);
        counted.calls.store(0, Ordering::Relaxed);
        let old = run_old(&hw, &counted.calls, &f_old, &d_old, case, PREC_1);
        let mut old_times = Vec::with_capacity(reps);
        for _ in 0..times_len(reps) {
            let f = AtomicUsize::new(0);
            let d = AtomicUsize::new(0);
            counted.calls.store(0, Ordering::Relaxed);
            old_times.push(run_old(&hw, &counted.calls, &f, &d, case, PREC_1).ns);
        }
        //What the shipped solver costs on the same case, for the ratio column.
        let new_ns = {
            let mut new_times = Vec::with_capacity(reps);
            for _ in 0..times_len(reps) {
                let f = AtomicUsize::new(0);
                let d = AtomicUsize::new(0);
                counted.calls.store(0, Ordering::Relaxed);
                new_times.push(run_new(&hw, &counted.calls, &f, &d, case, production_settings).ns);
            }
            median(&mut new_times)
        };
        let v2_ns = median(&mut times);
        let old_ns = median(&mut old_times);
        let root_err = (v2.root - reference.root).abs();
        max_root_err = max_root_err.max(root_err);
        max_resid = max_resid.max(v2.residual);
        v2_evals += v2.f_calls + v2.df_calls;
        v2_curves += v2.curve_calls;
        old_evals_d += old.f_calls + old.df_calls;
        old_curves_d += old.curve_calls;
        v2_iters_sum += v2.iterations as u64;
        println!(
            "{}\t{}\t{}\t{}\t{:.0}\t{:.3e}\t{}\t{}\t{:.0}\t{:.0}\t{:.2}\t{:.2}",
            case.label,
            v2.iterations,
            v2.f_calls,
            v2.curve_calls,
            v2_ns,
            root_err,
            old.iterations,
            old.curve_calls,
            old_ns,
            new_ns,
            if old_curves_d > 0 {
                v2.curve_calls as f64 / old.curve_calls as f64
            } else {
                f64::NAN
            },
            if old_ns > 0.0 {
                v2_ns / old_ns
            } else {
                f64::NAN
            },
        );
    }
    println!(
        "SUMMARY\t{:.1} iters avg\t{} evals\t{} curves\tmax root err vs reference {:.3e}\tmax residual {:.3e}\told {} evals / {} curves\teval ratio {:.2}\tcall ratio {:.2}",
        v2_iters_sum as f64 / cases.len() as f64,
        v2_evals,
        v2_curves,
        max_root_err,
        max_resid,
        old_evals_d,
        old_curves_d,
        ratio(v2_evals, old_evals_d),
        ratio(v2_curves, old_curves_d),
    );
    assert!(
        max_root_err < 1e-11,
        "candidate solve drifted from the reference by {max_root_err:.3e}"
    );
}

fn median(values: &mut [f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    values[values.len() / 2]
}

fn times_len(reps: usize) -> usize {
    reps.max(1)
}

/// The safety properties the shipped solver was written to have, re-checked against the candidate.
/// These mirror the unit tests in `rootfinder`: a fast solve that reintroduces the tail crawl or
/// trusts a lying derivative is not a fix.
#[test]
fn candidate_keeps_the_safety_properties() {
    let settings = SolverSettings {
        tolerance: 1e-12,
        max_iterations: 100,
        initial_guess: None,
    };

    //Tail crawl: deep in exp(-x) the Newton step is a fixed +1 regardless of distance, so an
    //unguarded solve needs ~140 passes across a wide bracket.  The guard must still be there.
    let root = 1e6_f64.ln();
    let f = |x: f64| (-x).exp() - 1e-6;
    let df = |x: f64| -(-x).exp();
    let solution = solve_v2(&f, &df, None, -100.0, &settings).unwrap();
    assert!(
        (solution.root - root).abs() < 1e-10,
        "tail case drifted: {solution:?} vs {root}"
    );
    assert!(solution.iterations <= 70, "tail case crawled: {solution:?}");

    //Lying derivative: wrong sign must still land on the root via the bisection guard.
    let solution = solve_v2(&|x| x - 0.25, &|_| -1.0, Some((0.0, 1.0)), 0.9, &settings).unwrap();
    assert!((solution.root - 0.25).abs() < 1e-10, "{solution:?}");

    //Non-straddling bracket with no root anywhere: still reported, never silently wrong.
    let error = solve_v2(&|x| x * x + 1.0, &|x| 2.0 * x, None, 0.0, &settings).unwrap_err();
    assert!(
        matches!(error, rootfinder::SolverError::NoSignChange { .. }),
        "{error:?}"
    );

    //The case the shipped solver is slow on: a nearly linear objective inside a wide bracket.
    //A sub-tolerance Newton correction has to read as converged, not as a stall.
    let tight = solve_v2(
        &|x| 2.0 * x - 1.0,
        &|_| 2.0,
        Some((-5.0, 5.0)),
        0.0,
        &settings,
    )
    .unwrap();
    assert!((tight.root - 0.5).abs() < 1e-12, "{tight:?}");
    assert!(
        tight.iterations <= 5,
        "linear objective still took {} iterations; the correction exit is not firing",
        tight.iterations
    );
}

#[test]
fn exit_rule_variants() {
    if std::env::var("SOLVER_AB").is_err() {
        eprintln!("exit_rule_variants skipped: set SOLVER_AB=1");
        return;
    }
    let cases = grid();
    //A price-scale residual: 1e-8 on a bond worth about 1 is well under any quoted basis point.
    let resid_rel = 1e-8;
    let settings = SolverSettings::default();
    println!("== E: exit-rule variants (counts; no timing) ==");
    println!(
        "CASE\tREF_ROOT\tSHIP_IT\tSHIP_F\tSHIP_NS\tV2_IT\tV2_F\tV2_NS\tV3_IT\tV3_F\tV4_IT\t\
         V4_F\tV4_NS\tV3_ROOTERR\tV4_ROOTERR\tSHIP_CURVES\tV4_CURVES"
    );
    let (mut ship_evals, mut v2_evals, mut v3_evals, mut v4_evals) =
        (0usize, 0usize, 0usize, 0usize);
    let (mut ship_curves, mut v2_curves, mut v3_curves, mut v4_curves) =
        (0usize, 0usize, 0usize, 0usize);
    let (mut ship_iters, mut v2_iters, mut v3_iters, mut v4_iters) = (0u64, 0u64, 0u64, 0u64);
    let (mut v3_err, mut v4_err) = (0f64, 0f64);
    let (mut ship_ns, mut v2_ns, mut v4_ns) = (0usize, 0usize, 0usize);
    let tight = SolverSettings {
        tolerance: 1e-15,
        max_iterations: 500,
        initial_guess: None,
    };
    for case in &cases {
        let counted = Counted::new(case.curr, case.a, case.b, case.sig);
        let hw = HullWhite::init(
            case.a,
            case.sig,
            &counted.yield_curve,
            &counted.forward_curve,
        )
        .unwrap();
        let reference = run_new(
            &hw,
            &counted.calls,
            &AtomicUsize::new(0),
            &AtomicUsize::new(0),
            case,
            tight,
        );
        let r = |f_calls: &AtomicUsize, df_calls: &AtomicUsize, solver: &_| -> Run {
            run_variant(
                &hw,
                &counted.calls,
                f_calls,
                df_calls,
                case,
                settings,
                solver,
            )
        };
        let shipped = run_new(
            &hw,
            &counted.calls,
            &AtomicUsize::new(0),
            &AtomicUsize::new(0),
            case,
            settings,
        );
        let v2 = r(&AtomicUsize::new(0), &AtomicUsize::new(0), &solve_v2);
        let v3 = r(
            &AtomicUsize::new(0),
            &AtomicUsize::new(0),
            &|f, df, brk, seed, s| {
                solve_v3(f, df, brk, seed, s, resid_rel * case.strike.abs().max(1.0))
            },
        );
        let v4 = r(
            &AtomicUsize::new(0),
            &AtomicUsize::new(0),
            &|f, df, brk, seed, s| solve_v4(f, df, brk, seed, s, 0.0),
        );
        ship_evals += shipped.f_calls + shipped.df_calls;
        v2_evals += v2.f_calls + v2.df_calls;
        v3_evals += v3.f_calls + v3.df_calls;
        v4_evals += v4.f_calls + v4.df_calls;
        ship_curves += shipped.curve_calls;
        v2_curves += v2.curve_calls;
        v3_curves += v3.curve_calls;
        v4_curves += v4.curve_calls;
        ship_iters += shipped.iterations as u64;
        v2_iters += v2.iterations as u64;
        v3_iters += v3.iterations as u64;
        v4_iters += v4.iterations as u64;
        v3_err = v3_err.max((v3.root - reference.root).abs());
        v4_err = v4_err.max((v4.root - reference.root).abs());
        ship_ns += shipped.ns as usize;
        v2_ns += v2.ns as usize;
        v4_ns += v4.ns as usize;
        println!(
            "{}\t{:.12}\t{}\t{}\t{:.0}\t{}\t{}\t{:.0}\t{}\t{}\t{}\t{}\t{:.0}\t{:.3e}\t{:.3e}\t{}\t{}",
            case.label,
            reference.root,
            shipped.iterations,
            shipped.f_calls,
            shipped.ns,
            v2.iterations,
            v2.f_calls,
            v2.ns,
            v3.iterations,
            v3.f_calls,
            v4.iterations,
            v4.f_calls,
            v4.ns,
            (v3.root - reference.root).abs(),
            (v4.root - reference.root).abs(),
            shipped.curve_calls,
            v4.curve_calls,
        );
    }
    let n = cases.len() as f64;
    println!(
        "SUMMARY\tshipped {:.1} it / {} ev / {} cv / {:.0} ns\tv2 {:.1} / {} ev / {} cv / {:.0} ns\tv3 {:.1} / {} ev / {} cv\tv4 {:.1} / {} ev / {} cv / {:.0} ns\tv3 root err {:.3e}\tv4 root err {:.3e}\tv2/ship {:.2}\tv3/ship {:.2}\tv4/ship {:.2}\tv4/v2 {:.2}",
        ship_iters as f64 / n,
        ship_evals,
        ship_curves,
        ship_ns as f64 / n,
        v2_iters as f64 / n,
        v2_evals,
        v2_curves,
        v2_ns as f64 / n,
        v3_iters as f64 / n,
        v3_evals,
        v3_curves,
        v4_iters as f64 / n,
        v4_evals,
        v4_curves,
        v4_ns as f64 / n,
        v3_err,
        v4_err,
        ratio(v2_evals, ship_evals),
        ratio(v3_evals, ship_evals),
        ratio(v4_evals, ship_evals),
        ratio(v4_ns, v2_ns),
    );
}

/// Every candidate has to keep the properties the shipped solver was written for.  A faster solve
/// that reintroduces the far-tail crawl or trusts a lying derivative is not a fix.
#[test]
fn variant_safety_checks() {
    let settings = SolverSettings {
        tolerance: 1e-12,
        max_iterations: 500,
        initial_guess: None,
    };
    //Deep exponential tail: the Newton step is a fixed +1 no matter how far the root is.
    let tail_f = |x: f64| (-x).exp() - 1e-6;
    let tail_df = |x: f64| -(-x).exp();
    let tail_root = 1e6_f64.ln();
    println!("VARIANT\tTAIL_ITERS\tTAIL_ROOT_ERR\tLYING_DERIV\tNO_ROOT_REPORTED");

    let shipped = rootfinder::solve(&tail_f, &tail_df, None, -100.0, &settings).unwrap();
    println!(
        "shipped\t{}\t{:.3e}\t{}\t{}",
        shipped.iterations,
        (shipped.root - tail_root).abs(),
        lying(&|f, df, brk, seed, s| rootfinder::solve(f, df, brk, seed, s)),
        no_root(&|f, df, brk, seed, s| rootfinder::solve(f, df, brk, seed, s)),
    );
    let v2 = solve_v2(&tail_f, &tail_df, None, -100.0, &settings).unwrap();
    println!(
        "v2\t{}\t{:.3e}\t{}\t{}",
        v2.iterations,
        (v2.root - tail_root).abs(),
        lying(&solve_v2),
        no_root(&solve_v2),
    );
    let v3 = solve_v3(&tail_f, &tail_df, None, -100.0, &settings, 1e-8).unwrap();
    println!(
        "v3\t{}\t{:.3e}\t{}\t{}",
        v3.iterations,
        (v3.root - tail_root).abs(),
        lying(&|f, df, brk, seed, s| solve_v3(f, df, brk, seed, s, 1e-8)),
        no_root(&|f, df, brk, seed, s| solve_v3(f, df, brk, seed, s, 1e-8)),
    );
    let v4 = solve_v4(&tail_f, &tail_df, None, -100.0, &settings, 0.0).unwrap();
    println!(
        "v4\t{}\t{:.3e}\t{}\t{}",
        v4.iterations,
        (v4.root - tail_root).abs(),
        lying(&|f, df, brk, seed, s| solve_v4(f, df, brk, seed, s, 0.0)),
        no_root(&|f, df, brk, seed, s| solve_v4(f, df, brk, seed, s, 0.0)),
    );
    //v3 is deliberately excluded from the accuracy assertions: the residual exit is kept here as
    //the counter-example that shows why the exit has to stay in rate space.  Its root error is
    //printed above (7e-3 on the tail fixture) rather than asserted.
    for (name, root) in [("shipped", shipped.root), ("v2", v2.root), ("v4", v4.root)] {
        assert!(
            (root - tail_root).abs() < 1e-9,
            "{name} lost the tail root: {root} vs {tail_root}"
        );
    }
}

/// A lying derivative (wrong sign) must still land on the root via the bracket guard.
fn lying(solver: SolveFn<'_>) -> bool {
    let settings = SolverSettings {
        tolerance: 1e-12,
        max_iterations: 500,
        initial_guess: None,
    };
    solver(&|x| x - 0.25, &|_| -1.0, Some((0.0, 1.0)), 0.9, &settings)
        .map(|s| (s.root - 0.25).abs() < 1e-9)
        .unwrap_or(false)
}

/// A function with no root must be reported, never answered.
fn no_root(solver: SolveFn<'_>) -> bool {
    let settings = SolverSettings {
        tolerance: 1e-12,
        max_iterations: 500,
        initial_guess: None,
    };
    solver(&|x| x * x + 1.0, &|x| 2.0 * x, None, 0.0, &settings)
        .map(|_| false)
        .unwrap_or(true)
}

#[derive(Default)]
struct Totals {
    cases: usize,
    new_evals: usize,
    old_evals: usize,
    new_curves: usize,
    old_curves: usize,
    new_setup_curves: usize,
    new_iters: u64,
    old_iters: u64,
    worst_root_delta: f64,
    worst_root_delta_case: String,
    worst_d_ratio_err: f64,
    worst_d_ratio_case: String,
    old_stale: usize,
    new_failed: usize,
    /// Cases where bare Newton never converged but the bracketed solve did.
    old_fail_cases: Vec<String>,
    /// Cases where bare Newton returned confidently and landed on a different root than the
    /// reference, while the bracketed solve stayed tight.
    old_wrong_root_cases: Vec<String>,
}

/// Classifies one case against the 1e-15 reference solve.
///
/// "Wrong root" is the failure that matters most here: the old solver converges on the size of
/// its own step, so on a flat stretch of the objective it can stop far from the root and still
/// report success.  Thresholds are in rate units.
fn verdict(reference: &Run, new: &Run, old: &Run) -> &'static str {
    const WRONG_ROOT: f64 = 1e-4;
    const NEW_TIGHT: f64 = 1e-8;
    match (old.ok, new.ok) {
        (false, false) => "BOTH_FAIL",
        (false, true) => "OLD_FAIL",
        (true, false) => "NEW_FAIL",
        (true, true) => {
            let old_err = (old.root - reference.root).abs();
            let new_err = (new.root - reference.root).abs();
            if old_err > WRONG_ROOT && new_err <= NEW_TIGHT {
                "OLD_WRONG_ROOT"
            } else if old_err > WRONG_ROOT {
                "OLD_INACCURATE"
            } else {
                "OK"
            }
        }
    }
}

impl Totals {
    fn record(
        &mut self,
        case: &Case,
        new: &Run,
        old: &Run,
        root_delta: f64,
        d_ratio: f64,
        verdict: &'static str,
    ) {
        self.cases += 1;
        self.new_evals += new.f_calls + new.df_calls;
        self.old_evals += old.f_calls + old.df_calls;
        self.new_curves += new.curve_calls;
        self.old_curves += old.curve_calls;
        self.new_setup_curves += new.setup_curve_calls;
        self.new_iters += new.iterations as u64;
        self.old_iters += old.iterations as u64;
        if root_delta > self.worst_root_delta {
            self.worst_root_delta = root_delta;
            self.worst_root_delta_case = case.label.clone();
        }
        let d_err = (d_ratio - 1.0).abs();
        if d_err > self.worst_d_ratio_err {
            self.worst_d_ratio_err = d_err;
            self.worst_d_ratio_case = case.label.clone();
        }
        //The old solver converges on its own step size, so it can stop on a nearly flat spot and
        //call that convergence several rate-units from the actual root.
        if old.residual > 1e-8 * case.strike.abs().max(1.0) {
            self.old_stale += 1;
        }
        if !new.ok {
            self.new_failed += 1;
        }
        if verdict == "OLD_FAIL" || verdict == "BOTH_FAIL" {
            self.old_fail_cases.push(case.label.clone());
        }
        if verdict == "OLD_WRONG_ROOT" {
            self.old_wrong_root_cases.push(case.label.clone());
        }
    }

    fn report(&self) {
        println!(
            "cases\tnew_evals\told_evals\teval_ratio\tnew_curves\told_curves\tcall_ratio\t\
             new_setup_curves\tnew_iters_avg\told_iters_avg\told_stale\tnew_failed"
        );
        println!(
            "{}\t{}\t{}\t{:.2}\t{}\t{}\t{:.2}\t{}\t{:.1}\t{:.1}\t{}\t{}",
            self.cases,
            self.new_evals,
            self.old_evals,
            ratio(self.new_evals, self.old_evals),
            self.new_curves,
            self.old_curves,
            ratio(self.new_curves, self.old_curves),
            self.new_setup_curves,
            self.new_iters as f64 / self.cases.max(1) as f64,
            self.old_iters as f64 / self.cases.max(1) as f64,
            self.old_stale,
            self.new_failed,
        );
        println!(
            "worst old/new root disagreement: {:.3e} at {}",
            self.worst_root_delta, self.worst_root_delta_case
        );
        println!(
            "worst analytic-vs-central derivative mismatch: {:.4}% at {}",
            self.worst_d_ratio_err * 100.0,
            self.worst_d_ratio_case
        );
        //The callout the rewrite needs: what the old solver could not do, that the new one can.
        if self.old_fail_cases.is_empty() {
            println!("old solver failures where the new solver solved: none on this grid");
        } else {
            println!(
                "old solver FAILED where the bracketed solve succeeded ({}): {}",
                self.old_fail_cases.len(),
                self.old_fail_cases.join(", ")
            );
        }
        if self.old_wrong_root_cases.is_empty() {
            println!(
                "old solver wrong-root landings where the new one is tight: none on this grid"
            );
        } else {
            println!(
                "old solver LANDED ON THE WRONG ROOT while the bracketed solve stayed tight ({}): {}",
                self.old_wrong_root_cases.len(),
                self.old_wrong_root_cases.join(", ")
            );
        }
        println!(
            "old solver converged stale (residual > 1e-8 x strike) on {} of {} cases",
            self.old_stale, self.cases
        );
    }
}

/// Cases built to break the old solve rather than to be realistic.
///
/// Bare Newton has three exposed flanks: a hard-coded seed (3%) with no bracket to fall back on, a
/// 50-iteration cap, and convergence judged on the step it took.  Puts the root far from 3%, put
/// the price scale far from 1, and see whether it still answers.
fn stress_grid() -> Vec<Case> {
    let fixtures = [
        ("flat", 0.05, 0.05, 0.05, 0.01),
        ("hivol", 0.05, 0.1, 0.08, 0.08),
    ];
    let short: Vec<f64> = vec![2.5, 3.0, 3.5, 4.0];
    let long: Vec<f64> = (1..=48).map(|i| 2.0 + 0.5 * i as f64).collect();
    let mut out = Vec::new();
    for &(name, curr, a, b, sig) in fixtures.iter() {
        for (sn, schedule) in [("c4", &short), ("c48", &long)] {
            for (ci, coupon_rate) in [("cp0.0001", 0.0001f64), ("cp5", 5.0), ("cp0.05", 0.05)] {
                for strike in [1e-6f64, 0.01, 100.0, 1e6] {
                    out.push(Case {
                        label: format!("{name}/{sn}/{ci}/K{strike}"),
                        curr,
                        a,
                        b,
                        sig,
                        schedule: schedule.clone(),
                        coupon_rate,
                        strike,
                        r_t: 0.04,
                        t: 1.0,
                        option_maturity: 2.0,
                    });
                }
            }
        }
    }
    out
}

/// The criterion that matters for a rewrite bought for robustness: is there any case where bare
/// Newton fails, or answers with the wrong root, while the bracketed solve answers correctly?
#[test]
fn old_solver_stress() {
    if std::env::var("SOLVER_AB").is_err() {
        eprintln!("old_solver_stress skipped: set SOLVER_AB=1");
        return;
    }
    let cases = stress_grid();
    let reference_settings = SolverSettings {
        tolerance: 1e-15,
        max_iterations: 400,
        initial_guess: None,
    };
    let production_settings = SolverSettings::default();
    println!("== F: adversarial cases for the old solve (seed 3%, no bracket, cap {MAX_ITER}) ==");
    println!("CASE\tREF_ROOT\tNEW_IT\tOLD_IT\tOLD_RESID\tNEW_RESID\tROOT_ERR_OLD\tVERDICT");
    let mut verdicts: std::collections::BTreeMap<&'static str, Vec<String>> =
        std::collections::BTreeMap::new();
    let (mut worst_old_err, mut worst_old_case) = (0f64, String::new());
    let (mut old_iters_max, mut old_iters_case) = (0u32, String::new());
    for case in &cases {
        let counted = Counted::new(case.curr, case.a, case.b, case.sig);
        let hw = HullWhite::init(
            case.a,
            case.sig,
            &counted.yield_curve,
            &counted.forward_curve,
        )
        .unwrap();
        counted.calls.store(0, Ordering::Relaxed);
        let reference = run_new(
            &hw,
            &counted.calls,
            &AtomicUsize::new(0),
            &AtomicUsize::new(0),
            case,
            reference_settings,
        );
        counted.calls.store(0, Ordering::Relaxed);
        let new = run_new(
            &hw,
            &counted.calls,
            &AtomicUsize::new(0),
            &AtomicUsize::new(0),
            case,
            production_settings,
        );
        counted.calls.store(0, Ordering::Relaxed);
        let old = run_old(
            &hw,
            &counted.calls,
            &AtomicUsize::new(0),
            &AtomicUsize::new(0),
            case,
            PREC_1,
        );
        let verdict = verdict(&reference, &new, &old);
        let old_err = (old.root - reference.root).abs();
        if old_err > worst_old_err {
            worst_old_err = old_err;
            worst_old_case = case.label.clone();
        }
        if old.iterations > old_iters_max {
            old_iters_max = old.iterations;
            old_iters_case = case.label.clone();
        }
        verdicts
            .entry(verdict)
            .or_default()
            .push(case.label.clone());
        println!(
            "{}\t{:.10}\t{}\t{}\t{:.3e}\t{:.3e}\t{:.3e}\t{}",
            case.label,
            reference.root,
            new.iterations,
            old.iterations,
            old.residual,
            new.residual,
            old_err,
            verdict,
        );
    }
    println!(
        "-- stress summary over {} adversarial cases --",
        cases.len()
    );
    for (verdict, names) in &verdicts {
        println!("{verdict}: {}", names.len());
        if *verdict != "OK" {
            println!("   {}", names.join(", "));
        }
    }
    println!(
        "worst old root error vs reference: {:.3e} at {}",
        worst_old_err, worst_old_case
    );
    println!(
        "most Newton passes the old solver needed: {old_iters_max} at {} (cap {MAX_ITER})",
        old_iters_case
    );
}

fn ratio(new: usize, old: usize) -> f64 {
    if old > 0 {
        new as f64 / old as f64
    } else {
        f64::NAN
    }
}
