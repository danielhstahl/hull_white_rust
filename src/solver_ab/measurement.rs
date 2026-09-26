//! The measurement plumbing: the counted objective and derivative the harness solves, and the
//! `Run` record every solver variant is normalised into.
//!
//! Every zero-coupon leg price funnels through the yield and forward curves, so a count of curve
//! calls (`Counted`) is a measure of model work no solver can game; the `f` / `df` counts are the
//! solver's own asks.  `run_new`, `run_old` and `run_variant` build the *same* closures and differ
//! only in the solver they hand them to -- the only independent variable is the solve.

use crate::bonds::coupon_bond_generic_t;
use crate::solver_ab::fixtures::{Case, MAX_ITER, R_INIT};
use crate::{HullWhite, SolverSettings, rootfinder};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

/// The shape shared by every solve routine measured here: objective, derivative, an optional
/// bracket, a seed, and the settings.  Four call sites take it, so it gets a name.
pub(super) type SolveFn<'a> = &'a dyn Fn(
    &dyn Fn(f64) -> f64,
    &dyn Fn(f64) -> f64,
    Option<(f64, f64)>,
    f64,
    &SolverSettings,
) -> Result<rootfinder::Solution, rootfinder::SolverError>;

pub(super) type Curve = Box<dyn Fn(f64) -> f64 + Sync>;

pub(super) type Model<'h> = HullWhite<'h, Curve, Curve>;

/// One solve's worth of measurements.
pub(super) struct Run {
    pub(super) ok: bool,
    pub(super) root: f64,
    pub(super) residual: f64,
    pub(super) iterations: u32,
    pub(super) f_calls: usize,
    pub(super) df_calls: usize,
    pub(super) curve_calls: usize,
    pub(super) setup_curve_calls: usize,
    pub(super) bracket: Option<(f64, f64)>,
    pub(super) ns: f64,
}

impl Run {
    pub(super) fn failed() -> Self {
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
    pub(super) fn width(&self) -> f64 {
        self.bracket
            .map(|(lo, hi)| (hi - lo).abs())
            .unwrap_or(f64::NAN)
    }
}

/// The production objective: the coupon bond's value at option expiry, less the strike.
pub(super) fn objective<'h>(
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
pub(super) fn derivative<'h>(
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
pub(super) fn run_new<'h>(
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
pub(super) fn run_old<'h>(
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

/// Runs a solve variant with the same counted objective/derivative plumbing.
pub(super) fn run_variant(
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

/// A straddle test that treats NaN as "does not straddle", matching the shipped solver.
pub(super) fn straddles(f_lower: f64, f_upper: f64) -> bool {
    (f_lower <= 0.0 && f_upper >= 0.0) || (f_lower >= 0.0 && f_upper <= 0.0)
}
