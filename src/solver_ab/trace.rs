//! Iteration-by-iteration tracing of a solve, reconstructed from the calls it made.
//!
//! The solver is not instrumented for this.  These loops evaluate `f(x)` then `df(x)` per iteration
//! and move the bracket by the sign of `f(x)`, and the ported loop additionally probes
//! `f(x_newton)` *after* `df(x)` for step acceptance.  Recording every call in order and replaying
//! the sign rule recovers the iterate, residual, Newton correction, probe and bracket per iteration
//! -- which is what shows where the passes go and which predicate accepted each step.

use crate::bonds::coupon_bond_generic_t;
use crate::solver_ab::fixtures::{Case, Counted, grid};
use crate::solver_ab::measurement::{Model, SolveFn};
use crate::solver_ab::variants::solve_v2;
use crate::{HullWhite, SolverSettings, rootfinder};
use std::sync::Arc;

/// One recorded call on the objective or its derivative, in the order the solver asked for it.
#[derive(Debug, Clone, Copy)]
pub(super) enum Event {
    F { x: f64, value: f64 },
    D { x: f64, value: f64 },
}

/// One loop iteration, recovered from the call log.
#[derive(Debug, Clone, Copy)]
pub(super) struct Step {
    x: f64,
    f: f64,
    slope: f64,
    /// The step-acceptance probe `f(x_newton)`, when the iteration evaluated one.
    probe: Option<(f64, f64)>,
}

/// Step-acceptance thresholds mirrored from `rootfinder`, for labelling the trace only:
/// `MIN_NEWTON_FRACTION` (a step must cross a tenth of the bracket to count as long-range work)
/// and `MAX_MODEL_ERROR_RATIO` (the tangent line's error at the probed point, relative to the
/// residual being corrected).  Kept here as literals because `rootfinder` keeps them private; if
/// they ever disagree with the shipped values the trace labels are wrong, so `trace_matches_the_
/// shipped_iteration_count` is the check on this pair.
pub(super) const TRACE_MIN_NEWTON_FRACTION: f64 = 0.1;

pub(super) const TRACE_MAX_MODEL_ERROR_RATIO: f64 = 0.25;

/// Reconstructs a solve's trajectory from the calls it made, without touching the solver.
///
/// These loops evaluate `f(x)` and then `df(x)` per iteration and move the bracket by the sign of
/// `f(x)`.  The ported loop additionally evaluates `f(x_newton)` *after* `df(x)`, as the
/// step-acceptance probe.  Recording the calls in the order they happened and replaying the sign
/// rule gives the iterate, residual, Newton correction, probe and bracket per iteration -- which
/// is what shows where the passes go and which predicate accepted each step.
///
/// The order matters and is the whole trick: an `f` call is a loop iterate exactly when the call
/// after it is a `df` at the same point.  A stray `f` is the probe.
pub(super) struct Trace {
    events: Vec<Event>,
}

pub(super) fn traced_solve(
    hw: &Model<'_>,
    case: &Case,
    settings: &SolverSettings,
    solver: SolveFn<'_>,
) -> (Trace, Result<rootfinder::Solution, rootfinder::SolverError>) {
    let events: Arc<std::sync::Mutex<Vec<Event>>> = Arc::new(std::sync::Mutex::new(Vec::new()));
    let objective = {
        let events = events.clone();
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
            events.lock().unwrap().push(Event::F { x: rate, value });
            value
        }
    };
    let derivative = {
        let events = events.clone();
        move |rate: f64| {
            let value = hw.coupon_bond_price_t_deriv(
                rate,
                case.option_maturity,
                &case.schedule,
                case.coupon_rate,
            );
            events.lock().unwrap().push(Event::D { x: rate, value });
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
        events: events.lock().unwrap().clone(),
    };
    (trace, outcome)
}

impl Trace {
    /// Iterations recovered from the call order, probes attached to the step that asked for them.
    ///
    /// `f` is called three times before the loop runs on the bracketed path -- the seed, then each
    /// end of the proposed bracket -- so those are skipped.  A trailing `f` with no `df` after it
    /// is the run returning on `f(x) == 0`, unless its point is the previous step's Newton point,
    /// in which case it was that step's probe and the run then exhausted.
    pub(super) fn steps(&self) -> Vec<Step> {
        let mut pre_loop = 0usize;
        let mut steps: Vec<Step> = Vec::new();
        for (i, event) in self.events.iter().enumerate() {
            let Event::F { x, value } = *event else {
                continue;
            };
            if pre_loop < 3 {
                pre_loop += 1;
                continue;
            }
            match self.events.get(i + 1) {
                Some(Event::D { x: dx, value: dv }) if *dx == x => {
                    //Is the call after the derivative the step probe, or the next iterate?  Only
                    //the probe is *not* followed by its own derivative at the same point.
                    let probe = match self.events.get(i + 2) {
                        Some(Event::F { x: px, value: pv }) => {
                            let next_is_head = matches!(
                                self.events.get(i + 3),
                                Some(Event::D { x: dx2, .. }) if *dx2 == *px
                            );
                            if next_is_head { None } else { Some((*px, *pv)) }
                        }
                        _ => None,
                    };
                    steps.push(Step {
                        x,
                        f: value,
                        slope: *dv,
                        probe,
                    });
                }
                _ => {
                    let is_previous_probe = steps.last().is_some_and(|s| {
                        s.slope.is_finite() && s.slope != 0.0 && s.x - s.f / s.slope == x
                    });
                    if !is_previous_probe {
                        steps.push(Step {
                            x,
                            f: value,
                            slope: f64::NAN,
                            probe: None,
                        });
                    }
                }
            }
        }
        steps
    }

    /// Prints one line per iteration: iterate, residual, Newton correction, the step probe and
    /// what the acceptance predicates made of it, and the bracket the replayed sign rule says the
    /// solver was holding.
    pub(super) fn print(&self, label: &str, lo: f64, hi: f64) {
        let f_pre = |index: usize| match self
            .events
            .iter()
            .filter(|e| matches!(e, Event::F { .. }))
            .nth(index)
        {
            Some(Event::F { value, .. }) => *value,
            _ => f64::NAN,
        };
        let mut a = lo;
        let mut b = hi;
        let mut f_a = f_pre(1);
        println!(
            "{label}\tITER\tX\tF(X)\tDF(X)\tCORRECTION\tPROBE_X\tPROBE_F\tTANGENT_ERR\tDECISION\t\
             BRACKET\tWIDTH"
        );
        for (iteration, step) in self.steps().iter().enumerate() {
            //Replay the bracket update exactly as the loop does it, before measuring the width the
            //acceptance test saw.
            if step.f.is_finite() && (step.f < 0.0) == (f_a < 0.0) {
                a = step.x;
                f_a = step.f;
            } else if step.f.is_finite() {
                b = step.x;
            }
            let (low, high) = (a.min(b), a.max(b));
            let width = high - low;
            let correction = if step.slope.is_finite() && step.slope != 0.0 {
                step.f / step.slope
            } else {
                f64::NAN
            };
            let newton = step.x - correction;
            let inside = correction.is_finite() && newton > low && newton < high;
            let makes_progress = inside && correction.abs() >= TRACE_MIN_NEWTON_FRACTION * width;
            let improves = inside
                && !makes_progress
                && step
                    .probe
                    .is_some_and(|(_, pv)| pv.abs() <= TRACE_MAX_MODEL_ERROR_RATIO * step.f.abs());
            let decision = if step.probe.is_some() && !inside {
                "bisect(probe outside)"
            } else if makes_progress {
                "newton(progress)"
            } else if improves {
                "newton(improves)"
            } else {
                "bisect"
            };
            let (probe_x, probe_f, tangent_err) = match step.probe {
                Some((px, pv)) => (
                    px,
                    pv,
                    if step.f.abs() > 0.0 {
                        pv.abs() / step.f.abs()
                    } else {
                        f64::NAN
                    },
                ),
                None => (f64::NAN, f64::NAN, f64::NAN),
            };
            println!(
                "{label}\t{}\t{x:.12}\t{f:.6e}\t{slope:.6e}\t{correction:.6e}\t{probe_x:.12}\t\
                 {probe_f:.6e}\t{tangent_err:.4}\t{decision}\t[{low:.6}, {high:.6}]\t{width:.6e}",
                iteration + 1,
                x = step.x,
                f = step.f,
                slope = step.slope,
            );
        }
    }
}

#[test]
pub(super) fn trace_one_solve() {
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
