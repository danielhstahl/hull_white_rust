//! Test fixtures for the solver A/B harness.
//!
//! A `Case` is one instrument -- a curve shape, a coupon schedule, a strike, a valuation point.
//! `grid()` is the 84-case production grid every table in the crate docs is computed from
//! (4 curve shapes x {4, 8, 48} coupons x 7 strikes); `stress_grid()` is the 48-case adversarial
//! grid, built to break bare Newton rather than to be realistic.
//!
//! The constants are the configuration of the solver this crate replaced (`nrfind`), kept verbatim
//! so "what would the old solve have cost, and what would it have answered" stays an apples-to-
//! apples question.

use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

// The pre-rootfinder configuration, verbatim from the implementation this replaced: bare Newton
// from a hard-coded 3% (`nrfind::find_root(&f, &df, R_INIT, PREC_1, MAX_ITER)`), converged on
// the size of the step it took rather than on the root.
pub(super) const R_INIT: f64 = 0.03;
pub(super) const PREC_1: f64 = 0.0000001;
pub(super) const MAX_ITER: i32 = 50;

/// Old-solver knobs pushed to roughly the new solver's default accuracy, to ask what the old
/// iteration would have cost if it had been asked for the same answer.
pub(super) const PREC_TIGHT: f64 = 1e-12;

/// New solver loosened to roughly where the old one stopped (old converged on a *step* of 1e-7,
/// about a bracket of that half-width).
pub(super) const TOL_MATCHED: f64 = 1e-7;

pub(super) fn reps() -> usize {
    std::env::var("AB_REPS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(25)
}

/// Counting stand-ins for the market curves.
///
/// Every zero-coupon leg price funnels through the yield and forward curves, so a count of curve
/// calls is a measure of model work that no solver can game.
pub(super) struct Counted {
    pub(super) yield_curve: Box<dyn Fn(f64) -> f64 + Sync>,
    pub(super) forward_curve: Box<dyn Fn(f64) -> f64 + Sync>,
    pub(super) calls: Arc<AtomicUsize>,
}

impl Counted {
    pub(super) fn new(curr: f64, a: f64, b: f64, sig: f64) -> Self {
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

#[derive(Clone)]
pub(super) struct Case {
    pub(super) label: String,
    pub(super) curr: f64,
    pub(super) a: f64,
    pub(super) b: f64,
    pub(super) sig: f64,
    pub(super) schedule: Vec<f64>,
    pub(super) coupon_rate: f64,
    pub(super) strike: f64,
    pub(super) r_t: f64,
    pub(super) t: f64,
    pub(super) option_maturity: f64,
}

pub(super) fn grid() -> Vec<Case> {
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

/// Cases built to break the old solve rather than to be realistic.
///
/// Bare Newton has three exposed flanks: a hard-coded seed (3%) with no bracket to fall back on, a
/// 50-iteration cap, and convergence judged on the step it took.  Puts the root far from 3%, put
/// the price scale far from 1, and see whether it still answers.
pub(super) fn stress_grid() -> Vec<Case> {
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
