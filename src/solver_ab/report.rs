//! Reducing `Run`s into the numbers and verdicts the harness prints: aggregate counts, medians,
//! the old-vs-new classification (`OK` / `OLD_FAIL` / `OLD_WRONG_ROOT` / ...), and the single-case
//! driver (`compare`) the always-on guards use.

use crate::solver_ab::fixtures::{Case, Counted, PREC_1};
use crate::solver_ab::measurement::{Run, SolveFn, run_new, run_old, run_variant};
use crate::solver_ab::reference_loop::solve_v0;
use crate::{HullWhite, SolverSettings};
use std::sync::atomic::AtomicUsize;

pub(super) fn median(values: &mut [f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    values[values.len() / 2]
}

pub(super) fn times_len(reps: usize) -> usize {
    reps.max(1)
}

pub(super) fn ratio(new: usize, old: usize) -> f64 {
    if old > 0 {
        new as f64 / old as f64
    } else {
        f64::NAN
    }
}

/// Classifies one case against the 1e-15 reference solve.
///
/// "Wrong root" is the failure that matters most here: the old solver converges on the size of
/// its own step, so on a flat stretch of the objective it can stop far from the root and still
/// report success.  Thresholds are in rate units.
pub(super) fn verdict(reference: &Run, new: &Run, old: &Run) -> &'static str {
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

#[derive(Default)]
pub(super) struct Totals {
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

impl Totals {
    pub(super) fn record(
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

    pub(super) fn report(&self) {
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

/// A lying derivative (wrong sign) must still land on the root via the bracket guard.
pub(super) fn lying(solver: SolveFn<'_>) -> bool {
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
pub(super) fn no_root(solver: SolveFn<'_>) -> bool {
    let settings = SolverSettings {
        tolerance: 1e-12,
        max_iterations: 500,
        initial_guess: None,
    };
    solver(&|x| x * x + 1.0, &|x| 2.0 * x, None, 0.0, &settings)
        .map(|_| false)
        .unwrap_or(true)
}

/// Runs one grid case through the ported solver, its tight reference and the old `nrfind` solve.
pub(super) struct Compared {
    pub(super) case: Case,
    pub(super) ported: Run,
    pub(super) reference: Run,
    pub(super) old: Run,
}

pub(super) fn compare(
    case: &Case,
    settings: SolverSettings,
    reference_settings: SolverSettings,
) -> Compared {
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
        &solve_v0,
    );
    let ported = run_new(
        &hw,
        &counted.calls,
        &AtomicUsize::new(0),
        &AtomicUsize::new(0),
        case,
        settings,
    );
    let old = run_old(
        &hw,
        &counted.calls,
        &AtomicUsize::new(0),
        &AtomicUsize::new(0),
        case,
        PREC_1,
    );
    Compared {
        case: case.clone(),
        ported,
        reference,
        old,
    }
}

/// The reference every guard below asks about: the pre-port loop driven far past any usable
/// accuracy.  1e-15 on a rate is above f64's own resolution so the loop runs to its bracket
/// floor; the budget is sized for that.
pub(super) fn reference_settings() -> SolverSettings {
    SolverSettings {
        tolerance: 1e-15,
        max_iterations: 500,
        initial_guess: None,
    }
}
