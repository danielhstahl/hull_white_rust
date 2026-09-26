//! Section F: the adversarial grid.  The criterion that matters for a rewrite bought for
//! robustness -- is there any case where bare Newton fails, or answers with the wrong root, while
//! the bracketed solve answers correctly?

use crate::solver_ab::fixtures::{Counted, MAX_ITER, PREC_1, stress_grid};
use crate::solver_ab::measurement::{run_new, run_old};
use crate::solver_ab::report::verdict;
use crate::{HullWhite, SolverSettings};
use std::sync::atomic::{AtomicUsize, Ordering};

/// The criterion that matters for a rewrite bought for robustness: is there any case where bare
/// Newton fails, or answers with the wrong root, while the bracketed solve answers correctly?
#[test]
pub(super) fn old_solver_stress() {
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
