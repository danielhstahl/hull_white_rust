//! Sections A-D of the A/B harness: per-case cost, aggregate counts, wall clock per solve, and the
//! candidate fix measured against both.  Env-gated on `SOLVER_AB=1`; the numbers it prints are the
//! ones quoted in the harness module docs.

use crate::solver_ab::fixtures::{Case, Counted, PREC_1, PREC_TIGHT, TOL_MATCHED, grid, reps};
use crate::solver_ab::measurement::{Run, derivative, objective, run_new, run_old, run_variant};
use crate::solver_ab::report::{Totals, median, ratio, times_len, verdict};
use crate::solver_ab::variants::solve_v2;
use crate::{HullWhite, SolverSettings};
use std::sync::atomic::{AtomicUsize, Ordering};

#[test]
pub(super) fn solver_ab_cost() {
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
