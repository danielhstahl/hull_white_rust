//! Section E: the exit-rule variants, plus the safety properties every variant has to keep.
//!
//! A faster solve that reintroduces the far-tail crawl or trusts a lying derivative is not a fix,
//! so each candidate is re-checked against the deep exponential tail, a wrong-sign derivative and a
//! rootless objective.

use crate::solver_ab::fixtures::{Counted, grid};
use crate::solver_ab::measurement::{Run, run_new, run_variant};
use crate::solver_ab::report::{lying, no_root, ratio};
use crate::solver_ab::variants::{solve_v2, solve_v3, solve_v4};
use crate::{HullWhite, SolverSettings, rootfinder};
use std::sync::atomic::AtomicUsize;

/// The safety properties the shipped solver was written to have, re-checked against the candidate.
/// These mirror the unit tests in `rootfinder`: a fast solve that reintroduces the tail crawl or
/// trusts a lying derivative is not a fix.
#[test]
pub(super) fn candidate_keeps_the_safety_properties() {
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
pub(super) fn exit_rule_variants() {
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
pub(super) fn variant_safety_checks() {
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
