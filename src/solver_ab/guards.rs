//! The always-on guards for the ported step rule (not env-gated).
//!
//! Three properties, pinned over the grids rather than asserted in prose: root agreement with the
//! pre-port reference, a cost ceiling relative to the old `nrfind` solve, and the far-tail crawl
//! staying shut.  84 solves cost milliseconds, so these run with every `cargo test`.

use crate::solver_ab::fixtures::{grid, stress_grid};
use crate::solver_ab::report::{Compared, compare, ratio, reference_settings};
use crate::solver_ab::variants::solve_v4;
use crate::{SolverSettings, rootfinder};

/// Guard on the step rule ported into `rootfinder::solve`, always on.
///
/// Two properties were regressions before the port and are pinned here over the production grid
/// (4 curve shapes x {4,8,48} coupons x 7 strikes = 84 cases):
///
/// * **Accuracy.**  At the shipped 1e-12 tolerance the ported loop lands on the same root as
///   `solve_v0` driven at 1e-15 -- worst disagreement `2.2e-16` measured, asserted `1e-12`.  A
///   step-acceptance rule that bought speed by trading the root would show up as this bound
///   breaking.
/// * **Cost.**  The evaluation count is capped at 2x what the old `nrfind` solve cost on the
///   same cases (measured 1,442 vs 832 = 1.73x, against the pre-port 9,092 = 10.9x) and the
///   mean iteration count at 10 (measured 5.4, against the pre-port 53.1).  Evaluation counts
///   are exact, not timed, so this assertion is machine-independent.
///
/// It is not env-gated like sections A-F: the cost of the shipped solve is the property the port
/// was for, and 84 solves cost milliseconds.
#[test]
pub(super) fn ported_step_rule_matches_a_tight_reference_and_stays_cheap() {
    let settings = SolverSettings::default();
    let cases: Vec<Compared> = grid()
        .iter()
        .map(|c| compare(c, settings, reference_settings()))
        .collect();

    let mut ported_evals = 0usize;
    let mut old_evals = 0usize;
    let mut iters_sum = 0u64;
    let mut max_iters = 0u32;
    let mut max_root_err = 0f64;
    let mut max_resid_rel = 0f64;
    let mut worst_root_case = String::new();
    let mut worst_iter_case = String::new();
    for c in &cases {
        assert!(c.ported.ok, "ported solve failed on {}", c.case.label);
        let err = (c.ported.root - c.reference.root).abs();
        if err > max_root_err {
            max_root_err = err;
            worst_root_case = c.case.label.clone();
        }
        max_resid_rel = max_resid_rel.max(c.ported.residual / c.case.strike.abs().max(1.0));
        if c.ported.iterations > max_iters {
            max_iters = c.ported.iterations;
            worst_iter_case = c.case.label.clone();
        }
        ported_evals += c.ported.f_calls + c.ported.df_calls;
        iters_sum += c.ported.iterations as u64;
        if c.old.ok {
            old_evals += c.old.f_calls + c.old.df_calls;
        }
    }
    let n = cases.len() as f64;
    println!(
        "GUARD grid n={} ported {:.1} it avg / max {} at {} / {} evals / max |f| rel {:.3e}",
        cases.len(),
        iters_sum as f64 / n,
        max_iters,
        worst_iter_case,
        ported_evals,
        max_resid_rel,
    );
    println!(
        "GUARD worst |root - v0@1e-15| = {:.3e} at {}; old evals {} (ratio {:.2})",
        max_root_err,
        worst_root_case,
        old_evals,
        ratio(ported_evals, old_evals)
    );
    assert!(
        max_root_err <= 1e-12,
        "ported solve drifted from the tight reference: {max_root_err:.3e} at {worst_root_case}"
    );
    assert!(
        max_resid_rel <= 1e-10,
        "ported solve left a {max_resid_rel:.3e} relative residual"
    );
    assert!(
        iters_sum as f64 / n <= 10.0,
        "ported solve averaged {:.1} iterations (pre-port 53.1); the step rule regressed",
        iters_sum as f64 / n
    );
    assert!(
        ported_evals <= old_evals * 2,
        "ported solve costs {ported_evals} evals vs the old solve's {old_evals} (ratio {}); \
         the port was supposed to remove an 11x premium, not keep a 2x one",
        ratio(ported_evals, old_evals)
    );
}

/// The same guard on the adversarial grid: K = 1e-6 ... 1e6 against a bond worth about 1, where
/// bare Newton fails 18 of 48 cases.  Cost is not asserted here -- the old solver burns its whole
/// iteration cap on most of these, so its evaluation count is not a baseline worth dividing by.
#[test]
pub(super) fn ported_step_rule_holds_on_the_adversarial_grid() {
    let settings = SolverSettings::default();
    let cases = stress_grid();
    let mut max_root_err = 0f64;
    let mut worst = String::new();
    let mut max_iters = 0u32;
    for case in &cases {
        let c = compare(case, settings, reference_settings());
        assert!(
            c.ported.ok,
            "ported solve failed on the adversarial case {}",
            case.label
        );
        let err = (c.ported.root - c.reference.root).abs();
        if err > max_root_err {
            max_root_err = err;
            worst = case.label.clone();
        }
        max_iters = max_iters.max(c.ported.iterations);
    }
    println!(
        "GUARD stress n={} max |root - v0@1e-15| = {:.3e} at {}, max iters {}",
        cases.len(),
        max_root_err,
        worst,
        max_iters
    );
    assert!(
        max_root_err <= 1e-12,
        "ported solve drifted from the tight reference on the stress grid: {max_root_err:.3e} at {worst}"
    );
}

/// The other half of the contract: accepting a residual-decreasing step must not reopen the
/// far-tail crawl the minimum-progress guard was written to close.  This is what fixes the
/// acceptance ratio from below -- see the rule scan in the module docs.
#[test]
pub(super) fn ported_step_rule_does_not_reopen_the_tail_crawl() {
    //Deep in `exp(-x)` the Newton step is a fixed +1 however far the root is, so an accepted
    //Newton walk across this bracket costs one unit per pass.  Bisection's own floor is ~48
    //passes to the 1e-12 rate tolerance; the ported solve exits early on its own Newton
    //correction and lands in 12.
    let tail_f = |x: f64| (-x).exp() - 1e-6;
    let tail_df = |x: f64| -(-x).exp();
    let tail_root = 1e6_f64.ln();
    let settings = SolverSettings {
        tolerance: 1e-12,
        max_iterations: 500,
        initial_guess: None,
    };
    let tail = rootfinder::solve(&tail_f, &tail_df, None, -100.0, &settings).unwrap();
    println!(
        "GUARD tail iters {} root err {:.3e}",
        tail.iterations,
        (tail.root - tail_root).abs()
    );
    assert!(tail.iterations <= 40, "tail crawl reopened: {tail:?}");
    assert!(
        (tail.root - tail_root).abs() < 1e-10,
        "tail root lost: {tail:?} vs {tail_root}"
    );
    //The strict-decrease variant (v4 as measured) takes 119 passes on the same fixture: the
    //shipped ratio is deliberately tighter than "any decrease" for exactly this reason.
    let v4 = solve_v4(&tail_f, &tail_df, None, -100.0, &settings, 0.0).unwrap();
    println!("GUARD v4 (strict decrease) tail iters {}", v4.iterations);
    assert!(
        tail.iterations < v4.iterations,
        "the shipped ratio should beat strict-decrease on the tail: {} vs {}",
        tail.iterations,
        v4.iterations
    );
}
