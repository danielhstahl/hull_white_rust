//! Unit tests for the bracketed rootfinder: the bracket is verified rather than trusted, every
//! exit is in rate units, and a solve that cannot converge says which stage failed.  These stay
//! deliberately independent of the instrument layer -- everything here is exercised on plain
//! scalar objectives (`2x - 1`, `exp(-x) - 1e-6`, `x^2 + 1`) so a change to a pricer cannot
//! hide a solver regression; the same properties are re-checked end-to-end on the production
//! Jamshidian objective in the `solver_ab` harness.

use super::*;

fn settings(tolerance: f64, max_iterations: u32) -> SolverSettings {
    SolverSettings {
        tolerance,
        max_iterations,
        initial_guess: None,
    }
}

#[test]
fn solves_a_simple_monotone_equation() {
    //f(x) = 2x - 1, root at 0.5, no bracket offered.
    let solution = solve(
        &|x| 2.0 * x - 1.0,
        &|_| 2.0,
        None,
        0.0,
        &settings(1e-12, 100),
    )
    .unwrap();
    assert!((solution.root - 0.5).abs() < 1e-12, "{solution:?}");
}

#[test]
fn a_seed_far_from_the_root_lands_on_the_same_answer() {
    //f(x) = exp(-x) - 0.3, root at ln(1/0.3) = 1.2039728043259361
    let expected = (1.0 / 0.3_f64).ln();
    let f = |x: f64| (-x).exp() - 0.3;
    let df = |x: f64| -(-x).exp();
    for seed in [-5.0, 0.0, 40.0, 500.0] {
        let solution = solve(&f, &df, None, seed, &settings(1e-12, 100)).unwrap();
        assert!(
            (solution.root - expected).abs() < 1e-10,
            "seed {seed} gave {solution:?}, expected {expected}"
        );
    }
}

#[test]
fn the_bracket_guard_stops_a_lying_derivative() {
    //A derivative of the wrong sign would send a plain Newton iteration off to infinity; the
    //bisection guard has to keep it on the root anyway.
    let solution = solve(
        &|x| x - 0.25,
        &|_| -1.0, //deliberately wrong sign
        Some((0.0, 1.0)),
        0.9,
        &settings(1e-12, 200),
    )
    .unwrap();
    assert!((solution.root - 0.25).abs() < 1e-10, "{solution:?}");
}

#[test]
fn newton_cannot_crawl_an_exponential_tail() {
    //f(x) = exp(-x) - 1e-6 has root x = ln(1e6) = 13.8155.  Deep in the tail the Newton step
    //is exactly +1 per pass regardless of how far the root is, so an iteration that is merely
    //"bracketed" still needs ~140 passes to cross a bracket this wide -- and a million for a
    //strike this extreme in the real model.  The minimum-progress guard hands the long distance
    //to bisection instead, which halves the bracket every pass.
    let root = 1e6_f64.ln();
    let f = |x: f64| (-x).exp() - 1e-6;
    let df = |x: f64| -(-x).exp();
    let solution = solve(&f, &df, None, -100.0, &settings(1e-12, 100)).unwrap();
    assert!(
        (solution.root - root).abs() < 1e-10,
        "{solution:?} vs root {root}"
    );
    assert!(solution.iterations <= 70, "{solution:?}");
}

#[test]
fn a_bracket_that_does_not_straddle_is_widened_not_trusted() {
    //The proposed bracket is entirely on the positive side of the root at -3.
    let solution = solve(
        &|x| x + 3.0,
        &|_| 1.0,
        Some((0.0, 1.0)),
        0.5,
        &settings(1e-12, 100),
    )
    .unwrap();
    assert!((solution.root + 3.0).abs() < 1e-9, "{solution:?}");
}

#[test]
fn running_out_of_iterations_says_so_with_context() {
    let error = solve(
        &|x| (-x).exp() - 1e-6,
        &|x| -(-x).exp(),
        None,
        0.0,
        &settings(1e-14, 2),
    )
    .unwrap_err();
    assert!(
        matches!(error, SolverError::Exhausted { iterations: 2, .. }),
        "{error:?}"
    );
    let text = error.to_string();
    assert!(text.contains("2 iterations"), "{text}");
    assert!(text.contains("residual"), "{text}");
}

#[test]
fn a_never_straddling_objective_reports_no_sign_change() {
    //Strictly positive: no root exists at any rate.
    let error = solve(
        &|x| x * x + 1.0,
        &|x| 2.0 * x,
        None,
        0.0,
        &settings(1e-12, 100),
    )
    .unwrap_err();
    assert!(
        matches!(error, SolverError::NoSignChange { .. }),
        "{error:?}"
    );
    assert!(error.to_string().contains("no sign change"), "{error}");
}

#[test]
fn a_nan_objective_is_reported_as_such() {
    let error = solve(
        &|x| if x > 0.0 { f64::NAN } else { x + 1.0 },
        &|_| 1.0,
        None,
        0.0,
        &settings(1e-12, 100),
    )
    .unwrap_err();
    assert!(
        matches!(error, SolverError::NonFiniteEvaluation { .. }),
        "{error:?}"
    );
}

#[test]
fn an_exact_endpoint_or_seed_root_short_circuits() {
    let settings = settings(1e-12, 100);
    let solution = solve(&|x| x - 1.0, &|_| 1.0, Some((1.0, 2.0)), 5.0, &settings).unwrap();
    assert_eq!(solution.root, 1.0);
    assert_eq!(solution.iterations, 0);
    let solution = solve(&|x| x - 7.0, &|_| 1.0, None, 7.0, &settings).unwrap();
    assert_eq!(solution.root, 7.0);
    assert_eq!(solution.iterations, 0);
}

#[test]
fn settings_are_validated_before_the_solve() {
    assert!(settings(0.0, 10).validate().is_err());
    assert!(settings(f64::NAN, 10).validate().is_err());
    assert!(settings(1e-9, 0).validate().is_err());
    assert!(
        SolverSettings {
            tolerance: 1e-9,
            max_iterations: 10,
            initial_guess: Some(f64::INFINITY),
        }
        .validate()
        .is_err()
    );
    assert!(settings(1e-9, 10).validate().is_ok());
}

#[test]
fn a_looser_tolerance_stops_looser() {
    //A nonlinear objective, so the number of iterations actually depends on the tolerance:
    //exp(x) = 2, root at ln 2.
    let f = |x: f64| x.exp() - 2.0;
    let df = |x: f64| x.exp();
    let root = 2.0_f64.ln();
    let loose = solve(&f, &df, Some((0.0, 2.0)), 0.0, &settings(1e-2, 100)).unwrap();
    let tight = solve(&f, &df, Some((0.0, 2.0)), 0.0, &settings(1e-14, 100)).unwrap();
    assert!(
        (loose.root - root).abs() <= 1e-2 && (tight.root - root).abs() <= 1e-13,
        "loose {loose:?} / tight {tight:?} vs root {root}"
    );
    assert!(loose.iterations < tight.iterations);
}
