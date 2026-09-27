//! Tests for the curve abstraction and the consistency relation the constructors enforce.
//!
//! The pair of closures this crate used to take could be mutually inconsistent and nothing said
//! so.  These are the tests that make that impossible: an accepted consistent pair, a rejected
//! inconsistent one, the boundary between them, and the two degenerate cases (a curve that cannot
//! be probed at all, and a curve whose forward is not a number).

use approx::*;

use crate::HullWhite;
use crate::curves::{
    CURVE_PROBE_TIMES, FORWARD_CONSISTENCY_RELATIVE_TOLERANCE, FORWARD_CONSISTENCY_TOLERANCE,
    ForwardInconsistency, YieldCurve, derivative, forward_consistency_tolerance, from_yield,
    from_yield_and_forward, max_forward_inconsistency, validate_curve,
};
use crate::test_support::ALL_SCENARIOS;

/// `y(t) = 0.05 t + 0.01 t^2`, so `f(0,t) = 0.05 + 0.02 t`.  The running example of a pair
/// that is actually one curve.
fn consistent_curve() -> impl YieldCurve {
    from_yield_and_forward(|t: f64| 0.05 * t + 0.01 * t * t, |t: f64| 0.05 + 0.02 * t)
}

/// The same yield with a forward that is *almost* right: `0.05 + 0.01 t` instead of `+ 0.02 t`.
/// The kind of slip the old API accepted silently at every price in the book.
fn half_slope_forward_curve() -> impl YieldCurve {
    from_yield_and_forward(|t: f64| 0.05 * t + 0.01 * t * t, |t: f64| 0.05 + 0.01 * t)
}

// --- accepted ---------------------------------------------------------------

/// A curve whose forward is the derivative of its own yield builds, on every probe.
#[test]
fn a_consistent_pair_is_accepted() {
    let curve = consistent_curve();
    validate_curve(&curve).unwrap();
    let model = HullWhite::new(0.15, 0.02, &curve);
    assert!(model.is_ok(), "consistent pair rejected: {model:?}");
}

/// Every fixture scenario is internally consistent, so every scenario builds a model.  This is
/// the same property the `test_support` tests pin algebraically, checked through the *constructor*
/// the tests actually use — if the check and the fixture ever drifted apart, this is where it shows.
#[test]
fn every_scenario_curve_is_accepted() {
    for s in ALL_SCENARIOS {
        let curve = s.curve();
        let bad = max_forward_inconsistency(&curve);
        assert!(
            bad.is_some_and(|b| b.error <= b.tolerance),
            "{} rejected by the consistency check: {bad:?}",
            s.name
        );
        assert!(
            HullWhite::new(s.a, s.sigma, &curve).is_ok(),
            "{} could not calibrate a model",
            s.name
        );
    }
}

/// A curve that implements only `zero_yield` has no pair to be inconsistent: its forward *is*
/// the numerical derivative the check computes, so the check is satisfied by construction.
#[test]
fn a_yield_only_curve_is_accepted_because_its_forward_is_derived() {
    // A curve that is not a polynomial, so the derived forward is not exactly the analytic one.
    let curve = from_yield(|t: f64| 0.05 * t + 0.005 * t * t + 0.001 * (2.0 * t).sin());
    validate_curve(&curve).unwrap();
    let model = HullWhite::new(0.2, 0.03, &curve);
    assert!(model.is_ok(), "yield-only curve rejected: {model:?}");

    // The derived forward tracks the analytic derivative (0.05 + 0.01 t + 0.002 cos 2t) to
    // well inside the tolerance -- this is the accuracy the trait doc promises for the default.
    for t in CURVE_PROBE_TIMES {
        let analytic = 0.05 + 0.01 * t + 0.002 * (2.0 * t).cos();
        assert_abs_diff_eq!(curve.forward(t), analytic, epsilon = 1e-8);
    }
}

/// The tolerance is a floor, not a cliff: an error *at* the tolerance is still a pass.  With a
/// `~0.05` forward the effective tolerance is `1e-4 * 0.05 = 5e-6`, so a sub-microstaff offset
/// must not fail.
#[test]
fn an_error_inside_the_tolerance_is_accepted() {
    let epsilon = FORWARD_CONSISTENCY_TOLERANCE * 0.5;
    let curve = from_yield_and_forward(
        |t: f64| 0.05 * t + 0.01 * t * t,
        |t: f64| 0.05 + 0.02 * t + epsilon,
    );
    let bad = max_forward_inconsistency(&curve).unwrap();
    assert!(
        bad.error <= bad.tolerance,
        "a {epsilon:e} offset should fit under a {} tolerance, got error {:e} vs tolerance {:e}",
        FORWARD_CONSISTENCY_TOLERANCE,
        bad.error,
        bad.tolerance
    );
    assert!(HullWhite::new(0.15, 0.02, &curve).is_ok());
}

// --- rejected ---------------------------------------------------------------

/// A pair that disagrees by more than the tolerance is refused, with a message that says what
/// disagreed, where, and by how much against what.
#[test]
fn an_inconsistent_pair_is_rejected() {
    let curve = half_slope_forward_curve();
    let err = HullWhite::new(0.15, 0.02, &curve).unwrap_err();
    let message = err.to_string();
    assert!(
        message.contains("inconsistent"),
        "expected an inconsistency error, got {message}"
    );
    //The message has to be actionable: which time, which values, which tolerance.
    assert!(
        message.contains("d/dt zero_yield"),
        "message lacks the derivative: {message}"
    );
    assert!(
        message.contains("forward("),
        "message lacks the supplied forward: {message}"
    );
    assert!(
        message.contains("tolerance"),
        "message lacks the tolerance: {message}"
    );
    //And the numbers in it are the real ones, not a generic complaint.
    let bad = max_forward_inconsistency(&curve).unwrap();
    assert!(
        message.contains(&format!("{}", bad.t)),
        "message should name the worst probe time {}; got {message}",
        bad.t
    );
    assert!(
        bad.error > bad.tolerance,
        "the pair the constructor rejected was actually inside tolerance: {bad:?}"
    );
}

/// `f = y'` for a *cumulative* yield.  A caller who hands over the annualised zero instead of
/// the integral fails the same check -- here `y(t) = 0.05` (constant, i.e. an annualised zero)
/// against the forward its cumulative self would have, which is not `0`.
#[test]
fn an_annualised_zero_supplied_as_the_cumulative_yield_is_rejected() {
    // y(t) = 0.05 for all t means P(0,t) = e^-0.05: a flat, one-instant curve whose derivative
    // is 0, not the 0.05 the caller's "forward" closure returns.
    let curve = from_yield_and_forward(|_t: f64| 0.05, |_t: f64| 0.05);
    let err = validate_curve(&curve).unwrap_err();
    assert!(
        err.to_string().contains("inconsistent"),
        "expected the annualised-convention slip to be caught, got {err}"
    );
}

/// The pair that every pre-0.9 example in this crate used: a linear cumulative yield and
/// `ln(t)` as the "forward".  They have nothing to do with each other, and it used to price.
#[test]
fn the_legacy_ln_forward_pair_is_rejected() {
    let curve = from_yield_and_forward(|t: f64| 0.05 * t, |t: f64| t.ln());
    let err = HullWhite::new(0.2, 0.03, &curve).unwrap_err();
    assert!(
        err.to_string().contains("inconsistent"),
        "the legacy ln-forward pair should not calibrate: {err}"
    );
}

/// The reported error is the *worst* probe, so a mismatch that grows with time (a forward that is
/// right at the front and wrong further out -- the usual way the two get edited apart) is caught
/// at the far end rather than averaged away.
#[test]
fn the_worst_probe_is_reported() {
    // error(t) = 0.001 * t: monotonic in t, so the worst probe is the largest one.
    let curve = from_yield_and_forward(
        |t: f64| 0.05 * t + 0.01 * t * t,
        |t: f64| 0.05 + 0.02 * t + 0.001 * t,
    );
    let bad: ForwardInconsistency = max_forward_inconsistency(&curve).unwrap();
    assert_abs_diff_eq!(
        bad.t,
        CURVE_PROBE_TIMES[CURVE_PROBE_TIMES.len() - 1],
        epsilon = 0.0
    );
    assert_abs_diff_eq!(bad.error, 0.001 * bad.t, epsilon = 1e-9);
}

/// A curve whose forward is not a number while its yield is fine is inconsistent by definition:
/// it claims a slope that does not exist.
#[test]
fn a_non_finite_forward_over_a_finite_yield_is_rejected() {
    struct NanForward;
    impl YieldCurve for NanForward {
        fn zero_yield(&self, t: f64) -> f64 {
            0.05 * t
        }
        fn forward(&self, _t: f64) -> f64 {
            f64::NAN
        }
    }
    let bad = max_forward_inconsistency(&NanForward).unwrap();
    assert!(
        bad.error.is_infinite(),
        "a NAN forward should read as an infinite error: {bad:?}"
    );
    assert!(validate_curve(&NanForward).is_err());
}

/// The same failure in its *infinite* form, which is the dangerous one.
///
/// A `NaN` propagates and shows up as a `NaN` price three frames later; an infinite forward is
/// worse, because a naive tolerance of the form `max(abs_floor, rel * |supplied|)` becomes `inf`
/// at an infinite `supplied`, and `inf > inf` is false -- so the check reports "consistent" for
/// a curve that claims every forward rate is infinite, and then `phi(t)` carries that `inf` into
/// `mu_r` / the Jamshidian bracket and produces a finite, wrong-looking number.
#[test]
fn an_infinite_forward_over_a_finite_yield_is_rejected() {
    for infinity in [f64::INFINITY, f64::NEG_INFINITY] {
        let curve = from_yield_and_forward(|t: f64| 0.05 * t, move |_t: f64| infinity);
        let bad = max_forward_inconsistency(&curve).unwrap_or_else(|| {
            panic!("{infinity} forward should still be probeable via its yield")
        });
        assert!(
            bad.error.is_infinite() && bad.error > bad.tolerance,
            "a {infinity} forward must be outside its tolerance, got {bad:?}"
        );
        assert!(
            bad.tolerance.is_finite(),
            "the tolerance may not become infinite just because the supplied forward is: {bad:?}"
        );
        let err = validate_curve(&curve).unwrap_err();
        assert!(
            err.to_string().contains("inconsistent"),
            "a {infinity} forward should be reported as inconsistent, got {err}"
        );
        assert!(
            HullWhite::new(0.2, 0.03, &curve).is_err(),
            "a model must not build on a {infinity} forward"
        );
    }
}

/// The `Sync`-trait-object route, for the same claim: an implementor that returns an infinite
/// forward is rejected exactly like a closure pair that does.
#[test]
fn an_infinite_forward_from_a_trait_implementor_is_rejected() {
    struct InfForward;
    impl YieldCurve for InfForward {
        fn zero_yield(&self, t: f64) -> f64 {
            0.05 * t
        }
        fn forward(&self, _t: f64) -> f64 {
            f64::INFINITY
        }
    }
    assert!(validate_curve(&InfForward).is_err());
    assert!(HullWhite::new(0.2, 0.03, &InfForward).is_err());
}

/// A curve that cannot be evaluated anywhere cannot be calibrated either.  Silence here would mean
/// "checked and fine" about a curve that was never looked at.
#[test]
fn an_unprobeable_curve_is_rejected_rather_than_trusted() {
    struct Undefined;
    impl YieldCurve for Undefined {
        fn zero_yield(&self, _t: f64) -> f64 {
            f64::NAN
        }
    }
    assert!(max_forward_inconsistency(&Undefined).is_none());
    let err = validate_curve(&Undefined).unwrap_err();
    assert!(
        err.to_string().contains("probe"),
        "expected an un-probeable curve error, got {err}"
    );
    assert!(HullWhite::new(0.2, 0.03, &Undefined).is_err());
}

/// The model's own parameters are still checked, and checked before the curve: a bad `a` is a
/// bad model whatever the curve says.
#[test]
fn bad_model_parameters_fail_before_the_curve_is_considered() {
    let curve = half_slope_forward_curve();
    assert!(HullWhite::new(0.0, 0.02, &curve).is_err());
    assert!(HullWhite::new(0.15, -1.0, &curve).is_err());
    assert!(HullWhite::new(f64::NAN, 0.02, &curve).is_err());
    //...and the curve's own complaint is a different one, so the two are not conflated.
    let curve_err = HullWhite::new(0.15, 0.02, &curve).unwrap_err().to_string();
    assert!(curve_err.contains("inconsistent"), "{curve_err}");
}

// --- the deprecated two-closure constructor ---------------------------------

/// The deprecated `init` keeps working for a consistent pair and gains the diagnostic the old one
/// never had for an inconsistent one.
#[test]
#[allow(deprecated)]
fn the_deprecated_constructor_runs_the_same_check() {
    let ok_yield = |t: f64| 0.05 * t + 0.01 * t * t;
    let ok_forward = |t: f64| 0.05 + 0.02 * t;
    let model = HullWhite::init(0.15, 0.02, &ok_yield, &ok_forward);
    assert!(model.is_ok(), "consistent closures rejected: {model:?}");
    //The model it builds behaves like the one built the new way.
    let new_curve = consistent_curve();
    let new_way = HullWhite::new(0.15, 0.02, &new_curve).unwrap();
    assert_abs_diff_eq!(
        model.unwrap().bond_price_now(3.0).unwrap(),
        new_way.bond_price_now(3.0).unwrap(),
        epsilon = 1e-15
    );

    let bad_forward = |t: f64| 0.05 + 0.01 * t;
    let err = HullWhite::init(0.15, 0.02, &ok_yield, &bad_forward).unwrap_err();
    assert!(
        err.to_string().contains("inconsistent"),
        "the deprecated path should reject an inconsistent pair too, got {err}"
    );
}

// --- the pieces -------------------------------------------------------------

/// `discount(t) = exp(-y(t))`, the identity every `_now` price rests on.
#[test]
fn discount_is_the_exponential_of_minus_the_cumulative_yield() {
    let curve = consistent_curve();
    for t in [0.0, 0.25, 1.0, 5.0, 20.0] {
        let y = curve.zero_yield(t);
        assert_abs_diff_eq!(curve.discount(t), (-y).exp(), epsilon = 1e-15);
        assert_abs_diff_eq!(
            HullWhite::new(0.15, 0.02, &curve)
                .unwrap()
                .bond_price_now(t)
                .unwrap(),
            curve.discount(t),
            epsilon = 1e-15
        );
    }
}

/// A curve can override `discount` directly (a bootstrapped curve stores prices, not yields); the
/// trait still describes one thing, so it must also be able to report the yield.
#[test]
fn a_curve_stored_as_discount_factors_agrees_with_its_yield() {
    struct DiscountStored;
    impl YieldCurve for DiscountStored {
        fn zero_yield(&self, t: f64) -> f64 {
            0.05 * t + 0.01 * t * t
        }
        fn discount(&self, t: f64) -> f64 {
            (-self.zero_yield(t)).exp()
        }
        fn forward(&self, t: f64) -> f64 {
            0.05 + 0.02 * t
        }
    }
    let stored = DiscountStored;
    let direct = consistent_curve();
    validate_curve(&stored).unwrap();
    for t in [0.5, 1.0, 4.0, 10.0] {
        assert_abs_diff_eq!(stored.discount(t), direct.discount(t), epsilon = 1e-15);
    }
}

/// The finite-difference derivative: exact for quartics (the truncation term is `f^(5)`, which a
/// quartic does not have), and within rounding for a real curve's smooth functions.
#[test]
fn the_finite_difference_derivative_is_exact_below_a_quintic() {
    // d/dt (a + b t + c t^2 + d t^3 + e t^4) at t = 1.7
    let f = |t: f64| 1.0 + 2.0 * t + 3.0 * t * t + 4.0 * t.powi(3) + 5.0 * t.powi(4);
    let exact = |t: f64| 2.0 + 6.0 * t + 12.0 * t * t + 20.0 * t.powi(3);
    for t in [0.0, 0.01, 0.25, 1.0, 1.7, 10.0] {
        assert_abs_diff_eq!(
            derivative(&f, t),
            exact(t),
            epsilon = (1e-6f64).max(1e-7 * exact(t).abs())
        );
    }
    //And on a function this crate actually differentiates: an exponential.
    let g = |t: f64| (-0.1 * t).exp();
    for t in [0.0, 1.0, 5.0] {
        assert_abs_diff_eq!(derivative(&g, t), -0.1 * (-0.1 * t).exp(), epsilon = 1e-9);
    }
}

/// The derivative never evaluates `f` at a negative time.  Every time in this crate is measured
/// from "now" (0), and a curve asked for `t < 0` is a curve asked to extrapolate into a domain
/// it does not have.
#[test]
fn the_finite_difference_derivative_stays_on_the_domain() {
    use std::sync::atomic::{AtomicBool, Ordering};
    let saw_negative = AtomicBool::new(false);
    let f = |t: f64| {
        if t < 0.0 {
            saw_negative.store(true, Ordering::Relaxed);
        }
        0.05 * t
    };
    let _ = derivative(&f, 0.0);
    let _ = derivative(&f, 1e-6);
    assert!(
        !saw_negative.load(Ordering::Relaxed),
        "the derivative stencil left the domain (t < 0)"
    );
}

/// The tolerance rule, stated where it is used: absolute floor, relative to the forward above it.
#[test]
fn the_tolerance_is_an_absolute_floor_that_grows_with_the_forward() {
    assert_abs_diff_eq!(
        forward_consistency_tolerance(0.0),
        FORWARD_CONSISTENCY_TOLERANCE,
        epsilon = 0.0
    );
    // A 5% forward: the relative term (1e-4 * 0.05 = 5e-6) is above the 1e-6 floor.
    assert_abs_diff_eq!(
        forward_consistency_tolerance(0.05),
        FORWARD_CONSISTENCY_RELATIVE_TOLERANCE * 0.05,
        epsilon = 0.0
    );
    assert_eq!(
        forward_consistency_tolerance(-0.05),
        forward_consistency_tolerance(0.05)
    );
}
