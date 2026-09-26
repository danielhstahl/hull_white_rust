//! Unit tests for the model primitives.
//!
//! `steep_fixture_actually_makes_phi_time_dependent` guards the fixture the time-coordinate
//! regressions elsewhere rely on: if `phi` were constant, a wrong clock would cancel out and
//! those tests would pass for the wrong reason.

use approx::*;

use crate::HullWhite;
use crate::error::HullWhiteError;
use crate::test_support::STEEP_CURVE;

#[test]
fn steep_fixture_actually_makes_phi_time_dependent() {
    //Protects the tests above: if phi were constant, a time-coordinate error would cancel out
    //and the regression tests would pass for the wrong reason.
    let (yield_curve, forward_curve) = STEEP_CURVE.curves();
    let hull_white = HullWhite::init(
        STEEP_CURVE.a,
        STEEP_CURVE.sigma,
        &yield_curve,
        &forward_curve,
    )
    .unwrap();
    let phi_0 = hull_white.phi_t(0.0);
    let phi_5 = hull_white.phi_t(5.0);
    assert!(
        (phi_5 - phi_0).abs() > 0.005,
        "steep fixture must make phi(t) time dependent: phi(0)={phi_0}, phi(5)={phi_5}"
    );
}

/// `r(0)` is the instantaneous forward at the front of the curve, not the cumulative yield.
///
/// The `hw_curves` fixture is a Vasicek short rate starting at `STEEP_CURVE.curr_rate`, so the right
/// answer is known independently of the model code: `short_rate_now()` has to come back with that
/// rate, and `yield_curve(0.0)` -- the tempting wrong answer -- is 0.
#[test]
fn short_rate_now_is_the_instantaneous_forward_at_zero() {
    let (yield_curve, forward_curve) = STEEP_CURVE.curves();
    let hull_white = HullWhite::init(
        STEEP_CURVE.a,
        STEEP_CURVE.sigma,
        &yield_curve,
        &forward_curve,
    )
    .unwrap();
    let r0 = hull_white.short_rate_now().unwrap();
    assert_abs_diff_eq!(r0, STEEP_CURVE.curr_rate, epsilon = 1e-12);
    assert_abs_diff_eq!(r0, forward_curve(0.0), epsilon = 1e-12);
    assert_abs_diff_eq!(r0, hull_white.phi_t(0.0), epsilon = 0.0);
    //...and emphatically not the cumulative yield at 0, which is zero for any integrated curve.
    assert!(
        (r0 - yield_curve(0.0)).abs() > 1e-3,
        "r(0) = {r0} came out at the yield curve's value at 0 ({})",
        yield_curve(0.0)
    );
}

/// The claim the whole `now`/`t` split rests on: `r(0)` is the rate that makes the affine bond
/// price at `t = 0` equal the calibrated `now` price.  If this drifted, every `*_now` variant that
/// routes through `r(0)` would price a different bond than the `*_now` bond pricers do.
#[test]
fn short_rate_now_makes_the_t_and_now_bond_prices_agree() {
    let (yield_curve, forward_curve) = STEEP_CURVE.curves();
    let hull_white = HullWhite::init(
        STEEP_CURVE.a,
        STEEP_CURVE.sigma,
        &yield_curve,
        &forward_curve,
    )
    .unwrap();
    let r0 = hull_white.short_rate_now().unwrap();
    for maturity in [0.0, 0.25, 1.0, 3.0, 7.5, 20.0] {
        let via_t = hull_white.bond_price_t(r0, 0.0, maturity).unwrap();
        let via_now = hull_white.bond_price_now(maturity).unwrap();
        assert_abs_diff_eq!(via_t, via_now, epsilon = 1e-14);
        //The coupon kernel has to agree too, since the Jamshidian legs discount with it.
        let times = [maturity + 0.25, maturity + 0.5, maturity + 0.75];
        let coupon_via_t = hull_white
            .coupon_bond_price_t(r0, 0.0, &times, 0.05)
            .unwrap();
        let coupon_via_now = hull_white.coupon_bond_price_now(&times, 0.05).unwrap();
        assert_abs_diff_eq!(coupon_via_t, coupon_via_now, epsilon = 1e-14);
    }
}

/// A curve with no finite front cannot supply a short rate.  `|t| t.ln()`, the shape used in older
/// examples in this crate, is `-inf` at 0; that surfaces as an error here rather than as a
/// `-inf`-priced instrument somewhere downstream.
#[test]
fn short_rate_now_rejects_a_curve_that_is_not_finite_at_zero() {
    let yield_curve = |t: f64| 0.05 * t;
    let forward_curve = |t: f64| t.ln();
    let hull_white = HullWhite::init(0.2, 0.03, &yield_curve, &forward_curve).unwrap();
    let err = hull_white.short_rate_now().unwrap_err();
    assert!(
        matches!(err, HullWhiteError::NumericalError(_)),
        "expected a numerical error for an -inf front curve, got {err:?}"
    );
    //The same model still prices the `now` zero coupon bond, which never touches a rate.
    assert!(hull_white.bond_price_now(2.0).unwrap().is_finite());
}
