//! Unit tests for the model primitives.
//!
//! `steep_fixture_actually_makes_phi_time_dependent` guards the fixture the time-coordinate
//! regressions elsewhere rely on: if `phi` were constant, a wrong clock would cancel out and
//! those tests would pass for the wrong reason.

use approx::*;

use crate::HullWhite;
use crate::curves::YieldCurve;
use crate::error::HullWhiteError;
use crate::test_support::STEEP_CURVE;

#[test]
fn steep_fixture_actually_makes_phi_time_dependent() {
    //Protects the tests above: if phi were constant, a time-coordinate error would cancel out
    //and the regression tests would pass for the wrong reason.
    let curve = STEEP_CURVE.curve();
    let hull_white = HullWhite::new(STEEP_CURVE.a, STEEP_CURVE.sigma, &curve).unwrap();
    let phi_0 = hull_white.phi_t(0.0);
    let phi_5 = hull_white.phi_t(5.0);
    assert!(
        (phi_5 - phi_0).abs() > 0.005,
        "steep fixture must make phi(t) time dependent: phi(0)={phi_0}, phi(5)={phi_5}"
    );
}

/// `r(0)` is the instantaneous forward at the front of the curve, not the cumulative yield.
///
/// The `HwCurve` fixture is a Vasicek short rate starting at `STEEP_CURVE.curr_rate`, so the right
/// answer is known independently of the model code: `short_rate_now()` has to come back with that
/// rate, and `zero_yield(0.0)` -- the tempting wrong answer -- is 0.
#[test]
fn short_rate_now_is_the_instantaneous_forward_at_zero() {
    let curve = STEEP_CURVE.curve();
    let hull_white = HullWhite::new(STEEP_CURVE.a, STEEP_CURVE.sigma, &curve).unwrap();
    let r0 = hull_white.short_rate_now().unwrap();
    assert_abs_diff_eq!(r0, STEEP_CURVE.curr_rate, epsilon = 1e-12);
    assert_abs_diff_eq!(r0, curve.forward(0.0), epsilon = 1e-12);
    assert_abs_diff_eq!(r0, hull_white.phi_t(0.0), epsilon = 0.0);
    //...and emphatically not the cumulative yield at 0, which is zero for any integrated curve.
    assert!(
        (r0 - curve.zero_yield(0.0)).abs() > 1e-3,
        "r(0) = {r0} came out at the yield curve's value at 0 ({})",
        curve.zero_yield(0.0)
    );
}

/// The claim the whole `now`/`t` split rests on: `r(0)` is the rate that makes the affine bond
/// price at `t = 0` equal the calibrated `now` price.  If this drifted, every `*_now` variant that
/// routes through `r(0)` would price a different bond than the `*_now` bond pricers do.
#[test]
fn short_rate_now_makes_the_t_and_now_bond_prices_agree() {
    let curve = STEEP_CURVE.curve();
    let hull_white = HullWhite::new(STEEP_CURVE.a, STEEP_CURVE.sigma, &curve).unwrap();
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

/// A curve with a *log-singular front*: perfectly self-consistent, no finite `f(0, 0)`.
///
/// `y(t) = t ln t - t` differentiates to `f(0,t) = ln t`, which is `-inf` at `0` and finite
/// everywhere after it.  This is the case the two checks deliberately split: the construction-time
/// check looks at `t > 0` only ([`CURVE_PROBE_TIMES`](crate::CURVE_PROBE_TIMES)), so this curve
/// builds, and the missing front is reported where it is actually needed — by
/// [`HullWhite::short_rate_now`], not by a `-inf`-priced instrument somewhere downstream.
struct LogFrontYield;

impl YieldCurve for LogFrontYield {
    fn zero_yield(&self, t: f64) -> f64 {
        if t <= 0.0 {
            0.0 // lim_{t->0+} t ln t = 0
        } else {
            t * t.ln() - t
        }
    }
    fn forward(&self, t: f64) -> f64 {
        if t <= 0.0 { f64::NEG_INFINITY } else { t.ln() }
    }
}

/// A curve with no finite front cannot supply a short rate.
///
/// A `|t| t.ln()`-shaped forward — the shape older examples in this crate used — is `-inf` at 0.
/// That is *not* an inconsistency (it is the derivative of the curve's own yield), so the model
/// builds; what it cannot do is answer `r(0)`, and that surfaces here rather than as a
/// `-inf`-priced instrument somewhere downstream.
#[test]
fn short_rate_now_rejects_a_curve_that_is_not_finite_at_zero() {
    let hull_white = HullWhite::new(0.2, 0.03, &LogFrontYield).unwrap();
    let err = hull_white.short_rate_now().unwrap_err();
    assert!(
        matches!(err, HullWhiteError::NumericalError(_)),
        "expected a numerical error for an -inf front curve, got {err:?}"
    );
    //The same model still prices the `now` zero coupon bond, which never touches a rate.
    assert!(hull_white.bond_price_now(2.0).unwrap().is_finite());
    //...and still prices a `t > 0` instrument, where the front is not needed either.
    assert!(hull_white.phi_t(0.25).is_finite());
}

/// The identity `t_forward_bond_vol`'s doc block claims, checked rather than asserted: the
/// number it returns is the *total* log-variance of the deliverable bond under the
/// option-expiry forward measure, i.e. the integral of that bond's instantaneous
/// forward-measure volatility `sigma * (B(s, T_b) - B(s, T))` from `s = t` to `s = T`.
///
/// Integrated here against [`HullWhite::bond_b`] rather than against another closed form, on
/// purpose.  The parameter names say which of the two dates is the expiry and which is the
/// deliverable; this says it again in numbers.  Swap `option_maturity` and `bond_maturity` in
/// the pricer, or read the doc backwards as it used to be written, and the two stop agreeing:
/// neither `B(s, .)` nor the integral's upper limit is symmetric in them.
#[test]
fn t_forward_bond_vol_is_the_integrated_forward_measure_variance() {
    let curve = STEEP_CURVE.curve();
    let sigma = STEEP_CURVE.sigma;
    let hull_white = HullWhite::new(STEEP_CURVE.a, sigma, &curve).unwrap();
    for &(t, option_maturity, bond_maturity) in &[
        (0.0, 1.0, 2.0),
        (0.25, 0.75, 10.0),
        (0.5, 2.0, 5.0),
        (1.0, 3.0, 3.25),
    ] {
        //Simpson's rule over the option's life; `steps` is even.
        let steps = 4000usize;
        let h = (option_maturity - t) / steps as f64;
        let squared_vol = |s: f64| {
            let instantaneous = sigma
                * (hull_white.bond_b(s, bond_maturity) - hull_white.bond_b(s, option_maturity));
            instantaneous * instantaneous
        };
        let mut sum = squared_vol(t) + squared_vol(option_maturity);
        for i in 1..steps {
            let s = t + i as f64 * h;
            let weight = if i % 2 == 1 { 4.0 } else { 2.0 };
            sum += weight * squared_vol(s);
        }
        let integrated_variance = sum * h / 3.0;
        let vol = hull_white
            .t_forward_bond_vol(t, option_maturity, bond_maturity)
            .unwrap();
        assert_relative_eq!(
            vol * vol,
            integrated_variance,
            max_relative = 1e-8,
            epsilon = 1e-18
        );
    }
}

/// The two affine coefficients are the price itself, not decoration beside it: reassembling
/// `exp(C - B * r)` from them reproduces [`HullWhite::bond_price_t`] bit for bit (it is the
/// same expression, which is what makes naming the two halves worth doing), and `B` is that
/// price's own duration `-dP/dr / P`, measured here by finite difference instead of taken on
/// trust from the algebra.
#[test]
fn the_affine_coefficients_rebuild_the_bond_price_and_its_duration() {
    let curve = STEEP_CURVE.curve();
    let hull_white = HullWhite::new(STEEP_CURVE.a, STEEP_CURVE.sigma, &curve).unwrap();
    let r = 0.037;
    for &(t, maturity) in &[(0.0, 1.0), (0.5, 2.0), (1.0, 5.0), (2.0, 2.0)] {
        let b = hull_white.bond_b(t, maturity);
        let c = hull_white.bond_c(t, maturity);
        let rebuilt = (c - b * r).exp();
        let price = hull_white.bond_price_t(r, t, maturity).unwrap();
        assert_abs_diff_eq!(rebuilt, price, epsilon = 0.0);
        let h = 1e-6;
        let d_price = (hull_white.bond_price_t(r + h, t, maturity).unwrap()
            - hull_white.bond_price_t(r - h, t, maturity).unwrap())
            / (2.0 * h);
        assert_abs_diff_eq!(-d_price / price, b, epsilon = 1e-8);
    }
    //On its own maturity date the bond is worth par and does not care what the rate is:
    //B(t, t) = 0 and C(t, t) = 0, exactly, not approximately.
    assert_abs_diff_eq!(hull_white.bond_b(2.0, 2.0), 0.0, epsilon = 0.0);
    assert_abs_diff_eq!(hull_white.bond_c(2.0, 2.0), 0.0, epsilon = 0.0);
}
