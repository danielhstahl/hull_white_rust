//! Unit tests for the bond pricers: the `now` and `t = 0` forms must agree, a bond at its own
//! maturity is at par, and the coupon kernel does not shift those identities.

use approx::*;

use crate::HullWhite;
use crate::schedules::get_coupon_times;

#[test]
fn test_bond_now_same_as_t_when_t_is_zero() {
    let curr_rate = 0.02;
    let sig: f64 = 0.02;
    let a: f64 = 0.3;
    let b = 0.04;
    let future_time = 0.0;
    let maturity = 1.5;
    let yield_curve = |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
    let bond_price_now = hull_white.bond_price_now(maturity).unwrap();
    let bond_price_t = hull_white
        .bond_price_t(curr_rate, future_time, maturity)
        .unwrap();
    assert_abs_diff_eq!(bond_price_now, bond_price_t, epsilon = 0.0000001);
}

#[test]
fn test_coupon_bond_now_same_as_t_when_t_is_zero() {
    let curr_rate = 0.02;
    let sig: f64 = 0.02;
    let a: f64 = 0.3;
    let b = 0.04;
    let delta = 0.25;
    let future_time = 0.0;
    let yield_curve = |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let coupon_times = get_coupon_times(6, future_time, delta).unwrap(); //this was 5, but made six since last payment is now included
    let coupon_rate = 0.05 * delta;
    let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
    let bond_price_now = hull_white
        .coupon_bond_price_now(&coupon_times, coupon_rate)
        .unwrap();
    let bond_price_t = hull_white
        .coupon_bond_price_t(curr_rate, future_time, &coupon_times, coupon_rate)
        .unwrap();
    assert_abs_diff_eq!(bond_price_now, bond_price_t, epsilon = 0.0000001);
}

#[test]
fn test_bond_price() {
    let curr_rate = 0.02;
    let sig: f64 = 0.02;
    let a: f64 = 0.3;
    let b = 0.04;
    let future_time = 0.5;
    let option_maturity = 1.5;
    let yield_curve = |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
    assert_eq!(
        hull_white
            .bond_price_t(curr_rate, future_time, option_maturity)
            .unwrap(),
        hull_white
            .bond_price_now(option_maturity - future_time)
            .unwrap()
    );
}

#[test]
fn test_bond_price_at_expiry() {
    let curr_rate = 0.02;
    let sig: f64 = 0.02;
    let a: f64 = 0.3;
    let b = 0.04;
    let future_time = 0.5;
    let yield_curve = |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
    assert_eq!(
        hull_white
            .bond_price_t(curr_rate, future_time, future_time)
            .unwrap(),
        1.0
    );
}
