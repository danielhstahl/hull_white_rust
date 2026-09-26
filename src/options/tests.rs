//! Unit tests for the bond-option pricers, against an external reference value
//! (bondoption_vasicek.html on quantcalc.net) and against the coupon-bond decomposition
//! degenerated to a single zero-coupon leg.

use approx::*;

use crate::HullWhite;

#[test]
fn zero_coupon_reference() {
    //http://www.quantcalc.net/BondOption_Vasicek.html
    let curr_rate = 0.01;
    let sig: f64 = 0.03;
    let a = 0.05;
    let b = 0.04;
    let strike = 0.96;
    let future_time = 0.0;
    let bond_maturity = 3.0;
    let option_maturity = 2.0;
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
    let bond_call = hull_white
        .bond_call_t(
            curr_rate,
            future_time,
            option_maturity,
            bond_maturity,
            strike,
        )
        .unwrap();
    assert_abs_diff_eq!(bond_call, 0.033282, epsilon = 0.0001)
}

#[test]
fn zero_coupon_to_coupon() {
    let curr_rate = 0.01;
    let sig: f64 = 0.03;
    let a = 0.05;
    let b = 0.04;
    let strike = 0.96;
    let future_time = 0.0;
    let bond_maturity = 3.0;
    let option_maturity = 2.0;
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
    let bond_call = hull_white
        .bond_call_t(
            curr_rate,
            future_time,
            option_maturity,
            bond_maturity,
            strike,
        )
        .unwrap();
    let coupon_rate = 0.0;
    let coupon_bond_call = hull_white
        .coupon_bond_call_t(
            curr_rate,
            future_time,
            option_maturity,
            &[bond_maturity],
            coupon_rate,
            strike,
        )
        .unwrap();

    assert_abs_diff_eq!(bond_call, coupon_bond_call, epsilon = 0.0001)
}
