//! Unit tests for the bond pricers: the `now` and `t = 0` forms must agree, a bond at its own
//! maturity is at par, and the coupon kernel does not shift those identities.

use approx::*;

use crate::HullWhite;
use crate::schedules::get_coupon_times;
use crate::test_support::BASELINE;

#[test]
fn test_bond_now_same_as_t_when_t_is_zero() {
    let fixture = BASELINE;
    let curr_rate = fixture.curr_rate;
    let future_time = 0.0;
    let maturity = 1.5;
    let curve = fixture.curve();
    let hull_white = HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();
    let bond_price_now = hull_white.bond_price_now(maturity).unwrap();
    let bond_price_t = hull_white
        .bond_price_t(curr_rate, future_time, maturity)
        .unwrap();
    assert_abs_diff_eq!(bond_price_now, bond_price_t, epsilon = 0.0000001);
}

#[test]
fn test_coupon_bond_now_same_as_t_when_t_is_zero() {
    let fixture = BASELINE;
    let curr_rate = fixture.curr_rate;
    let delta = fixture.delta;
    let future_time = 0.0;
    let curve = fixture.curve();
    let coupon_times = get_coupon_times(6, future_time, delta).unwrap(); //this was 5, but made six since last payment is now included
    let coupon_rate = 0.05 * delta;
    let hull_white = HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();
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
    let fixture = BASELINE;
    let curr_rate = fixture.curr_rate;
    let future_time = 0.5;
    let option_maturity = 1.5;
    let curve = fixture.curve();
    let hull_white = HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();
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
    let fixture = BASELINE;
    let curr_rate = fixture.curr_rate;
    let future_time = 0.5;
    let curve = fixture.curve();
    let hull_white = HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();
    assert_eq!(
        hull_white
            .bond_price_t(curr_rate, future_time, future_time)
            .unwrap(),
        1.0
    );
}
