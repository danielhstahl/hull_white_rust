//! Unit tests for the bond-option pricers, against an external reference value
//! (bondoption_vasicek.html on quantcalc.net) and against the coupon-bond decomposition
//! degenerated to a single zero-coupon leg.

use approx::*;

use crate::HullWhite;
use crate::error::HullWhiteError;
use crate::test_support::{STEEP_CURVE, VASICEK_REFERENCE};

#[test]
fn zero_coupon_reference() {
    //http://www.quantcalc.net/BondOption_Vasicek.html
    let fixture = VASICEK_REFERENCE;
    let curr_rate = fixture.curr_rate;
    let strike = 0.96;
    let future_time = 0.0;
    let bond_maturity = 3.0;
    let option_maturity = 2.0;
    let curve = fixture.curve();
    let hull_white = HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();
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
    let fixture = VASICEK_REFERENCE;
    let curr_rate = fixture.curr_rate;
    let strike = 0.96;
    let future_time = 0.0;
    let bond_maturity = 3.0;
    let option_maturity = 2.0;
    let curve = fixture.curve();
    let hull_white = HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();
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

/// The `now` twins of the Jamshidian entry points are the `t` form at `(r(0), 0)`, exactly.
///
/// This is the property the crate docs claim for every `_now` variant of a state-dependent
/// instrument: `r(0)` from `short_rate_now` is the state, and nothing else about the price changes.
#[test]
fn coupon_bond_option_now_matches_the_t_form_at_zero() {
    let curve = STEEP_CURVE.curve();
    let hull_white = HullWhite::new(STEEP_CURVE.a, STEEP_CURVE.sigma, &curve).unwrap();
    let r0 = hull_white.short_rate_now().unwrap();
    let coupon_times = [1.75, 2.0, 2.25, 2.5];
    let coupon_rate = 0.05;
    for option_maturity in [1.0, 1.5, 2.25] {
        for strike in [0.85, 1.0, 1.15] {
            let call_now = hull_white
                .coupon_bond_call_now(option_maturity, &coupon_times, coupon_rate, strike)
                .unwrap();
            let call_t = hull_white
                .coupon_bond_call_t(r0, 0.0, option_maturity, &coupon_times, coupon_rate, strike)
                .unwrap();
            assert_abs_diff_eq!(call_now, call_t, epsilon = 1e-15);
            let put_now = hull_white
                .coupon_bond_put_now(option_maturity, &coupon_times, coupon_rate, strike)
                .unwrap();
            let put_t = hull_white
                .coupon_bond_put_t(r0, 0.0, option_maturity, &coupon_times, coupon_rate, strike)
                .unwrap();
            assert_abs_diff_eq!(put_now, put_t, epsilon = 1e-15);
        }
    }
}

/// Put-call parity across a whole Jamshidian decomposition, priced from `now`.
///
/// The decomposition replaces a coupon-bond option with a basket of leg options; what it must not do
/// is create or destroy value, so `call - put` has to be the forward bond leg,
/// `P_c - K * P(0, option_maturity)`, with no option value left over.  Each strike used below puts
/// the option on a different side of the money, so a one-sided error in the leg strikes would show.
#[test]
fn coupon_bond_option_parity_at_now() {
    let curve = STEEP_CURVE.curve();
    let hull_white = HullWhite::new(STEEP_CURVE.a, STEEP_CURVE.sigma, &curve).unwrap();
    let coupon_times = [1.75, 2.0, 2.25, 2.5];
    let coupon_rate = 0.05;
    let option_maturity = 1.5;
    let bond = hull_white
        .coupon_bond_price_now(&coupon_times, coupon_rate)
        .unwrap();
    let discount = hull_white.bond_price_now(option_maturity).unwrap();
    for strike in [0.85, 0.95, 1.0, 1.05, 1.2] {
        let call = hull_white
            .coupon_bond_call_now(option_maturity, &coupon_times, coupon_rate, strike)
            .unwrap();
        let put = hull_white
            .coupon_bond_put_now(option_maturity, &coupon_times, coupon_rate, strike)
            .unwrap();
        let leg = bond - strike * discount;
        assert_abs_diff_eq!(call - put, leg, epsilon = 1e-11);
        assert!(call >= 0.0 && put >= 0.0);
    }
    //Deep enough in the money on the call side that the strike is below the cash the bond throws off,
    //the call goes to its exercise value and the put to zero -- parity still holds.
    let deep_strike = 0.0;
    let call = hull_white
        .coupon_bond_call_now(option_maturity, &coupon_times, coupon_rate, deep_strike)
        .unwrap();
    let put = hull_white
        .coupon_bond_put_now(option_maturity, &coupon_times, coupon_rate, deep_strike)
        .unwrap();
    assert_abs_diff_eq!(call, bond, epsilon = 1e-12);
    assert_abs_diff_eq!(put, 0.0, epsilon = 1e-12);
}

/// The straddle conventions survive the trip to `now`: the pre-expiry coupon is not part of the
/// underlying, and the coupon settling on the expiry date is a strike reduction on the residual.
#[test]
fn coupon_bond_option_now_keeps_the_straddle_convention() {
    let curve = STEEP_CURVE.curve();
    let hull_white = HullWhite::new(STEEP_CURVE.a, STEEP_CURVE.sigma, &curve).unwrap();
    let option_maturity = 1.5;
    let coupon_rate = 0.05;
    let full = [1.25, 1.5, 1.75, 2.0, 2.25, 2.5];
    let tail = [1.5, 1.75, 2.0, 2.25, 2.5];
    let residual = [1.75, 2.0, 2.25, 2.5];
    for strike in [0.9, 1.0, 1.1] {
        let call = hull_white
            .coupon_bond_call_now(option_maturity, &full, coupon_rate, strike)
            .unwrap();
        //Dropping the coupon paid strictly before expiry cannot change the price.
        let dropped = hull_white
            .coupon_bond_call_now(option_maturity, &tail, coupon_rate, strike)
            .unwrap();
        assert_abs_diff_eq!(call, dropped, epsilon = 1e-15);
        //The coupon paid exactly on the expiry date is cash there, so on the residual bond it is a
        //strike reduction.
        let folded = hull_white
            .coupon_bond_call_now(
                option_maturity,
                &residual,
                coupon_rate,
                strike - coupon_rate,
            )
            .unwrap();
        assert_abs_diff_eq!(call, folded, epsilon = 1e-15);
    }
    //Nothing deliverable after expiry is not an underlying, on either side.
    assert!(matches!(
        hull_white.coupon_bond_call_now(1.5, &[1.25, 1.5], 0.05, 1.0),
        Err(HullWhiteError::InvalidInput(_))
    ));
    assert!(matches!(
        hull_white.coupon_bond_put_now(1.5, &[1.25, 1.5], 0.05, 1.0),
        Err(HullWhiteError::InvalidInput(_))
    ));
}
