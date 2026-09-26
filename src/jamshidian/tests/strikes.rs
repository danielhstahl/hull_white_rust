//! Strike-edge behaviour of the decomposition: zero strike, put-call parity, and strikes so far
//! out that the bracket degenerates but the price is still reported rather than refused.

use crate::HullWhite;
use crate::testutil::{STEEP_A, STEEP_B, STEEP_CURR_RATE, STEEP_SIG, hw_curves};

use super::FIXTURES;

#[test]
fn a_zero_strike_call_is_the_underlying_and_a_zero_strike_put_is_worthless() {
    //A zero strike puts the critical rate at the far end of the bond's range, which used to
    //be a guaranteed RootFindingError.  What the answer has to be: the call is the present
    //value of the underlying, the put is nothing.
    let times = [2.5, 3.0, 3.5, 4.0];
    for &(curr, a, b, sigma) in FIXTURES.iter() {
        let (yield_curve, forward_curve) = hw_curves(curr, a, b, sigma);
        let hull_white = HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
        let underlying = hull_white
            .coupon_bond_price_t(0.04, 1.0, &times, 0.05)
            .unwrap();
        let call = hull_white
            .coupon_bond_call_t(0.04, 1.0, 2.0, &times, 0.05, 0.0)
            .unwrap();
        let put = hull_white
            .coupon_bond_put_t(0.04, 1.0, 2.0, &times, 0.05, 0.0)
            .unwrap();
        assert!(
            (call - underlying).abs() <= underlying.abs() * 1e-12,
            "fixture {curr},{a},{b},{sigma}: call {call} vs underlying {underlying}"
        );
        assert_eq!(put, 0.0, "fixture {curr},{a},{b},{sigma}: put {put}");
    }
}

#[test]
fn put_call_parity_holds_whatever_the_strike() {
    //C - P = PV(underlying) - K * P(t,U), to machine precision, including at zero strike.
    //A solve that lands on the wrong critical rate breaks this immediately.
    let times = [2.5, 3.0, 3.5, 4.0];
    let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
    let hull_white = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve).unwrap();
    let (r_t, t, u) = (0.04, 1.0, 2.0);
    let underlying = hull_white
        .coupon_bond_price_t(r_t, t, &times, 0.05)
        .unwrap();
    let discount = hull_white.bond_price_t(r_t, t, u).unwrap();
    for strike in [0.0f64, 0.3, 0.7, 0.95, 1.0, 1.2, 2.0, 1e3, 1e6] {
        let call = hull_white
            .coupon_bond_call_t(r_t, t, u, &times, 0.05, strike)
            .unwrap();
        let put = hull_white
            .coupon_bond_put_t(r_t, t, u, &times, 0.05, strike)
            .unwrap();
        let parity = underlying - strike * discount;
        //At a strike of 1e6 the put is a difference of huge cancelling leg values, so the
        //residual is a few parts in 1e11 of the numbers involved rather of the option.
        assert!(
            (call - put - parity).abs() <= 1e-11f64.max(parity.abs() * 1e-10),
            "strike {strike}: C-P {} vs parity {parity}",
            call - put
        );
    }
}

#[test]
fn an_extreme_strike_prices_rather_than_erroring() {
    //A strike orders of magnitude away from anywhere the bond can reach has no optionality
    //left: the answer is zero (or parity), not an error.  This used to come back as
    //RootFindingError("NaN") because the old Newton iterate blew up on the flat of the
    //exponential.
    let times = [2.5, 3.0, 3.5, 4.0];
    let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
    let hull_white = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve).unwrap();
    for strike in [1e3f64, 1e8, 1e12] {
        let call = hull_white
            .coupon_bond_call_t(0.04, 1.0, 2.0, &times, 0.05, strike)
            .unwrap();
        assert_eq!(call, 0.0, "strike {strike} call {call}");
        let put = hull_white
            .coupon_bond_put_t(0.04, 1.0, 2.0, &times, 0.05, strike)
            .unwrap();
        assert!(put.is_finite() && put > 0.0, "strike {strike} put {put}");
    }
}
