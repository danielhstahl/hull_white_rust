//! The other side of the input contract: nothing that used to be valid got caught by the guards.
//! Walks every public entry point once; the module-specific numeric tests are the real net.

use approx::*;

use crate::schedules::get_coupon_times;
use crate::testutil::yvf_setup;

#[test]
fn valid_instruments_still_price() {
    //Smoke test for the other side of the contract: nothing that used to be valid got caught by
    //the guards.  The 32 pre-existing numerical tests (including the bit-exact swaption pins)
    //are the real regression net; this walks every public entry point once.
    yvf_setup!(hull_white);
    let coupon_times = get_coupon_times(4, 1.0, 0.25).unwrap();
    assert_eq!(coupon_times.len(), 4);
    //bond schedules for the option entry points must sit entirely beyond expiry
    let post_expiry = get_coupon_times(4, 1.5, 0.25).unwrap();

    assert!(
        hull_white
            .t_forward_bond_vol(1.0, 2.0, 3.0)
            .unwrap()
            .is_finite()
    );
    assert!(hull_white.mu_r(0.05, 1.0, 2.0).unwrap().is_finite());
    assert!(hull_white.variance_r(1.0, 2.0).unwrap().is_finite());
    //a bond priced at its own maturity is at par
    assert_abs_diff_eq!(
        hull_white.bond_price_t(0.05, 2.0, 2.0).unwrap(),
        1.0,
        epsilon = 1e-12
    );
    assert!(hull_white.bond_price_now(2.0).unwrap() > 0.0);
    //t == 0 entry points are valid at spot
    assert!(hull_white.bond_price_t(0.05, 0.0, 2.0).unwrap().is_finite());
    assert!(
        hull_white
            .coupon_bond_price_t(0.05, 1.0, &coupon_times, 0.05)
            .unwrap()
            > 0.0
    );
    assert!(
        hull_white
            .coupon_bond_price_now(&coupon_times, 0.05)
            .unwrap()
            > 0.0
    );
    assert!(hull_white.bond_call_t(0.05, 1.0, 1.5, 2.0, 0.98).unwrap() >= 0.0);
    assert!(hull_white.bond_call_now(1.5, 2.0, 0.98).unwrap() >= 0.0);
    assert!(hull_white.bond_put_t(0.05, 1.0, 1.5, 2.0, 0.98).unwrap() >= 0.0);
    assert!(hull_white.bond_put_now(1.5, 2.0, 0.98).unwrap() >= 0.0);
    assert!(
        hull_white
            .coupon_bond_call_t(0.05, 1.0, 1.5, &post_expiry, 0.05, 1.0)
            .unwrap()
            >= 0.0
    );
    //This put is deep out of the money (forward ~0.95 against a strike of 1.0), so the honest
    //answer is ~0; the ~1e-17 residual is floating-point noise and is identical to what the
    //pre-change code returned, so assert the magnitude rather than a non-negative bound.
    assert!(
        hull_white
            .coupon_bond_put_t(0.05, 1.0, 1.5, &post_expiry, 0.05, 1.0)
            .unwrap()
            .abs()
            < 1e-12
    );
    assert!(hull_white.caplet_now(1.5, 0.25, 0.04).unwrap() >= 0.0);
    assert!(hull_white.caplet_t(0.05, 1.0, 1.5, 0.25, 0.04).unwrap() >= 0.0);
    assert!(
        hull_white
            .euro_dollar_future_t(0.05, 1.0, 1.5, 0.25)
            .unwrap()
            .is_finite()
    );
    assert!(
        hull_white
            .euro_dollar_future_now(1.5, 0.25)
            .unwrap()
            .is_finite()
    );
    //spot libor fixing at t is legal (maturity == t)
    assert!(
        hull_white
            .forward_libor_rate_t(0.05, 1.0, 1.0, 0.25)
            .unwrap()
            .is_finite()
    );
    assert!(
        hull_white
            .forward_libor_rate_now(1.5, 0.25)
            .unwrap()
            .is_finite()
    );
    assert!(
        hull_white
            .libor_rate_t(0.05, 1.0, 0.25)
            .unwrap()
            .is_finite()
    );
    assert!(
        hull_white
            .forward_swap_rate_t(0.05, 1.0, 1.5, 8, 0.25)
            .unwrap()
            .is_finite()
    );
    assert!(
        hull_white
            .swap_rate_t(0.05, 1.0, 8, 0.25)
            .unwrap()
            .is_finite()
    );
    assert!(
        hull_white
            .swap_price_t(0.05, 1.0, 3.0, 0.25, 0.04)
            .unwrap()
            .is_finite()
    );
    assert!(
        hull_white
            .swap_price_t_init(0.05, 1.0, 1.0, 8, 0.25, 0.04)
            .unwrap()
            .is_finite()
    );
    assert!(
        hull_white
            .european_payer_swaption_t(0.05, 1.0, 1.5, 8, 0.25, 0.04)
            .unwrap()
            >= 0.0
    );
    assert!(
        hull_white
            .european_receiver_swaption_t(0.05, 1.0, 1.5, 8, 0.25, 0.04)
            .unwrap()
            >= 0.0
    );
    assert!(
        hull_white
            .american_payer_swaption_t(0.05, 1.0, 1.5, 8, 0.25, 0.04, 100)
            .unwrap()
            >= 0.0
    );
    assert!(
        hull_white
            .american_receiver_swaption_t(0.05, 1.0, 1.5, 8, 0.25, 0.04, 100)
            .unwrap()
            >= 0.0
    );
}
