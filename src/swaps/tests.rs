//! Unit tests for swap pricing: pricing at the forward swap rate is worth zero, the anchored and
//! derived forms of the swap price agree, and a whole-but-not-binary schedule anchors at `t`.

use approx::*;

use crate::HullWhite;
use crate::schedules::get_coupon_times;
use crate::test_support::{BASELINE, STEEP_CURVE};

/// The float comparison was not cosmetic: `swap_price_t` selects `swap_start` from `is_exact`, so a
/// whole-but-not-binary schedule priced the swap from the wrong anchor.  With an integral number of
/// periods remaining, `swap_price_t` must be the same call as `swap_price_t_init(..., t, n, ...).unwrap()`.
#[test]
fn swap_price_t_anchors_at_t_for_whole_non_binary_schedules() {
    //Needs a curve that moves: on a flat curve a one-period anchor shift is nearly free, which is
    //why the legacy fixture never showed this.
    let curve = STEEP_CURVE.curve();
    let hull_white = HullWhite::new(STEEP_CURVE.a, STEEP_CURVE.sigma, &curve).unwrap();
    let r_t = STEEP_CURVE.curr_rate;
    let swap_rate = 0.045;
    for (t, swap_maturity, delta, num_payments) in [
        (0.1, 1.1, 0.1, 10),
        (0.1, 2.0, 0.1, 19),
        (0.3, 2.3, 0.2, 10),
        (0.25, 2.25, 0.5, 4),
    ] {
        let derived = hull_white
            .swap_price_t(r_t, t, swap_maturity, delta, swap_rate)
            .unwrap();
        let anchored_at_t = hull_white
            .swap_price_t_init(r_t, t, t, num_payments, delta, swap_rate)
            .unwrap();
        let anchored_one_period_late = hull_white
            .swap_price_t_init(r_t, t, t + delta, num_payments, delta, swap_rate)
            .unwrap();
        assert_abs_diff_eq!(derived, anchored_at_t, epsilon = 1e-12);
        //Sanity: the wrong anchor really is a materially different number, so the assertion above
        //has teeth on this fixture (measured gap is ~1e-2 on a ~1e-1 swap price).
        assert!(
            (anchored_at_t - anchored_one_period_late).abs() > 1e-4,
            "fixture too flat to detect an anchor shift at t={t}: {anchored_at_t} vs \
                 {anchored_one_period_late}"
        );
    }
}

#[test]
fn test_swap() {
    let fixture = BASELINE;
    let curr_rate = fixture.curr_rate;
    let delta = fixture.delta;
    let future_time = 0.5;
    let swap_maturity = 5.5;
    let num_swap_payments = 20;
    let curve = fixture.curve();
    let hull_white = HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();
    assert_abs_diff_eq!(
        hull_white
            .swap_price_t(
                curr_rate,
                future_time,
                swap_maturity,
                delta,
                hull_white
                    .swap_rate_t(curr_rate, future_time, num_swap_payments, delta)
                    .unwrap()
            )
            .unwrap(),
        0.0,
        epsilon = 0.000000001
    );
}

#[test]
fn test_swap_init() {
    let fixture = BASELINE;
    let curr_rate = fixture.curr_rate;
    let delta = fixture.delta;
    let future_time = 0.5;
    let swap_maturity = 5.5;
    let num_swap_payments = 20;
    let swap_rate = 0.03;
    let curve = fixture.curve();
    let hull_white = HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();
    let sp_init = hull_white
        .swap_price_t_init(
            curr_rate,
            future_time,
            future_time,
            num_swap_payments,
            delta,
            swap_rate,
        )
        .unwrap();
    let sp = hull_white
        .swap_price_t(curr_rate, future_time, swap_maturity, delta, swap_rate)
        .unwrap();
    assert_eq!(sp_init, sp);
}

/// Every swap `*_now` variant is the `_t` form at `(r(0), 0)` — the state that a calibration leaves
/// you with, rather than a rate the caller has to invent.
#[test]
fn swap_now_variants_match_the_t_form_at_zero() {
    let curve = STEEP_CURVE.curve();
    let hull_white = HullWhite::new(STEEP_CURVE.a, STEEP_CURVE.sigma, &curve).unwrap();
    let r0 = hull_white.short_rate_now().unwrap();
    let delta = 0.25;
    let num_payments = 16;
    let swap_rate = 0.045;

    assert_abs_diff_eq!(
        hull_white
            .forward_swap_rate_now(1.5, num_payments, delta)
            .unwrap(),
        hull_white
            .forward_swap_rate_t(r0, 0.0, 1.5, num_payments, delta)
            .unwrap(),
        epsilon = 1e-15
    );
    assert_abs_diff_eq!(
        hull_white.swap_rate_now(num_payments, delta).unwrap(),
        hull_white
            .swap_rate_t(r0, 0.0, num_payments, delta)
            .unwrap(),
        epsilon = 1e-15
    );
    //A swap already part-run: 4 payments left, maturing at 1.0.
    assert_abs_diff_eq!(
        hull_white.swap_price_now(1.0, delta, swap_rate).unwrap(),
        hull_white
            .swap_price_t(r0, 0.0, 1.0, delta, swap_rate)
            .unwrap(),
        epsilon = 1e-15
    );
    assert_abs_diff_eq!(
        hull_white
            .swap_price_now_init(num_payments, delta, swap_rate)
            .unwrap(),
        hull_white
            .swap_price_t_init(r0, 0.0, 0.0, num_payments, delta, swap_rate)
            .unwrap(),
        epsilon = 1e-15
    );
    assert_abs_diff_eq!(
        hull_white
            .european_payer_swaption_now(2.0, num_payments, delta, swap_rate)
            .unwrap(),
        hull_white
            .european_payer_swaption_t(r0, 0.0, 2.0, num_payments, delta, swap_rate)
            .unwrap(),
        epsilon = 1e-15
    );
    assert_abs_diff_eq!(
        hull_white
            .european_receiver_swaption_now(2.0, num_payments, delta, swap_rate)
            .unwrap(),
        hull_white
            .european_receiver_swaption_t(r0, 0.0, 2.0, num_payments, delta, swap_rate)
            .unwrap(),
        epsilon = 1e-15
    );
}

/// Pricing today's forward swap at today's forward rate has to leave nothing on the table, on either
/// of the two ways of stating the instrument (derived maturity, or anchored payment count).
#[test]
fn swaps_are_at_par_at_the_now_forward_rate() {
    let curve = STEEP_CURVE.curve();
    let hull_white = HullWhite::new(STEEP_CURVE.a, STEEP_CURVE.sigma, &curve).unwrap();
    let delta = 0.25;
    for num_payments in [1usize, 4, 12, 20] {
        let forward = hull_white.swap_rate_now(num_payments, delta).unwrap();
        let maturity = num_payments as f64 * delta;
        assert_abs_diff_eq!(
            hull_white.swap_price_now(maturity, delta, forward).unwrap(),
            0.0,
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            hull_white
                .swap_price_now_init(num_payments, delta, forward)
                .unwrap(),
            0.0,
            epsilon = 1e-12
        );
        //Moving the fixed rate off the forward moves the price off zero, in the payer's favour when
        //the fixed rate falls.
        assert!(
            hull_white
                .swap_price_now(maturity, delta, forward - 0.01)
                .unwrap()
                > 0.0
        );
        assert!(
            hull_white
                .swap_price_now(maturity, delta, forward + 0.01)
                .unwrap()
                < 0.0
        );
    }
}

/// Payer and receiver swaptions priced from `now` differ by the swap leg they are written on, not by
/// any option value: `receiver - payer = fixed_leg - P(0, option_maturity)`, which collapses to
/// zero exactly when the strike is the forward swap rate.
#[test]
fn swaption_payer_receiver_parity_at_now() {
    let curve = STEEP_CURVE.curve();
    let hull_white = HullWhite::new(STEEP_CURVE.a, STEEP_CURVE.sigma, &curve).unwrap();
    let delta = 0.25;
    let option_maturity = 2.0;
    let num_payments = 16;
    let coupon_times = get_coupon_times(num_payments, option_maturity, delta).unwrap();
    let discount = hull_white.bond_price_now(option_maturity).unwrap();
    let atm = hull_white
        .forward_swap_rate_now(option_maturity, num_payments, delta)
        .unwrap();
    for strike in [atm - 0.02, atm, atm + 0.02, 0.03, 0.06] {
        let payer = hull_white
            .european_payer_swaption_now(option_maturity, num_payments, delta, strike)
            .unwrap();
        let receiver = hull_white
            .european_receiver_swaption_now(option_maturity, num_payments, delta, strike)
            .unwrap();
        let fixed_leg = hull_white
            .coupon_bond_price_now(&coupon_times, strike * delta)
            .unwrap();
        assert_abs_diff_eq!(
            receiver - payer,
            fixed_leg - 1.0 * discount,
            epsilon = 1e-11
        );
        assert!(payer > 0.0 && receiver > 0.0);
    }
    //At the forward strike the fixed leg and the discounted par coincide, so the two swaptions are
    //worth the same: the parity above is not hiding a constant offset.
    let payer_atm = hull_white
        .european_payer_swaption_now(option_maturity, num_payments, delta, atm)
        .unwrap();
    let receiver_atm = hull_white
        .european_receiver_swaption_now(option_maturity, num_payments, delta, atm)
        .unwrap();
    assert_abs_diff_eq!(payer_atm, receiver_atm, epsilon = 1e-11);
}
