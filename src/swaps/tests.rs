//! Unit tests for swap pricing: pricing at the forward swap rate is worth zero, the anchored and
//! derived forms of the swap price agree, a whole-but-not-binary schedule anchors at `t`, and the
//! payment convention the shared annuity kernel documents is what the prices are actually built on
//! (`annuity_t_is_the_discounted_weight_of_the_payment_dates_it_documents`,
//! `one_more_payment_is_a_different_swap_not_the_principal_being_counted`).

use approx::*;

use crate::HullWhite;
use crate::schedules::get_coupon_times;
use crate::test_support::{ALL_SCENARIOS, BASELINE, STEEP_CURVE};

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

/// The acceptance test for the annuity extraction, on a **second schedule**: 8 payments at
/// `delta = 0.5` (semi-annual), so `swap_maturity = t + 4.0` and the whole remaining life is a
/// whole number of periods from the valuation date, which is the case `swap_price_t` derives the
/// payment count from.  Walked over every named calibration and four valuation times, and asserted
/// **exactly** zero, not to a tolerance: the price is `P(t,start) - K * A - P(t, start + n*delta)`
/// and the rate is `(P(t,start) - P(t, start + n*delta)) / A`, so at `K =` that rate the product
/// `K * A` reproduces the numerator exactly on this grid.  Zero is the identity here, which is the
/// point of having both functions read the same `annuity_t`.
#[test]
fn forward_swap_rate_is_the_exact_zero_of_the_swap_price_on_the_semiannual_schedule() {
    let num_swap_payments = 8usize;
    let delta = 0.5;
    let bp = 1e-4; //one basis point, for the "the zero is unique" half of the claim
    for scenario in ALL_SCENARIOS {
        let curve = scenario.curve();
        let hull_white = HullWhite::new(scenario.a, scenario.sigma, &curve).unwrap();
        let r_t = scenario.curr_rate;
        for t in [0.0, 0.25, 0.5, 1.0] {
            let swap_maturity = t + num_swap_payments as f64 * delta;
            let forward = hull_white
                .forward_swap_rate_t(r_t, t, t, num_swap_payments, delta)
                .unwrap();
            let annuity = hull_white.annuity_t(r_t, t, t, num_swap_payments, delta);
            assert_eq!(
                hull_white
                    .swap_price_t(r_t, t, swap_maturity, delta, forward)
                    .unwrap(),
                0.0,
                "{}: {num_swap_payments} semiannual payments priced at their own forward                  {forward} are not flat",
                scenario.name
            );
            //And it is the *only* flat rate: the price is linear in the strike with slope
            //`-annuity_t`, so a bp of strike is `annuity * bp` of value, and paying above the
            //forward is a loss to the payer.  This is the annuity showing up as the swap's DV01,
            //which is the other way of checking that the helper sums the right cash flows.
            assert_abs_diff_eq!(
                hull_white
                    .swap_price_t(r_t, t, swap_maturity, delta, forward + bp)
                    .unwrap(),
                -annuity * bp,
                epsilon = 1e-14
            );
            assert_abs_diff_eq!(
                hull_white
                    .swap_price_t(r_t, t, swap_maturity, delta, forward - bp)
                    .unwrap(),
                annuity * bp,
                epsilon = 1e-14
            );
            assert!(
                hull_white
                    .swap_price_t(r_t, t, swap_maturity, delta, forward + bp)
                    .unwrap()
                    < 0.0
            );
        }
    }
}

/// `annuity_t` is the sum its doc comment enumerates, on the dates that comment lists, and the
/// forward rate and the price are the two halves of one balance sheet read off it.
///
/// The doc's worked example: `start = 1.0`, `delta = 0.25`, `n = 4` — payment dates
/// 1.25, 1.50, 1.75, 2.00; every one carries a coupon and **2.00, the last, also carries the
/// principal**.  Priced on a curve that moves (the fixture is `steep_curve`) so a date used twice
/// or a date skipped is not invisible.
#[test]
fn annuity_t_is_the_discounted_weight_of_the_payment_dates_it_documents() {
    let scenario = STEEP_CURVE;
    let curve = scenario.curve();
    let hull_white = HullWhite::new(scenario.a, scenario.sigma, &curve).unwrap();
    let r_t = scenario.curr_rate;
    let t = 0.0;
    let start = 1.0;
    let delta = 0.25;
    let n = 4usize;
    //The four dates the convention generates, spelled out rather than generated, so the helper is
    //checked against the list and not against itself.
    let payment_dates = [1.25, 1.5, 1.75, 2.0];
    let by_hand = delta
        * payment_dates
            .iter()
            .map(|&date| hull_white.bond_price_t(r_t, t, date).unwrap())
            .sum::<f64>();
    let annuity = hull_white.annuity_t(r_t, t, start, n, delta);
    assert_abs_diff_eq!(annuity, by_hand, epsilon = 1e-15);
    //The smallest schedule that exists: one coupon, and the principal on that same date.
    assert_abs_diff_eq!(
        hull_white.annuity_t(r_t, t, start, 1, delta),
        delta * hull_white.bond_price_t(r_t, t, 1.25).unwrap(),
        epsilon = 1e-15
    );
    //One more payment adds exactly that payment's weight -- the sum runs out to `start + n * delta`,
    //and that last date is the maturity the convention names.
    assert_abs_diff_eq!(
        hull_white.annuity_t(r_t, t, start, n + 1, delta),
        annuity + delta * hull_white.bond_price_t(r_t, t, 2.25).unwrap(),
        epsilon = 1e-15
    );
    //The balance the convention is built on, at the forward rate: the fixed leg (its coupons, i.e.
    //K * A, plus the principal on the last payment date) equals the floating leg, which is par
    //discounted from the reset date and so is just `P(t, start)`.
    let forward = hull_white
        .forward_swap_rate_t(r_t, t, start, n, delta)
        .unwrap();
    let principal_at_maturity = hull_white
        .bond_price_t(r_t, t, start + n as f64 * delta)
        .unwrap();
    assert_abs_diff_eq!(
        forward * annuity + principal_at_maturity,
        hull_white.bond_price_t(r_t, t, start).unwrap(),
        epsilon = 1e-15
    );
    //Same statement from the pricer's side.  Note `swap_price_t_init`, not `swap_price_t`: this
    //swap starts at 1.0 and `t = 0`, and `swap_price_t` derives its schedule from *its own* `t`,
    //so the same maturity would describe a swap that had been running since 0 with 8 payments.
    assert_abs_diff_eq!(
        hull_white
            .swap_price_t_init(r_t, t, start, n, delta, forward)
            .unwrap(),
        0.0,
        epsilon = 1e-15
    );
}

/// The `//open question, should num_swap_payments be num_swap_payments+1??` that
/// `swap_price_t_init_raw` used to carry, answered with numbers: the `+1` reading is a different
/// instrument.  `n` counts coupons; the principal rides on the last coupon's date rather than
/// adding one, so `n + 1` means an extra period of coupon *and* the principal paid one period
/// past the maturity.
#[test]
fn one_more_payment_is_a_different_swap_not_the_principal_being_counted() {
    let scenario = STEEP_CURVE;
    let curve = scenario.curve();
    let hull_white = HullWhite::new(scenario.a, scenario.sigma, &curve).unwrap();
    let r_t = scenario.curr_rate;
    let t = 0.0;
    let start = 1.0;
    let delta = 0.25;
    let forward_4 = hull_white
        .forward_swap_rate_t(r_t, t, start, 4, delta)
        .unwrap();
    let at_4 = hull_white
        .swap_price_t_init(r_t, t, start, 4, delta, forward_4)
        .unwrap();
    let at_5 = hull_white
        .swap_price_t_init(r_t, t, start, 5, delta, forward_4)
        .unwrap();
    assert_abs_diff_eq!(at_4, 0.0, epsilon = 1e-15);
    //Read as 5 payments, the same strike leaves the payer holding an extra period's cash flow:
    //the 2.25 coupon (`delta * forward_4 * P(2.25)`) plus principal repaid at 2.25 instead of
    //2.00.  Which is `P(2.0) - P(2.25) * (1 + delta * forward_4)`: 7.0e-4 on this fixture, about
    //twelve orders of magnitude above the 4.4e-16 worst move this refactor made.  A convention
    //error, not rounding.
    let p_2_00 = hull_white.bond_price_t(r_t, t, 2.0).unwrap();
    let p_2_25 = hull_white.bond_price_t(r_t, t, 2.25).unwrap();
    assert_abs_diff_eq!(
        at_5,
        at_4 + p_2_00 - (1.0 + delta * forward_4) * p_2_25,
        epsilon = 1e-15
    );
    assert!(
        (at_5 - at_4).abs() > 1e-4,
        "the 5-payment misreading should be a visible instrument change, got {at_5}"
    );
    //And the 5-payment swap is only flat at the 5-payment forward: the count is part of the
    //instrument, and every count has its own rate that makes it worth zero.
    let forward_5 = hull_white
        .forward_swap_rate_t(r_t, t, start, 5, delta)
        .unwrap();
    assert_abs_diff_eq!(
        hull_white
            .swap_price_t_init(r_t, t, start, 5, delta, forward_5)
            .unwrap(),
        0.0,
        epsilon = 1e-15
    );
    assert!(
        (forward_5 - forward_4).abs() > 1e-6,
        "{forward_5} vs {forward_4}: a longer swap should not quote the identical rate here"
    );
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
