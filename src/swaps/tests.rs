//! Unit tests for swap pricing: pricing at the forward swap rate is worth zero, the anchored and
//! derived forms of the swap price agree, and a whole-but-not-binary schedule anchors at `t`.

use approx::*;

use crate::HullWhite;
use crate::testutil::{STEEP_A, STEEP_B, STEEP_CURR_RATE, STEEP_SIG, hw_curves};

/// The float comparison was not cosmetic: `swap_price_t` selects `swap_start` from `is_exact`, so a
/// whole-but-not-binary schedule priced the swap from the wrong anchor.  With an integral number of
/// periods remaining, `swap_price_t` must be the same call as `swap_price_t_init(..., t, n, ...).unwrap()`.
#[test]
fn swap_price_t_anchors_at_t_for_whole_non_binary_schedules() {
    //Needs a curve that moves: on a flat curve a one-period anchor shift is nearly free, which is
    //why the legacy fixture never showed this.
    let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
    let hull_white = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve).unwrap();
    let r_t = STEEP_CURR_RATE;
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
    let curr_rate = 0.02;
    let sig: f64 = 0.02;
    let a: f64 = 0.3;
    let b = 0.04;
    let delta = 0.25;
    let future_time = 0.5;
    let swap_maturity = 5.5;
    let num_swap_payments = 20;
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
    let curr_rate = 0.02;
    let sig: f64 = 0.02;
    let a: f64 = 0.3;
    let b = 0.04;
    let delta = 0.25;
    let future_time = 0.5;
    let swap_maturity = 5.5;
    let num_swap_payments = 20;
    let swap_rate = 0.03;
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
