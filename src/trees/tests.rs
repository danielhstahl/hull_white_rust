//! Unit tests for the short-rate tree.
//!
//! The tree is where the shifted-time/absolute-time bug lived, so these tests pin: the payoff
//! helpers, the European tree against the analytic swaption at `t = 0` and `t > 0`, a
//! bit-identical `t = 0` result versus the pre-fix pricers, and a positive early-exercise
//! premium for the American swaptions at `t > 0`.

use approx::*;

use super::{max_or_zero, payoff_swaption};
use crate::HullWhite;
use crate::testutil::{STEEP_A, STEEP_B, STEEP_CURR_RATE, STEEP_SIG, hw_curves};

#[test]
fn test_max_or_zero() {
    let v = 1.0;
    assert_eq!(max_or_zero(v), 1.0);
    assert_eq!(max_or_zero(-v), 0.0);
}

#[test]
fn test_payoff_swaption() {
    let v = 1.0;
    assert_eq!(payoff_swaption(true, v), 1.0);
    assert_eq!(payoff_swaption(false, v), 0.0);
    assert_eq!(payoff_swaption(true, -v), 0.0);
    assert_eq!(payoff_swaption(false, -v), 1.0);
}

#[test]
fn european_swaption_tree_matches_analytic_when_t_is_zero() {
    //Guard: when the valuation time is 0 the shifted tree clock and the absolute model clock
    //coincide, so this pins the (already correct) behaviour across the time-coordinate fix.
    let curr_rate = STEEP_CURR_RATE;
    let sig = STEEP_SIG;
    let a = STEEP_A;
    let b = STEEP_B;
    let delta = 0.25;
    let future_time = 0.0;
    let option_maturity = 1.5;
    let num_swap_payments = 20;
    let (yield_curve, forward_curve) = hw_curves(curr_rate, a, b, sig);
    let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
    let swap_rate = hull_white
        .forward_swap_rate_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
        )
        .unwrap();
    //Tree convergence for this fixture: |tree - analytic| is ~1.4e-4 at 50 steps, ~1.1e-4 at
    //100, ~7.0e-5 at 200, ~3.1e-5 at 400, ~4e-7 at 800.  1e-4 at 400 steps gives >3x headroom
    //over the measured discretisation noise while still being ~40x tighter than the time-shift bug.
    let steps = 400;
    let payer = hull_white
        .european_payer_swaption_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
        )
        .unwrap();
    let tree_payer = hull_white.european_swaption_tree(
        curr_rate,
        future_time,
        option_maturity,
        num_swap_payments,
        delta,
        swap_rate,
        true,
        steps,
    );
    assert_abs_diff_eq!(payer, tree_payer, epsilon = 0.0001);
}

/// Bit-exact guard on the `t = 0` path.  The values below were produced by the pre-fix pricers at
/// commit 254d089 (`option_maturity` 1.5, `num_swap_payments` 20, `delta` 0.25, ATM-forward
/// strike, 400 tree steps), captured with a round-trip-checked `{:?}` print -- note `{:.17}` is
/// 17 digits *after the decimal point*, which for these magnitudes loses the last bit.
/// At `t = 0` the shifted tree clock and the absolute model clock coincide, so the
/// time-coordinate fix must leave every one of them bit-identical -- pinning bits rather than a
/// tolerance means a real behaviour change at `t = 0` cannot slip through unnoticed, and cannot
/// be "absorbed" by widening an epsilon later.
#[test]
fn swaption_tree_at_t_is_bit_identical_to_pre_fix() {
    let option_maturity = 1.5;
    let num_swap_payments = 20;
    let delta = 0.25;
    let steps = 400;
    //             name     r0     a     b     sig   eur payer        eur receiver     amer payer   amer receiver
    let golden: [(&str, f64, f64, f64, f64, f64, f64, f64, f64); 2] = [
        (
            "legacy",
            0.05,
            0.05,
            0.05,
            0.01,
            0.017330477644662997,
            0.017329771203617984,
            0.01834265355592532,
            0.017797483448434452,
        ),
        (
            "steep",
            0.02,
            0.2,
            0.06,
            0.03,
            0.03597512274296273,
            0.03597334589011368,
            0.03781406402902323,
            0.04452822113093644,
        ),
    ];
    for (name, r0, a, b, sig, eur_p, eur_r, amer_p, amer_r) in golden {
        let (yield_curve, forward_curve) = hw_curves(r0, a, b, sig);
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        let swap_rate = hull_white
            .forward_swap_rate_t(r0, 0.0, option_maturity, num_swap_payments, delta)
            .unwrap();
        assert_bits_eq(
            hull_white.european_swaption_tree(
                r0,
                0.0,
                option_maturity,
                num_swap_payments,
                delta,
                swap_rate,
                true,
                steps,
            ),
            eur_p,
            &format!("{name} european payer tree"),
        );
        assert_bits_eq(
            hull_white.european_swaption_tree(
                r0,
                0.0,
                option_maturity,
                num_swap_payments,
                delta,
                swap_rate,
                false,
                steps,
            ),
            eur_r,
            &format!("{name} european receiver tree"),
        );
        assert_bits_eq(
            hull_white
                .american_payer_swaption_t(
                    r0,
                    0.0,
                    option_maturity,
                    num_swap_payments,
                    delta,
                    swap_rate,
                    steps,
                )
                .unwrap(),
            amer_p,
            &format!("{name} american payer"),
        );
        assert_bits_eq(
            hull_white
                .american_receiver_swaption_t(
                    r0,
                    0.0,
                    option_maturity,
                    num_swap_payments,
                    delta,
                    swap_rate,
                    steps,
                )
                .unwrap(),
            amer_r,
            &format!("{name} american receiver"),
        );
    }
}

/// Golden comparison helper: asserts the two f64 have identical bit patterns, so "unchanged" here
/// really means unchanged, not "within a tolerance that can be nudged".
#[allow(clippy::float_cmp)]
fn assert_bits_eq(actual: f64, expected: f64, what: &str) {
    assert_eq!(
        actual.to_bits(),
        expected.to_bits(),
        "{what} moved at t = 0: actual={actual:.17} expected={expected:.17}"
    );
}

#[test]
fn european_swaption_tree_matches_analytic_when_t_is_nonzero() {
    //Regression for the shifted-vs-absolute time bug: the tree runs on the clock
    //`option_maturity - t` but phi() and the swap legs need absolute time from "now" (0).
    let curr_rate = STEEP_CURR_RATE;
    let sig = STEEP_SIG;
    let a = STEEP_A;
    let b = STEEP_B;
    let delta = 0.25;
    let future_time = 0.5;
    let option_maturity = 1.5;
    let num_swap_payments = 20;
    let (yield_curve, forward_curve) = hw_curves(curr_rate, a, b, sig);
    let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
    let swap_rate = hull_white
        .forward_swap_rate_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
        )
        .unwrap();
    //Same budget as the t = 0 guard above.  With the tree clocking the shifted time into phi and
    //into the swap legs, this diff plateaus at ~4.0e-3 (payer) / ~5.0e-3 (receiver) no matter how
    //many steps are used -- a systematic pricing error of ~13-16% of the option value, not
    //discretisation noise.  After the fix the residual is ~2e-5 at 400 steps.
    let steps = 400;
    let payer = hull_white
        .european_payer_swaption_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
        )
        .unwrap();
    let receiver = hull_white
        .european_receiver_swaption_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
        )
        .unwrap();
    let tree_payer = hull_white.european_swaption_tree(
        curr_rate,
        future_time,
        option_maturity,
        num_swap_payments,
        delta,
        swap_rate,
        true,
        steps,
    );
    let tree_receiver = hull_white.european_swaption_tree(
        curr_rate,
        future_time,
        option_maturity,
        num_swap_payments,
        delta,
        swap_rate,
        false,
        steps,
    );
    assert_abs_diff_eq!(payer, tree_payer, epsilon = 0.0001);
    assert_abs_diff_eq!(receiver, tree_receiver, epsilon = 0.0001);
}

#[test]
fn american_swaption_tree_at_nonzero_t_carries_a_positive_early_exercise_premium() {
    //The American tree shares the time-coordinate fix with the European tree; this pins that the
    //American price at t > 0 is a sane number above the European analytic price rather than a
    //shifted-clock artefact.
    let delta = 0.25;
    let future_time = 0.5;
    let option_maturity = 1.5;
    let num_swap_payments = 20;
    let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
    let hull_white = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve).unwrap();
    let swap_rate = hull_white
        .forward_swap_rate_t(
            STEEP_CURR_RATE,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
        )
        .unwrap();
    let american = |is_payer: bool, steps: usize| {
        if is_payer {
            hull_white
                .american_payer_swaption_t(
                    STEEP_CURR_RATE,
                    future_time,
                    option_maturity,
                    num_swap_payments,
                    delta,
                    swap_rate,
                    steps,
                )
                .unwrap()
        } else {
            hull_white
                .american_receiver_swaption_t(
                    STEEP_CURR_RATE,
                    future_time,
                    option_maturity,
                    num_swap_payments,
                    delta,
                    swap_rate,
                    steps,
                )
                .unwrap()
        }
    };
    for is_payer in [true, false] {
        let side = if is_payer { "payer" } else { "receiver" };
        let european = if is_payer {
            hull_white
                .european_payer_swaption_t(
                    STEEP_CURR_RATE,
                    future_time,
                    option_maturity,
                    num_swap_payments,
                    delta,
                    swap_rate,
                )
                .unwrap()
        } else {
            hull_white
                .european_receiver_swaption_t(
                    STEEP_CURR_RATE,
                    future_time,
                    option_maturity,
                    num_swap_payments,
                    delta,
                    swap_rate,
                )
                .unwrap()
        };
        let american_200 = american(is_payer, 200);
        let american_400 = american(is_payer, 400);
        assert!(
            american_200.is_finite() && american_400.is_finite(),
            "{side} swaption: non-finite price: 200={american_200}, 400={american_400}"
        );
        assert!(
            american_400 > european,
            "{side} swaption: early exercise premium must be positive: european={european}, \
                 american(400)={american_400}"
        );
        //Converging, not oscillating: refining the tree moves the price by less than the premium.
        assert!(
            (american_400 - american_200).abs() < (american_400 - european),
            "{side} swaption: american tree not converging: 200={american_200}, \
                 400={american_400}, european={european}"
        );
    }
}

mod legacy_swaptions;
