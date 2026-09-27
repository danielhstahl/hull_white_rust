//! Unit tests for the short-rate tree.
//!
//! The tree is where the shifted-time/absolute-time bug lived, so these tests pin: the payoff
//! helpers, the European tree against the analytic swaption at `t = 0` and `t > 0`, a
//! bit-identical `t = 0` result versus the pre-fix pricers, and a positive early-exercise
//! premium for the American swaptions at `t > 0`.
//!
//! Every tree price here is taken through the **public** `european_*_swaption_tree` /
//! `american_*_swaption_t` entry points rather than through the private [`HullWhite::swaption_tree`]
//! kernel: both sides of the parity checks are then calls a consumer could make, which is the
//! whole reason the European tree helper is public.  The kernel itself has no behaviour beyond
//! what those entry points reach.

use approx::*;

use super::{max_or_zero, payoff_swaption};
use crate::HullWhite;
use crate::test_support::{ALL_SCENARIOS, STEEP_CURVE, hw_curve};

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
    let curr_rate = STEEP_CURVE.curr_rate;
    let sig = STEEP_CURVE.sigma;
    let a = STEEP_CURVE.a;
    let b = STEEP_CURVE.b;
    let delta = 0.25;
    let future_time = 0.0;
    let option_maturity = 1.5;
    let num_swap_payments = 20;
    let curve = hw_curve(curr_rate, a, b, sig);
    let hull_white = HullWhite::new(a, sig, &curve).unwrap();
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
    let tree_payer = hull_white
        .european_payer_swaption_tree(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            steps,
        )
        .unwrap();
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
    /// One row of captured pre-fix values.  Named fields rather than a nine-slot
    /// tuple: the tuple version tripped `clippy::type_complexity`, and an
    /// unlabeled `(f64, f64, f64, f64, f64, f64, f64, f64)` is exactly the shape
    /// in which two values can be swapped without anything noticing.
    #[derive(Clone, Copy)]
    struct GoldenRow {
        name: &'static str,
        r0: f64,
        a: f64,
        b: f64,
        sig: f64,
        eur_payer: f64,
        eur_receiver: f64,
        amer_payer: f64,
        amer_receiver: f64,
    }
    let golden = [
        GoldenRow {
            name: "legacy",
            r0: 0.05,
            a: 0.05,
            b: 0.05,
            sig: 0.01,
            eur_payer: 0.017330477644662997,
            eur_receiver: 0.017329771203617984,
            amer_payer: 0.01834265355592532,
            amer_receiver: 0.017797483448434452,
        },
        GoldenRow {
            name: "steep",
            r0: 0.02,
            a: 0.2,
            b: 0.06,
            sig: 0.03,
            eur_payer: 0.03597512274296273,
            eur_receiver: 0.03597334589011368,
            amer_payer: 0.03781406402902323,
            amer_receiver: 0.04452822113093644,
        },
    ];
    for row in &golden {
        let GoldenRow {
            name,
            r0,
            a,
            b,
            sig,
            eur_payer: eur_p,
            eur_receiver: eur_r,
            amer_payer: amer_p,
            amer_receiver: amer_r,
        } = *row;
        let curve = hw_curve(r0, a, b, sig);
        let hull_white = HullWhite::new(a, sig, &curve).unwrap();
        let swap_rate = hull_white
            .forward_swap_rate_t(r0, 0.0, option_maturity, num_swap_payments, delta)
            .unwrap();
        assert_bits_eq(
            hull_white
                .european_payer_swaption_tree(
                    r0,
                    0.0,
                    option_maturity,
                    num_swap_payments,
                    delta,
                    swap_rate,
                    steps,
                )
                .unwrap(),
            eur_p,
            &format!("{name} european payer tree"),
        );
        assert_bits_eq(
            hull_white
                .european_receiver_swaption_tree(
                    r0,
                    0.0,
                    option_maturity,
                    num_swap_payments,
                    delta,
                    swap_rate,
                    steps,
                )
                .unwrap(),
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
    let curr_rate = STEEP_CURVE.curr_rate;
    let sig = STEEP_CURVE.sigma;
    let a = STEEP_CURVE.a;
    let b = STEEP_CURVE.b;
    let delta = 0.25;
    let future_time = 0.5;
    let option_maturity = 1.5;
    let num_swap_payments = 20;
    let curve = hw_curve(curr_rate, a, b, sig);
    let hull_white = HullWhite::new(a, sig, &curve).unwrap();
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
    let tree_payer = hull_white
        .european_payer_swaption_tree(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            steps,
        )
        .unwrap();
    let tree_receiver = hull_white
        .european_receiver_swaption_tree(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            steps,
        )
        .unwrap();
    assert_abs_diff_eq!(payer, tree_payer, epsilon = 0.0001);
    assert_abs_diff_eq!(receiver, tree_receiver, epsilon = 0.0001);
}

/// The same parity check the two tests above run on one fixture, run over every calibration and
/// every valuation time this crate names.
///
/// Both products come out of the *shared* tree wiring ([`HullWhite::swaption_tree`]), so the
/// cross-check has to hold everywhere the wiring is used, not just at the two points that were
/// written down when the time-coordinate fix was made.  The tolerance is relative rather than
/// absolute because the option value varies by an order of magnitude across the scenarios: a flat
/// `1e-4` would be tight against `low_vol`'s tiny prices and loose against `high_vol`'s big ones,
/// while the tree's error scales with the price.
///
/// Worst measured relative residual at 400 steps, over all 42 (scenario, t, side) points:
/// `7.0e-4` (`high_vol` receiver).  The bound below is `2.5e-3`, ~3.5x that, and ~5 orders
/// tighter than the shifted-clock error it is guarding against (13-16% at the same point).
#[test]
fn european_tree_matches_analytic_across_every_scenario_and_valuation_time() {
    let num_swap_payments = 8;
    let steps = 400;
    //Relative slack on the tree-vs-analytic residual; see the doc comment for the measured worst case.
    let relative_slack = 2.5e-3;
    for scenario in ALL_SCENARIOS {
        let curve = scenario.curve();
        let hull_white = HullWhite::new(scenario.a, scenario.sigma, &curve).unwrap();
        for t in [0.0, 0.5, 1.0] {
            let option_maturity = t + 1.0;
            let swap_rate = hull_white
                .forward_swap_rate_t(
                    scenario.curr_rate,
                    t,
                    option_maturity,
                    num_swap_payments,
                    scenario.delta,
                )
                .unwrap();
            for is_payer in [true, false] {
                let side = if is_payer { "payer" } else { "receiver" };
                let analytic = if is_payer {
                    hull_white
                        .european_payer_swaption_t(
                            scenario.curr_rate,
                            t,
                            option_maturity,
                            num_swap_payments,
                            scenario.delta,
                            swap_rate,
                        )
                        .unwrap()
                } else {
                    hull_white
                        .european_receiver_swaption_t(
                            scenario.curr_rate,
                            t,
                            option_maturity,
                            num_swap_payments,
                            scenario.delta,
                            swap_rate,
                        )
                        .unwrap()
                };
                let tree = if is_payer {
                    hull_white
                        .european_payer_swaption_tree(
                            scenario.curr_rate,
                            t,
                            option_maturity,
                            num_swap_payments,
                            scenario.delta,
                            swap_rate,
                            steps,
                        )
                        .unwrap()
                } else {
                    hull_white
                        .european_receiver_swaption_tree(
                            scenario.curr_rate,
                            t,
                            option_maturity,
                            num_swap_payments,
                            scenario.delta,
                            swap_rate,
                            steps,
                        )
                        .unwrap()
                };
                let relative = (analytic - tree).abs() / analytic;
                assert!(
                    relative < relative_slack,
                    "{} at t = {t}, {side} swaption: tree {tree} vs analytic {analytic} \
                     (relative {relative:.3e}, bound {relative_slack:.1e})",
                    scenario.name,
                );
            }
        }
    }
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
    let curve = STEEP_CURVE.curve();
    let hull_white = HullWhite::new(STEEP_CURVE.a, STEEP_CURVE.sigma, &curve).unwrap();
    let swap_rate = hull_white
        .forward_swap_rate_t(
            STEEP_CURVE.curr_rate,
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
                    STEEP_CURVE.curr_rate,
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
                    STEEP_CURVE.curr_rate,
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
                    STEEP_CURVE.curr_rate,
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
                    STEEP_CURVE.curr_rate,
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

/// The tree entry points validate as one set, because they share one validator.  A call that the
/// American side refuses has to be refused by the European tree side too -- otherwise the two
/// halves of the cross-check accept different instruments and stop being comparable.
#[test]
fn every_tree_entry_point_rejects_the_same_bad_instruments() {
    let curve = STEEP_CURVE.curve();
    let hull_white = HullWhite::new(STEEP_CURVE.a, STEEP_CURVE.sigma, &curve).unwrap();
    //num_swap_payments = 0, on both sides and both styles.
    assert!(
        hull_white
            .european_payer_swaption_tree(0.05, 1.0, 1.5, 0, 0.25, 0.04, 50)
            .is_err()
    );
    assert!(
        hull_white
            .european_receiver_swaption_tree(0.05, 1.0, 1.5, 0, 0.25, 0.04, 50)
            .is_err()
    );
    //num_steps = 0 (the tree would have no nodes).
    assert!(
        hull_white
            .european_payer_swaption_tree(0.05, 1.0, 1.5, 8, 0.25, 0.04, 0)
            .is_err()
    );
    assert!(
        hull_white
            .european_receiver_swaption_tree(0.05, 1.0, 1.5, 8, 0.25, 0.04, 0)
            .is_err()
    );
    //An option that has already expired relative to the valuation time.
    assert!(
        hull_white
            .european_payer_swaption_tree(0.05, 1.5, 1.5, 8, 0.25, 0.04, 50)
            .is_err()
    );
    //A non-finite strike, and a non-positive tenor.
    let nan = f64::NAN;
    assert!(
        hull_white
            .european_receiver_swaption_tree(0.05, 1.0, 1.5, 8, 0.25, nan, 50)
            .is_err()
    );
    assert!(
        hull_white
            .european_payer_swaption_tree(0.05, 1.0, 1.5, 8, 0.0, 0.04, 50)
            .is_err()
    );
    //A negative valuation time.
    assert!(
        hull_white
            .european_receiver_swaption_tree(0.05, -1.0, 1.5, 8, 0.25, 0.04, 50)
            .is_err()
    );
}

mod legacy_swaptions;
