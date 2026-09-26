//! Behaviour of the solve itself: the bracket really brackets, the answer does not depend on the
//! seed, a starved solver says which stage failed, and the configured tolerance bounds the search
//! rather than capping the accuracy of the price.

use crate::HullWhite;
use crate::SolverSettings;
use crate::testutil::{STEEP_A, STEEP_B, STEEP_CURR_RATE, STEEP_SIG, hw_curves};

use super::{FIXTURES, direct_payoff_price};

#[test]
fn the_price_does_not_depend_on_where_the_solver_starts() {
    //The old solver started from a hard-coded 3% and its answer moved with that guess.  With
    //a real bracket and a bisection guard, any seed reaches the same root.
    let times = [2.5, 3.0, 3.5, 4.0];
    let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
    let reference = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve)
        .unwrap()
        .coupon_bond_call_t(0.04, 1.0, 2.0, &times, 0.05, 0.95)
        .unwrap();
    for seed in [-1.0f64, 0.0, 0.03, 0.5, 5.0, 50.0, 1e3, 1e6] {
        let hull_white = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve)
            .unwrap()
            .with_solver(SolverSettings {
                initial_guess: Some(seed),
                ..SolverSettings::default()
            })
            .unwrap();
        let priced = hull_white
            .coupon_bond_call_t(0.04, 1.0, 2.0, &times, 0.05, 0.95)
            .unwrap();
        assert!(
            (priced - reference).abs() <= 1e-10f64.max(reference.abs() * 1e-9),
            "seed {seed}: {priced} vs {reference}"
        );
    }
}

#[test]
fn the_bracket_actually_brackets() {
    //The analytic bracket has to straddle the critical rate: bond value above the strike at
    //the low end, below it at the high end, and the solved root in between.
    let times = [2.5, 3.0, 3.5, 4.0];
    for &(curr, a, b, sigma) in FIXTURES.iter() {
        let (yield_curve, forward_curve) = hw_curves(curr, a, b, sigma);
        let hull_white = HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
        let (r_t, t, u) = (0.04, 1.0, 2.0);
        for strike in [0.5f64, 0.8, 0.95, 1.0, 1.05, 1.3, 3.0] {
            let (lower, upper) = hull_white
                .critical_rate_bracket(u, &times, 0.05, strike)
                .unwrap_or_else(|| panic!("no bracket for {curr},{a},{b},{sigma} strike {strike}"));
            let gap = |rate: f64| {
                hull_white
                    .coupon_bond_price_t(rate, u, &times, 0.05)
                    .unwrap()
                    - strike
            };
            assert!(
                gap(lower) >= 0.0,
                "{curr},{a},{b},{sigma} strike {strike}: low end {lower} gives {}",
                gap(lower)
            );
            assert!(
                gap(upper) <= 0.0,
                "{curr},{a},{b},{sigma} strike {strike}: high end {upper} gives {}",
                gap(upper)
            );
            let solution = hull_white
                .solve_critical_rate(r_t, t, u, &times, 0.05, strike)
                .unwrap();
            assert!(
                solution.root >= lower && solution.root <= upper,
                "{curr},{a},{b},{sigma} strike {strike}: root {} outside [{lower}, {upper}]",
                solution.root
            );
        }
    }
}

#[test]
fn a_single_coupon_schedule_reduces_to_the_zero_coupon_option() {
    //The degenerate schedule: one payment date, carrying coupon plus redemption.  Jamshidian
    //has exactly one leg, so the answer has to be that leg weighted by (1 + coupon) -- which
    //is the plain zero-coupon bond option, priced here without any root find at all.
    //
    //  (1 + c) * call_discount(P(t,T), K / (1 + c), P(t,U), sigma_leg)
    let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
    let hull_white = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve).unwrap();
    let (r_t, t, u) = (0.04, 1.0, 2.0);
    for bond_maturity in [2.5f64, 4.0, 10.0] {
        for coupon_rate in [0.05f64, 0.01, -0.02] {
            let weight = 1.0 + coupon_rate;
            for strike in [0.5f64, 0.9, 0.95, 1.0, 1.3, 3.0] {
                for is_call in [true, false] {
                    let schedule = [bond_maturity];
                    let priced = if is_call {
                        hull_white
                            .coupon_bond_call_t(r_t, t, u, &schedule, coupon_rate, strike)
                            .unwrap()
                    } else {
                        hull_white
                            .coupon_bond_put_t(r_t, t, u, &schedule, coupon_rate, strike)
                            .unwrap()
                    };
                    //The same option through the zero-coupon pricer: one leg, struck at the
                    //strike the single leg would carry.  No decomposition, no solver.
                    let leg = if is_call {
                        hull_white
                            .bond_call_t(r_t, t, u, bond_maturity, strike / weight)
                            .unwrap()
                    } else {
                        hull_white
                            .bond_put_t(r_t, t, u, bond_maturity, strike / weight)
                            .unwrap()
                    };
                    let reference = weight * leg;
                    assert!(
                        (priced - reference).abs() <= 1e-12f64.max(reference.abs() * 1e-11),
                        "T {bond_maturity} c {coupon_rate} {} strike {strike}: {priced} vs {reference}",
                        if is_call { "call" } else { "put" }
                    );
                    //...and the direct integral agrees too.
                    let integral = direct_payoff_price(
                        &hull_white,
                        r_t,
                        t,
                        u,
                        &schedule,
                        coupon_rate,
                        strike,
                        is_call,
                    );
                    assert!(
                        (priced - integral).abs() <= 1e-10f64.max(integral.abs() * 1e-9),
                        "T {bond_maturity} c {coupon_rate} {} strike {strike}: {priced} vs integral {integral}",
                        if is_call { "call" } else { "put" }
                    );
                }
            }
        }
    }
}

#[test]
fn a_forty_eight_coupon_schedule_prices_against_the_integral() {
    //A long schedule is where the solve travels furthest and where the leg sum has the most
    //terms to get wrong: 48 coupons, maturities out to 26 years, strikes from deep ITM to
    //deep OTM.  Still agrees with the direct payoff integral.
    let schedule: Vec<f64> = (1..=48).map(|i| 2.0 + 0.5 * i as f64).collect();
    assert_eq!(schedule.len(), 48);
    for &(curr, a, b, sigma) in FIXTURES.iter() {
        let (yield_curve, forward_curve) = hw_curves(curr, a, b, sigma);
        let hull_white = HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
        let (r_t, t, u) = (0.04, 1.0, 2.0);
        let underlying = hull_white
            .coupon_bond_price_t(r_t, t, &schedule, 0.05)
            .unwrap();
        for factor in [0.2f64, 0.6, 0.9, 0.99, 1.0, 1.01, 1.4, 2.5] {
            let strike = underlying * factor;
            for is_call in [true, false] {
                let priced = if is_call {
                    hull_white
                        .coupon_bond_call_t(r_t, t, u, &schedule, 0.05, strike)
                        .unwrap()
                } else {
                    hull_white
                        .coupon_bond_put_t(r_t, t, u, &schedule, 0.05, strike)
                        .unwrap()
                };
                let reference =
                    direct_payoff_price(&hull_white, r_t, t, u, &schedule, 0.05, strike, is_call);
                assert!(
                    (priced - reference).abs() <= 1e-9f64.max(reference.abs() * 1e-9),
                    "fixture {curr},{a},{b},{sigma} {side} strike {strike}: {priced} vs {reference}",
                    side = if is_call { "call" } else { "put" }
                );
                //Parity on the same instrument, for a second independent angle.
                let other = if is_call {
                    hull_white
                        .coupon_bond_put_t(r_t, t, u, &schedule, 0.05, strike)
                        .unwrap()
                } else {
                    hull_white
                        .coupon_bond_call_t(r_t, t, u, &schedule, 0.05, strike)
                        .unwrap()
                };
                let parity = underlying - strike * hull_white.bond_price_t(r_t, t, u).unwrap();
                //C - P = parity, so a put's difference runs the other way.
                let expected = if is_call { parity } else { -parity };
                assert!(
                    (priced - other - expected).abs() <= 1e-9f64.max(parity.abs() * 1e-10),
                    "fixture {curr},{a},{b},{sigma} strike {strike}: {} vs {}",
                    priced - other,
                    expected
                );
            }
        }
    }
}

#[test]
fn a_starved_solver_says_which_stage_failed() {
    //A failure names the stage that failed instead of handing back a bare number: too few
    //iterations to converge, in the context of the instrument being priced.
    let times = [2.5, 3.0, 3.5, 4.0];
    let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
    let starved = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve)
        .unwrap()
        .with_solver(SolverSettings {
            max_iterations: 1,
            ..SolverSettings::default()
        })
        .unwrap();
    let error = starved
        .coupon_bond_call_t(0.04, 1.0, 2.0, &times, 0.05, 0.95)
        .unwrap_err();
    let message = error.to_string();
    assert!(
        message.contains("iterations"),
        "expected an iteration failure, got {message}"
    );
    assert!(
        message.contains("Jamshidian"),
        "expected instrument context, got {message}"
    );
}

/// Worst price error over a few strikes, against a 1e-15-tolerance reference, for one solve
/// tolerance.
fn worst_price_err(
    tolerance: f64,
    yield_curve: &(impl Fn(f64) -> f64 + Sync),
    forward_curve: &(impl Fn(f64) -> f64 + Sync),
    times: &[f64],
) -> f64 {
    let reference_model = HullWhite::init(STEEP_A, STEEP_SIG, yield_curve, forward_curve)
        .unwrap()
        .with_solver(SolverSettings {
            tolerance: 1e-15,
            max_iterations: 500,
            initial_guess: None,
        })
        .unwrap();
    let model = HullWhite::init(STEEP_A, STEEP_SIG, yield_curve, forward_curve)
        .unwrap()
        .with_solver(SolverSettings {
            tolerance,
            ..SolverSettings::default()
        })
        .unwrap();
    [0.8f64, 0.95, 1.0, 1.05]
        .into_iter()
        .map(|strike| {
            let reference =
                direct_payoff_price(&reference_model, 0.04, 1.0, 2.0, times, 0.05, strike, true);
            let price = model
                .coupon_bond_call_t(0.04, 1.0, 2.0, times, 0.05, strike)
                .unwrap();
            (price - reference).abs()
        })
        .fold(0.0f64, f64::max)
}

#[test]
fn a_looser_root_tolerance_is_visible_in_the_price() {
    //The reason the tolerance is configurable at all: the old hard-coded 1e-7 capped how
    //accurate a price anyone could get, and now that cap is measurable instead of invisible.
    let times = [2.5, 3.0, 3.5, 4.0];
    let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
    let loose = worst_price_err(1e-4, &yield_curve, &forward_curve, &times);
    let tight = worst_price_err(1e-14, &yield_curve, &forward_curve, &times);
    println!("tolerance 1e-4 worst price err {loose:e}; 1e-14 worst price err {tight:e}");
    assert!(tight < loose, "tight {tight:e} should beat loose {loose:e}");
    assert!(tight < 1e-11, "tight tolerance left {tight:e}");
}

#[test]
fn the_root_tolerance_bounds_the_search_not_the_achieved_error() {
    //The tolerance is a statement about how hard the solver keeps hunting (bracket width in rate
    //units), not about how far off the answer it returns.  When the solve exits on the Newton
    //correction the returned root is Newton's *prediction* of the root, which is good to far
    //better than the correction's size, so the old 1e-7 cap is no longer a cap on accuracy:
    //1e-7 now lands on the same price as 1e-14 (both sit on the model's round-off floor).
    let times = [2.5, 3.0, 3.5, 4.0];
    let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
    let legacy_cap = worst_price_err(1e-7, &yield_curve, &forward_curve, &times);
    let tight = worst_price_err(1e-14, &yield_curve, &forward_curve, &times);
    println!("tolerance 1e-7 worst price err {legacy_cap:e}; 1e-14 worst price err {tight:e}");
    //A 1e-7 rate error on this instrument is worth ~1e-9 of price; the achieved error is three
    //orders below that, so the legacy tolerance no longer costs accuracy -- only iterations.
    assert!(
        legacy_cap < 1e-11 && legacy_cap <= tight * 1.5 + 1e-15,
        "tolerance 1e-7 left {legacy_cap:e}, tight 1e-14 left {tight:e}"
    );
}
