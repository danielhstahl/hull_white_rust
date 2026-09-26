//! Tests for the Jamshidian decomposition, checked against a quadrature of the payoff that
//! shares none of its machinery.
//!
//! [`direct_payoff_price`] integrates `(P_c(r_u) - K)^+` against the option-maturity-forward
//! distribution of the short rate ([`expiry_rate_moments`]), finding the payoff kink by its own
//! bisection so the kink is never inside a Simpson panel.  Agreement between that and the
//! decomposition pins the bracket, the solve and the leg weighting at once.

use crate::HullWhite;
use crate::testutil::{STEEP_A, STEEP_B, STEEP_CURR_RATE, STEEP_SIG, hw_curves};

fn simpson(f: &dyn Fn(f64) -> f64, a: f64, b: f64, panels: usize) -> f64 {
    let panels = panels.max(2) + (panels % 2); //Simpson needs an even panel count
    let h = (b - a) / panels as f64;
    let mut sum = f(a) + f(b);
    for i in 1..panels {
        let x = a + i as f64 * h;
        sum += if i % 2 == 1 { 4.0 } else { 2.0 } * f(x);
    }
    sum * h / 3.0
}

/// Mean and standard deviation of the short rate at `option_maturity`, under the
/// *option-maturity-forward* measure.
///
/// A claim settled on the option's expiry date is priced with the zero-coupon bond maturing
/// then as numeraire, so that -- not the risk-neutral money-market measure -- is the law the
/// option's payoff has to be integrated against.  Under this model `r_U` is Gaussian either
/// way with the same variance, and the drift moves by
/// `sigma^2 * int_t^U exp(-a(U-s)) B(s,U) ds`, which closes to the form below.
fn expiry_rate_moments(
    hull_white: &HullWhite<impl Fn(f64) -> f64 + Sync, impl Fn(f64) -> f64 + Sync>,
    r_t: f64,
    t: f64,
    option_maturity: f64,
) -> (f64, f64) {
    let tau = option_maturity - t;
    let a = hull_white.a;
    let sigma = hull_white.sigma;
    let risk_neutral_mean = hull_white.mu_r(r_t, t, option_maturity).unwrap();
    let variance = hull_white.variance_r(t, option_maturity).unwrap();
    let numeraire_drift = (sigma * sigma / (a * a))
        * ((1.0 - (-a * tau).exp()) - 0.5 * (1.0 - (-2.0 * a * tau).exp()));
    (risk_neutral_mean - numeraire_drift, variance.sqrt())
}

/// Prices the payoff directly by integrating over the expiry-date distribution of the short
/// rate, as an independent check on the Jamshidian machinery.
///
/// ```text
/// price = P(t,U) * E^U[(P_c(mu + sd*Z) - strike)^+]
/// ```
///
/// Nothing in here uses Jamshidian's decomposition, the analytic bracket or the root solver,
/// so agreement pins those down rather than echoing them.  The payoff kinks where
/// `P_c(r) = strike`; that kink is found by this test's own plain bisection and each side is
/// integrated separately, so the kink is never interior to a panel.  `Z` is truncated at
/// +/-12 sigma, where `phi` underflows, putting the tail error far below any price worth
/// quoting.
#[allow(clippy::too_many_arguments)] //a test helper with an instrument's whole description
fn direct_payoff_price(
    hull_white: &HullWhite<impl Fn(f64) -> f64 + Sync, impl Fn(f64) -> f64 + Sync>,
    r_t: f64,
    t: f64,
    option_maturity: f64,
    coupon_times: &[f64],
    coupon_rate: f64,
    strike: f64,
    is_call: bool,
) -> f64 {
    let coupon_bond_at = |rate: f64| {
        hull_white
            .coupon_bond_price_t(rate, option_maturity, coupon_times, coupon_rate)
            .unwrap()
    };
    let (mean, sd) = expiry_rate_moments(hull_white, r_t, t, option_maturity);
    assert!(
        sd > 0.0,
        "the reference needs a non-degenerate expiry distribution"
    );

    //Kink: P_c(r) = strike.  P_c falls strictly with the rate, so widen until the ends
    //disagree and then halve.
    let g = |rate: f64| coupon_bond_at(rate) - strike;
    let mut lo = mean - sd;
    let mut hi = mean + sd;
    while g(lo) < 0.0 {
        lo = mean - 2.0 * (mean - lo);
        assert!(lo > -1e6, "no kink found below the mean");
    }
    while g(hi) > 0.0 {
        hi = mean + 2.0 * (hi - mean);
        assert!(hi < 1e6, "no kink found above the mean");
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if g(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let kink_z = (0.5 * (lo + hi) - mean) / sd;

    let pdf = |z: f64| (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt();
    let integrand = |z: f64| {
        let density = pdf(z);
        if density == 0.0 {
            //Stops inf * 0 = NaN in a deep tail where the payoff itself overflows.
            return 0.0;
        }
        let value = coupon_bond_at(mean + sd * z);
        let payoff = if is_call {
            (value - strike).max(0.0)
        } else {
            (strike - value).max(0.0)
        };
        payoff * density
    };
    //A call is in the money below the kink, a put above it.
    let (a, b) = if is_call {
        (-12.0, kink_z.min(12.0))
    } else {
        (kink_z.max(-12.0), 12.0)
    };
    if a >= b {
        return 0.0;
    }
    let coarse = simpson(&integrand, a, b, 800);
    let fine = simpson(&integrand, a, b, 1600);
    assert!(
        (coarse - fine).abs() < 1e-8f64.max(fine.abs() * 1e-9),
        "reference quadrature not converged: {coarse} vs {fine}"
    );
    hull_white.bond_price_t(r_t, t, option_maturity).unwrap() * fine
}

const FIXTURES: [(f64, f64, f64, f64); 4] = [
    //curr_rate, a, b, sigma
    (0.05, 0.05, 0.05, 0.01),
    (STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG),
    (0.05, 0.1, 0.08, 0.08),
    (0.03, 0.4, 0.045, 0.005),
];

#[test]
fn the_forward_measure_drift_reproduces_the_bond_price() {
    //Guards the reference itself: repricing the coupon bond through the expiry-date
    //distribution of the rate must return today's bond price exactly.  Without the numeraire
    //drift this is out by ~1e-4, which is exactly the size of error that was showing up in
    //every option comparison below until the measure was fixed.
    let times = [2.5, 3.0, 3.5, 4.0];
    for &(curr, a, b, sigma) in FIXTURES.iter() {
        let (yield_curve, forward_curve) = hw_curves(curr, a, b, sigma);
        let hull_white = HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
        let (r_t, t, u) = (0.04, 1.0, 2.0);
        let bond = hull_white
            .coupon_bond_price_t(r_t, t, &times, 0.05)
            .unwrap();
        let (mean, sd) = expiry_rate_moments(&hull_white, r_t, t, u);
        let expected_value = simpson(
            &|z: f64| {
                let density = (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt();
                hull_white
                    .coupon_bond_price_t(mean + sd * z, u, &times, 0.05)
                    .unwrap()
                    * density
            },
            -12.0,
            12.0,
            4000,
        );
        let repriced = hull_white.bond_price_t(r_t, t, u).unwrap() * expected_value;
        assert!(
            (repriced - bond).abs() < 1e-12,
            "fixture {curr},{a},{b},{sigma}: repriced {repriced} vs bond {bond}"
        );
    }
}

#[test]
fn jamshidian_matches_the_direct_payoff_integral() {
    //The headline check: the decomposition, its bracket and its solver, against a quadrature
    //over the payoff that shares none of them.  The worst deviation across this whole grid is
    //about 2e-12.
    let times = [2.5, 3.0, 3.5, 4.0];
    for &(curr, a, b, sigma) in FIXTURES.iter() {
        let (yield_curve, forward_curve) = hw_curves(curr, a, b, sigma);
        let hull_white = HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
        for strike in [0.5f64, 0.8, 0.95, 1.0, 1.05, 1.3, 3.0] {
            for is_call in [true, false] {
                let priced = if is_call {
                    hull_white
                        .coupon_bond_call_t(0.04, 1.0, 2.0, &times, 0.05, strike)
                        .unwrap()
                } else {
                    hull_white
                        .coupon_bond_put_t(0.04, 1.0, 2.0, &times, 0.05, strike)
                        .unwrap()
                };
                let reference =
                    direct_payoff_price(&hull_white, 0.04, 1.0, 2.0, &times, 0.05, strike, is_call);
                assert!(
                    (priced - reference).abs() <= 1e-10f64.max(reference.abs() * 1e-9),
                    "fixture {curr},{a},{b},{sigma} {side} strike {strike}: {priced} vs {reference}",
                    side = if is_call { "call" } else { "put" }
                );
            }
        }
    }
}

#[test]
fn a_fine_strike_ladder_around_the_money_matches_the_integral() {
    //At the money the price is most sensitive to the critical rate, so a fine ladder through
    //the ATM region is where a sloppy root shows up first.
    let times = [2.5, 3.0, 3.5, 4.0];
    for &(curr, a, b, sigma) in FIXTURES.iter() {
        let (yield_curve, forward_curve) = hw_curves(curr, a, b, sigma);
        let hull_white = HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
        let (r_t, t, u) = (0.04, 1.0, 2.0);
        let underlying = hull_white
            .coupon_bond_price_t(r_t, t, &times, 0.05)
            .unwrap();
        for step in 0..21 {
            let strike = underlying * (0.9 + step as f64 * 0.01);
            for is_call in [true, false] {
                let priced = if is_call {
                    hull_white
                        .coupon_bond_call_t(r_t, t, u, &times, 0.05, strike)
                        .unwrap()
                } else {
                    hull_white
                        .coupon_bond_put_t(r_t, t, u, &times, 0.05, strike)
                        .unwrap()
                };
                let reference =
                    direct_payoff_price(&hull_white, r_t, t, u, &times, 0.05, strike, is_call);
                assert!(
                    (priced - reference).abs() <= 1e-10f64.max(reference.abs() * 1e-9),
                    "fixture {curr},{a},{b},{sigma} {side} strike {strike}: {priced} vs {reference}",
                    side = if is_call { "call" } else { "put" }
                );
            }
        }
    }
}

#[test]
fn a_negative_coupon_schedule_still_prices() {
    //A coupon below zero is not nonsense -- deep-discount instruments and some cross-currency
    //legs pay one -- and nothing in the decomposition needs the coupon to be positive as
    //long as the weights keep their sign.  The reference checks it independently.
    let times = [2.5, 3.0, 3.5, 4.0];
    let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
    let hull_white = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve).unwrap();
    let (r_t, t, u) = (0.04, 1.0, 2.0);
    for coupon_rate in [-0.02f64, -0.005, -0.0001] {
        for strike in [0.5f64, 0.9, 0.95, 1.0, 1.2] {
            for is_call in [true, false] {
                let priced = if is_call {
                    hull_white
                        .coupon_bond_call_t(r_t, t, u, &times, coupon_rate, strike)
                        .unwrap()
                } else {
                    hull_white
                        .coupon_bond_put_t(r_t, t, u, &times, coupon_rate, strike)
                        .unwrap()
                };
                let reference = direct_payoff_price(
                    &hull_white,
                    r_t,
                    t,
                    u,
                    &times,
                    coupon_rate,
                    strike,
                    is_call,
                );
                assert!(
                    (priced - reference).abs() <= 1e-10f64.max(reference.abs() * 1e-9),
                    "coupon {coupon_rate} {side} strike {strike}: {priced} vs {reference}",
                    side = if is_call { "call" } else { "put" }
                );
            }
        }
    }
}

mod solver_behaviour;
mod strikes;
