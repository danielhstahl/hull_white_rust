//! Unit tests for the simple-rate instruments: `now`/`t` agreement for caplets and forward
//! Libor, the cap/floor identities (caplet-floorlet parity, a cap or floor summing its periods,
//! the forward-strike coincidence, schedule validation), and the Monte-Carlo checks of a caplet,
//! a floorlet and the Eurodollar-futures convexity against the closed forms.
//!
//! The Monte-Carlo checks go through [`crate::mc`]: a fixed seed per test, antithetic variates,
//! one integer-step time-grid convention, and a tolerance computed as
//! `k * standard_error + estimated discretisation residual` rather than a literal.  Sample size
//! and grid resolution come from [`crate::mc::scale`], tiered on the `slow` feature: a plain
//! `cargo test` runs the light tier (~0.25s for the three checks, bands of order `1e-4`) and
//! `cargo test --features slow` runs the heavy one (~18s, bands ~4x tighter); CI runs both.
//! The conventions, the measured variance reductions and the bias-study numbers are written up in
//! [`crate::mc`].

use approx::*;

use crate::HullWhite;
use crate::error::HullWhiteError;
use crate::mc::{self, Estimate, Grid};
use crate::test_support::{BASELINE, STEEP_CURVE};

/// Seed slot for the caplet / floorlet path estimator.  One slot per test so the tests do not
/// share a stretch of random numbers and cannot covary — see `mc::rng_for`.
const SIMPLE_RATE_OPTION_SEED: u8 = 11;

/// Seed slot for the Eurodollar-future estimator.
const EURODOLLAR_SEED: u8 = 12;

/// The two routes to a floorlet have to be the same number: the bond-call twin of the caplet's bond
/// put, and the `t` form run at `(r(0), 0)`.
#[test]
fn floorlet_now_is_the_t_form_at_zero() {
    let (yield_curve, forward_curve) = STEEP_CURVE.curves();
    let hull_white = HullWhite::init(
        STEEP_CURVE.a,
        STEEP_CURVE.sigma,
        &yield_curve,
        &forward_curve,
    )
    .unwrap();
    let r0 = hull_white.short_rate_now().unwrap();
    let delta = 0.25;
    for option_maturity in [0.25, 1.0, 1.5, 5.0] {
        for strike in [0.0, 0.02, 0.05, 0.10] {
            let now = hull_white
                .floorlet_now(option_maturity, delta, strike)
                .unwrap();
            let via_t = hull_white
                .floorlet_t(r0, 0.0, option_maturity, delta, strike)
                .unwrap();
            assert_abs_diff_eq!(now, via_t, epsilon = 1e-15);
            assert!(now > 0.0, "floorlet({option_maturity}, {strike}) = {now}");
        }
    }
}

/// Cap-floor parity, per period and on both time shapes.
///
/// A caplet and a floorlet on the same period are the call and the put of one bond option scaled by
/// `1 + delta * K`, so their difference is not an approximation shared between two independent
/// pricers -- it is the forward leg itself:
///
/// ```text
/// caplet - floorlet = P(t, T) - (1 + delta * K) * P(t, T + delta)
/// ```
///
/// which is `delta * (forward Libor - K)` discounted, and vanishes when the strike is the forward.
#[test]
fn cap_floor_parity_is_the_forward_leg() {
    let (yield_curve, forward_curve) = STEEP_CURVE.curves();
    let hull_white = HullWhite::init(
        STEEP_CURVE.a,
        STEEP_CURVE.sigma,
        &yield_curve,
        &forward_curve,
    )
    .unwrap();
    let delta = 0.25;
    for (t, option_maturity) in [(0.0, 1.0), (0.5, 2.0), (1.0, 1.25), (2.0, 7.0)] {
        for strike in [-0.02, 0.0, 0.02, 0.05, 0.20] {
            let caplet = hull_white
                .caplet_t(STEEP_CURVE.curr_rate, t, option_maturity, delta, strike)
                .unwrap();
            let floorlet = hull_white
                .floorlet_t(STEEP_CURVE.curr_rate, t, option_maturity, delta, strike)
                .unwrap();
            let near = hull_white
                .bond_price_t(STEEP_CURVE.curr_rate, t, option_maturity)
                .unwrap();
            let far = hull_white
                .bond_price_t(STEEP_CURVE.curr_rate, t, option_maturity + delta)
                .unwrap();
            let leg = near - (1.0 + delta * strike) * far;
            assert_abs_diff_eq!(caplet - floorlet, leg, epsilon = 1e-13);
            //Both option sides are non-negative, and neither is ever exactly zero here.
            assert!(caplet >= 0.0 && floorlet >= 0.0);
            if t == 0.0 {
                //The `now` twins are the same instrument at `r_t = r(0)` = the fixture's start rate.
                assert_abs_diff_eq!(
                    hull_white
                        .caplet_now(option_maturity, delta, strike)
                        .unwrap(),
                    caplet,
                    epsilon = 1e-15
                );
                assert_abs_diff_eq!(
                    hull_white
                        .floorlet_now(option_maturity, delta, strike)
                        .unwrap(),
                    floorlet,
                    epsilon = 1e-15
                );
            }
        }
    }
}

/// At the forward strike a cap and a floor on the same schedule differ only by the (zero) forward
/// leg, so the swap of one for the other is a swap at par -- the check that the strike convention on
/// the two sides is the same one.
#[test]
fn floor_and_cap_agree_at_the_forward_strike() {
    let (yield_curve, forward_curve) = STEEP_CURVE.curves();
    let hull_white = HullWhite::init(
        STEEP_CURVE.a,
        STEEP_CURVE.sigma,
        &yield_curve,
        &forward_curve,
    )
    .unwrap();
    let delta = 0.25;
    for option_maturity in [1.0, 2.5, 4.0] {
        let forward = hull_white
            .forward_libor_rate_now(option_maturity, delta)
            .unwrap();
        let periods = [(option_maturity, forward)];
        let cap = hull_white.cap_now(&periods, delta).unwrap();
        let floor = hull_white.floor_now(&periods, delta).unwrap();
        assert_abs_diff_eq!(cap, floor, epsilon = 1e-12);
        //Away from the forward the two separate in the expected direction: a cap is dear when the
        //forward is above the strike, a floor when it is below.
        let below = forward - 0.02;
        assert!(
            hull_white
                .cap_now(&[(option_maturity, below)], delta)
                .unwrap()
                > hull_white
                    .floor_now(&[(option_maturity, below)], delta)
                    .unwrap(),
            "cap should be dearer than the floor with the forward above the strike"
        );
        let above = forward + 0.02;
        assert!(
            hull_white
                .floor_now(&[(option_maturity, above)], delta)
                .unwrap()
                > hull_white
                    .cap_now(&[(option_maturity, above)], delta)
                    .unwrap(),
            "floor should be dearer than the cap with the forward below the strike"
        );
    }
}

/// A cap (floor) is the sum of its caplets (floorlets), and a 20 period instrument prices in one
/// call rather than in a caller-managed loop.
#[test]
fn cap_and_floor_are_the_sum_of_their_periods() {
    let (yield_curve, forward_curve) = STEEP_CURVE.curves();
    let hull_white = HullWhite::init(
        STEEP_CURVE.a,
        STEEP_CURVE.sigma,
        &yield_curve,
        &forward_curve,
    )
    .unwrap();
    let delta = 0.25;
    let strikes = [0.02, 0.03, 0.04, 0.05, 0.06];
    let periods: Vec<(f64, f64)> = (1..=20)
        .map(|index| (index as f64 * delta, strikes[(index - 1) % strikes.len()]))
        .collect();

    let cap = hull_white.cap_now(&periods, delta).unwrap();
    let floor = hull_white.floor_now(&periods, delta).unwrap();
    let caplets: f64 = periods
        .iter()
        .map(|&(option_maturity, strike)| {
            hull_white
                .caplet_now(option_maturity, delta, strike)
                .unwrap()
        })
        .sum();
    let floorlets: f64 = periods
        .iter()
        .map(|&(option_maturity, strike)| {
            hull_white
                .floorlet_now(option_maturity, delta, strike)
                .unwrap()
        })
        .sum();
    assert!(cap > 0.0 && floor > 0.0, "cap {cap}, floor {floor}");
    assert_abs_diff_eq!(cap, caplets, epsilon = 1e-12);
    assert_abs_diff_eq!(floor, floorlets, epsilon = 1e-12);

    //Schedule-level cap-floor parity: the gap is the sum of the forward legs.
    let legs: f64 = periods
        .iter()
        .map(|&(option_maturity, strike)| {
            let near = hull_white.bond_price_now(option_maturity).unwrap();
            let far = hull_white.bond_price_now(option_maturity + delta).unwrap();
            near - (1.0 + delta * strike) * far
        })
        .sum();
    assert_abs_diff_eq!(cap - floor, legs, epsilon = 1e-11);

    //Adding a period adds exactly that caplet: aggregation never re-weights what was already there.
    let mut longer = periods.clone();
    longer.push((21.0 * delta, 0.045));
    assert_abs_diff_eq!(
        hull_white.cap_now(&longer, delta).unwrap() - cap,
        hull_white.caplet_now(21.0 * delta, delta, 0.045).unwrap(),
        epsilon = 1e-12
    );
    assert_abs_diff_eq!(
        hull_white.floor_now(&longer, delta).unwrap() - floor,
        hull_white.floorlet_now(21.0 * delta, delta, 0.045).unwrap(),
        epsilon = 1e-12
    );

    //The `t` shape does the same thing from a state.
    let future_periods: Vec<(f64, f64)> = (1..=20)
        .map(|index| {
            (
                1.0 + index as f64 * delta,
                strikes[(index - 1) % strikes.len()],
            )
        })
        .collect();
    let cap_t = hull_white.cap_t(0.04, 1.0, &future_periods, delta).unwrap();
    let caplets_t: f64 = future_periods
        .iter()
        .map(|&(option_maturity, strike)| {
            hull_white
                .caplet_t(0.04, 1.0, option_maturity, delta, strike)
                .unwrap()
        })
        .sum();
    assert_abs_diff_eq!(cap_t, caplets_t, epsilon = 1e-12);
    let floor_t = hull_white
        .floor_t(0.04, 1.0, &future_periods, delta)
        .unwrap();
    let floorlets_t: f64 = future_periods
        .iter()
        .map(|&(option_maturity, strike)| {
            hull_white
                .floorlet_t(0.04, 1.0, option_maturity, delta, strike)
                .unwrap()
        })
        .sum();
    assert_abs_diff_eq!(floor_t, floorlets_t, epsilon = 1e-12);
}

/// An empty schedule is not a free cap, and one bad period is not a dropped leg.
#[test]
fn cap_and_floor_schedule_errors() {
    let (yield_curve, forward_curve) = STEEP_CURVE.curves();
    let hull_white = HullWhite::init(
        STEEP_CURVE.a,
        STEEP_CURVE.sigma,
        &yield_curve,
        &forward_curve,
    )
    .unwrap();
    let delta = STEEP_CURVE.delta;
    for empty in [
        hull_white.cap_now(&[], delta),
        hull_white.floor_now(&[], delta),
        hull_white.cap_t(0.04, 1.0, &[], delta),
        hull_white.floor_t(0.04, 1.0, &[], delta),
    ] {
        assert!(
            matches!(empty, Err(HullWhiteError::InvalidInput(_))),
            "an empty schedule should be InvalidInput, got {empty:?}"
        );
    }
    //A period whose expiry is not strictly after the valuation date is refused, not silently skipped:
    //dropping it would return a price for a shorter cap than was asked for.
    let bad = hull_white
        .cap_t(0.04, 1.0, &[(1.25, 0.04), (0.5, 0.04)], delta)
        .unwrap_err();
    assert!(matches!(bad, HullWhiteError::InvalidInput(_)), "{bad:?}");
    //A strike at or beyond the `1 + delta * K = 0` singularity is refused on every period.
    let singular = -1.0 / delta;
    assert!(
        hull_white
            .floor_now(&[(1.0, 0.04), (2.0, singular)], delta)
            .is_err()
    );
    assert!(hull_white.cap_now(&[(1.0, 0.04)], 0.0).is_err());
    //Order carries no meaning here: the same periods shuffled price identically.
    let one = [(1.0, 0.05), (2.0, 0.04), (3.0, 0.03)];
    let other = [(3.0, 0.03), (1.0, 0.05), (2.0, 0.04)];
    assert_abs_diff_eq!(
        hull_white.cap_now(&one, delta).unwrap(),
        hull_white.cap_now(&other, delta).unwrap(),
        epsilon = 1e-15
    );
    assert_abs_diff_eq!(
        hull_white.floor_now(&one, delta).unwrap(),
        hull_white.floor_now(&other, delta).unwrap(),
        epsilon = 1e-15
    );
}

#[test]
fn compare_caplet() {
    let fixture = BASELINE;
    let curr_rate = fixture.curr_rate;
    let future_time = 0.0;
    let option_maturity = 1.5;
    let strike = 0.02;
    let (yield_curve, forward_curve) = fixture.curves();
    let delta = fixture.delta;
    let hull_white =
        HullWhite::init(fixture.a, fixture.sigma, &yield_curve, &forward_curve).unwrap();
    let caplet_n = hull_white
        .caplet_now(option_maturity, delta, strike)
        .unwrap();
    let caplet = hull_white
        .caplet_t(curr_rate, future_time, option_maturity, delta, strike)
        .unwrap();
    assert_abs_diff_eq!(caplet_n, caplet, epsilon = 0.00001);
}

#[test]
fn compare_libor() {
    let fixture = BASELINE;
    let curr_rate = fixture.curr_rate;
    let future_time = 0.0;
    let maturity = 1.5;
    let (yield_curve, forward_curve) = fixture.curves();
    let delta = fixture.delta;
    let hull_white =
        HullWhite::init(fixture.a, fixture.sigma, &yield_curve, &forward_curve).unwrap();
    let libor_n = hull_white.forward_libor_rate_now(maturity, delta).unwrap();
    let libor_t = hull_white
        .forward_libor_rate_t(curr_rate, future_time, maturity, delta)
        .unwrap();
    assert_abs_diff_eq!(libor_n, libor_t, epsilon = 0.0001);
}

// ---------------------------------------------------------------------------------------
// Monte-Carlo checks against the closed forms.
//
// These use `crate::mc`: fixed per-test seeds, antithetic variates, a uniform integer-step time
// grid, and a tolerance computed as `k * standard_error + estimated discretisation residual`
// instead of a literal.  Size comes from `mc::scale()`, which is tiered on the `slow` feature:
// plain `cargo test` runs a light tier (~0.25s for the three checks, bands of order 1e-4) and
// `cargo test --features slow` runs the heavy one (~18s, bands ~4x tighter at a few times 1e-5);
// CI runs both.  See `src/mc.rs` for the conventions and the measured numbers.
// ---------------------------------------------------------------------------------------

/// One period of a cap or a floor, described the way the Monte-Carlo simulates it rather than the
/// way the closed form is written.
///
/// The simulated claim is the *arrears* one: `delta * (L_T - K)^+` (or `(K - L_T)^+`) paid at
/// `T + delta`, discounted by the sampled `exp(-integral r ds)` over `[0, T + delta]`.  That is
/// a genuinely independent route to the closed form's number, which arrives through the
/// deposit-bond identity instead:
///
/// ```text
///   E[ exp(-int_0^{T+delta} r) * delta * (L_T - K)^+ ]
/// = E[ exp(-int_0^T r) * P(T, T+delta) * delta * (L_T - K)^+ ]   tower, P(T,T+d) = E[e^-int | F_T]
/// = E[ exp(-int_0^T r) * (1 - (1 + delta*K) * P(T, T+delta))^+ ]  since delta*L_T = 1/B_T - 1
/// = (1 + delta*K) * put_on_B(P(., T+delta); K' = 1/(1 + delta*K)) = caplet_now
/// ```
///
/// with the call standing in for the put, and the sign flipped, for the floorlet.  Putting the
/// payoff direction in a field means the caplet and the floorlet walk exactly the same path and
/// cannot drift apart in convention.
struct ArrearsOption {
    /// The short rate the simulated path starts from.
    start_rate: f64,
    /// `T`: the fixing date, which is also the option's expiry.
    option_maturity: f64,
    /// The accrual period: the claim fixes over `[T, T + delta]` and pays at `T + delta`.
    delta: f64,
    strike: f64,
    /// `true` for a caplet (`(L - K)^+`), `false` for a floorlet (`(K - L)^+`).
    cap: bool,
}

impl ArrearsOption {
    /// The two legs of this instrument's life, from the valuation date to the fixing and from the
    /// fixing to the payment date, each stepped uniformly at `steps_per_year`
    /// ([`Grid::for_legs`], which never resolves a coarser `dt` than asked for).
    fn legs(&self, valuation_time: f64, steps_per_year: usize) -> Grid {
        Grid::for_legs(
            valuation_time,
            &[self.option_maturity, self.option_maturity + self.delta],
            steps_per_year,
        )
    }

    /// This path's discounted payoff: the fixing read off the rate *at* the expiry date, the
    /// period's simple-rate payoff, discounted by this path's own short-rate integral to the
    /// payment date.
    fn payoff(
        &self,
        hull_white: &HullWhite<impl Fn(f64) -> f64 + Sync, impl Fn(f64) -> f64 + Sync>,
        state: &mc::PathState,
    ) -> f64 {
        //The rate at the end of the first leg is the rate at the fixing date — not one step past
        //it, which is what the old `dt = span / (n - 1)` convention did.
        let fixing_rate = state.rate_at_end_of_leg(0);
        let fixing = hull_white
            .libor_rate_t(fixing_rate, self.option_maturity, self.delta)
            .unwrap();
        let payoff = if self.cap {
            self.delta * (fixing - self.strike).max(0.0)
        } else {
            self.delta * (self.strike - fixing).max(0.0)
        };
        payoff * (-state.integral).exp()
    }

    /// The closed form this estimator is checked against.
    fn analytic(
        &self,
        hull_white: &HullWhite<impl Fn(f64) -> f64 + Sync, impl Fn(f64) -> f64 + Sync>,
    ) -> f64 {
        if self.cap {
            hull_white
                .caplet_now(self.option_maturity, self.delta, self.strike)
                .unwrap()
        } else {
            hull_white
                .floorlet_now(self.option_maturity, self.delta, self.strike)
                .unwrap()
        }
    }
}

/// Simulate the arrears option on `grid` with antithetic pairs.
fn simulate(
    hull_white: &HullWhite<impl Fn(f64) -> f64 + Sync, impl Fn(f64) -> f64 + Sync>,
    option: &ArrearsOption,
    grid: &Grid,
    pairs: usize,
) -> Estimate {
    mc::run_paths(
        hull_white,
        mc::rng_for(SIMPLE_RATE_OPTION_SEED),
        grid,
        option.start_rate,
        pairs,
        |state: &mc::PathState| option.payoff(hull_white, state),
    )
}

/// The whole caplet/floorlet check: coarse grid, fine grid at exactly half the step size, the
/// Richardson residual from those two, and the `k * SE + residual` assertion — plus the report on
/// stdout, which is what `--nocapture` is for.
fn assert_arrears_option_within_budget(
    label: &str,
    hull_white: &HullWhite<impl Fn(f64) -> f64 + Sync, impl Fn(f64) -> f64 + Sync>,
    option: &ArrearsOption,
) -> Estimate {
    let scale = mc::scale();
    let coarse = option.legs(0.0, scale.steps_per_year / 2);
    let fine = coarse.refined(2);
    let coarse_estimate = simulate(hull_white, option, &coarse, scale.pairs);
    let fine_estimate = simulate(hull_white, option, &fine, scale.pairs);
    //Refined by an exact factor of 2 on every leg, so the Richardson factor is exactly 2 rather
    //than something that fell out of how each leg happened to round.
    let ratio = coarse.step_ratio(&fine);
    let bias = mc::discretisation_bias(&coarse_estimate, &fine_estimate, ratio);
    mc::Check::with_bias_study(label, fine_estimate, option.analytic(hull_white), bias)
        .assert_within_budget();
    //The variance reduction has to actually be there.  If it ever comes out at ~1 the pairing is
    //not doing its job, and the tight heavy-tier bands would be resting on nothing.
    assert!(
        fine_estimate.variance_reduction > 1.2,
        "{label}: antithetic pairing bought only {:.2}x",
        fine_estimate.variance_reduction
    );
    fine_estimate
}

/// Monte-Carlo caplet vs the closed form, with the error budget stated rather than asserted.
#[test]
fn monte_carlo_caplet_matches_the_closed_form() {
    let fixture = BASELINE;
    let (yield_curve, forward_curve) = fixture.curves();
    let hull_white =
        HullWhite::init(fixture.a, fixture.sigma, &yield_curve, &forward_curve).unwrap();
    assert_arrears_option_within_budget(
        "caplet(T=1.5, K=0.02) simulated arrears, antithetic, fine dt = 1/steps_per_year",
        &hull_white,
        &ArrearsOption {
            start_rate: fixture.curr_rate,
            option_maturity: 1.5,
            delta: fixture.delta,
            strike: 0.02,
            cap: true,
        },
    );
}

/// Monte-Carlo floorlet vs the closed form.  The strike sits well above the forward, so the floor
/// carries the value and the caplet on the same period is worth a fraction of it: a pricer with the
/// cap and floor orientations swapped could not reproduce that asymmetry, and would also fail
/// `cap_floor_parity_is_the_forward_leg` in the other direction.
#[test]
fn monte_carlo_floorlet_matches_the_closed_form() {
    let fixture = BASELINE;
    let (yield_curve, forward_curve) = fixture.curves();
    let hull_white =
        HullWhite::init(fixture.a, fixture.sigma, &yield_curve, &forward_curve).unwrap();
    let strike = 0.05;
    assert_arrears_option_within_budget(
        "floorlet(T=1.5, K=0.05) simulated arrears, antithetic, fine dt = 1/steps_per_year",
        &hull_white,
        &ArrearsOption {
            start_rate: fixture.curr_rate,
            option_maturity: 1.5,
            delta: fixture.delta,
            strike,
            cap: false,
        },
    );
    let caplet = hull_white.caplet_now(1.5, fixture.delta, strike).unwrap();
    let floorlet = hull_white.floorlet_now(1.5, fixture.delta, strike).unwrap();
    assert!(
        floorlet > 10.0 * caplet,
        "floor {floorlet} vs cap {caplet} at a strike far above the forward"
    );
}

/// Monte-Carlo check of the Eurodollar-future convexity against the closed form.
///
/// Unlike the caplet and floorlet checks this estimator has **no** discretisation at all: the
/// terminal rate is drawn from its exact conditional distribution `N(mu_r, var_r)` over
/// `[0, T]`, and the futures is a deterministic function of that one draw
/// (`(E[1 / B(T, T + delta)] - 1) / delta`), so there is no time step whose coarseness could
/// bias the answer.  The budget is therefore sampling error and nothing else — `k * SE`, with the
/// bias term declared as [`mc::NO_DISCRETISATION`] rather than quietly left out.
///
/// Because each path costs a single normal draw rather than hundreds of steps, the same wall-clock
/// budget buys about ten times as many pairs as the path estimators get, so this test runs at
/// `10 x scale().pairs` and lands on a tighter band than the caplet check for less time.
#[test]
fn monte_carlo_euro_dollar_future_matches_the_closed_form() {
    let fixture = BASELINE;
    let (yield_curve, forward_curve) = fixture.curves();
    let hull_white =
        HullWhite::init(fixture.a, fixture.sigma, &yield_curve, &forward_curve).unwrap();
    let option_maturity = 1.5;
    let delta = fixture.delta;
    let mu = hull_white
        .mu_r(fixture.curr_rate, 0.0, option_maturity)
        .unwrap();
    let vol = hull_white.variance_r(0.0, option_maturity).unwrap().sqrt();
    let pairs = mc::scale().pairs * 10;
    let estimate = mc::run(mc::rng_for(EURODOLLAR_SEED), 1, pairs, |normals| {
        let rate = mu + vol * normals[0];
        1.0 / hull_white
            .bond_price_t(rate, option_maturity, option_maturity + delta)
            .unwrap()
    })
    .rescale(1.0 / delta, -1.0 / delta);
    let analytic = hull_white
        .euro_dollar_future_t(fixture.curr_rate, 0.0, option_maturity, delta)
        .unwrap();
    mc::Check::undiscretised(
        "Eurodollar future(T=1.5, exact terminal draw, antithetic, no time stepping)",
        estimate,
        analytic,
    )
    .assert_within_budget();
    //The convexity itself: the future sits above the forward it quotes against, and only because of
    //volatility.  A zero-variance limit collapses the two together.
    let forward = hull_white
        .forward_libor_rate_now(option_maturity, delta)
        .unwrap();
    assert!(
        analytic > forward,
        "future {analytic} should exceed the forward {forward} by the convexity term"
    );
}
