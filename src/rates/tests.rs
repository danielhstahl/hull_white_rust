//! Unit tests for the simple-rate instruments: `now`/`t` agreement for caplets and forward
//! Libor, Monte-Carlo checks of a caplet and of the Eurodollar-futures convexity against
//! the closed forms, and the cap/floor identities: caplet-floorlet parity, cap and floor summing
//! their periods, and the Monte-Carlo check that the floorlet really is `delta * max(K - L, 0)`.

use approx::*;
use rand::distributions::{Distribution, StandardNormal};

use crate::HullWhite;
use crate::error::HullWhiteError;
use crate::test_support::{BASELINE, STEEP_CURVE, get_rng_seed};

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

#[test]
fn test_caplet() {
    let fixture = BASELINE;
    let curr_rate = fixture.curr_rate;
    let future_time = 0.0;
    let option_maturity = 1.5;
    let strike = 0.02;
    let (yield_curve, forward_curve) = fixture.curves();
    let delta = fixture.delta;
    let seed: [u8; 32] = [2; 32];
    let mut rng_seed = get_rng_seed(seed);
    let normal = StandardNormal;
    let num_sims: usize = 1000; //hopefully accurate
    let num_discrete_steps: usize = 1000;
    let hull_white =
        HullWhite::init(fixture.a, fixture.sigma, &yield_curve, &forward_curve).unwrap();
    let total_sum = (0..num_sims).fold(0.0, |accum, _sample_index| {
        let mut sum_r = 0.0;
        let mut running_r = curr_rate;
        let dt = (option_maturity - future_time) / (num_discrete_steps as f64 - 1.0);
        (0..num_discrete_steps).for_each(|t_index| {
            let norm = normal.sample(&mut rng_seed);
            let curr_t = dt * (t_index as f64) + future_time;
            let curr_vol = hull_white.variance_r(curr_t, curr_t + dt).unwrap().sqrt();
            let curr_mu = hull_white.mu_r(running_r, curr_t, curr_t + dt).unwrap();
            running_r = curr_mu + curr_vol * norm;
            sum_r = sum_r + running_r * dt;
        });
        let libor_at_option_maturity = hull_white
            .libor_rate_t(running_r, option_maturity, delta)
            .unwrap();
        //And some more steps since discounted in arrears
        let more_steps = (delta / dt).floor() as usize;
        let new_dt = delta / (more_steps as f64 - 1.0);
        (1..more_steps).for_each(|t_index| {
            let norm = normal.sample(&mut rng_seed);
            let curr_t = new_dt * (t_index as f64) + option_maturity;
            let curr_vol = hull_white
                .variance_r(curr_t, curr_t + new_dt)
                .unwrap()
                .sqrt();
            let curr_mu = hull_white.mu_r(running_r, curr_t, curr_t + new_dt).unwrap();
            running_r = curr_mu + curr_vol * norm;
            sum_r = sum_r + running_r * new_dt;
        });

        if libor_at_option_maturity > strike {
            accum + (libor_at_option_maturity - strike) * ((-sum_r).exp()) //discount
        } else {
            accum
        }
    });
    let average_caplet = delta * (total_sum / (num_sims as f64));
    let analytical_caplet = hull_white
        .caplet_now(option_maturity, delta, strike)
        .unwrap();
    assert_abs_diff_eq!(average_caplet, analytical_caplet, epsilon = 0.0001);
}

#[test]
fn test_edf() {
    let fixture = BASELINE;
    let curr_rate = fixture.curr_rate;
    let future_time = 0.0;
    let option_maturity = 1.5;
    let (yield_curve, forward_curve) = fixture.curves();
    let delta = fixture.delta;
    let seed: [u8; 32] = [2; 32];
    let mut rng_seed = get_rng_seed(seed);
    let normal = StandardNormal;
    let num_sims: usize = 1000000; //hopefully accurate
    let hull_white =
        HullWhite::init(fixture.a, fixture.sigma, &yield_curve, &forward_curve).unwrap();
    let mu = hull_white
        .mu_r(curr_rate, future_time, option_maturity)
        .unwrap();
    let vol = hull_white
        .variance_r(future_time, option_maturity)
        .unwrap()
        .sqrt();
    let total_sum = (0..num_sims).fold(0.0, |accum, _sample_index| {
        let norm = normal.sample(&mut rng_seed);
        let final_r = mu + vol * norm;
        let final_bond = hull_white
            .bond_price_t(final_r, option_maturity, option_maturity + delta)
            .unwrap();
        accum + 1.0 / final_bond
    });
    let average_edf = ((total_sum / (num_sims as f64)) - 1.0) / delta;

    let analytical_edf = hull_white
        .euro_dollar_future_t(curr_rate, future_time, option_maturity, delta)
        .unwrap();
    assert_abs_diff_eq!(average_edf, analytical_edf, epsilon = 0.0001);
}

/// Monte-Carlo check of the floorlet against the closed form, mirroring `test_caplet`.  The
/// simulated payoff is `delta * max(strike - Libor, 0)` discounted along the sampled short-rate
/// path, which is the orientation the bond-call twin in `floorlet_now` is supposed to reproduce.
///
/// The strike sits well above the forward here, so the floor carries the value and the caplet on the
/// same period is worth about a twentieth as much: a pricer that had the cap and floor payoff
/// orientations swapped could not produce that asymmetry, and would also fail the parity check in
/// `cap_floor_parity_is_the_forward_leg` in the other direction.
#[test]
fn test_floorlet() {
    let fixture = BASELINE;
    let curr_rate = fixture.curr_rate;
    let future_time = 0.0;
    let option_maturity = 1.5;
    let strike = 0.05;
    let (yield_curve, forward_curve) = fixture.curves();
    let delta = 0.25;
    let seed: [u8; 32] = [2; 32];
    let mut rng_seed = get_rng_seed(seed);
    let normal = StandardNormal;
    let num_sims: usize = 1000; //hopefully accurate
    let num_discrete_steps: usize = 1000;
    let hull_white =
        HullWhite::init(fixture.a, fixture.sigma, &yield_curve, &forward_curve).unwrap();
    let total_sum = (0..num_sims).fold(0.0, |accum, _sample_index| {
        let mut sum_r = 0.0;
        let mut running_r = curr_rate;
        let dt = (option_maturity - future_time) / (num_discrete_steps as f64 - 1.0);
        (0..num_discrete_steps).for_each(|t_index| {
            let norm = normal.sample(&mut rng_seed);
            let curr_t = dt * (t_index as f64) + future_time;
            let curr_vol = hull_white.variance_r(curr_t, curr_t + dt).unwrap().sqrt();
            let curr_mu = hull_white.mu_r(running_r, curr_t, curr_t + dt).unwrap();
            running_r = curr_mu + curr_vol * norm;
            sum_r += running_r * dt;
        });
        let libor_at_option_maturity = hull_white
            .libor_rate_t(running_r, option_maturity, delta)
            .unwrap();
        //And some more steps since discounted in arrears
        let more_steps = (delta / dt).floor() as usize;
        let new_dt = delta / (more_steps as f64 - 1.0);
        (1..more_steps).for_each(|t_index| {
            let norm = normal.sample(&mut rng_seed);
            let curr_t = new_dt * (t_index as f64) + option_maturity;
            let curr_vol = hull_white
                .variance_r(curr_t, curr_t + new_dt)
                .unwrap()
                .sqrt();
            let curr_mu = hull_white.mu_r(running_r, curr_t, curr_t + new_dt).unwrap();
            running_r = curr_mu + curr_vol * norm;
            sum_r += running_r * new_dt;
        });

        if libor_at_option_maturity < strike {
            accum + (strike - libor_at_option_maturity) * ((-sum_r).exp()) //discount
        } else {
            accum
        }
    });
    let average_floorlet = delta * (total_sum / (num_sims as f64));
    let analytical_floorlet = hull_white
        .floorlet_now(option_maturity, delta, strike)
        .unwrap();
    let analytical_caplet = hull_white
        .caplet_now(option_maturity, delta, strike)
        .unwrap();
    assert_abs_diff_eq!(average_floorlet, analytical_floorlet, epsilon = 0.0001);
    assert!(
        analytical_floorlet > 10.0 * analytical_caplet,
        "floor {analytical_floorlet} vs cap {analytical_caplet} at a strike far above the forward"
    );
}
