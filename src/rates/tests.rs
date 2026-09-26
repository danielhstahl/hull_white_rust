//! Unit tests for the simple-rate instruments: `now`/`t` agreement for caplets and forward
//! Libor, and Monte-Carlo checks of a caplet and of the Eurodollar-futures convexity against
//! the closed forms.

use approx::*;
use rand::distributions::{Distribution, StandardNormal};

use crate::HullWhite;
use crate::testutil::get_rng_seed;

#[test]
fn compare_caplet() {
    let curr_rate = 0.02;
    let sig: f64 = 0.02;
    let a: f64 = 0.3;
    let b = 0.04;
    let future_time = 0.0;
    let option_maturity = 1.5;
    let strike = 0.02;
    let yield_curve = |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let delta = 0.25;
    let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
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
    let curr_rate = 0.02;
    let sig: f64 = 0.02;
    let a: f64 = 0.3;
    let b = 0.04;
    let future_time = 0.0;
    let maturity = 1.5;
    let yield_curve = |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let delta = 0.25;
    let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
    let libor_n = hull_white.forward_libor_rate_now(maturity, delta).unwrap();
    let libor_t = hull_white
        .forward_libor_rate_t(curr_rate, future_time, maturity, delta)
        .unwrap();
    assert_abs_diff_eq!(libor_n, libor_t, epsilon = 0.0001);
}

#[test]
fn test_caplet() {
    let curr_rate = 0.02;
    let sig: f64 = 0.02;
    let a: f64 = 0.3;
    let b = 0.04;
    let future_time = 0.0;
    let option_maturity = 1.5;
    let strike = 0.02;
    let yield_curve = |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let delta = 0.25;
    let seed: [u8; 32] = [2; 32];
    let mut rng_seed = get_rng_seed(seed);
    let normal = StandardNormal;
    let num_sims: usize = 1000; //hopefully accurate
    let num_discrete_steps: usize = 1000;
    let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
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
    let curr_rate = 0.02;
    let sig: f64 = 0.02;
    let a: f64 = 0.3;
    let b = 0.04;
    let future_time = 0.0;
    let option_maturity = 1.5;
    let yield_curve = |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let delta = 0.25;
    let seed: [u8; 32] = [2; 32];
    let mut rng_seed = get_rng_seed(seed);
    let normal = StandardNormal;
    let num_sims: usize = 1000000; //hopefully accurate
    let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
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
