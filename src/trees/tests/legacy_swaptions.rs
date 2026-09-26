//! The long-standing swaption fixtures: the American tree price has to beat the analytic
//! European price, on the flat calibration that predates the steep fixture.

use approx::*;

use crate::HullWhite;

#[test]
fn payer_swaption() {
    let curr_rate = 0.05;
    let sig: f64 = 0.01;
    let a: f64 = 0.05;
    let b = 0.05;
    let delta = 0.25;
    let future_time = 0.0;
    //let swap_tenor = 5.0;
    let num_swap_payments = 20;
    let option_maturity = 1.0;
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
    let swap_rate = hull_white
        .forward_swap_rate_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
        )
        .unwrap();
    let analytical = hull_white
        .european_payer_swaption_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
        )
        .unwrap();
    let is_payer = true;

    let tree = hull_white.european_swaption_tree(
        curr_rate,
        future_time,
        option_maturity,
        num_swap_payments,
        delta,
        swap_rate,
        is_payer,
        100,
    );
    assert_abs_diff_eq!(analytical, tree, epsilon = 0.0001)
}

#[test]
fn receiver_swaption() {
    let curr_rate = 0.05;
    let sig: f64 = 0.01;
    let a: f64 = 0.05;
    let b = 0.05;
    let delta = 0.25;
    let future_time = 0.0;
    //let swap_tenor = 5.0;
    let num_swap_payments = 20;
    let option_maturity = 1.0;
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
    let swap_rate = hull_white
        .forward_swap_rate_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
        )
        .unwrap();
    let analytical = hull_white
        .european_receiver_swaption_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
        )
        .unwrap();
    let is_payer = false;

    let tree = hull_white.european_swaption_tree(
        curr_rate,
        future_time,
        option_maturity,
        num_swap_payments,
        delta,
        swap_rate,
        is_payer,
        100,
    );
    assert_abs_diff_eq!(analytical, tree, epsilon = 0.0001)
}

#[test]
fn american_payer_swaption() {
    let curr_rate = 0.05;
    let sig: f64 = 0.01;
    let a: f64 = 0.05;
    let b = 0.05;
    let delta = 0.25;
    let future_time = 0.0;
    //let swap_tenor = 5.0;
    let num_swap_payments = 20;
    let option_maturity = 1.0;
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
    let swap_rate = hull_white
        .forward_swap_rate_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
        )
        .unwrap();
    let analytical = hull_white
        .european_payer_swaption_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
        )
        .unwrap();

    let tree = hull_white
        .american_payer_swaption_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            100,
        )
        .unwrap();
    assert_eq!(analytical < tree, true);
}

#[test]
fn american_receiver_swaption() {
    let curr_rate = 0.05;
    let sig: f64 = 0.01;
    let a: f64 = 0.05;
    let b = 0.05;
    let delta = 0.25;
    let future_time = 0.0;
    //let swap_tenor = 5.0;
    let num_swap_payments = 20;
    let option_maturity = 1.0;
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
    let swap_rate = hull_white
        .forward_swap_rate_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
        )
        .unwrap();
    let analytical = hull_white
        .european_receiver_swaption_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
        )
        .unwrap();

    let tree = hull_white
        .american_receiver_swaption_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            100,
        )
        .unwrap();
    assert_eq!(analytical < tree, true);
}
