//! The long-standing swaption fixtures: the American tree price has to beat the analytic
//! European price, on the flat calibration that predates the steep fixture.
//!
//! The European side of each pair is priced with the public
//! [`european_payer_swaption_tree`](crate::HullWhite::european_payer_swaption_tree) /
//! [`european_receiver_swaption_tree`](crate::HullWhite::european_receiver_swaption_tree)
//! helpers, so this file is the same "tree vs Jamshidian" cross-check from the crate's own tests
//! that a downstream user can now write.

use approx::*;

use crate::HullWhite;
use crate::test_support::FLAT_5PCT;

#[test]
fn payer_swaption() {
    let fixture = FLAT_5PCT;
    let curr_rate = fixture.curr_rate;
    let delta = fixture.delta;
    let future_time = 0.0;
    //let swap_tenor = 5.0;
    let num_swap_payments = 20;
    let option_maturity = 1.0;
    let curve = fixture.curve();
    let hull_white = HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();
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
        .european_payer_swaption_tree(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            100,
        )
        .unwrap();
    assert_abs_diff_eq!(analytical, tree, epsilon = 0.0001)
}

#[test]
fn receiver_swaption() {
    let fixture = FLAT_5PCT;
    let curr_rate = fixture.curr_rate;
    let delta = fixture.delta;
    let future_time = 0.0;
    //let swap_tenor = 5.0;
    let num_swap_payments = 20;
    let option_maturity = 1.0;
    let curve = fixture.curve();
    let hull_white = HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();
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
        .european_receiver_swaption_tree(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            100,
        )
        .unwrap();
    assert_abs_diff_eq!(analytical, tree, epsilon = 0.0001)
}

#[test]
fn american_payer_swaption() {
    let fixture = FLAT_5PCT;
    let curr_rate = fixture.curr_rate;
    let delta = fixture.delta;
    let future_time = 0.0;
    //let swap_tenor = 5.0;
    let num_swap_payments = 20;
    let option_maturity = 1.0;
    let curve = fixture.curve();
    let hull_white = HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();
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
    assert!(analytical < tree);
}

#[test]
fn american_receiver_swaption() {
    let fixture = FLAT_5PCT;
    let curr_rate = fixture.curr_rate;
    let delta = fixture.delta;
    let future_time = 0.0;
    //let swap_tenor = 5.0;
    let num_swap_payments = 20;
    let option_maturity = 1.0;
    let curve = fixture.curve();
    let hull_white = HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();
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
    assert!(analytical < tree);
}
