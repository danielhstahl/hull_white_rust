//! Tests for the input contract: an invalid instrument is an error that names the argument,
//! never a panic and never a plausible-looking number.
//!
//! Before these guards three shapes of failure were live:
//!   * `coupon_times.len() - 1` on an empty schedule -> usize underflow panic;
//!   * an expired swap / `delta <= 0` -> the float-to-usize casts saturated at 0, the payment
//!     loop vanished, and a bond leg came back out looking like a swap price;
//!   * NaN/inf travelled straight through `exp`/`/` into the returned price.
//!
//! The guards under `validation::` now turn each of those into `InvalidInput` naming the
//! argument, and a non-finite *computed* result into `NumericalError`.

use crate::error::HullWhiteError;
use crate::schedules::get_coupon_times;
use crate::testutil::yvf_setup;

fn expect_invalid<T: std::fmt::Debug>(result: Result<T, HullWhiteError>, needle: &str) {
    match result {
        Err(HullWhiteError::InvalidInput(msg)) => assert!(
            msg.contains(needle),
            "expected the error to name {needle:?}, got: {msg}"
        ),
        Err(other) => panic!("expected InvalidInput naming {needle:?}, got {other:?}"),
        Ok(value) => panic!("expected InvalidInput naming {needle:?}, got Ok({value:?})"),
    }
}

#[test]
fn empty_coupon_times_is_invalid_input() {
    yvf_setup!(hull_white);
    expect_invalid(
        hull_white.coupon_bond_price_t(0.05, 1.0, &[], 0.05),
        "coupon_times is empty",
    );
    expect_invalid(
        hull_white.coupon_bond_price_now(&[], 0.05),
        "coupon_times is empty",
    );
    expect_invalid(
        hull_white.coupon_bond_call_t(0.05, 1.0, 1.5, &[], 0.05, 1.0),
        "coupon_times is empty",
    );
    expect_invalid(
        hull_white.coupon_bond_put_t(0.05, 1.0, 1.5, &[], 0.05, 1.0),
        "coupon_times is empty",
    );
}

#[test]
fn unascending_coupon_times_is_invalid_input() {
    yvf_setup!(hull_white);
    //descending step at index 1
    expect_invalid(
        hull_white.coupon_bond_price_t(0.05, 1.0, &[1.25, 1.0, 1.5], 0.05),
        "coupon_times[1]",
    );
    //a duplicated date is not ascending either; it silently double-weights a leg
    expect_invalid(
        hull_white.coupon_bond_price_t(0.05, 1.0, &[1.25, 1.5, 1.5], 0.05),
        "coupon_times[2]",
    );
}

#[test]
fn coupon_time_at_or_before_valuation_time_is_invalid_input() {
    yvf_setup!(hull_white);
    //exactly on the valuation date: already paid
    expect_invalid(
        hull_white.coupon_bond_price_t(0.05, 1.0, &[1.25, 1.0], 0.05),
        "coupon_times[1]",
    );
    //in the past
    expect_invalid(
        hull_white.coupon_bond_price_t(0.05, 1.0, &[0.5, 2.0], 0.05),
        "coupon_times[0]",
    );
    //for the `now` entry point the valuation time is 0, so a payment at 0 is already made
    expect_invalid(
        hull_white.coupon_bond_price_now(&[0.0, 2.0], 0.05),
        "coupon_times[0]",
    );
}

#[test]
fn expired_swap_is_invalid_input() {
    yvf_setup!(hull_white);
    //matured before the valuation date
    expect_invalid(
        hull_white.swap_price_t(0.05, 1.0, 0.75, 0.25, 0.04),
        "swap_maturity",
    );
    //matures exactly on the valuation date: nothing left to price
    expect_invalid(
        hull_white.swap_price_t(0.05, 1.0, 1.0, 0.25, 0.04),
        "swap_maturity",
    );
}

#[test]
fn non_positive_delta_is_invalid_input() {
    yvf_setup!(hull_white);
    for delta in [0.0, -0.25] {
        expect_invalid(get_coupon_times(4, 1.0, delta), "delta");
        expect_invalid(
            hull_white.swap_price_t(0.05, 1.0, 3.0, delta, 0.04),
            "delta",
        );
        expect_invalid(
            hull_white.swap_price_t_init(0.05, 1.0, 1.0, 8, delta, 0.04),
            "delta",
        );
        expect_invalid(hull_white.swap_rate_t(0.05, 1.0, 8, delta), "delta");
        expect_invalid(
            hull_white.forward_swap_rate_t(0.05, 1.0, 1.5, 8, delta),
            "delta",
        );
        expect_invalid(hull_white.caplet_t(0.05, 1.0, 1.5, delta, 0.04), "delta");
        expect_invalid(
            hull_white.euro_dollar_future_t(0.05, 1.0, 1.5, delta),
            "delta",
        );
        expect_invalid(hull_white.libor_rate_t(0.05, 1.0, delta), "delta");
        expect_invalid(
            hull_white.european_payer_swaption_t(0.05, 1.0, 1.5, 8, delta, 0.04),
            "delta",
        );
    }
}

#[test]
fn zero_period_instrument_is_invalid_input() {
    yvf_setup!(hull_white);
    expect_invalid(
        hull_white.swap_price_t_init(0.05, 1.0, 1.0, 0, 0.25, 0.04),
        "num_swap_payments",
    );
    expect_invalid(
        hull_white.forward_swap_rate_t(0.05, 1.0, 1.5, 0, 0.25),
        "num_swap_payments",
    );
    expect_invalid(
        hull_white.swap_rate_t(0.05, 1.0, 0, 0.25),
        "num_swap_payments",
    );
    expect_invalid(
        hull_white.american_payer_swaption_t(0.05, 1.0, 1.5, 0, 0.25, 0.04, 50),
        "num_swap_payments",
    );
    expect_invalid(
        hull_white.american_receiver_swaption_t(0.05, 1.0, 1.5, 8, 0.25, 0.04, 0),
        "num_steps",
    );
}

#[test]
fn non_finite_input_is_invalid_input_and_is_named() {
    yvf_setup!(hull_white);
    let nan = f64::NAN;
    let inf = f64::INFINITY;
    expect_invalid(hull_white.bond_price_t(nan, 0.0, 2.0), "r_t");
    expect_invalid(hull_white.bond_price_t(0.05, inf, 2.0), "t");
    expect_invalid(hull_white.bond_price_now(nan), "bond_maturity");
    expect_invalid(
        hull_white.swap_price_t(0.05, 1.0, 3.0, 0.25, inf),
        "swap_rate",
    );
    expect_invalid(
        hull_white.coupon_bond_price_t(0.05, 1.0, &[1.5, 2.0], nan),
        "coupon_rate",
    );
    expect_invalid(
        hull_white.american_payer_swaption_t(0.05, 1.0, 1.5, 8, 0.25, nan, 50),
        "swap_rate",
    );
    expect_invalid(hull_white.variance_r(nan, 2.0), "t");
    expect_invalid(hull_white.t_forward_bond_vol(1.0, 2.0, nan), "t_f");
}

#[test]
fn negative_valuation_time_is_invalid_input() {
    yvf_setup!(hull_white);
    expect_invalid(hull_white.bond_price_t(0.05, -1.0, 2.0), "t");
    expect_invalid(hull_white.swap_price_t(0.05, -1.0, 3.0, 0.25, 0.04), "t");
    expect_invalid(
        hull_white.swap_price_t_init(0.05, 1.0, 0.5, 8, 0.25, 0.04),
        "swap_start",
    );
    expect_invalid(
        hull_white.forward_swap_rate_t(0.05, 2.0, 1.5, 8, 0.25),
        "swap_initiation",
    );
}

#[test]
fn bond_option_needs_the_underlying_to_outlive_the_option() {
    yvf_setup!(hull_white);
    //bond maturing at expiry leaves nothing to deliver; before expiry is not deliverable at all
    expect_invalid(
        hull_white.bond_call_t(0.05, 0.5, 1.5, 1.5, 0.98),
        "bond_maturity",
    );
    expect_invalid(hull_white.bond_put_now(1.5, 1.25, 0.98), "bond_maturity");
    //an option with no time left is not an option (zero vol also blows up Black-Scholes)
    expect_invalid(
        hull_white.bond_call_t(0.05, 1.0, 1.0, 2.0, 0.98),
        "option_maturity",
    );
}

#[test]
fn jamshidian_rejects_coupons_paid_before_option_expiry() {
    yvf_setup!(hull_white);
    //The underlying of the decomposed option is the bond *at* expiry.  A coupon paid before
    //expiry is not part of that bond, so pricing it would value a leg that does not exist.
    expect_invalid(
        hull_white.coupon_bond_call_t(0.05, 1.0, 1.5, &[1.25, 1.75, 2.0], 0.05, 1.0),
        "coupon_times[0]",
    );
    expect_invalid(
        hull_white.coupon_bond_put_t(0.05, 1.0, 1.5, &[1.5, 2.0], 0.05, 1.0),
        "coupon_times[0]",
    );
}

#[test]
fn caplet_strike_that_breaks_the_bond_put_transform_is_invalid() {
    yvf_setup!(hull_white);
    //1 + delta * strike == 0 makes the transformed strike 1/0
    expect_invalid(
        hull_white.caplet_t(0.05, 1.0, 1.5, 0.25, -4.0),
        "1 + delta * strike",
    );
    expect_invalid(hull_white.caplet_now(1.5, 0.25, -5.0), "1 + delta * strike");
    //a negative strike deeper than -1/delta is not a caplet
    expect_invalid(
        hull_white.caplet_t(0.05, 1.0, 1.5, 0.25, -40.0),
        "1 + delta * strike",
    );
}

#[test]
fn non_finite_computed_result_surfaces_as_numerical_error() {
    yvf_setup!(hull_white);
    //Inputs are all finite and ordered, but exp() overflows: that is a numerical failure and
    //must not be reported as a price.
    let err = hull_white.bond_price_t(-1e3, 0.0, 10.0).unwrap_err();
    assert!(
        matches!(err, HullWhiteError::NumericalError(_)),
        "expected NumericalError for an overflowing price, got {err:?}"
    );
    assert!(err.to_string().contains("bond_price_t"), "{err}");
}

mod valid_instruments;
