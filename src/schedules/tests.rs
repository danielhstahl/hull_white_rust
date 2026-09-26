//! Unit tests for schedule generation and the remaining-payment count.
//!
//! The payment-count tests are the interesting ones: the count is derived from a float quotient,
//! so both directions of error (a whole schedule read as off-schedule, and an off-schedule date
//! snapped onto a period) are pinned here, plus the downstream consequence for the swap anchor.

use approx::*;

use super::{get_coupon_times, get_num_remaining_payments, get_time_from_t_index};

#[test]
fn test_get_num_payments_if_exact_integer() {
    let t = 0.5;
    let maturity = 2.0;
    let delta = 0.25;
    let (num_payments, is_exact) = get_num_remaining_payments(t, maturity, delta);
    assert_eq!(num_payments, 6);
    assert!(is_exact);
    let next_exchange_date = maturity - (num_payments as f64 - 1.0) * delta;
    assert_abs_diff_eq!(next_exchange_date, t + delta, epsilon = 0.0000001);
}

#[test]
fn test_get_num_payments_not_even() {
    let t = 0.5;
    let maturity = 2.0;
    let delta = 0.4;
    let (num_payments, is_exact) = get_num_remaining_payments(t, maturity, delta);
    assert_eq!(num_payments, 4);
    assert_eq!(is_exact, false);
    let next_exchange_date = maturity - (num_payments as f64 - 1.0) * delta;
    assert_abs_diff_eq!(next_exchange_date, 0.8, epsilon = 0.0000001);
}

/// Whole-number-of-period schedules that are NOT binary-representable.  The f64 quotient drifts to
/// either side of the integer, and each direction broke differently under `== trunc()`:
///   low  e.g. (1.1 - 0.1) / 0.1   = 9.999999999999998  -> is_exact flipped false, so the swap
///        got anchored at maturity - (n-1)*delta = 0.2 = t + delta, one period late
///   high e.g. (1.3000000000000003 - 0.1) / 0.1 = 12.000000000000002 -> floor+1 gave n+1,
///        one spurious payment period
#[test]
fn num_remaining_payments_tolerates_whole_schedules_that_are_not_binary_exact() {
    let whole_schedules: [(f64, f64, f64, usize); 8] = [
        (0.1, 1.1, 0.1, 10), //  9.999999999999998
        (0.1, 2.0, 0.1, 19), // 18.999999999999996
        (0.3, 2.3, 0.2, 10),
        (0.25, 2.25, 0.5, 4),
        (0.1, 0.30000000000000004, 0.1, 2), //  2.0000000000000004
        (0.1, 0.4, 0.1, 3),                 //  3.0000000000000004
        (0.1, 0.7000000000000001, 0.1, 6),  //  6.000000000000001
        (0.1, 1.3000000000000003, 0.1, 12), // 12.000000000000002
    ];
    for (t, maturity, delta, expected_payments) in whole_schedules {
        let raw = (maturity - t) / delta;
        let (num_payments, is_exact) = get_num_remaining_payments(t, maturity, delta);
        assert_eq!(
            num_payments, expected_payments,
            "({maturity} - {t}) / {delta} = {raw:?} should be {expected_payments} payments"
        );
        assert!(
            is_exact,
            "({maturity} - {t}) / {delta} = {raw:?} is a whole schedule and must read as exact"
        );
        //Whole schedule => the anchor the swap derives must be t itself (payments at t+delta ...).
        let next_exchange_date = maturity - (num_payments as f64 - 1.0) * delta;
        assert_abs_diff_eq!(next_exchange_date, t + delta, epsilon = 1e-9);
    }
}

/// Guard the other way: a mid-life schedule that really is between payment dates must NOT be
/// snapped onto a period boundary by the new tolerance.
#[test]
fn num_remaining_payments_does_not_snap_genuinely_off_schedule_dates() {
    //      t, maturity, delta, expected payments, expected first remaining payment
    let mid_life: [(f64, f64, f64, usize, f64); 4] = [
        (0.5, 2.0, 0.4, 4, 0.8),    //  3.75
        (0.1, 1.15, 0.1, 11, 0.15), // 10.5
        (0.0, 1.0, 0.3, 4, 0.1),    //  3.3333333333333335
        (0.05, 1.0, 0.3, 4, 0.1),   //  3.1666666666666665
    ];
    for (t, maturity, delta, expected_payments, expected_first_payment) in mid_life {
        let raw = (maturity - t) / delta;
        let (num_payments, is_exact) = get_num_remaining_payments(t, maturity, delta);
        assert_eq!(
            num_payments, expected_payments,
            "({maturity} - {t}) / {delta} = {raw:?} should be {expected_payments} payments"
        );
        assert!(
            !is_exact,
            "({maturity} - {t}) / {delta} = {raw:?} is off-schedule and must not be exact"
        );
        let next_exchange_date = maturity - (num_payments as f64 - 1.0) * delta;
        assert_abs_diff_eq!(next_exchange_date, expected_first_payment, epsilon = 1e-9);
    }
}

#[test]
fn test_get_time_from_t_index() {
    let t = 0.5;
    let delta = 0.4;
    let index = 3;
    let time = get_time_from_t_index(index, t, delta);
    assert_abs_diff_eq!(time, 1.7, epsilon = 0.0000001);
}

#[test]
fn test_get_time_from_t_index_with_zero_index() {
    let t = 0.5;
    let delta = 0.4;
    let index = 0;
    let time = get_time_from_t_index(index, t, delta);
    assert_eq!(time, t);
}

#[test]
fn test_get_coupon_times() {
    let num_payments = 5;
    let t = 1.0;
    let delta = 0.25;
    let coupon_times = get_coupon_times(num_payments, t, delta).unwrap();
    let expected_coupon_times = vec![1.25, 1.5, 1.75, 2.0, 2.25];
    coupon_times
        .iter()
        .zip(expected_coupon_times.iter())
        .for_each(|(actual, expected)| assert_eq!(actual, expected))
}

#[test]
fn test_get_coupon_times_no_payments() {
    let num_payments = 0;
    let t = 1.0;
    let delta = 0.25;
    let coupon_times = get_coupon_times(num_payments, t, delta).unwrap();
    assert_eq!(coupon_times.len(), 0);
}
