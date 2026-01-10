//! Utility functions for the Hull-White model

/// Helper function to return the maximum of a value and zero
pub fn max_or_zero(v: f64) -> f64 {
    if v > 0.0 {
        v
    } else {
        0.0
    }
}

/// Calculate the payoff for a swaption based on whether it's a payer or receiver
pub fn payoff_swaption(is_payer: bool, swap_value: f64) -> f64 {
    match is_payer {
        true => max_or_zero(swap_value),
        false => max_or_zero(-swap_value),
    }
}

/// Compute the LIBOR rate from two bond prices
pub fn compute_libor_rate(nearest_bond: f64, farthest_bond: f64, tenor: f64) -> f64 {
    (nearest_bond - farthest_bond) / (farthest_bond * tenor)
}

/// Determine the number of remaining payments for a contract
/// 
/// # Arguments
/// * `current_time` - Current time in the contract lifecycle
/// * `maturity` - Maturity time of the contract
/// * `payment_interval` - Time interval between payments
/// 
/// # Returns
/// A tuple containing (number of remaining payments, whether the next payment is exactly at the interval)
pub fn get_num_remaining_payments(current_time: f64, maturity: f64, payment_interval: f64) -> (usize, bool) {
    let raw_payments = (maturity - current_time) / payment_interval;
    
    // If exactly integer, then "next payment" happened already
    if raw_payments == raw_payments.trunc() {
        (raw_payments as usize, true)
    } else {
        ((raw_payments.floor() + 1.0) as usize, false)
    }
}

/// Generate coupon payment times
/// 
/// # Arguments
/// * `num_payments` - Number of coupon payments
/// * `start_time` - Time when coupon payments start
/// * `payment_interval` - Time interval between payments
pub fn get_coupon_times(num_payments: usize, start_time: f64, payment_interval: f64) -> Vec<f64> {
    (1..=(num_payments))
        .map(|index| start_time + (index as f64) * payment_interval)
        .collect()
}

/// Helper function to calculate time from index
pub fn get_time_from_index(index: usize, start_time: f64, payment_interval: f64) -> f64 {
    start_time + (index as f64) * payment_interval
}