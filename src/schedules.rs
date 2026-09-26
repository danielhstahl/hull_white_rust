//! Payment schedules: which dates an instrument still has left, and how they are indexed.
//!
//! A schedule is generated from a valuation time `t`, a period `delta` and a payment count, as
//! `t + delta, t + 2*delta, ..., t + n*delta`.  [`get_num_remaining_payments`] goes the other
//! way, deriving the count of remaining periods from a maturity — the one place where the
//! answer has to survive the fact that a realistic schedule is not binary-exact.

use crate::error::HullWhiteError;
use crate::validation;

//This is used for when the number of remaining payments must be derived
//eg when the contract already is in the middle of its life.  If the contract
// start is at time 0.0, and current time is 1.1, and delta is 0.25, then the
//next payment is in .15 units of time.

/// Tolerance (relative to the payment count) for deciding that `(maturity - t) / delta` is a whole
/// number of payment periods.  IEEE-754 rounds a realistic schedule's quotient by a few ulps, i.e.
/// ~1e-15 relative: `(1.1 - 0.1) / 0.1 == 9.999999999999998`, not 10.  1e-9 sits ~6 orders of
/// magnitude above that noise while staying far below the smallest gap (~1e-1 periods) that separates
/// two genuinely different schedules, so it never snaps a real off-schedule date onto a period.
const PAYMENT_PERIOD_TOLERANCE: f64 = 1e-9;

pub(crate) fn get_num_remaining_payments(t: f64, maturity: f64, delta: f64) -> (usize, bool) {
    let raw_payments = (maturity - t) / delta;
    //Never test an f64 quotient for exact integrality: `raw == raw.trunc()` only fires when the
    //schedule happens to be binary-representable, so on real schedules `is_exact` flipped and the
    //swap got anchored a whole coupon period early/late (see `swap_price_t`).  Compare against the
    //nearest integer within tolerance instead; below tolerance means the "next payment" is now.
    let whole_payments = raw_payments.round();
    let tolerance = PAYMENT_PERIOD_TOLERANCE * raw_payments.abs().max(1.0);
    if (raw_payments - whole_payments).abs() <= tolerance {
        (whole_payments as usize, true)
    } else {
        ((raw_payments.floor() + 1.0) as usize, false)
    }
}

/// Builds a payment schedule `t + delta, t + 2*delta, ..., t + num_payments*delta`.
///
/// # Errors
///
/// [`HullWhiteError::InvalidInput`] if `t` is negative or `delta` is not strictly positive.  With
/// `delta == 0` every payment landed on `t`, which read as a schedule and priced as a pile of
/// coincident (and, downstream, divided-by) cash flows.
pub fn get_coupon_times(
    num_payments: usize,
    t: f64,
    delta: f64,
) -> Result<Vec<f64>, HullWhiteError> {
    validation::valuation_time(t)?;
    validation::positive("delta", delta)?;
    Ok((1..(num_payments + 1))
        .map(|index| get_time_from_t_index(index, t, delta))
        .collect())
}

pub(crate) fn get_time_from_t_index(index: usize, t: f64, delta: f64) -> f64 {
    t + (index as f64) * delta
}

#[cfg(test)]
mod tests;
