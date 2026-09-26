//! Input validation for instrument time coordinates and payment schedules.
//!
//! Every public pricing entry point runs its arguments through these checks before any arithmetic
//! happens.  The contract is deliberately narrow: a bad instrument — an expired contract, an empty
//! payment schedule, a non-positive tenor, a non-finite argument — is reported as
//! [`HullWhiteError::InvalidInput`](crate::HullWhiteError::InvalidInput) naming the offending
//! argument, instead of panicking on a `usize` underflow, dividing by zero, or returning a number
//! that looks like a price.
//!
//! The checks are O(1) (O(number of coupons) for a schedule) and run once per public call, so the
//! inner pricing kernels keep their original shape and cost: they call the `*_raw`/private helpers
//! and never re-validate.
use crate::error::HullWhiteError;

fn invalid(msg: String) -> HullWhiteError {
    HullWhiteError::InvalidInput(msg)
}

/// Anything used as a time, rate or notional has to be an actual number.
pub fn finite(name: &str, value: f64) -> Result<(), HullWhiteError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(invalid(format!("{name} = {value} is not finite")))
    }
}

/// `value` must be strictly positive (payment tenor, tree steps, ...).
pub fn positive(name: &str, value: f64) -> Result<(), HullWhiteError> {
    finite(name, value)?;
    if value > 0.0 {
        Ok(())
    } else {
        Err(invalid(format!("{name} = {value} must be > 0")))
    }
}

/// `value` must be zero or positive.
pub fn non_negative(name: &str, value: f64) -> Result<(), HullWhiteError> {
    finite(name, value)?;
    if value >= 0.0 {
        Ok(())
    } else {
        Err(invalid(format!("{name} = {value} must be >= 0")))
    }
}

/// All times in this crate are measured from "now" (0), so the valuation time cannot be negative.
pub fn valuation_time(t: f64) -> Result<(), HullWhiteError> {
    non_negative("t", t)
}

/// `later` must be at or after `earlier`.  Used where the degenerate (equal) case is legitimate —
/// e.g. pricing a zero coupon bond exactly at its maturity, which is worth par.
pub fn not_before(
    later_name: &str,
    later: f64,
    earlier_name: &str,
    earlier: f64,
) -> Result<(), HullWhiteError> {
    finite(later_name, later)?;
    finite(earlier_name, earlier)?;
    if later >= earlier {
        Ok(())
    } else {
        Err(invalid(format!(
            "{later_name} = {later} must be >= {earlier_name} = {earlier}"
        )))
    }
}

/// `later` must be strictly after `earlier`.  Used where equality means the contract (or a leg of
/// it) has already matured/been paid and so has nothing left to price.
pub fn strictly_after(
    later_name: &str,
    later: f64,
    earlier_name: &str,
    earlier: f64,
) -> Result<(), HullWhiteError> {
    finite(later_name, later)?;
    finite(earlier_name, earlier)?;
    if later > earlier {
        Ok(())
    } else {
        Err(invalid(format!(
            "{later_name} = {later} must be > {earlier_name} = {earlier}"
        )))
    }
}

/// A live payment schedule: non-empty, finite, strictly ascending, and every payment strictly after
/// `t`.  A payment at or before the valuation date has already been made and cannot be represented
/// by the remaining-payment machinery, and an unordered schedule silently re-weights the legs.
pub fn payment_schedule(coupon_times: &[f64], t: f64) -> Result<(), HullWhiteError> {
    if coupon_times.is_empty() {
        return Err(invalid(
            "coupon_times is empty; an instrument with no remaining payments has no price here"
                .to_string(),
        ));
    }
    finite("t", t)?;
    for (index, &coupon_time) in coupon_times.iter().enumerate() {
        let name = format!("coupon_times[{index}]");
        strictly_after(&name, coupon_time, "t", t)?;
        if index > 0 {
            strictly_after(
                &name,
                coupon_time,
                &format!("coupon_times[{}]", index - 1),
                coupon_times[index - 1],
            )?;
        }
    }
    Ok(())
}

/// There is no instrument with zero periods.
pub fn at_least_one(name: &str, count: usize) -> Result<(), HullWhiteError> {
    if count >= 1 {
        Ok(())
    } else {
        Err(invalid(format!("{name} = {count} must be >= 1")))
    }
}

/// An option on a bond (or the unit strike of a swaption) cannot have a negative strike.
pub fn strike(value: f64) -> Result<(), HullWhiteError> {
    non_negative("strike", value)
}

/// The caplet <-> bond-put transformation divides by `1 + delta * strike`, so that quantity has to
/// be strictly positive: a strike at or below `-1/delta` is not a caplet.
pub fn caplet_strike(delta: f64, strike_value: f64) -> Result<(), HullWhiteError> {
    positive("delta", delta)?;
    finite("strike", strike_value)?;
    let one_plus = 1.0 + delta * strike_value;
    if one_plus > 0.0 {
        Ok(())
    } else {
        Err(invalid(format!(
            "1 + delta * strike = {one_plus} must be > 0 (strike = {strike_value}, delta = {delta})"
        )))
    }
}

/// A computed price that is not finite is a numerical failure, not a price of zero.
pub fn finish(what: &str, value: f64) -> Result<f64, HullWhiteError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(HullWhiteError::NumericalError(format!(
            "{what} produced a non-finite result: {value}"
        )))
    }
}

#[cfg(test)]
mod tests;
