//! Bond pricers: zero coupon and coupon-bearing, at a future date and at "now".
//!
//! [`coupon_bond_generic_t`] / [`coupon_bond_generic_now`] are the cash-flow kernels: they walk a
//! coupon schedule and hand each date to a supplied discounting function, adding the par payment
//! on the last one.  Zero-coupon and coupon-bond (and, through the kernels, option) pricers are
//! the same code with a different kernel.
//!
//! Each public pricer has an infallible `*_raw` twin.  The public entry points validate and the
//! kernels call the raw form, so a schedule is validated once per call rather than once per
//! leg — which is also what lets the Jamshidian closures stay non-`Result`.

use crate::HullWhite;
use crate::curves::{at_t, ct_t};
use crate::error::HullWhiteError;
use crate::validation;

//The coupon-sum kernels below are infallible: an empty schedule sums to 0.0 instead of overflowing
//an index, and `is_last` is derived with `index + 1 == len` so there is no `len() - 1` anywhere to
//underflow.  Emptiness and every other bad-instrument case is the public boundary's job (see
//`validation`), which is what lets the root-finding closures below be infallible too.

pub(crate) fn coupon_bond_generic_t(
    r_t: f64,
    t: f64,
    coupon_times: &[f64], //includes bond_maturity
    coupon_rate: f64,
    generic_fn: &impl Fn(f64, f64, f64) -> f64,
) -> f64 {
    let par_value = 1.0; //without loss of generality
    let last_index = coupon_times.len();
    coupon_times
        .iter()
        .enumerate()
        .map(|(index, coupon_time)| {
            let is_last = index + 1 == last_index;
            (coupon_rate + if is_last { par_value } else { 0.0 }) * generic_fn(r_t, t, *coupon_time)
        })
        .sum()
}
pub(crate) fn coupon_bond_generic_now(
    coupon_times: &[f64], //includes bond_maturity
    coupon_rate: f64,
    generic_fn: &impl Fn(f64) -> f64,
) -> f64 {
    let par_value = 1.0; //without loss of generality
    let last_index = coupon_times.len();
    coupon_times
        .iter()
        .enumerate()
        .map(|(index, coupon_time)| {
            let is_last = index + 1 == last_index;
            (coupon_rate + if is_last { par_value } else { 0.0 }) * generic_fn(*coupon_time)
        })
        .sum()
}
impl<'a, T, U> HullWhite<'a, T, U>
where
    T: Fn(f64) -> f64 + std::marker::Sync,
    U: Fn(f64) -> f64 + std::marker::Sync,
{
    /// Returns price of a zero coupon bond at some future date
    /// given the interest rate at that future date
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //instantaneous rate at date t
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let bond_maturity = 2.0;
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let bond_price = hull_white.bond_price_t(r_t, t, bond_maturity).unwrap();
    /// ```
    pub fn bond_price_t(
        &self,
        r_t: f64,
        t: f64,
        bond_maturity: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        //bond_maturity == t is legitimate: the bond is at par on its maturity date.
        validation::not_before("bond_maturity", bond_maturity, "t", t)?;
        validation::finish("bond_price_t", self.bond_price_t_raw(r_t, t, bond_maturity))
    }
    /// Unvalidated bond price.  The coupon/option kernels price a whole schedule per node and are
    /// fed arguments already checked at the public boundary, so they use this instead of paying for
    /// (and having to propagate) the checks again.
    pub(crate) fn bond_price_t_raw(&self, r_t: f64, t: f64, bond_maturity: f64) -> f64 {
        (-r_t * at_t(self.a, t, bond_maturity)
            + ct_t(
                self.a,
                self.sigma,
                t,
                bond_maturity,
                self.yield_curve,
                self.forward_curve,
            ))
        .exp()
    }
    //used for newton's method
    pub(crate) fn bond_price_t_deriv(&self, r_t: f64, t: f64, bond_maturity: f64) -> f64 {
        let at_t_c = at_t(self.a, t, bond_maturity);
        -(-r_t * at_t_c
            + ct_t(
                self.a,
                self.sigma,
                t,
                bond_maturity,
                self.yield_curve,
                self.forward_curve,
            ))
        .exp()
            * at_t_c
    }
    /// Returns price of a zero coupon bond at current date
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let bond_maturity = 2.0;
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let bond_price = hull_white.bond_price_now(bond_maturity).unwrap();
    /// ```
    pub fn bond_price_now(&self, bond_maturity: f64) -> Result<f64, HullWhiteError> {
        validation::non_negative("bond_maturity", bond_maturity)?;
        validation::finish("bond_price_now", self.bond_price_now_raw(bond_maturity))
    }
    /// Unvalidated counterpart of [`HullWhite::bond_price_now`], for internals that have already
    /// checked their arguments.
    pub(crate) fn bond_price_now_raw(&self, bond_maturity: f64) -> f64 {
        (-(self.yield_curve)(bond_maturity)).exp()
    }
    /// Returns price of a coupon bond at some future date
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //instantaneous rate at date t
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let coupon_times = vec![1.25, 1.5, 1.75, 2.0]; //measure time from now (0), all should be greater than t.  Final coupon is the bond maturity
    /// let coupon_rate = 0.05;
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let bond_price = hull_white.coupon_bond_price_t(r_t, t, &coupon_times, coupon_rate).unwrap();
    /// ```
    pub fn coupon_bond_price_t(
        &self,
        r_t: f64,
        t: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        validation::finite("coupon_rate", coupon_rate)?;
        validation::payment_schedule(coupon_times, t)?;
        validation::finish(
            "coupon_bond_price_t",
            coupon_bond_generic_t(
                r_t,
                t,
                coupon_times,
                coupon_rate,
                &|r_t: f64, t: f64, bond_maturity: f64| {
                    self.bond_price_t_raw(r_t, t, bond_maturity)
                },
            ),
        )
    }
    pub(crate) fn coupon_bond_price_t_deriv(
        &self,
        r_t: f64,
        t: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
    ) -> f64 {
        coupon_bond_generic_t(
            r_t,
            t,
            coupon_times,
            coupon_rate,
            &|r_t: f64, t: f64, bond_maturity: f64| self.bond_price_t_deriv(r_t, t, bond_maturity),
        )
    }
    /// Returns price of a coupon bond at current date
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let coupon_times = vec![1.25, 1.5, 1.75, 2.0]; //measure time from now (0), all should be greater than t.  Final coupon is the bond maturity
    /// let coupon_rate = 0.05;
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let bond_price = hull_white.coupon_bond_price_now(&coupon_times, coupon_rate).unwrap();
    /// ```
    pub fn coupon_bond_price_now(
        &self,
        coupon_times: &[f64], //does not include the bond_maturity, but the function does check for that
        coupon_rate: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::finite("coupon_rate", coupon_rate)?;
        //"now" is t = 0 for this entry point.
        validation::payment_schedule(coupon_times, 0.0)?;
        validation::finish(
            "coupon_bond_price_now",
            coupon_bond_generic_now(coupon_times, coupon_rate, &|bond_maturity: f64| {
                self.bond_price_now_raw(bond_maturity)
            }),
        )
    }
}

#[cfg(test)]
mod tests;
