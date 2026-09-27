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
impl<'a> HullWhite<'a> {
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
        (-r_t * self.bond_b(t, bond_maturity) + self.bond_c(t, bond_maturity)).exp()
    }
    /// `dP(t, T)/dr_t = -P(t, T) * B(t, T)`, the Newton step for the critical-rate solve.
    pub(crate) fn bond_price_t_deriv(&self, r_t: f64, t: f64, bond_maturity: f64) -> f64 {
        let b = self.bond_b(t, bond_maturity);
        -self.bond_price_t_raw(r_t, t, bond_maturity) * b
    }
    /// Returns price of a zero coupon bond at current date
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let bond_maturity = 2.0;
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let bond_price = hull_white.bond_price_now(bond_maturity).unwrap();
    /// ```
    pub fn bond_price_now(&self, bond_maturity: f64) -> Result<f64, HullWhiteError> {
        validation::non_negative("bond_maturity", bond_maturity)?;
        validation::finish("bond_price_now", self.bond_price_now_raw(bond_maturity))
    }
    /// Unvalidated counterpart of [`HullWhite::bond_price_now`], for internals that have already
    /// checked their arguments.
    pub(crate) fn bond_price_now_raw(&self, bond_maturity: f64) -> f64 {
        self.curve().discount(bond_maturity)
    }
    /// Price of a coupon bond at some future date `t`, given the short rate observed there.
    ///
    /// Same schedule convention as [`HullWhite::coupon_bond_price_now`]: `coupon_times` holds every
    /// payment date still to come, in strictly increasing order, **with the bond's maturity as the
    /// last element** — that last date is where the par value is added, so a schedule that stops one
    /// payment short repays principal a coupon period early.
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
    /// Price of a coupon bond at the current date.
    ///
    /// # The schedule convention: the last coupon time *is* the bond maturity
    ///
    /// `coupon_times` must contain every remaining payment date **including** the maturity date, in
    /// strictly increasing order.  The kernel adds the `1.0` of par value to the **last** element
    /// of the slice and coupons to all of them, so the maturity is not a separate argument — it is
    /// `coupon_times[coupon_times.len() - 1]`, and that final payment is `coupon_rate + 1.0`.
    ///
    /// Two consequences worth naming, because both used to be documented backwards (the parameter
    /// comment claimed the schedule "does not include the bond_maturity, but the function does
    /// check for that"; it neither excludes it nor checks it):
    ///
    /// * leave the maturity date out and you get a bond whose principal is repaid on the *last
    ///   coupon* date, i.e. one coupon period early, with no diagnostic;
    /// * repeat the maturity date and you are asking to be paid the principal twice — which the
    ///   strictly-increasing requirement on the schedule refuses outright rather than pricing.
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let coupon_rate = 0.05;
    /// //Three dates, quarterly-ish; the last one is the maturity, not "one before it".
    /// let coupon_times = vec![1.0, 1.5, 2.0];
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let price = hull_white
    ///     .coupon_bond_price_now(&coupon_times, coupon_rate)
    ///     .unwrap();
    /// //Discount every coupon off the curve; only the final date also carries the principal.
    /// let coupons_only = coupon_rate
    ///     * (hull_white.bond_price_now(1.0).unwrap()
    ///         + hull_white.bond_price_now(1.5).unwrap()
    ///         + hull_white.bond_price_now(2.0).unwrap());
    /// let expected = coupons_only + 1.0 * hull_white.bond_price_now(2.0).unwrap();
    /// assert!((price - expected).abs() < 1e-12, "{price} vs {expected}");
    /// //So a bond whose schedule is just its maturity is a zero coupon bond, priced identically.
    /// let zero = hull_white.coupon_bond_price_now(&[2.0], 0.0).unwrap();
    /// assert!((zero - hull_white.bond_price_now(2.0).unwrap()).abs() < 1e-12);
    /// //And a repeated date is refused, not silently paid twice.
    /// assert!(
    ///     hull_white
    ///         .coupon_bond_price_now(&[1.0, 1.5, 2.0, 2.0], coupon_rate)
    ///         .is_err()
    /// );
    /// ```
    pub fn coupon_bond_price_now(
        &self,
        //Every remaining payment date, strictly increasing, with the bond maturity as the LAST
        //element: `coupon_bond_generic_now` adds the par value to that last date, so the maturity
        //is a member of the schedule rather than a separate argument (see the doc comment above).
        coupon_times: &[f64],
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
