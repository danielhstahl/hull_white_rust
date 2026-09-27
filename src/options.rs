//! Option pricers that are plain Black-Scholes on a zero coupon bond, plus the coupon-bond
//! option entry points that delegate to the Jamshidian decomposition in
//! [`crate::jamshidian`].
//!
//! A zero-coupon bond option needs no root finding: under the `option_maturity`-forward measure
//! the bond is lognormal, so [`HullWhite::bond_call_t`] / [`HullWhite::bond_put_t`] are one
//! Black-Scholes call with the bond as underlying, the `option_maturity` bond as discount factor
//! and `t_forward_bond_vol` as volatility.  A coupon bond is a *sum* of such bonds and is not
//! lognormal, which is what the decomposition in [`crate::jamshidian`] exists to handle.

use crate::HullWhite;
use crate::error::HullWhiteError;
use crate::jamshidian::Side;
use crate::validation;

impl<'a> HullWhite<'a> {
    /// Returns price of a call option on zero coupon bond at some future time
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //rate at time t
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let option_maturity = 1.5;
    /// let bond_maturity = 2.0;
    /// let strike = 0.98;
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let bond_call = hull_white.bond_call_t(r_t, t, option_maturity, bond_maturity, strike).unwrap();
    /// ```
    pub fn bond_call_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        bond_maturity: f64,
        strike: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        validation::strictly_after("option_maturity", option_maturity, "t", t)?;
        validation::strictly_after(
            "bond_maturity",
            bond_maturity,
            "option_maturity",
            option_maturity,
        )?;
        validation::strike(strike)?;
        let price = black_scholes::call_discount(
            self.bond_price_t_raw(r_t, t, bond_maturity), //underlying
            strike,
            self.bond_price_t_raw(r_t, t, option_maturity), //discount
            self.t_forward_bond_vol(t, option_maturity, bond_maturity)?, //volatility of the deliverable bond
        );
        validation::finish("bond_call_t", price)
    }
    /// Returns price of a call option on zero coupon bond at current time
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let option_maturity = 1.5;
    /// let bond_maturity = 2.0;
    /// let strike = 0.98;
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let bond_call = hull_white.bond_call_now(option_maturity, bond_maturity, strike).unwrap();
    /// ```
    pub fn bond_call_now(
        &self,
        option_maturity: f64,
        bond_maturity: f64,
        strike: f64,
    ) -> Result<f64, HullWhiteError> {
        let t = 0.0; //since "now"
        validation::strictly_after("option_maturity", option_maturity, "t", t)?;
        validation::strictly_after(
            "bond_maturity",
            bond_maturity,
            "option_maturity",
            option_maturity,
        )?;
        validation::strike(strike)?;
        let price = black_scholes::call_discount(
            self.bond_price_now_raw(bond_maturity), //underlying
            strike,
            self.bond_price_now_raw(option_maturity), //discount
            self.t_forward_bond_vol(t, option_maturity, bond_maturity)?, //volatility of the deliverable bond
        );
        validation::finish("bond_call_now", price)
    }
    /// Returns price of a call option on a coupon bond at some future time
    ///
    /// `coupon_times` is the whole payment schedule of the bond (the final entry is the maturity
    /// payment) and it may straddle the option's expiry: the underlying of a coupon-bond option is
    /// the bond *as it stands on the expiry date*, so
    ///
    /// * coupons paid **strictly before** `option_maturity` are not part of the underlying and are
    ///   dropped -- the option holder never receives them;
    /// * a coupon falling **exactly on** `option_maturity` is cash the holder does receive on the
    ///   expiry date, and is subtracted from the strike (it is a constant at expiry, not a
    ///   function of the rate);
    /// * coupons **strictly after** `option_maturity` are the residual bond that gets decomposed.
    ///
    /// A schedule with nothing strictly after `option_maturity` is refused: the bond is settled
    /// before the option can be exercised, so there is nothing to deliver.  Dropping the
    /// pre-expiry coupons cannot change the price, so the full schedule and the post-expiry tail
    /// alone price identically (see the example).  The full reasoning is in the `jamshidian`
    /// module (`src/jamshidian.rs`).
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //rate at time t
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let option_maturity = 1.5;
    /// //a 2.5y bond with coupons before, on, and after the 1.5y expiry
    /// let coupon_times = vec![1.25, 1.5, 1.75, 2.0, 2.25, 2.5];
    /// let coupon_rate = 0.05;
    /// let strike = 1.0;
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let call = hull_white.coupon_bond_call_t(r_t, t, option_maturity, &coupon_times, coupon_rate, strike).unwrap();
    /// //Dropping the one coupon paid strictly before expiry (1.25) cannot change the price: it is
    /// //not part of the bond the holder gets on the expiry date.
    /// let tail = vec![1.5, 1.75, 2.0, 2.25, 2.5];
    /// let dropped = hull_white.coupon_bond_call_t(r_t, t, option_maturity, &tail, coupon_rate, strike).unwrap();
    /// assert_eq!(call, dropped);
    /// //The coupon paid *on* the expiry date is cash worth `coupon_rate` there, so on the residual
    /// //bond it shows up as a strike reduction.
    /// let residual = vec![1.75, 2.0, 2.25, 2.5];
    /// let folded = hull_white.coupon_bond_call_t(r_t, t, option_maturity, &residual, coupon_rate, strike - coupon_rate).unwrap();
    /// assert_eq!(call, folded);
    /// ```
    pub fn coupon_bond_call_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
    ) -> Result<f64, HullWhiteError> {
        self.coupon_bond_option_generic_t(
            r_t,
            t,
            option_maturity,
            coupon_times,
            coupon_rate,
            strike,
            Side::Call,
            &|r_t: f64, t: f64, option_maturity: f64, bond_maturity: f64, strike: f64| {
                self.bond_call_t(r_t, t, option_maturity, bond_maturity, strike)
            },
        )
    }
    /// Returns price of a call option on a coupon bond at current time
    ///
    /// Exactly [`HullWhite::coupon_bond_call_t`] run at `t = 0` with `r_t = r(0)`, the model's own
    /// initial short rate ([`HullWhite::short_rate_now`]).  The Jamshidian solve needs a rate to
    /// stand at the option's expiry, and that rate is conditioned on the state today; at `t = 0` the
    /// only state there is is the calibrated curve, so the argument disappears rather than becoming
    /// a guess.  Schedule conventions — pre-expiry coupons dropped, coupons on the expiry date folded
    /// into the strike — are unchanged from [`HullWhite::coupon_bond_call_t`].
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let option_maturity = 1.5;
    /// //a 2.5y bond whose first coupon falls after the 1.5y expiry
    /// let coupon_times = vec![1.75, 2.0, 2.25, 2.5];
    /// let coupon_rate = 0.05;
    /// let strike = 1.0;
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let call = hull_white.coupon_bond_call_now(option_maturity, &coupon_times, coupon_rate, strike).unwrap();
    /// //Same number as the `t` form seeded with the initial short rate.
    /// let r0 = hull_white.short_rate_now().unwrap();
    /// let via_t = hull_white
    ///     .coupon_bond_call_t(r0, 0.0, option_maturity, &coupon_times, coupon_rate, strike)
    ///     .unwrap();
    /// assert!((call - via_t).abs() < 1e-12, "{call} vs {via_t}");
    /// assert!(call > 0.0, "call {call}");
    /// ```
    pub fn coupon_bond_call_now(
        &self,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
    ) -> Result<f64, HullWhiteError> {
        let t = 0.0; //since "now"
        let r_t = self.short_rate_now()?;
        self.coupon_bond_call_t(r_t, t, option_maturity, coupon_times, coupon_rate, strike)
    }
    /// Returns price of a put option on zero coupon bond at some future time
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //rate at time t
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let option_maturity = 1.5;
    /// let bond_maturity = 2.0;
    /// let strike = 0.98;
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let bond_put = hull_white.bond_put_t(r_t, t, option_maturity, bond_maturity, strike).unwrap();
    /// ```
    pub fn bond_put_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        bond_maturity: f64,
        strike: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        validation::strictly_after("option_maturity", option_maturity, "t", t)?;
        validation::strictly_after(
            "bond_maturity",
            bond_maturity,
            "option_maturity",
            option_maturity,
        )?;
        validation::strike(strike)?;
        let price = black_scholes::put_discount(
            self.bond_price_t_raw(r_t, t, bond_maturity), //underlying
            strike,
            self.bond_price_t_raw(r_t, t, option_maturity), //discount
            self.t_forward_bond_vol(t, option_maturity, bond_maturity)?, //volatility of the deliverable bond
        );
        validation::finish("bond_put_t", price)
    }
    /// Returns price of a put option on zero coupon bond at current time
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let option_maturity = 1.5;
    /// let bond_maturity = 2.0;
    /// let strike = 0.98;
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let bond_put = hull_white.bond_put_now(option_maturity, bond_maturity, strike).unwrap();
    /// ```
    pub fn bond_put_now(
        &self,
        option_maturity: f64,
        bond_maturity: f64,
        strike: f64,
    ) -> Result<f64, HullWhiteError> {
        let t = 0.0; //since "now"
        validation::strictly_after("option_maturity", option_maturity, "t", t)?;
        validation::strictly_after(
            "bond_maturity",
            bond_maturity,
            "option_maturity",
            option_maturity,
        )?;
        validation::strike(strike)?;
        let price = black_scholes::put_discount(
            self.bond_price_now_raw(bond_maturity), //underlying
            strike,
            self.bond_price_now_raw(option_maturity), //discount
            self.t_forward_bond_vol(t, option_maturity, bond_maturity)?, //volatility of the deliverable bond
        );
        validation::finish("bond_put_now", price)
    }
    /// Returns price of a put option on a coupon bond at some future time
    ///
    /// Same schedule convention as [`HullWhite::coupon_bond_call_t`]: the underlying is the bond
    /// *as it stands on the option's expiry date*.  Coupons strictly before `option_maturity` are
    /// dropped (the holder never receives them), a coupon falling exactly on `option_maturity` is
    /// cash on the expiry date and is subtracted from the strike, and only coupons strictly after
    /// `option_maturity` are decomposed.  A schedule with nothing strictly after the expiry is
    /// refused.  See the `jamshidian` module (`src/jamshidian.rs`) for the reasoning.
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //rate at time t
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let option_maturity = 1.5;
    /// //a 2.5y bond with coupons before, on, and after the 1.5y expiry
    /// let coupon_times = vec![1.25, 1.5, 1.75, 2.0, 2.25, 2.5];
    /// let coupon_rate = 0.05;
    /// let strike = 1.0;
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let bond_put = hull_white.coupon_bond_put_t(r_t, t, option_maturity, &coupon_times, coupon_rate, strike).unwrap();
    /// //The pre-expiry coupon is not part of the underlying: dropping it leaves the price alone.
    /// let tail = vec![1.5, 1.75, 2.0, 2.25, 2.5];
    /// let dropped = hull_white.coupon_bond_put_t(r_t, t, option_maturity, &tail, coupon_rate, strike).unwrap();
    /// assert_eq!(bond_put, dropped);
    /// ```
    pub fn coupon_bond_put_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
    ) -> Result<f64, HullWhiteError> {
        self.coupon_bond_option_generic_t(
            r_t,
            t,
            option_maturity,
            coupon_times,
            coupon_rate,
            strike,
            Side::Put,
            &|r_t: f64, t: f64, option_maturity: f64, bond_maturity: f64, strike: f64| {
                self.bond_put_t(r_t, t, option_maturity, bond_maturity, strike)
            },
        )
    }
    /// Returns price of a put option on a coupon bond at current time
    ///
    /// Exactly [`HullWhite::coupon_bond_put_t`] run at `t = 0` with `r_t = r(0)`
    /// ([`HullWhite::short_rate_now`]); see [`HullWhite::coupon_bond_call_now`] for why the rate
    /// argument drops out and the same schedule conventions as
    /// [`HullWhite::coupon_bond_put_t`].
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let option_maturity = 1.5;
    /// let coupon_times = vec![1.75, 2.0, 2.25, 2.5];
    /// let coupon_rate = 0.05;
    /// let strike = 1.0;
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let put = hull_white.coupon_bond_put_now(option_maturity, &coupon_times, coupon_rate, strike).unwrap();
    /// let r0 = hull_white.short_rate_now().unwrap();
    /// let via_t = hull_white
    ///     .coupon_bond_put_t(r0, 0.0, option_maturity, &coupon_times, coupon_rate, strike)
    ///     .unwrap();
    /// assert!((put - via_t).abs() < 1e-12, "{put} vs {via_t}");
    /// //Put-call parity on the Jamshidian decomposition: C - P is the forward bond leg, with no
    /// //option value left in it.  Everything settles after expiry, so the deliverable is the whole
    /// //bond and nothing was folded out of the strike.
    /// let bond = hull_white.coupon_bond_price_now(&coupon_times, coupon_rate).unwrap();
    /// let discount = hull_white.bond_price_now(option_maturity).unwrap();
    /// let call = hull_white
    ///     .coupon_bond_call_now(option_maturity, &coupon_times, coupon_rate, strike)
    ///     .unwrap();
    /// assert!(
    ///     (call - put - (bond - strike * discount)).abs() < 1e-11,
    ///     "{} - {} vs {} - {} * {}",
    ///     call,
    ///     put,
    ///     bond,
    ///     strike,
    ///     discount
    /// );
    /// ```
    pub fn coupon_bond_put_now(
        &self,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
    ) -> Result<f64, HullWhiteError> {
        let t = 0.0; //since "now"
        let r_t = self.short_rate_now()?;
        self.coupon_bond_put_t(r_t, t, option_maturity, coupon_times, coupon_rate, strike)
    }
}

#[cfg(test)]
mod tests;
