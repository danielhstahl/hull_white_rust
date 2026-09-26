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
use crate::validation;

impl<'a, T, U> HullWhite<'a, T, U>
where
    T: Fn(f64) -> f64 + std::marker::Sync,
    U: Fn(f64) -> f64 + std::marker::Sync,
{
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
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
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
            self.t_forward_bond_vol(t, option_maturity, bond_maturity)?, //volatility with maturity
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
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
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
            self.t_forward_bond_vol(t, option_maturity, bond_maturity)?, //volatility with maturity
        );
        validation::finish("bond_call_now", price)
    }
    /// Returns price of a call option on a coupon bond at some future time
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //rate at time t
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let option_maturity = 1.5;
    /// let coupon_times = vec![1.75, 2.0, 2.25, 2.5, 2.75, 3.0];
    /// let coupon_rate = 0.05;
    /// let strike = 1.0;
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let bond_call = hull_white.coupon_bond_call_t(r_t, t, option_maturity, &coupon_times, coupon_rate, strike).unwrap();
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
            &|r_t: f64, t: f64, option_maturity: f64, bond_maturity: f64, strike: f64| {
                self.bond_call_t(r_t, t, option_maturity, bond_maturity, strike)
            },
        )
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
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
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
            self.t_forward_bond_vol(t, option_maturity, bond_maturity)?, //volatility with maturity
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
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
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
            self.t_forward_bond_vol(t, option_maturity, bond_maturity)?, //volatility with maturity
        );
        validation::finish("bond_put_now", price)
    }
    /// Returns price of a put option on a coupon bond at some future time
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //rate at time t
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let option_maturity = 1.5;
    /// let coupon_times = vec![1.75, 2.0, 2.25, 2.5, 2.75, 3.0];
    /// let coupon_rate = 0.05;
    /// let strike = 1.0;
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let bond_put = hull_white.coupon_bond_put_t(r_t, t, option_maturity, &coupon_times, coupon_rate, strike).unwrap();
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
            &|r_t: f64, t: f64, option_maturity: f64, bond_maturity: f64, strike: f64| {
                self.bond_put_t(r_t, t, option_maturity, bond_maturity, strike)
            },
        )
    }
}

#[cfg(test)]
mod tests;
