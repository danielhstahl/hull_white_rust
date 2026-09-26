//! Simple-rate instruments: caplets, Eurodollar futures and forward spot Libor.
//!
//! These are all translations of bond prices into a simple (non-compounding) rate over a tenor
//! `delta`, and the option on such a rate is priced by converting the caplet into a put on the
//! underlying deposit bond — `1 + delta * K` has to stay positive for that map to exist, which is
//! what [`crate::validation::caplet_strike`] guards.  The futures price differs from the forward
//! rate by the convexity term built from [`crate::curves::gamma_edf`].

use crate::HullWhite;
use crate::curves::{compute_libor_rate, edf_compute, gamma_edf};
use crate::error::HullWhiteError;
use crate::validation;

impl<'a, T, U> HullWhite<'a, T, U>
where
    T: Fn(f64) -> f64 + std::marker::Sync,
    U: Fn(f64) -> f64 + std::marker::Sync,
{
    /// Returns price of a caplet at current time
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let option_maturity = 1.5;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let strike = 0.04;
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let caplet = hull_white.caplet_now(option_maturity, delta, strike).unwrap();
    /// ```
    pub fn caplet_now(
        &self,
        option_maturity: f64,
        delta: f64,
        strike: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::strictly_after("option_maturity", option_maturity, "t", 0.0)?;
        validation::caplet_strike(delta, strike)?;
        self.bond_put_now(
            option_maturity,
            option_maturity + delta,
            1.0 / (delta * strike + 1.0),
        )
        .map(|put| (strike * delta + 1.0) * put)
    }
    /// Returns price of a caplet at some future time
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //rate at time t
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let option_maturity = 1.5;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let strike = 0.04;
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let caplet = hull_white.caplet_t(r_t, t, option_maturity, delta, strike).unwrap();
    /// ```
    pub fn caplet_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        delta: f64,
        strike: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        validation::strictly_after("option_maturity", option_maturity, "t", t)?;
        validation::caplet_strike(delta, strike)?;
        self.bond_put_t(
            r_t,
            t,
            option_maturity,
            option_maturity + delta,
            1.0 / (delta * strike + 1.0),
        )
        .map(|put| (strike * delta + 1.0) * put)
    }
    /// Returns price of a Euro Dollar Future at some future time
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; // rate at time t
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let option_maturity = 1.5;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let edf = hull_white.euro_dollar_future_t(r_t,  t, option_maturity, delta).unwrap();
    /// ```
    pub fn euro_dollar_future_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        delta: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        validation::strictly_after("option_maturity", option_maturity, "t", t)?;
        validation::positive("delta", delta)?;
        let gamma = gamma_edf(self.a, self.sigma, t, option_maturity, delta);
        validation::finish(
            "euro_dollar_future_t",
            edf_compute(
                self.bond_price_t_raw(r_t, t, option_maturity),
                self.bond_price_t_raw(r_t, t, option_maturity + delta),
                gamma,
                delta,
            ),
        )
    }
    /// Returns price of a Euro Dollar Future at some future time
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let option_maturity = 1.5;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let edf = hull_white.euro_dollar_future_now(option_maturity, delta).unwrap();
    /// ```
    pub fn euro_dollar_future_now(
        &self,
        option_maturity: f64,
        delta: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::strictly_after("option_maturity", option_maturity, "t", 0.0)?;
        validation::positive("delta", delta)?;
        let gamma = gamma_edf(self.a, self.sigma, 0.0, option_maturity, delta);
        validation::finish(
            "euro_dollar_future_now",
            edf_compute(
                self.bond_price_now_raw(option_maturity),
                self.bond_price_now_raw(option_maturity + delta),
                gamma,
                delta,
            ),
        )
    }
    /// Returns forward Libor rate at some future time
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //current rate
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let maturity = 1.5;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let forward_libor = hull_white.forward_libor_rate_t(r_t, t, maturity, delta).unwrap();
    /// ```
    pub fn forward_libor_rate_t(
        &self,
        r_t: f64,
        t: f64,
        maturity: f64,
        delta: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        //maturity == t is the spot fixing (what `libor_rate_t` asks for).
        validation::not_before("maturity", maturity, "t", t)?;
        validation::positive("delta", delta)?;
        let nearest_bond = self.bond_price_t_raw(r_t, t, maturity);
        let farthest_bond = self.bond_price_t_raw(r_t, t, maturity + delta);
        validation::finish(
            "forward_libor_rate_t",
            compute_libor_rate(nearest_bond, farthest_bond, delta),
        )
    }

    /// Returns forward Libor rate at current time
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let maturity = 1.5;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let forward_libor = hull_white.forward_libor_rate_now(maturity, delta).unwrap();
    /// ```
    pub fn forward_libor_rate_now(&self, maturity: f64, delta: f64) -> Result<f64, HullWhiteError> {
        validation::non_negative("maturity", maturity)?;
        validation::positive("delta", delta)?;
        let nearest_bond = self.bond_price_now_raw(maturity);
        let farthest_bond = self.bond_price_now_raw(maturity + delta);
        validation::finish(
            "forward_libor_rate_now",
            compute_libor_rate(nearest_bond, farthest_bond, delta),
        )
    }
    /// Returns Libor rate at some future time
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //current rate
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let libor = hull_white.libor_rate_t(r_t, t, delta).unwrap();
    /// ```
    pub fn libor_rate_t(&self, r_t: f64, t: f64, delta: f64) -> Result<f64, HullWhiteError> {
        validation::valuation_time(t)?;
        validation::positive("delta", delta)?;
        self.forward_libor_rate_t(r_t, t, t, delta)
    }
}

#[cfg(test)]
mod tests;
