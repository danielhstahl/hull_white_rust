//! Vanilla swaps: forward swap rate, swap price, and the European (analytic) swaptions.
//!
//! A swap is priced as the difference of two bond legs — the fixed leg is a coupon bond, the
//! floating leg par at reset — so `swap_price_t(swap_rate)` is zero exactly at the forward swap
//! rate.  [`HullWhite::swap_price_t`] derives the remaining-payment anchor from the maturity
//! (via [`crate::schedules::get_num_remaining_payments`]); [`HullWhite::swap_price_t_init`]
//! takes the anchor as an argument instead, which is the form the tree and the swaptions use.
//!
//! European swaptions here are Black-style on the forward swap rate (Jamshidian/Chen-Holden
//! single-curve form).  American, tree-based swaptions live in [`crate::trees`].

use crate::HullWhite;
use crate::error::HullWhiteError;
use crate::schedules::{get_coupon_times, get_num_remaining_payments};
use crate::validation;

impl<'a, T, U> HullWhite<'a, T, U>
where
    T: Fn(f64) -> f64 + std::marker::Sync,
    U: Fn(f64) -> f64 + std::marker::Sync,
{
    /// Returns forward swap rate at some future time
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //rate at t
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let swap_initiation = 1.5;
    /// let num_swap_payments = 14;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let forward_swap = hull_white.forward_swap_rate_t(r_t,  t, swap_initiation, num_swap_payments, delta).unwrap();
    /// ```
    pub fn forward_swap_rate_t(
        &self,
        r_t: f64,
        t: f64,
        swap_initiation: f64, //must be greater than or equl to t
        num_swap_payments: usize,
        delta: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        validation::not_before("swap_initiation", swap_initiation, "t", t)?;
        validation::at_least_one("num_swap_payments", num_swap_payments)?;
        validation::positive("delta", delta)?;
        let denominator_swap: f64 = (1..(num_swap_payments + 1))
            .map(|curr| self.bond_price_t_raw(r_t, t, swap_initiation + delta * (curr as f64)))
            .sum::<f64>()
            * delta;
        validation::finish(
            "forward_swap_rate_t",
            (self.bond_price_t_raw(r_t, t, swap_initiation)
                - self.bond_price_t_raw(
                    r_t,
                    t,
                    swap_initiation + (num_swap_payments as f64) * delta,
                )) //swap_initiation + (num_swap_payments as f64) * delta)=swap_maturity+delta
                / denominator_swap,
        )
    }
    /// Returns forward swap rate at current time
    ///
    /// Exactly [`HullWhite::forward_swap_rate_t`] at `t = 0` with `r_t = r(0)`
    /// ([`HullWhite::short_rate_now`]).  The two bond legs that make the forward swap rate are
    /// priced with `bond_price_now`, so on this side of the trade the curve — not a state variable
    /// — does all the work.
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let swap_initiation = 1.5;
    /// let num_swap_payments = 14;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let yield_curve = |t:f64|0.05*t + 0.01*t*t; //cumulative yield: the integral of the forward curve
    /// let forward_curve = |t:f64|0.05 + 0.02*t; //instantaneous forward, finite at 0
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let forward_swap = hull_white
    ///     .forward_swap_rate_now(swap_initiation, num_swap_payments, delta)
    ///     .unwrap();
    /// let r0 = hull_white.short_rate_now().unwrap();
    /// let via_t = hull_white
    ///     .forward_swap_rate_t(r0, 0.0, swap_initiation, num_swap_payments, delta)
    ///     .unwrap();
    /// assert!((forward_swap - via_t).abs() < 1e-12, "{forward_swap} vs {via_t}");
    /// assert!(forward_swap.is_finite(), "forward swap rate {forward_swap}");
    /// ```
    pub fn forward_swap_rate_now(
        &self,
        swap_initiation: f64,
        num_swap_payments: usize,
        delta: f64,
    ) -> Result<f64, HullWhiteError> {
        let t = 0.0; //since "now"
        let r_t = self.short_rate_now()?;
        self.forward_swap_rate_t(r_t, t, swap_initiation, num_swap_payments, delta)
    }
    /// Returns swap rate at some future time
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //rate at t
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let num_swap_payments = 16;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let swap_rate = hull_white.swap_rate_t(r_t, t, num_swap_payments, delta).unwrap();
    /// ```
    pub fn swap_rate_t(
        &self,
        r_t: f64,
        t: f64,
        num_swap_payments: usize,
        delta: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::valuation_time(t)?;
        validation::at_least_one("num_swap_payments", num_swap_payments)?;
        validation::positive("delta", delta)?;
        self.forward_swap_rate_t(
            r_t,
            t,
            t, //swap is a forward swap starting at time "0" (t-t=0)
            num_swap_payments,
            delta,
        )
    }
    /// Returns swap rate at current time
    ///
    /// The swap starting today, priced with [`HullWhite::forward_swap_rate_now`] — the `now` twin of
    /// [`HullWhite::swap_rate_t`].
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let num_swap_payments = 16;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let yield_curve = |t:f64|0.05*t + 0.01*t*t; //cumulative yield: the integral of the forward curve
    /// let forward_curve = |t:f64|0.05 + 0.02*t; //instantaneous forward, finite at 0
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let swap_rate = hull_white.swap_rate_now(num_swap_payments, delta).unwrap();
    /// let forward = hull_white.forward_swap_rate_now(0.0, num_swap_payments, delta).unwrap();
    /// assert!((swap_rate - forward).abs() < 1e-12, "{swap_rate} vs {forward}");
    /// ```
    pub fn swap_rate_now(
        &self,
        num_swap_payments: usize,
        delta: f64,
    ) -> Result<f64, HullWhiteError> {
        self.forward_swap_rate_now(0.0, num_swap_payments, delta)
    }
    /// Returns price of a swap at some future time, not necessarily at initiation of the swap
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //current rate
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let swap_maturity = 5.0;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let swap_rate = 0.04; //at initiation, the swap rate is such that the swap has zero value
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let swap = hull_white.swap_price_t(r_t, t, swap_maturity, delta, swap_rate).unwrap();
    /// ```
    pub fn swap_price_t(
        &self,
        r_t: f64,
        t: f64,
        swap_maturity: f64,
        delta: f64,
        swap_rate: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        validation::finite("swap_rate", swap_rate)?;
        validation::positive("delta", delta)?;
        //An expired swap is not a swap.  Without this the derived payment count came out at 0 (or
        //negative-and-saturated), the payment loop vanished, and the function returned the value of
        //a bond leg as if it were a swap price.
        validation::strictly_after("swap_maturity", swap_maturity, "t", t)?;
        let (num_payments, is_exact) = get_num_remaining_payments(t, swap_maturity, delta);
        let swap_start = if is_exact {
            t
        } else {
            swap_maturity - (num_payments as f64 - 1.0) * delta
        };
        self.swap_price_t_init(r_t, t, swap_start, num_payments, delta, swap_rate)
    }
    /// Returns price of a swap at current time, not necessarily at initiation of the swap
    ///
    /// Exactly [`HullWhite::swap_price_t`] at `t = 0` with `r_t = r(0)`
    /// ([`HullWhite::short_rate_now`]), including that function's derivation of the remaining
    /// payment count from `swap_maturity`.  At the current forward swap rate the swap prices at
    /// zero; above it the payer loses.
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let swap_maturity = 5.0;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let yield_curve = |t:f64|0.05*t + 0.01*t*t; //cumulative yield: the integral of the forward curve
    /// let forward_curve = |t:f64|0.05 + 0.02*t; //instantaneous forward, finite at 0
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// //A swap already running: 4 payments left, maturing in 1 year.
    /// let atm = hull_white.swap_rate_now(4, delta).unwrap();
    /// assert!(hull_white.swap_price_now(1.0, delta, atm).unwrap().abs() < 1e-12);
    /// //Paying above the forward is a loss to the payer.
    /// assert!(hull_white.swap_price_now(1.0, delta, atm + 0.01).unwrap() < 0.0);
    /// //Already matured is not a swap.
    /// assert!(hull_white.swap_price_now(0.0, delta, atm).is_err());
    /// ```
    pub fn swap_price_now(
        &self,
        swap_maturity: f64,
        delta: f64,
        swap_rate: f64,
    ) -> Result<f64, HullWhiteError> {
        let t = 0.0; //since "now"
        let r_t = self.short_rate_now()?;
        self.swap_price_t(r_t, t, swap_maturity, delta, swap_rate)
    }
    /// Returns price of a swap at the start of the swap
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //current rate
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let num_swap_payments = 16;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let swap_rate = 0.04;
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let swap = hull_white.swap_price_t_init(r_t, t, t, num_swap_payments, delta, swap_rate).unwrap();
    /// ```
    pub fn swap_price_t_init(
        &self,
        r_t: f64,
        t: f64,
        swap_start: f64, //must be greater than or equal to t, and less than t+delta
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        validation::finite("swap_rate", swap_rate)?;
        validation::not_before("swap_start", swap_start, "t", t)?;
        validation::at_least_one("num_swap_payments", num_swap_payments)?;
        validation::positive("delta", delta)?;
        validation::finish(
            "swap_price_t_init",
            self.swap_price_t_init_raw(r_t, t, swap_start, num_swap_payments, delta, swap_rate),
        )
    }
    /// Returns price of a swap at the start of the swap, as of now
    ///
    /// Exactly [`HullWhite::swap_price_t_init`] at `t = 0` with `r_t = r(0)`
    /// ([`HullWhite::short_rate_now`]), for a swap whose first payment is `delta` away.  Since the
    /// valuation date and the swap start coincide, `swap_start` is `0.0`.
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let num_swap_payments = 16;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let swap_rate = 0.04;
    /// let yield_curve = |t:f64|0.05*t + 0.01*t*t; //cumulative yield: the integral of the forward curve
    /// let forward_curve = |t:f64|0.05 + 0.02*t; //instantaneous forward, finite at 0
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let swap = hull_white.swap_price_now_init(num_swap_payments, delta, swap_rate).unwrap();
    /// let via_t = hull_white
    ///     .swap_price_t_init(
    ///         hull_white.short_rate_now().unwrap(),
    ///         0.0,
    ///         0.0,
    ///         num_swap_payments,
    ///         delta,
    ///         swap_rate,
    ///     )
    ///     .unwrap();
    /// assert!((swap - via_t).abs() < 1e-12, "{swap} vs {via_t}");
    /// //At the rate that par-replaces the float leg, the swap is worth nothing.
    /// let atm = hull_white.swap_rate_now(num_swap_payments, delta).unwrap();
    /// assert!(hull_white.swap_price_now_init(num_swap_payments, delta, atm).unwrap().abs() < 1e-12);
    /// ```
    pub fn swap_price_now_init(
        &self,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
    ) -> Result<f64, HullWhiteError> {
        let t = 0.0; //since "now"
        let r_t = self.short_rate_now()?;
        self.swap_price_t_init(r_t, t, t, num_swap_payments, delta, swap_rate)
    }
    /// Unvalidated swap pricing, used by the swaption tree: the tree callback has to return an `f64`
    /// (`binomial_tree`'s contract), so it cannot propagate a `Result`, and its inputs are checked
    /// once at the public boundary instead of once per node.
    pub(crate) fn swap_price_t_init_raw(
        &self,
        r_t: f64,
        t: f64,
        swap_start: f64,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
    ) -> f64 {
        //open question, should num_swap_payments be num_swap_payments+1??
        let sm_bond: f64 = (1..num_swap_payments)
            .map(|curr| {
                self.bond_price_t_raw(r_t, t, swap_start + delta * (curr as f64))
                    * swap_rate
                    * delta
            })
            .sum();
        self.bond_price_t_raw(r_t, t, swap_start)
            - sm_bond
            - (1.0 + swap_rate * delta)
                * self.bond_price_t_raw(r_t, t, swap_start + delta * (num_swap_payments as f64))
    }
    /// Returns price of a payer swaption at some future time t
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //current rate
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let option_maturity = 2.0;
    /// let num_swap_payments = 16;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let swap_rate = 0.04; //the swap rate is what the payer agrees to pay if option is exercised
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let swaption = hull_white.european_payer_swaption_t(r_t, t, option_maturity, num_swap_payments, delta, swap_rate).unwrap();
    /// ```
    pub fn european_payer_swaption_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
    ) -> Result<f64, HullWhiteError> {
        let coupon_times = get_coupon_times(num_swap_payments, option_maturity, delta)?;
        let strike = 1.0;
        self.coupon_bond_put_t(
            r_t,
            t,
            option_maturity,
            &coupon_times,
            swap_rate * delta,
            strike,
        ) //swaption is equal to put on coupon bond with coupon=swaption swapRate*delta and strike 1.
    }
    /// Returns price of a payer swaption at current time
    ///
    /// Exactly [`HullWhite::european_payer_swaption_t`] at `t = 0` with `r_t = r(0)`
    /// ([`HullWhite::short_rate_now`]): the Jamshidian decomposition of the underlying coupon bond
    /// runs off the initial short rate, which is the only state that exists today.
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let option_maturity = 2.0;
    /// let num_swap_payments = 16;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let yield_curve = |t:f64|0.05*t + 0.01*t*t; //cumulative yield: the integral of the forward curve
    /// let forward_curve = |t:f64|0.05 + 0.02*t; //instantaneous forward, finite at 0
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// //Strike at today's forward swap rate for the same swap.
    /// let atm = hull_white.swap_rate_now(num_swap_payments, delta).unwrap();
    /// let swaption = hull_white
    ///     .european_payer_swaption_now(option_maturity, num_swap_payments, delta, atm)
    ///     .unwrap();
    /// assert!(swaption > 0.0, "at the money payer swaption {swaption}");
    /// //A payer wants to pay below the forward, so a lower strike is worth more.
    /// let cheaper = hull_white
    ///     .european_payer_swaption_now(option_maturity, num_swap_payments, delta, atm - 0.01)
    ///     .unwrap();
    /// assert!(cheaper > swaption, "{cheaper} vs {swaption}");
    /// ```
    pub fn european_payer_swaption_now(
        &self,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
    ) -> Result<f64, HullWhiteError> {
        let t = 0.0; //since "now"
        let r_t = self.short_rate_now()?;
        self.european_payer_swaption_t(r_t, t, option_maturity, num_swap_payments, delta, swap_rate)
    }
    /// Returns price of a payer swaption at some future time t
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //current rate
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let option_maturity = 2.0;
    /// let num_swap_payments = 16;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let swap_rate = 0.04; //the swap rate is what the payer agrees to pay if option is exercised
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let swaption = hull_white.european_receiver_swaption_t(r_t, t, option_maturity, num_swap_payments, delta, swap_rate).unwrap();
    /// ```
    pub fn european_receiver_swaption_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
    ) -> Result<f64, HullWhiteError> {
        let coupon_times = get_coupon_times(num_swap_payments, option_maturity, delta)?;
        let strike = 1.0;
        self.coupon_bond_call_t(
            r_t,
            t,
            option_maturity,
            &coupon_times,
            swap_rate * delta,
            strike,
        ) //swaption is equal to call on coupon bond with coupon=swapRate*delta and strike 1.
    }
    /// Returns price of a receiver swaption at current time
    ///
    /// Exactly [`HullWhite::european_receiver_swaption_t`] at `t = 0` with `r_t = r(0)`
    /// ([`HullWhite::short_rate_now`]).
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let option_maturity = 2.0;
    /// let num_swap_payments = 16;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let yield_curve = |t:f64|0.05*t + 0.01*t*t; //cumulative yield: the integral of the forward curve
    /// let forward_curve = |t:f64|0.05 + 0.02*t; //instantaneous forward, finite at 0
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let atm = hull_white.swap_rate_now(num_swap_payments, delta).unwrap();
    /// let swaption = hull_white
    ///     .european_receiver_swaption_now(option_maturity, num_swap_payments, delta, atm)
    ///     .unwrap();
    /// assert!(swaption > 0.0, "at the money receiver swaption {swaption}");
    /// //A receiver wants to receive above the forward, so a higher strike is worth more.
    /// let richer = hull_white
    ///     .european_receiver_swaption_now(option_maturity, num_swap_payments, delta, atm + 0.01)
    ///     .unwrap();
    /// assert!(richer > swaption, "{richer} vs {swaption}");
    /// ```
    pub fn european_receiver_swaption_now(
        &self,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
    ) -> Result<f64, HullWhiteError> {
        let t = 0.0; //since "now"
        let r_t = self.short_rate_now()?;
        self.european_receiver_swaption_t(
            r_t,
            t,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
        )
    }
}

#[cfg(test)]
mod tests;
