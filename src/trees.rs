//! Short-rate tree plumbing: the European tree check and the American (early-exercise) swaptions.
//!
//! Both walk a trinomial-style Black-Vasicek tree from `binomial_tree`, with the model's `phi(t)`
//! cached per time step so the state is the shifted rate `y = r - phi`.  Tree time `tau` means
//! *absolute* model time `t + tau` everywhere: the option is priced at `t`, the tree runs for
//! `option_maturity - t`, and every call into `phi` or a swap leg adds `t` back.
//!
//! [`HullWhite::european_swaption_tree`] is `#[cfg(test)]`: it exists to check the analytic
//! swaption in [`crate::swaps`] against the same tree machinery the American pricer uses.

use crate::HullWhite;
use crate::error::HullWhiteError;
use crate::validation;

fn max_or_zero(v: f64) -> f64 {
    if v > 0.0 { v } else { 0.0 }
}
fn payoff_swaption(is_payer: bool, swp: f64) -> f64 {
    match is_payer {
        true => max_or_zero(swp),
        false => max_or_zero(-swp),
    }
}
impl<'a, T, U> HullWhite<'a, T, U>
where
    T: Fn(f64) -> f64 + std::marker::Sync,
    U: Fn(f64) -> f64 + std::marker::Sync,
{
    fn american_swaption(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64, //tenor of simple yield
        swap_rate: f64,
        is_payer: bool,
        num_steps: usize,
    ) -> f64 {
        let alpha_div_sigma = |_t_step: f64, curr_val: f64, _dt: f64, _width: usize| {
            -(self.a * curr_val) / self.sigma
        };
        let sigma_prime = |_t_step: f64, _curr_val: f64, _dt: f64, _j: usize| 0.0;
        let sigma_inv = |_t_step: f64, y: f64, _dt: f64, _j: usize| self.sigma * y;
        let t_of_option = option_maturity - t;
        //The tree runs on a clock shifted by the valuation time `t` (it spans `option_maturity - t`),
        //but phi and every bond/swap leg are measured from "now" (0).  A node at tree time `tau` is
        //therefore at absolute time `t + tau`; using the shifted time directly misprices the swap
        //legs whenever `t > 0` (and whenever phi is not constant).
        let mut phi_cache: Vec<f64> = binomial_tree::get_all_t(t_of_option, num_steps)
            .map(|tau| self.phi_t(t + tau))
            .collect();
        phi_cache.push(self.phi_t(t + t_of_option));
        let payoff = |t_step: f64, curr_val: f64, _dt: f64, j: usize| {
            let t_abs = t + t_step;
            let swp = self.swap_price_t_init_raw(
                curr_val + phi_cache[j],
                t_abs,
                t_abs,
                num_swap_payments,
                delta,
                swap_rate,
            );
            payoff_swaption(is_payer, swp)
        };
        let discount = |_t_step: f64, curr_val: f64, dt: f64, j: usize| {
            (-(curr_val + phi_cache[j]) * dt).exp()
        };
        binomial_tree::compute_price_american(
            &alpha_div_sigma,
            &sigma_prime,
            &sigma_inv,
            &payoff,
            &discount,
            (r_t - self.phi_t(t)) / self.sigma, //initial "y"
            t_of_option,
            num_steps,
        )
    }
    /// Returns price of an American payer swaption at some future time t
    ///
    /// # Comments
    ///
    /// This function uses a tree to solve and will take longer to compute
    /// than other pricing functions.
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
    /// let num_tree_steps = 100;
    /// let swaption = hull_white.american_payer_swaption_t(
    ///     r_t, t, option_maturity, num_swap_payments, delta, swap_rate, num_tree_steps
    /// ).unwrap();
    /// ```
    pub fn american_payer_swaption_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64, //tenor of simple yield
        swap_rate: f64,
        num_steps: usize,
    ) -> Result<f64, HullWhiteError> {
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        validation::strictly_after("option_maturity", option_maturity, "t", t)?;
        validation::at_least_one("num_swap_payments", num_swap_payments)?;
        validation::positive("delta", delta)?;
        validation::finite("swap_rate", swap_rate)?;
        validation::at_least_one("num_steps", num_steps)?;
        Ok(self.american_swaption(
            r_t,
            t,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            true,
            num_steps,
        ))
    }
    /// Returns price of an American payer swaption at some future time t
    ///
    /// # Comments
    ///
    /// This function uses a tree to solve and will take longer to compute
    /// than other pricing functions.
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
    /// let num_tree_steps = 100;
    /// let swaption = hull_white.american_receiver_swaption_t(
    ///     r_t, t, option_maturity, num_swap_payments, delta, swap_rate, num_tree_steps
    /// ).unwrap();
    /// ```
    pub fn american_receiver_swaption_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64, //tenor of simple yield
        swap_rate: f64,
        num_steps: usize,
    ) -> Result<f64, HullWhiteError> {
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        validation::strictly_after("option_maturity", option_maturity, "t", t)?;
        validation::at_least_one("num_swap_payments", num_swap_payments)?;
        validation::positive("delta", delta)?;
        validation::finite("swap_rate", swap_rate)?;
        validation::at_least_one("num_steps", num_steps)?;
        Ok(self.american_swaption(
            r_t,
            t,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            false,
            num_steps,
        ))
    }
    #[cfg(test)]
    fn european_swaption_tree(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64, //tenor of simple yield
        swap_rate: f64,
        is_payer: bool,
        num_steps: usize,
    ) -> f64 {
        let alpha_div_sigma = |_t_step: f64, curr_val: f64, _dt: f64, _width: usize| {
            -(self.a * curr_val) / self.sigma
        };
        let sigma_prime = |_t_step: f64, _curr_val: f64, _dt: f64, _j: usize| 0.0;
        let sigma_inv = |_t_step: f64, y: f64, _dt: f64, _j: usize| self.sigma * y;
        let t_of_option = option_maturity - t;
        //Same time-coordinate convention as `american_swaption`: tree time `tau` means absolute
        //time `t + tau` for phi and for the swap legs.
        let mut phi_cache: Vec<f64> = binomial_tree::get_all_t(t_of_option, num_steps)
            .map(|tau| self.phi_t(t + tau))
            .collect();
        phi_cache.push(self.phi_t(t + t_of_option));
        let payoff = |t_step: f64, curr_val: f64, _dt: f64, j: usize| {
            let t_abs = t + t_step;
            let swp = self.swap_price_t_init_raw(
                curr_val + phi_cache[j],
                t_abs,
                t_abs,
                num_swap_payments,
                delta,
                swap_rate,
            );
            payoff_swaption(is_payer, swp)
        };
        let discount =
            |_t_a: f64, curr_val: f64, dt: f64, j: usize| (-(curr_val + phi_cache[j]) * dt).exp();
        binomial_tree::compute_price_raw(
            &alpha_div_sigma,
            &sigma_prime,
            &sigma_inv,
            &payoff,
            &discount,
            (r_t - self.phi_t(t)) / self.sigma, //initial "y"
            option_maturity - t,
            num_steps,
            false,
        )
    }
}

#[cfg(test)]
mod tests;
