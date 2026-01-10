//! The Hull-White interest rate model implementation

use crate::constants::{INITIAL_GUESS, MAX_ITERATIONS, PRECISION};
use crate::error::HullWhiteError;
use crate::utils;
use crate::Result;

/// The Hull-White interest rate model
///
/// This struct encapsulates the parameters and functions of the Hull-White model.
/// The model assumes that the short rate follows the stochastic differential equation:
/// dr(t) = [θ(t) - a*r(t)]dt + σ*dW(t)
///
/// where:
/// - a is the speed of mean reversion
/// - σ is the volatility parameter
/// - θ(t) is a time-dependent function calibrated to the initial term structure
///
/// # Example
///
/// ```
/// use hull_white::HullWhite;
///
/// // Define yield and forward curves
/// let yield_curve = |t: f64| 0.05 * t;  // Simple linear yield curve
/// let forward_curve = |t: f64| t.ln();  // Natural log forward curve
///
/// // Create a Hull-White model with parameters
/// let a = 0.1;      // Mean reversion speed
/// let sigma = 0.01; // Volatility parameter
/// let hull_white = HullWhite::new_panicking(a, sigma, &yield_curve, &forward_curve);
///
/// // Price a zero-coupon bond maturing in 2 years
/// let bond_price = hull_white.bond_price_now(2.0);
/// ```
pub struct HullWhite<'a, T, U>
where
    T: Fn(f64) -> f64 + Sync,
    U: Fn(f64) -> f64 + Sync,
{
    /// Speed of mean reversion parameter
    pub a: f64,
    /// Volatility parameter
    pub sigma: f64,
    /// Yield curve function (returns the "raw" yield, not divided by time)
    pub yield_curve: &'a T,
    /// Forward curve function
    pub forward_curve: &'a U,
}

impl<'a, T, U> HullWhite<'a, T, U>
where
    T: Fn(f64) -> f64 + Sync,
    U: Fn(f64) -> f64 + Sync,
{
    /// Create a new Hull-White model instance
    ///
    /// # Arguments
    /// * `a` - Speed of mean reversion parameter
    /// * `sigma` - Volatility parameter
    /// * `yield_curve` - Function representing the yield curve
    /// * `forward_curve` - Function representing the forward curve
    ///
    /// # Returns
    /// Result<HullWhite, HullWhiteError> indicating success or validation error
    pub fn new(a: f64, sigma: f64, yield_curve: &'a T, forward_curve: &'a U) -> Result<Self> {
        Self::validate_parameters(a, sigma)?;
        Ok(HullWhite {
            a,
            sigma,
            yield_curve,
            forward_curve,
        })
    }

    /// Create a new Hull-White model instance (panicking version)
    ///
    /// # Arguments
    /// * `a` - Speed of mean reversion parameter
    /// * `sigma` - Volatility parameter
    /// * `yield_curve` - Function representing the yield curve
    /// * `forward_curve` - Function representing the forward curve
    ///
    /// # Panics
    /// Panics if `a <= 0` or `sigma <= 0`
    pub fn new_panicking(a: f64, sigma: f64, yield_curve: &'a T, forward_curve: &'a U) -> Self {
        Self::validate_parameters(a, sigma).expect("Invalid Hull-White parameters");
        HullWhite {
            a,
            sigma,
            yield_curve,
            forward_curve,
        }
    }

    /// Validate Hull-White model parameters
    ///
    /// # Arguments
    /// * `a` - Speed of mean reversion parameter
    /// * `sigma` - Volatility parameter
    ///
    /// # Returns
    /// Result<(), HullWhiteError> indicating whether parameters are valid
    fn validate_parameters(a: f64, sigma: f64) -> Result<()> {
        if a <= 0.0 {
            return Err(HullWhiteError::InvalidInput(
                "Mean reversion parameter 'a' must be positive".to_string()
            ));
        }
        if sigma <= 0.0 {
            return Err(HullWhiteError::InvalidInput(
                "Volatility parameter 'sigma' must be positive".to_string()
            ));
        }
        if !a.is_finite() {
            return Err(HullWhiteError::InvalidInput(
                "Mean reversion parameter 'a' must be finite".to_string()
            ));
        }
        if !sigma.is_finite() {
            return Err(HullWhiteError::InvalidInput(
                "Volatility parameter 'sigma' must be finite".to_string()
            ));
        }
        Ok(())
    }


    /// Calculate the A(t) function used in the Hull-White model
    /// 
    /// # Arguments
    /// * `t_diff` - Time difference (T - t)
    fn a_t(&self, t_diff: f64) -> f64 {
        (1.0 - (-self.a * t_diff).exp()) / self.a
    }

    /// Calculate the A(t, T) function
    /// 
    /// # Arguments
    /// * `t` - Start time
    /// * `t_m` - End time (maturity)
    fn at_t(&self, t: f64, t_m: f64) -> f64 {
        self.a_t(t_m - t)
    }

    /// Calculate the C(t, T) function in the Hull-White model
    /// 
    /// # Arguments
    /// * `t` - Start time
    /// * `t_m` - End time (maturity)
    fn ct_t(&self, t: f64, t_m: f64) -> f64 {
        let sqr = (-self.a * t_m).exp() - (-self.a * t).exp();
        (self.yield_curve)(t) - (self.yield_curve)(t_m) + (self.forward_curve)(t) * self.at_t(t, t_m)
            - (self.sigma * sqr).powi(2) * ((2.0 * self.a * t).exp() - 1.0) / (4.0 * self.a.powi(3))
    }

    /// Calculate the φ(t) function in the Hull-White model
    /// 
    /// # Arguments
    /// * `t` - Time point
    fn phi_t(&self, t: f64) -> f64 {
        let exp_t = 1.0 - (-self.a * t).exp();
        (self.forward_curve)(t) + (self.sigma * exp_t).powi(2) / (2.0 * self.a.powi(2))
    }

    /// Calculate the volatility of a bond under the t-forward measure
    /// 
    /// # Arguments
    /// * `t` - Start time
    /// * `t_m` - Bond maturity
    /// * `t_f` - Forward measure maturity
    pub fn t_forward_bond_vol(&self, t: f64, t_m: f64, t_f: f64) -> f64 {
        let exp_d = 1.0 - (-self.a * (t_f - t_m)).exp();
        let exp_t = 1.0 - (-2.0 * self.a * (t_m - t)).exp();
        self.sigma * (exp_t / (2.0 * self.a.powi(3))).sqrt() * exp_d
    }

    /// Calculate the mean of the interest rate process
    /// 
    /// # Arguments
    /// * `r_t` - Current interest rate at time t
    /// * `t` - Start time
    /// * `t_m` - End time
    pub fn mu_r(&self, r_t: f64, t: f64, t_m: f64) -> f64 {
        self.phi_t(t_m) + (r_t - self.phi_t(t)) * (-self.a * (t_m - t)).exp()
    }

    /// Calculate the variance of the interest rate process
    /// 
    /// # Arguments
    /// * `t` - Start time
    /// * `t_m` - End time
    pub fn variance_r(&self, t: f64, t_m: f64) -> f64 {
        self.sigma.powi(2) * (1.0 - (-2.0 * self.a * (t_m - t)).exp()) / (2.0 * self.a)
    }

    /// Calculate the price of a zero-coupon bond at a future time t
    ///
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when bond is valued
    /// * `bond_maturity` - Maturity of the bond
    pub fn bond_price_t(&self, r_t: f64, t: f64, bond_maturity: f64) -> f64 {
        // Note: In a production environment, we might want to validate inputs here too
        // For performance reasons, we're skipping validation in the hot path
        (-r_t * self.at_t(t, bond_maturity) + self.ct_t(t, bond_maturity)).exp()
    }

    /// Calculate the derivative of bond price with respect to the interest rate
    fn bond_price_t_deriv(&self, r_t: f64, t: f64, bond_maturity: f64) -> f64 {
        let at_t_val = self.at_t(t, bond_maturity);
        -self.bond_price_t(r_t, t, bond_maturity) * at_t_val
    }

    /// Calculate the price of a zero-coupon bond at current time (t=0)
    ///
    /// # Arguments
    /// * `bond_maturity` - Maturity of the bond
    pub fn bond_price_now(&self, bond_maturity: f64) -> f64 {
        // Validate maturity
        if bond_maturity <= 0.0 || !bond_maturity.is_finite() {
            return 0.0; // or handle error appropriately
        }
        (-(self.yield_curve)(bond_maturity)).exp()
    }

    /// Calculate the price of a coupon bond at a future time t
    /// 
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when bond is valued
    /// * `coupon_times` - Vector of coupon payment times (including maturity)
    /// * `coupon_rate` - Rate of coupon payments
    pub fn coupon_bond_price_t(&self, r_t: f64, t: f64, coupon_times: &[f64], coupon_rate: f64) -> f64 {
        let par_value = 1.0; // Without loss of generality
        let last_index_coupon = coupon_times.len() - 1;
        
        coupon_times
            .iter()
            .enumerate()
            .map(|(index, coupon_time)| {
                let is_last = index == last_index_coupon;
                (coupon_rate + if is_last { par_value } else { 0.0 }) * self.bond_price_t(r_t, t, *coupon_time)
            })
            .sum()
    }

    /// Calculate the derivative of coupon bond price with respect to the interest rate
    fn coupon_bond_price_t_deriv(&self, r_t: f64, t: f64, coupon_times: &[f64], coupon_rate: f64) -> f64 {
        let par_value = 1.0; // Without loss of generality
        let last_index_coupon = coupon_times.len() - 1;
        
        coupon_times
            .iter()
            .enumerate()
            .map(|(index, coupon_time)| {
                let is_last = index == last_index_coupon;
                (coupon_rate + if is_last { par_value } else { 0.0 }) * self.bond_price_t_deriv(r_t, t, *coupon_time)
            })
            .sum()
    }

    /// Calculate the price of a coupon bond at current time (t=0)
    /// 
    /// # Arguments
    /// * `coupon_times` - Vector of coupon payment times (including maturity)
    /// * `coupon_rate` - Rate of coupon payments
    pub fn coupon_bond_price_now(&self, coupon_times: &[f64], coupon_rate: f64) -> f64 {
        let par_value = 1.0; // Without loss of generality
        let last_index_coupon = coupon_times.len() - 1;
        
        coupon_times
            .iter()
            .enumerate()
            .map(|(index, coupon_time)| {
                let is_last = index == last_index_coupon;
                (coupon_rate + if is_last { par_value } else { 0.0 }) * self.bond_price_now(*coupon_time)
            })
            .sum()
    }

    /// Calculate the price of a call option on a zero-coupon bond at a future time t
    /// 
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when option is valued
    /// * `option_maturity` - Maturity of the option
    /// * `bond_maturity` - Maturity of the underlying bond
    /// * `strike` - Strike price of the option
    pub fn bond_call_t(&self, r_t: f64, t: f64, option_maturity: f64, bond_maturity: f64, strike: f64) -> f64 {
        black_scholes::call_discount(
            self.bond_price_t(r_t, t, bond_maturity), // underlying
            strike,
            self.bond_price_t(r_t, t, option_maturity), // discount
            self.t_forward_bond_vol(t, option_maturity, bond_maturity), // volatility
        )
    }

    /// Calculate the price of a call option on a zero-coupon bond at current time (t=0)
    /// 
    /// # Arguments
    /// * `option_maturity` - Maturity of the option
    /// * `bond_maturity` - Maturity of the underlying bond
    /// * `strike` - Strike price of the option
    pub fn bond_call_now(&self, option_maturity: f64, bond_maturity: f64, strike: f64) -> f64 {
        black_scholes::call_discount(
            self.bond_price_now(bond_maturity), // underlying
            strike,
            self.bond_price_now(option_maturity), // discount
            self.t_forward_bond_vol(0.0, option_maturity, bond_maturity), // volatility
        )
    }

    /// Calculate the price of a put option on a zero-coupon bond at a future time t
    /// 
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when option is valued
    /// * `option_maturity` - Maturity of the option
    /// * `bond_maturity` - Maturity of the underlying bond
    /// * `strike` - Strike price of the option
    pub fn bond_put_t(&self, r_t: f64, t: f64, option_maturity: f64, bond_maturity: f64, strike: f64) -> f64 {
        black_scholes::put_discount(
            self.bond_price_t(r_t, t, bond_maturity), // underlying
            strike,
            self.bond_price_t(r_t, t, option_maturity), // discount
            self.t_forward_bond_vol(t, option_maturity, bond_maturity), // volatility
        )
    }

    /// Calculate the price of a put option on a zero-coupon bond at current time (t=0)
    /// 
    /// # Arguments
    /// * `option_maturity` - Maturity of the option
    /// * `bond_maturity` - Maturity of the underlying bond
    /// * `strike` - Strike price of the option
    pub fn bond_put_now(&self, option_maturity: f64, bond_maturity: f64, strike: f64) -> f64 {
        black_scholes::put_discount(
            self.bond_price_now(bond_maturity), // underlying
            strike,
            self.bond_price_now(option_maturity), // discount
            self.t_forward_bond_vol(0.0, option_maturity, bond_maturity), // volatility
        )
    }

    /// Calculate the price of a caplet at current time (t=0)
    /// 
    /// # Arguments
    /// * `option_maturity` - Maturity of the caplet
    /// * `delta` - Tenor of the LIBOR rate
    /// * `strike` - Strike rate of the caplet
    pub fn caplet_now(&self, option_maturity: f64, delta: f64, strike: f64) -> f64 {
        (strike * delta + 1.0)
            * self.bond_put_now(
                option_maturity,
                option_maturity + delta,
                1.0 / (delta * strike + 1.0),
            )
    }

    /// Calculate the price of a caplet at a future time t
    /// 
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when caplet is valued
    /// * `option_maturity` - Maturity of the caplet
    /// * `delta` - Tenor of the LIBOR rate
    /// * `strike` - Strike rate of the caplet
    pub fn caplet_t(&self, r_t: f64, t: f64, option_maturity: f64, delta: f64, strike: f64) -> f64 {
        (strike * delta + 1.0)
            * self.bond_put_t(
                r_t,
                t,
                option_maturity,
                option_maturity + delta,
                1.0 / (delta * strike + 1.0),
            )
    }

    /// Calculate the price of a Eurodollar futures contract at a future time t
    /// 
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when futures is valued
    /// * `option_maturity` - Maturity of the futures
    /// * `delta` - Tenor of the LIBOR rate
    pub fn euro_dollar_future_t(&self, r_t: f64, t: f64, option_maturity: f64, delta: f64) -> f64 {
        let gamma = self.gamma_edf(t, option_maturity, delta);
        self.edf_compute(
            self.bond_price_t(r_t, t, option_maturity),
            self.bond_price_t(r_t, t, option_maturity + delta),
            gamma,
            delta,
        )
    }

    /// Calculate the price of a Eurodollar futures contract at current time (t=0)
    /// 
    /// # Arguments
    /// * `option_maturity` - Maturity of the futures
    /// * `delta` - Tenor of the LIBOR rate
    pub fn euro_dollar_future_now(&self, option_maturity: f64, delta: f64) -> f64 {
        let gamma = self.gamma_edf(0.0, option_maturity, delta);
        self.edf_compute(
            self.bond_price_now(option_maturity),
            self.bond_price_now(option_maturity + delta),
            gamma,
            delta,
        )
    }

    /// Calculate the forward LIBOR rate at a future time t
    /// 
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when rate is calculated
    /// * `maturity` - Maturity of the forward rate
    /// * `delta` - Tenor of the LIBOR rate
    pub fn forward_libor_rate_t(&self, r_t: f64, t: f64, maturity: f64, delta: f64) -> f64 {
        let nearest_bond = self.bond_price_t(r_t, t, maturity);
        let farthest_bond = self.bond_price_t(r_t, t, maturity + delta);
        utils::compute_libor_rate(nearest_bond, farthest_bond, delta)
    }

    /// Calculate the forward LIBOR rate at current time (t=0)
    /// 
    /// # Arguments
    /// * `maturity` - Maturity of the forward rate
    /// * `delta` - Tenor of the LIBOR rate
    pub fn forward_libor_rate_now(&self, maturity: f64, delta: f64) -> f64 {
        let nearest_bond = self.bond_price_now(maturity);
        let farthest_bond = self.bond_price_now(maturity + delta);
        utils::compute_libor_rate(nearest_bond, farthest_bond, delta)
    }

    /// Calculate the LIBOR rate at a future time t
    /// 
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when rate is calculated
    /// * `delta` - Tenor of the LIBOR rate
    pub fn libor_rate_t(&self, r_t: f64, t: f64, delta: f64) -> f64 {
        self.forward_libor_rate_t(r_t, t, t, delta)
    }

    /// Calculate the forward swap rate at a future time t
    /// 
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when swap rate is calculated
    /// * `swap_initiation` - Time when the swap starts
    /// * `num_swap_payments` - Number of swap payments
    /// * `delta` - Tenor of the swap payments
    pub fn forward_swap_rate_t(
        &self,
        r_t: f64,
        t: f64,
        swap_initiation: f64, // must be greater than or equal to t
        num_swap_payments: usize,
        delta: f64,
    ) -> f64 {
        let denominator_swap: f64 = (1..=(num_swap_payments))
            .map(|curr| self.bond_price_t(r_t, t, swap_initiation + delta * (curr as f64)))
            .sum::<f64>()
            * delta;
        (self.bond_price_t(r_t, t, swap_initiation)
            - self.bond_price_t(r_t, t, swap_initiation + (num_swap_payments as f64) * delta))
            / denominator_swap
    }

    /// Calculate the swap rate at a future time t
    /// 
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when swap rate is calculated
    /// * `num_swap_payments` - Number of swap payments
    /// * `delta` - Tenor of the swap payments
    pub fn swap_rate_t(&self, r_t: f64, t: f64, num_swap_payments: usize, delta: f64) -> f64 {
        self.forward_swap_rate_t(
            r_t,
            t,
            t, // swap is a forward swap starting at time "0" (t-t=0)
            num_swap_payments,
            delta,
        )
    }

    /// Calculate the price of a swap at a future time t (not necessarily at initiation)
    /// 
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when swap is valued
    /// * `swap_maturity` - Maturity of the swap
    /// * `delta` - Tenor of the swap payments
    /// * `swap_rate` - Fixed rate of the swap
    pub fn swap_price_t(
        &self,
        r_t: f64,
        t: f64,
        swap_maturity: f64,
        delta: f64,
        swap_rate: f64,
    ) -> f64 {
        let (num_payments, is_exact) = utils::get_num_remaining_payments(t, swap_maturity, delta);
        let swap_start = if is_exact {
            t
        } else {
            swap_maturity - (num_payments as f64 - 1.0) * delta
        };
        self.swap_price_t_init(r_t, t, swap_start, num_payments, delta, swap_rate)
    }

    /// Calculate the price of a swap at the start of the swap
    /// 
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when swap is valued
    /// * `swap_start` - Time when the swap starts
    /// * `num_swap_payments` - Number of swap payments
    /// * `delta` - Tenor of the swap payments
    /// * `swap_rate` - Fixed rate of the swap
    pub fn swap_price_t_init(
        &self,
        r_t: f64,
        t: f64,
        swap_start: f64, // must be greater than or equal to t, and less than t+delta
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
    ) -> f64 {
        let sm_bond: f64 = (1..num_swap_payments)
            .map(|curr| {
                self.bond_price_t(r_t, t, swap_start + delta * (curr as f64)) * swap_rate * delta
            })
            .sum();
        self.bond_price_t(r_t, t, swap_start)
            - sm_bond
            - (1.0 + swap_rate * delta)
                * self.bond_price_t(r_t, t, swap_start + delta * (num_swap_payments as f64))
    }

    /// Calculate the price of a European payer swaption at a future time t
    /// 
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when swaption is valued
    /// * `option_maturity` - Maturity of the swaption
    /// * `num_swap_payments` - Number of swap payments
    /// * `delta` - Tenor of the swap payments
    /// * `swap_rate` - Strike rate of the swaption
    pub fn european_payer_swaption_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
    ) -> Result<f64> {
        let coupon_times = utils::get_coupon_times(num_swap_payments, option_maturity, delta);
        let strike = 1.0;
        self.coupon_bond_option_generic_t(
            r_t,
            t,
            option_maturity,
            &coupon_times,
            swap_rate * delta,
            strike,
            false, // is_put = true for payer swaption
        )
    }

    /// Calculate the price of a European receiver swaption at a future time t
    /// 
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when swaption is valued
    /// * `option_maturity` - Maturity of the swaption
    /// * `num_swap_payments` - Number of swap payments
    /// * `delta` - Tenor of the swap payments
    /// * `swap_rate` - Strike rate of the swaption
    pub fn european_receiver_swaption_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
    ) -> Result<f64> {
        let coupon_times = utils::get_coupon_times(num_swap_payments, option_maturity, delta);
        let strike = 1.0;
        self.coupon_bond_option_generic_t(
            r_t,
            t,
            option_maturity,
            &coupon_times,
            swap_rate * delta,
            strike,
            true, // is_call = true for receiver swaption
        )
    }

    /// Generic function for pricing options on coupon bonds using Jamshidian's trick
    fn coupon_bond_option_generic_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
        is_put: bool,
    ) -> Result<f64> {
        let par_value = 1.0;
        let final_coupon_index = coupon_times.len() - 1;

        let fn_to_optimize = |r| {
            self.coupon_bond_price_t(r, option_maturity, coupon_times, coupon_rate) - strike
        };
        let fn_deriv = |r| {
            self.coupon_bond_price_t_deriv(r, option_maturity, coupon_times, coupon_rate)
        };

        let r_optimal = nrfind::find_root(&fn_to_optimize, &fn_deriv, INITIAL_GUESS, PRECISION, MAX_ITERATIONS)
            .map_err(|_| HullWhiteError::RootFindingError("Failed to find optimal rate".to_string()))?;

        let result = coupon_times
            .iter()
            .enumerate()
            .map(|(index, coupon_time)| {
                let is_last = final_coupon_index == index;
                let bond_price_at_optimal = self.bond_price_t(r_optimal, option_maturity, *coupon_time);

                let option_price = if !is_put {  // is_call = !is_put
                    self.bond_call_t(r_t, t, option_maturity, *coupon_time, bond_price_at_optimal)
                } else {
                    self.bond_put_t(r_t, t, option_maturity, *coupon_time, bond_price_at_optimal)
                };

                option_price * (coupon_rate + if is_last { par_value } else { 0.0 })
            })
            .sum();

        Ok(result)
    }

    /// Calculate the price of a call option on a coupon bond at a future time t
    ///
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when option is valued
    /// * `option_maturity` - Maturity of the option
    /// * `coupon_times` - Vector of coupon payment times (including maturity)
    /// * `coupon_rate` - Rate of coupon payments
    /// * `strike` - Strike price of the option
    pub fn coupon_bond_call_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
    ) -> Result<f64> {
        self.coupon_bond_option_generic_t(
            r_t,
            t,
            option_maturity,
            coupon_times,
            coupon_rate,
            strike,
            false, // is_put = false means it's a call
        )
    }

    /// Calculate the price of a put option on a coupon bond at a future time t
    ///
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when option is valued
    /// * `option_maturity` - Maturity of the option
    /// * `coupon_times` - Vector of coupon payment times (including maturity)
    /// * `coupon_rate` - Rate of coupon payments
    /// * `strike` - Strike price of the option
    pub fn coupon_bond_put_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
    ) -> Result<f64> {
        self.coupon_bond_option_generic_t(
            r_t,
            t,
            option_maturity,
            coupon_times,
            coupon_rate,
            strike,
            true, // is_put = true means it's a put
        )
    }

    /// Calculate the gamma for Eurodollar futures
    fn gamma_edf(&self, t: f64, option_maturity: f64, delta: f64) -> f64 {
        let exp_t = (-self.a * (option_maturity - t)).exp();
        let exp_d = (-self.a * delta).exp();
        (self.sigma.powi(2) / self.a.powi(3))
            * (1.0 - exp_d)
            * ((1.0 - exp_t) - exp_d * 0.5 * (1.0 - exp_t.powi(2)))
    }

    /// Calculate the Eurodollar futures price
    fn edf_compute(&self, bond_num: f64, bond_den: f64, gamma: f64, delta: f64) -> f64 {
        ((bond_num / bond_den) * gamma.exp() - 1.0) / delta
    }

    /// Calculate the price of an American payer swaption at a future time t
    ///
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when swaption is valued
    /// * `option_maturity` - Maturity of the swaption
    /// * `num_swap_payments` - Number of swap payments
    /// * `delta` - Tenor of the swap payments
    /// * `swap_rate` - Strike rate of the swaption
    /// * `num_steps` - Number of steps in the binomial tree
    pub fn american_payer_swaption_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64, // tenor of simple yield
        swap_rate: f64,
        num_steps: usize,
    ) -> f64 {
        self.american_swaption(
            r_t,
            t,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            true,  // is_payer
            num_steps,
        )
    }

    /// Calculate the price of an American receiver swaption at a future time t
    ///
    /// # Arguments
    /// * `r_t` - Interest rate at time t
    /// * `t` - Time when swaption is valued
    /// * `option_maturity` - Maturity of the swaption
    /// * `num_swap_payments` - Number of swap payments
    /// * `delta` - Tenor of the swap payments
    /// * `swap_rate` - Strike rate of the swaption
    /// * `num_steps` - Number of steps in the binomial tree
    pub fn american_receiver_swaption_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64, // tenor of simple yield
        swap_rate: f64,
        num_steps: usize,
    ) -> f64 {
        self.american_swaption(
            r_t,
            t,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            false, // is_payer = false means receiver
            num_steps,
        )
    }

    /// Internal function for American swaption pricing using binomial trees
    fn american_swaption(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64, // tenor of simple yield
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
        let mut phi_cache: Vec<f64> = binomial_tree::get_all_t(t_of_option, num_steps)
            .map(|t_a| self.phi_t(t_a))
            .collect();
        phi_cache.push(self.phi_t(t_of_option));
        let payoff = |t_step: f64, curr_val: f64, _dt: f64, j: usize| {
            let swp = self.swap_price_t_init(
                curr_val + phi_cache[j],
                t_step,
                t_step,
                num_swap_payments,
                delta,
                swap_rate,
            );
            utils::payoff_swaption(is_payer, swp)
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
            (r_t - self.phi_t(t)) / self.sigma, // initial "y"
            t_of_option,
            num_steps,
        )
    }
}