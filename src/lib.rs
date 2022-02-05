//! This is a library of fixed income pricers using a
//! Hull White process as the underlying process.  The
//! fundamental times here are (0, t, T, TM).  0 is the
//! current time (and is  reflective of the current
//! yield  curve) while t is some future time  that we may
//! want to price options at given the underlying at
//! that time.  T and TM are shorthands for a variety of
//! asset times.  For example, an option on a bond requires
//! an option maturity and a bond maturity.  The option
//! maturity should be before the bond maturity, but after
//! the future time t.  Note that ALL TIMES ARE WITH
//! RESPECT TO 0!
//!

use core::num;

const PREC_1: f64 = 0.0000001;
const R_INIT: f64 = 0.03;
const MAX_ITER: i32 = 50;

//tdiff=T-t
fn a_t(a: f64, t_diff: f64) -> f64 {
    (1.0 - (-a * t_diff).exp()) / a
}

//t is first future time
//t_m is second future time
fn at_t(a: f64, t: f64, t_m: f64) -> f64 {
    a_t(a, t_m - t)
}

fn ct_t(
    a: f64,
    sigma: f64,
    t: f64,
    t_m: f64,
    yield_curve: &dyn Fn(f64) -> f64,
    forward_curve: &dyn Fn(f64) -> f64,
) -> f64 {
    let sqr = (-a * t_m).exp() - (-a * t).exp();
    yield_curve(t) - yield_curve(t_m) + forward_curve(t) * at_t(a, t, t_m)
        - (sigma * sqr).powi(2) * ((2.0 * a * t).exp() - 1.0) / (4.0 * a.powi(3))
}

//https://www.math.nyu.edu/~alberts/spring07/Lecture5.pdf
//https://developers.opengamma.com/quantitative-research/Hull-White-One-Factor-Model-OpenGamma.pdf (note that in the open gamma derivation, t0=option_maturity)
fn gamma_edf(a: f64, sigma: f64, t: f64, option_maturity: f64, delta: f64) -> f64 {
    let exp_t = (-a * (option_maturity - t)).exp();
    let exp_d = (-a * delta).exp();
    (sigma.powi(2) / a.powi(3))
        * (1.0 - exp_d)
        * ((1.0 - exp_t) - exp_d * 0.5 * (1.0 - exp_t.powi(2)))
}
fn edf_compute(bond_num: f64, bond_den: f64, gamma: f64, delta: f64) -> f64 {
    ((bond_num / bond_den) * gamma.exp() - 1.0) / delta
}

fn compute_libor_rate(nearest_bond: f64, farthest_bond: f64, tenor: f64) -> f64 {
    (nearest_bond - farthest_bond) / (farthest_bond * tenor)
}

//payments should be an integer!
fn get_num_payments(t: f64, maturity: f64, delta: f64) -> Result<usize, f64> {
    let num_payments = (maturity - t) / delta + 1.0;
    if num_payments.trunc() == num_payments {
        Ok(num_payments as usize)
    } else {
        Err(num_payments)
    }
}

pub fn get_coupon_times(num_payments: usize, t: f64, delta: f64) -> Vec<f64> {
    (1..(num_payments + 1))
        .map(|index| get_time_from_t_index(index, t, delta))
        .collect()
}

fn get_time_from_t_index(index: usize, t: f64, delta: f64) -> f64 {
    t + (index as f64) * delta
}

fn max_or_zero(v: f64) -> f64 {
    if v > 0.0 {
        v
    } else {
        0.0
    }
}
fn payoff_swaption(is_payer: bool, swp: f64) -> f64 {
    match is_payer {
        true => max_or_zero(swp),
        false => max_or_zero(-swp),
    }
}

fn coupon_bond_generic_t(
    r_t: f64,
    t: f64,
    coupon_times: &[f64], //includes bond_maturity
    coupon_rate: f64,
    generic_fn: &impl Fn(f64, f64, f64) -> f64,
) -> f64 {
    let par_value = 1.0; //without loss of generality
    let last_index_coupon = coupon_times.len() - 1;
    coupon_times
        .iter()
        .enumerate()
        .map(|(index, coupon_time)| {
            let is_last = index == last_index_coupon;
            (coupon_rate + if is_last { par_value } else { 0.0 }) * generic_fn(r_t, t, *coupon_time)
        })
        .sum()
}
fn coupon_bond_generic_now(
    coupon_times: &[f64], //includes bond_maturity
    coupon_rate: f64,
    generic_fn: &impl Fn(f64) -> f64,
) -> f64 {
    let par_value = 1.0; //without loss of generality
    let last_index_coupon = coupon_times.len() - 1;
    coupon_times
        .iter()
        .enumerate()
        .map(|(index, coupon_time)| {
            let is_last = index == last_index_coupon;
            (coupon_rate + if is_last { par_value } else { 0.0 }) * generic_fn(*coupon_time)
        })
        .sum()
}
pub struct HullWhite<'a> {
    a: f64,
    sigma: f64,
    yield_curve: &'a dyn Fn(f64) -> f64, //this is not divided by time, so this gets perpetually larger (unless rates are negative)
    forward_curve: &'a dyn Fn(f64) -> f64,
}

impl HullWhite<'_> {
    pub fn init<'a>(
        a: f64,
        sigma: f64,
        yield_curve: &'a dyn Fn(f64) -> f64,
        forward_curve: &'a dyn Fn(f64) -> f64,
    ) -> HullWhite<'a> {
        HullWhite {
            a,
            sigma,
            yield_curve,
            forward_curve,
        }
    }
    /// Returns volality of bond under the t-forward measure.
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0;
    /// let t_m = 2.0;
    /// let t_f = 3.0;
    /// let yield_curve = |t:f64|0.05*t;
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white= hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let bond_vol = hull_white.t_forward_bond_vol(
    ///     t, t_m, t_f
    /// );
    /// ```
    pub fn t_forward_bond_vol(&self, t: f64, t_m: f64, t_f: f64) -> f64 {
        let exp_d = 1.0 - (-self.a * (t_f - t_m)).exp();
        let exp_t = 1.0 - (-2.0 * self.a * (t_m - t)).exp();
        self.sigma * (exp_t / (2.0 * self.a.powi(3))).sqrt() * exp_d
    }
    fn phi_t(&self, t: f64) -> f64 {
        let exp_t = 1.0 - (-self.a * t).exp();
        (self.forward_curve)(t) + (self.sigma * exp_t).powi(2) / (2.0 * self.a.powi(2))
    }
    /// Returns volality of bond under the t-forward measure.
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start taking the expectation
    /// let t_m = 2.0; //horizon of the expectation
    /// let r_t = 0.04; //rate at t
    /// let yield_curve = |t:f64|0.05*t;
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white= hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let bond_vol = hull_white.mu_r(r_t, t, t_m);
    /// ```
    pub fn mu_r(&self, r_t: f64, t: f64, t_m: f64) -> f64 {
        self.phi_t(t_m) + (r_t - self.phi_t(t)) * (-self.a * (t_m - t)).exp()
    }
    /// Returns variance of the interest rate process
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start taking the variance
    /// let t_m = 2.0; //horizon of the variance
    /// let yield_curve = |t:f64|0.05*t;
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let variance = hull_white.variance_r(t, t_m);
    /// ```
    pub fn variance_r(&self, t: f64, t_m: f64) -> f64 {
        self.sigma.powi(2) * (1.0 - (-2.0 * self.a * (t_m - t)).exp()) / (2.0 * self.a)
    }
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let bond_price = hull_white.bond_price_t(r_t, t, bond_maturity);
    /// ```
    pub fn bond_price_t(&self, r_t: f64, t: f64, bond_maturity: f64) -> f64 {
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
    fn bond_price_t_deriv(&self, r_t: f64, t: f64, bond_maturity: f64) -> f64 {
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let bond_price = hull_white.bond_price_now(bond_maturity);
    /// ```
    pub fn bond_price_now(&self, bond_maturity: f64) -> f64 {
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let bond_price = hull_white.coupon_bond_price_t(r_t, t, &coupon_times, coupon_rate);
    /// ```
    pub fn coupon_bond_price_t(
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
            &|r_t: f64, t: f64, bond_maturity: f64| self.bond_price_t(r_t, t, bond_maturity),
        )
    }
    fn coupon_bond_price_t_deriv(
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let bond_price = hull_white.coupon_bond_price_now(&coupon_times, coupon_rate);
    /// ```
    pub fn coupon_bond_price_now(
        &self,
        coupon_times: &[f64], //does not include the bond_maturity, but the function does check for that
        coupon_rate: f64,
    ) -> f64 {
        coupon_bond_generic_now(coupon_times, coupon_rate, &|bond_maturity: f64| {
            self.bond_price_now(bond_maturity)
        })
    }
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let bond_call = hull_white.bond_call_t(r_t, t, option_maturity, bond_maturity, strike);
    /// ```
    pub fn bond_call_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        bond_maturity: f64,
        strike: f64,
    ) -> f64 {
        black_scholes::call_discount(
            self.bond_price_t(r_t, t, bond_maturity), //underlying
            strike,
            self.bond_price_t(r_t, t, option_maturity), //discount
            self.t_forward_bond_vol(t, option_maturity, bond_maturity), //volatility with maturity
        )
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let bond_call = hull_white.bond_call_now(option_maturity, bond_maturity, strike);
    /// ```
    pub fn bond_call_now(&self, option_maturity: f64, bond_maturity: f64, strike: f64) -> f64 {
        let t = 0.0; //since "now"
        black_scholes::call_discount(
            self.bond_price_now(bond_maturity), //underlying
            strike,
            self.bond_price_now(option_maturity), //discount
            self.t_forward_bond_vol(t, option_maturity, bond_maturity), //volatility with maturity
        )
    }
    //The price of a call option on coupon bond under Hull White...uses jamshidian's trick*
    fn coupon_bond_option_generic_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
        generic_fn: &impl Fn(f64, f64, f64, f64, f64) -> f64,
    ) -> Result<f64, f64> {
        let par_value = 1.0;
        let final_coupon_index = coupon_times.len() - 1;
        let fn_to_optimize =
            |r| self.coupon_bond_price_t(r, option_maturity, coupon_times, coupon_rate) - strike;
        let fn_derv =
            |r| self.coupon_bond_price_t_deriv(r, option_maturity, coupon_times, coupon_rate);
        let r_optimal = nrfind::find_root(&fn_to_optimize, &fn_derv, R_INIT, PREC_1, MAX_ITER)?;
        Ok(coupon_times
            .iter()
            .enumerate()
            .map(|(index, coupon_time)| {
                let is_last = final_coupon_index == index;
                generic_fn(
                    r_t,
                    t,
                    option_maturity,
                    *coupon_time,
                    self.bond_price_t(r_optimal, option_maturity, *coupon_time),
                ) * (coupon_rate + if is_last { par_value } else { 0.0 })
            })
            .sum())
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
    /// let coupon_times = vec![1.25, 1.5, 1.75, 2.0, 2.5, 3.0];
    /// let coupon_rate = 0.05;
    /// let strike = 1.0;
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let bond_call = hull_white.coupon_bond_call_t(r_t, t, option_maturity, &coupon_times, coupon_rate, strike);
    /// ```
    pub fn coupon_bond_call_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
    ) -> Result<f64, f64> {
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let bond_put = hull_white.bond_put_t(r_t, t, option_maturity, bond_maturity, strike);
    /// ```
    pub fn bond_put_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        bond_maturity: f64,
        strike: f64,
    ) -> f64 {
        black_scholes::put_discount(
            self.bond_price_t(r_t, t, bond_maturity), //underlying
            strike,
            self.bond_price_t(r_t, t, option_maturity), //discount
            self.t_forward_bond_vol(t, option_maturity, bond_maturity), //volatility with maturity
        )
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let bond_put = hull_white.bond_put_now(option_maturity, bond_maturity, strike);
    /// ```
    pub fn bond_put_now(&self, option_maturity: f64, bond_maturity: f64, strike: f64) -> f64 {
        let t = 0.0; //since "now"
        black_scholes::put_discount(
            self.bond_price_now(bond_maturity), //underlying
            strike,
            self.bond_price_now(option_maturity), //discount
            self.t_forward_bond_vol(t, option_maturity, bond_maturity), //volatility with maturity
        )
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
    /// let coupon_times = vec![1.25, 1.5, 1.75, 2.0, 2.5, 3.0];
    /// let coupon_rate = 0.05;
    /// let strike = 1.0;
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let bond_put = hull_white.coupon_bond_put_t(r_t, t, option_maturity, &coupon_times, coupon_rate, strike);
    /// ```
    pub fn coupon_bond_put_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
    ) -> Result<f64, f64> {
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let caplet = hull_white.caplet_now(option_maturity, delta, strike);
    /// ```
    pub fn caplet_now(&self, option_maturity: f64, delta: f64, strike: f64) -> f64 {
        (strike * delta + 1.0)
            * self.bond_put_now(
                option_maturity,
                option_maturity + delta,
                1.0 / (delta * strike + 1.0),
            )
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let caplet = hull_white.caplet_t(r_t, t, option_maturity, delta, strike);
    /// ```
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let edf = hull_white.euro_dollar_future_t(r_t,  t, option_maturity, delta);
    /// ```
    pub fn euro_dollar_future_t(&self, r_t: f64, t: f64, option_maturity: f64, delta: f64) -> f64 {
        let gamma = gamma_edf(self.a, self.sigma, t, option_maturity, delta);
        edf_compute(
            self.bond_price_t(r_t, t, option_maturity),
            self.bond_price_t(r_t, t, option_maturity + delta),
            gamma,
            delta,
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let edf = hull_white.euro_dollar_future_now(option_maturity, delta);
    /// ```
    pub fn euro_dollar_future_now(&self, option_maturity: f64, delta: f64) -> f64 {
        let gamma = gamma_edf(self.a, self.sigma, 0.0, option_maturity, delta);
        edf_compute(
            self.bond_price_now(option_maturity),
            self.bond_price_now(option_maturity + delta),
            gamma,
            delta,
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let forward_libor = hull_white.forward_libor_rate_t(r_t, t, maturity, delta);
    /// ```
    pub fn forward_libor_rate_t(&self, r_t: f64, t: f64, maturity: f64, delta: f64) -> f64 {
        let nearest_bond = self.bond_price_t(r_t, t, maturity);
        let farthest_bond = self.bond_price_t(r_t, t, maturity + delta);
        compute_libor_rate(nearest_bond, farthest_bond, delta)
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let forward_libor = hull_white.forward_libor_rate_now(maturity, delta);
    /// ```
    pub fn forward_libor_rate_now(&self, maturity: f64, delta: f64) -> f64 {
        let nearest_bond = self.bond_price_now(maturity);
        let farthest_bond = self.bond_price_now(maturity + delta);
        compute_libor_rate(nearest_bond, farthest_bond, delta)
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let libor = hull_white.libor_rate_t(r_t, t, delta);
    /// ```
    pub fn libor_rate_t(&self, r_t: f64, t: f64, delta: f64) -> f64 {
        self.forward_libor_rate_t(r_t, t, t, delta)
    }

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
    /// let swap_maturity = 5.0;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let forward_swap = hull_white.forward_swap_rate_t(r_t,  t, swap_initiation, swap_maturity, delta).unwrap();
    /// ```
    pub fn forward_swap_rate_t(
        &self,
        r_t: f64,
        t: f64,
        swap_initiation: f64,
        swap_maturity: f64,
        delta: f64,
    ) -> Result<f64, f64> {
        let num_payments = get_num_payments(swap_initiation, swap_maturity, delta)?; //this should be an integer!  remember, T-t is the total swap length
        let denominator_swap: f64 = (1..(num_payments + 1))
            .map(|curr| self.bond_price_t(r_t, t, swap_initiation + delta * (curr as f64)))
            .sum::<f64>()
            * delta;
        Ok((self.bond_price_t(r_t, t, swap_initiation)
            - self.bond_price_t(r_t, t, swap_maturity + delta))
            / denominator_swap)
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
    /// let swap_maturity = 5.0;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let swap_rate = hull_white.swap_rate_t(r_t, t, swap_maturity, delta).unwrap();
    /// ```
    pub fn swap_rate_t(
        &self,
        r_t: f64,
        t: f64,
        swap_maturity: f64,
        delta: f64,
    ) -> Result<f64, f64> {
        self.forward_swap_rate_t(
            r_t,
            t,
            t, //swap is a forward swap starting at time "0" (t-t=0)
            swap_maturity,
            delta,
        )
    }
    /// Returns price of a swap at some future time
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let swap = hull_white.swap_price_t(r_t, t, swap_maturity, delta, swap_rate).unwrap();
    /// ```
    pub fn swap_price_t(
        &self,
        r_t: f64,
        t: f64,
        swap_maturity: f64,
        delta: f64,
        swap_rate: f64,
    ) -> Result<f64, f64> {
        let num_payments = get_num_payments(t, swap_maturity, delta)?; //this should be an integer!  remember, T-t is the total swap length
        let first_exchange_date = swap_maturity - (num_payments as f64 - 1.0) * delta;
        let sm_bond: f64 = (1..num_payments)
            .map(|curr| {
                self.bond_price_t(r_t, t, first_exchange_date + delta * (curr as f64))
                    * swap_rate
                    * delta
            })
            .sum();
        Ok(self.bond_price_t(r_t, t, first_exchange_date)
            - sm_bond
            - (1.0 + swap_rate * delta)
                * self.bond_price_t(r_t, t, first_exchange_date + delta * (num_payments as f64)))
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
    /// let swap_tenor = 4.0; //swap_tenor is how long the swap will be once entered in
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let swap_rate = 0.04; //the swap rate is what the payer agrees to pay if option is exercised
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let swaption = hull_white.european_payer_swaption_t(r_t, t, swap_tenor, option_maturity, delta, swap_rate).unwrap();
    /// ```
    pub fn european_payer_swaption_t(
        &self,
        r_t: f64,
        t: f64,
        swap_tenor: f64,
        option_maturity: f64,
        delta: f64,
        swap_rate: f64,
    ) -> Result<f64, f64> {
        let num_payments = get_num_payments(t, swap_tenor, delta)?;
        let coupon_times = get_coupon_times(num_payments, option_maturity, delta);
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
    /// let swap_tenor = 4.0; //swap_tenor is how long the swap will be once entered in
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let swap_rate = 0.04; //the swap rate is what the payer agrees to pay if option is exercised
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let swaption = hull_white.european_receiver_swaption_t(r_t, t, swap_tenor, option_maturity, delta, swap_rate).unwrap();
    /// ```
    pub fn european_receiver_swaption_t(
        &self,
        r_t: f64,
        t: f64,
        swap_tenor: f64,
        option_maturity: f64,
        delta: f64,
        swap_rate: f64,
    ) -> Result<f64, f64> {
        let num_payments = get_num_payments(t, swap_tenor, delta)?;
        let coupon_times = get_coupon_times(num_payments, option_maturity, delta);
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

    fn american_swaption(
        &self,
        r_t: f64,
        t: f64,
        swap_tenor: f64,
        option_maturity: f64,
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
        let mut phi_cache: Vec<f64> = binomial_tree::get_all_t(t_of_option, num_steps)
            .map(|t_a| self.phi_t(t_a))
            .collect();
        phi_cache.push(self.phi_t(t_of_option));
        let payoff = |t_step: f64, curr_val: f64, _dt: f64, j: usize| {
            let swp = self
                .swap_price_t(
                    curr_val + phi_cache[j],
                    t_step,
                    swap_tenor + t_step,
                    delta,
                    swap_rate,
                )
                .expect("must be integer number periods"); //ugh
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
    /// let swap_tenor = 4.0; //swap_tenor is how long the swap will be once entered in
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let swap_rate = 0.04; //the swap rate is what the payer agrees to pay if option is exercised
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let num_tree_steps = 100;
    /// let swaption = hull_white.american_payer_swaption_t(
    ///     r_t, t, swap_tenor, option_maturity, delta, swap_rate, num_tree_steps
    /// );
    /// ```
    pub fn american_payer_swaption_t(
        &self,
        r_t: f64,
        t: f64,
        swap_tenor: f64,
        option_maturity: f64,
        delta: f64, //tenor of simple yield
        swap_rate: f64,
        num_steps: usize,
    ) -> f64 {
        self.american_swaption(
            r_t,
            t,
            swap_tenor,
            option_maturity,
            delta,
            swap_rate,
            true,
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
    /// let swap_tenor = 4.0; //swap_tenor is how long the swap will be once entered in
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let swap_rate = 0.04; //the swap rate is what the payer agrees to pay if option is exercised
    /// let yield_curve = |t:f64|0.05*t; //yield curve returns the "raw" yield (not divided by maturity)
    /// let forward_curve = |t:f64|t.ln();
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve);
    /// let num_tree_steps = 100;
    /// let swaption = hull_white.american_receiver_swaption_t(
    ///     r_t, t, swap_tenor, option_maturity, delta, swap_rate, num_tree_steps
    /// );
    /// ```
    pub fn american_receiver_swaption_t(
        &self,
        r_t: f64,
        t: f64,
        swap_tenor: f64,
        option_maturity: f64,
        delta: f64, //tenor of simple yield
        swap_rate: f64,
        num_steps: usize,
    ) -> f64 {
        self.american_swaption(
            r_t,
            t,
            swap_tenor,
            option_maturity,
            delta,
            swap_rate,
            false,
            num_steps,
        )
    }
    #[cfg(test)]
    fn european_swaption_tree(
        &self,
        r_t: f64,
        t: f64,
        swap_tenor: f64,
        option_maturity: f64,
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
        let mut phi_cache: Vec<f64> = binomial_tree::get_all_t(option_maturity - t, num_steps)
            .map(|t_a| self.phi_t(t_a))
            .collect();
        phi_cache.push(self.phi_t(option_maturity - t));
        let payoff = |t_step: f64, curr_val: f64, _dt: f64, j: usize| {
            let swp = self
                .swap_price_t(
                    curr_val + phi_cache[j],
                    t_step,
                    swap_tenor + t_step,
                    delta,
                    swap_rate,
                )
                .unwrap();
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
mod tests {
    use super::*;
    use approx::*;
    use rand::distributions::{Distribution, StandardNormal};
    use rand::{SeedableRng, StdRng};
    fn get_rng_seed(seed: [u8; 32]) -> StdRng {
        SeedableRng::from_seed(seed)
    }
    #[test]
    fn test_get_num_payments() {
        let t = 0.5;
        let maturity = 2.0;
        let delta = 0.25;
        let num_payments = get_num_payments(t, maturity, delta).unwrap();
        assert_eq!(num_payments, 7);
    }

    #[test]
    fn test_get_num_payments_not_even() {
        let t = 0.5;
        let maturity = 2.0;
        let delta = 0.4;
        assert_eq!(get_num_payments(t, maturity, delta).is_err(), true);
    }

    #[test]
    fn test_get_time_from_t_index() {
        let t = 0.5;
        let delta = 0.4;
        let index = 3;
        let time = get_time_from_t_index(index, t, delta);
        assert_abs_diff_eq!(time, 1.7, epsilon = 0.0000001);
    }
    #[test]
    fn test_get_time_from_t_index_with_zero_index() {
        let t = 0.5;
        let delta = 0.4;
        let index = 0;
        let time = get_time_from_t_index(index, t, delta);
        assert_eq!(time, t);
    }
    #[test]
    fn test_get_coupon_times() {
        let num_payments = 5;
        let t = 1.0;
        let delta = 0.25;
        let coupon_times = get_coupon_times(num_payments, t, delta);
        let expected_coupon_times = vec![1.25, 1.5, 1.75, 2.0, 2.25];
        coupon_times
            .iter()
            .zip(expected_coupon_times.iter())
            .for_each(|(actual, expected)| assert_eq!(actual, expected))
    }
    #[test]
    fn test_get_coupon_times_no_payments() {
        let num_payments = 0;
        let t = 1.0;
        let delta = 0.25;
        let coupon_times = get_coupon_times(num_payments, t, delta);
        assert_eq!(coupon_times.len(), 0);
    }

    #[test]
    fn test_max_or_zero() {
        let v = 1.0;
        assert_eq!(max_or_zero(v), 1.0);
        assert_eq!(max_or_zero(-v), 0.0);
    }

    #[test]
    fn test_payoff_swaption() {
        let v = 1.0;
        assert_eq!(payoff_swaption(true, v), 1.0);
        assert_eq!(payoff_swaption(false, v), 0.0);
        assert_eq!(payoff_swaption(true, -v), 0.0);
        assert_eq!(payoff_swaption(false, -v), 1.0);
    }
    #[test]
    fn compare_caplet() {
        let curr_rate = 0.02;
        let sig: f64 = 0.02;
        let a: f64 = 0.3;
        let b = 0.04;
        let future_time = 0.0;
        let option_maturity = 1.5;
        let strike = 0.02;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct =
                (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let delta = 0.25;
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve);
        let caplet_n = hull_white.caplet_now(option_maturity, delta, strike);
        let caplet = hull_white.caplet_t(curr_rate, future_time, option_maturity, delta, strike);
        assert_abs_diff_eq!(caplet_n, caplet, epsilon = 0.00001);
    }
    #[test]
    fn compare_libor() {
        let curr_rate = 0.02;
        let sig: f64 = 0.02;
        let a: f64 = 0.3;
        let b = 0.04;
        let future_time = 0.0;
        let maturity = 1.5;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct =
                (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let delta = 0.25;
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve);
        let libor_n = hull_white.forward_libor_rate_now(maturity, delta);
        let libor_t = hull_white.forward_libor_rate_t(curr_rate, future_time, maturity, delta);
        assert_abs_diff_eq!(libor_n, libor_t, epsilon = 0.0001);
    }
    #[test]
    fn test_caplet() {
        let curr_rate = 0.02;
        let sig: f64 = 0.02;
        let a: f64 = 0.3;
        let b = 0.04;
        let future_time = 0.0;
        let option_maturity = 1.5;
        let strike = 0.02;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct =
                (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let delta = 0.25;
        let seed: [u8; 32] = [2; 32];
        let mut rng_seed = get_rng_seed(seed);
        let normal = StandardNormal;
        let num_sims: usize = 1000; //hopefully accurate
        let num_discrete_steps: usize = 1000;
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve);
        let total_sum = (0..num_sims).fold(0.0, move |accum, _sample_index| {
            let mut sum_r = 0.0;
            let mut running_r = curr_rate;
            let dt = (option_maturity - future_time) / (num_discrete_steps as f64 - 1.0);
            (0..num_discrete_steps).for_each(|t_index| {
                let norm = normal.sample(&mut rng_seed);
                let curr_t = dt * (t_index as f64) + future_time;
                let curr_vol = hull_white.variance_r(curr_t, curr_t + dt).sqrt();
                let curr_mu = hull_white.mu_r(running_r, curr_t, curr_t + dt);
                running_r = curr_mu + curr_vol * norm;
                sum_r = sum_r + running_r * dt;
            });
            let libor_at_option_maturity =
                hull_white.libor_rate_t(running_r, option_maturity, delta);
            //And some more steps since discounted in arrears
            let more_steps = (delta / dt).floor() as usize;
            let new_dt = delta / (more_steps as f64 - 1.0);
            (1..more_steps).for_each(|t_index| {
                let norm = normal.sample(&mut rng_seed);
                let curr_t = new_dt * (t_index as f64) + option_maturity;
                let curr_vol = hull_white.variance_r(curr_t, curr_t + new_dt).sqrt();
                let curr_mu = hull_white.mu_r(running_r, curr_t, curr_t + new_dt);
                running_r = curr_mu + curr_vol * norm;
                sum_r = sum_r + running_r * new_dt;
            });

            if libor_at_option_maturity > strike {
                accum + (libor_at_option_maturity - strike) * ((-sum_r).exp()) //discount
            } else {
                accum
            }
        });
        let average_caplet = delta * (total_sum / (num_sims as f64));
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve);

        let analytical_caplet = hull_white.caplet_now(option_maturity, delta, strike);
        assert_abs_diff_eq!(average_caplet, analytical_caplet, epsilon = 0.0001);
    }

    #[test]
    fn test_edf() {
        let curr_rate = 0.02;
        let sig: f64 = 0.02;
        let a: f64 = 0.3;
        let b = 0.04;
        let future_time = 0.0;
        let option_maturity = 1.5;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct =
                (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let delta = 0.25;
        let seed: [u8; 32] = [2; 32];
        let mut rng_seed = get_rng_seed(seed);
        let normal = StandardNormal;
        let num_sims: usize = 1000000; //hopefully accurate
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve);
        let mu = hull_white.mu_r(curr_rate, future_time, option_maturity);
        let vol = hull_white.variance_r(future_time, option_maturity).sqrt();
        let total_sum = (0..num_sims).fold(0.0, move |accum, _sample_index| {
            let norm = normal.sample(&mut rng_seed);
            let final_r = mu + vol * norm;
            let final_bond =
                hull_white.bond_price_t(final_r, option_maturity, option_maturity + delta);
            accum + 1.0 / final_bond
        });
        let average_edf = ((total_sum / (num_sims as f64)) - 1.0) / delta;

        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve);

        let analytical_edf =
            hull_white.euro_dollar_future_t(curr_rate, future_time, option_maturity, delta);
        assert_abs_diff_eq!(average_edf, analytical_edf, epsilon = 0.0001);
    }

    #[test]
    fn test_bond_now_same_as_t_when_t_is_zero() {
        let curr_rate = 0.02;
        let sig: f64 = 0.02;
        let a: f64 = 0.3;
        let b = 0.04;
        let future_time = 0.0;
        let maturity = 1.5;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct =
                (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve);
        let bond_price_now = hull_white.bond_price_now(maturity);
        let bond_price_t = hull_white.bond_price_t(curr_rate, future_time, maturity);
        assert_abs_diff_eq!(bond_price_now, bond_price_t, epsilon = 0.0000001);
    }
    #[test]
    fn test_coupon_bond_now_same_as_t_when_t_is_zero() {
        let curr_rate = 0.02;
        let sig: f64 = 0.02;
        let a: f64 = 0.3;
        let b = 0.04;
        let delta = 0.25;
        let future_time = 0.0;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct =
                (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let coupon_times = get_coupon_times(6, future_time, delta); //this was 5, but made six since last payment is now included
        let coupon_rate = 0.05 * delta;
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve);
        let bond_price_now = hull_white.coupon_bond_price_now(&coupon_times, coupon_rate);
        let bond_price_t =
            hull_white.coupon_bond_price_t(curr_rate, future_time, &coupon_times, coupon_rate);
        assert_abs_diff_eq!(bond_price_now, bond_price_t, epsilon = 0.0000001);
    }

    #[test]
    fn test_bond_price() {
        let curr_rate = 0.02;
        let sig: f64 = 0.02;
        let a: f64 = 0.3;
        let b = 0.04;
        let future_time = 0.5;
        let option_maturity = 1.5;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct =
                (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve);
        assert_eq!(
            hull_white.bond_price_t(curr_rate, future_time, option_maturity),
            hull_white.bond_price_now(option_maturity - future_time)
        );
    }
    #[test]
    fn test_bond_price_at_expiry() {
        let curr_rate = 0.02;
        let sig: f64 = 0.02;
        let a: f64 = 0.3;
        let b = 0.04;
        let future_time = 0.5;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct =
                (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve);
        assert_eq!(
            hull_white.bond_price_t(curr_rate, future_time, future_time),
            1.0
        );
    }
    #[test]
    fn test_swap() {
        let curr_rate = 0.02;
        let sig: f64 = 0.02;
        let a: f64 = 0.3;
        let b = 0.04;
        let delta = 0.25;
        let future_time = 0.5;
        let swap_maturity = 5.5;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct =
                (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve);
        assert_abs_diff_eq!(
            hull_white
                .swap_price_t(
                    curr_rate,
                    future_time,
                    swap_maturity,
                    delta,
                    hull_white
                        .swap_rate_t(curr_rate, future_time, swap_maturity, delta)
                        .unwrap()
                )
                .unwrap(),
            0.0,
            epsilon = 0.000000001
        );
    }
    #[test]
    fn payer_swaption() {
        let curr_rate = 0.05;
        let sig: f64 = 0.01;
        let a: f64 = 0.05;
        let b = 0.05;
        let delta = 0.25;
        let future_time = 0.0;
        let swap_tenor = 5.0;
        let option_maturity = 1.0;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct =
                (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve);
        let swap_rate = hull_white
            .forward_swap_rate_t(
                curr_rate,
                future_time,
                option_maturity,
                swap_tenor + option_maturity,
                delta,
            )
            .unwrap();
        let analytical = hull_white
            .european_payer_swaption_t(
                curr_rate,
                future_time,
                swap_tenor,
                option_maturity,
                delta,
                swap_rate,
            )
            .unwrap();
        let is_payer = true;

        let tree = hull_white.european_swaption_tree(
            curr_rate,
            future_time,
            swap_tenor,
            option_maturity,
            delta,
            swap_rate,
            is_payer,
            100,
        );
        assert_abs_diff_eq!(analytical, tree, epsilon = 0.0001)
    }
    #[test]
    fn receiver_swaption() {
        let curr_rate = 0.05;
        let sig: f64 = 0.01;
        let a: f64 = 0.05;
        let b = 0.05;
        let delta = 0.25;
        let future_time = 0.0;
        let swap_tenor = 5.0;
        let option_maturity = 1.0;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct =
                (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve);
        let swap_rate = hull_white
            .forward_swap_rate_t(
                curr_rate,
                future_time,
                option_maturity,
                swap_tenor + option_maturity,
                delta,
            )
            .unwrap();
        let analytical = hull_white
            .european_receiver_swaption_t(
                curr_rate,
                future_time,
                swap_tenor,
                option_maturity,
                delta,
                swap_rate,
            )
            .unwrap();
        let is_payer = false;

        let tree = hull_white.european_swaption_tree(
            curr_rate,
            future_time,
            swap_tenor,
            option_maturity,
            delta,
            swap_rate,
            is_payer,
            100,
        );
        assert_abs_diff_eq!(analytical, tree, epsilon = 0.0001)
    }
    #[test]
    fn american_payer_swaption() {
        let curr_rate = 0.05;
        let sig: f64 = 0.01;
        let a: f64 = 0.05;
        let b = 0.05;
        let delta = 0.25;
        let future_time = 0.0;
        let swap_tenor = 5.0;
        let option_maturity = 1.0;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct =
                (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve);
        let swap_rate = hull_white
            .forward_swap_rate_t(
                curr_rate,
                future_time,
                option_maturity,
                swap_tenor + option_maturity,
                delta,
            )
            .unwrap();
        let analytical = hull_white
            .european_payer_swaption_t(
                curr_rate,
                future_time,
                swap_tenor,
                option_maturity,
                delta,
                swap_rate,
            )
            .unwrap();

        let tree = hull_white.american_payer_swaption_t(
            curr_rate,
            future_time,
            swap_tenor,
            option_maturity,
            delta,
            swap_rate,
            100,
        );
        assert_eq!(analytical < tree, true);
    }
    #[test]
    fn american_receiver_swaption() {
        let curr_rate = 0.05;
        let sig: f64 = 0.01;
        let a: f64 = 0.05;
        let b = 0.05;
        let delta = 0.25;
        let future_time = 0.0;
        let swap_tenor = 5.0;
        let option_maturity = 1.0;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct =
                (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve);
        let swap_rate = hull_white
            .forward_swap_rate_t(
                curr_rate,
                future_time,
                option_maturity,
                swap_tenor + option_maturity,
                delta,
            )
            .unwrap();
        let analytical = hull_white
            .european_receiver_swaption_t(
                curr_rate,
                future_time,
                swap_tenor,
                option_maturity,
                delta,
                swap_rate,
            )
            .unwrap();

        let tree = hull_white.american_receiver_swaption_t(
            curr_rate,
            future_time,
            swap_tenor,
            option_maturity,
            delta,
            swap_rate,
            100,
        );
        assert_eq!(analytical < tree, true);
    }
    #[test]
    fn zero_coupon_reference() {
        //http://www.quantcalc.net/BondOption_Vasicek.html
        let curr_rate = 0.01;
        let sig: f64 = 0.03;
        let a = 0.05;
        let b = 0.04;
        let strike = 0.96;
        let future_time = 0.0;
        let bond_maturity = 3.0;
        let option_maturity = 2.0;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct =
                (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve);
        let bond_call = hull_white.bond_call_t(
            curr_rate,
            future_time,
            option_maturity,
            bond_maturity,
            strike,
        );
        assert_abs_diff_eq!(bond_call, 0.033282, epsilon = 0.0001)
    }
    #[test]
    fn zero_coupon_to_coupon() {
        let curr_rate = 0.01;
        let sig: f64 = 0.03;
        let a = 0.05;
        let b = 0.04;
        let strike = 0.96;
        let future_time = 0.0;
        let bond_maturity = 3.0;
        let option_maturity = 2.0;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct =
                (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve);
        let bond_call = hull_white.bond_call_t(
            curr_rate,
            future_time,
            option_maturity,
            bond_maturity,
            strike,
        );
        let coupon_rate = 0.0;
        let coupon_bond_call = hull_white
            .coupon_bond_call_t(
                curr_rate,
                future_time,
                option_maturity,
                &[bond_maturity],
                coupon_rate,
                strike,
            )
            .unwrap();

        assert_abs_diff_eq!(bond_call, coupon_bond_call, epsilon = 0.0001)
    }
}
