//! # Hull-White Interest Rate Model Library
//!
//! This library implements pricing functions for fixed income products using a Hull-White interest rate model.
//! The Hull-White model is a mathematical model describing the evolution of interest rates, commonly used
//! in financial mathematics to price interest rate derivatives.
//!
//! ## Overview
//!
//! The library provides functions to price various fixed income instruments including:
//! - Zero-coupon and coupon bonds
//! - Bond options (calls and puts)
//! - Interest rate caps and floors
//! - Swaps and swaptions
//! - Eurodollar futures
//!
//! ## Key Concepts
//!
//! The fundamental time points in this library are (0, t, T, TM), where:
//! - 0 is the current time (reflective of the current yield curve)
//! - t is some future time for pricing options given the underlying at that time
//! - T and TM represent various asset times (option maturity, bond maturity, etc.)
//!
//! All times are measured with respect to time 0.
//!
//! ## Example Usage
//!
//! ```rust
//! use hull_white::HullWhite;
//!
//! // Define yield and forward curves
//! let yield_curve = |t: f64| 0.05 * t;  // Simple linear yield curve
//! let forward_curve = |t: f64| t.ln();  // Natural log forward curve
//!
//! // Create a Hull-White model with parameters
//! let a = 0.1;      // Mean reversion speed
//! let sigma = 0.01; // Volatility parameter
//! let hull_white = HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
//!
//! // Price a zero-coupon bond maturing in 2 years
//! let bond_price = hull_white.bond_price_now(2.0).unwrap();
//! println!("Bond price: {}", bond_price);
//!
//! // Price a call option on a bond with 2-year maturity, expiring in 1 year, with strike 0.95
//! let option_price = hull_white.bond_call_now(1.0, 2.0, 0.95).unwrap();
//! println!("Call option price: {}", option_price);
//! ```
//!
//! ## Error handling
//!
//! Every pricing function returns `Result<f64, HullWhiteError>`.  An instrument that cannot be
//! priced — an empty or unordered coupon schedule, a coupon already paid at the valuation date, an
//! expired swap, a non-positive tenor, a non-finite argument — is reported as
//! `HullWhiteError::InvalidInput` naming the offending argument, and a computation that overflows
//! to a non-finite value is reported as `HullWhiteError::NumericalError`.  For a bad instrument
//! neither a panic nor a silent `0.0` is a valid answer.
//!
//! ## The Jamshidian solve
//!
//! A coupon-bond option (and so a swaption) is priced with Jamshidian's decomposition, which
//! needs the **critical rate**: the rate at which the underlying bond, valued on the option's
//! expiry date, is worth exactly the strike.  That solve is bracketed analytically from the
//! bond's own leg structure and then safeguarded with bisection, rather than run as the
//! open-ended Newton iteration this crate used to use.
//!
//! Boundary cases return the economically correct answer instead of an error:
//!
//! * `strike = 0` needs no solve at all: the call is the underlying, the put is worthless.
//! * A strike too far away for any rate to reach prices to zero (call) or to parity (put).
//!
//! Tolerance, iteration budget and starting point are configurable with
//! [`HullWhite::with_solver`] / [`SolverSettings`]; the default is a `1e-12` tolerance inside a
//! 100-iteration budget, and prices agree with a direct quadrature of the payoff to about
//! `1e-12`.  When a solve genuinely cannot converge it is `HullWhiteError::RootFindingError`,
//! and the message carries the option maturity, the strike and the seed alongside the bracket,
//! its width and the residual at the point the solver gave up.
//!
//! Every exit from that solve is in **rate** space: the bracket width, and the Newton correction
//! once it is below tolerance, are read as distance-to-root in rate units.  A price residual is
//! never an exit — `f` is a price and can be steep enough that a `1e-8` price residual hides a
//! `7e-3` error in the critical rate, which would put every leg strike of the decomposition
//! wrong while the objective reported "zero".  Because the exit is a statement about how hard the
//! search keeps hunting rather than about how far off the answer is, a loosened tolerance costs
//! iterations rather than capping accuracy: at `1e-7` the price lands on the same value as at
//! `1e-14`.
//!
//! ## Mathematical Foundation
//!
//! The Hull-White model assumes that the short rate follows the stochastic differential equation:
//! `dr(t) = [θ(t) - a*r(t)]dt + σ*dW(t)`
//!
//! where:
//! - `a` is the speed of mean reversion
//! - `σ` is the volatility parameter
//! - `θ(t)` is a time-dependent function calibrated to the initial term structure
//! - `W(t)` is a Wiener process (Brownian motion)

pub mod error;
use error::HullWhiteError;
mod rootfinder;
pub use rootfinder::{Solution, SolverError, SolverSettings};
#[cfg(test)]
mod solver_ab;
mod validation;

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

//This is used for when the number of remaining payments must be derived
//eg when the contract already is in the middle of its life.  If the contract
// start is at time 0.0, and current time is 1.1, and delta is 0.25, then the
//next payment is in .15 units of time.

/// Tolerance (relative to the payment count) for deciding that `(maturity - t) / delta` is a whole
/// number of payment periods.  IEEE-754 rounds a realistic schedule's quotient by a few ulps, i.e.
/// ~1e-15 relative: `(1.1 - 0.1) / 0.1 == 9.999999999999998`, not 10.  1e-9 sits ~6 orders of
/// magnitude above that noise while staying far below the smallest gap (~1e-1 periods) that separates
/// two genuinely different schedules, so it never snaps a real off-schedule date onto a period.
const PAYMENT_PERIOD_TOLERANCE: f64 = 1e-9;

fn get_num_remaining_payments(t: f64, maturity: f64, delta: f64) -> (usize, bool) {
    let raw_payments = (maturity - t) / delta;
    //Never test an f64 quotient for exact integrality: `raw == raw.trunc()` only fires when the
    //schedule happens to be binary-representable, so on real schedules `is_exact` flipped and the
    //swap got anchored a whole coupon period early/late (see `swap_price_t`).  Compare against the
    //nearest integer within tolerance instead; below tolerance means the "next payment" is now.
    let whole_payments = raw_payments.round();
    let tolerance = PAYMENT_PERIOD_TOLERANCE * raw_payments.abs().max(1.0);
    if (raw_payments - whole_payments).abs() <= tolerance {
        (whole_payments as usize, true)
    } else {
        ((raw_payments.floor() + 1.0) as usize, false)
    }
}

/// Builds a payment schedule `t + delta, t + 2*delta, ..., t + num_payments*delta`.
///
/// # Errors
///
/// [`HullWhiteError::InvalidInput`] if `t` is negative or `delta` is not strictly positive.  With
/// `delta == 0` every payment landed on `t`, which read as a schedule and priced as a pile of
/// coincident (and, downstream, divided-by) cash flows.
pub fn get_coupon_times(
    num_payments: usize,
    t: f64,
    delta: f64,
) -> Result<Vec<f64>, HullWhiteError> {
    validation::valuation_time(t)?;
    validation::positive("delta", delta)?;
    Ok((1..(num_payments + 1))
        .map(|index| get_time_from_t_index(index, t, delta))
        .collect())
}

fn get_time_from_t_index(index: usize, t: f64, delta: f64) -> f64 {
    t + (index as f64) * delta
}

fn max_or_zero(v: f64) -> f64 {
    if v > 0.0 { v } else { 0.0 }
}
fn payoff_swaption(is_payer: bool, swp: f64) -> f64 {
    match is_payer {
        true => max_or_zero(swp),
        false => max_or_zero(-swp),
    }
}

//The coupon-sum kernels below are infallible: an empty schedule sums to 0.0 instead of overflowing
//an index, and `is_last` is derived with `index + 1 == len` so there is no `len() - 1` anywhere to
//underflow.  Emptiness and every other bad-instrument case is the public boundary's job (see
//`validation`), which is what lets the root-finding closures below be infallible too.

fn coupon_bond_generic_t(
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
fn coupon_bond_generic_now(
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
pub struct HullWhite<'a, T, U>
where
    T: Fn(f64) -> f64 + std::marker::Sync,
    U: Fn(f64) -> f64 + std::marker::Sync,
{
    a: f64,
    sigma: f64,
    //yield_curve is not divided by time, so this gets perpetually larger (unless rates are negative)
    yield_curve: &'a T,
    forward_curve: &'a U,
    /// Tolerance / iteration budget for the Jamshidian critical-rate solve.  [`HullWhite::init`]
    /// uses [`SolverSettings::default`]; [`HullWhite::with_solver`] changes it.
    solver: SolverSettings,
}

impl<'a, T, U> HullWhite<'a, T, U>
where
    T: Fn(f64) -> f64 + std::marker::Sync,
    U: Fn(f64) -> f64 + std::marker::Sync,
{
    pub fn init(
        a: f64,
        sigma: f64,
        yield_curve: &'a T,
        forward_curve: &'a U,
    ) -> Result<Self, HullWhiteError> {
        Self::validate_parameters(a, sigma)?;
        Ok(Self {
            a,
            sigma,
            yield_curve,
            forward_curve,
            solver: SolverSettings::default(),
        })
    }
    /// Returns a model with a different root-finding budget.
    ///
    /// Deep in-the-money coupon-bond options with long coupon schedules are where the Jamshidian
    /// solve has to travel furthest, so that is where a tighter `tolerance` or a larger
    /// `max_iterations` earns its keep.  A bad configuration is rejected up front rather than
    /// showing up later as an unconverged price.
    ///
    /// # Examples
    ///
    /// ```
    /// use hull_white::{HullWhite, SolverSettings};
    ///
    /// let yield_curve = |t: f64| 0.05 * t;
    /// let forward_curve = |t: f64| t.ln();
    /// let hull_white = HullWhite::init(0.2, 0.3, &yield_curve, &forward_curve).unwrap();
    /// let tuned = hull_white
    ///     .with_solver(SolverSettings {
    ///         tolerance: 1e-14,
    ///         max_iterations: 200,
    ///         initial_guess: None,
    ///     })
    ///     .unwrap();
    /// let price = tuned
    ///     .coupon_bond_call_t(0.04, 1.0, 1.5, &[1.75, 2.0, 2.25], 0.05, 1.0)
    ///     .unwrap();
    /// assert!(price > 0.0);
    /// ```
    #[must_use = "with_solver returns a new model carrying the settings"]
    pub fn with_solver(self, solver: SolverSettings) -> Result<Self, HullWhiteError> {
        solver
            .validate()
            .map_err(|reason| HullWhiteError::InvalidInput(format!("solver settings: {reason}")))?;
        Ok(Self { solver, ..self })
    }
    /// The root-finding configuration currently in use.
    pub fn solver_settings(&self) -> SolverSettings {
        self.solver
    }
    fn validate_parameters(a: f64, sigma: f64) -> Result<(), HullWhiteError> {
        if a <= 0.0 {
            return Err(HullWhiteError::InvalidInput(
                "Mean reversion parameter 'a' must be positive".to_string(),
            ));
        }
        if sigma <= 0.0 {
            return Err(HullWhiteError::InvalidInput(
                "Volatility parameter 'sigma' must be positive".to_string(),
            ));
        }
        if !a.is_finite() {
            return Err(HullWhiteError::InvalidInput(
                "Mean reversion parameter 'a' must be finite".to_string(),
            ));
        }
        if !sigma.is_finite() {
            return Err(HullWhiteError::InvalidInput(
                "Volatility parameter 'sigma' must be finite".to_string(),
            ));
        }
        Ok(())
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
    /// let hull_white= hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let bond_vol = hull_white.t_forward_bond_vol(
    ///     t, t_m, t_f
    /// ).unwrap();
    /// ```
    pub fn t_forward_bond_vol(&self, t: f64, t_m: f64, t_f: f64) -> Result<f64, HullWhiteError> {
        validation::valuation_time(t)?;
        //t_m == t collapses the variance term to zero (and the Black term to a division by zero),
        //and t_f <= t_m flips the sign of the vol, so neither is a volatility at all.
        validation::strictly_after("t_m", t_m, "t", t)?;
        validation::strictly_after("t_f", t_f, "t_m", t_m)?;
        let exp_d = 1.0 - (-self.a * (t_f - t_m)).exp();
        let exp_t = 1.0 - (-2.0 * self.a * (t_m - t)).exp();
        validation::finish(
            "t_forward_bond_vol",
            self.sigma * (exp_t / (2.0 * self.a.powi(3))).sqrt() * exp_d,
        )
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
    /// let hull_white= hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let bond_vol = hull_white.mu_r(r_t, t, t_m).unwrap();
    /// ```
    pub fn mu_r(&self, r_t: f64, t: f64, t_m: f64) -> Result<f64, HullWhiteError> {
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        validation::not_before("t_m", t_m, "t", t)?;
        validation::finish(
            "mu_r",
            self.phi_t(t_m) + (r_t - self.phi_t(t)) * (-self.a * (t_m - t)).exp(),
        )
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let variance = hull_white.variance_r(t, t_m).unwrap();
    /// ```
    pub fn variance_r(&self, t: f64, t_m: f64) -> Result<f64, HullWhiteError> {
        validation::valuation_time(t)?;
        validation::not_before("t_m", t_m, "t", t)?;
        validation::finish(
            "variance_r",
            self.sigma.powi(2) * (1.0 - (-2.0 * self.a * (t_m - t)).exp()) / (2.0 * self.a),
        )
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
    fn bond_price_t_raw(&self, r_t: f64, t: f64, bond_maturity: f64) -> f64 {
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
    /// let hull_white = hull_white::HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
    /// let bond_price = hull_white.bond_price_now(bond_maturity).unwrap();
    /// ```
    pub fn bond_price_now(&self, bond_maturity: f64) -> Result<f64, HullWhiteError> {
        validation::non_negative("bond_maturity", bond_maturity)?;
        validation::finish("bond_price_now", self.bond_price_now_raw(bond_maturity))
    }
    /// Unvalidated counterpart of [`HullWhite::bond_price_now`], for internals that have already
    /// checked their arguments.
    fn bond_price_now_raw(&self, bond_maturity: f64) -> f64 {
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
    //The price of a call option on coupon bond under Hull White...uses jamshidian's trick*
    /// A rigorous bracket on the Jamshidian critical rate, or `None` when it cannot be built (in
    /// which case [`rootfinder::solve`] widens around the seed instead).
    ///
    /// Write the bond's value on the option's expiry date `u` as a sum of zero-coupon legs.  With
    /// `B_i = B(u, T_i)` and `C_i = C(u, T_i)`, and the affine zero-coupon price
    /// `P_i(r) = exp(C_i - B_i r)`, that sum is
    ///
    /// ```text
    /// P_c(r) = sum_i w_i * exp(C_i - B_i * r),   w_i = coupon_rate, w_last = 1 + coupon_rate
    /// ```
    ///
    /// and the objective is `f(r) = P_c(r) - strike`.  Each end of the bracket comes from a
    /// one-sided bound:
    ///
    /// * **low end** — one term of a non-negative sum bounds the whole thing from below, so
    ///   `P_c(r) >= w_last * exp(C_last - B_last r)`.  Solving that for the strike gives an `r`
    ///   where `f >= 0`, valid at any rate.
    /// * **high end** — pulling the largest `C` and a single `B` out of the exponent bounds the
    ///   whole sum by `W * exp(C_max - B r)`, where `W = sum_i w_i` — *provided* `B` is chosen
    ///   for the sign of the rate, `B_min` when `r >= 0` and `B_max` when `r <= 0`.  Solving for
    ///   the strike gives an `r` where `f <= 0`, and taking the `B` whose regime the answer really
    ///   lands in keeps the inequality live.
    ///
    /// This costs one pass over the schedule and no function evaluations at all, which is what
    /// replaces the "iterate blindly from 3% and hope" behaviour.
    fn critical_rate_bracket(
        &self,
        u: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
    ) -> Option<(f64, f64)> {
        if strike <= 0.0 || coupon_times.is_empty() {
            return None;
        }
        let last_index = coupon_times.len() - 1;
        let mut total_weight = 1.0 + coupon_rate; //the redemption-bearing leg
        let mut c_max = f64::NEG_INFINITY;
        let mut b_min = f64::INFINITY;
        let mut b_max: f64 = 0.0;
        for (index, coupon_time) in coupon_times.iter().enumerate() {
            if index != last_index {
                total_weight += coupon_rate;
            }
            let b = at_t(self.a, u, *coupon_time);
            let c = ct_t(
                self.a,
                self.sigma,
                u,
                *coupon_time,
                self.yield_curve,
                self.forward_curve,
            );
            c_max = c_max.max(c);
            b_min = b_min.min(b);
            b_max = b_max.max(b);
        }
        let w_last = 1.0 + coupon_rate;
        let last_time = coupon_times[last_index];
        let b_last = at_t(self.a, u, last_time);
        let c_last = ct_t(
            self.a,
            self.sigma,
            u,
            last_time,
            self.yield_curve,
            self.forward_curve,
        );
        //Both bounds need non-negative weights.  A coupon_rate at or below -100% breaks them, so
        //hand that nonsense to the widening path instead of trusting a bracket built on sand.
        if !(b_min > 0.0 && b_max > 0.0 && total_weight > 0.0 && w_last > 0.0) {
            return None;
        }
        let log_strike = strike.ln();
        //Differences of logs, not a log of a ratio: a strike of 1e-320 would overflow 1/ratio.
        let lower = (c_last + (w_last.ln() - log_strike)) / b_last;
        let numerator = c_max + total_weight.ln() - log_strike;
        let b_for_regime = if numerator >= 0.0 { b_min } else { b_max };
        let upper = numerator / b_for_regime;
        //exp/log rounding leaves each end a few ulp off `f = 0` rather than safely beyond it, so
        //shove both ends outwards.  Outward is always safe: `P_c` falls with the rate at the low
        //end, and at the high end it is already comfortably under the strike.
        let nudge = |x: f64| 1e-7 * (1.0 + x.abs());
        let (lower, upper) = (lower - nudge(lower), upper + nudge(upper));
        if lower.is_finite() && upper.is_finite() && lower <= upper {
            Some((lower, upper))
        } else {
            None
        }
    }

    /// Solves for the critical rate of the Jamshidian decomposition.
    ///
    /// The critical rate `r*` is the rate at which the coupon bond, valued at the option's
    /// expiry date, is worth exactly the strike.  Striking every zero-coupon leg at its own price
    /// at that rate makes the leg portfolio level with the whole bond, which is what lets a sum
    /// of leg options equal the bond option.
    ///
    /// The seed, tolerance and iteration budget come from `SolverSettings`; the bracket comes
    /// from `critical_rate_bracket` where the schedule allows it, and from geometric widening
    /// otherwise.
    fn solve_critical_rate(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
    ) -> Result<rootfinder::Solution, HullWhiteError> {
        let objective = |rate: f64| {
            coupon_bond_generic_t(
                rate,
                option_maturity,
                coupon_times,
                coupon_rate,
                &|leg_rate: f64, leg_t: f64, bond_maturity: f64| {
                    self.bond_price_t_raw(leg_rate, leg_t, bond_maturity)
                },
            ) - strike
        };
        let derivative = |rate: f64| {
            self.coupon_bond_price_t_deriv(rate, option_maturity, coupon_times, coupon_rate)
        };
        //Seed with the model's own expected short rate at expiry: same units as the unknown,
        //curve aware and state aware, where a fixed guess can sit a very long way from the answer.
        let seed = match self.solver.initial_guess {
            Some(guess) => guess,
            None => self.mu_r(r_t, t, option_maturity)?,
        };
        let bracket =
            self.critical_rate_bracket(option_maturity, coupon_times, coupon_rate, strike);
        rootfinder::solve(&objective, &derivative, bracket, seed, &self.solver).map_err(|error| {
            HullWhiteError::RootFindingError(format!(
                "critical rate for the Jamshidian decomposition of a coupon-bond option \
                 (option_maturity = {option_maturity}, strike = {strike}, seed = {seed}): {error}"
            ))
        })
    }

    fn coupon_bond_option_generic_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
        generic_fn: &impl Fn(f64, f64, f64, f64, f64) -> Result<f64, HullWhiteError>,
    ) -> Result<f64, HullWhiteError> {
        let par_value = 1.0;
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        validation::strictly_after("option_maturity", option_maturity, "t", t)?;
        validation::finite("coupon_rate", coupon_rate)?;
        validation::strike(strike)?;
        //Jamshidian's decomposition prices the option against the bond *as it stands at expiry*, so
        //every coupon in the schedule has to be paid after `option_maturity`.  A coupon paid before
        //expiry is not part of the underlying at expiry; taking it through here priced a leg that no
        //longer exists (see the follow-up issue for supporting such bonds properly).
        validation::payment_schedule(coupon_times, option_maturity)?;
        //Each leg is struck at the zero-coupon price that puts the *whole* bond level with the
        //strike at the critical rate.  A strike of zero means that critical rate sits at
        //+infinity, where every modified strike is zero -- and Black-Scholes is exact there (a call
        //struck at nothing on nothing is worth the underlying, the put is worth nothing), so the
        //solve is skipped rather than chased out to infinity.
        let critical_rate = if strike == 0.0 {
            None
        } else {
            Some(self.solve_critical_rate(
                r_t,
                t,
                option_maturity,
                coupon_times,
                coupon_rate,
                strike,
            )?)
        };
        let last_index = coupon_times.len();
        coupon_times
            .iter()
            .enumerate()
            .map(|(index, coupon_time)| {
                let is_last = index + 1 == last_index;
                let strike_leg = match critical_rate {
                    Some(solution) => {
                        self.bond_price_t_raw(solution.root, option_maturity, *coupon_time)
                    }
                    None => 0.0,
                };
                generic_fn(r_t, t, option_maturity, *coupon_time, strike_leg)
                    .map(|leg| leg * (coupon_rate + if is_last { par_value } else { 0.0 }))
            })
            .sum::<Result<f64, HullWhiteError>>()
            .and_then(|price| {
                validation::finish("Jamshidian decomposition of a coupon-bond option", price)
            })
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
    /// Unvalidated swap pricing, used by the swaption tree: the tree callback has to return an `f64`
    /// (`binomial_tree`'s contract), so it cannot propagate a `Result`, and its inputs are checked
    /// once at the public boundary instead of once per node.
    fn swap_price_t_init_raw(
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
mod tests {
    use super::*;
    use approx::*;
    use rand::distributions::{Distribution, StandardNormal};
    use rand::{SeedableRng, StdRng};
    fn get_rng_seed(seed: [u8; 32]) -> StdRng {
        SeedableRng::from_seed(seed)
    }
    #[test]
    fn test_get_num_payments_if_exact_integer() {
        let t = 0.5;
        let maturity = 2.0;
        let delta = 0.25;
        let (num_payments, is_exact) = get_num_remaining_payments(t, maturity, delta);
        assert_eq!(num_payments, 6);
        assert!(is_exact);
        let next_exchange_date = maturity - (num_payments as f64 - 1.0) * delta;
        assert_abs_diff_eq!(next_exchange_date, t + delta, epsilon = 0.0000001);
    }

    #[test]
    fn test_get_num_payments_not_even() {
        let t = 0.5;
        let maturity = 2.0;
        let delta = 0.4;
        let (num_payments, is_exact) = get_num_remaining_payments(t, maturity, delta);
        assert_eq!(num_payments, 4);
        assert_eq!(is_exact, false);
        let next_exchange_date = maturity - (num_payments as f64 - 1.0) * delta;
        assert_abs_diff_eq!(next_exchange_date, 0.8, epsilon = 0.0000001);
    }

    /// Whole-number-of-period schedules that are NOT binary-representable.  The f64 quotient drifts to
    /// either side of the integer, and each direction broke differently under `== trunc()`:
    ///   low  e.g. (1.1 - 0.1) / 0.1   = 9.999999999999998  -> is_exact flipped false, so the swap
    ///        got anchored at maturity - (n-1)*delta = 0.2 = t + delta, one period late
    ///   high e.g. (1.3000000000000003 - 0.1) / 0.1 = 12.000000000000002 -> floor+1 gave n+1,
    ///        one spurious payment period
    #[test]
    fn num_remaining_payments_tolerates_whole_schedules_that_are_not_binary_exact() {
        let whole_schedules: [(f64, f64, f64, usize); 8] = [
            (0.1, 1.1, 0.1, 10), //  9.999999999999998
            (0.1, 2.0, 0.1, 19), // 18.999999999999996
            (0.3, 2.3, 0.2, 10),
            (0.25, 2.25, 0.5, 4),
            (0.1, 0.30000000000000004, 0.1, 2), //  2.0000000000000004
            (0.1, 0.4, 0.1, 3),                 //  3.0000000000000004
            (0.1, 0.7000000000000001, 0.1, 6),  //  6.000000000000001
            (0.1, 1.3000000000000003, 0.1, 12), // 12.000000000000002
        ];
        for (t, maturity, delta, expected_payments) in whole_schedules {
            let raw = (maturity - t) / delta;
            let (num_payments, is_exact) = get_num_remaining_payments(t, maturity, delta);
            assert_eq!(
                num_payments, expected_payments,
                "({maturity} - {t}) / {delta} = {raw:?} should be {expected_payments} payments"
            );
            assert!(
                is_exact,
                "({maturity} - {t}) / {delta} = {raw:?} is a whole schedule and must read as exact"
            );
            //Whole schedule => the anchor the swap derives must be t itself (payments at t+delta ...).
            let next_exchange_date = maturity - (num_payments as f64 - 1.0) * delta;
            assert_abs_diff_eq!(next_exchange_date, t + delta, epsilon = 1e-9);
        }
    }

    /// Guard the other way: a mid-life schedule that really is between payment dates must NOT be
    /// snapped onto a period boundary by the new tolerance.
    #[test]
    fn num_remaining_payments_does_not_snap_genuinely_off_schedule_dates() {
        //      t, maturity, delta, expected payments, expected first remaining payment
        let mid_life: [(f64, f64, f64, usize, f64); 4] = [
            (0.5, 2.0, 0.4, 4, 0.8),    //  3.75
            (0.1, 1.15, 0.1, 11, 0.15), // 10.5
            (0.0, 1.0, 0.3, 4, 0.1),    //  3.3333333333333335
            (0.05, 1.0, 0.3, 4, 0.1),   //  3.1666666666666665
        ];
        for (t, maturity, delta, expected_payments, expected_first_payment) in mid_life {
            let raw = (maturity - t) / delta;
            let (num_payments, is_exact) = get_num_remaining_payments(t, maturity, delta);
            assert_eq!(
                num_payments, expected_payments,
                "({maturity} - {t}) / {delta} = {raw:?} should be {expected_payments} payments"
            );
            assert!(
                !is_exact,
                "({maturity} - {t}) / {delta} = {raw:?} is off-schedule and must not be exact"
            );
            let next_exchange_date = maturity - (num_payments as f64 - 1.0) * delta;
            assert_abs_diff_eq!(next_exchange_date, expected_first_payment, epsilon = 1e-9);
        }
    }

    /// The float comparison was not cosmetic: `swap_price_t` selects `swap_start` from `is_exact`, so a
    /// whole-but-not-binary schedule priced the swap from the wrong anchor.  With an integral number of
    /// periods remaining, `swap_price_t` must be the same call as `swap_price_t_init(..., t, n, ...).unwrap()`.
    #[test]
    fn swap_price_t_anchors_at_t_for_whole_non_binary_schedules() {
        //Needs a curve that moves: on a flat curve a one-period anchor shift is nearly free, which is
        //why the legacy fixture never showed this.
        let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
        let hull_white = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve).unwrap();
        let r_t = STEEP_CURR_RATE;
        let swap_rate = 0.045;
        for (t, swap_maturity, delta, num_payments) in [
            (0.1, 1.1, 0.1, 10),
            (0.1, 2.0, 0.1, 19),
            (0.3, 2.3, 0.2, 10),
            (0.25, 2.25, 0.5, 4),
        ] {
            let derived = hull_white
                .swap_price_t(r_t, t, swap_maturity, delta, swap_rate)
                .unwrap();
            let anchored_at_t = hull_white
                .swap_price_t_init(r_t, t, t, num_payments, delta, swap_rate)
                .unwrap();
            let anchored_one_period_late = hull_white
                .swap_price_t_init(r_t, t, t + delta, num_payments, delta, swap_rate)
                .unwrap();
            assert_abs_diff_eq!(derived, anchored_at_t, epsilon = 1e-12);
            //Sanity: the wrong anchor really is a materially different number, so the assertion above
            //has teeth on this fixture (measured gap is ~1e-2 on a ~1e-1 swap price).
            assert!(
                (anchored_at_t - anchored_one_period_late).abs() > 1e-4,
                "fixture too flat to detect an anchor shift at t={t}: {anchored_at_t} vs \
                 {anchored_one_period_late}"
            );
        }
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
        let coupon_times = get_coupon_times(num_payments, t, delta).unwrap();
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
        let coupon_times = get_coupon_times(num_payments, t, delta).unwrap();
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
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        let caplet_n = hull_white
            .caplet_now(option_maturity, delta, strike)
            .unwrap();
        let caplet = hull_white
            .caplet_t(curr_rate, future_time, option_maturity, delta, strike)
            .unwrap();
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
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        let libor_n = hull_white.forward_libor_rate_now(maturity, delta).unwrap();
        let libor_t = hull_white
            .forward_libor_rate_t(curr_rate, future_time, maturity, delta)
            .unwrap();
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
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        let total_sum = (0..num_sims).fold(0.0, |accum, _sample_index| {
            let mut sum_r = 0.0;
            let mut running_r = curr_rate;
            let dt = (option_maturity - future_time) / (num_discrete_steps as f64 - 1.0);
            (0..num_discrete_steps).for_each(|t_index| {
                let norm = normal.sample(&mut rng_seed);
                let curr_t = dt * (t_index as f64) + future_time;
                let curr_vol = hull_white.variance_r(curr_t, curr_t + dt).unwrap().sqrt();
                let curr_mu = hull_white.mu_r(running_r, curr_t, curr_t + dt).unwrap();
                running_r = curr_mu + curr_vol * norm;
                sum_r = sum_r + running_r * dt;
            });
            let libor_at_option_maturity = hull_white
                .libor_rate_t(running_r, option_maturity, delta)
                .unwrap();
            //And some more steps since discounted in arrears
            let more_steps = (delta / dt).floor() as usize;
            let new_dt = delta / (more_steps as f64 - 1.0);
            (1..more_steps).for_each(|t_index| {
                let norm = normal.sample(&mut rng_seed);
                let curr_t = new_dt * (t_index as f64) + option_maturity;
                let curr_vol = hull_white
                    .variance_r(curr_t, curr_t + new_dt)
                    .unwrap()
                    .sqrt();
                let curr_mu = hull_white.mu_r(running_r, curr_t, curr_t + new_dt).unwrap();
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
        let analytical_caplet = hull_white
            .caplet_now(option_maturity, delta, strike)
            .unwrap();
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
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        let mu = hull_white
            .mu_r(curr_rate, future_time, option_maturity)
            .unwrap();
        let vol = hull_white
            .variance_r(future_time, option_maturity)
            .unwrap()
            .sqrt();
        let total_sum = (0..num_sims).fold(0.0, |accum, _sample_index| {
            let norm = normal.sample(&mut rng_seed);
            let final_r = mu + vol * norm;
            let final_bond = hull_white
                .bond_price_t(final_r, option_maturity, option_maturity + delta)
                .unwrap();
            accum + 1.0 / final_bond
        });
        let average_edf = ((total_sum / (num_sims as f64)) - 1.0) / delta;

        let analytical_edf = hull_white
            .euro_dollar_future_t(curr_rate, future_time, option_maturity, delta)
            .unwrap();
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
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        let bond_price_now = hull_white.bond_price_now(maturity).unwrap();
        let bond_price_t = hull_white
            .bond_price_t(curr_rate, future_time, maturity)
            .unwrap();
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
        let coupon_times = get_coupon_times(6, future_time, delta).unwrap(); //this was 5, but made six since last payment is now included
        let coupon_rate = 0.05 * delta;
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        let bond_price_now = hull_white
            .coupon_bond_price_now(&coupon_times, coupon_rate)
            .unwrap();
        let bond_price_t = hull_white
            .coupon_bond_price_t(curr_rate, future_time, &coupon_times, coupon_rate)
            .unwrap();
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
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        assert_eq!(
            hull_white
                .bond_price_t(curr_rate, future_time, option_maturity)
                .unwrap(),
            hull_white
                .bond_price_now(option_maturity - future_time)
                .unwrap()
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
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        assert_eq!(
            hull_white
                .bond_price_t(curr_rate, future_time, future_time)
                .unwrap(),
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
        let num_swap_payments = 20;
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
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        assert_abs_diff_eq!(
            hull_white
                .swap_price_t(
                    curr_rate,
                    future_time,
                    swap_maturity,
                    delta,
                    hull_white
                        .swap_rate_t(curr_rate, future_time, num_swap_payments, delta)
                        .unwrap()
                )
                .unwrap(),
            0.0,
            epsilon = 0.000000001
        );
    }
    #[test]
    fn test_swap_init() {
        let curr_rate = 0.02;
        let sig: f64 = 0.02;
        let a: f64 = 0.3;
        let b = 0.04;
        let delta = 0.25;
        let future_time = 0.5;
        let swap_maturity = 5.5;
        let num_swap_payments = 20;
        let swap_rate = 0.03;
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
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        let sp_init = hull_white
            .swap_price_t_init(
                curr_rate,
                future_time,
                future_time,
                num_swap_payments,
                delta,
                swap_rate,
            )
            .unwrap();
        let sp = hull_white
            .swap_price_t(curr_rate, future_time, swap_maturity, delta, swap_rate)
            .unwrap();
        assert_eq!(sp_init, sp);
    }

    /// Curves consistent with a Hull-White process (Vasicek short rate with mean reversion `a`
    /// towards `b`), so that `bond_price_now`/`bond_price_t` are exact for the same model.
    fn hw_curves(
        curr_rate: f64,
        a: f64,
        b: f64,
        sig: f64,
    ) -> (impl Fn(f64) -> f64, impl Fn(f64) -> f64) {
        let yield_curve = move |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct =
                (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = move |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        (yield_curve, forward_curve)
    }

    /// A deliberately steep (strongly time-inhomogeneous) calibration: the short rate sits far below
    /// the long-run mean, so phi(t) moves a lot over the option's life. The existing fixture uses
    /// curr_rate == b, which makes phi(t) exactly constant and hides any time-coordinate error.
    const STEEP_CURR_RATE: f64 = 0.02;
    const STEEP_SIG: f64 = 0.03;
    const STEEP_A: f64 = 0.2;
    const STEEP_B: f64 = 0.06;

    #[test]
    fn steep_fixture_actually_makes_phi_time_dependent() {
        //Protects the tests above: if phi were constant, a time-coordinate error would cancel out
        //and the regression tests would pass for the wrong reason.
        let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
        let hull_white = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve).unwrap();
        let phi_0 = hull_white.phi_t(0.0);
        let phi_5 = hull_white.phi_t(5.0);
        assert!(
            (phi_5 - phi_0).abs() > 0.005,
            "steep fixture must make phi(t) time dependent: phi(0)={phi_0}, phi(5)={phi_5}"
        );
    }

    #[test]
    fn european_swaption_tree_matches_analytic_when_t_is_zero() {
        //Guard: when the valuation time is 0 the shifted tree clock and the absolute model clock
        //coincide, so this pins the (already correct) behaviour across the time-coordinate fix.
        let curr_rate = STEEP_CURR_RATE;
        let sig = STEEP_SIG;
        let a = STEEP_A;
        let b = STEEP_B;
        let delta = 0.25;
        let future_time = 0.0;
        let option_maturity = 1.5;
        let num_swap_payments = 20;
        let (yield_curve, forward_curve) = hw_curves(curr_rate, a, b, sig);
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        let swap_rate = hull_white
            .forward_swap_rate_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
            )
            .unwrap();
        //Tree convergence for this fixture: |tree - analytic| is ~1.4e-4 at 50 steps, ~1.1e-4 at
        //100, ~7.0e-5 at 200, ~3.1e-5 at 400, ~4e-7 at 800.  1e-4 at 400 steps gives >3x headroom
        //over the measured discretisation noise while still being ~40x tighter than the time-shift bug.
        let steps = 400;
        let payer = hull_white
            .european_payer_swaption_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
                swap_rate,
            )
            .unwrap();
        let tree_payer = hull_white.european_swaption_tree(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            true,
            steps,
        );
        assert_abs_diff_eq!(payer, tree_payer, epsilon = 0.0001);
    }

    /// Bit-exact guard on the `t = 0` path.  The values below were produced by the pre-fix pricers at
    /// commit 254d089 (`option_maturity` 1.5, `num_swap_payments` 20, `delta` 0.25, ATM-forward
    /// strike, 400 tree steps), captured with a round-trip-checked `{:?}` print -- note `{:.17}` is
    /// 17 digits *after the decimal point*, which for these magnitudes loses the last bit.
    /// At `t = 0` the shifted tree clock and the absolute model clock coincide, so the
    /// time-coordinate fix must leave every one of them bit-identical -- pinning bits rather than a
    /// tolerance means a real behaviour change at `t = 0` cannot slip through unnoticed, and cannot
    /// be "absorbed" by widening an epsilon later.
    #[test]
    fn swaption_tree_at_t_is_bit_identical_to_pre_fix() {
        let option_maturity = 1.5;
        let num_swap_payments = 20;
        let delta = 0.25;
        let steps = 400;
        //             name     r0     a     b     sig   eur payer        eur receiver     amer payer   amer receiver
        let golden: [(&str, f64, f64, f64, f64, f64, f64, f64, f64); 2] = [
            (
                "legacy",
                0.05,
                0.05,
                0.05,
                0.01,
                0.017330477644662997,
                0.017329771203617984,
                0.01834265355592532,
                0.017797483448434452,
            ),
            (
                "steep",
                0.02,
                0.2,
                0.06,
                0.03,
                0.03597512274296273,
                0.03597334589011368,
                0.03781406402902323,
                0.04452822113093644,
            ),
        ];
        for (name, r0, a, b, sig, eur_p, eur_r, amer_p, amer_r) in golden {
            let (yield_curve, forward_curve) = hw_curves(r0, a, b, sig);
            let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
            let swap_rate = hull_white
                .forward_swap_rate_t(r0, 0.0, option_maturity, num_swap_payments, delta)
                .unwrap();
            assert_bits_eq(
                hull_white.european_swaption_tree(
                    r0,
                    0.0,
                    option_maturity,
                    num_swap_payments,
                    delta,
                    swap_rate,
                    true,
                    steps,
                ),
                eur_p,
                &format!("{name} european payer tree"),
            );
            assert_bits_eq(
                hull_white.european_swaption_tree(
                    r0,
                    0.0,
                    option_maturity,
                    num_swap_payments,
                    delta,
                    swap_rate,
                    false,
                    steps,
                ),
                eur_r,
                &format!("{name} european receiver tree"),
            );
            assert_bits_eq(
                hull_white
                    .american_payer_swaption_t(
                        r0,
                        0.0,
                        option_maturity,
                        num_swap_payments,
                        delta,
                        swap_rate,
                        steps,
                    )
                    .unwrap(),
                amer_p,
                &format!("{name} american payer"),
            );
            assert_bits_eq(
                hull_white
                    .american_receiver_swaption_t(
                        r0,
                        0.0,
                        option_maturity,
                        num_swap_payments,
                        delta,
                        swap_rate,
                        steps,
                    )
                    .unwrap(),
                amer_r,
                &format!("{name} american receiver"),
            );
        }
    }

    /// Golden comparison helper: asserts the two f64 have identical bit patterns, so "unchanged" here
    /// really means unchanged, not "within a tolerance that can be nudged".
    #[allow(clippy::float_cmp)]
    fn assert_bits_eq(actual: f64, expected: f64, what: &str) {
        assert_eq!(
            actual.to_bits(),
            expected.to_bits(),
            "{what} moved at t = 0: actual={actual:.17} expected={expected:.17}"
        );
    }

    #[test]
    fn european_swaption_tree_matches_analytic_when_t_is_nonzero() {
        //Regression for the shifted-vs-absolute time bug: the tree runs on the clock
        //`option_maturity - t` but phi() and the swap legs need absolute time from "now" (0).
        let curr_rate = STEEP_CURR_RATE;
        let sig = STEEP_SIG;
        let a = STEEP_A;
        let b = STEEP_B;
        let delta = 0.25;
        let future_time = 0.5;
        let option_maturity = 1.5;
        let num_swap_payments = 20;
        let (yield_curve, forward_curve) = hw_curves(curr_rate, a, b, sig);
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        let swap_rate = hull_white
            .forward_swap_rate_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
            )
            .unwrap();
        //Same budget as the t = 0 guard above.  With the tree clocking the shifted time into phi and
        //into the swap legs, this diff plateaus at ~4.0e-3 (payer) / ~5.0e-3 (receiver) no matter how
        //many steps are used -- a systematic pricing error of ~13-16% of the option value, not
        //discretisation noise.  After the fix the residual is ~2e-5 at 400 steps.
        let steps = 400;
        let payer = hull_white
            .european_payer_swaption_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
                swap_rate,
            )
            .unwrap();
        let receiver = hull_white
            .european_receiver_swaption_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
                swap_rate,
            )
            .unwrap();
        let tree_payer = hull_white.european_swaption_tree(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            true,
            steps,
        );
        let tree_receiver = hull_white.european_swaption_tree(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            false,
            steps,
        );
        assert_abs_diff_eq!(payer, tree_payer, epsilon = 0.0001);
        assert_abs_diff_eq!(receiver, tree_receiver, epsilon = 0.0001);
    }

    #[test]
    fn american_swaption_tree_at_nonzero_t_carries_a_positive_early_exercise_premium() {
        //The American tree shares the time-coordinate fix with the European tree; this pins that the
        //American price at t > 0 is a sane number above the European analytic price rather than a
        //shifted-clock artefact.
        let delta = 0.25;
        let future_time = 0.5;
        let option_maturity = 1.5;
        let num_swap_payments = 20;
        let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
        let hull_white = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve).unwrap();
        let swap_rate = hull_white
            .forward_swap_rate_t(
                STEEP_CURR_RATE,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
            )
            .unwrap();
        let american = |is_payer: bool, steps: usize| {
            if is_payer {
                hull_white
                    .american_payer_swaption_t(
                        STEEP_CURR_RATE,
                        future_time,
                        option_maturity,
                        num_swap_payments,
                        delta,
                        swap_rate,
                        steps,
                    )
                    .unwrap()
            } else {
                hull_white
                    .american_receiver_swaption_t(
                        STEEP_CURR_RATE,
                        future_time,
                        option_maturity,
                        num_swap_payments,
                        delta,
                        swap_rate,
                        steps,
                    )
                    .unwrap()
            }
        };
        for is_payer in [true, false] {
            let side = if is_payer { "payer" } else { "receiver" };
            let european = if is_payer {
                hull_white
                    .european_payer_swaption_t(
                        STEEP_CURR_RATE,
                        future_time,
                        option_maturity,
                        num_swap_payments,
                        delta,
                        swap_rate,
                    )
                    .unwrap()
            } else {
                hull_white
                    .european_receiver_swaption_t(
                        STEEP_CURR_RATE,
                        future_time,
                        option_maturity,
                        num_swap_payments,
                        delta,
                        swap_rate,
                    )
                    .unwrap()
            };
            let american_200 = american(is_payer, 200);
            let american_400 = american(is_payer, 400);
            assert!(
                american_200.is_finite() && american_400.is_finite(),
                "{side} swaption: non-finite price: 200={american_200}, 400={american_400}"
            );
            assert!(
                american_400 > european,
                "{side} swaption: early exercise premium must be positive: european={european}, \
                 american(400)={american_400}"
            );
            //Converging, not oscillating: refining the tree moves the price by less than the premium.
            assert!(
                (american_400 - american_200).abs() < (american_400 - european),
                "{side} swaption: american tree not converging: 200={american_200}, \
                 400={american_400}, european={european}"
            );
        }
    }

    #[test]
    fn payer_swaption() {
        let curr_rate = 0.05;
        let sig: f64 = 0.01;
        let a: f64 = 0.05;
        let b = 0.05;
        let delta = 0.25;
        let future_time = 0.0;
        //let swap_tenor = 5.0;
        let num_swap_payments = 20;
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
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        let swap_rate = hull_white
            .forward_swap_rate_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
            )
            .unwrap();
        let analytical = hull_white
            .european_payer_swaption_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
                swap_rate,
            )
            .unwrap();
        let is_payer = true;

        let tree = hull_white.european_swaption_tree(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
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
        //let swap_tenor = 5.0;
        let num_swap_payments = 20;
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
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        let swap_rate = hull_white
            .forward_swap_rate_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
            )
            .unwrap();
        let analytical = hull_white
            .european_receiver_swaption_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
                swap_rate,
            )
            .unwrap();
        let is_payer = false;

        let tree = hull_white.european_swaption_tree(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
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
        //let swap_tenor = 5.0;
        let num_swap_payments = 20;
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
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        let swap_rate = hull_white
            .forward_swap_rate_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
            )
            .unwrap();
        let analytical = hull_white
            .european_payer_swaption_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
                swap_rate,
            )
            .unwrap();

        let tree = hull_white
            .american_payer_swaption_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
                swap_rate,
                100,
            )
            .unwrap();
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
        //let swap_tenor = 5.0;
        let num_swap_payments = 20;
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
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        let swap_rate = hull_white
            .forward_swap_rate_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
            )
            .unwrap();
        let analytical = hull_white
            .european_receiver_swaption_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
                swap_rate,
            )
            .unwrap();

        let tree = hull_white
            .american_receiver_swaption_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
                swap_rate,
                100,
            )
            .unwrap();
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
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        let bond_call = hull_white
            .bond_call_t(
                curr_rate,
                future_time,
                option_maturity,
                bond_maturity,
                strike,
            )
            .unwrap();
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
        let hull_white = HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
        let bond_call = hull_white
            .bond_call_t(
                curr_rate,
                future_time,
                option_maturity,
                bond_maturity,
                strike,
            )
            .unwrap();
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

    // ---------------------------------------------------------------------------
    // Invalid / empty instrument inputs must error, not panic and not return a number.
    //
    // Before this change three shapes of failure were live:
    //   * `coupon_times.len() - 1` on an empty schedule -> usize underflow panic;
    //   * an expired swap / `delta <= 0` -> the float-to-usize casts saturated at 0, the payment
    //     loop vanished, and a bond leg came back out looking like a swap price;
    //   * NaN/inf travelled straight through `exp`/`/` into the returned price.
    // The guards under `validation::` now turn each of those into `InvalidInput` naming the
    // argument, and a non-finite *computed* result into `NumericalError`.
    // ---------------------------------------------------------------------------

    macro_rules! yvf_setup {
        ($model:ident) => {
            let (yield_curve, forward_curve) =
                hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
            let $model = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve).unwrap();
        };
    }

    fn expect_invalid<T: std::fmt::Debug>(result: Result<T, HullWhiteError>, needle: &str) {
        match result {
            Err(HullWhiteError::InvalidInput(msg)) => assert!(
                msg.contains(needle),
                "expected the error to name {needle:?}, got: {msg}"
            ),
            Err(other) => panic!("expected InvalidInput naming {needle:?}, got {other:?}"),
            Ok(value) => panic!("expected InvalidInput naming {needle:?}, got Ok({value:?})"),
        }
    }

    #[test]
    fn empty_coupon_times_is_invalid_input() {
        yvf_setup!(hull_white);
        expect_invalid(
            hull_white.coupon_bond_price_t(0.05, 1.0, &[], 0.05),
            "coupon_times is empty",
        );
        expect_invalid(
            hull_white.coupon_bond_price_now(&[], 0.05),
            "coupon_times is empty",
        );
        expect_invalid(
            hull_white.coupon_bond_call_t(0.05, 1.0, 1.5, &[], 0.05, 1.0),
            "coupon_times is empty",
        );
        expect_invalid(
            hull_white.coupon_bond_put_t(0.05, 1.0, 1.5, &[], 0.05, 1.0),
            "coupon_times is empty",
        );
    }

    #[test]
    fn unascending_coupon_times_is_invalid_input() {
        yvf_setup!(hull_white);
        //descending step at index 1
        expect_invalid(
            hull_white.coupon_bond_price_t(0.05, 1.0, &[1.25, 1.0, 1.5], 0.05),
            "coupon_times[1]",
        );
        //a duplicated date is not ascending either; it silently double-weights a leg
        expect_invalid(
            hull_white.coupon_bond_price_t(0.05, 1.0, &[1.25, 1.5, 1.5], 0.05),
            "coupon_times[2]",
        );
    }

    #[test]
    fn coupon_time_at_or_before_valuation_time_is_invalid_input() {
        yvf_setup!(hull_white);
        //exactly on the valuation date: already paid
        expect_invalid(
            hull_white.coupon_bond_price_t(0.05, 1.0, &[1.25, 1.0], 0.05),
            "coupon_times[1]",
        );
        //in the past
        expect_invalid(
            hull_white.coupon_bond_price_t(0.05, 1.0, &[0.5, 2.0], 0.05),
            "coupon_times[0]",
        );
        //for the `now` entry point the valuation time is 0, so a payment at 0 is already made
        expect_invalid(
            hull_white.coupon_bond_price_now(&[0.0, 2.0], 0.05),
            "coupon_times[0]",
        );
    }

    #[test]
    fn expired_swap_is_invalid_input() {
        yvf_setup!(hull_white);
        //matured before the valuation date
        expect_invalid(
            hull_white.swap_price_t(0.05, 1.0, 0.75, 0.25, 0.04),
            "swap_maturity",
        );
        //matures exactly on the valuation date: nothing left to price
        expect_invalid(
            hull_white.swap_price_t(0.05, 1.0, 1.0, 0.25, 0.04),
            "swap_maturity",
        );
    }

    #[test]
    fn non_positive_delta_is_invalid_input() {
        yvf_setup!(hull_white);
        for delta in [0.0, -0.25] {
            expect_invalid(get_coupon_times(4, 1.0, delta), "delta");
            expect_invalid(
                hull_white.swap_price_t(0.05, 1.0, 3.0, delta, 0.04),
                "delta",
            );
            expect_invalid(
                hull_white.swap_price_t_init(0.05, 1.0, 1.0, 8, delta, 0.04),
                "delta",
            );
            expect_invalid(hull_white.swap_rate_t(0.05, 1.0, 8, delta), "delta");
            expect_invalid(
                hull_white.forward_swap_rate_t(0.05, 1.0, 1.5, 8, delta),
                "delta",
            );
            expect_invalid(hull_white.caplet_t(0.05, 1.0, 1.5, delta, 0.04), "delta");
            expect_invalid(
                hull_white.euro_dollar_future_t(0.05, 1.0, 1.5, delta),
                "delta",
            );
            expect_invalid(hull_white.libor_rate_t(0.05, 1.0, delta), "delta");
            expect_invalid(
                hull_white.european_payer_swaption_t(0.05, 1.0, 1.5, 8, delta, 0.04),
                "delta",
            );
        }
    }

    #[test]
    fn zero_period_instrument_is_invalid_input() {
        yvf_setup!(hull_white);
        expect_invalid(
            hull_white.swap_price_t_init(0.05, 1.0, 1.0, 0, 0.25, 0.04),
            "num_swap_payments",
        );
        expect_invalid(
            hull_white.forward_swap_rate_t(0.05, 1.0, 1.5, 0, 0.25),
            "num_swap_payments",
        );
        expect_invalid(
            hull_white.swap_rate_t(0.05, 1.0, 0, 0.25),
            "num_swap_payments",
        );
        expect_invalid(
            hull_white.american_payer_swaption_t(0.05, 1.0, 1.5, 0, 0.25, 0.04, 50),
            "num_swap_payments",
        );
        expect_invalid(
            hull_white.american_receiver_swaption_t(0.05, 1.0, 1.5, 8, 0.25, 0.04, 0),
            "num_steps",
        );
    }

    #[test]
    fn non_finite_input_is_invalid_input_and_is_named() {
        yvf_setup!(hull_white);
        let nan = f64::NAN;
        let inf = f64::INFINITY;
        expect_invalid(hull_white.bond_price_t(nan, 0.0, 2.0), "r_t");
        expect_invalid(hull_white.bond_price_t(0.05, inf, 2.0), "t");
        expect_invalid(hull_white.bond_price_now(nan), "bond_maturity");
        expect_invalid(
            hull_white.swap_price_t(0.05, 1.0, 3.0, 0.25, inf),
            "swap_rate",
        );
        expect_invalid(
            hull_white.coupon_bond_price_t(0.05, 1.0, &[1.5, 2.0], nan),
            "coupon_rate",
        );
        expect_invalid(
            hull_white.american_payer_swaption_t(0.05, 1.0, 1.5, 8, 0.25, nan, 50),
            "swap_rate",
        );
        expect_invalid(hull_white.variance_r(nan, 2.0), "t");
        expect_invalid(hull_white.t_forward_bond_vol(1.0, 2.0, nan), "t_f");
    }

    #[test]
    fn negative_valuation_time_is_invalid_input() {
        yvf_setup!(hull_white);
        expect_invalid(hull_white.bond_price_t(0.05, -1.0, 2.0), "t");
        expect_invalid(hull_white.swap_price_t(0.05, -1.0, 3.0, 0.25, 0.04), "t");
        expect_invalid(
            hull_white.swap_price_t_init(0.05, 1.0, 0.5, 8, 0.25, 0.04),
            "swap_start",
        );
        expect_invalid(
            hull_white.forward_swap_rate_t(0.05, 2.0, 1.5, 8, 0.25),
            "swap_initiation",
        );
    }

    #[test]
    fn bond_option_needs_the_underlying_to_outlive_the_option() {
        yvf_setup!(hull_white);
        //bond maturing at expiry leaves nothing to deliver; before expiry is not deliverable at all
        expect_invalid(
            hull_white.bond_call_t(0.05, 0.5, 1.5, 1.5, 0.98),
            "bond_maturity",
        );
        expect_invalid(hull_white.bond_put_now(1.5, 1.25, 0.98), "bond_maturity");
        //an option with no time left is not an option (zero vol also blows up Black-Scholes)
        expect_invalid(
            hull_white.bond_call_t(0.05, 1.0, 1.0, 2.0, 0.98),
            "option_maturity",
        );
    }

    #[test]
    fn jamshidian_rejects_coupons_paid_before_option_expiry() {
        yvf_setup!(hull_white);
        //The underlying of the decomposed option is the bond *at* expiry.  A coupon paid before
        //expiry is not part of that bond, so pricing it would value a leg that does not exist.
        expect_invalid(
            hull_white.coupon_bond_call_t(0.05, 1.0, 1.5, &[1.25, 1.75, 2.0], 0.05, 1.0),
            "coupon_times[0]",
        );
        expect_invalid(
            hull_white.coupon_bond_put_t(0.05, 1.0, 1.5, &[1.5, 2.0], 0.05, 1.0),
            "coupon_times[0]",
        );
    }

    #[test]
    fn caplet_strike_that_breaks_the_bond_put_transform_is_invalid() {
        yvf_setup!(hull_white);
        //1 + delta * strike == 0 makes the transformed strike 1/0
        expect_invalid(
            hull_white.caplet_t(0.05, 1.0, 1.5, 0.25, -4.0),
            "1 + delta * strike",
        );
        expect_invalid(hull_white.caplet_now(1.5, 0.25, -5.0), "1 + delta * strike");
        //a negative strike deeper than -1/delta is not a caplet
        expect_invalid(
            hull_white.caplet_t(0.05, 1.0, 1.5, 0.25, -40.0),
            "1 + delta * strike",
        );
    }

    #[test]
    fn non_finite_computed_result_surfaces_as_numerical_error() {
        yvf_setup!(hull_white);
        //Inputs are all finite and ordered, but exp() overflows: that is a numerical failure and
        //must not be reported as a price.
        let err = hull_white.bond_price_t(-1e3, 0.0, 10.0).unwrap_err();
        assert!(
            matches!(err, HullWhiteError::NumericalError(_)),
            "expected NumericalError for an overflowing price, got {err:?}"
        );
        assert!(err.to_string().contains("bond_price_t"), "{err}");
    }

    #[test]
    fn valid_instruments_still_price() {
        //Smoke test for the other side of the contract: nothing that used to be valid got caught by
        //the guards.  The 32 pre-existing numerical tests (including the bit-exact swaption pins)
        //are the real regression net; this walks every public entry point once.
        yvf_setup!(hull_white);
        let coupon_times = get_coupon_times(4, 1.0, 0.25).unwrap();
        assert_eq!(coupon_times.len(), 4);
        //bond schedules for the option entry points must sit entirely beyond expiry
        let post_expiry = get_coupon_times(4, 1.5, 0.25).unwrap();

        assert!(
            hull_white
                .t_forward_bond_vol(1.0, 2.0, 3.0)
                .unwrap()
                .is_finite()
        );
        assert!(hull_white.mu_r(0.05, 1.0, 2.0).unwrap().is_finite());
        assert!(hull_white.variance_r(1.0, 2.0).unwrap().is_finite());
        //a bond priced at its own maturity is at par
        assert_abs_diff_eq!(
            hull_white.bond_price_t(0.05, 2.0, 2.0).unwrap(),
            1.0,
            epsilon = 1e-12
        );
        assert!(hull_white.bond_price_now(2.0).unwrap() > 0.0);
        //t == 0 entry points are valid at spot
        assert!(hull_white.bond_price_t(0.05, 0.0, 2.0).unwrap().is_finite());
        assert!(
            hull_white
                .coupon_bond_price_t(0.05, 1.0, &coupon_times, 0.05)
                .unwrap()
                > 0.0
        );
        assert!(
            hull_white
                .coupon_bond_price_now(&coupon_times, 0.05)
                .unwrap()
                > 0.0
        );
        assert!(hull_white.bond_call_t(0.05, 1.0, 1.5, 2.0, 0.98).unwrap() >= 0.0);
        assert!(hull_white.bond_call_now(1.5, 2.0, 0.98).unwrap() >= 0.0);
        assert!(hull_white.bond_put_t(0.05, 1.0, 1.5, 2.0, 0.98).unwrap() >= 0.0);
        assert!(hull_white.bond_put_now(1.5, 2.0, 0.98).unwrap() >= 0.0);
        assert!(
            hull_white
                .coupon_bond_call_t(0.05, 1.0, 1.5, &post_expiry, 0.05, 1.0)
                .unwrap()
                >= 0.0
        );
        //This put is deep out of the money (forward ~0.95 against a strike of 1.0), so the honest
        //answer is ~0; the ~1e-17 residual is floating-point noise and is identical to what the
        //pre-change code returned, so assert the magnitude rather than a non-negative bound.
        assert!(
            hull_white
                .coupon_bond_put_t(0.05, 1.0, 1.5, &post_expiry, 0.05, 1.0)
                .unwrap()
                .abs()
                < 1e-12
        );
        assert!(hull_white.caplet_now(1.5, 0.25, 0.04).unwrap() >= 0.0);
        assert!(hull_white.caplet_t(0.05, 1.0, 1.5, 0.25, 0.04).unwrap() >= 0.0);
        assert!(
            hull_white
                .euro_dollar_future_t(0.05, 1.0, 1.5, 0.25)
                .unwrap()
                .is_finite()
        );
        assert!(
            hull_white
                .euro_dollar_future_now(1.5, 0.25)
                .unwrap()
                .is_finite()
        );
        //spot libor fixing at t is legal (maturity == t)
        assert!(
            hull_white
                .forward_libor_rate_t(0.05, 1.0, 1.0, 0.25)
                .unwrap()
                .is_finite()
        );
        assert!(
            hull_white
                .forward_libor_rate_now(1.5, 0.25)
                .unwrap()
                .is_finite()
        );
        assert!(
            hull_white
                .libor_rate_t(0.05, 1.0, 0.25)
                .unwrap()
                .is_finite()
        );
        assert!(
            hull_white
                .forward_swap_rate_t(0.05, 1.0, 1.5, 8, 0.25)
                .unwrap()
                .is_finite()
        );
        assert!(
            hull_white
                .swap_rate_t(0.05, 1.0, 8, 0.25)
                .unwrap()
                .is_finite()
        );
        assert!(
            hull_white
                .swap_price_t(0.05, 1.0, 3.0, 0.25, 0.04)
                .unwrap()
                .is_finite()
        );
        assert!(
            hull_white
                .swap_price_t_init(0.05, 1.0, 1.0, 8, 0.25, 0.04)
                .unwrap()
                .is_finite()
        );
        assert!(
            hull_white
                .european_payer_swaption_t(0.05, 1.0, 1.5, 8, 0.25, 0.04)
                .unwrap()
                >= 0.0
        );
        assert!(
            hull_white
                .european_receiver_swaption_t(0.05, 1.0, 1.5, 8, 0.25, 0.04)
                .unwrap()
                >= 0.0
        );
        assert!(
            hull_white
                .american_payer_swaption_t(0.05, 1.0, 1.5, 8, 0.25, 0.04, 100)
                .unwrap()
                >= 0.0
        );
        assert!(
            hull_white
                .american_receiver_swaption_t(0.05, 1.0, 1.5, 8, 0.25, 0.04, 100)
                .unwrap()
                >= 0.0
        );
    }
    // ---------------------------------------------------------------------------
    // Jamshidian robustness (workspace-76t)
    // ---------------------------------------------------------------------------

    fn simpson(f: &dyn Fn(f64) -> f64, a: f64, b: f64, panels: usize) -> f64 {
        let panels = panels.max(2) + (panels % 2); //Simpson needs an even panel count
        let h = (b - a) / panels as f64;
        let mut sum = f(a) + f(b);
        for i in 1..panels {
            let x = a + i as f64 * h;
            sum += if i % 2 == 1 { 4.0 } else { 2.0 } * f(x);
        }
        sum * h / 3.0
    }

    /// Mean and standard deviation of the short rate at `option_maturity`, under the
    /// *option-maturity-forward* measure.
    ///
    /// A claim settled on the option's expiry date is priced with the zero-coupon bond maturing
    /// then as numeraire, so that -- not the risk-neutral money-market measure -- is the law the
    /// option's payoff has to be integrated against.  Under this model `r_U` is Gaussian either
    /// way with the same variance, and the drift moves by
    /// `sigma^2 * int_t^U exp(-a(U-s)) B(s,U) ds`, which closes to the form below.
    fn expiry_rate_moments(
        hull_white: &HullWhite<impl Fn(f64) -> f64 + Sync, impl Fn(f64) -> f64 + Sync>,
        r_t: f64,
        t: f64,
        option_maturity: f64,
    ) -> (f64, f64) {
        let tau = option_maturity - t;
        let a = hull_white.a;
        let sigma = hull_white.sigma;
        let risk_neutral_mean = hull_white.mu_r(r_t, t, option_maturity).unwrap();
        let variance = hull_white.variance_r(t, option_maturity).unwrap();
        let numeraire_drift = (sigma * sigma / (a * a))
            * ((1.0 - (-a * tau).exp()) - 0.5 * (1.0 - (-2.0 * a * tau).exp()));
        (risk_neutral_mean - numeraire_drift, variance.sqrt())
    }

    /// Prices the payoff directly by integrating over the expiry-date distribution of the short
    /// rate, as an independent check on the Jamshidian machinery.
    ///
    /// ```text
    /// price = P(t,U) * E^U[(P_c(mu + sd*Z) - strike)^+]
    /// ```
    ///
    /// Nothing in here uses Jamshidian's decomposition, the analytic bracket or the root solver,
    /// so agreement pins those down rather than echoing them.  The payoff kinks where
    /// `P_c(r) = strike`; that kink is found by this test's own plain bisection and each side is
    /// integrated separately, so the kink is never interior to a panel.  `Z` is truncated at
    /// +/-12 sigma, where `phi` underflows, putting the tail error far below any price worth
    /// quoting.
    #[allow(clippy::too_many_arguments)] //a test helper with an instrument's whole description
    fn direct_payoff_price(
        hull_white: &HullWhite<impl Fn(f64) -> f64 + Sync, impl Fn(f64) -> f64 + Sync>,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
        is_call: bool,
    ) -> f64 {
        let coupon_bond_at = |rate: f64| {
            hull_white
                .coupon_bond_price_t(rate, option_maturity, coupon_times, coupon_rate)
                .unwrap()
        };
        let (mean, sd) = expiry_rate_moments(hull_white, r_t, t, option_maturity);
        assert!(
            sd > 0.0,
            "the reference needs a non-degenerate expiry distribution"
        );

        //Kink: P_c(r) = strike.  P_c falls strictly with the rate, so widen until the ends
        //disagree and then halve.
        let g = |rate: f64| coupon_bond_at(rate) - strike;
        let mut lo = mean - sd;
        let mut hi = mean + sd;
        while g(lo) < 0.0 {
            lo = mean - 2.0 * (mean - lo);
            assert!(lo > -1e6, "no kink found below the mean");
        }
        while g(hi) > 0.0 {
            hi = mean + 2.0 * (hi - mean);
            assert!(hi < 1e6, "no kink found above the mean");
        }
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if g(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let kink_z = (0.5 * (lo + hi) - mean) / sd;

        let pdf = |z: f64| (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt();
        let integrand = |z: f64| {
            let density = pdf(z);
            if density == 0.0 {
                //Stops inf * 0 = NaN in a deep tail where the payoff itself overflows.
                return 0.0;
            }
            let value = coupon_bond_at(mean + sd * z);
            let payoff = if is_call {
                (value - strike).max(0.0)
            } else {
                (strike - value).max(0.0)
            };
            payoff * density
        };
        //A call is in the money below the kink, a put above it.
        let (a, b) = if is_call {
            (-12.0, kink_z.min(12.0))
        } else {
            (kink_z.max(-12.0), 12.0)
        };
        if a >= b {
            return 0.0;
        }
        let coarse = simpson(&integrand, a, b, 800);
        let fine = simpson(&integrand, a, b, 1600);
        assert!(
            (coarse - fine).abs() < 1e-8f64.max(fine.abs() * 1e-9),
            "reference quadrature not converged: {coarse} vs {fine}"
        );
        hull_white.bond_price_t(r_t, t, option_maturity).unwrap() * fine
    }

    const FIXTURES: [(f64, f64, f64, f64); 4] = [
        //curr_rate, a, b, sigma
        (0.05, 0.05, 0.05, 0.01),
        (STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG),
        (0.05, 0.1, 0.08, 0.08),
        (0.03, 0.4, 0.045, 0.005),
    ];

    #[test]
    fn the_forward_measure_drift_reproduces_the_bond_price() {
        //Guards the reference itself: repricing the coupon bond through the expiry-date
        //distribution of the rate must return today's bond price exactly.  Without the numeraire
        //drift this is out by ~1e-4, which is exactly the size of error that was showing up in
        //every option comparison below until the measure was fixed.
        let times = [2.5, 3.0, 3.5, 4.0];
        for &(curr, a, b, sigma) in FIXTURES.iter() {
            let (yield_curve, forward_curve) = hw_curves(curr, a, b, sigma);
            let hull_white = HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
            let (r_t, t, u) = (0.04, 1.0, 2.0);
            let bond = hull_white
                .coupon_bond_price_t(r_t, t, &times, 0.05)
                .unwrap();
            let (mean, sd) = expiry_rate_moments(&hull_white, r_t, t, u);
            let expected_value = simpson(
                &|z: f64| {
                    let density = (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt();
                    hull_white
                        .coupon_bond_price_t(mean + sd * z, u, &times, 0.05)
                        .unwrap()
                        * density
                },
                -12.0,
                12.0,
                4000,
            );
            let repriced = hull_white.bond_price_t(r_t, t, u).unwrap() * expected_value;
            assert!(
                (repriced - bond).abs() < 1e-12,
                "fixture {curr},{a},{b},{sigma}: repriced {repriced} vs bond {bond}"
            );
        }
    }

    #[test]
    fn jamshidian_matches_the_direct_payoff_integral() {
        //The headline check: the decomposition, its bracket and its solver, against a quadrature
        //over the payoff that shares none of them.  The worst deviation across this whole grid is
        //about 2e-12.
        let times = [2.5, 3.0, 3.5, 4.0];
        for &(curr, a, b, sigma) in FIXTURES.iter() {
            let (yield_curve, forward_curve) = hw_curves(curr, a, b, sigma);
            let hull_white = HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
            for strike in [0.5f64, 0.8, 0.95, 1.0, 1.05, 1.3, 3.0] {
                for is_call in [true, false] {
                    let priced = if is_call {
                        hull_white
                            .coupon_bond_call_t(0.04, 1.0, 2.0, &times, 0.05, strike)
                            .unwrap()
                    } else {
                        hull_white
                            .coupon_bond_put_t(0.04, 1.0, 2.0, &times, 0.05, strike)
                            .unwrap()
                    };
                    let reference = direct_payoff_price(
                        &hull_white,
                        0.04,
                        1.0,
                        2.0,
                        &times,
                        0.05,
                        strike,
                        is_call,
                    );
                    assert!(
                        (priced - reference).abs() <= 1e-10f64.max(reference.abs() * 1e-9),
                        "fixture {curr},{a},{b},{sigma} {side} strike {strike}: {priced} vs {reference}",
                        side = if is_call { "call" } else { "put" }
                    );
                }
            }
        }
    }

    #[test]
    fn a_fine_strike_ladder_around_the_money_matches_the_integral() {
        //At the money the price is most sensitive to the critical rate, so a fine ladder through
        //the ATM region is where a sloppy root shows up first.
        let times = [2.5, 3.0, 3.5, 4.0];
        for &(curr, a, b, sigma) in FIXTURES.iter() {
            let (yield_curve, forward_curve) = hw_curves(curr, a, b, sigma);
            let hull_white = HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
            let (r_t, t, u) = (0.04, 1.0, 2.0);
            let underlying = hull_white
                .coupon_bond_price_t(r_t, t, &times, 0.05)
                .unwrap();
            for step in 0..21 {
                let strike = underlying * (0.9 + step as f64 * 0.01);
                for is_call in [true, false] {
                    let priced = if is_call {
                        hull_white
                            .coupon_bond_call_t(r_t, t, u, &times, 0.05, strike)
                            .unwrap()
                    } else {
                        hull_white
                            .coupon_bond_put_t(r_t, t, u, &times, 0.05, strike)
                            .unwrap()
                    };
                    let reference =
                        direct_payoff_price(&hull_white, r_t, t, u, &times, 0.05, strike, is_call);
                    assert!(
                        (priced - reference).abs() <= 1e-10f64.max(reference.abs() * 1e-9),
                        "fixture {curr},{a},{b},{sigma} {side} strike {strike}: {priced} vs {reference}",
                        side = if is_call { "call" } else { "put" }
                    );
                }
            }
        }
    }

    #[test]
    fn a_negative_coupon_schedule_still_prices() {
        //A coupon below zero is not nonsense -- deep-discount instruments and some cross-currency
        //legs pay one -- and nothing in the decomposition needs the coupon to be positive as
        //long as the weights keep their sign.  The reference checks it independently.
        let times = [2.5, 3.0, 3.5, 4.0];
        let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
        let hull_white = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve).unwrap();
        let (r_t, t, u) = (0.04, 1.0, 2.0);
        for coupon_rate in [-0.02f64, -0.005, -0.0001] {
            for strike in [0.5f64, 0.9, 0.95, 1.0, 1.2] {
                for is_call in [true, false] {
                    let priced = if is_call {
                        hull_white
                            .coupon_bond_call_t(r_t, t, u, &times, coupon_rate, strike)
                            .unwrap()
                    } else {
                        hull_white
                            .coupon_bond_put_t(r_t, t, u, &times, coupon_rate, strike)
                            .unwrap()
                    };
                    let reference = direct_payoff_price(
                        &hull_white,
                        r_t,
                        t,
                        u,
                        &times,
                        coupon_rate,
                        strike,
                        is_call,
                    );
                    assert!(
                        (priced - reference).abs() <= 1e-10f64.max(reference.abs() * 1e-9),
                        "coupon {coupon_rate} {side} strike {strike}: {priced} vs {reference}",
                        side = if is_call { "call" } else { "put" }
                    );
                }
            }
        }
    }

    #[test]
    fn a_zero_strike_call_is_the_underlying_and_a_zero_strike_put_is_worthless() {
        //A zero strike puts the critical rate at the far end of the bond's range, which used to
        //be a guaranteed RootFindingError.  What the answer has to be: the call is the present
        //value of the underlying, the put is nothing.
        let times = [2.5, 3.0, 3.5, 4.0];
        for &(curr, a, b, sigma) in FIXTURES.iter() {
            let (yield_curve, forward_curve) = hw_curves(curr, a, b, sigma);
            let hull_white = HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
            let underlying = hull_white
                .coupon_bond_price_t(0.04, 1.0, &times, 0.05)
                .unwrap();
            let call = hull_white
                .coupon_bond_call_t(0.04, 1.0, 2.0, &times, 0.05, 0.0)
                .unwrap();
            let put = hull_white
                .coupon_bond_put_t(0.04, 1.0, 2.0, &times, 0.05, 0.0)
                .unwrap();
            assert!(
                (call - underlying).abs() <= underlying.abs() * 1e-12,
                "fixture {curr},{a},{b},{sigma}: call {call} vs underlying {underlying}"
            );
            assert_eq!(put, 0.0, "fixture {curr},{a},{b},{sigma}: put {put}");
        }
    }

    #[test]
    fn put_call_parity_holds_whatever_the_strike() {
        //C - P = PV(underlying) - K * P(t,U), to machine precision, including at zero strike.
        //A solve that lands on the wrong critical rate breaks this immediately.
        let times = [2.5, 3.0, 3.5, 4.0];
        let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
        let hull_white = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve).unwrap();
        let (r_t, t, u) = (0.04, 1.0, 2.0);
        let underlying = hull_white
            .coupon_bond_price_t(r_t, t, &times, 0.05)
            .unwrap();
        let discount = hull_white.bond_price_t(r_t, t, u).unwrap();
        for strike in [0.0f64, 0.3, 0.7, 0.95, 1.0, 1.2, 2.0, 1e3, 1e6] {
            let call = hull_white
                .coupon_bond_call_t(r_t, t, u, &times, 0.05, strike)
                .unwrap();
            let put = hull_white
                .coupon_bond_put_t(r_t, t, u, &times, 0.05, strike)
                .unwrap();
            let parity = underlying - strike * discount;
            //At a strike of 1e6 the put is a difference of huge cancelling leg values, so the
            //residual is a few parts in 1e11 of the numbers involved rather of the option.
            assert!(
                (call - put - parity).abs() <= 1e-11f64.max(parity.abs() * 1e-10),
                "strike {strike}: C-P {} vs parity {parity}",
                call - put
            );
        }
    }

    #[test]
    fn an_extreme_strike_prices_rather_than_erroring() {
        //A strike orders of magnitude away from anywhere the bond can reach has no optionality
        //left: the answer is zero (or parity), not an error.  This used to come back as
        //RootFindingError("NaN") because the old Newton iterate blew up on the flat of the
        //exponential.
        let times = [2.5, 3.0, 3.5, 4.0];
        let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
        let hull_white = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve).unwrap();
        for strike in [1e3f64, 1e8, 1e12] {
            let call = hull_white
                .coupon_bond_call_t(0.04, 1.0, 2.0, &times, 0.05, strike)
                .unwrap();
            assert_eq!(call, 0.0, "strike {strike} call {call}");
            let put = hull_white
                .coupon_bond_put_t(0.04, 1.0, 2.0, &times, 0.05, strike)
                .unwrap();
            assert!(put.is_finite() && put > 0.0, "strike {strike} put {put}");
        }
    }

    #[test]
    fn the_price_does_not_depend_on_where_the_solver_starts() {
        //The old solver started from a hard-coded 3% and its answer moved with that guess.  With
        //a real bracket and a bisection guard, any seed reaches the same root.
        let times = [2.5, 3.0, 3.5, 4.0];
        let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
        let reference = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve)
            .unwrap()
            .coupon_bond_call_t(0.04, 1.0, 2.0, &times, 0.05, 0.95)
            .unwrap();
        for seed in [-1.0f64, 0.0, 0.03, 0.5, 5.0, 50.0, 1e3, 1e6] {
            let hull_white = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve)
                .unwrap()
                .with_solver(SolverSettings {
                    initial_guess: Some(seed),
                    ..SolverSettings::default()
                })
                .unwrap();
            let priced = hull_white
                .coupon_bond_call_t(0.04, 1.0, 2.0, &times, 0.05, 0.95)
                .unwrap();
            assert!(
                (priced - reference).abs() <= 1e-10f64.max(reference.abs() * 1e-9),
                "seed {seed}: {priced} vs {reference}"
            );
        }
    }

    #[test]
    fn the_bracket_actually_brackets() {
        //The analytic bracket has to straddle the critical rate: bond value above the strike at
        //the low end, below it at the high end, and the solved root in between.
        let times = [2.5, 3.0, 3.5, 4.0];
        for &(curr, a, b, sigma) in FIXTURES.iter() {
            let (yield_curve, forward_curve) = hw_curves(curr, a, b, sigma);
            let hull_white = HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
            let (r_t, t, u) = (0.04, 1.0, 2.0);
            for strike in [0.5f64, 0.8, 0.95, 1.0, 1.05, 1.3, 3.0] {
                let (lower, upper) = hull_white
                    .critical_rate_bracket(u, &times, 0.05, strike)
                    .unwrap_or_else(|| {
                        panic!("no bracket for {curr},{a},{b},{sigma} strike {strike}")
                    });
                let gap = |rate: f64| {
                    hull_white
                        .coupon_bond_price_t(rate, u, &times, 0.05)
                        .unwrap()
                        - strike
                };
                assert!(
                    gap(lower) >= 0.0,
                    "{curr},{a},{b},{sigma} strike {strike}: low end {lower} gives {}",
                    gap(lower)
                );
                assert!(
                    gap(upper) <= 0.0,
                    "{curr},{a},{b},{sigma} strike {strike}: high end {upper} gives {}",
                    gap(upper)
                );
                let solution = hull_white
                    .solve_critical_rate(r_t, t, u, &times, 0.05, strike)
                    .unwrap();
                assert!(
                    solution.root >= lower && solution.root <= upper,
                    "{curr},{a},{b},{sigma} strike {strike}: root {} outside [{lower}, {upper}]",
                    solution.root
                );
            }
        }
    }

    #[test]
    fn a_single_coupon_schedule_reduces_to_the_zero_coupon_option() {
        //The degenerate schedule: one payment date, carrying coupon plus redemption.  Jamshidian
        //has exactly one leg, so the answer has to be that leg weighted by (1 + coupon) -- which
        //is the plain zero-coupon bond option, priced here without any root find at all.
        //
        //  (1 + c) * call_discount(P(t,T), K / (1 + c), P(t,U), sigma_leg)
        let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
        let hull_white = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve).unwrap();
        let (r_t, t, u) = (0.04, 1.0, 2.0);
        for bond_maturity in [2.5f64, 4.0, 10.0] {
            for coupon_rate in [0.05f64, 0.01, -0.02] {
                let weight = 1.0 + coupon_rate;
                for strike in [0.5f64, 0.9, 0.95, 1.0, 1.3, 3.0] {
                    for is_call in [true, false] {
                        let schedule = [bond_maturity];
                        let priced = if is_call {
                            hull_white
                                .coupon_bond_call_t(r_t, t, u, &schedule, coupon_rate, strike)
                                .unwrap()
                        } else {
                            hull_white
                                .coupon_bond_put_t(r_t, t, u, &schedule, coupon_rate, strike)
                                .unwrap()
                        };
                        //The same option through the zero-coupon pricer: one leg, struck at the
                        //strike the single leg would carry.  No decomposition, no solver.
                        let leg = if is_call {
                            hull_white
                                .bond_call_t(r_t, t, u, bond_maturity, strike / weight)
                                .unwrap()
                        } else {
                            hull_white
                                .bond_put_t(r_t, t, u, bond_maturity, strike / weight)
                                .unwrap()
                        };
                        let reference = weight * leg;
                        assert!(
                            (priced - reference).abs() <= 1e-12f64.max(reference.abs() * 1e-11),
                            "T {bond_maturity} c {coupon_rate} {} strike {strike}: {priced} vs {reference}",
                            if is_call { "call" } else { "put" }
                        );
                        //...and the direct integral agrees too.
                        let integral = direct_payoff_price(
                            &hull_white,
                            r_t,
                            t,
                            u,
                            &schedule,
                            coupon_rate,
                            strike,
                            is_call,
                        );
                        assert!(
                            (priced - integral).abs() <= 1e-10f64.max(integral.abs() * 1e-9),
                            "T {bond_maturity} c {coupon_rate} {} strike {strike}: {priced} vs integral {integral}",
                            if is_call { "call" } else { "put" }
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn a_forty_eight_coupon_schedule_prices_against_the_integral() {
        //A long schedule is where the solve travels furthest and where the leg sum has the most
        //terms to get wrong: 48 coupons, maturities out to 26 years, strikes from deep ITM to
        //deep OTM.  Still agrees with the direct payoff integral.
        let schedule: Vec<f64> = (1..=48).map(|i| 2.0 + 0.5 * i as f64).collect();
        assert_eq!(schedule.len(), 48);
        for &(curr, a, b, sigma) in FIXTURES.iter() {
            let (yield_curve, forward_curve) = hw_curves(curr, a, b, sigma);
            let hull_white = HullWhite::init(a, sigma, &yield_curve, &forward_curve).unwrap();
            let (r_t, t, u) = (0.04, 1.0, 2.0);
            let underlying = hull_white
                .coupon_bond_price_t(r_t, t, &schedule, 0.05)
                .unwrap();
            for factor in [0.2f64, 0.6, 0.9, 0.99, 1.0, 1.01, 1.4, 2.5] {
                let strike = underlying * factor;
                for is_call in [true, false] {
                    let priced = if is_call {
                        hull_white
                            .coupon_bond_call_t(r_t, t, u, &schedule, 0.05, strike)
                            .unwrap()
                    } else {
                        hull_white
                            .coupon_bond_put_t(r_t, t, u, &schedule, 0.05, strike)
                            .unwrap()
                    };
                    let reference = direct_payoff_price(
                        &hull_white,
                        r_t,
                        t,
                        u,
                        &schedule,
                        0.05,
                        strike,
                        is_call,
                    );
                    assert!(
                        (priced - reference).abs() <= 1e-9f64.max(reference.abs() * 1e-9),
                        "fixture {curr},{a},{b},{sigma} {side} strike {strike}: {priced} vs {reference}",
                        side = if is_call { "call" } else { "put" }
                    );
                    //Parity on the same instrument, for a second independent angle.
                    let other = if is_call {
                        hull_white
                            .coupon_bond_put_t(r_t, t, u, &schedule, 0.05, strike)
                            .unwrap()
                    } else {
                        hull_white
                            .coupon_bond_call_t(r_t, t, u, &schedule, 0.05, strike)
                            .unwrap()
                    };
                    let parity = underlying - strike * hull_white.bond_price_t(r_t, t, u).unwrap();
                    //C - P = parity, so a put's difference runs the other way.
                    let expected = if is_call { parity } else { -parity };
                    assert!(
                        (priced - other - expected).abs() <= 1e-9f64.max(parity.abs() * 1e-10),
                        "fixture {curr},{a},{b},{sigma} strike {strike}: {} vs {}",
                        priced - other,
                        expected
                    );
                }
            }
        }
    }

    #[test]
    fn a_starved_solver_says_which_stage_failed() {
        //A failure names the stage that failed instead of handing back a bare number: too few
        //iterations to converge, in the context of the instrument being priced.
        let times = [2.5, 3.0, 3.5, 4.0];
        let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
        let starved = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve)
            .unwrap()
            .with_solver(SolverSettings {
                max_iterations: 1,
                ..SolverSettings::default()
            })
            .unwrap();
        let error = starved
            .coupon_bond_call_t(0.04, 1.0, 2.0, &times, 0.05, 0.95)
            .unwrap_err();
        let message = error.to_string();
        assert!(
            message.contains("iterations"),
            "expected an iteration failure, got {message}"
        );
        assert!(
            message.contains("Jamshidian"),
            "expected instrument context, got {message}"
        );
    }

    /// Worst price error over a few strikes, against a 1e-15-tolerance reference, for one solve
    /// tolerance.
    fn worst_price_err(
        tolerance: f64,
        yield_curve: &(impl Fn(f64) -> f64 + Sync),
        forward_curve: &(impl Fn(f64) -> f64 + Sync),
        times: &[f64],
    ) -> f64 {
        let reference_model = HullWhite::init(STEEP_A, STEEP_SIG, yield_curve, forward_curve)
            .unwrap()
            .with_solver(SolverSettings {
                tolerance: 1e-15,
                max_iterations: 500,
                initial_guess: None,
            })
            .unwrap();
        let model = HullWhite::init(STEEP_A, STEEP_SIG, yield_curve, forward_curve)
            .unwrap()
            .with_solver(SolverSettings {
                tolerance,
                ..SolverSettings::default()
            })
            .unwrap();
        [0.8f64, 0.95, 1.0, 1.05]
            .into_iter()
            .map(|strike| {
                let reference = direct_payoff_price(
                    &reference_model,
                    0.04,
                    1.0,
                    2.0,
                    times,
                    0.05,
                    strike,
                    true,
                );
                let price = model
                    .coupon_bond_call_t(0.04, 1.0, 2.0, times, 0.05, strike)
                    .unwrap();
                (price - reference).abs()
            })
            .fold(0.0f64, f64::max)
    }

    #[test]
    fn a_looser_root_tolerance_is_visible_in_the_price() {
        //The reason the tolerance is configurable at all: the old hard-coded 1e-7 capped how
        //accurate a price anyone could get, and now that cap is measurable instead of invisible.
        let times = [2.5, 3.0, 3.5, 4.0];
        let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
        let loose = worst_price_err(1e-4, &yield_curve, &forward_curve, &times);
        let tight = worst_price_err(1e-14, &yield_curve, &forward_curve, &times);
        println!("tolerance 1e-4 worst price err {loose:e}; 1e-14 worst price err {tight:e}");
        assert!(tight < loose, "tight {tight:e} should beat loose {loose:e}");
        assert!(tight < 1e-11, "tight tolerance left {tight:e}");
    }

    #[test]
    fn the_root_tolerance_bounds_the_search_not_the_achieved_error() {
        //The tolerance is a statement about how hard the solver keeps hunting (bracket width in rate
        //units), not about how far off the answer it returns.  When the solve exits on the Newton
        //correction the returned root is Newton's *prediction* of the root, which is good to far
        //better than the correction's size, so the old 1e-7 cap is no longer a cap on accuracy:
        //1e-7 now lands on the same price as 1e-14 (both sit on the model's round-off floor).
        let times = [2.5, 3.0, 3.5, 4.0];
        let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
        let legacy_cap = worst_price_err(1e-7, &yield_curve, &forward_curve, &times);
        let tight = worst_price_err(1e-14, &yield_curve, &forward_curve, &times);
        println!("tolerance 1e-7 worst price err {legacy_cap:e}; 1e-14 worst price err {tight:e}");
        //A 1e-7 rate error on this instrument is worth ~1e-9 of price; the achieved error is three
        //orders below that, so the legacy tolerance no longer costs accuracy -- only iterations.
        assert!(
            legacy_cap < 1e-11 && legacy_cap <= tight * 1.5 + 1e-15,
            "tolerance 1e-7 left {legacy_cap:e}, tight 1e-14 left {tight:e}"
        );
    }
}
