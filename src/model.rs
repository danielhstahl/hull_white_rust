//! The model itself: [`HullWhite`], its construction and the primitive short-rate moments.
//!
//! Everything else in the crate is built from these four things: the calibrated drift `phi(t)`,
//! the conditional mean `mu_r`, the conditional variance `variance_r` and the zero-coupon bond
//! volatility under the forward measure (`t_forward_bond_vol`).
//!
//! Time is always measured from "now" (0).  A path runs `0 -> t -> T`: `t` is the date the rate
//! `r_t` is observed at, `T` (written `t_m` here) the date a moment is taken at, and instruments
//! add their own coordinates on top (`option_maturity`, `bond_maturity`, `swap_start`, ...).
//! `phi(t)` is the deterministic function the short rate mean-reverts towards; it is fixed by
//! requiring the model to fit the initial term structure given by the yield and forward curves.

use crate::error::HullWhiteError;
use crate::rootfinder::SolverSettings;
use crate::validation;

pub struct HullWhite<'a, T, U>
where
    T: Fn(f64) -> f64 + std::marker::Sync,
    U: Fn(f64) -> f64 + std::marker::Sync,
{
    pub(crate) a: f64,
    pub(crate) sigma: f64,
    //yield_curve is not divided by time, so this gets perpetually larger (unless rates are negative)
    pub(crate) yield_curve: &'a T,
    pub(crate) forward_curve: &'a U,
    /// Tolerance / iteration budget for the Jamshidian critical-rate solve.  [`HullWhite::init`]
    /// uses [`SolverSettings::default`]; [`HullWhite::with_solver`] changes it.
    pub(crate) solver: SolverSettings,
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
    pub(crate) fn phi_t(&self, t: f64) -> f64 {
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
}

#[cfg(test)]
mod tests;
