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
//! requiring the model to fit the initial term structure the [`YieldCurve`] it was built from
//! describes.

use crate::curves::{CurveRef, YieldCurve, validate_curve};
use crate::error::HullWhiteError;
use crate::rootfinder::SolverSettings;
use crate::validation;

/// A one-factor Hull-White model, calibrated to one initial term structure.
///
/// Built with [`HullWhite::new`] from a single [`YieldCurve`], which supplies the cumulative
/// yield the model discounts with and the instantaneous forward it mean-reverts towards; see
/// [`crate::curves`] for the trait, the invariant between those two, and the builders that make
/// one out of a closure ([`from_yield`](crate::from_yield)) or a closure pair
/// ([`from_yield_and_forward`](crate::from_yield_and_forward)).  The curve is borrowed, not
/// owned, so a model is only good for as long as the curve it was calibrated to.
///
/// The struct carries the two model parameters (`a`, `sigma`), the curve, and the root-finding
/// budget used by the Jamshidian solve; see [`HullWhite::with_solver`].  It has no closure type
/// parameters: `HullWhite<'a>` is the whole type, so a function that takes a model says
/// `&HullWhite` rather than naming two curve closures it does not care about.
pub struct HullWhite<'a> {
    pub(crate) a: f64,
    pub(crate) sigma: f64,
    //The initial term structure: cumulative yield and instantaneous forward, as one object.
    //`zero_yield` is not divided by time, so it gets perpetually larger (unless rates are negative).
    curve: CurveRef<'a>,
    /// Tolerance / iteration budget for the Jamshidian critical-rate solve.  [`HullWhite::new`]
    /// uses [`SolverSettings::default`]; [`HullWhite::with_solver`] changes it.
    pub(crate) solver: SolverSettings,
}

/// A `Debug` that shows the model's own numbers.
///
/// The curve is a `dyn` object — a trait object has no `Debug` of its own, and printing the
/// term structure would mean printing a function.  What is worth seeing in a failure message is
/// `a`, `sigma` and the solver budget, so those are printed and the curve is named rather than
/// dumped.
impl core::fmt::Debug for HullWhite<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("HullWhite")
            .field("a", &self.a)
            .field("sigma", &self.sigma)
            .field("solver", &self.solver)
            .field("curve", &"(dyn YieldCurve)")
            .finish()
    }
}

impl<'a> HullWhite<'a> {
    /// Calibrate a model: a mean-reversion speed, a volatility, and one term structure.
    ///
    /// The two parameters are the model's own dynamics (`dr = [theta(t) - a r] dt + sigma dW`);
    /// the curve is what pins `theta(t)`, so a `HullWhite` is meaningless without one and the
    /// calibration is a constructor argument rather than a later setter.
    ///
    /// `a` and `sigma` are checked for positivity and finiteness, and the curve is checked for
    /// internal consistency: `forward(t)` must equal `d/dt zero_yield(t)` at every time in
    /// [`CURVE_PROBE_TIMES`](crate::CURVE_PROBE_TIMES), to within
    /// [`forward_consistency_tolerance`](crate::forward_consistency_tolerance).  A curve that
    /// fails is [`HullWhiteError::InvalidInput`] naming the worst probe time, the supplied
    /// forward, the derived one and the tolerance — because a model built on a cumulative yield
    /// and a forward that disagree does not misprice loudly, it misprices quietly, and the
    /// disagreement is far cheaper to report here than to hunt down in a book.
    ///
    /// # Examples
    ///
    /// ```
    /// use hull_white::{HullWhite, from_yield_and_forward};
    ///
    /// // Cumulative yield 0.05 t + 0.01 t^2, whose derivative is the forward 0.05 + 0.02 t.
    /// let curve = from_yield_and_forward(
    ///     |t: f64| 0.05 * t + 0.01 * t * t,
    ///     |t: f64| 0.05 + 0.02 * t,
    /// );
    /// let hull_white = HullWhite::new(0.15, 0.02, &curve).unwrap();
    /// assert!(hull_white.bond_price_now(2.0).unwrap() > 0.0);
    /// ```
    pub fn new(a: f64, sigma: f64, curve: &'a dyn YieldCurve) -> Result<Self, HullWhiteError> {
        Self::validate_parameters(a, sigma)?;
        validate_curve(curve)?;
        Ok(Self {
            a,
            sigma,
            curve: CurveRef::Borrowed(curve),
            solver: SolverSettings::default(),
        })
    }

    /// The initial term structure this model was calibrated to, as a borrow.
    ///
    /// Read-only access to the term structure the prices come off: `discount(t)` for `P(0,t)`,
    /// `zero_yield(t)` for the cumulative yield, `forward(t)` for the instantaneous forward
    /// `phi(t)` is built on.
    ///
    /// ```
    /// use hull_white::{HullWhite, from_yield};
    ///
    /// let curve = from_yield(|t: f64| 0.05 * t);
    /// let hull_white = HullWhite::new(0.2, 0.03, &curve).unwrap();
    /// let price = hull_white.bond_price_now(2.0).unwrap();
    /// assert!((price - hull_white.curve().discount(2.0)).abs() < 1e-15);
    /// ```
    #[must_use = "curve() is a borrow of the model's term structure"]
    pub fn curve(&self) -> &dyn YieldCurve {
        self.curve.get()
    }

    /// Calibrate a model from a cumulative-yield closure and a forward closure.
    ///
    /// Deprecated: those two closures are one term structure, and passing them separately is what
    /// let a caller pass two that disagreed.  Build a [`YieldCurve`] instead —
    /// [`from_yield`](crate::from_yield) when only the cumulative yield is at hand (the forward
    /// is then derived), [`from_yield_and_forward`](crate::from_yield_and_forward) when both are
    /// known in closed form, or your own implementor — and pass it to [`HullWhite::new`]:
    ///
    /// ```text
    /// // 0.8
    /// HullWhite::init(a, sigma, &yield_curve, &forward_curve)
    /// // 0.9
    /// HullWhite::new(a, sigma, &from_yield_and_forward(yield_curve, forward_curve))
    /// ```
    ///
    /// The old call still compiles, and now runs the same consistency check: a pair that disagrees
    /// by more than [`forward_consistency_tolerance`](crate::forward_consistency_tolerance) is an
    /// error rather than a mispricing.  Note that this is a *behaviour* change as well as a type
    /// change — pairs that used to price (wrongly) are rejected; the replacement for "I only have
    /// the yield curve" is [`from_yield`](crate::from_yield).
    ///
    /// Removal: this ships deprecated through the whole of `0.9.x` and is removed in `0.10.0`.
    /// The version bump for the breaking part landed with `0.9.0` (the type collapse); the
    /// destructor step is deliberately one release behind so a caller can move to the new
    /// constructor and still compile against the old one while doing it.
    #[deprecated(
        since = "0.9.0",
        note = "two closures for one term structure; use HullWhite::new with a YieldCurve (from_yield / from_yield_and_forward). Removed in 0.10.0"
    )]
    pub fn init<F, G>(
        a: f64,
        sigma: f64,
        yield_curve: &'a F,
        forward_curve: &'a G,
    ) -> Result<Self, HullWhiteError>
    where
        F: Fn(f64) -> f64 + std::marker::Sync,
        G: Fn(f64) -> f64 + std::marker::Sync,
    {
        Self::validate_parameters(a, sigma)?;
        //`&F` is itself `Fn + Sync` when `F` is, so the wrapper borrows the caller's closures
        //exactly as the 0.8 struct did; nothing is cloned and nothing is copied out.
        let curve = Box::new(crate::curves::from_yield_and_forward(
            yield_curve,
            forward_curve,
        ));
        validate_curve(curve.as_ref())?;
        Ok(Self {
            a,
            sigma,
            curve: CurveRef::Owned(curve),
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = HullWhite::new(0.2, 0.3, &curve).unwrap();
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
    /// Volatility of the `t_m`-maturity bond price under the `t_f`-forward measure.
    ///
    /// This is the Black volatility of the underlying of a bond option: valued at `t`, the bond that
    /// matures at `t_m` is the asset, the option expires at `t_f`, and pricing that option under the
    /// `t_f`-maturity numeraire makes the bond price a lognormal martingale whose sigma is the
    /// number returned here.  It is what feeds the Black / Jamshidian leg strike, not a moment of
    /// the short rate — for those see [`HullWhite::mu_r`] (mean) and [`HullWhite::variance_r`]
    /// (variance).
    ///
    /// The two time gaps do different things to it: `t_m - t` sets how volatile the bond price is
    /// while the option is alive, and `t_f - t_m` scales that down as the option expiry approaches
    /// the bond's own maturity — a bond that matures on (or before) the expiry date has a known
    /// price at expiry, so its forward-measure volatility goes to zero there.  Hence both
    /// `t_m > t` and `t_f > t_m` are required: outside them the expression is not a volatility at
    /// all, and it is refused rather than returned as a negative or zero number that Black will
    /// still happily swallow.
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //valuation date, measured from "now" (0)
    /// let t_m = 2.0; //maturity of the bond that is the option's underlying
    /// let t_f = 3.0; //option expiry == the numeraire's maturity
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let bond_vol = hull_white
    ///     .t_forward_bond_vol(t, t_m, t_f)
    ///     .unwrap();
    /// assert!(bond_vol > 0.0, "a live option has positive bond vol: {bond_vol}");
    /// //Pushing the bond maturity towards the expiry shrinks the vol towards zero: at expiry the
    /// //deliverable's price is already settled, so there is nothing left to be uncertain about.
    /// let near_expiry = hull_white.t_forward_bond_vol(t, 2.9, 3.0).unwrap();
    /// assert!(near_expiry < bond_vol, "{near_expiry} vs {bond_vol}");
    /// //The numeraire has to be the later date: a bond maturing after the option is not priced by
    /// //this measure at all, and asking for it is an error, not a negative volatility.
    /// assert!(hull_white.t_forward_bond_vol(t, 3.0, 2.0).is_err());
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
        self.curve().forward(t) + (self.sigma * exp_t).powi(2) / (2.0 * self.a.powi(2))
    }
    /// Conditional mean of the short rate: `E[r(t_m) | r(t) = r_t]`.
    ///
    /// A Hull-White short rate is Gaussian, so its conditional distribution is fully described by
    /// this mean and [`HullWhite::variance_r`], and the mean is an exponential pull of the observed
    /// rate `r_t` towards the model's own mean-reversion target `phi`:
    ///
    /// ```text
    /// E[r(t_m) | r(t) = r_t] = phi(t_m) + (r_t - phi(t)) * exp(-a * (t_m - t))
    /// ```
    ///
    /// `phi` is the deterministic function the rate reverts towards, fixed by requiring the model
    /// to fit the initial term structure (the same `phi(t)` the crate docs describe, computed by
    /// the crate-internal `phi_t`).  Nothing about this number is a volatility: the spread around
    /// this mean is [`HullWhite::variance_r`], and the *bond*-under-forward-measure volatility is
    /// [`HullWhite::t_forward_bond_vol`].  (This function's doc comment used to read exactly like
    /// that one's, copied and never changed; if you arrived here looking for a vol, that is the
    /// other method.)
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) that r_t is observed at
    /// let t_m = 2.0; //horizon the mean is taken at
    /// let r_t = 0.04; //observed short rate at t
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let mean_rate = hull_white.mu_r(r_t, t, t_m).unwrap();
    /// //No time to move: the horizon is the observation date, so the mean is the observed rate.
    /// assert!((hull_white.mu_r(r_t, t, t).unwrap() - r_t).abs() < 1e-12);
    /// //Here the long-run level is above r_t, so time pulls the mean up from r_t ...
    /// assert!(mean_rate > r_t, "{mean_rate} vs {r_t}");
    /// // ... by less than the full distance, because mean reversion is exponential, not instant.
    /// let further = hull_white.mu_r(r_t, t, 6.0).unwrap();
    /// assert!(further > mean_rate, "{further} vs {mean_rate}");
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
    /// The model's short rate right now, `r(0)`.
    ///
    /// Every `*_now` pricing function that needs a state variable takes it from here rather than
    /// from an argument: at the valuation date `t = 0` there is nothing to condition on except the
    /// initial curve the model was calibrated to, so `r(0)` is not an input, it is a fact about
    /// the calibration.  Concretely `r(0) = phi(0)`, and the volatility term of `phi`
    /// (`sigma^2 (1 - e^{-a t})^2 / (2 a^2)`) vanishes at `t = 0`, so
    ///
    /// ```text
    /// r(0) = phi(0) = curve.forward(0.0)
    /// ```
    ///
    /// i.e. the *instantaneous* forward rate at the front of the curve.  Note which half of the
    /// curve that is: `YieldCurve::forward` is the instantaneous forward `f(0, t)`, while
    /// `YieldCurve::zero_yield` is cumulative (`zero_yield(T)` is the integral of `f(0, .)` over
    /// `[0, T]`, which is why `bond_price_now` is `exp(-zero_yield(T))`).  So `r(0)` comes off
    /// the forward, *not* off `zero_yield(0.0)`, which is `0` for any curve that is an integral.
    ///
    /// This is the value to hand a `*_t` function as `r_t` alongside `t = 0.0` when the state has
    /// to be passed explicitly; every `*_now` variant in this crate is exactly that call, and the
    /// two routes agree to machine precision because the same `r(0)` is what makes
    /// `bond_price_t(r(0), 0, T) == bond_price_now(T)`.
    ///
    /// A curve that is not finite at `0` — a `|t| t.ln()` forward, say — cannot give a short rate
    /// and is an error rather than a `-inf` price.  (It is also a curve the construction check
    /// tolerates, deliberately: it is consistent with its own yield, it just has no front.  See
    /// [`CURVE_PROBE_TIMES`](crate::CURVE_PROBE_TIMES).)
    ///
    /// # Examples
    ///
    /// ```
    /// use hull_white::{HullWhite, YieldCurve, from_yield_and_forward};
    ///
    /// // Cumulative yield: 0.05*T + 0.01*T^2 is the integral of 0.05 + 0.02*t.
    /// let curve = from_yield_and_forward(
    ///     |t: f64| 0.05 * t + 0.01 * t * t,
    ///     // Instantaneous forward; finite at 0, as `r(0)` requires.
    ///     |t: f64| 0.05 + 0.02 * t,
    /// );
    /// let hull_white = HullWhite::new(0.15, 0.02, &curve).unwrap();
    ///
    /// let r0 = hull_white.short_rate_now().unwrap();
    /// assert!((r0 - 0.05).abs() < 1e-12, "r(0) = {r0}");
    /// assert!((r0 - curve.forward(0.0)).abs() < 1e-12, "r(0) is the forward curve at 0");
    ///
    /// // Feeding that rate back into the `t` form at `t = 0` reproduces the `now` bond price.
    /// for maturity in [1.0, 3.0, 7.5] {
    ///     let via_t = hull_white.bond_price_t(r0, 0.0, maturity).unwrap();
    ///     let via_now = hull_white.bond_price_now(maturity).unwrap();
    ///     assert!((via_t - via_now).abs() < 1e-12, "maturity {maturity}: {via_t} vs {via_now}");
    /// }
    /// ```
    pub fn short_rate_now(&self) -> Result<f64, HullWhiteError> {
        validation::finish("short_rate_now", self.short_rate_now_raw())
    }
    /// Unvalidated [`HullWhite::short_rate_now`], for internals that have already checked the
    /// curve or do not need to propagate its failure.
    pub(crate) fn short_rate_now_raw(&self) -> f64 {
        self.phi_t(0.0)
    }
}

#[cfg(test)]
mod tests;
