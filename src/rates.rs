//! Simple-rate instruments: caplets, caps and floors, Eurodollar futures and forward spot Libor.
//!
//! These are all translations of bond prices into a simple (non-compounding) rate over a tenor
//! `delta`, and the option on such a rate is priced by converting the caplet into a put on the
//! underlying deposit bond — `1 + delta * K` has to stay positive for that map to exist, which is
//! what [`crate::validation::caplet_strike`] guards.  The floorlet is the same map with the bond
//! *call* in place of the bond put, so a cap and a floor on identical periods differ by exactly
//! the forward-rate leg they are written on, not by anything numerical.  The futures price differs
//! from the forward rate by the convexity term built from [`HullWhite::gamma_edf`].
//!
//! A cap (or floor) is a *schedule* of caplets (floorlets), each with its own expiry and strike.
//! [`HullWhite::cap_now`] / [`HullWhite::cap_t`] and [`HullWhite::floor_now`] /
//! [`HullWhite::floor_t`] take that schedule — `&[(option_maturity, strike)]` — and price it in
//! one call rather than leaving the caller to loop.

use crate::HullWhite;
use crate::error::HullWhiteError;
use crate::validation;

/// The simple (non-compounding, non-convexity-adjusted) rate over `tenor` implied by two zero
/// coupon bond prices: `(P_near - P_far) / (P_far * tenor)`.
///
/// Stays a free function rather than a method because it reads no model state at all — no `a`,
/// no `sigma`, no curve — just the two bond prices the caller already has.  A `&self` here would
/// advertise a dependency that is not in the arithmetic.  It used to sit in
/// [`crate::curves`] with the curve maths that does need a curve.
fn compute_libor_rate(nearest_bond: f64, farthest_bond: f64, tenor: f64) -> f64 {
    (nearest_bond - farthest_bond) / (farthest_bond * tenor)
}

/// The Eurodollar futures rate out of the two bracketing bond prices and the convexity gamma:
/// `(P_near / P_far * exp(gamma) - 1) / delta`.
///
/// Same reason [`compute_libor_rate`] stays free: the model's `a`, `sigma` and curve are all
/// already spent, into `bond_num`/`bond_den` and into `gamma`, by the time this runs.  The
/// difference from a plain forward Libor is the `exp(gamma)` — the daily-marked-to-market
/// settlement of the future, which [`HullWhite::gamma_edf`] prices.
fn edf_compute(bond_num: f64, bond_den: f64, gamma: f64, delta: f64) -> f64 {
    ((bond_num / bond_den) * gamma.exp() - 1.0) / delta
}

impl<'a> HullWhite<'a> {
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
    /// Returns price of a floorlet at current time
    ///
    /// A floorlet pays `delta * max(strike - Libor, 0)` at `option_maturity + delta` — the mirror
    /// image of [`HullWhite::caplet_now`], which pays `delta * max(Libor - strike, 0)` there.
    /// It uses the same deposit-bond map as the caplet, with the bond *call* where the caplet uses
    /// the bond put:
    ///
    /// ```text
    /// floorlet = (1 + delta*K) * call(P(., T+delta); K' = 1/(1 + delta*K), expiry T)
    /// caplet   = (1 + delta*K) * put (P(., T+delta); K' = 1/(1 + delta*K), expiry T)
    /// ```
    ///
    /// Subtracting the two gives the parity a floor has with a cap on the same period:
    ///
    /// ```text
    /// caplet - floorlet = P(0, T) - (1 + delta*K) * P(0, T + delta)
    ///                 = delta * (forward Libor(T, T + delta) - K) * P(0, T + delta)
    /// ```
    ///
    /// which is the value of the forward-rate leg itself: positive when the forward sits above the
    /// strike (the cap is then the expensive side), zero at the forward strike.
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let option_maturity = 1.5;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let strike = 0.04;
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let floorlet = hull_white.floorlet_now(option_maturity, delta, strike).unwrap();
    /// let caplet = hull_white.caplet_now(option_maturity, delta, strike).unwrap();
    /// assert!(floorlet > 0.0, "floorlet {floorlet}");
    /// //Cap-floor parity: the gap is the forward leg, priced off the two bond prices.
    /// let near = hull_white.bond_price_now(option_maturity).unwrap();
    /// let far = hull_white.bond_price_now(option_maturity + delta).unwrap();
    /// let leg = near - (1.0 + delta * strike) * far;
    /// assert!((caplet - floorlet - leg).abs() < 1e-12, "{caplet} - {floorlet} vs {leg}");
    /// ```
    pub fn floorlet_now(
        &self,
        option_maturity: f64,
        delta: f64,
        strike: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::strictly_after("option_maturity", option_maturity, "t", 0.0)?;
        validation::caplet_strike(delta, strike)?;
        self.bond_call_now(
            option_maturity,
            option_maturity + delta,
            1.0 / (delta * strike + 1.0),
        )
        .map(|call| (strike * delta + 1.0) * call)
    }
    /// Returns price of a floorlet at some future time
    ///
    /// Same map as [`HullWhite::floorlet_now`], from the state `(r_t, t)`; see that function for
    /// the cap-floor parity this sits on.
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let floorlet = hull_white.floorlet_t(r_t, t, option_maturity, delta, strike).unwrap();
    /// let caplet = hull_white.caplet_t(r_t, t, option_maturity, delta, strike).unwrap();
    /// let near = hull_white.bond_price_t(r_t, t, option_maturity).unwrap();
    /// let far = hull_white.bond_price_t(r_t, t, option_maturity + delta).unwrap();
    /// let leg = near - (1.0 + delta * strike) * far;
    /// assert!((caplet - floorlet - leg).abs() < 1e-12, "{caplet} - {floorlet} vs {leg}");
    /// ```
    pub fn floorlet_t(
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
        self.bond_call_t(
            r_t,
            t,
            option_maturity,
            option_maturity + delta,
            1.0 / (delta * strike + 1.0),
        )
        .map(|call| (strike * delta + 1.0) * call)
    }
    /// Returns price of a whole cap at current time.
    ///
    /// `periods` is the cap's schedule as `(option_maturity, strike)` pairs — one per caplet, so a
    /// 20 period cap is 20 pairs — and every period shares the tenor `delta`, the way one curve
    /// carries one Libor tenor.  The price is the sum of the individual caplets, each priced with
    /// [`HullWhite::caplet_now`]: a cap is a plain sum of its parts and aggregating them changes
    /// no single period's value.
    ///
    /// No ordering is assumed among the periods: they are priced independently, so unsorted or
    /// duplicated maturities are simply more caplets.  An empty schedule is an error rather than a
    /// `0.0` — a cap with no periods is not a free cap.
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// //A three period cap, each period with its own strike.
    /// let periods = [(1.0, 0.04), (1.25, 0.045), (1.5, 0.05)];
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let cap = hull_white.cap_now(&periods, delta).unwrap();
    /// let caplets: f64 = periods
    ///     .iter()
    ///     .map(|&(option_maturity, strike)| hull_white.caplet_now(option_maturity, delta, strike).unwrap())
    ///     .sum();
    /// assert!(cap > 0.0, "cap {cap}");
    /// assert!((cap - caplets).abs() <= 1e-12 * caplets.abs(), "{cap} vs {caplets}");
    /// //An empty schedule is refused, not priced at zero.
    /// assert!(hull_white.cap_now(&[], delta).is_err());
    /// ```
    pub fn cap_now(&self, periods: &[(f64, f64)], delta: f64) -> Result<f64, HullWhiteError> {
        validation::caplet_schedule(periods)?;
        validation::positive("delta", delta)?;
        let price = periods
            .iter()
            .map(|&(option_maturity, strike)| self.caplet_now(option_maturity, delta, strike))
            .sum::<Result<f64, HullWhiteError>>()?;
        validation::finish("cap_now", price)
    }
    /// Returns price of a whole cap at some future time.
    ///
    /// Same schedule as [`HullWhite::cap_now`] — `(option_maturity, strike)` per period, common
    /// tenor `delta` — priced from the state `(r_t, t)` with [`HullWhite::caplet_t`] per period.
    /// Every period's expiry has to be strictly after `t`; a period that already expired is an
    /// error rather than a silently dropped leg.
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //rate at time t
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let periods = [(1.5, 0.04), (1.75, 0.045), (2.0, 0.05)];
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let cap = hull_white.cap_t(r_t, t, &periods, delta).unwrap();
    /// let caplets: f64 = periods
    ///     .iter()
    ///     .map(|&(option_maturity, strike)| hull_white.caplet_t(r_t, t, option_maturity, delta, strike).unwrap())
    ///     .sum();
    /// assert!((cap - caplets).abs() <= 1e-12 * caplets.abs(), "{cap} vs {caplets}");
    /// //A period that expired before the valuation date is an error, not a dropped leg.
    /// assert!(hull_white.cap_t(r_t, t, &[(0.5, 0.04)], delta).is_err());
    /// ```
    pub fn cap_t(
        &self,
        r_t: f64,
        t: f64,
        periods: &[(f64, f64)],
        delta: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        validation::positive("delta", delta)?;
        validation::caplet_schedule(periods)?;
        let price = periods
            .iter()
            .map(|&(option_maturity, strike)| self.caplet_t(r_t, t, option_maturity, delta, strike))
            .sum::<Result<f64, HullWhiteError>>()?;
        validation::finish("cap_t", price)
    }
    /// Returns price of a whole floor at current time.
    ///
    /// The floor-side twin of [`HullWhite::cap_now`]: `periods` is `(option_maturity, strike)` per
    /// floorlet, all at tenor `delta`, and the price sums [`HullWhite::floorlet_now`].  On the
    /// *same* schedule and strikes the cap/floor gap telescopes into the forward legs:
    ///
    /// ```text
    /// cap - floor = sum_i [ P(0, T_i) - (1 + delta*K_i) * P(0, T_i + delta) ]
    /// ```
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let periods = [(1.0, 0.04), (1.25, 0.045), (1.5, 0.05)];
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let cap = hull_white.cap_now(&periods, delta).unwrap();
    /// let floor = hull_white.floor_now(&periods, delta).unwrap();
    /// let legs: f64 = periods
    ///     .iter()
    ///     .map(|&(option_maturity, strike)| {
    ///         let near = hull_white.bond_price_now(option_maturity).unwrap();
    ///         let far = hull_white.bond_price_now(option_maturity + delta).unwrap();
    ///         near - (1.0 + delta * strike) * far
    ///     })
    ///     .sum();
    /// assert!(floor > 0.0, "floor {floor}");
    /// let tol = 1e-11 * legs.abs().max(1.0);
    /// assert!((cap - floor - legs).abs() <= tol, "{cap} - {floor} vs {legs}");
    /// ```
    pub fn floor_now(&self, periods: &[(f64, f64)], delta: f64) -> Result<f64, HullWhiteError> {
        validation::caplet_schedule(periods)?;
        validation::positive("delta", delta)?;
        let price = periods
            .iter()
            .map(|&(option_maturity, strike)| self.floorlet_now(option_maturity, delta, strike))
            .sum::<Result<f64, HullWhiteError>>()?;
        validation::finish("floor_now", price)
    }
    /// Returns price of a whole floor at some future time.
    ///
    /// Same schedule as [`HullWhite::floor_now`], priced from the state `(r_t, t)` with
    /// [`HullWhite::floorlet_t`] per period.
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //rate at time t
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let periods = [(1.5, 0.04), (1.75, 0.045), (2.0, 0.05)];
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let floor = hull_white.floor_t(r_t, t, &periods, delta).unwrap();
    /// let floorlets: f64 = periods
    ///     .iter()
    ///     .map(|&(option_maturity, strike)| hull_white.floorlet_t(r_t, t, option_maturity, delta, strike).unwrap())
    ///     .sum();
    /// assert!(floor > 0.0, "floor {floor}");
    /// assert!((floor - floorlets).abs() <= 1e-12 * floorlets.abs(), "{floor} vs {floorlets}");
    /// ```
    pub fn floor_t(
        &self,
        r_t: f64,
        t: f64,
        periods: &[(f64, f64)],
        delta: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        validation::positive("delta", delta)?;
        validation::caplet_schedule(periods)?;
        let price = periods
            .iter()
            .map(|&(option_maturity, strike)| {
                self.floorlet_t(r_t, t, option_maturity, delta, strike)
            })
            .sum::<Result<f64, HullWhiteError>>()?;
        validation::finish("floor_t", price)
    }
    /// The Eurodollar convexity adjustment: the log-normal convexity of the fixing that makes a
    /// future worth more than the forward rate it is written on.
    ///
    /// A forward rate is paid at the *end* of its period and is discounted at the rate being
    /// fixed; a future is marked to market daily against that same fixing, so the margin is
    /// funded at the rate it is being marked by.  The futures rate is the forward rate times
    /// `exp(gamma)` — the `exp()` is applied by [`edf_compute`], not here.
    ///
    /// ```text
    /// gamma = (sigma^2 / a^3) * (1 - exp(-a * delta))
    ///         * [ (1 - exp(-a (option_maturity - t)))
    ///             - exp(-a * delta) * 0.5 * (1 - exp(-2a (option_maturity - t))) ]
    /// ```
    ///
    /// Derivations: <https://www.math.nyu.edu/~alberts/spring07/Lecture5.pdf> and the
    /// `OpenGamma` note
    /// <https://developers.opengamma.com/quantitative-research/Hull-White-One-Factor-Model-OpenGamma.pdf>,
    /// in whose notation the `t0` of the adjustment is this function's `option_maturity`.
    pub(crate) fn gamma_edf(&self, t: f64, option_maturity: f64, delta: f64) -> f64 {
        let exp_t = (-self.a * (option_maturity - t)).exp();
        let exp_d = (-self.a * delta).exp();
        (self.sigma.powi(2) / self.a.powi(3))
            * (1.0 - exp_d)
            * ((1.0 - exp_t) - exp_d * 0.5 * (1.0 - exp_t.powi(2)))
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
        let gamma = self.gamma_edf(t, option_maturity, delta);
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
    /// Price of a Eurodollar future at the current date (today, `t = 0`).
    ///
    /// Same contract as [`HullWhite::euro_dollar_future_t`] with the valuation collapsed to today:
    /// the state is not an argument, it is the calibration, so `r_t = r(0)`
    /// ([`HullWhite::short_rate_now`]) and the delivery date `option_maturity` is measured from now.
    /// The future is quoted on the Libor rate for the period starting at delivery, so the
    /// instrument is priced off the two zero-coupon prices that bracket that period
    /// (`option_maturity` and `option_maturity + delta`), which is why nothing here needs a rate at
    /// all: `P(0, T)` and `P(0, T + delta)` come straight off the curve.
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let option_maturity = 1.5; //delivery date, measured from today
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let edf = hull_white.euro_dollar_future_now(option_maturity, delta).unwrap();
    /// //The number returned is the futures rate itself (the convexity-adjusted Libor), not the
    /// //100-based index quote: it is a rate, so small and positive on a positive-rate curve.
    /// assert!(edf > 0.0 && edf < 0.5, "{edf}");
    /// //Futures settle daily against a falling-then-rising margin, which is worth a positive
    /// //convexity adjustment over the plain forward: the future is higher than forward Libor.
    /// let fwd = hull_white
    ///     .forward_libor_rate_now(option_maturity, delta)
    ///     .unwrap();
    /// assert!(edf > fwd, "{edf} vs {fwd}");
    /// //"Now" is exactly the `t = 0` case of the state-ful twin.
    /// let via_t = hull_white
    ///     .euro_dollar_future_t(hull_white.short_rate_now().unwrap(), 0.0, option_maturity, delta)
    ///     .unwrap();
    /// assert!((edf - via_t).abs() < 1e-12, "{edf} vs {via_t}");
    /// ```
    pub fn euro_dollar_future_now(
        &self,
        option_maturity: f64,
        delta: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::strictly_after("option_maturity", option_maturity, "t", 0.0)?;
        validation::positive("delta", delta)?;
        let gamma = self.gamma_edf(0.0, option_maturity, delta);
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
    /// Spot Libor rate at `t`: the rate fixed at `t` for borrowing over `[t, t + delta]`.
    ///
    /// The fixing itself, not a forecast of it: this is
    /// [`forward_libor_rate_t`](HullWhite::forward_libor_rate_t) with the period starting on the
    /// valuation date, so the near leg is `P(t, t) = 1` and the rate reduces to
    /// `(1 / delta) * (1 / P(t, t + delta) - 1)` under the state `r(t) = r_t`.  The `now` case of
    /// the same fixing — no state argument, read straight off the initial curve — is
    /// [`HullWhite::libor_rate_now`], and a period that starts later than the valuation date is the
    /// forward, not this.
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //current rate
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let libor = hull_white.libor_rate_t(r_t, t, delta).unwrap();
    /// ```
    pub fn libor_rate_t(&self, r_t: f64, t: f64, delta: f64) -> Result<f64, HullWhiteError> {
        validation::valuation_time(t)?;
        validation::positive("delta", delta)?;
        self.forward_libor_rate_t(r_t, t, t, delta)
    }
    /// Returns today's spot Libor rate: the rate fixed now for borrowing over `delta`.
    ///
    /// The `now` twin of [`HullWhite::libor_rate_t`] — the fixing whose settlement date is the
    /// valuation date — and so the same thing as `forward_libor_rate_now(0.0, delta)`.  Because
    /// both bond legs are priced with `bond_price_now`, this needs no short-rate argument at all:
    /// the spot fixing is read straight off the initial curve.
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let spot = hull_white.libor_rate_now(delta).unwrap();
    /// //Same instrument as the forward fixing starting today, and as the `t` form run at `t = 0`
    /// //with the model's own initial short rate.
    /// let forward = hull_white.forward_libor_rate_now(0.0, delta).unwrap();
    /// let r0 = hull_white.short_rate_now().unwrap();
    /// let via_t = hull_white.libor_rate_t(r0, 0.0, delta).unwrap();
    /// assert!((spot - forward).abs() < 1e-12, "{spot} vs {forward}");
    /// assert!((spot - via_t).abs() < 1e-12, "{spot} vs {via_t}");
    /// ```
    pub fn libor_rate_now(&self, delta: f64) -> Result<f64, HullWhiteError> {
        validation::positive("delta", delta)?;
        self.forward_libor_rate_now(0.0, delta)
    }
}

#[cfg(test)]
mod tests;
