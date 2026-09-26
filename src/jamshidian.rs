//! Jamshidian's decomposition of a coupon-bond option, and the critical-rate solve it needs.
//!
//! A coupon bond option pays `(P_c(r_u) - K)^+` at expiry `u`, where the bond is a sum of
//! zero-coupon legs.  Jamshidian's trick rewrites that payoff as a sum of bond-option payoffs
//! `(x_i - K/N)^+` that are each priced as if their own bond leg were the underlying, which
//! requires the **critical rate** `r*`: the (unique, since the bond price strictly falls in the
//! rate) rate at which `P_c(r*) = K`.
//!
//! [`HullWhite::critical_rate_bracket`] derives a bracket for `r*` analytically from the leg
//! structure — one-sided bounds on the sum of exponentials, costing no function evaluations —
//! and [`HullWhite::solve_critical_rate`] runs the safeguarded solve inside it.  Both exits of
//! that solve are in *rate* space; see the crate docs for why a price residual is not an exit.

use crate::HullWhite;
use crate::bonds::coupon_bond_generic_t;
use crate::curves::{at_t, ct_t};
use crate::error::HullWhiteError;
use crate::rootfinder::{self, Solution};
use crate::validation;

impl<'a, T, U> HullWhite<'a, T, U>
where
    T: Fn(f64) -> f64 + std::marker::Sync,
    U: Fn(f64) -> f64 + std::marker::Sync,
{
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
    pub(crate) fn critical_rate_bracket(
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
    ) -> Result<Solution, HullWhiteError> {
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

    pub(crate) fn coupon_bond_option_generic_t(
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
}

#[cfg(test)]
mod tests;
