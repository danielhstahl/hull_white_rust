//! Jamshidian's decomposition of a coupon-bond option, and the critical-rate solve it needs.
//!
//! A coupon bond option pays `(P_c(r_u) - K)^+` at expiry `u`, where the bond is a sum of
//! zero-coupon legs.  Jamshidian's trick rewrites that payoff as a sum of bond-option payoffs
//! `(x_i - K/N)^+` that are each priced as if their own bond leg were the underlying, which
//! requires the **critical rate** `r*`: the (unique, since the bond price strictly falls in the
//! rate) rate at which `P_c(r*) = K`.
//!
//! ## Which bond is the underlying: the deliverable convention
//!
//! `coupon_times` is the *whole* payment schedule of the bond, and `option_maturity` can fall
//! anywhere inside it: a 5y bond with a 2y option is the ordinary case, and until now that
//! schedule had to be rewritten to start after expiry to be priceable at all.  The payoff decides
//! what the underlying is.  `(P_c(r_u) - K)^+` is a claim on the bond *valued on the expiry
//! date*, so the underlying is whatever the holder has on that date, and the schedule splits into
//! three contiguous pieces ([`Deliverable::at_expiry`] does the split):
//!
//! * **payments strictly before `option_maturity` are dropped.**  They are cash the option holder
//!   never receives: exercising on the expiry date hands over a bond that has already paid them,
//!   and not exercising hands over nothing.  They are not part of `P_c(r_u)` for any `r_u`.
//! * **payments falling exactly on `option_maturity` are cash on the expiry date.**  The holder
//!   does get those, and their expiry-date value is their weight at `P(u, u) = 1` -- a constant,
//!   not a function of `r_u` -- so they move to the strike side exactly:
//!   `(cash + R(r_u) - K)^+ == (R(r_u) - (K - cash))^+`.
//! * **payments strictly after `option_maturity` are the residual bond** `R`: the legs the
//!   decomposition strikes.
//!
//! A schedule with nothing strictly after the expiry date is refused as
//! [`InvalidInput`](crate::HullWhiteError::InvalidInput): a bond fully settled by the expiry date
//! is not something an option can be written on.
//!
//! Why *drop* rather than net the pre-expiry coupons against the strike, the other candidate
//! convention.  The payoff is a claim on the bond at expiry, and a coupon paid before expiry is
//! worth nothing there, so no adjustment of either side of `(P_c(r_u) - K)^+` can put it back --
//! shrinking the strike by its present value would pay the holder for a cash flow they never got,
//! which is only right for an instrument whose holder *does* receive the interim coupons.  Netting
//! against the underlying instead is what you do when the input is the cum-coupon bond *price* at
//! `t`; here the model builds the underlying from the schedule, and that netting is already true
//! of it by no-arbitrage,
//!
//! ```text
//! full_bond_price(t) - PV(pre-expiry coupons) = P(t, u) * E^u[R(r_u)]
//! ```
//!
//! so that version of the convention collapses onto this one rather than producing a different
//! number.  Two consequences, both tested in the `straddling` tests:
//!
//! * dropping the pre-expiry payments cannot change the price, so passing the full schedule and
//!   passing only the tail from the first post-expiry payment onwards agree;
//! * the price is right-continuous in a coupon's date as that date slides across expiry from above
//!   (at `u + eps` the leg is worth about its weight, which is what the fold hands over at `u`).
//!   The jump the other way is the instrument's, not the convention's: a coupon one day before
//!   expiry really is not the same as one one day after.
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

/// Which side of the coupon-bond option is being priced.
///
/// The decomposition itself does not care: every leg is handed to the caller's Black-Scholes
/// closure.  The side is needed for one degenerate case only.  When the deliverable is worth more
/// than the strike in *every* rate state the option is exercised whatever happens, which makes the
/// call a parity value and the put nothing -- and that difference is not readable off the leg
/// closures, both of which return their own zero-strike answer.
pub(crate) enum Side {
    Call,
    Put,
}

/// A coupon schedule split at the option's expiry date, into what an exercising holder has.
///
/// See the [module docs](self) for the convention and why it is the only one the payoff supports.
/// `at_expiry` relies on [`validation::payment_schedule`] having already put `coupon_times` in
/// strictly ascending order, which makes the three pieces contiguous and lets `partition_point`
/// find both cuts in one pass.
pub(crate) struct Deliverable<'a> {
    /// Payments strictly after expiry: the residual bond, and the only legs the decomposition
    /// strikes.  Empty means the option has nothing to deliver.
    after_expiry: &'a [f64],
    /// Expiry-date value of the payments falling exactly on the expiry date -- cash, worth its
    /// weight at `P(u, u) = 1` -- which is subtracted from the strike.
    cash_at_expiry: f64,
    /// Payments made strictly before expiry, and therefore dropped.  Counted so the "nothing to
    /// deliver" error can say where the schedule went.
    paid_before_expiry: usize,
    /// The schedule's final payment, i.e. the bond's maturity, for that error message.
    last_payment: f64,
}

impl<'a> Deliverable<'a> {
    /// Split the whole schedule `coupon_times` at the expiry date `u`.
    fn at_expiry(coupon_times: &'a [f64], coupon_rate: f64, u: f64) -> Self {
        //Par rides on the schedule's final payment, so a payment's weight is a property of where it
        //sits in the whole schedule rather than of where it lands in a slice of it -- hence the
        //index rather than the slice.  Written as `index + 1 == len` so that there is no
        //`len() - 1` anywhere in here to underflow on an empty schedule.
        let len = coupon_times.len();
        let weight = |index: usize| coupon_rate + if index + 1 == len { 1.0 } else { 0.0 };
        let first_at_expiry = coupon_times.partition_point(|&time| time < u);
        let first_after_expiry = coupon_times.partition_point(|&time| time <= u);
        let cash_at_expiry = (first_at_expiry..first_after_expiry).map(weight).sum();
        Deliverable {
            after_expiry: &coupon_times[first_after_expiry..],
            cash_at_expiry,
            paid_before_expiry: first_at_expiry,
            last_payment: coupon_times.last().copied().unwrap_or(f64::NAN),
        }
    }
}

impl<'a> HullWhite<'a> {
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
            let c = ct_t(self.a, self.sigma, u, *coupon_time, self.curve());
            c_max = c_max.max(c);
            b_min = b_min.min(b);
            b_max = b_max.max(b);
        }
        let w_last = 1.0 + coupon_rate;
        let last_time = coupon_times[last_index];
        let b_last = at_t(self.a, u, last_time);
        let c_last = ct_t(self.a, self.sigma, u, last_time, self.curve());
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
    /// `coupon_times` is the *deliverable* schedule -- payments strictly after `option_maturity`,
    /// with anything settled on the expiry date already folded out of the strike.  See the
    /// [module docs](crate::jamshidian) for that split.
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

    /// Price a European option on a coupon bond by Jamshidian's decomposition.
    ///
    /// `coupon_times` is the whole schedule of the bond; the underlying is the bond *as it stands
    /// on the option's expiry date*, split by [`Deliverable::at_expiry`] as documented on the
    /// module.  Only the payments strictly after `option_maturity` reach the decomposition, and a
    /// payment landing exactly on `option_maturity` is folded into the strike, so everything below
    /// -- the bracket, the solve, the leg weights -- works on the residual bond with the
    /// cash-adjusted strike.
    ///
    /// `generic_fn` prices one zero-coupon leg option (a call or a put, matching `side`) and is
    /// what makes this one routine serve both `coupon_bond_call_t` and `coupon_bond_put_t`.
    #[allow(clippy::too_many_arguments)] //an instrument's whole description, already the crate's shape
    pub(crate) fn coupon_bond_option_generic_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
        side: Side,
        generic_fn: &impl Fn(f64, f64, f64, f64, f64) -> Result<f64, HullWhiteError>,
    ) -> Result<f64, HullWhiteError> {
        let par_value = 1.0;
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        validation::strictly_after("option_maturity", option_maturity, "t", t)?;
        validation::finite("coupon_rate", coupon_rate)?;
        validation::strike(strike)?;
        //The schedule has to be a live, ordered set of payments measured from `t`.  How it sits
        //relative to the expiry date is not a validity question -- see the module docs -- so it is
        //not checked here.
        validation::payment_schedule(coupon_times, t)?;
        let deliverable = Deliverable::at_expiry(coupon_times, coupon_rate, option_maturity);
        if deliverable.after_expiry.is_empty() {
            return Err(HullWhiteError::InvalidInput(format!(
                "coupon_times has no payment strictly after option_maturity = {option_maturity}: \
                 {paid} paid before it, {at_cash} on it, last payment at {last}.  A bond settled by \
                 the expiry date is not a deliverable underlying.",
                paid = deliverable.paid_before_expiry,
                at_cash = deliverable.cash_at_expiry,
                last = deliverable.last_payment,
            )));
        }
        //Cash settled on the expiry date is worth one per unit of weight there, whatever the rate
        //does, so folding it into the strike is an identity rather than an approximation:
        //  (cash + R(r_u) - K)^+  ==  (R(r_u) - (K - cash))^+.
        let strike = strike - deliverable.cash_at_expiry;
        let discount = self.bond_price_t_raw(r_t, t, option_maturity);
        //A strike at or below the folded cash leaves no state the option is not exercised in, so
        //there is no critical rate to find and the answer is the exercise value: for the call, the
        //deliverable less the strike discounted; for the put, nothing.  Note that `strike` is
        //already net of the expiry cash here, so the residual bond is the only thing left to value
        //-- adding the cash to the underlying as well would pay it out twice.
        //
        //Arriving here with a strictly negative strike needs `cash_at_expiry > strike >= 0`, hence
        //a positive coupon rate, hence every residual weight positive (the par-bearing weight is
        //1 + coupon_rate), hence the residual really is strictly positive in every expiry state and
        //the exercise value really is the whole answer rather than a floor.  A zero strike with no
        //folded cash is the case that has always been here: the call is the underlying, the put is
        //nothing.
        if strike <= 0.0 {
            let residual_at_t = coupon_bond_generic_t(
                r_t,
                t,
                deliverable.after_expiry,
                coupon_rate,
                &|rate: f64, leg_t: f64, maturity: f64| {
                    self.bond_price_t_raw(rate, leg_t, maturity)
                },
            );
            let price = match side {
                Side::Call => residual_at_t - strike * discount,
                Side::Put => 0.0,
            };
            return validation::finish("Jamshidian decomposition of a coupon-bond option", price);
        }
        //Each leg is struck at the zero-coupon price that puts the *whole* deliverable level with
        //the strike at the critical rate.
        let solution = self.solve_critical_rate(
            r_t,
            t,
            option_maturity,
            deliverable.after_expiry,
            coupon_rate,
            strike,
        )?;
        let last_index = deliverable.after_expiry.len();
        deliverable
            .after_expiry
            .iter()
            .enumerate()
            .map(|(index, coupon_time)| {
                //The residual slice is a suffix of the schedule, so its last payment is the
                //schedule's last payment and carries the par.
                let is_last = index + 1 == last_index;
                let strike_leg =
                    self.bond_price_t_raw(solution.root, option_maturity, *coupon_time);
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
