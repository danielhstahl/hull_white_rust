//! Vanilla swaps: forward swap rate, swap price, and the European (analytic) swaptions.
//!
//! The payment convention every price here is built on — `n` coupon dates, `start + delta` through
//! `start + n * delta`, with the principal riding on the last of them, so `swap_maturity` *is*
//! that last date — is stated once, on the module-private `annuity_t` that both
//! [`HullWhite::forward_swap_rate_t`] and [`HullWhite::swap_price_t_init`] call, with the worked
//! example and the algebra that shows the two consumers agree.
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

impl<'a> HullWhite<'a> {
    /// The annuity of a swap's remaining payments:
    ///
    /// ```text
    /// A(t) = delta * sum_{i=1..n} P(t, start + i * delta)
    /// ```
    ///
    /// where `P` is [`HullWhite::bond_price_t_raw`].  Both things this module prices are read off
    /// this one sum — the forward rate is `(P(t, start) - P(t, start + n*delta)) / A(t)` and the
    /// swap price is `P(t, start) - K * A(t) - P(t, start + n*delta)` — and those two only agree
    /// with each other if the schedule behind the sum is counted the same way in both.  So the
    /// payment convention lives here, once, rather than being implicit in two hand-rolled loops.
    ///
    /// # The payment convention
    ///
    /// A schedule is three numbers: `start` (when the swap begins and the floating leg resets to
    /// par), `delta` (one payment period) and `n` (the caller's `num_swap_payments`, how many
    /// coupons are left).  They generate the `n` payment dates
    /// `start + delta, start + 2*delta, ..., start + n*delta`, and on those dates:
    ///
    /// * **every one of the `n` dates carries a coupon** — `K * delta` on the fixed leg, the
    ///   period's reset Libor times `delta` on the floating leg;
    /// * **the principal — 1.0 of notional — is carried by the last of those same `n` dates**,
    ///   together with that date's coupon.  Exchanging principal does not append a date.
    ///
    /// Worked example: `start = 1.0`, `delta = 0.25`, `n = 4` (quarterly, starting in a year):
    ///
    /// ```text
    /// i                  1        2        3        4
    /// date              1.25     1.50     1.75  2.00  = swap_maturity
    /// coupon            K*d      K*d      K*d     K*d
    /// principal         --       --       --      + 1.00
    /// fixed cash flow   K*d      K*d      K*d     K*d + 1.00
    ///
    /// A = d * ( P(1.25) + P(1.50) + P(1.75) + P(2.00) )
    /// fixed leg value   = K * A + P(2.00)
    /// ```
    ///
    /// and `swap_maturity` is `start + n * delta` — 2.00 above, i.e. `1.0 + 4 * 0.25`.  Note the
    /// maturity date is *in* the schedule (it is the `i = n` member), and principal sits on it on
    /// top of that date's coupon.
    ///
    /// # Resolved: `num_swap_payments` is `n`, not `n + 1`
    ///
    /// This function used to carry `//open question, should num_swap_payments be
    /// num_swap_payments+1??`.  It was written against the old body, which summed the fixed
    /// coupons over `1..n` — visibly stopping one short of `n` — and then paid `(1 + K*delta)` at
    /// the `n`-th date.  The "missing" `n`-th coupon is not missing: the `K*delta` inside that
    /// `(1 + K*delta)` *is* the `i = n` coupon, paid alongside the principal on the maturity
    /// date.  Extending the coupon sum to `n` and folding the principal back out is the same
    /// number, term by term (writing `P_i` for `P(t, start + i * delta)`):
    ///
    /// ```text
    ///   P(t,start) - sum_{i=1..n-1} K*d*P_i - (1 + K*d)*P_n
    /// = P(t,start) - K*d*sum_{i=1..n}   P_i - P_n        (K*d*P_n added to each sum)
    /// = P(t,start) - K*A(t)             - P_n
    /// ```
    ///
    /// which is what [`HullWhite::swap_price_t_init_raw`] now returns.  Passing `n + 1` instead
    /// would not just add a payment, it would move one: the schedule would run to `start +
    /// (n+1)*delta`, so the worked example would pay a fifth coupon — and the principal — at 2.25
    /// on a swap that matured at 2.00.
    ///
    /// The check that ties the two consumers together: setting fixed leg equal to floating leg
    /// gives `K * A(t) + P(t, start + n*delta) = P(t, start)`, i.e.
    /// `K = (P(t,start) - P(t,start+n*delta)) / A(t)`, which is exactly what
    /// [`HullWhite::forward_swap_rate_t`] computes.  Feed that `K` back into the price and every
    /// term cancels, which is another way of saying the two functions are counting the same four
    /// dates.
    ///
    /// # Why the floating leg needs no annuity
    ///
    /// A floating-rate note is worth par on a reset date, so the whole remaining float leg — all
    /// `n` of its reset coupons and its principal — is worth `1.0` at `start`, and `P(t, start)`
    /// from `t`.  That single bond is the entire floating side; no-arbitrage is what makes its
    /// sum of reset coupons plus principal collapse to the par it resets to.  The fixed leg has
    /// no such shortcut, because its coupon does not reset: that sum is [`HullWhite::annuity_t`].
    fn annuity_t(&self, r_t: f64, t: f64, start: f64, n: usize, delta: f64) -> f64 {
        (1..=n)
            .map(|i| self.bond_price_t_raw(r_t, t, start + delta * (i as f64)))
            .sum::<f64>()
            * delta
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let forward_swap = hull_white.forward_swap_rate_t(r_t,  t, swap_initiation, num_swap_payments, delta).unwrap();
    /// ```
    pub fn forward_swap_rate_t(
        &self,
        r_t: f64,
        t: f64,
        swap_initiation: f64, //must be greater than or equal to t
        num_swap_payments: usize,
        delta: f64,
    ) -> Result<f64, HullWhiteError> {
        validation::finite("r_t", r_t)?;
        validation::valuation_time(t)?;
        validation::not_before("swap_initiation", swap_initiation, "t", t)?;
        validation::at_least_one("num_swap_payments", num_swap_payments)?;
        validation::positive("delta", delta)?;
        let denominator_swap = self.annuity_t(r_t, t, swap_initiation, num_swap_payments, delta);
        //The two zero-coupon legs the swap collapses to: par discounted from the reset date (the
        //floating side, worth par at `swap_initiation`) and the principal discounted from the last
        //payment date.  Under the convention stated on `annuity_t`,
        //`swap_initiation + num_swap_payments * delta` *is* the swap's maturity -- the `i = n`
        //payment date, the one that carries the principal on top of its coupon -- not maturity
        //plus a period.
        //One cast of the payment count into the time domain, named once instead of inline in the
        //bond call: `swap_initiation` plus a full set of `delta` periods *is* the last payment
        //date, i.e. the swap's maturity under the convention stated on `annuity_t`.
        let swap_maturity = swap_initiation + num_swap_payments as f64 * delta;
        let par_at_start = self.bond_price_t_raw(r_t, t, swap_initiation);
        let principal_at_maturity = self.bond_price_t_raw(r_t, t, swap_maturity);
        validation::finish(
            "forward_swap_rate_t",
            (par_at_start - principal_at_maturity) / denominator_swap,
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
    /// Swap rate at time `t` for the swap that starts at `t`, given the short rate `r(t) = r_t`.
    ///
    /// The spot swap rate as seen from `t`: [`HullWhite::forward_swap_rate_t`] with
    /// `swap_initiation == t`, so the fixed leg is the coupon bond running from `t` and the rate is
    /// the one that makes the swap worth zero at that moment.  Rate, not price — for the value of a
    /// swap already struck at some other rate use [`HullWhite::swap_price_t`].  The same quantity
    /// from today, with the state coming from the calibration rather than an argument, is
    /// [`HullWhite::swap_rate_now`].
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
    /// Today's spot swap rate: the fixed rate that values a swap starting today at zero.
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
    /// The remaining-payment count and the `swap_price_t_init` anchor are derived from
    /// `swap_maturity`.  When maturity is a whole number of periods from `t` the anchor is `t`
    /// itself, which is exactly what this module's `annuity_t` convention wants: payments
    /// at `t + delta, ..., t + n * delta = swap_maturity`, principal on the last one.
    ///
    /// **Known gap, deliberately not fixed here.**  When maturity is *not* a whole number of
    /// periods past `t`, the branch below anchors on `swap_maturity - (n - 1) * delta`, which is
    /// the **next payment date** (`schedules`' own tests name that expression
    /// `next_exchange_date`), not the reset date one period before it.  Under the convention the
    /// anchor is `swap_maturity - n * delta`, so the off-schedule branch drops the first remaining
    /// payment and repays principal a period past the stated maturity: with `t = 0.5`,
    /// `swap_maturity = 2.0`, `delta = 0.4` it prices coupons at 1.2 / 1.6 / 2.0 / **2.4** with
    /// principal at 2.4, where the real remaining schedule is 0.8 / 1.2 / 1.6 / 2.0 with principal
    /// at 2.0.  The payment *count* is right; only the anchor is a period late.  Correcting it
    /// moves prices by a whole period, not by rounding, so it belongs in its own change rather than
    /// riding along inside a refactor that otherwise moves nothing past the last ulp.
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
        //Floating leg at par, less the fixed leg (K * annuity), less the principal repaid on the
        //last payment date -- which is also where the last coupon of the annuity falls, see
        //`annuity_t` for the convention and for why `num_swap_payments` needs no +1 here.
        let swap_maturity = swap_start + num_swap_payments as f64 * delta;
        self.bond_price_t_raw(r_t, t, swap_start)
            - swap_rate * self.annuity_t(r_t, t, swap_start, num_swap_payments, delta)
            - self.bond_price_t_raw(r_t, t, swap_maturity)
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
    /// Price of a **receiver** swaption at some future time `t`.
    ///
    /// The holder has the right, at `option_maturity`, to enter the swap as the **fixed-rate
    /// receiver** (and floating payer).  Worth having when the fixed rate they will lock in beats
    /// the forward swap rate, i.e. when the fixed-coupon bond they are about to be long is worth
    /// more than par — which is exactly why this is a **call** on that coupon bond:
    /// `max(Bond(swap_rate * delta) - 1, 0)`, struck at par, on the schedule
    /// `get_coupon_times(num_swap_payments, option_maturity, delta)`.  Mirror image of
    /// [`HullWhite::european_payer_swaption_t`], which is the corresponding put.
    ///
    /// # Examples
    ///
    /// ```
    /// let r_t = 0.04; //short rate observed at t
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let t = 1.0; //time from "now" (0) to start valuing the bond
    /// let option_maturity = 2.0;
    /// let num_swap_payments = 16;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// //The rate the receiver is indifferent at: the forward swap rate of the swap the option
    /// //delivers.  Everything below is stated against that rather than a hard-coded strike, so it
    /// //stays true whatever curve the model is calibrated to.
    /// let atm = hull_white
    ///     .forward_swap_rate_t(r_t, t, option_maturity, num_swap_payments, delta)
    ///     .unwrap();
    /// let receiver = |strike: f64| {
    ///     hull_white
    ///         .european_receiver_swaption_t(
    ///             r_t, t, option_maturity, num_swap_payments, delta, strike,
    ///         )
    ///         .unwrap()
    /// };
    /// let payer = |strike: f64| {
    ///     hull_white
    ///         .european_payer_swaption_t(r_t, t, option_maturity, num_swap_payments, delta, strike)
    ///         .unwrap()
    /// };
    /// //Above the forward the receiver is in the money: they lock in a fixed rate richer than the
    /// //market's, and the fixed-coupon bond they receive is worth more than the par they pay.
    /// assert!(receiver(atm + 0.01) > payer(atm + 0.01));
    /// //Below it the trade is worth more to the payer.
    /// assert!(receiver(atm - 0.01) < payer(atm - 0.01));
    /// //At the forward the two sides are worth the same: `receiver - payer` is the call-put spread
    /// //on the deliverable, which is its forward price minus the par strike, and that is zero here.
    /// assert!((receiver(atm) - payer(atm)).abs() < 1e-12);
    /// //Being paid a higher fixed rate is monotonically better for a receiver, worse for a payer.
    /// assert!(receiver(atm + 0.01) > receiver(atm - 0.01));
    /// assert!(payer(atm + 0.01) < payer(atm - 0.01));
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
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
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
