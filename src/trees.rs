//! Short-rate tree plumbing: the European tree check and the American (early-exercise) swaptions.
//!
//! Both walk a trinomial-style Black-Vasicek tree from `binomial_tree`, with the model's `phi(t)`
//! cached per time step so the state is the shifted rate `y = r - phi`.  Tree time `tau` means
//! *absolute* model time `t + tau` everywhere: the option is priced at `t`, the tree runs for
//! `option_maturity - t`, and every call into `phi` or a swap leg adds `t` back.
//!
//! The wiring exists once: [`HullWhite::tree_price`] is the engine -- drift, volatility, the
//! per-step `phi` cache, the node discount, the exercise style -- and takes the payoff as a
//! closure written against the short rate at the node.  [`HullWhite::swaption_tree`] turns a
//! [`SwaptionSpec`] into such a payoff, and the `american_*` entry points hand it [`AMERICAN`]
//! while the `european_*_tree` entry points hand it [`EUROPEAN`].  Before this was folded there
//! were two hand-maintained copies of the same thirty-five lines -- and a third in the Jamshidian
//! straddling test -- which meant a fix to the time-coordinate convention above, the one bug this
//! module has actually had to carry, had to be made in every one of them and could be made in one
//! and leave the others wrong.
//!
//! The European tree methods are public rather than `cfg(test)`, deliberately.  A tree and the
//! analytic Jamshidian decomposition reach the same swaption price by independent means, so the
//! tree is the cross-check a consumer can run against
//! [`european_payer_swaption_t`](crate::HullWhite::european_payer_swaption_t) without taking
//! either route on faith; while the helper was test-only that check could not be run from outside
//! the crate at all.  They are verification helpers, not the preferred pricer: the closed-form
//! European price in [`crate::swaps`] is exact within the model's own approximation and far
//! cheaper than a tree.  No `american_*` public signature changed.

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

/// One swaption as the tree sees it: the state, the swap leg, the side and the resolution.
///
/// The fields are exactly the arguments the public entry points take, bundled so the shared
/// pricing helper [`HullWhite::swaption_tree`] is called as `swaption_tree(&spec, is_american)`
/// rather than as eight positional parameters whose types are five `f64`s, two `usize`s and a
/// `bool` -- the shape in which two adjacent same-typed arguments can be swapped without the
/// compiler noticing, which is exactly how a tree-wiring bug hides.  It stays private: the public
/// API keeps its named, positional entry points (`american_payer_swaption_t` and friends), so no
/// existing caller changes.
struct SwaptionSpec {
    /// Short rate observed at the valuation time `t`.
    r_t: f64,
    /// Valuation time, measured from "now" (0).
    t: f64,
    /// Option maturity, measured from "now".  The tree spans `option_maturity - t`.
    option_maturity: f64,
    /// Number of remaining payments on the underlying swap.
    num_swap_payments: usize,
    /// Tenor of the simple yield: the period length of each swap payment.
    delta: f64,
    /// The fixed rate the payer agrees to pay on exercise.
    swap_rate: f64,
    /// `true` for the payer side, `false` for the receiver side.
    is_payer: bool,
    /// Number of tree steps across the option's life.
    num_steps: usize,
}

/// Named values for the `is_american` flag on [`HullWhite::swaption_tree`].
///
/// That flag is the entire difference between the two products this module prices: early exercise
/// compared against the continuation value at every node, or plain discounted expectation at
/// expiry.  Four entry points pass it, and a bare `true`/`false` typo in one of them is silent --
/// the tests below caught exactly that while this file was being folded -- so the call sites say
/// which style they mean.
const AMERICAN: bool = true;
/// See [`AMERICAN`]: no early exercise, the payoff is taken at expiry only.
const EUROPEAN: bool = false;

/// The input contract every swaption entry point on the tree is checked against.
///
/// The same seven rules, in the same order and with the same messages, that each of the four
/// `american_*` entry points used to carry its own copy of.  The European tree entry points
/// inherit them rather than re-typing them, so an argument rejected on one side of the tree is
/// rejected on the other.
fn validate_swaption(spec: &SwaptionSpec) -> Result<(), HullWhiteError> {
    validation::finite("r_t", spec.r_t)?;
    validation::valuation_time(spec.t)?;
    validation::strictly_after("option_maturity", spec.option_maturity, "t", spec.t)?;
    validation::at_least_one("num_swap_payments", spec.num_swap_payments)?;
    validation::positive("delta", spec.delta)?;
    validation::finite("swap_rate", spec.swap_rate)?;
    validation::at_least_one("num_steps", spec.num_steps)?;
    Ok(())
}

impl<'a> HullWhite<'a> {
    /// The engine every tree-priced instrument in this crate rolls its payoff back through: the
    /// drift, the (zero) drift derivative and the inverse volatility of the shifted short rate,
    /// the per-step `phi` cache, the one-step discount at each node, and the exercise style.
    ///
    /// `payoff(tau, r, dt, j)` is written against **`r`, the short rate at the node**, on tree
    /// time `tau`; the engine un-shifts the tree state (`y = r - phi`) so that callers never have
    /// to touch the cache themselves.  `is_american` switches the early-exercise comparison on at
    /// every node; with it off the payoff is taken at expiry and rolled back as a plain discounted
    /// expectation.
    ///
    /// Keeping the time-coordinate convention in this one function is the point.  The tree runs on
    /// a clock shifted by the valuation time `t` (it spans `horizon`, not `option_maturity`), but
    /// `phi` and every leg of every underlying are measured from "now" (0).  A node at tree time
    /// `tau` is therefore at absolute time `t + tau`, and using the shifted time directly
    /// misprices the legs whenever `t > 0` (or whenever `phi` is not constant).  That was one bug
    /// paid for three times -- here, in the European tree copy, and in the reference tree in
    /// `jamshidian::tests::straddling` -- before this was the only place it is written down.
    pub(crate) fn tree_price(
        &self,
        r_t: f64,
        t: f64,
        horizon: f64,
        num_steps: usize,
        is_american: bool,
        payoff: &dyn Fn(f64, f64, f64, usize) -> f64,
    ) -> f64 {
        //The three callbacks that describe the diffusion, in this model's terms.  Note that the
        //second argument is not the same kind of number in each of them: `binomial_tree` takes a
        //diffusion `dX = alpha dt + sigma dW` and builds its lattice on a pure-Brownian
        //coordinate `w`, recovering the state from `w` through the `sigma_inverse` callback.
        //The state diffused here is not the short rate but the *shifted* rate `y = r - phi(t)`,
        //which obeys `dy = -a y dt + sigma dW` -- drift `-a y`, a constant volatility.  So:
        //
        //  * `alpha_div_sigma` and `d_sigma_d_state` are handed the state `y`;
        //  * `state_from_lattice_coord` -- the engine's `sigma_inverse` -- is handed the
        //    lattice coordinate `w` and *returns* the state: the primitive of `1 / sigma` is
        //    `y -> y / sigma`, so its inverse is `w -> sigma * w`.  It used to be named
        //    `sigma_inv` with its argument named `y`, which read as though the engine handed
        //    over a rate that then got multiplied by sigma, instead of being the `w -> y` map.
        //    And `sigma_prime` is `d sigma / d(state)`: identically zero here, because a
        //    constant volatility does not vary with the state.
        let alpha_div_sigma =
            |_t_step: f64, y: f64, _dt: f64, _width: usize| -(self.a * y) / self.sigma;
        let d_sigma_d_state = |_t_step: f64, _state: f64, _dt: f64, _j: usize| 0.0;
        let state_from_lattice_coord = |_t_step: f64, w: f64, _dt: f64, _j: usize| self.sigma * w;
        let mut phi_cache: Vec<f64> = binomial_tree::get_all_t(horizon, num_steps)
            .map(|tau| self.phi_t(t + tau))
            .collect();
        phi_cache.push(self.phi_t(t + horizon));
        //The tree's state is the shifted rate `y`; the payoff wants the rate.  Un-shift here rather
        //than repeating `y + phi_cache[j]` at every payoff site, which is how the index `j` and the
        //time base get out of step.
        let node_payoff =
            |t_step: f64, y: f64, dt: f64, j: usize| payoff(t_step, y + phi_cache[j], dt, j);
        let discount = |_t_step: f64, curr_val: f64, dt: f64, j: usize| {
            (-(curr_val + phi_cache[j]) * dt).exp()
        };
        binomial_tree::compute_price_raw(
            &alpha_div_sigma,
            &d_sigma_d_state,
            &state_from_lattice_coord,
            &node_payoff,
            &discount,
            //The starting lattice coordinate `w`: the value that `state_from_lattice_coord`
            //sends to the initial state `y = r_t - phi(t)`.
            (r_t - self.phi_t(t)) / self.sigma,
            horizon,
            num_steps,
            is_american,
        )
    }
    /// The swaption on the tree: the payoff of the underlying swap leg at each node, handed to
    /// [`HullWhite::tree_price`] with the exercise style.
    ///
    /// Upstream, `compute_price_american` *is* `compute_price_raw(..., true)`, so routing both
    /// sides through the one engine leaves the American numbers bit-identical to what they were
    /// before the two copies were folded; `swaption_tree_at_t_is_bit_identical_to_pinned_values`
    /// pins that with bits rather than with a tolerance.
    ///
    /// The caller validates (see [`validate_swaption`]); this is the kernel and does not re-check
    /// its arguments.
    fn swaption_tree(&self, spec: &SwaptionSpec, is_american: bool) -> f64 {
        self.tree_price(
            spec.r_t,
            spec.t,
            spec.option_maturity - spec.t,
            spec.num_steps,
            is_american,
            &|t_step: f64, rate: f64, _dt: f64, _j: usize| {
                //Priced at the node's absolute time on the rate observed at that node: the swap
                //leg is worth what it is worth then, not what it was worth at the valuation date.
                let t_abs = spec.t + t_step;
                payoff_swaption(
                    spec.is_payer,
                    self.swap_price_t_init_raw(
                        rate,
                        t_abs,
                        t_abs,
                        spec.num_swap_payments,
                        spec.delta,
                        spec.swap_rate,
                    ),
                )
            },
        )
    }
    /// Returns the price of an American payer swaption at some future time t, by tree.
    ///
    /// # Comments
    ///
    /// This function uses a tree to solve and will take longer to compute
    /// than other pricing functions.
    ///
    /// The parameter list mirrors the swap leg being priced -- current rate, clock,
    /// option maturity, number of payments, payment tenor, strike, tree resolution.
    /// Folding it into a parameter struct would change the public signature, so the
    /// `too_many_arguments` allow on this item is deliberate, not unfinished work.
    /// (Inside the module the bundling has already happened: this builds a `SwaptionSpec`
    /// and hands it to the shared `swaption_tree` kernel -- the same one the European
    /// tree prices run through.)
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
    /// let num_tree_steps = 100;
    /// let swaption = hull_white.american_payer_swaption_t(
    ///     r_t, t, option_maturity, num_swap_payments, delta, swap_rate, num_tree_steps
    /// ).unwrap();
    /// ```
    #[allow(clippy::too_many_arguments)]
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
        let spec = SwaptionSpec {
            r_t,
            t,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            is_payer: true,
            num_steps,
        };
        validate_swaption(&spec)?;
        Ok(self.swaption_tree(&spec, AMERICAN))
    }
    /// Returns price of an American payer swaption at current time
    ///
    /// Exactly [`HullWhite::american_payer_swaption_t`] at `t = 0` with `r_t = r(0)`
    /// ([`HullWhite::short_rate_now`]), so the tree starts on the model's initial short rate.
    ///
    /// # Comments
    ///
    /// Tree-based like [`HullWhite::american_payer_swaption_t`]: slower than the analytic pricers.
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let option_maturity = 2.0;
    /// let num_swap_payments = 16;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let swap_rate = 0.04; //the swap rate is what the payer agrees to pay if option is exercised
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let num_tree_steps = 100;
    /// let swaption = hull_white.american_payer_swaption_now(
    ///     option_maturity, num_swap_payments, delta, swap_rate, num_tree_steps
    /// ).unwrap();
    /// let r0 = hull_white.short_rate_now().unwrap();
    /// let via_t = hull_white.american_payer_swaption_t(
    ///     r0, 0.0, option_maturity, num_swap_payments, delta, swap_rate, num_tree_steps
    /// ).unwrap();
    /// assert!((swaption - via_t).abs() < 1e-12, "{swaption} vs {via_t}");
    /// assert!(swaption > 0.0, "american payer swaption {swaption}");
    /// ```
    pub fn american_payer_swaption_now(
        &self,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
        num_steps: usize,
    ) -> Result<f64, HullWhiteError> {
        let t = 0.0; //since "now"
        let r_t = self.short_rate_now()?;
        self.american_payer_swaption_t(
            r_t,
            t,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            num_steps,
        )
    }
    /// Returns price of an American receiver swaption at some future time t
    ///
    /// # Comments
    ///
    /// This function uses a tree to solve and will take longer to compute
    /// than other pricing functions.
    ///
    /// The parameter list mirrors the swap leg being priced -- current rate, clock,
    /// option maturity, number of payments, payment tenor, strike, tree resolution.
    /// Folding it into a parameter struct would change the public signature, so the
    /// `too_many_arguments` allow on this item is deliberate, not unfinished work.
    /// (Inside the module the bundling has already happened: this builds a `SwaptionSpec`
    /// and hands it to the shared `swaption_tree` kernel -- the same one the European
    /// tree prices run through.)
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
    /// let swap_rate = 0.04; //the swap rate is what the receiver receives if option is exercised
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let num_tree_steps = 100;
    /// let swaption = hull_white.american_receiver_swaption_t(
    ///     r_t, t, option_maturity, num_swap_payments, delta, swap_rate, num_tree_steps
    /// ).unwrap();
    /// ```
    #[allow(clippy::too_many_arguments)]
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
        let spec = SwaptionSpec {
            r_t,
            t,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            is_payer: false,
            num_steps,
        };
        validate_swaption(&spec)?;
        Ok(self.swaption_tree(&spec, AMERICAN))
    }
    /// Returns price of an American receiver swaption at current time
    ///
    /// Exactly [`HullWhite::american_receiver_swaption_t`] at `t = 0` with `r_t = r(0)`
    /// ([`HullWhite::short_rate_now`]), so the tree starts on the model's initial short rate.
    ///
    /// # Comments
    ///
    /// Tree-based like [`HullWhite::american_receiver_swaption_t`]: slower than the analytic
    /// pricers.
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2; //speed of mean reversion for underlying Hull White process
    /// let sigma = 0.3; //volatility of underlying Hull White process
    /// let option_maturity = 2.0;
    /// let num_swap_payments = 16;
    /// let delta = 0.25; //delta is the tenor of the Libor rate
    /// let swap_rate = 0.04; //the swap rate is what the receiver receives if option is exercised
    /// // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose
    /// // derivative f(0,t) = 0.05 + 0.02 t is the instantaneous forward.
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let num_tree_steps = 100;
    /// let swaption = hull_white.american_receiver_swaption_now(
    ///     option_maturity, num_swap_payments, delta, swap_rate, num_tree_steps
    /// ).unwrap();
    /// let r0 = hull_white.short_rate_now().unwrap();
    /// let via_t = hull_white.american_receiver_swaption_t(
    ///     r0, 0.0, option_maturity, num_swap_payments, delta, swap_rate, num_tree_steps
    /// ).unwrap();
    /// assert!((swaption - via_t).abs() < 1e-12, "{swaption} vs {via_t}");
    /// assert!(swaption > 0.0, "american receiver swaption {swaption}");
    /// ```
    pub fn american_receiver_swaption_now(
        &self,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
        num_steps: usize,
    ) -> Result<f64, HullWhiteError> {
        let t = 0.0; //since "now"
        let r_t = self.short_rate_now()?;
        self.american_receiver_swaption_t(
            r_t,
            t,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            num_steps,
        )
    }
    /// Prices a European payer swaption on the short-rate tree, as an independent check on
    /// [`HullWhite::european_payer_swaption_t`].
    ///
    /// # Comments
    ///
    /// **This is a verification helper, not the recommended way to price a European swaption.**
    /// [`HullWhite::european_payer_swaption_t`] is the closed-form (Jamshidian) price and is what
    /// the rest of the crate uses; this walks the same trinomial tree as
    /// [`HullWhite::american_payer_swaption_t`] but with early exercise switched off, so the two
    /// routes agree only to the tree's discretisation error.  That is the point: run them against
    /// each other and a wiring bug in either one shows up as a gap, with no third pricer required.
    ///
    /// Same instrument and same argument list as the analytic pricer, plus `num_steps` -- the tree
    /// resolution.  Refining it shrinks the gap (`O(1/num_steps)` in practice here), so the
    /// tolerance in a cross-check should be read as a statement about the tree, never widened to
    /// absorb a real disagreement.
    ///
    /// `american_*` signatures are unchanged by the existence of this method.
    ///
    /// # Examples
    ///
    /// ```
    /// // Steep calibration: phi(t) travels, so a time-coordinate mistake cannot hide.
    /// let curr_rate = 0.02;
    /// let a = 0.2;
    /// let sigma = 0.03;
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let t = 0.5; //valuation time, strictly after "now"
    /// let option_maturity = 1.5;
    /// let num_swap_payments = 20;
    /// let delta = 0.25;
    /// let swap_rate = hull_white
    ///     .forward_swap_rate_t(curr_rate, t, option_maturity, num_swap_payments, delta)
    ///     .unwrap();
    /// let analytic = hull_white
    ///     .european_payer_swaption_t(
    ///         curr_rate, t, option_maturity, num_swap_payments, delta, swap_rate
    ///     )
    ///     .unwrap();
    /// let tree = hull_white
    ///     .european_payer_swaption_tree(
    ///         curr_rate, t, option_maturity, num_swap_payments, delta, swap_rate, 400
    ///     )
    ///     .unwrap();
    /// //Measured residual here: 7.1e-6 (payer) / 1.2e-5 (receiver) at 400 steps.
    /// //1e-4 leaves ~8x headroom without being so loose that a real disagreement hides in it.
    /// assert!((analytic - tree).abs() < 1e-4, "{analytic} vs {tree}");
    /// ```
    #[allow(clippy::too_many_arguments)]
    pub fn european_payer_swaption_tree(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64, //tenor of simple yield
        swap_rate: f64,
        num_steps: usize,
    ) -> Result<f64, HullWhiteError> {
        let spec = SwaptionSpec {
            r_t,
            t,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            is_payer: true,
            num_steps,
        };
        validate_swaption(&spec)?;
        Ok(self.swaption_tree(&spec, EUROPEAN))
    }
    /// Returns the price of a European payer swaption at current time from the short-rate tree.
    ///
    /// Exactly [`HullWhite::european_payer_swaption_tree`] at `t = 0` with `r_t = r(0)`
    /// ([`HullWhite::short_rate_now`]), so the tree starts on the model's initial short rate.
    ///
    /// # Comments
    ///
    /// A verification helper like [`HullWhite::european_payer_swaption_tree`]: prefer the
    /// closed-form [`HullWhite::european_payer_swaption_now`], and use this to check it against
    /// the tree.
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2;
    /// let sigma = 0.03;
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let option_maturity = 1.5;
    /// let num_swap_payments = 20;
    /// let delta = 0.25;
    /// let r0 = hull_white.short_rate_now().unwrap();
    /// let swap_rate = hull_white
    ///     .forward_swap_rate_t(r0, 0.0, option_maturity, num_swap_payments, delta)
    ///     .unwrap();
    /// let now = hull_white
    ///     .european_payer_swaption_tree_now(
    ///         option_maturity, num_swap_payments, delta, swap_rate, 200
    ///     )
    ///     .unwrap();
    /// let via_t = hull_white
    ///     .european_payer_swaption_tree(
    ///         r0, 0.0, option_maturity, num_swap_payments, delta, swap_rate, 200
    ///     )
    ///     .unwrap();
    /// assert!((now - via_t).abs() < 1e-12, "{now} vs {via_t}");
    /// assert!(now > 0.0, "european payer swaption tree {now}");
    /// ```
    pub fn european_payer_swaption_tree_now(
        &self,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
        num_steps: usize,
    ) -> Result<f64, HullWhiteError> {
        let t = 0.0; //since "now"
        let r_t = self.short_rate_now()?;
        self.european_payer_swaption_tree(
            r_t,
            t,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            num_steps,
        )
    }
    /// Prices a European receiver swaption on the short-rate tree, as an independent check on
    /// [`HullWhite::european_receiver_swaption_t`].
    ///
    /// # Comments
    ///
    /// The receiver-side twin of [`HullWhite::european_payer_swaption_tree`], and the same
    /// proposition: a verification helper walking the tree with early exercise off, to be
    /// cross-checked against the closed-form Jamshidian price rather than used in its place.
    ///
    /// # Examples
    ///
    /// ```
    /// let curr_rate = 0.02;
    /// let a = 0.2;
    /// let sigma = 0.03;
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let t = 0.5;
    /// let option_maturity = 1.5;
    /// let num_swap_payments = 20;
    /// let delta = 0.25;
    /// let swap_rate = hull_white
    ///     .forward_swap_rate_t(curr_rate, t, option_maturity, num_swap_payments, delta)
    ///     .unwrap();
    /// let analytic = hull_white
    ///     .european_receiver_swaption_t(
    ///         curr_rate, t, option_maturity, num_swap_payments, delta, swap_rate
    ///     )
    ///     .unwrap();
    /// let tree = hull_white
    ///     .european_receiver_swaption_tree(
    ///         curr_rate, t, option_maturity, num_swap_payments, delta, swap_rate, 400
    ///     )
    ///     .unwrap();
    /// assert!((analytic - tree).abs() < 1e-4, "{analytic} vs {tree}");
    /// ```
    #[allow(clippy::too_many_arguments)]
    pub fn european_receiver_swaption_tree(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64, //tenor of simple yield
        swap_rate: f64,
        num_steps: usize,
    ) -> Result<f64, HullWhiteError> {
        let spec = SwaptionSpec {
            r_t,
            t,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            is_payer: false,
            num_steps,
        };
        validate_swaption(&spec)?;
        Ok(self.swaption_tree(&spec, EUROPEAN))
    }
    /// Returns the price of a European receiver swaption at current time from the short-rate tree.
    ///
    /// Exactly [`HullWhite::european_receiver_swaption_tree`] at `t = 0` with `r_t = r(0)`
    /// ([`HullWhite::short_rate_now`]), so the tree starts on the model's initial short rate.
    ///
    /// # Comments
    ///
    /// A verification helper like [`HullWhite::european_receiver_swaption_tree`]: prefer the
    /// closed-form [`HullWhite::european_receiver_swaption_now`], and use this to check it
    /// against the tree.
    ///
    /// # Examples
    ///
    /// ```
    /// let a = 0.2;
    /// let sigma = 0.03;
    /// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
    /// let hull_white = hull_white::HullWhite::new(a, sigma, &curve).unwrap();
    /// let option_maturity = 1.5;
    /// let num_swap_payments = 20;
    /// let delta = 0.25;
    /// let r0 = hull_white.short_rate_now().unwrap();
    /// let swap_rate = hull_white
    ///     .forward_swap_rate_t(r0, 0.0, option_maturity, num_swap_payments, delta)
    ///     .unwrap();
    /// let now = hull_white
    ///     .european_receiver_swaption_tree_now(
    ///         option_maturity, num_swap_payments, delta, swap_rate, 200
    ///     )
    ///     .unwrap();
    /// let via_t = hull_white
    ///     .european_receiver_swaption_tree(
    ///         r0, 0.0, option_maturity, num_swap_payments, delta, swap_rate, 200
    ///     )
    ///     .unwrap();
    /// assert!((now - via_t).abs() < 1e-12, "{now} vs {via_t}");
    /// assert!(now > 0.0, "european receiver swaption tree {now}");
    /// ```
    pub fn european_receiver_swaption_tree_now(
        &self,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
        num_steps: usize,
    ) -> Result<f64, HullWhiteError> {
        let t = 0.0; //since "now"
        let r_t = self.short_rate_now()?;
        self.european_receiver_swaption_tree(
            r_t,
            t,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
            num_steps,
        )
    }
}

#[cfg(test)]
mod tests;
