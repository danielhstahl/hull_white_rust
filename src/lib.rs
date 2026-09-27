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
//! All times are measured with respect to time 0 ("now"), and every time parameter is named for
//! the role it plays rather than for a letter of the model notation:
//!
//! - `t` — the valuation date, the time the state `r_t` is observed at; the `*_now` variants
//!   are the `t = 0` case of the same instrument, with the state read off the calibration
//!   instead of passed in
//! - `option_maturity` — when an option is exercised, and the maturity of the numeraire its
//!   price is written under
//! - `bond_maturity` — when the bond an option delivers is repaid; always later than the
//!   option's own maturity, or there is nothing left to deliver
//! - `delta` — the tenor of a simple (Libor) rate, the length of one payment period
//!
//! Textbook symbols (`T`, `T_b`, `t0`) show up only inside formulas, where the doc comment
//! spells each one against the parameter name it stands for.  The older `(0, t, T, TM)`
//! shorthand is deliberately not used in signatures: `T` and `TM` never said which maturity was
//! the option's and which was the deliverable's, which is exactly how `t_forward_bond_vol` came
//! to carry documentation the reverse of its own formula.
//!
//! ## Example Usage
//!
//! ```rust
//! use hull_white::{HullWhite, from_yield};
//!
//! // One curve object: the cumulative yield y(t) = 0.05 t + 0.01 t^2, whose derivative
//! // f(0,t) = 0.05 + 0.02 t is the instantaneous forward the model drifts to.
//! let curve = from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
//!
//! // Create a Hull-White model with parameters
//! let a = 0.1;      // Mean reversion speed
//! let sigma = 0.01; // Volatility parameter
//! let hull_white = HullWhite::new(a, sigma, &curve).unwrap();
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
//! ## The curve
//!
//! A model is calibrated to one [`YieldCurve`], not to a pair of closures.  The trait's one
//! required method is the **cumulative** yield `y(t)` (so `P(0,t) = exp(-y(t))`), and the
//! instantaneous forward `f(0,t) = y'(t)` is derived from it — or supplied in closed form with
//! [`from_yield_and_forward`] when a finite difference is not wanted.  Because the two halves are
//! one object, the relation between them is checkable, and [`HullWhite::new`] checks it: a curve
//! whose forward is not the derivative of its own yield is refused at construction instead of
//! mispricing quietly.  See [`curves`] for the trait, the tolerance, and why the convention here
//! is cumulative rather than annualised.
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
//! A coupon schedule may straddle the option's expiry date, because real bonds do.  The underlying
//! is the bond *as it stands on that date*, so payments made strictly before expiry are dropped
//! (the holder never receives them), a payment falling exactly on expiry is cash that folds into
//! the strike, and only payments strictly after expiry are decomposed.  A schedule with nothing
//! left after the expiry date is refused, since a bond settled before exercise is not a
//! deliverable underlying.  See the `jamshidian` module (`src/jamshidian.rs`) for the convention and
//! the reasoning behind it.
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
//! ## Pricing at `now`: the `_now` variants and `r(0)`
//!
//! Every pricing function comes in two shapes.  The `_t` shape takes the state explicitly — `r_t`,
//! the short rate observed at the valuation time `t` — and prices the instrument from there.  The
//! `_now` shape prices the same instrument at `t = 0`, where the state is not an input but a
//! consequence of the calibration, so it drops the `r_t` argument entirely.  The rate that stands
//! in for it is [`HullWhite::short_rate_now`], and each `now` variant of a state-dependent
//! instrument is exactly:
//!
//! ```rust
//! # use hull_white::{HullWhite, from_yield_and_forward};
//! # let curve = from_yield_and_forward(
//! #     |t: f64| 0.05 * t + 0.01 * t * t,
//! #     |t: f64| 0.05 + 0.02 * t,
//! # );
//! # let hull_white = HullWhite::new(0.2, 0.03, &curve).unwrap();
//! # let periods = [(1.0, 0.04), (1.25, 0.045), (1.5, 0.05)];
//! let r0 = hull_white.short_rate_now().unwrap();
//! let now = hull_white.cap_now(&periods, 0.25).unwrap();
//! let via_t = hull_white.cap_t(r0, 0.0, &periods, 0.25).unwrap();
//! assert!((now - via_t).abs() < 1e-12, "{now} vs {via_t}");
//! ```
//!
//! `r(0)` is the **instantaneous** forward at the front of the curve, `curve.forward(0.0)` — which
//! is the same thing as `phi(0)`, because the volatility term of `phi`,
//! `sigma^2 (1 - e^{-a t})^2 / (2 a^2)`, vanishes at `t = 0`.  It is *not* `curve.zero_yield(0.0)`:
//! the yield is cumulative (`bond_price_now(T) = exp(-zero_yield(T))` is the integral of the
//! forward), so `zero_yield(0.0)` is `0` for any curve worth the name.  A curve whose forward is
//! `-inf` at `0` — `|t| t.ln()`, which appears in the older examples in this crate — cannot price
//! anything on the `now` side, and says so at [`short_rate_now`](HullWhite::short_rate_now).
//!
//! The zero-coupon `_now` prices never look at a rate at all — `bond_price_now` is a closed form off
//! the yield curve — but a Jamshidian-priced option does, and not by accident: the decomposition
//! needs the rate at which the deliverable equals the strike, and the distribution of that rate is
//! anchored at `r(0)`.  Getting it from the curve rather than from a guess is what makes
//! `coupon_bond_call_now(r0-dependent args)` agree with `bond_call_now`'s discounting.
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
//! ## Crate layout
//!
//! The crate is split by instrument family so a numeric change touches one file, not all of them:
//!
//! | Module | What lives there |
//! |---|---|
//! | `model` | the [`HullWhite`] struct, calibration entry, `phi_t` / `mu_r` / `variance_r` / [`short_rate_now`](HullWhite::short_rate_now) / `t_forward_bond_vol` |
//! | [`curves`] | the [`YieldCurve`] trait and its builders, the construction-time consistency check, and the [`HullWhite`] affine bond coefficients `bond_b` (the duration) and `bond_c` (the price constant) |
//! | `schedules` | coupon/payment schedules and the remaining-payment count |
//! | `bonds` | zero coupon and coupon bond prices, and the coupon-sum kernels they share |
//! | `jamshidian` | the critical-rate bracket and solve, and the decomposition that consumes them |
//! | `options` | bond options and coupon-bond option entry points |
//! | `rates` | caplets, whole caps and floors over a period schedule, Eurodollar futures (with the convexity `gamma_edf`), forward and spot Libor |
//! | `swaps` | forward swap rate, swap price, European swaptions |
//! | `trees` | the short-rate tree: the shared lattice engine (`tree_price`), the European tree cross-check (`european_*_swaption_tree`) and the American swaptions |
//! | [`error`] | [`error::HullWhiteError`] |
//! | `validation` | the input contract every public entry point is checked against |
//! | `rootfinder` | the bracketed, safeguarded scalar root solver |
//! | `test_support` | *test / bench only, `#[doc(hidden)]`* — the shared HW-consistent `HwCurve` implementation of [`YieldCurve`], and the named calibrations (`flat_5pct`, `steep_curve`, `low_vol`, ...) every test and bench prices off |
//! | `mc` | *test only* — the Monte-Carlo harness: per-test seeds, antithetic variates, the one time-grid convention, and the `k * standard_error + discretisation bias` error budget every simulated price is asserted against |
//!
//! All public items are re-exported here, so `hull_white::HullWhite`, `hull_white::get_coupon_times`
//! and the rest of the surface resolve from the crate root regardless of which module defines them.
//! Each module keeps its tests in a sibling `tests.rs` (`src/bonds/tests.rs`, ...); the shared
//! curve fixtures every module prices off live in `test_support`, compiled under `cfg(test)` and
//! under the hidden `test-support` feature so `benches/` — a separate crate, which cannot see
//! `cfg(test)` — builds against the same fixture rather than a copy of it.  The Monte-Carlo
//! harness in `mc` is `cfg(test)` only, because it needs the `rand` dev-dependency, and its
//! sample size is tiered on the `slow` feature so that plain `cargo test` stays quick; `src/mc.rs`
//! documents what `--features slow` buys and how to run it.  Outside the crate,
//! `tests/public_api.rs` walks the whole public surface the way a consumer does: that is the check
//! that the split moved code without changing what `hull_white::` resolves to.

pub mod error;
// Monte-Carlo test harness (seeds, antithetic pairing, error budgets).  Test-only: it needs the
// `rand` dev-dependency, so it is not reachable from a non-test build of the crate.
#[cfg(test)]
mod mc;
mod rootfinder;
pub use rootfinder::{Solution, SolverError, SolverSettings};
#[cfg(test)]
mod solver_ab;
mod validation;

mod bonds;
pub mod curves;
mod jamshidian;
mod model;
mod options;
mod rates;
mod schedules;
mod swaps;
mod trees;

// Shared fixtures for the tests and the benches: `#[doc(hidden)]` and behind a non-default
// feature, so it is a build-time helper and not part of `hull_white::`'s public surface.
#[doc(hidden)]
#[cfg(any(test, feature = "test-support"))]
pub mod test_support;

pub use curves::{
    CURVE_PROBE_TIMES, FORWARD_CONSISTENCY_RELATIVE_TOLERANCE, FORWARD_CONSISTENCY_TOLERANCE,
    ForwardInconsistency, FromYield, YieldAndForward, YieldCurve, forward_consistency_tolerance,
    from_yield, from_yield_and_forward, max_forward_inconsistency, validate_curve,
};
pub use model::HullWhite;
pub use schedules::get_coupon_times;
