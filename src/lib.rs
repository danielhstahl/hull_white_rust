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
//! ## Crate layout
//!
//! The crate is split by instrument family so a numeric change touches one file, not all of them:
//!
//! | Module | What lives there |
//! |---|---|
//! | `model` | the [`HullWhite`] struct, calibration entry, `phi_t` / `mu_r` / `variance_r` / `t_forward_bond_vol` |
//! | `curves` | `a_t` (bond duration), `ct_t` (bond price constant), the Eurodollar variance integral |
//! | `schedules` | coupon/payment schedules and the remaining-payment count |
//! | `bonds` | zero coupon and coupon bond prices, and the coupon-sum kernels they share |
//! | `jamshidian` | the critical-rate bracket and solve, and the decomposition that consumes them |
//! | `options` | bond options and coupon-bond option entry points |
//! | `rates` | caplets, Eurodollar futures, forward and spot Libor |
//! | `swaps` | forward swap rate, swap price, European swaptions |
//! | `trees` | the short-rate tree: European tree check and American swaptions |
//! | [`error`] | [`error::HullWhiteError`] |
//! | `validation` | the input contract every public entry point is checked against |
//! | `rootfinder` | the bracketed, safeguarded scalar root solver |
//!
//! All public items are re-exported here, so `hull_white::HullWhite`, `hull_white::get_coupon_times`
//! and the rest of the surface resolve from the crate root regardless of which module defines them.
//! Each module keeps its tests in a sibling `tests.rs` (`src/bonds/tests.rs`, ...); the shared
//! fixtures used across modules are in `testutil`, compiled under `cfg(test)`.

pub mod error;
mod rootfinder;
pub use rootfinder::{Solution, SolverError, SolverSettings};
#[cfg(test)]
mod solver_ab;
mod validation;

mod bonds;
mod curves;
mod jamshidian;
mod model;
mod options;
mod rates;
mod schedules;
mod swaps;
#[cfg(test)]
mod testutil;
mod trees;

pub use model::HullWhite;
pub use schedules::get_coupon_times;
