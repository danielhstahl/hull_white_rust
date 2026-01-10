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
//! use hull_white::{HullWhite, Result};
//!
//! // Define yield and forward curves
//! let yield_curve = |t: f64| 0.05 * t;  // Simple linear yield curve
//! let forward_curve = |t: f64| t.ln();  // Natural log forward curve
//!
//! // Create a Hull-White model with parameters
//! let a = 0.1;      // Mean reversion speed
//! let sigma = 0.01; // Volatility parameter
//! let hull_white = HullWhite::new_panicking(a, sigma, &yield_curve, &forward_curve);
//!
//! // Price a zero-coupon bond maturing in 2 years
//! let bond_price = hull_white.bond_price_now(2.0);
//! println!("Bond price: {}", bond_price);
//!
//! // Price a call option on a bond with 2-year maturity, expiring in 1 year, with strike 0.95
//! let option_price = hull_white.bond_call_now(1.0, 2.0, 0.95);
//! println!("Call option price: {}", option_price);
//! ```
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

pub mod constants;
pub mod error;
pub mod model;
pub mod products;
pub mod utils;

pub use error::HullWhiteError;
pub use model::HullWhite;

/// Type alias for result with HullWhiteError
pub type Result<T> = std::result::Result<T, HullWhiteError>;

#[cfg(test)]
mod test_utils;