//! Consumer-level check of the crate's public surface.
//!
//! The crate is split into one module per instrument family (see the crate docs in `src/lib.rs`).
//! A split is only safe if nothing a downstream user wrote has to change, and that cannot be
//! checked from inside the crate: an item can be `pub` in `src/swaps.rs` and still be invisible
//! from outside if `src/lib.rs` forgets to re-export it.  This file compiles as a separate crate
//! against `hull_white`, so every path below is exactly what a consumer writes, on the stable
//! toolchain — which `benches/bench.rs`, the other end-to-end compile check, is not (it needs
//! `#![feature(test)]`).
//!
//! Two things are pinned here:
//!
//! 1. Every public item that existed before the split still resolves from the crate root:
//!    `hull_white::HullWhite`, `hull_white::get_coupon_times`, the root-finding types
//!    (`SolverSettings`, `Solution`, `SolverError`) and `hull_white::error::HullWhiteError`.
//! 2. Every public method is callable with its pre-split argument list and returns a sane value.
//!    Argument *order* is part of the API, and a swap of two `f64`s compiles happily while
//!    pricing the wrong thing — so each call below uses distinct values per argument and checks
//!    the result against the economics it is supposed to describe.
//!
//! These are shape checks, not pricing checks: the numerics are covered by the per-module tests
//! under `src/`.

use hull_white::error::HullWhiteError;
use hull_white::{HullWhite, Solution, SolverError, SolverSettings, get_coupon_times};

type Model = HullWhite<'static, fn(f64) -> f64, fn(f64) -> f64>;

fn yield_curve(t: f64) -> f64 {
    0.05 * t
}
fn forward_curve(t: f64) -> f64 {
    t.ln()
}

/// `HullWhite<'a, T, U>` *borrows* its curves, so a consumer that wants a `'static` model has to
/// hand it something with `'static` type.  `fn` items are distinct from `fn` pointers, hence the
/// statics rather than a pair of locals.
static YIELD_CURVE: fn(f64) -> f64 = yield_curve;
static FORWARD_CURVE: fn(f64) -> f64 = forward_curve;

fn model() -> Model {
    HullWhite::init(0.1, 0.01, &YIELD_CURVE, &FORWARD_CURVE).unwrap()
}

/// A schedule far enough away that nothing has been paid yet at `t = 1.0`.
fn coupon_times() -> Vec<f64> {
    get_coupon_times(4, 1.0, 0.25).unwrap()
}

/// A schedule whose first coupon is strictly after the 1.5y expiry, as an option on a coupon bond
/// requires: a coupon paid before expiry is not part of the exercise decision.
fn post_expiry_coupon_times() -> Vec<f64> {
    get_coupon_times(4, 1.5, 0.25).unwrap()
}

// ---- type-level surface -------------------------------------------------------------

#[test]
fn public_types_resolve_from_crate_root() {
    fn takes_model<T: Fn(f64) -> f64 + Sync, U: Fn(f64) -> f64 + Sync>(_: &HullWhite<'_, T, U>) {}
    fn takes_error(_: HullWhiteError) {}
    fn takes_solution(_: Solution) {}
    fn takes_solver_error(_: SolverError) {}
    fn takes_settings(_: SolverSettings) {}

    let hw = model();
    takes_model::<fn(f64) -> f64, fn(f64) -> f64>(&hw);
    takes_error(HullWhiteError::InvalidInput("shape check".to_string()));
    takes_settings(hw.solver_settings());
    takes_solution(Solution {
        root: 0.05,
        iterations: 4,
        residual: 1e-14,
    });
    takes_solver_error(SolverError::Exhausted {
        iterations: 100,
        lower: 0.0,
        upper: 1.0,
        residual: 1e-6,
    });

    let _: fn(usize, f64, f64) -> Result<Vec<f64>, HullWhiteError> = get_coupon_times;
}

#[test]
fn solver_settings_are_public_and_round_trip() {
    let settings = SolverSettings {
        tolerance: 1e-11,
        max_iterations: 250,
        initial_guess: None,
    };
    let tuned = model().with_solver(settings).unwrap();
    assert_eq!(tuned.solver_settings().tolerance, 1e-11);
    assert_eq!(tuned.solver_settings().max_iterations, 250);
    // The default is still reachable and still tight.
    assert!(model().solver_settings().tolerance <= 1e-10);
    // A bad budget is rejected at the boundary, not at pricing time.
    assert!(
        model()
            .with_solver(SolverSettings {
                tolerance: 0.0,
                ..Default::default()
            })
            .is_err()
    );
}

// ---- schedules ----------------------------------------------------------------------

#[test]
fn get_coupon_times_generates_the_schedule() {
    assert_eq!(
        get_coupon_times(4, 1.0, 0.25).unwrap(),
        vec![1.25, 1.5, 1.75, 2.0]
    );
    // No payments is an empty schedule, not an error; a zero period is.
    assert_eq!(get_coupon_times(0, 1.0, 0.25).unwrap(), Vec::<f64>::new());
    assert!(get_coupon_times(4, 1.0, 0.0).is_err());
    assert!(get_coupon_times(4, -1.0, 0.25).is_err());
}

// ---- model primitives ---------------------------------------------------------------

#[test]
fn model_primitive_calls() {
    let hw = model();
    // sigma * sqrt((1 - exp(-2a(t_m - t))) / (2a)) * (1 - exp(-a(t_f - t_m))) > 0
    let vol = hw.t_forward_bond_vol(1.0, 2.0, 3.0).unwrap();
    assert!(vol > 0.0 && vol.is_finite(), "bond vol {vol}");
    assert!(hw.mu_r(0.05, 1.0, 2.0).unwrap().is_finite());
    assert!(hw.variance_r(1.0, 2.0).unwrap() > 0.0);
}

// ---- bonds --------------------------------------------------------------------------

#[test]
fn bond_price_calls() {
    let hw = model();
    let p = hw.bond_price_t(0.05, 1.0, 3.0).unwrap();
    assert!(p > 0.0 && p < 1.0, "zero coupon price {p}");
    let p_now = hw.bond_price_now(3.0).unwrap();
    assert!(p_now > 0.0 && p_now < 1.0, "zero coupon now price {p_now}");

    let times = coupon_times();
    let cp = hw.coupon_bond_price_t(0.05, 1.0, &times, 0.05).unwrap();
    // Coupon bond = coupon stream + par at maturity, so it is worth more than the bare discount bond.
    assert!(cp > p, "coupon bond {cp} vs discount bond {p}");
    let cp_now = hw.coupon_bond_price_now(&times, 0.05).unwrap();
    assert!(
        cp_now > p_now,
        "coupon bond now {cp_now} vs discount bond now {p_now}"
    );
}

// ---- zero coupon bond options -------------------------------------------------------

#[test]
fn bond_option_calls() {
    let hw = model();
    let call = hw.bond_call_t(0.05, 1.0, 1.5, 3.0, 0.9).unwrap();
    let put = hw.bond_put_t(0.05, 1.0, 1.5, 3.0, 0.9).unwrap();
    assert!(call > 0.0 && put > 0.0);
    // Put-call parity on the same strike: C - K*P(0; u, T) == Put? use monotonicity instead:
    // a deeper strike can never be worth more for a call, less for a put.
    let call_low_k = hw.bond_call_t(0.05, 1.0, 1.5, 3.0, 0.8).unwrap();
    let put_high_k = hw.bond_put_t(0.05, 1.0, 1.5, 3.0, 1.0).unwrap();
    assert!(
        call_low_k > call,
        "call is increasing in strike? {call_low_k} vs {call}"
    );
    assert!(
        put_high_k > put,
        "put is decreasing in strike? {put_high_k} vs {put}"
    );

    assert!(hw.bond_call_now(1.5, 3.0, 0.9).unwrap() > 0.0);
    assert!(hw.bond_put_now(1.5, 3.0, 0.9).unwrap() > 0.0);
}

#[test]
fn coupon_bond_option_calls() {
    let hw = model();
    let times = post_expiry_coupon_times();
    // The underlying coupon bond is worth ~1.045 here, so K = 0.95 is in the money and K = 1.1 is
    // out of the money for the call / in the money for the put.
    let call = |k: f64| {
        hw.coupon_bond_call_t(0.05, 1.0, 1.5, &times, 0.05, k)
            .unwrap()
    };
    let put = |k: f64| {
        hw.coupon_bond_put_t(0.05, 1.0, 1.5, &times, 0.05, k)
            .unwrap()
    };
    assert!(call(0.95) > 0.0, "ITM coupon bond call {}", call(0.95));
    assert!(put(1.1) > 0.0, "OTM-strike put {}", put(1.1));
    // Strike monotonicity, across the whole Jamshidian decomposition: a cheaper strike is a
    // better call, a richer strike is a better put.
    assert!(call(0.85) > call(0.95), "call decreasing in strike");
    assert!(put(1.2) > put(1.1), "put increasing in strike");
    // Jamshidian is exact, so the decomposition has to agree with the plain call/put relation
    // C - P = P_c(r) - K on the same underlying (checked to model precision by the module tests;
    // here only that both legs stay finite and non-negative).
    assert!(call(1.0) >= 0.0 && put(1.0) >= 0.0);
}

// ---- simple-rate instruments -------------------------------------------------------

#[test]
fn rate_instrument_calls() {
    let hw = model();
    assert!(hw.caplet_now(1.0, 0.25, 0.05).unwrap() > 0.0);
    assert!(hw.caplet_t(0.05, 0.5, 1.0, 0.25, 0.05).unwrap() > 0.0);
    assert!(
        hw.euro_dollar_future_t(0.05, 0.5, 1.0, 0.25)
            .unwrap()
            .is_finite()
    );
    assert!(hw.euro_dollar_future_now(1.0, 0.25).unwrap().is_finite());
    let fwd = hw.forward_libor_rate_t(0.05, 0.5, 1.0, 0.25).unwrap();
    let spot = hw.libor_rate_t(0.05, 0.5, 0.25).unwrap();
    assert!(fwd.is_finite() && spot.is_finite());
    assert!(hw.forward_libor_rate_now(1.0, 0.25).unwrap().is_finite());
}

// ---- swaps and swaptions -----------------------------------------------------------

#[test]
fn swap_calls() {
    let hw = model();
    let rate = hw.forward_swap_rate_t(0.05, 1.0, 1.0, 4, 0.25).unwrap();
    assert!(rate.is_finite());
    assert!(hw.swap_rate_t(0.05, 1.0, 4, 0.25).unwrap().is_finite());
    // At the forward swap rate the swap is worth zero; above it the payer pays.
    assert!(hw.swap_price_t(0.05, 1.0, 2.0, 0.25, rate).unwrap().abs() < 1e-6);
    assert!(hw.swap_price_t(0.05, 1.0, 2.0, 0.25, rate + 0.01).unwrap() < 0.0);
    assert!(
        hw.swap_price_t_init(0.05, 1.0, 1.0, 4, 0.25, rate)
            .unwrap()
            .abs()
            < 1e-6
    );
}

#[test]
fn swaption_calls() {
    let hw = model();
    // Anchor on the forward swap rate *at the valuation time*; a strike is only in or out of the
    // money relative to that, so anchoring anywhere else makes the assertions meaningless.
    let atm = hw.forward_swap_rate_t(0.05, 0.5, 1.0, 4, 0.25).unwrap();
    let eur_payer = |k: f64| {
        hw.european_payer_swaption_t(0.05, 0.5, 1.0, 4, 0.25, k)
            .unwrap()
    };
    let eur_receiver = |k: f64| {
        hw.european_receiver_swaption_t(0.05, 0.5, 1.0, 4, 0.25, k)
            .unwrap()
    };
    let am_payer = |k: f64| {
        hw.american_payer_swaption_t(0.05, 0.5, 1.0, 4, 0.25, k, 40)
            .unwrap()
    };
    let am_receiver = |k: f64| {
        hw.american_receiver_swaption_t(0.05, 0.5, 1.0, 4, 0.25, k, 40)
            .unwrap()
    };
    assert!(eur_payer(atm) > 0.0 && eur_receiver(atm) > 0.0);
    // A payer wants to pay below the forward; a receiver wants to receive above it.
    assert!(eur_payer(atm - 0.01) > eur_payer(atm));
    assert!(eur_receiver(atm + 0.01) > eur_receiver(atm));
    assert!(am_payer(atm) > 0.0 && am_receiver(atm) > 0.0);
    // American prices stay monotone in the strike too.  They are *not* compared against the
    // European numbers here: the tree discretises the same payoff on `num_steps` buckets while the
    // analytic pricer integrates a Black approximation, so the two only agree in the limit.
    // That comparison is the tree's own job — see `european_swaption_tree` in `src/trees.rs`.
    assert!(am_payer(atm - 0.01) > am_payer(atm));
    assert!(am_receiver(atm + 0.01) > am_receiver(atm));
}

// ---- errors reach the consumer -----------------------------------------------------

#[test]
fn invalid_instruments_are_errors_not_panics() {
    let hw = model();
    // Valuation time before now, an expired option, an empty schedule: each is InvalidInput.
    let bad = hw.bond_price_now(-1.0).unwrap_err();
    assert!(matches!(bad, HullWhiteError::InvalidInput(_)), "{bad:?}");
    let bad = hw.bond_call_now(0.0, 1.0, 0.9).unwrap_err();
    assert!(matches!(bad, HullWhiteError::InvalidInput(_)), "{bad:?}");
    let bad = hw.coupon_bond_price_now(&[], 0.05).unwrap_err();
    assert!(matches!(bad, HullWhiteError::InvalidInput(_)), "{bad:?}");
}
