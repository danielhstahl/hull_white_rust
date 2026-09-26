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
//! Three things are pinned here:
//!
//! 1. Every public item that existed before the split still resolves from the crate root:
//!    `hull_white::HullWhite`, `hull_white::get_coupon_times`, the root-finding types
//!    (`SolverSettings`, `Solution`, `SolverError`) and `hull_white::error::HullWhiteError`.
//! 2. Every public method is callable with its pre-split argument list and returns a sane value.
//!    Argument *order* is part of the API, and a swap of two `f64`s compiles happily while
//!    pricing the wrong thing — so each call below uses distinct values per argument and checks
//!    the result against the economics it is supposed to describe.
//! 3. The `now` side behaves as documented: `HullWhite::short_rate_now` is the curve-derived `r(0)`,
//!    every `*_now` variant is its `*_t` twin at `(short_rate_now(), 0.0)`, and `cap` / `floor`
//!    aggregate a whole `(option_maturity, strike)` schedule into the sum of their caplets /
//!    floorlets.  See `every_now_variant_is_its_t_variant_at_zero` and
//!    `cap_and_floor_aggregate_periods`.
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

/// A curve pair whose *front* is finite.  The `now` side of the crate needs this: `r(0)` is
/// `forward_curve(0.0)`, and the `t.ln()` shape used above is `-inf` there.  It is also the shape
/// a real curve has — the instantaneous forward at the front of the curve is a number.
fn finite_yield_curve(t: f64) -> f64 {
    0.05 * t + 0.01 * t * t //the integral of 0.05 + 0.02 * t
}
fn finite_forward_curve(t: f64) -> f64 {
    0.05 + 0.02 * t
}

/// `HullWhite<'a, T, U>` *borrows* its curves, so a consumer that wants a `'static` model has to
/// hand it something with `'static` type.  `fn` items are distinct from `fn` pointers, hence the
/// statics rather than a pair of locals.
static YIELD_CURVE: fn(f64) -> f64 = yield_curve;
static FORWARD_CURVE: fn(f64) -> f64 = forward_curve;
static FINITE_YIELD_CURVE: fn(f64) -> f64 = finite_yield_curve;
static FINITE_FORWARD_CURVE: fn(f64) -> f64 = finite_forward_curve;

fn model() -> Model {
    HullWhite::init(0.1, 0.01, &YIELD_CURVE, &FORWARD_CURVE).unwrap()
}

/// The model to price the `now` side with: same shape, curves that are finite at `0`.
fn finite_model() -> Model {
    HullWhite::init(0.1, 0.01, &FINITE_YIELD_CURVE, &FINITE_FORWARD_CURVE).unwrap()
}

/// A schedule far enough away that nothing has been paid yet at `t = 1.0`.
fn coupon_times() -> Vec<f64> {
    get_coupon_times(4, 1.0, 0.25).unwrap()
}

/// A schedule starting strictly after the 1.5y expiry, so an option on it has nothing dropped.
fn post_expiry_coupon_times() -> Vec<f64> {
    get_coupon_times(4, 1.5, 0.25).unwrap()
}

/// The same bond as written at `t = 1.0`: coupons paid before the 1.5y expiry and one exactly on
/// it, so the deliverable underlying is the residual that starts at 1.75.
fn straddling_coupon_times() -> Vec<f64> {
    get_coupon_times(6, 1.0, 0.25).unwrap()
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

    // A bond's schedule may straddle the expiry: the underlying is the bond as it stands on the
    // expiry date, so the two coupons paid before it (1.25, 1.5 ...) drop out and the coupon
    // falling exactly on it arrives as a reduction of the strike.
    let straddle_call = |k: f64| {
        hw.coupon_bond_call_t(0.05, 1.0, 1.5, &straddling_coupon_times(), 0.05, k)
            .unwrap()
    };
    let residual_call = |k: f64| {
        hw.coupon_bond_call_t(0.05, 1.0, 1.5, &post_expiry_coupon_times(), 0.05, k)
            .unwrap()
    };
    assert!(straddle_call(0.95) > 0.0, "straddling call");
    assert_eq!(
        straddle_call(1.0),
        residual_call(1.0 - 0.05),
        "the coupon on the expiry date is a strike reduction on the residual bond"
    );
    // A bond with nothing left to deliver after the expiry date is not an underlying.
    assert!(matches!(
        hw.coupon_bond_call_t(0.05, 1.0, 3.0, &[1.25, 1.5, 2.5], 0.05, 1.0),
        Err(HullWhiteError::InvalidInput(_))
    ));
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

// ---- the `now` side ---------------------------------------------------------------

/// `r(0)` is public, and it is derived rather than supplied: the `now` side of the crate takes no
/// rate argument because the calibration already fixes what the short rate is today.
#[test]
fn short_rate_now_is_public_and_curve_derived() {
    let hw = finite_model();
    let r0 = hw.short_rate_now().unwrap();
    assert!(r0.is_finite(), "r(0) = {r0}");
    //The instantaneous forward at the front of the curve, not the cumulative yield there.
    assert!((r0 - finite_forward_curve(0.0)).abs() < 1e-12, "{r0}");
    assert!((r0 - 0.05).abs() < 1e-12, "{r0}");
    assert!(
        (r0 - finite_yield_curve(0.0)).abs() > 1e-3,
        "r(0) must not be read off the cumulative yield"
    );
}

/// The contract every `now`/`t` pair keeps, checked from outside the crate: the `now` variant is the
/// `t` variant at `(short_rate_now(), 0.0)`.  This pins shape *and* semantics — a pair wired to
/// the wrong argument still compiles, so each call below carries distinct values and compares the
/// two routes.
#[test]
fn every_now_variant_is_its_t_variant_at_zero() {
    let hw = finite_model();
    let r0 = hw.short_rate_now().unwrap();
    let close = |a: f64, b: f64, what: &str| {
        assert!(
            (a - b).abs() <= 1e-12 * a.abs().max(1.0),
            "{what}: now {a} vs t {b}"
        );
    };

    //Bonds.
    close(
        hw.bond_price_now(3.0).unwrap(),
        hw.bond_price_t(r0, 0.0, 3.0).unwrap(),
        "bond_price",
    );
    let coupon_times = get_coupon_times(4, 1.75, 0.25).unwrap();
    close(
        hw.coupon_bond_price_now(&coupon_times, 0.05).unwrap(),
        hw.coupon_bond_price_t(r0, 0.0, &coupon_times, 0.05)
            .unwrap(),
        "coupon_bond_price",
    );
    //Zero coupon bond options.
    close(
        hw.bond_call_now(1.5, 3.0, 0.9).unwrap(),
        hw.bond_call_t(r0, 0.0, 1.5, 3.0, 0.9).unwrap(),
        "bond_call",
    );
    close(
        hw.bond_put_now(1.5, 3.0, 0.9).unwrap(),
        hw.bond_put_t(r0, 0.0, 1.5, 3.0, 0.9).unwrap(),
        "bond_put",
    );
    //Coupon bond options (the Jamshidian side).
    close(
        hw.coupon_bond_call_now(1.5, &coupon_times, 0.05, 1.0)
            .unwrap(),
        hw.coupon_bond_call_t(r0, 0.0, 1.5, &coupon_times, 0.05, 1.0)
            .unwrap(),
        "coupon_bond_call",
    );
    close(
        hw.coupon_bond_put_now(1.5, &coupon_times, 0.05, 1.0)
            .unwrap(),
        hw.coupon_bond_put_t(r0, 0.0, 1.5, &coupon_times, 0.05, 1.0)
            .unwrap(),
        "coupon_bond_put",
    );
    //Caplets, floorlets and their aggregates.
    let delta = 0.25;
    let periods = [(1.0, 0.04), (1.25, 0.045), (1.5, 0.05)];
    close(
        hw.caplet_now(1.0, delta, 0.04).unwrap(),
        hw.caplet_t(r0, 0.0, 1.0, delta, 0.04).unwrap(),
        "caplet",
    );
    close(
        hw.floorlet_now(1.0, delta, 0.04).unwrap(),
        hw.floorlet_t(r0, 0.0, 1.0, delta, 0.04).unwrap(),
        "floorlet",
    );
    close(
        hw.cap_now(&periods, delta).unwrap(),
        hw.cap_t(r0, 0.0, &periods, delta).unwrap(),
        "cap",
    );
    close(
        hw.floor_now(&periods, delta).unwrap(),
        hw.floor_t(r0, 0.0, &periods, delta).unwrap(),
        "floor",
    );
    //Simple rates.
    close(
        hw.euro_dollar_future_now(1.0, delta).unwrap(),
        hw.euro_dollar_future_t(r0, 0.0, 1.0, delta).unwrap(),
        "euro_dollar_future",
    );
    close(
        hw.libor_rate_now(delta).unwrap(),
        hw.libor_rate_t(r0, 0.0, delta).unwrap(),
        "libor_rate",
    );
    //Swaps.
    close(
        hw.forward_swap_rate_now(1.0, 4, delta).unwrap(),
        hw.forward_swap_rate_t(r0, 0.0, 1.0, 4, delta).unwrap(),
        "forward_swap_rate",
    );
    close(
        hw.swap_rate_now(4, delta).unwrap(),
        hw.swap_rate_t(r0, 0.0, 4, delta).unwrap(),
        "swap_rate",
    );
    close(
        hw.swap_price_now(1.0, delta, 0.045).unwrap(),
        hw.swap_price_t(r0, 0.0, 1.0, delta, 0.045).unwrap(),
        "swap_price",
    );
    close(
        hw.swap_price_now_init(4, delta, 0.045).unwrap(),
        hw.swap_price_t_init(r0, 0.0, 0.0, 4, delta, 0.045).unwrap(),
        "swap_price_init",
    );
    //European swaptions.
    close(
        hw.european_payer_swaption_now(1.0, 4, delta, 0.045)
            .unwrap(),
        hw.european_payer_swaption_t(r0, 0.0, 1.0, 4, delta, 0.045)
            .unwrap(),
        "european_payer_swaption",
    );
    close(
        hw.european_receiver_swaption_now(1.0, 4, delta, 0.045)
            .unwrap(),
        hw.european_receiver_swaption_t(r0, 0.0, 1.0, 4, delta, 0.045)
            .unwrap(),
        "european_receiver_swaption",
    );
    //American swaptions, tree priced at the initial short rate.
    close(
        hw.american_payer_swaption_now(1.0, 4, delta, 0.045, 40)
            .unwrap(),
        hw.american_payer_swaption_t(r0, 0.0, 1.0, 4, delta, 0.045, 40)
            .unwrap(),
        "american_payer_swaption",
    );
    close(
        hw.american_receiver_swaption_now(1.0, 4, delta, 0.045, 40)
            .unwrap(),
        hw.american_receiver_swaption_t(r0, 0.0, 1.0, 4, delta, 0.045, 40)
            .unwrap(),
        "american_receiver_swaption",
    );
}

/// A cap or floor is one call for a whole schedule, and the schedule-level identities hold from
/// outside the crate: the aggregate equals the sum of its caplets (floorlets), and the cap/floor
/// gap is the forward legs.
#[test]
fn cap_and_floor_aggregate_periods() {
    let hw = finite_model();
    let delta = 0.25;
    let periods: Vec<(f64, f64)> = (1..=20)
        .map(|index| (index as f64 * delta, 0.04 + 0.001 * (index % 5) as f64))
        .collect();
    let cap = hw.cap_now(&periods, delta).unwrap();
    let floor = hw.floor_now(&periods, delta).unwrap();
    let caplets: f64 = periods
        .iter()
        .map(|&(option_maturity, strike)| hw.caplet_now(option_maturity, delta, strike).unwrap())
        .sum();
    let floorlets: f64 = periods
        .iter()
        .map(|&(option_maturity, strike)| hw.floorlet_now(option_maturity, delta, strike).unwrap())
        .sum();
    assert!(cap > 0.0 && floor > 0.0, "cap {cap}, floor {floor}");
    assert!(
        (cap - caplets).abs() <= 1e-12 * caplets.abs(),
        "cap {cap} vs sum of caplets {caplets}"
    );
    assert!(
        (floor - floorlets).abs() <= 1e-12 * floorlets.abs(),
        "floor {floor} vs sum of floorlets {floorlets}"
    );
    let legs: f64 = periods
        .iter()
        .map(|&(option_maturity, strike)| {
            let near = hw.bond_price_now(option_maturity).unwrap();
            let far = hw.bond_price_now(option_maturity + delta).unwrap();
            near - (1.0 + delta * strike) * far
        })
        .sum();
    assert!(
        (cap - floor - legs).abs() <= 1e-11 * legs.abs().max(1.0),
        "cap {cap} - floor {floor} vs forward legs {legs}"
    );
}

/// A curve with no finite front cannot supply `r(0)`, and says so instead of pricing a nonsense
/// number.  This is the trap the `|t| t.ln()` fixture in the older examples walks into: that curve
/// is `-inf` at `0`, so nothing on the state-dependent `now` side can be priced with it.
#[test]
fn an_infinite_front_curve_reports_r0_as_an_error() {
    let hw = model();
    assert!(matches!(
        hw.short_rate_now(),
        Err(HullWhiteError::NumericalError(_))
    ));
    assert!(
        hw.coupon_bond_call_now(1.5, &[1.75, 2.0], 0.05, 1.0)
            .is_err()
    );
    assert!(hw.swap_price_now(1.0, 0.25, 0.045).is_err());
    //The zero coupon `now` prices never consult a rate, so they still work.
    assert!(hw.bond_price_now(2.0).unwrap().is_finite());
    assert!(hw.caplet_now(1.0, 0.25, 0.04).unwrap() > 0.0);
}
