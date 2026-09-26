//! Shared fixtures for the crate's unit tests.
//!
//! [`hw_curves`] builds a yield/forward-curve pair that is *internal* to a Hull-White model
//! (a Vasicek short rate reverting to `b`), so prices computed from it have an independent
//! closed form to check against.  The `STEEP_*` fixture is the deliberately time-inhomogeneous
//! calibration: the short rate sits far below the long-run mean, so `phi(t)` moves a lot and a
//! time-coordinate error cannot cancel out.
//!
//! Lives at the crate root (rather than in one module) because every module's tests use it.

use rand::SeedableRng;
use rand::StdRng;

/// A deliberately steep (strongly time-inhomogeneous) calibration: the short rate sits far below
/// the long-run mean, so phi(t) moves a lot over the option's life. The existing fixture uses
/// curr_rate == b, which makes phi(t) exactly constant and hides any time-coordinate error.
pub(crate) const STEEP_CURR_RATE: f64 = 0.02;
pub(crate) const STEEP_SIG: f64 = 0.03;
pub(crate) const STEEP_A: f64 = 0.2;
pub(crate) const STEEP_B: f64 = 0.06;

/// Curves consistent with a Hull-White process (Vasicek short rate with mean reversion `a`
/// towards `b`), so that `bond_price_now`/`bond_price_t` are exact for the same model.
pub(crate) fn hw_curves(
    curr_rate: f64,
    a: f64,
    b: f64,
    sig: f64,
) -> (impl Fn(f64) -> f64, impl Fn(f64) -> f64) {
    let yield_curve = move |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = move |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    (yield_curve, forward_curve)
}

pub(crate) fn get_rng_seed(seed: [u8; 32]) -> StdRng {
    SeedableRng::from_seed(seed)
}

/// The steeply calibrated model every invalid-input test hangs its assertions off.
///
/// Every path inside is `$crate`-qualified so a caller only needs `use
/// $crate::testutil::yvf_setup;` -- the fixture stays in one place and cannot be
/// shadowed at the call site.
macro_rules! yvf_setup {
    ($model:ident) => {
        let (yield_curve, forward_curve) = $crate::testutil::hw_curves(
            $crate::testutil::STEEP_CURR_RATE,
            $crate::testutil::STEEP_A,
            $crate::testutil::STEEP_B,
            $crate::testutil::STEEP_SIG,
        );
        let $model = $crate::HullWhite::init(
            $crate::testutil::STEEP_A,
            $crate::testutil::STEEP_SIG,
            &yield_curve,
            &forward_curve,
        )
        .unwrap();
    };
}

pub(crate) use yvf_setup;
