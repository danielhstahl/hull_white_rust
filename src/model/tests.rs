//! Unit tests for the model primitives.
//!
//! `steep_fixture_actually_makes_phi_time_dependent` guards the fixture the time-coordinate
//! regressions elsewhere rely on: if `phi` were constant, a wrong clock would cancel out and
//! those tests would pass for the wrong reason.

use crate::HullWhite;
use crate::testutil::{STEEP_A, STEEP_B, STEEP_CURR_RATE, STEEP_SIG, hw_curves};

#[test]
fn steep_fixture_actually_makes_phi_time_dependent() {
    //Protects the tests above: if phi were constant, a time-coordinate error would cancel out
    //and the regression tests would pass for the wrong reason.
    let (yield_curve, forward_curve) = hw_curves(STEEP_CURR_RATE, STEEP_A, STEEP_B, STEEP_SIG);
    let hull_white = HullWhite::init(STEEP_A, STEEP_SIG, &yield_curve, &forward_curve).unwrap();
    let phi_0 = hull_white.phi_t(0.0);
    let phi_5 = hull_white.phi_t(5.0);
    assert!(
        (phi_5 - phi_0).abs() > 0.005,
        "steep fixture must make phi(t) time dependent: phi(0)={phi_0}, phi(5)={phi_5}"
    );
}
