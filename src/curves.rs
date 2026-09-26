//! Curve maths: the building blocks every pricer integrates, differentiates or discounts with.
//!
//! * [`a_t`] / [`at_t`] — the Hull-White (Vasicek) bond "duration" `B(t, T)`, the coefficient
//!   that multiplies the rate in an exponent.
//! * [`ct_t`] — the affine constant `C(t, T)` of the zero-coupon price
//!   `P(t, T) = exp(C(t, T) - B(t, T) * r_t)`, fitted to the initial yield and forward curves.
//! * [`gamma_edf`] — the variance integral used by the Eurodollar-futures convexity term.
//! * [`edf_compute`] / [`compute_libor_rate`] — the small algebraic transforms from bond prices
//!   to a futures price and a simple (Libor) rate.

//tdiff=T-t
pub(crate) fn a_t(a: f64, t_diff: f64) -> f64 {
    (1.0 - (-a * t_diff).exp()) / a
}

//t is first future time
//t_m is second future time
pub(crate) fn at_t(a: f64, t: f64, t_m: f64) -> f64 {
    a_t(a, t_m - t)
}

pub(crate) fn ct_t(
    a: f64,
    sigma: f64,
    t: f64,
    t_m: f64,
    yield_curve: &dyn Fn(f64) -> f64,
    forward_curve: &dyn Fn(f64) -> f64,
) -> f64 {
    let sqr = (-a * t_m).exp() - (-a * t).exp();
    yield_curve(t) - yield_curve(t_m) + forward_curve(t) * at_t(a, t, t_m)
        - (sigma * sqr).powi(2) * ((2.0 * a * t).exp() - 1.0) / (4.0 * a.powi(3))
}

//https://www.math.nyu.edu/~alberts/spring07/Lecture5.pdf
//https://developers.opengamma.com/quantitative-research/Hull-White-One-Factor-Model-OpenGamma.pdf (note that in the open gamma derivation, t0=option_maturity)
pub(crate) fn gamma_edf(a: f64, sigma: f64, t: f64, option_maturity: f64, delta: f64) -> f64 {
    let exp_t = (-a * (option_maturity - t)).exp();
    let exp_d = (-a * delta).exp();
    (sigma.powi(2) / a.powi(3))
        * (1.0 - exp_d)
        * ((1.0 - exp_t) - exp_d * 0.5 * (1.0 - exp_t.powi(2)))
}
pub(crate) fn edf_compute(bond_num: f64, bond_den: f64, gamma: f64, delta: f64) -> f64 {
    ((bond_num / bond_den) * gamma.exp() - 1.0) / delta
}

pub(crate) fn compute_libor_rate(nearest_bond: f64, farthest_bond: f64, tenor: f64) -> f64 {
    (nearest_bond - farthest_bond) / (farthest_bond * tenor)
}
