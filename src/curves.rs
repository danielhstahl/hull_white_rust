//! The curve the model calibrates to, and the affine bond maths built on it.
//!
//! ## One curve object, one invariant
//!
//! A Hull-White model is calibrated to *one* initial term structure.  That object is
//! [`YieldCurve`]: a thing that answers two questions about the initial curve —
//! [`zero_yield`](YieldCurve::zero_yield), the **cumulative** yield `y(t)`, and
//! [`forward`](YieldCurve::forward), the **instantaneous** forward `f(0, t)`.  Only
//! `zero_yield` has to be implemented; everything else follows from it:
//!
//! ```text
//! discount(t) = exp(-zero_yield(t))          (provided: the zero coupon bond price)
//! forward(t)  = d/dt zero_yield(t)           (provided: a 5-point finite difference)
//! ```
//!
//! Overriding `forward` with the closed-form derivative is optional and cheap, and it is what
//! [`from_yield_and_forward`] exists for: it wraps a `(yield, forward)` pair as the single object
//! the model wants, so a caller who knows the derivative never pays for a finite difference.
//! [`from_yield`] is the one-closure case, where the forward really is derived.
//!
//! The relation between them — **`forward(t) == d/dt zero_yield(t)`** — is not a style
//! preference, it is the statement that the two numbers describe the *same* term structure.  Up
//! to 0.8.0 the model took the cumulative yield and the forward as two unrelated closures
//! (`HullWhite<'a, T, U>`), so a caller could hand it a pair that disagreed and every price came
//! out wrong with nothing said about it.  [`YieldCurve`] collapses the pair into one object and
//! [`HullWhite::new`](crate::HullWhite::new) checks the relation at construction, over the
//! probe grid in [`CURVE_PROBE_TIMES`] and to the tolerance in
//! [`forward_consistency_tolerance`].  An inconsistent curve is an
//! [`HullWhiteError::InvalidInput`] naming the time
//! it disagreed worst, the supplied value and the derived one.
//!
//! ### Which "yield" this is
//!
//! `zero_yield` is *cumulative*: `y(t) = integral of f(0, .) over [0, t]`, which is why
//! `discount(t)` is `exp(-y(t))` with no division by maturity.  It is **not** the annualised
//! zero rate `z(t) = y(t) / t` that a quoting screen prints.  Under the cumulative convention
//! the consistency relation is the plain `f = y'`; under the annualised one the same curve
//! satisfies `f(t) = z(t) + t z'(t)`.  Same curve, different function of it — so a
//! `YieldCurve` implemented against the annualised zero fails the check, loudly, rather than
//! pricing a curve that is a factor of `t` away from the one that was meant.
//!
//! ```
//! use hull_white::{HullWhite, YieldCurve};
//!
//! // A curve whose only primitive is the cumulative yield.  The forward is derived, and the
//! // model accepts it because the derived forward *is* its own derivative.
//! struct LinearYield(f64);
//! impl YieldCurve for LinearYield {
//!     fn zero_yield(&self, t: f64) -> f64 {
//!         self.0 * t
//!     }
//! }
//! let model = HullWhite::new(0.2, 0.03, &LinearYield(0.05)).unwrap();
//! assert!((model.bond_price_now(2.0).unwrap() - (-0.1f64).exp()).abs() < 1e-12);
//! // f(0,t) = d/dt (0.05 t) = 0.05 everywhere.
//! assert!((model.curve().forward(3.0) - 0.05).abs() < 1e-9);
//! ```

use crate::error::HullWhiteError;

/// The initial term structure a [`HullWhite`](crate::HullWhite) model is calibrated to.
///
/// Two accessors matter to the model, both at a time `t` measured from "now":
///
/// * [`zero_yield`](YieldCurve::zero_yield) — the **cumulative** yield `y(t)`, the integral of
///   the instantaneous forward over `[0, t]`.  This is the required method and the primitive:
///   the curve is fully described by it.
/// * [`forward`](YieldCurve::forward) — the **instantaneous** forward `f(0, t) = y'(t)`.  The
///   model's drift `phi(t)` is this plus a volatility term, so it is the number the short rate
///   mean-reverts towards, and it has to be the actual derivative of `zero_yield`.
///
/// [`discount`](YieldCurve::discount) — `P(0, t) = exp(-y(t))` — is provided in terms of
/// `zero_yield`, and `forward` is provided as a finite difference of it.  Overriding `forward`
/// with the analytic derivative costs one method and removes the finite difference from the
/// pricing path; do that for anything calibration-derived or speed-sensitive.
///
/// The invariant `forward == zero_yield'` is verified when a model is built — see
/// [`max_forward_inconsistency`] and [`HullWhite::new`](crate::HullWhite::new).
///
/// # Implementing it
///
/// Implement `zero_yield` and, if the derivative is known, `forward`:
///
/// ```
/// use hull_white::YieldCurve;
///
/// /// y(t) = 0.05 t + 0.01 t^2, so f(0,t) = 0.05 + 0.02 t.
/// struct QuadraticYield;
/// impl YieldCurve for QuadraticYield {
///     fn zero_yield(&self, t: f64) -> f64 {
///         0.05 * t + 0.01 * t * t
///     }
///     fn forward(&self, t: f64) -> f64 {
///         0.05 + 0.02 * t
///     }
/// }
///
/// let curve = QuadraticYield;
/// // y(5) = 0.25 + 0.25 = 0.5, so P(0,5) = e^-0.5
/// assert!((curve.discount(5.0) - (-0.5f64).exp()).abs() < 1e-15);
/// ```
///
/// Anything that is `Sync` and can answer for a `t` is a curve — a bootstrapped curve over a
/// vector of discount factors, a piecewise-linear forward, a function, a constant.  A curve that
/// cannot supply a finite `zero_yield` at any of the [`CURVE_PROBE_TIMES`] cannot calibrate a
/// model either, and is rejected at construction for that reason.
pub trait YieldCurve: Sync {
    /// The cumulative yield `y(t)` to time `t`, i.e. the integral of the instantaneous forward
    /// over `[0, t]`.  The single required piece of a curve.
    fn zero_yield(&self, t: f64) -> f64;

    /// The zero coupon bond price `P(0, t) = exp(-y(t))`.
    ///
    /// Override this (and `zero_yield`) if the curve is stored as discount factors, so the two
    /// accessors do not round-trip through a `ln`/`exp` pair.
    fn discount(&self, t: f64) -> f64 {
        (-self.zero_yield(t)).exp()
    }

    /// The instantaneous forward `f(0, t) = d/dt zero_yield(t)`.
    ///
    /// The default is a 5-point finite difference of [`zero_yield`](YieldCurve::zero_yield)
    /// (accurate to roughly `1e-11` absolute on rate-curve-shaped inputs — comfortably inside
    /// the [`forward_consistency_tolerance`] the construction check applies, and inside what a
    /// price is sensitive to).  Override with the analytic derivative when it is available.
    fn forward(&self, t: f64) -> f64 {
        derivative(&|s: f64| self.zero_yield(s), t)
    }
}

/// Absolute tolerance for the construction-time forward-consistency check, in annualised rate
/// units: `1e-6` is one hundredth of a basis point.
pub const FORWARD_CONSISTENCY_TOLERANCE: f64 = 1e-6;

/// Relative tolerance for the same check, applied to the magnitude of the supplied forward so a
/// 20% curve is not failed on a difference that is a rounding error at that scale.
pub const FORWARD_CONSISTENCY_RELATIVE_TOLERANCE: f64 = 1e-4;

/// The tolerance the consistency check allows at a probe point: the larger of
/// [`FORWARD_CONSISTENCY_TOLERANCE`] and [`FORWARD_CONSISTENCY_RELATIVE_TOLERANCE`] times the
/// magnitude of the supplied forward.
pub fn forward_consistency_tolerance(supplied: f64) -> f64 {
    FORWARD_CONSISTENCY_TOLERANCE.max(FORWARD_CONSISTENCY_RELATIVE_TOLERANCE * supplied.abs())
}

/// The times, in years from "now", the construction-time check probes the curve.
///
/// Six points from a quarter to a decade: enough span that a forward curve which is right at the
/// front and wrong further out — the usual way the two get edited apart — cannot pass, and few
/// enough that construction stays cheap on a curve whose evaluation is not.
///
/// `0.0` is deliberately absent.  The derivative there is one-sided, and the front of the curve
/// already has its own gate: [`HullWhite::short_rate_now`](crate::HullWhite::short_rate_now)
/// reports a curve that is not finite at `0` as an error rather than pricing `-inf` off it.
pub const CURVE_PROBE_TIMES: [f64; 6] = [0.25, 0.5, 1.0, 2.0, 5.0, 10.0];

/// The worst disagreement between a curve's own `forward` and the numerical derivative of its
/// `zero_yield`, over [`CURVE_PROBE_TIMES`].
///
/// A curve that supplies `forward` as the analytic derivative of its own `zero_yield` shows an
/// error at the finite-difference noise floor (`1e-12` or better on the fixtures in this crate);
/// a pair of closures that were written against different assumptions shows `1e-2` and up.
///
/// `None` means the check could not be run at all: the curve returned no finite
/// `zero_yield`+`forward` pair at any probe time.  That is treated as a failure to calibrate,
/// not as a pass.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ForwardInconsistency {
    /// The probe time that disagreed worst.
    pub t: f64,
    /// `curve.forward(t)`.
    pub supplied: f64,
    /// `d/dt curve.zero_yield(t)`, finite-differenced.
    pub derived: f64,
    /// `|supplied - derived|`.
    pub error: f64,
    /// [`forward_consistency_tolerance`] evaluated at `supplied`, or
    /// [`FORWARD_CONSISTENCY_TOLERANCE`] when `supplied` is not finite — the relative term is
    /// deliberately not applied to an infinite forward, which would make the tolerance infinite
    /// and the comparison `inf > inf` false.  See
    /// [`max_forward_inconsistency`].
    pub tolerance: f64,
}

/// Walk [`CURVE_PROBE_TIMES`] and report the worst `forward` vs `d/dt zero_yield` disagreement.
///
/// See [`ForwardInconsistency`]; [`HullWhite::new`](crate::HullWhite::new) calls this and
/// rejects anything over tolerance.
pub fn max_forward_inconsistency(curve: &dyn YieldCurve) -> Option<ForwardInconsistency> {
    let yield_fn = |s: f64| curve.zero_yield(s);
    let mut worst: Option<ForwardInconsistency> = None;
    for &t in CURVE_PROBE_TIMES.iter() {
        //A non-finite derivative means the curve is not defined (or not smooth) at this
        //probe, which says nothing about the pair; skip rather than call it inconsistent.
        let derived = derivative(&yield_fn, t);
        if !derived.is_finite() {
            continue;
        }
        let supplied = curve.forward(t);
        //A supplied forward that is not finite while the derivative is, *is* an inconsistency:
        //the curve claims a slope that does not exist.  Such a value never lands inside
        //tolerance, and the tolerance deliberately does not follow the relative term here:
        //`1e-4 * |inf|` is `inf`, and `inf > inf` is false, so scaling the tolerance by a
        //non-finite forward would let `forward(t) = inf` pass the check that exists to fail it.
        //`phi(t)` would then carry that infinity into `mu_r` and the Jamshidian bracket, and
        //hand back a finite, wrong, unexplained price.
        let (error, tolerance) = if supplied.is_finite() {
            (
                (supplied - derived).abs(),
                forward_consistency_tolerance(supplied),
            )
        } else {
            (f64::INFINITY, FORWARD_CONSISTENCY_TOLERANCE)
        };
        let bad = ForwardInconsistency {
            t,
            supplied,
            derived,
            error,
            tolerance,
        };
        if worst.is_none_or(|w| error > w.error) {
            worst = Some(bad);
        }
    }
    worst
}

/// Turn a [`ForwardInconsistency`] into the error a constructor returns.
pub(crate) fn inconsistency_error(bad: &ForwardInconsistency) -> HullWhiteError {
    HullWhiteError::InvalidInput(format!(
        "forward curve is inconsistent with the yield curve: d/dt zero_yield({t}) = {derived} \
         but forward({t}) = {supplied} (|diff| = {error}, tolerance = {tolerance}).  A \
         Hull-White model calibrated to a cumulative yield must be given the instantaneous \
         forward f(0,t) = d/dt zero_yield(t); see hull_white::YieldCurve.",
        t = bad.t,
        derived = bad.derived,
        supplied = bad.supplied,
        error = bad.error,
        tolerance = bad.tolerance,
    ))
}

/// The error returned when a curve cannot be probed at all.
pub(crate) fn unprobeable_error() -> HullWhiteError {
    HullWhiteError::InvalidInput(format!(
        "yield curve returned no finite value at any probe time {CURVE_PROBE_TIMES:?}; \
         a model cannot be calibrated to a curve that is not defined"
    ))
}

/// Validate a curve for use as a model's initial term structure.
///
/// Checks [`max_forward_inconsistency`] against the tolerance and reports both failure modes
/// (nothing probeable, and inconsistent) as
/// [`HullWhiteError::InvalidInput`].  A `forward`
/// that is not finite over a finite `zero_yield` is always in the second group: it is a slope
/// the curve does not have.
pub fn validate_curve(curve: &dyn YieldCurve) -> Result<(), HullWhiteError> {
    match max_forward_inconsistency(curve) {
        None => Err(unprobeable_error()),
        Some(bad) if bad.error > bad.tolerance => Err(inconsistency_error(&bad)),
        Some(_) => Ok(()),
    }
}

/// A [`YieldCurve`] from one closure: the cumulative yield.
///
/// The forward is derived by `derivative()`, so there is no pair to be inconsistent — which is
/// the point of this shape.  Supply the derivative instead (implement the trait, or use
/// [`from_yield_and_forward`]) when the finite difference's ~`1e-11` is not enough or when the
/// curve is evaluated in a hot loop.
///
/// ```
/// use hull_white::HullWhite;
///
/// // y(t) = 0.05 t + 0.01 t^2 -- one closure, forward derived as 0.05 + 0.02 t.
/// let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
/// let model = HullWhite::new(0.2, 0.03, &curve).unwrap();
/// assert!((model.short_rate_now().unwrap() - 0.05).abs() < 1e-9);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct FromYield<F>
where
    F: Fn(f64) -> f64 + Sync,
{
    yield_fn: F,
}

impl<F> YieldCurve for FromYield<F>
where
    F: Fn(f64) -> f64 + Sync,
{
    fn zero_yield(&self, t: f64) -> f64 {
        (self.yield_fn)(t)
    }
    // `forward` and `discount` come from the trait defaults: the forward is the finite
    // difference of the yield, the price is its exponential.
}

/// Wrap a cumulative-yield closure as a [`YieldCurve`], deriving the forward.
///
/// See [`FromYield`].
pub fn from_yield<F>(yield_fn: F) -> FromYield<F>
where
    F: Fn(f64) -> f64 + Sync,
{
    FromYield { yield_fn }
}

/// A [`YieldCurve`] from the pair of closures this crate used to demand, carried as one object.
///
/// Both the cumulative yield and the instantaneous forward are supplied, so nothing is
/// finite-differenced; the construction check then *verifies* the pair instead of trusting it.
/// This is the direct migration from `HullWhite::init(a, sigma, &yield_curve, &forward_curve)` —
/// same closures, one object, plus the diagnostic the old call never gave.
///
/// ```
/// use hull_white::HullWhite;
///
/// // Consistent: the forward is the derivative of the cumulative yield.
/// let curve = hull_white::from_yield_and_forward(
///     |t: f64| 0.05 * t + 0.01 * t * t,
///     |t: f64| 0.05 + 0.02 * t,
/// );
/// assert!(HullWhite::new(0.2, 0.03, &curve).is_ok());
///
/// // Not consistent: 0.05 + 0.01 t is not the derivative of 0.05 t + 0.01 t^2.
/// let broken = hull_white::from_yield_and_forward(
///     |t: f64| 0.05 * t + 0.01 * t * t,
///     |t: f64| 0.05 + 0.01 * t,
/// );
/// let err = HullWhite::new(0.2, 0.03, &broken).unwrap_err();
/// assert!(format!("{err}").contains("inconsistent"));
/// ```
#[derive(Debug, Clone, Copy)]
pub struct YieldAndForward<F, G>
where
    F: Fn(f64) -> f64 + Sync,
    G: Fn(f64) -> f64 + Sync,
{
    yield_fn: F,
    forward_fn: G,
}

impl<F, G> YieldCurve for YieldAndForward<F, G>
where
    F: Fn(f64) -> f64 + Sync,
    G: Fn(f64) -> f64 + Sync,
{
    fn zero_yield(&self, t: f64) -> f64 {
        (self.yield_fn)(t)
    }
    fn forward(&self, t: f64) -> f64 {
        (self.forward_fn)(t)
    }
}

/// Wrap a `(cumulative yield, instantaneous forward)` closure pair as a single [`YieldCurve`].
///
/// See [`YieldAndForward`].
pub fn from_yield_and_forward<F, G>(yield_fn: F, forward_fn: G) -> YieldAndForward<F, G>
where
    F: Fn(f64) -> f64 + Sync,
    G: Fn(f64) -> f64 + Sync,
{
    YieldAndForward {
        yield_fn,
        forward_fn,
    }
}

/// How the finite-difference step in [`derivative`] is chosen: `1e-4` years per unit of `|t|`,
/// at least `1e-4` years.
pub(crate) const DERIVATIVE_STEP: f64 = 1e-4;

/// Five-point first derivative of `f` at `t`.
///
/// Central when the stencil stays on the domain (`t >= 2h`), forward otherwise: every time in
/// this crate is measured from "now" (0) and there is nothing to the left of it.
///
/// At `h = 1e-4 * max(1, |t|)` the two error terms balance: rounding is `~eps * |f| / h`
/// (`~1e-12` for a ten-year cumulative yield) and truncation is `h^4 f^(5) / 30`, which is
/// below that for anything a rate curve looks like.  Exact for polynomials up to degree four, so
/// the example curves in the docs come back to their analytic derivatives to the last bit.
pub(crate) fn derivative(f: &dyn Fn(f64) -> f64, t: f64) -> f64 {
    let h = DERIVATIVE_STEP * t.abs().max(1.0);
    if t - 2.0 * h >= 0.0 {
        (f(t - 2.0 * h) - 8.0 * f(t - h) + 8.0 * f(t + h) - f(t + 2.0 * h)) / (12.0 * h)
    } else {
        (-25.0 * f(t) + 48.0 * f(t + h) - 36.0 * f(t + 2.0 * h) + 16.0 * f(t + 3.0 * h)
            - 3.0 * f(t + 4.0 * h))
            / (12.0 * h)
    }
}

//tdiff=T-t
pub(crate) fn a_t(a: f64, t_diff: f64) -> f64 {
    (1.0 - (-a * t_diff).exp()) / a
}

//t is first future time
//t_m is second future time
pub(crate) fn at_t(a: f64, t: f64, t_m: f64) -> f64 {
    a_t(a, t_m - t)
}

/// The affine constant `C(t, t_m)` of the zero-coupon price
/// `P(t, t_m) = exp(C(t, t_m) - B(t, t_m) * r_t)`, fitted to the initial curve.
///
/// Takes the curve as one object because both halves of it are needed here and they are only
/// meaningful together: `y(t) - y(t_m)` moves along the cumulative yield and `f(t)` is the
/// slope that yield is built from.  A model whose `forward` is not the derivative of its
/// `zero_yield` makes this the difference of two unrelated curves — which is exactly why
/// [`validate_curve`] runs before anything calls this.
pub(crate) fn ct_t(a: f64, sigma: f64, t: f64, t_m: f64, curve: &dyn YieldCurve) -> f64 {
    let sqr = (-a * t_m).exp() - (-a * t).exp();
    curve.zero_yield(t) - curve.zero_yield(t_m) + curve.forward(t) * at_t(a, t, t_m)
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

/// What a [`HullWhite`](crate::HullWhite) holds onto: the curve, borrowed or built.
///
/// [`HullWhite::new`](crate::HullWhite::new) borrows the caller's curve object;
/// [`HullWhite::init`](crate::HullWhite::init) (deprecated) has no curve object to borrow, so
/// it builds the [`YieldAndForward`] wrapper for the two closures and owns it.  Both read back
/// through [`CurveRef::get`] as `&dyn YieldCurve`, so nothing downstream cares which it is.
pub(crate) enum CurveRef<'a> {
    /// Borrowed from the caller, the zero-cost case.
    Borrowed(&'a dyn YieldCurve),
    /// Built by a constructor from parts that were not a curve to begin with.
    Owned(Box<dyn YieldCurve + 'a>),
}

impl CurveRef<'_> {
    pub(crate) fn get(&self) -> &dyn YieldCurve {
        match self {
            CurveRef::Borrowed(curve) => *curve,
            CurveRef::Owned(curve) => curve.as_ref(),
        }
    }
}

#[cfg(test)]
mod tests;
