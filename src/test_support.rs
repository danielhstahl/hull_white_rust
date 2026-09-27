//! Shared, self-consistent market-data fixtures for this crate's tests and benches.
//!
//! Every model here is calibrated to one curve: a [`YieldCurve`] whose `zero_yield` is the
//! cumulative yield (`bond_price_now(T) = exp(-zero_yield(T))`) and whose `forward` is the
//! instantaneous forward that yield is the integral of.  Tests and benches need a curve that is
//! *internal* to the model it calibrates, so a fixture is built from a Vasicek short rate
//! reverting to `b` at speed `a` with volatility `sigma`, started at `curr_rate`.  For that
//! curve the model's closed forms are exact, so a test can put `x_now` against `x_t`, or a
//! Monte-Carlo average against an analytic price, without a second pricer to disagree with — and
//! a time-coordinate error cannot hide, because the yield and the forward carry the same
//! convexity the model carries.
//!
//! That last property is also what makes the fixture clear the construction-time consistency
//! check with room to spare.  Measured worst disagreement over the seven scenarios, `forward(t)`
//! against the finite-difference derivative of `zero_yield(t)`:
//!
//! ```text
//! scenario            worst |diff|   at t     tolerance
//! baseline            5.0e-14      10.0     3.7e-6
//! vasicek_reference   9.5e-13       0.25    1.0e-6
//! flat_5pct           1.3e-13       0.25    5.0e-6
//! steep_curve         1.2e-13       1.0     2.7e-6
//! high_vol            2.0e-12       0.5     5.1e-6
//! quick_reversion     4.8e-14       5.0     4.3e-6
//! low_vol             6.1e-14       5.0     5.0e-6
//! ```
//!
//! Six orders of headroom on the loosest scenario, all of it finite-difference rounding on the
//! check's side rather than sloppiness in the fixture.  That is not true of the ad-hoc curve
//! pairs that used to be pasted around this suite: the crate's own long-standing example pair
//! (`0.05 * t` cumulative yield against `t.ln()` forward) is off by order unity and would not
//! build a model at all today.
//!
//! This maths used to be copy-pasted into every test and bench.  It is written once, in
//! [`HwCurve`], and nothing else repeats it:
//!
//! * [`hw_curve`] is the maths, for the rare test that needs an off-scenario calibration.
//! * [`Scenario`] names one calibration of it — `curr_rate`, `a`, `b`, `sigma`, plus `delta`,
//!   the period length the instruments built on it are quoted on.  A test says which part of the
//!   parameter space it covers by *name*, instead of by four unexplained literals next to a
//!   re-typed curve.
//! * [`ALL_SCENARIOS`] is the list of them, and `test_support::tests` walks it, so the coverage
//!   is a greppable set rather than an accident of how many times the block was pasted.
//!
//! The scenarios are spread over the shape space on purpose.  `curr_rate` relative to `b` decides
//! whether `phi(t)` moves at all (`curr_rate == b` freezes it), `a` decides how fast it moves,
//! `sigma` decides how much the stochastic part is worth.  A calibration that is flat *and*
//! low-vol cannot expose a time error, and one that is steep cannot separate a drift error from a
//! vol error, so both ends are kept rather than averaged into a single "default" curve.
//!
//! ```text
//! scenario            curr   a     b     sigma  what it exercises
//! baseline            0.02   0.30  0.04  0.020  fast reversion, 2s spot under a 4s mean
//! vasicek_reference   0.01   0.05  0.04  0.030  slow reversion; matches the quantcalc.net
//!                                              BondOption_Vasicek worked example
//! flat_5pct           0.05   0.05  0.05  0.010  curr == b: phi(t) is (nearly) constant,
//!                                              the swaption/bench calibration
//! steep_curve         0.02   0.20  0.06  0.030  curr far below a high mean, high vol:
//!                                              phi(t) travels, so time errors show
//! high_vol            0.05   0.10  0.08  0.080  the far-vol corner: the rate wanders, the
//!                                              Jamshidian bracket has to travel with it
//! quick_reversion     0.03   0.40  0.045 0.005  pulled to its mean in a couple of years,
//!                                              barely stochastic
//! low_vol             0.05   0.10  0.05  0.002  near-deterministic: convexity and option
//!                                              value are small by construction
//! ```
//!
//! Compiled into the crate for tests (`cfg(test)`) and for anything that turns on the hidden
//! `test-support` feature — which is how `benches/` reaches it, since a bench is a separate crate
//! and cannot see `cfg(test)` items.  It is `#[doc(hidden)]`: it is not part of the public
//! surface of this crate.

use crate::curves::YieldCurve;

/// The Hull-White-consistent initial curve for one calibration, as a single [`YieldCurve`].
///
/// This is the whole fixture maths, and the only place it appears:
///
/// ```text
/// a(t)   = (1 - e^{-a t}) / a                                  (bond duration)
/// c(t)   = (b - sigma^2 / 2a^2) (a(t) - t) - (sigma a(t))^2 / 4a
/// y(t)   = a(t) * curr_rate - c(t)                             cumulative yield
/// F(t)   = b + e^{-a t} (curr_rate - b)
///          - (sigma^2 / 2a^2) (1 - e^{-a t})^2                 instantaneous forward
/// ```
///
/// `F(t)` is the model's own `phi(t)` minus its volatility term, and `y` is its integral, which
/// is what makes `exp(-y(T))` the exact `T`-maturity bond price under the model these curves
/// calibrate.  Both accessors take `t` as *time from now*, the same convention as every pricer in
/// the crate, and both are closed form — the construction-time consistency check compares
/// `forward` against a finite difference of `zero_yield` and finds them agreeing to about
/// `1e-12` (worst measured case, `high_vol`) where the tolerance is `1e-6`, so nothing here
/// pays for a finite difference in the price either.
///
/// The shape is one object on purpose: a fixture that hands back a `(yield, forward)` *tuple* is
/// a fixture that can be edited apart, which is the exact failure mode the 52 copies of this
/// maths used to invite.  [`hw_curve`] builds it; [`Scenario::curve`] is the same thing with a
/// name attached.
///
/// Reach for a [`Scenario`] rather than [`hw_curve`] unless the calibration is genuinely a
/// one-off; a bare `hw_curve(0.017, 0.13, ...)` call tells a reader nothing about which corner
/// of the parameter space the test covers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HwCurve {
    /// Short rate at `t = 0`, which is also the front of the forward curve (`forward(0.0)`).
    pub curr_rate: f64,
    /// Mean-reversion speed of the underlying Vasicek/Hull-White short rate.
    pub a: f64,
    /// Long-run mean the short rate reverts to.
    pub b: f64,
    /// Short-rate volatility.
    pub sigma: f64,
}

impl HwCurve {
    fn duration(&self, t: f64) -> f64 {
        (1.0 - (-self.a * t).exp()) / self.a
    }
}

impl YieldCurve for HwCurve {
    /// `y(t) = a(t) * curr_rate - c(t)`, the integral of [`HwCurve::forward`].
    fn zero_yield(&self, t: f64) -> f64 {
        let at = self.duration(t);
        let ct = (self.b - self.sigma.powi(2) / (2.0 * self.a.powi(2))) * (at - t)
            - (self.sigma * at).powi(2) / (4.0 * self.a);
        at * self.curr_rate - ct
    }

    /// `F(t) = b + e^{-a t} (curr_rate - b) - (sigma^2 / 2a^2) (1 - e^{-a t})^2`, the
    /// derivative of [`HwCurve::zero_yield`], and `phi(t)` with its volatility term taken back off.
    fn forward(&self, t: f64) -> f64 {
        self.b + (-self.a * t).exp() * (self.curr_rate - self.b)
            - (self.sigma.powi(2) / (2.0 * self.a.powi(2))) * (1.0 - (-self.a * t).exp()).powi(2)
    }
}

/// The [`HwCurve`] for `curr_rate` / `a` / `b` / `sigma`.
///
/// Prefer [`Scenario::curve`] so the calibration has a name.
pub fn hw_curve(curr_rate: f64, a: f64, b: f64, sigma: f64) -> HwCurve {
    HwCurve {
        curr_rate,
        a,
        b,
        sigma,
    }
}

/// One named calibration of [`hw_curve`].
///
/// The four model numbers plus `delta`, the accrual period the instruments priced off the
/// scenario are quoted on.  Cheap `Copy` value: take what a test needs
/// (`let curr_rate = FLAT_5PCT.curr_rate;`) or the curves straight off it
/// (`let curve = FLAT_5PCT.curve();`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Scenario {
    /// The scenario's name, for assertion messages and for grepping this file.
    pub name: &'static str,
    /// Short rate at `t = 0`, which is also the front of the forward curve (`forward(0.0)`).
    pub curr_rate: f64,
    /// Mean-reversion speed of the underlying Vasicek/Hull-White short rate.
    pub a: f64,
    /// Long-run mean the short rate reverts to.
    pub b: f64,
    /// Short-rate volatility.
    pub sigma: f64,
    /// Period length of the instruments priced off this scenario (quarterly by convention here).
    pub delta: f64,
}

impl Scenario {
    /// This scenario's curve, as the single [`YieldCurve`] a [`HullWhite`](crate::HullWhite) is
    /// built from.
    pub fn curve(self) -> HwCurve {
        hw_curve(self.curr_rate, self.a, self.b, self.sigma)
    }
}

/// The crate's default legacy calibration: a 2s short rate under a 4s mean, reverting fast
/// (`a = 0.3`), moderate vol.
///
/// `curr_rate < b` with a fast `a` means the curve bends early and then flattens onto the long-run
/// mean, which is what most of the bond, caplet and swap tests were written against.
pub const BASELINE: Scenario = Scenario {
    name: "baseline",
    curr_rate: 0.02,
    a: 0.3,
    b: 0.04,
    sigma: 0.02,
    delta: 0.25,
};

/// Slow reversion (`a = 0.05`) from a 1s short rate towards a 4s mean, high vol.
///
/// This is the calibration the closed-form bond-option values in
/// `src/options/tests.rs` were read off quantcalc.net's `BondOption_Vasicek` page with; changing
/// any of the four numbers invalidates that external reference value.
pub const VASICEK_REFERENCE: Scenario = Scenario {
    name: "vasicek_reference",
    curr_rate: 0.01,
    a: 0.05,
    b: 0.04,
    sigma: 0.03,
    delta: 0.25,
};

/// Flat at 5%: `curr_rate == b`, low vol, slow reversion.
///
/// The term structure barely moves — `phi(t)` is constant to first order — so this is the fixture
/// for instruments whose *shape* is the interesting part (swaptions, the bench calibration) and
/// the wrong fixture for anything checking that time is being carried through correctly.  Use
/// [`STEEP_CURVE`] for that.
pub const FLAT_5PCT: Scenario = Scenario {
    name: "flat_5pct",
    curr_rate: 0.05,
    a: 0.05,
    b: 0.05,
    sigma: 0.01,
    delta: 0.25,
};

/// Deliberately steep: a 2s short rate far below a 6s mean, `a = 0.2`, `sigma = 0.03`.
///
/// `phi(t)` travels a long way over an option's life here, so a price that used the wrong time
/// coordinate cannot come out right by accident.  This is the fixture behind the invalid-input and
/// Jamshidian tests; the old `flat_5pct` made `phi(t)` constant and hid exactly that class of
/// bug, which is why this one exists.
pub const STEEP_CURVE: Scenario = Scenario {
    name: "steep_curve",
    curr_rate: 0.02,
    a: 0.2,
    b: 0.06,
    sigma: 0.03,
    delta: 0.25,
};

/// Near-deterministic: flat at 5% with `sigma = 0.002`.
///
/// Convexity is `O(sigma^2)`, so on this scenario the gap between a Eurodollar future and the
/// forward Libor it quotes against is under a tenth of a basis point, and an at-the-money caplet
/// is a rounding error on the period's notional.  That is the point: this scenario lets a test
/// assert the *limit* the others approximate, so the vol-sensitive terms can be shown to be
/// small rather than merely asserted absent.
pub const LOW_VOL: Scenario = Scenario {
    name: "low_vol",
    curr_rate: 0.05,
    a: 0.1,
    b: 0.05,
    sigma: 0.002,
    delta: 0.25,
};

/// The far-volatility corner: `sigma = 0.08` against a 300bp gap between spot and the long-run
/// mean.
///
/// The short rate genuinely wanders here, so the critical rate a coupon-bond option solves for can
/// sit a long way from where it starts and the tree's tails matter.  This is the scenario that
/// breaks a solver with no bracket, not one with a bad price.
pub const HIGH_VOL: Scenario = Scenario {
    name: "high_vol",
    curr_rate: 0.05,
    a: 0.1,
    b: 0.08,
    sigma: 0.08,
    delta: 0.25,
};

/// A fast, quiet calibration: `a = 0.4` pulls the rate onto its 4.5% mean in about two years and
/// `sigma = 0.005` barely moves it after that.
///
/// The mirror image of [`HIGH_VOL`]: the critical rate sits close to the mean, the objective is
/// nearly linear in it, and anything that is only right "by a little" looks exactly right here.
pub const QUICK_REVERSION: Scenario = Scenario {
    name: "quick_reversion",
    curr_rate: 0.03,
    a: 0.4,
    b: 0.045,
    sigma: 0.005,
    delta: 0.25,
};

/// Every named scenario, so coverage can be walked and asserted over rather than listed by hand.
///
/// A test that takes `[Scenario]` from this list silently inherits any scenario added to it.
pub const ALL_SCENARIOS: [Scenario; 7] = [
    BASELINE,
    VASICEK_REFERENCE,
    FLAT_5PCT,
    STEEP_CURVE,
    HIGH_VOL,
    QUICK_REVERSION,
    LOW_VOL,
];

/// Set `model` up on a named scenario.
///
/// Every path inside is `$crate`-qualified so a caller only needs `use
/// $crate::test_support::hw_setup;` — the fixture stays in one place and cannot be shadowed at
/// the call site.
///
/// ```ignore
/// hw_setup!(hull_white);                      // steep_curve, the default
/// hw_setup!(hull_white, FLAT_5PCT);          // an explicit scenario
/// ```
#[allow(unused_macros)]
macro_rules! hw_setup {
    ($model:ident) => {
        $crate::test_support::hw_setup!($model, $crate::test_support::STEEP_CURVE);
    };
    ($model:ident, $scenario:expr) => {
        let __scenario = $scenario;
        let __curve = __scenario.curve();
        let $model = $crate::HullWhite::new(__scenario.a, __scenario.sigma, &__curve).unwrap();
    };
}

#[cfg(test)]
pub(crate) use hw_setup;

#[cfg(test)]
mod tests;
