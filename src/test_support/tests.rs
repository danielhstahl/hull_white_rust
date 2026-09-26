//! Tests for the fixture itself.
//!
//! These exist because the shared curves replaced ~50 hand-copied ones: if the replacement were
//! not the same maths, every test that prices off it would drift quietly.  So the pair is pinned
//! three ways — to its own algebra (`y` is the integral of `F`, `F(0)` is the spot), to golden
//! values computed independently of this file, and to the model's own bond pricer.

use approx::*;

use crate::HullWhite;
use crate::test_support::*;

/// `y(0) = 0` (nothing accrued), `F(0) = curr_rate` (the volatility term of `phi` vanishes at
/// zero), and `y` integrates `F`: `dy/dt == F(t)`.
///
/// The last one is the property the whole fixture rests on.  `exp(-y(T))` is the discount factor
/// only because `y` is the accumulated forward curve, so if the two closures were ever edited
/// apart — the exact failure mode of copy-pasting them — this fails on every scenario.
#[test]
fn every_scenario_is_internally_consistent() {
    for s in ALL_SCENARIOS {
        let (yield_curve, forward_curve) = s.curves();
        assert_abs_diff_eq!(yield_curve(0.0), 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(forward_curve(0.0), s.curr_rate, epsilon = 1e-15);
        // Central difference: truncation error is O(h^2) and rounding O(eps/h), which at
        // h = 1e-5 leaves ~1e-11 on the sloppiest of these calibrations.
        let h = 1e-5;
        for t in [0.25, 0.5, 1.0, 2.0, 5.0, 10.0] {
            let slope = (yield_curve(t + h) - yield_curve(t - h)) / (2.0 * h);
            assert_abs_diff_eq!(slope, forward_curve(t), epsilon = 1e-9);
        }
    }
}

/// Golden values, computed from the same closed form outside this crate (double-precision
/// arithmetic on the formulas in [`hw_curves`]' doc), and so a description of what the
/// copy-pasted fixture used to produce rather than of what this code happens to compute.
#[allow(clippy::excessive_precision)]
#[test]
fn fixture_matches_golden_curve_values() {
    // (scenario, t, yield_curve(t), forward_curve(t))
    let goldens: [(Scenario, f64, f64, f64); 15] = [
        (
            BASELINE,
            1.0,
            2.266_765_453_940_028_4e-2,
            2.503_435_737_585_322_0e-2,
        ),
        (
            BASELINE,
            5.0,
            1.450_874_418_905_402_9e-1,
            3.419_622_624_576_250_4e-2,
        ),
        (
            BASELINE,
            10.0,
            3.248_129_544_789_879_5e-1,
            3.699_780_393_166_286_1e-2,
        ),
        (
            VASICEK_REFERENCE,
            1.0,
            1.059_315_076_001_475_2e-2,
            1.103_497_483_876_290_2e-2,
        ),
        (
            VASICEK_REFERENCE,
            5.0,
            5.167_001_921_146_794_5e-2,
            7.828_739_665_289_595_5e-3,
        ),
        (
            VASICEK_REFERENCE,
            10.0,
            5.908_064_000_521_569_3e-2,
            -6.063_181_705_690_583_2e-3,
        ),
        (
            FLAT_5PCT,
            1.0,
            4.998_394_400_662_070_6e-2,
            4.995_242_861_930_937_5e-2,
        ),
        (
            FLAT_5PCT,
            5.0,
            2.482_655_054_854_028_0e-1,
            4.902_141_812_860_352_8e-2,
        ),
        (
            FLAT_5PCT,
            10.0,
            4.883_513_604_641_817_1e-1,
            4.690_363_756_507_649_1e-2,
        ),
        (
            STEEP_CURVE,
            1.0,
            2.361_669_218_907_576_4e-2,
            2.688_111_130_323_437_8e-2,
        ),
        (
            STEEP_CURVE,
            5.0,
            1.641_207_559_435_309_3e-1,
            4.078_958_784_308_787_3e-2,
        ),
        (
            STEEP_CURVE,
            10.0,
            3.842_319_646_269_493_1e-1,
            4.617_558_160_586_102_1e-2,
        ),
        (
            LOW_VOL,
            1.0,
            4.999_938_108_093_414_4e-2,
            4.999_818_881_659_878_8e-2,
        ),
        (
            LOW_VOL,
            5.0,
            2.499_417_568_023_209_2e-1,
            4.996_903_637_565_076_9e-2,
        ),
        (
            LOW_VOL,
            10.0,
            4.996_638_175_185_508_2e-1,
            4.992_008_471_982_125_5e-2,
        ),
    ];
    for (s, t, y, f) in goldens {
        let (yield_curve, forward_curve) = s.curves();
        assert_abs_diff_eq!(yield_curve(t), y, epsilon = 1e-12);
        assert_abs_diff_eq!(forward_curve(t), f, epsilon = 1e-12);
    }
}

/// The fixture is *consistent* with the model rather than merely smooth: `exp(-y(T))` has to be
/// the model's own `T`-maturity bond price, and `F(t)` the model's own `phi(t)` at the front.
#[test]
fn fixture_prices_are_what_the_model_says_they_are() {
    for s in ALL_SCENARIOS {
        let (yield_curve, forward_curve) = s.curves();
        let model = HullWhite::init(s.a, s.sigma, &yield_curve, &forward_curve).unwrap();
        for t in [0.25, 1.0, 3.0, 7.5] {
            assert_abs_diff_eq!(
                model.bond_price_now(t).unwrap(),
                (-yield_curve(t)).exp(),
                epsilon = 1e-15
            );
        }
        // `short_rate_now` is the instantaneous forward at the front of the curve, so it is
        // `phi(0)`, which is `forward_curve(0.0)`.
        assert_abs_diff_eq!(
            model.short_rate_now().unwrap(),
            forward_curve(0.0),
            epsilon = 1e-15
        );
    }
}

/// The scenarios have to stay *different*.
///
/// If two scenarios collapsed onto the same calibration, the suite would silently lose a
/// dimension of coverage while still counting every test.
#[test]
fn scenarios_are_distinct() {
    let names: Vec<&str> = ALL_SCENARIOS.iter().map(|s| s.name).collect();
    let mut sorted = names.clone();
    sorted.sort_unstable();
    sorted.dedup();
    assert_eq!(sorted.len(), ALL_SCENARIOS.len(), "duplicate scenario name");
    for (i, s) in ALL_SCENARIOS.iter().enumerate() {
        for o in ALL_SCENARIOS.iter().skip(i + 1) {
            assert_ne!(
                (s.curr_rate, s.a, s.b, s.sigma),
                (o.curr_rate, o.a, o.b, o.sigma),
                "{} and {} are the same calibration",
                s.name,
                o.name
            );
        }
    }
}

/// What the names claim about the shape of the curve is actually true.
///
/// `flat_5pct` really is nearly flat and `steep_curve` really travels; if a scenario were edited
/// into the other one, the "coverage" the list advertises would be a lie, and the tests written
/// against the flat fixture would stop being the cheap ones they claim to be.
#[test]
fn the_scenario_names_describe_the_curves() {
    let forward_at = |s: Scenario, t: f64| s.curves().1(t);

    // Flat: 5% spot, and the forward barely moves over ten years.
    assert!(
        (forward_at(FLAT_5PCT, 10.0) - FLAT_5PCT.curr_rate).abs() < 5e-3,
        "flat_5pct drifts too much: {}",
        forward_at(FLAT_5PCT, 10.0)
    );
    // Steep: the same 10y drift is an order of magnitude larger, and it moves monotonically up.
    let steep_drift = forward_at(STEEP_CURVE, 10.0) - STEEP_CURVE.curr_rate;
    assert!(
        steep_drift > 2.0e-2,
        "steep_curve barely moves: {steep_drift}"
    );
    assert!(
        steep_drift > 5.0 * (forward_at(FLAT_5PCT, 10.0) - FLAT_5PCT.curr_rate).abs(),
        "steep_curve is not distinguishable from flat_5pct"
    );
    assert!(forward_at(STEEP_CURVE, 1.0) > STEEP_CURVE.curr_rate);
    assert!(forward_at(STEEP_CURVE, 10.0) > forward_at(STEEP_CURVE, 1.0));

    // Low vol: `sigma` is the smallest of any scenario's, and comfortably below the next one,
    // so "low vol" is not just a name on a curve that is a shade quieter than its neighbours.
    let lowest = ALL_SCENARIOS
        .iter()
        .map(|s| s.sigma)
        .fold(1.0_f64, f64::min);
    assert_eq!(lowest, LOW_VOL.sigma);
    assert!(
        ALL_SCENARIOS
            .iter()
            .filter(|s| s.name != LOW_VOL.name)
            .all(|s| s.sigma > 2.0 * LOW_VOL.sigma),
        "low_vol is not actually the low-volatility scenario"
    );
}

/// `LOW_VOL` earns its place instead of merely sounding quiet.
///
/// With `sigma = 0.002` the convexity that separates a Eurodollar future from the forward it is
/// quoted against is a fiftieth of a basis point, and an at-the-money caplet is under 2% of the
/// period's rate notional.  On `BASELINE` both terms are an order of magnitude larger, so the
/// comparison below measures volatility rather than a fixture that is uniformly tiny.
#[test]
fn low_vol_collapses_the_convexity_terms() {
    // (edf-vs-forward gap, at-the-money caplet as a share of `delta * forward`)
    let measure = |s: Scenario| {
        let (yield_curve, forward_curve) = s.curves();
        let model = HullWhite::init(s.a, s.sigma, &yield_curve, &forward_curve).unwrap();
        let option_maturity = 1.5;
        let forward = model
            .forward_libor_rate_now(option_maturity, s.delta)
            .unwrap();
        let edf = model
            .euro_dollar_future_now(option_maturity, s.delta)
            .unwrap();
        let caplet = model.caplet_now(option_maturity, s.delta, forward).unwrap();
        ((forward - edf).abs(), caplet / (forward * s.delta))
    };
    let low = measure(LOW_VOL);
    let base = measure(BASELINE);
    assert!(
        low.0 < 1e-5,
        "{}: edf convexity {} is not a rounding error",
        LOW_VOL.name,
        low.0
    );
    assert!(
        base.0 > 50.0 * low.0,
        "baseline convexity {} is not clearly bigger than {}'s {}",
        base.0,
        LOW_VOL.name,
        low.0
    );
    assert!(
        low.1 < 0.02,
        "{}: atm caplet is {} of the period notional, not a rounding error",
        LOW_VOL.name,
        low.1
    );
    assert!(
        base.1 > 10.0 * low.1,
        "baseline atm caplet share {} is not clearly bigger than {}'s {}",
        base.1,
        LOW_VOL.name,
        low.1
    );
}

/// `Scenario::curves` is `hw_curves` with the scenario's numbers, not a second implementation.
#[test]
fn scenario_curves_are_just_hw_curves() {
    for s in ALL_SCENARIOS {
        let (a_y, a_f) = s.curves();
        let (b_y, b_f) = hw_curves(s.curr_rate, s.a, s.b, s.sigma);
        for t in [0.0, 0.5, 2.5, 20.0] {
            assert_eq!(a_y(t), b_y(t));
            assert_eq!(a_f(t), b_f(t));
        }
    }
}
