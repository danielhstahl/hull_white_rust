//! Harness self-tests: the estimator statistics themselves, the grid convention, and the bias
//! arithmetic.  These check the harness, not any instrument — the instrument-level checks that
//! consume it are in `crate::rates::tests`.

use super::*;

/// Welford in the harness has to agree with a two-pass computation over the very same payoffs.
///
/// The payoff closure records every value it is handed, so the reference mean and variance are
/// computed from exactly the numbers the accumulator saw, and a one-pass/online bug shows up as a
/// mismatch rather than as a slightly-off price somewhere downstream.
#[test]
fn pair_statistics_match_a_two_pass_computation() {
    //Every payoff the accumulator sees is recorded, in order, so the reference statistics below
    //are computed over exactly the same numbers (round index included, to make each pair's value
    //depend on something other than the increments alone).
    let mut recorder: Vec<f64> = Vec::new();
    let mut mc = AntitheticMonteCarlo::new(rng_for(1), 3);
    for round in 0..500usize {
        mc.push_pair(|normals| {
            let value =
                normals.iter().sum::<f64>() + 0.5 * normals[0] * normals[0] + 0.01 * round as f64;
            recorder.push(value);
            value
        });
    }
    let est = mc.finish();
    assert_eq!(recorder.len(), 1000, "two paths per pair");

    let pairs: Vec<f64> = recorder
        .chunks(2)
        .map(|pair| 0.5 * (pair[0] + pair[1]))
        .collect();
    let n = pairs.len() as f64;
    let mean = pairs.iter().sum::<f64>() / n;
    let var = pairs.iter().map(|p| (p - mean).powi(2)).sum::<f64>() / (n - 1.0);
    let two_pass_se = (var / n).sqrt();

    assert!(
        (est.mean - mean).abs() <= 1e-9 * mean.abs().max(1.0),
        "online mean {} vs two-pass {mean}",
        est.mean
    );
    assert!(
        (est.standard_error - two_pass_se).abs() <= 1e-9 * two_pass_se,
        "online SE {} vs two-pass {two_pass_se}",
        est.standard_error
    );
    assert_eq!(est.pairs, 500);
    assert_eq!(est.paths, 1000);

    //And the unpaired half of the report: same numbers, treated as independent draws.
    let all_mean = recorder.iter().sum::<f64>() / 1000.0;
    let all_var = recorder.iter().map(|v| (v - all_mean).powi(2)).sum::<f64>() / 999.0;
    let expected_naive = (all_var / 1000.0).sqrt();
    assert!(
        (est.naive_standard_error - expected_naive).abs() <= 1e-9 * expected_naive,
        "online unpaired SE {} vs two-pass {expected_naive}",
        est.naive_standard_error
    );
    //The payoff has an odd part, so pairing has to have bought something.
    assert!(est.variance_reduction > 1.0, "{}", est.variance_reduction);
}

/// An odd payoff is cancelled exactly by the mirror path: the antithetic estimator has zero
/// variance while the unpaired one has the full `sqrt(k / n_paths)`.
#[test]
fn an_odd_payoff_is_cancelled_by_the_mirror_path() {
    let steps = 9;
    let est = run(rng_for(2), steps, 200, |normals| {
        normals.iter().sum::<f64>()
    });
    assert!(
        est.standard_error < 1e-12,
        "antithetic SE for an odd payoff should vanish, got {}",
        est.standard_error
    );
    //Var(sum of 9 standard normals) = 9, 400 independent paths -> SE = sqrt(9 / 400).
    let expected_naive = (9.0f64 / 400.0).sqrt();
    assert!(
        (est.naive_standard_error - expected_naive).abs() <= 0.1 * expected_naive,
        "naive SE {} vs theory {expected_naive}",
        est.naive_standard_error
    );
    assert!(est.variance_reduction > 1e6);
}

/// A convex-but-monotone payoff is where pairing is supposed to pay: the mirrored leg sits on the
/// other side of the mean, so the pair average varies less than either path.
#[test]
fn pairing_reduces_the_error_of_a_monotone_payoff() {
    let est = run(rng_for(3), 1, 4000, |normals| {
        let z = normals[0];
        //Strictly increasing in z, so cov(z, f) > 0 and the pair variance is below the path
        //variance by construction.
        0.05 + 0.02 * z + 0.001 * z * z
    });
    assert!(
        est.variance_reduction > 1.5,
        "expected a meaningful variance reduction, got {}",
        est.variance_reduction
    );
    //E[f] = 0.05 + 0.001.
    assert!(
        (est.mean - 0.051).abs() <= 5.0 * est.standard_error,
        "mean {} is not within 5 SE of 0.051",
        est.mean
    );
}

/// The grid rule: integer steps, exact endpoints, `dt = span / steps`.
#[test]
fn grids_land_on_every_leg_endpoint() {
    let grid = Grid::for_legs(0.0, &[1.5, 1.75], 100);
    assert_eq!(grid.legs().len(), 2);
    assert_eq!(grid.increments(), 175);
    for leg in grid.legs() {
        assert_eq!(leg.time(0), leg.start);
        assert!(
            (leg.time(leg.steps) - leg.end).abs() <= 1e-12,
            "leg {} -> {} ended at {}",
            leg.start,
            leg.end,
            leg.time(leg.steps)
        );
        assert!(leg.dt() <= 1.0 / 100.0 + 1e-12, "dt overshot the request");
    }
    assert!(
        (grid.legs()[0].dt() - 0.01).abs() <= 1e-12,
        "{}",
        grid.legs()[0].dt()
    );
    assert!(
        (grid.legs()[1].dt() - 0.01).abs() <= 1e-12,
        "{}",
        grid.legs()[1].dt()
    );
}

/// A doubling of the resolution is an exact halving of `dt` on every leg, which is what makes the
/// Richardson factor in the bias estimate meaningful.
#[test]
fn refining_the_resolution_gives_a_uniform_step_ratio() {
    let coarse = Grid::for_legs(0.0, &[1.5, 1.75], 100);
    let fine = Grid::for_legs(0.0, &[1.5, 1.75], 200);
    let ratio = coarse.step_ratio(&fine);
    assert!(
        (ratio - 2.0).abs() <= 1e-9,
        "expected a factor of exactly 2, got {ratio}"
    );
}

/// The bias arithmetic, on numbers with a known answer: a bias of `c * dt` with the fine grid at
/// half the step size leaves a residual of exactly the observed gap on the fine run.
#[test]
fn richardson_residual_matches_first_order_bias() {
    //Suppose the true (continuous) value is 1.0 and the bias is -0.02 * dt.
    let continuous = 1.0;
    let dt_coarse = 0.01;
    let dt_fine = 0.005;
    let coarse = Estimate {
        mean: continuous - 0.02 * dt_coarse,
        standard_error: 1e-6,
        naive_standard_error: 2e-4,
        variance_reduction: 2.0,
        pairs: 1000,
        paths: 2000,
    };
    let fine = Estimate {
        mean: continuous - 0.02 * dt_fine,
        standard_error: 1e-6,
        naive_standard_error: 2e-4,
        variance_reduction: 2.0,
        pairs: 1000,
        paths: 2000,
    };
    let bias = discretisation_bias(&coarse, &fine, dt_coarse / dt_fine);
    assert_eq!(bias.ratio, 2.0);
    //gap = 0.02 * (dt_coarse - dt_fine) = 1e-4, residual = gap / (ratio - 1) = 1e-4.
    assert!((bias.gap - 1e-4).abs() <= 1e-12, "{}", bias.gap);
    assert!(
        (bias.residual - (continuous - fine.mean).abs()).abs() <= 1e-12,
        "residual {} should equal the true fine-run bias {}",
        bias.residual,
        continuous - fine.mean
    );
    assert!(
        (bias.extrapolated - continuous).abs() <= 1e-12,
        "extrapolation {} should land on the continuous value {continuous}",
        bias.extrapolated
    );
    assert!(bias.resolvable());
}

/// When the refinement response is below the noise in the difference, the residual is floored at
/// that noise, not at the (smaller) point estimate.
#[test]
fn the_residual_is_floored_at_the_noise_floor_of_the_gap() {
    let coarse = sample_estimate(1.0, 1e-3);
    let fine = sample_estimate(1.0 + 1e-5, 1e-3);
    let bias = discretisation_bias(&coarse, &fine, 2.0);
    let gap_se = (2.0f64.sqrt()) * 1e-3;
    assert!(
        (bias.residual - gap_se).abs() <= 1e-12,
        "residual {} should be the gap SE {gap_se}, not the 1e-5 point estimate",
        bias.residual
    );
}

/// A refinement response smaller than the noise in the difference is reported as unresolved, not
/// as a measured bias.
#[test]
fn an_unresolved_refinement_is_flagged_not_silently_used() {
    let coarse = sample_estimate(1.0, 1e-3);
    let fine = sample_estimate(1.0 + 1e-5, 1e-3);
    let bias = discretisation_bias(&coarse, &fine, 2.0);
    assert!(
        !bias.resolvable(),
        "a 1e-5 gap with a 1.4e-3 gap SE is not resolved: {bias:?}"
    );
    //Still usable as a bound; it just is not a measurement.
    assert!(bias.residual > 0.0);
}

fn sample_estimate(mean: f64, standard_error: f64) -> Estimate {
    Estimate {
        mean,
        standard_error,
        naive_standard_error: standard_error * 2.0,
        variance_reduction: 2.0,
        pairs: 1000,
        paths: 2000,
    }
}

/// The budget is the sum of its stated parts, and the check's deviation is measured against that
/// sum rather than against a literal.
#[test]
fn the_budget_is_variance_plus_the_stated_bias() {
    let est = sample_estimate(0.025, 1e-4);
    let check = Check {
        label: "test check",
        estimate: est,
        analytic: 0.0255,
        k: K_SIGMA,
        bias: 3e-4,
        discretisation: None,
    };
    assert!(
        (check.budget() - (K_SIGMA * 1e-4 + 3e-4)).abs() <= 1e-15,
        "{}",
        check.budget()
    );
    assert!(
        (check.deviation() - 5e-4).abs() <= 1e-15,
        "{}",
        check.deviation()
    );
    //5e-4 <= 7e-4: inside budget even though it is 5 SE of the analytic value, because the
    //discretisation term was declared.
    check.assert_within_budget();
}

#[test]
#[should_panic(expected = "budget exceeded")]
fn a_deviation_outside_the_budget_fails() {
    let est = sample_estimate(0.03, 1e-4);
    Check {
        label: "test check",
        estimate: est,
        analytic: 0.025,
        k: K_SIGMA,
        bias: 3e-4,
        discretisation: None,
    }
    .assert_within_budget();
}

/// Every test gets its own stream: the same slot reproduces, different slots do not.
#[test]
fn seeds_are_reproducible_and_slots_are_independent() {
    let normal = StandardNormal;
    let draw = |slot: u8| {
        let mut rng = rng_for(slot);
        (0..4)
            .map(|_| normal.sample(&mut rng))
            .collect::<Vec<f64>>()
    };
    assert_eq!(draw(1), draw(1), "a slot must reproduce exactly");
    assert_ne!(draw(1), draw(2), "slots must not share their stream");
    assert_ne!(draw(2), draw(3), "and neighbouring slots must not either");
}

/// The increments themselves really are standard normal: mean `0`, variance `1`, each at four
/// standard errors of that moment's own estimate.
///
/// Worth having next to the seed test because it covers the failure the price checks cannot.  Every
/// instrument-level band is sized from the dispersion the run just measured, so a draw stream that
/// is *not* the model's `N(0, 1)` — a distribution swapped in at an import site, a `Normal`
/// built with a variance of 2, a `sample` call wired to the wrong generator — would not show up
/// as a blown band on some price.  It would show up as every band being self-consistently wrong,
/// with each simulated price drifting from its closed form by an amount the reported SE happily
/// absorbs.  This is the one place the primitive is checked rather than the estimator built on it.
///
/// For a unit-normal sample, `SE(mean) = 1 / sqrt(N)` and `SE(var) = sqrt(2 / N)` (the variance
/// estimate's own standard error, `sqrt((kurtosis - 1) / N)` at kurtosis 3).  Four of each is
/// the same `k` the price checks run at, so each moment is a ~1-in-15,000 event under the null
/// that the draws are what they claim to be.  Measured here, at `N = 200_000` on slot 4: `mean
/// = -1.4e-3` and `variance = 1.0055`, against tolerances of `8.9e-3` and `1.3e-2`
/// respectively — each moment using about an eighth of the band it is given.
#[test]
fn the_increment_draws_are_standard_normal() {
    const N: usize = 200_000;
    //A dedicated slot: nothing else in the crate draws from these numbers, so the moments below are
    //a fixed function of this test, reproducible run to run.
    let mut rng = rng_for(4);
    let normal = StandardNormal;
    //Online (Welford) moments, so the check does not hold 200k draws alive.
    let mut mean = 0.0;
    let mut m2 = 0.0;
    for i in 1..=N {
        let z: f64 = normal.sample(&mut rng);
        assert!(z.is_finite(), "draw {i} of {N} was not finite: {z}");
        let delta = z - mean;
        mean += delta / i as f64;
        m2 += delta * (z - mean);
    }
    let variance = m2 / (N - 1) as f64;
    let se_mean = 1.0 / (N as f64).sqrt();
    let se_variance = (2.0 / N as f64).sqrt();

    assert!(
        mean.abs() <= K_SIGMA * se_mean,
        "mean {mean:+.5} is more than {K_SIGMA} SE ({:.5}) from 0: the draws are not centred",
        K_SIGMA * se_mean,
    );
    assert!(
        (variance - 1.0).abs() <= K_SIGMA * se_variance,
        "variance {variance:.5} is more than {K_SIGMA} SE ({:.5}) from 1: the draws are not unit \
         scale (a std-dev of {} would produce this)",
        K_SIGMA * se_variance,
        variance.sqrt(),
    );
}

/// With zero noise the simulated path is the model's own mean path, which makes the walk checkable
/// exactly: the rate at each leg boundary is the chained conditional mean *at that date* (not one
/// step past it, which is the bug the grid convention exists to prevent), and splitting a span into
/// legs does not move the grid, drop a step, or change the integral.
#[test]
fn the_walk_lands_on_its_leg_boundaries_and_covers_its_whole_span() {
    use crate::HullWhite;
    use crate::test_support::STEEP_CURVE;
    let curve = STEEP_CURVE.curve();
    let model = HullWhite::new(STEEP_CURVE.a, STEEP_CURVE.sigma, &curve).unwrap();
    let start = 0.02f64;
    let grid = Grid::for_legs(0.0, &[1.5, 1.75], 100);
    let zeros = vec![0.0; grid.increments()];
    let mut state = PathState::new(&grid);
    walk(&model, &grid, start, &zeros, &mut state);

    assert_eq!(state.rates[0], start);
    let mut chained = start;
    for (index, leg) in grid.legs().iter().enumerate() {
        chained = model.mu_r(chained, leg.start, leg.end).unwrap();
        assert!(
            (state.rate_at_end_of_leg(index) - chained).abs() <= 1e-12 * chained.abs(),
            "leg {index} ended at {} but the model's mean at {} is {chained}",
            state.rate_at_end_of_leg(index),
            leg.end
        );
    }

    //The same span as one leg at the same resolution is the same grid, so it must integrate to the
    //same number: no step gained, lost or shifted at the join between legs.
    let one_leg = Grid::for_legs(0.0, &[1.75], 100);
    assert_eq!(one_leg.increments(), grid.increments());
    assert!(
        (one_leg.legs()[0].dt() - grid.legs()[0].dt()).abs() <= 1e-15,
        "same steps_per_year should give the same dt"
    );
    let mut single_state = PathState::new(&one_leg);
    walk(
        &model,
        &one_leg,
        start,
        &vec![0.0; one_leg.increments()],
        &mut single_state,
    );
    assert!(
        (single_state.integral - state.integral).abs() <= 1e-12 * state.integral.abs(),
        "two-leg integral {} vs one-leg {}",
        state.integral,
        single_state.integral
    );

    //And the integral is a right-endpoint sum over the full span: recomputed independently here,
    //it has to agree to rounding.
    let mut r = start;
    let mut expected = 0.0;
    for leg in grid.legs() {
        let dt = leg.dt();
        for step in 0..leg.steps {
            let t = leg.time(step);
            r = model.mu_r(r, t, t + dt).unwrap();
            expected += r * dt;
        }
    }
    assert!(
        (state.integral - expected).abs() <= 1e-12 * expected.abs(),
        "walk integral {} vs right-endpoint sum {expected}",
        state.integral
    );
}
