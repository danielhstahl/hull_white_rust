//! Monte-Carlo test harness: fixed seeds, antithetic variates, and an explicit error budget.
//!
//! Every Monte-Carlo test in this crate answers the same question — "does the simulated price of
//! this instrument agree with the closed form?" — and that question has two separate sources of
//! disagreement, which need to be budgeted separately:
//!
//! * **Sampling error.** The estimate is an average over a finite number of paths, so it carries
//!   a standard error that shrinks like `1 / sqrt(n_paths)`.  Nothing about the model changes it;
//!   it is the irreducible noise of the estimator.
//! * **Discretisation bias.** Where a path is built from Euler steps, the simulated integral of
//!   the short rate is not the real integral, so the simulated price converges to something
//!   *near* the closed form rather than *to* it.  Halving `dt` halves this error; at a finite
//!   step count it never vanishes.
//!
//! A hard-coded `epsilon = 0.0001` conflates the two, and it is a number that can be *smaller*
//! than the sampling error of the run it guards.  Measured on
//! [`crate::test_support::BASELINE`]: the caplet's per-path discounted payoff has a standard
//! deviation of `3.2e-3` against a mean of `2.9e-3`, so the pre-harness check — which averaged
//! 1,000 such paths — had a standard error of about `1.0e-4`, i.e. **one standard error of
//! sampling noise was the size of the whole tolerance**.  A reseed, a `rand` bump, or a codegen
//! change that moved the draw stream could have failed that assertion with every price in the
//! crate unchanged.  This module replaces the constant with a budget computed from the run:
//!
//! ```text
//! | mc_mean - analytic |  <=  k * standard_error             (sampling, from this very run)
//!                       +     residual_discretisation_bias    (from a coarse/fine pair; 0 if none)
//! ```
//!
//! with [`K_SIGMA`] documented below.  A failure then means one of two things, and the message
//! says which: the estimate moved further than its own noise can explain, or the discretisation
//! residual grew.  Neither can be made to go away by editing a constant, because the constant is
//! not there — tightening the band means buying paths or refining the grid, and both are visible in
//! the report the assertion prints.
//!
//! # Antithetic variates
//!
//! [`AntitheticMonteCarlo`] draws the standard normal increments for a whole path and evaluates
//! the payoff twice: once on the drawn increments `z`, once on the mirrored `-z`.  The estimator
//! is the mean over the `n_pairs` *pair averages*, so its standard error comes from the dispersion
//! of the pair averages, which is exactly where the negative correlation between the two legs
//! shows up as a discount.  The cost is one extra payoff evaluation per pair — free next to the
//! path simulation — and the same pairing also makes the estimator's sampling distribution much
//! more symmetric than the per-path payoff distribution, which is what makes a Gaussian
//! `k * SE` band worth asserting on at a few thousand pairs instead of a few hundred thousand.
//!
//! Measured on [`crate::test_support::BASELINE`], `T = 1.5`, `delta = 0.25`, against the
//! unpaired estimator at the same path count (`mc::tests` asserts the direction of each):
//!
//! ```text
//! instrument                                  variance reduction (unpaired SE / paired SE)
//! Eurodollar future, one draw per path        ~208x  (payoff nearly affine in z over the range)
//! caplet,   K = 0.02  (in the money)          ~1.9x
//! floorlet, K = 0.05  (the valuable side)     ~3.8x
//! ```
//!
//! The Eurodollar number is large because that estimator draws a single normal per path and the
//! payoff `1 / B(T, T + delta)` is close to affine in it; the caplet and floorlet gain less but
//! still buy nearly a factor of two and four in effective sample size.
//!
//! # Step-size convention
//!
//! The single rule, used by every simulated path in this crate:
//!
//! > A simulated interval `[t_start, t_end]` is split into an **integer number of steps** `n`,
//! > with `dt = (t_end - t_start) / n`.  `n` always counts *steps*, never grid points.
//!
//! The two older conventions these tests grew up with, in the same test file, were
//!
//! * `dt = span / (n - 1)` with `n` draws, which walks one whole step *past* the interval's end —
//!   the caplet and floorlet checks read their Libor fixing off the rate at `T + dt`, one step
//!   after the fixing date; and
//! * the arrears leg's `dt = span / (n - 1)` with the *first* step skipped, which left the
//!   discount integral covering `[T + dt, T + delta + dt]` instead of `[T, T + delta]`: right in
//!   length, wrong in position, by one step at each end.
//!
//! Both are off by one in a way that is invisible at small `dt`, that differs between tests, and
//! that nobody can check by reading the test.  The rule above is stated once here and enforced by
//! [`Grid`] rather than retyped per test: a grid endpoint *is* the event date.
//!
//! A path that has to read a rate at an intermediate model event — a caplet fixes Libor at its
//! own expiry and then keeps discounting to the payment date — is simulated **leg by leg**, each
//! leg getting its own integer step count, so no event ever falls between two grid points.
//! [`Grid::for_legs`] builds that from one `steps_per_year` figure, and
//! [`Grid::refined`] produces the finer grid for the bias study, so the two resolutions differ by
//! an exact factor in `dt` *on every leg*.
//!
//! The state update inside a step is the same backward (right-hand) Euler rule the pre-existing
//! tests used — step the rate with the model's own exact conditional moments over the step, then
//! hold that *new* rate across the step for the running `integral r ds`:
//!
//! ```text
//! r_{i+1} = mu_r(r_i, t_i, t_{i+1}) + sqrt(var_r(t_i, t_{i+1})) * z_i
//! integral += r_{i+1} * dt
//! ```
//!
//! Note that the transition itself is *exact* — `mu_r` and `variance_r` are the model's own
//! conditional mean and variance over `[t_i, t_{i+1}]`, so the sampled rate is right at the grid
//! points and the only discretisation lives in the integral.  That is why the residual in the
//! table below is a couple of parts in `1e-5` on a price rather than the much larger number a
//! hand-rolled `r += a(b - r)dt + sigma sqrt(dt) z` Euler loop produces at the same `dt`.
//!
//! Measured, at the heavy tier (fine `dt = 1/400`, coarse `dt = 1/200`, `BASELINE`,
//! `T = 1.5`, `delta = 0.25`; 30,000 antithetic pairs for the two path checks, `10 x` that for
//! the Eurodollar check).  Re-measured after the `rand` 0.5 -> 0.9 bump, which replaced the
//! whole draw stream (`ChaCha20` in `rand 0.5` -> `ChaCha12` in `rand_chacha 0.9`, normals
//! from `rand_distr::StandardNormal` instead of `rand::distributions::StandardNormal`) and not
//! one tolerance: every number below is printed by the run it describes.
//!
//! ```text
//! instrument     |mc - analytic|   k * SE      residual     budget     budget used
//! caplet          3.7e-6          2.7e-5      9.7e-6       3.7e-5      9.9%
//! floorlet        8.1e-7          1.8e-5      6.5e-6       2.5e-5      3.3%
//! EDF (exact)     1.9e-7          4.8e-7      0            4.8e-7     40.4%
//! ```
//!
//! That the same bands absorb a completely different stream — same orders of magnitude, same
//! residuals, nothing retuned — is the property this harness exists to provide.
//!
//! The two path-based checks use roughly a tenth of their budget; the discretisation residual is a
//! quarter to a third of the budget, which is the sign that the grid is fine enough that sampling
//! error, not discretisation, is what the band is really made of.  The Eurodollar check has no
//! discretisation term at all and is the tightest of the three.
//!
//! # Bias estimation
//!
//! [`discretisation_bias`] runs the Richardson argument on a coarse/fine pair under the stated
//! first-order assumption: with step ratio `ratio = dt_coarse / dt_fine > 1` and a bias `c * dt`,
//!
//! ```text
//! residual bias of the fine run = |fine_mean - coarse_mean| / (ratio - 1)
//! extrapolated (bias-free) value = fine_mean + (fine_mean - coarse_mean) / (ratio - 1)
//! ```
//!
//! The residual goes into the assertion's budget next to `k * SE`, and the whole breakdown is
//! reported, so the size of the discretisation term relative to the noise term is visible.
//!
//! The residual is `max(|gap|, gap standard error) / (ratio - 1)`, not the bare point estimate:
//! the Richardson number is itself noisy, and a run whose refinement response happens to come out
//! small must not be allowed to claim a bias allowance tighter than its own noise permits.  A run
//! whose two step sizes differ by less than the noise in that difference
//! ([`Discretisation::resolvable`]) has a bias this pair cannot actually see; the report says so
//! rather than pretending the number was measured.
//!
//! Tests whose estimator has **no** discretisation at all pass [`NO_DISCRETISATION`] and get a
//! pure `k * SE` band.  Saying so explicitly is the point: it forces the question of which term
//! belongs in the budget for which instrument.  The Eurodollar-future check is the case in point —
//! it draws the terminal rate from its exact conditional distribution, so the expectation of
//! `1 / B(T, T + delta)` over that one draw *is* the model expectation and there is no time
//! discretisation to bound: its budget is sampling error and nothing else.
//!
//! # Choosing `k`
//!
//! [`K_SIGMA`] is 4.  The pair-averaged payoff is close enough to Gaussian at these sample sizes
//! that `k = 4` is a two-sided roughly 1-in-15,000 event per assertion under the null that the
//! simulator and the closed form agree.  It is deliberately not 2 (~1-in-20 per assertion, which
//! fails on ordinary noise) and deliberately not 8 (~1-in-10^15, a band so wide that a real sign
//! error or a wrong-leg payoff would pass inside it).  A test that needs a different `k` puts the
//! reason next to it.
//!
//! # Run cost, and how to run this
//!
//! Sizes are tiered by the `slow` feature ([`scale`]), so the heavy end of the error-budget curve
//! is opt-in and `cargo test` stays quick:
//!
//! ```text
//! cargo test                                  # quick tier: ~0.3s of Monte Carlo in total
//! cargo test --features slow                  # heavy tier: ~18s, bands ~4x tighter
//! cargo test --features slow -- --nocapture   # the full per-run error-budget report on stdout
//! ```
//!
//! CI runs both tiers on every push: `.github/workflows/test.yml` keeps plain `cargo test` for the
//! fast signal and adds a `cargo test --features slow` step, and the nightly coverage run uses
//! `--all-features`, which includes the heavy tier.  The reason for tiering rather than one size
//! is that the sample size worth having under a `1e-5` budget is four orders of magnitude above
//! the one worth having under a `1e-3` budget, and paying that on every save is what the
//! pre-harness version of these three tests cost — about a third of the whole `cargo test --lib`
//! wall clock, for a tolerance no tighter than the quick tier's.

use rand::SeedableRng;
use rand::rngs::StdRng;
use rand_distr::{Distribution, StandardNormal};

use crate::HullWhite;

/// Two-sided deviation multiplier for Monte-Carlo assertions.  See the `k` section of the module
/// docs: 4 is ~1-in-15,000 per assertion under the null hypothesis that the simulator and the
/// closed form agree.
pub const K_SIGMA: f64 = 4.0;

/// The budget entry for an estimator that has no discretisation bias because it does not step
/// through time at all.
pub const NO_DISCRETISATION: f64 = 0.0;

/// The two sample-size tiers, selected by the `slow` feature.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Scale {
    /// Antithetic pairs per Monte-Carlo estimate.  Paths simulated = `2 * pairs`.
    pub pairs: usize,
    /// Euler steps per year for the *fine* run of a bias study.  The coarse run is exactly half.
    pub steps_per_year: usize,
}

/// The tier this build was compiled for.
///
/// Quick: a `k * SE` band of order `1e-4` on a rate-price — the same order the hand-written
/// `epsilon = 0.0001` these tests replace was aiming at, but *earned* rather than asserted, and
/// about a tenth of the cost.  Heavy: ~20x the pairs on a 4x finer grid, which tightens the band
/// to `~1e-5` and takes the discretisation residual from "noise-dominated bound" to "measured".
///
/// ```text
///                   pairs   steps/yr    dt         total budget            all-MC wall clock
///                                              caplet    floorlet   EDF
/// quick             1 500     100      1/100 y   1.7e-4    1.1e-4   2.2e-6   ~0.3s
/// heavy (`slow`)   30 000     400      1/400 y   3.7e-5    2.5e-5   4.8e-7   ~18s
/// ```
///
/// Measured on `BASELINE`, `T = 1.5`, `delta = 0.25`, debug profile, the run the reports print:
/// the heavy tier buys about 4x the tightness of the band for 60x the Monte-Carlo time, which is
/// precisely why it is the tier that CI runs rather than the tier a save waits on.
pub const fn scale() -> Scale {
    if cfg!(feature = "slow") {
        Scale {
            pairs: 30_000,
            steps_per_year: 400,
        }
    } else {
        Scale {
            pairs: 1_500,
            steps_per_year: 100,
        }
    }
}

/// A deterministic RNG for one Monte-Carlo test.
///
/// Each test passes its own `slot`, so no two tests draw from the same stretch of the same stream:
/// adding, removing or changing one test cannot shift another's numbers, and every test is
/// bit-reproducible regardless of `--test-threads`.  The stream for a slot is a splitmix64 fill of
/// the slot index, so neighbouring slots are not trivially-related seeds.
pub fn rng_for(slot: u8) -> StdRng {
    let mut seed = [0u8; 32];
    let mut state = (slot as u64) ^ 0x9E37_79B9_7F4A_7C15;
    for byte in &mut seed {
        state ^= state >> 30;
        state = state.wrapping_mul(0xBF58_476D_1CE4_E5B9);
        state ^= state >> 27;
        state = state.wrapping_mul(0x94D0_49BB_1331_11EB);
        *byte = (state >> 33) as u8;
    }
    <StdRng as SeedableRng>::from_seed(seed)
}

/// The result of averaging a payoff over a finite sample of antithetic pairs.
#[derive(Debug, Clone, Copy)]
pub struct Estimate {
    /// Mean of the pair averages: the estimator proper.
    pub mean: f64,
    /// Standard error of [`Estimate::mean`] — the sampling half of the error budget.
    pub standard_error: f64,
    /// Standard error the same number of *independent* paths would have had, i.e. without the
    /// variance reduction.  Only used to report what the pairing bought.
    pub naive_standard_error: f64,
    /// `naive_standard_error / standard_error`: how many times as many paths an unpaired estimator
    /// would need for the same precision.
    pub variance_reduction: f64,
    pub pairs: usize,
    pub paths: usize,
}

impl Estimate {
    /// The sampling half of the budget: `k` standard errors.
    pub fn band(&self, k: f64) -> f64 {
        k * self.standard_error
    }

    /// Push the estimator through an affine map `x -> scale * x + shift`.
    ///
    /// Needed because a simulated quantity is often not the quoted one: the Eurodollar-future
    /// check simulates `E[1 / B(T, T + delta)]` and the future is `(E[1 / B] - 1) / delta`.  An
    /// affine map moves the mean by `shift`, scales the standard error by `|scale|`, and leaves
    /// the variance-reduction ratio (a ratio of two equally-scaled things) alone.
    pub fn rescale(&self, scale: f64, shift: f64) -> Estimate {
        Estimate {
            mean: scale * self.mean + shift,
            standard_error: self.standard_error * scale.abs(),
            naive_standard_error: self.naive_standard_error * scale.abs(),
            variance_reduction: self.variance_reduction,
            pairs: self.pairs,
            paths: self.paths,
        }
    }
}

/// Accumulator for antithetic pairs.
///
/// Build one with a seeded RNG and the number of normal increments a single path consumes, then
/// [`AntitheticMonteCarlo::push_pair`] repeatedly.  Statistics are accumulated online (Welford)
/// over the pair averages and, separately, over the individual path payoffs, so nothing but the
/// final numbers is retained and the pair values never have to be held in memory.
pub struct AntitheticMonteCarlo {
    rng: StdRng,
    normals: Vec<f64>,
    pairs: usize,
    pair_mean: f64,
    pair_m2: f64,
    path_mean: f64,
    path_m2: f64,
    paths: usize,
}

impl AntitheticMonteCarlo {
    /// `increments` is how many standard normal draws one complete path consumes: the total step
    /// count of the [`Grid`] the payoff walks.
    pub fn new(rng: StdRng, increments: usize) -> Self {
        assert!(increments > 0, "a path needs at least one increment");
        Self {
            rng,
            normals: vec![0.0; increments],
            pairs: 0,
            pair_mean: 0.0,
            pair_m2: 0.0,
            path_mean: 0.0,
            path_m2: 0.0,
            paths: 0,
        }
    }

    /// Push one antithetic pair: `payoff` is called with the drawn increments `z`, then with every
    /// increment negated, and must return that path's (already discounted, already scaled) payoff.
    ///
    /// The payoff closure is handed the *whole* path's increments rather than being stepped one at
    /// a time, because the mirror path is only available if the increments are all materialised:
    /// an increment-at-a-time API could not reuse the same draws as their own negation.
    pub fn push_pair(&mut self, mut payoff: impl FnMut(&[f64]) -> f64) {
        for slot in &mut self.normals {
            *slot = StandardNormal.sample(&mut self.rng);
        }
        let first = payoff(&self.normals);
        for slot in &mut self.normals {
            *slot = -*slot;
        }
        let second = payoff(&self.normals);
        assert!(
            first.is_finite() && second.is_finite(),
            "payoff produced a non-finite value ({first}, {second})"
        );

        let pair = 0.5 * (first + second);
        self.pairs += 1;
        let delta = pair - self.pair_mean;
        self.pair_mean += delta / self.pairs as f64;
        self.pair_m2 += delta * (pair - self.pair_mean);

        for value in [first, second] {
            self.paths += 1;
            let delta = value - self.path_mean;
            self.path_mean += delta / self.paths as f64;
            self.path_m2 += delta * (value - self.path_mean);
        }
    }

    /// Finish and reduce to an [`Estimate`].
    pub fn finish(self) -> Estimate {
        assert!(
            self.pairs >= 2,
            "need at least two pairs for a sample standard error"
        );
        let pair_variance = self.pair_m2 / (self.pairs - 1) as f64;
        let path_variance = self.path_m2 / (self.paths - 1) as f64;
        let standard_error = (pair_variance / self.pairs as f64).sqrt();
        let naive_standard_error = (path_variance / self.paths as f64).sqrt();
        Estimate {
            mean: self.pair_mean,
            standard_error,
            naive_standard_error,
            variance_reduction: if standard_error > 0.0 {
                naive_standard_error / standard_error
            } else {
                f64::INFINITY
            },
            pairs: self.pairs,
            paths: self.paths,
        }
    }
}

/// Run `pairs` antithetic pairs of an `increments`-draw path through `payoff`.
pub fn run(
    rng: StdRng,
    increments: usize,
    pairs: usize,
    mut payoff: impl FnMut(&[f64]) -> f64,
) -> Estimate {
    let mut mc = AntitheticMonteCarlo::new(rng, increments);
    for _ in 0..pairs {
        mc.push_pair(&mut payoff);
    }
    mc.finish()
}

/// A uniform time grid over one interval, in the convention stated in the module docs: an integer
/// step count, `dt = span / steps`, both endpoints exactly on the grid.
#[derive(Debug, Clone, Copy)]
pub struct Leg {
    /// First time of the leg (a grid point).
    pub start: f64,
    /// Last time of the leg (a grid point).
    pub end: f64,
    /// Number of steps in the leg, not of points.
    pub steps: usize,
}

impl Leg {
    /// The leg's step size.
    pub fn dt(&self) -> f64 {
        (self.end - self.start) / self.steps as f64
    }

    /// Time of grid point `i` (0-based, `0 <= i <= steps`).
    pub fn time(&self, i: usize) -> f64 {
        self.start + self.dt() * i as f64
    }
}

/// A path's worth of legs, each uniformly stepped.
///
/// Built from a `steps_per_year` figure so that two resolutions of the same instrument differ by an
/// exact factor in `dt` on every leg — which is what makes the Richardson factor in
/// [`discretisation_bias`] exactly that factor rather than something per-leg.  Use
/// [`Grid::refined`] rather than a second [`Grid::for_legs`] call for the fine grid of a bias study.
#[derive(Debug, Clone)]
pub struct Grid {
    legs: Vec<Leg>,
}

impl Grid {
    /// Step through `ends` in order, starting at `start`, at `steps_per_year`: each leg gets
    /// `ceil(len * steps_per_year)` steps, and at least one.
    ///
    /// Rounding *up* rather than to nearest is deliberate: it guarantees the resolved `dt` is no
    /// larger than the requested one, so a test asking for a resolution never silently gets a
    /// coarser path than the figure it quoted, and the bias bound it computes is a bound on a path
    /// at least as fine as the one it names.
    pub fn for_legs(start: f64, ends: &[f64], steps_per_year: usize) -> Self {
        assert!(
            steps_per_year > 0 && !ends.is_empty(),
            "a grid needs a positive resolution and at least one leg end"
        );
        let mut cursor = start;
        let mut legs = Vec::with_capacity(ends.len());
        for &end in ends {
            assert!(
                end > cursor,
                "grid legs must advance: {end} is not after {cursor}"
            );
            let len = end - cursor;
            let steps = ((len * steps_per_year as f64).ceil() as usize).max(1);
            legs.push(Leg {
                start: cursor,
                end,
                steps,
            });
            cursor = end;
        }
        Self { legs }
    }

    /// The legs, in time order.
    pub fn legs(&self) -> &[Leg] {
        &self.legs
    }

    /// Total steps over all legs, i.e. the number of normal increments one path through this grid
    /// consumes.
    pub fn increments(&self) -> usize {
        self.legs.iter().map(|leg| leg.steps).sum()
    }

    /// This grid with every leg's step count multiplied by `factor`, endpoints unchanged.
    ///
    /// This, not a second [`Grid::for_legs`] at a different `steps_per_year`, is how a bias study
    /// gets its fine grid: the refinement is exact on every leg, so the step ratio is the factor
    /// exactly rather than something that depends on how each leg's length happened to round.  Two
    /// independently rounded grids over the same legs generally do **not** share a leg-to-leg
    /// ratio, which makes the Richardson factor meaningless — [`Grid::step_ratio`] refuses them.
    #[must_use = "refined returns the finer grid"]
    pub fn refined(&self, factor: usize) -> Grid {
        assert!(factor >= 2, "a refinement factor must be at least 2");
        Grid {
            legs: self
                .legs
                .iter()
                .map(|leg| Leg {
                    start: leg.start,
                    end: leg.end,
                    steps: leg.steps * factor,
                })
                .collect(),
        }
    }

    /// `dt_coarse / dt_fine` between this grid and a finer one over the same legs.
    ///
    /// Requires the ratio to be uniform across legs to within a rounding of the step counts; a
    /// non-uniform ratio would make the Richardson factor ill-defined, so it is refused rather
    /// than averaged into something misleading.
    pub fn step_ratio(&self, fine: &Grid) -> f64 {
        assert_eq!(
            self.legs.len(),
            fine.legs.len(),
            "step ratio needs the same leg structure"
        );
        let mut ratio = None;
        for (coarse, fine_leg) in self.legs.iter().zip(&fine.legs) {
            assert_eq!(
                coarse.start, fine_leg.start,
                "step ratio needs the same leg boundaries"
            );
            assert_eq!(
                coarse.end, fine_leg.end,
                "step ratio needs the same leg boundaries"
            );
            let this = coarse.dt() / fine_leg.dt();
            match ratio {
                None => ratio = Some(this),
                Some(previous) => assert!(
                    (this - previous).abs() <= 1e-6 * previous.abs(),
                    "step ratio is not uniform across legs: {previous} then {this}"
                ),
            }
        }
        ratio.unwrap_or(1.0)
    }
}

/// What a simulated path leaves behind: the short rate at every leg boundary, and the running
/// integral of the short rate over the whole grid.
///
/// `rates[0]` is the rate the path started from, `rates[i + 1]` the rate at the end of leg `i`.
/// For a simple-rate option built as `[t, T]` then `[T, T + delta]`, the fixing rate is
/// `rates[1]`, which is the rate *at* `T` rather than one step past it.
#[derive(Debug, Clone)]
pub struct PathState {
    /// Short rate at each leg boundary, start first.
    pub rates: Vec<f64>,
    /// Right-endpoint Riemann integral of the short rate over the whole grid.
    pub integral: f64,
}

impl PathState {
    /// An empty state sized for `grid`, reusable across paths.
    pub fn new(grid: &Grid) -> Self {
        Self {
            rates: vec![0.0; grid.legs().len() + 1],
            integral: 0.0,
        }
    }

    /// The rate at the end of leg `index`.
    pub fn rate_at_end_of_leg(&self, index: usize) -> f64 {
        self.rates[index + 1]
    }
}

/// Walk `grid` from `start_rate` with the model's own exact conditional moments, consuming one
/// normal per step and accumulating `integral r ds` by the right-endpoint rule.
///
/// The rate at each grid point is exact (the transition kernel is the model's, not a hand-rolled
/// Euler step), so the only discretised quantity is the integral.  See the module docs for the
/// convention and for why the *fixing* lands exactly on the leg boundary.
pub fn walk(
    model: &HullWhite,
    grid: &Grid,
    start_rate: f64,
    normals: &[f64],
    state: &mut PathState,
) {
    assert_eq!(
        normals.len(),
        grid.increments(),
        "the walk needs exactly one normal per grid step"
    );
    state.rates[0] = start_rate;
    let mut rate = start_rate;
    let mut integral = 0.0;
    let mut offset = 0usize;
    for (index, leg) in grid.legs().iter().enumerate() {
        let dt = leg.dt();
        for step in 0..leg.steps {
            let t = leg.time(step);
            //Each step conditions on the rate the previous step produced, which is what makes the
            //rate at a grid point exact; only the integral below is a Riemann sum.
            rate = model.mu_r(rate, t, t + dt).unwrap()
                + model.variance_r(t, t + dt).unwrap().sqrt() * normals[offset];
            integral += rate * dt;
            offset += 1;
        }
        //The rate at the leg boundary is the rate *at* that date, not one step past it.
        state.rates[index + 1] = rate;
    }
    state.integral = integral;
}

/// Run `pairs` antithetic pairs of paths through `grid`, handing each completed [`PathState`] to
/// `payoff`.
///
/// This is the entry point an instrument-level Monte-Carlo check uses: the walk convention lives in
/// the harness, so a test says what is paid and when and cannot accidentally invent its own
/// `dt`.
pub fn run_paths(
    model: &HullWhite,
    rng: StdRng,
    grid: &Grid,
    start_rate: f64,
    pairs: usize,
    mut payoff: impl FnMut(&PathState) -> f64,
) -> Estimate {
    let mut mc = AntitheticMonteCarlo::new(rng, grid.increments());
    //One scratch state, re-filled on each of the two mirrored evaluations, so the hot loop does
    //not allocate per path.
    let mut state = PathState::new(grid);
    for _ in 0..pairs {
        mc.push_pair(|normals| {
            walk(model, grid, start_rate, normals, &mut state);
            payoff(&state)
        });
    }
    mc.finish()
}

/// What a coarse/fine pair says about the discretisation bias of the fine run.
#[derive(Debug, Clone, Copy)]
pub struct Discretisation {
    /// `dt_coarse / dt_fine`.
    pub ratio: f64,
    /// `fine_mean - coarse_mean`: the observed refinement response.
    pub gap: f64,
    /// Independent-run noise in [`Discretisation::gap`].
    pub gap_standard_error: f64,
    /// Richardson residual bias of the fine run: `max(|gap|, gap standard error) / (ratio - 1)`
    /// under a first-order (`O(dt)`) bias.  This is the number that goes into the assertion's
    /// budget.
    ///
    /// The floor at the gap's own standard error is deliberate.  `|gap| / (ratio - 1)` is a
    /// *point estimate* of the residual and it is noisy: when the true bias is small relative to
    /// the noise in the difference, the point estimate can come out anywhere near zero and would
    /// license a bias allowance far tighter than the run is entitled to claim.  Flooring at
    /// `gap_standard_error` means the admitted bias is never smaller than what the difference
    /// could have hidden, which is the honest upper bound of the Richardson estimate at this
    /// sample size.  It widens the budget by at most one gap-standard-error over the plain
    /// Richardson number, and only when the plain number was already below the noise.
    pub residual: f64,
    /// Richardson-extrapolated, bias-free value: `fine + gap / (ratio - 1)`.
    pub extrapolated: f64,
}

impl Discretisation {
    /// Was the refinement response actually resolved, or is it inside the noise of the difference?
    ///
    /// `false` does not mean "no bias"; it means this pair could not see it, and the honest
    /// statement about the bias is then the noise floor rather than the Richardson number.
    pub fn resolvable(&self) -> bool {
        self.gap.abs() > self.gap_standard_error
    }
}

/// Estimate the residual discretisation bias of `fine` from `coarse`, assuming first-order Euler
/// convergence (see the bias section of the module docs).
pub fn discretisation_bias(coarse: &Estimate, fine: &Estimate, ratio: f64) -> Discretisation {
    assert!(
        ratio > 1.0 + 1e-12,
        "a bias study needs the coarse grid to be coarser than the fine one, got ratio {ratio}"
    );
    let gap = fine.mean - coarse.mean;
    let gap_standard_error = coarse
        .standard_error
        .hypot(fine.standard_error)
        .max(f64::EPSILON);
    let shift = gap / (ratio - 1.0);
    Discretisation {
        ratio,
        gap,
        gap_standard_error,
        residual: shift.abs().max(gap_standard_error / (ratio - 1.0)),
        extrapolated: fine.mean + shift,
    }
}

/// Everything a Monte-Carlo check is asserted on, assembled so the failure message carries the
/// whole budget rather than just the two numbers that failed to match.
#[derive(Debug)]
pub struct Check<'a> {
    /// What is being checked, e.g. `"caplet(T=1.5, K=0.02, antithetic, dt=1/400)"`.
    pub label: &'a str,
    /// The Monte-Carlo estimate.
    pub estimate: Estimate,
    /// The closed form it is checked against.
    pub analytic: f64,
    /// Deviation multiplier applied to [`Estimate::standard_error`].
    pub k: f64,
    /// Residual discretisation bias admitted into the budget ([`NO_DISCRETISATION`] where the
    /// estimator does not step through time).
    pub bias: f64,
    /// Optional: what the two-grid study said, for the report.
    pub discretisation: Option<Discretisation>,
}

impl<'a> Check<'a> {
    /// A check with no bias term: the estimator does not step through time.
    pub fn undiscretised(label: &'a str, estimate: Estimate, analytic: f64) -> Self {
        Self {
            label,
            estimate,
            analytic,
            k: K_SIGMA,
            bias: NO_DISCRETISATION,
            discretisation: None,
        }
    }

    /// A check against the bias study of a coarse/fine pair: the residual is the budget's bias
    /// term and the study itself is carried into the report.
    pub fn with_bias_study(
        label: &'a str,
        estimate: Estimate,
        analytic: f64,
        discretisation: Discretisation,
    ) -> Self {
        Self {
            label,
            estimate,
            analytic,
            k: K_SIGMA,
            bias: discretisation.residual,
            discretisation: Some(discretisation),
        }
    }
}

impl Check<'_> {
    /// `k * standard_error + bias`.
    pub fn budget(&self) -> f64 {
        self.estimate.band(self.k) + self.bias
    }

    /// `|mc - analytic|`.
    pub fn deviation(&self) -> f64 {
        (self.estimate.mean - self.analytic).abs()
    }

    /// Assert `|mc - analytic| <= budget`.
    ///
    /// The report is printed first, so a passing run under `--nocapture` and a failing run print
    /// the same text: the numbers behind a pass are recoverable from CI logs, and a failure
    /// arrives with the whole budget already spelled out.
    pub fn assert_within_budget(&self) {
        let deviation = self.deviation();
        let budget = self.budget();
        println!("{}", self.report());
        assert!(
            deviation <= budget,
            "{}\nbudget exceeded: |{:.3e} - {:.3e}| = {:.3e} > {:.3e}  (k = {} x SE {:.3e} = {:.3e}, \
             + discretisation residual {:.3e})\nEither the estimator moved by more than {} standard \
             errors can explain, or the discretisation residual above now understates the bias: refine \
             the grid before touching k.",
            self.label,
            self.estimate.mean,
            self.analytic,
            deviation,
            budget,
            self.k,
            self.estimate.standard_error,
            self.estimate.band(self.k),
            self.bias,
            self.k,
        );
    }

    /// Human-readable breakdown of the run and of how much of its budget it used.
    pub fn report(&self) -> String {
        let est = &self.estimate;
        let mut out = format!(
            "[mc] {label}\n\
             \x20     mc            = {mean:.12}\n\
             \x20     analytic      = {analytic:.12}\n\
             \x20     |difference|  = {deviation:.3e}\n\
             \x20     pairs         = {pairs} (paths = {paths}, antithetic)\n\
             \x20     SE            = {se:.3e}  (unpaired: {naive:.3e}, variance reduction {vr:.2}x)\n\
             \x20     k * SE        = {band:.3e}  with k = {k}\n\
             \x20     bias residual = {bias:.3e}\n\
             \x20     budget        = {budget:.3e}, used {used:.1}% of it",
            label = self.label,
            mean = est.mean,
            analytic = self.analytic,
            deviation = self.deviation(),
            pairs = est.pairs,
            paths = est.paths,
            se = est.standard_error,
            naive = est.naive_standard_error,
            vr = est.variance_reduction,
            k = self.k,
            band = est.band(self.k),
            bias = self.bias,
            budget = self.budget(),
            used = if self.budget() > 0.0 {
                100.0 * self.deviation() / self.budget()
            } else {
                0.0
            },
        );
        if let Some(d) = &self.discretisation {
            out.push_str(&format!(
                "\n\x20     bias study    = ratio {} | gap {gap:+.3e} (gap SE {gap_se:.3e}, {resolved})\n\
                 \x20                   Richardson extrapolation {extrapolated:.12}",
                if d.ratio == 2.0 {
                    "2 (coarse = 2x dt)".to_string()
                } else {
                    format!("{:.3}", d.ratio)
                },
                gap = d.gap,
                gap_se = d.gap_standard_error,
                resolved = if d.resolvable() {
                    "resolved"
                } else {
                    "below the noise floor: the bound is loose"
                },
                extrapolated = d.extrapolated,
            ));
        }
        out
    }
}

#[cfg(test)]
mod tests;
