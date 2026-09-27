# Changelog

All notable changes to this crate are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html). The per-release
procedure — what gets renamed, bumped and re-linked, in order — is in
[README.md](README.md#doc-and-version-bookkeeping), and
[`tests/release_sync.rs`](tests/release_sync.rs) fails the build if a step is
skipped.

## [Unreleased]

Intended next release: **0.10.0**. The block below adds public methods, which is a
minor bump rather than the **0.9.1** patch first sketched here. No price moves: every
number this crate returns is the number it was in 0.9.0, and the tree results are
pinned bit-for-bit (`trees::tests::swaption_tree_at_t_is_bit_identical_to_pre_fix`), so
the extra surface is not carrying a numeric change in disguise.

### Added — the European swaption tree, public, so it can be cross-checked

* `HullWhite::european_payer_swaption_tree`, `HullWhite::european_receiver_swaption_tree`,
  and the `..._tree_now` variants. Same instrument and same argument list as the
  closed-form `european_*_swaption_t` / `european_*_swaption_now` pricers, plus
  `num_steps` for the tree resolution.
* The tree already existed — as a `#[cfg(test)]` private function. So the parity check it
  was written for could not be run by anyone outside this crate, even though it is the one
  check of the analytic Jamshidian swaption that needs no third pricer: a lattice and a
  decomposition reach the same price by independent means, and a wiring fault in either one
  shows up as a gap. The time-coordinate bug this module carried was a 13–16% error at
  `t > 0` that neither route could see on its own.
* They are documented as **verification helpers**, not as a way to price a European
  swaption. `european_*_swaption_t` is exact within the model's own approximation and
  costs a tree nothing; reach for the tree to check the closed form, not to replace it.
* `HullWhite::tree_price` (`pub(crate)`): the lattice engine — drift, volatility, the
  per-step `phi` cache, the node discount, the exercise style — taking the payoff as a
  closure written against the short rate at the node.

### Changed — one copy of the tree wiring

* `american_swaption` and `european_swaption_tree` were the same ~35 lines, hand-maintained
  twice: identical `alpha_div_sigma` / `sigma_prime` / `sigma_inv`, identical `phi_cache`
  construction including its extra tail element, identical payoff and discount closures,
  differing only in the final call. They are now one `SwaptionSpec` (state, swap leg, side,
  resolution) handed to one engine with an `is_american` flag.
* A third near-copy of the same wiring, `tree_option_price` in
  `jamshidian/tests/straddling.rs`, runs on that engine too. The tree-time
  (`tau`) versus absolute-model-time (`t + tau`) convention — the thing that had the bug —
  is now written down once, in `HullWhite::tree_price`, instead of three times in three
  files that had to be edited together.
* `binomial_tree::compute_price_american` is no longer called. Both sides go through
  `compute_price_raw(..., is_american)`; upstream the former *is* the latter with `true`, so
  the American numbers are unchanged to the bit — which is why the existing golden test
  still passes with its exact bit patterns rather than a tolerance.
* The four `american_*` public signatures are unchanged.
* The four copies of the swaption argument-validation block are one `validate_swaption`, so
  a bad instrument is refused with the same message by every entry point on the tree, on
  either side and in either exercise style.
* Tests: the tree tests now price through the **public** entry points rather than the private
  helper, so what they cover is what a consumer calls. Added
  `european_tree_matches_analytic_across_every_scenario_and_valuation_time` (all seven named
  calibrations, `t` = 0 / 0.5 / 1.0, both sides, 400 steps, relative bound 2.5e-3 against a
  measured worst case of 7.0e-4),
  `every_tree_entry_point_rejects_the_same_bad_instruments` (the shared validator), and
  `tests/public_api.rs::european_tree_prices_cross_check_the_analytic_swaptions`, which
  runs the cross-check from outside the crate — the thing the `cfg(test)` helper made
  impossible.

### Fixed — documentation that contradicted the code

Copy/paste drift in doc comments: the sentences were present, and wrong, so nothing
failed. Each is now stated as what the function returns.

* `mu_r` carried `t_forward_bond_vol`'s summary line ("Returns volality of bond
  under the t-forward measure"). It returns the **conditional mean of the short
  rate**, `E[r(t_m) | r(t) = r_t]`; the doc now gives the mean-reversion formula,
  names `phi`, and cross-links the two things a reader might actually have been
  after (`variance_r`, `t_forward_bond_vol`).
* The misspelling "volality" (twice, both in that shared sentence) is gone;
  `grep -rn volality src` is now empty and
  `tests/release_sync.rs::no_volality_typos_in_src` keeps it empty.
* `euro_dollar_future_now` promised a price "at some future time" on a function
  whose whole point is that the valuation is today and the state comes from the
  calibration. Restated, with the two-bond-bracket construction of the futures
  rate and an assertion that it equals the `t = 0` case of `euro_dollar_future_t`.
* `european_receiver_swaption_t` opened with "Returns price of a **payer**
  swaption". It is the receiver side: a call, struck at par, on the deliverable
  coupon bond. The summary now says receiver, and the example asserts the
  payer/receiver spread against the forward swap rate so the sign is not a matter
  of trusting the prose.
* `coupon_bond_price_now`'s `coupon_times` parameter comment said the schedule
  "does not include the bond_maturity, but the function does check for that" — it
  does neither. The kernel adds the par value to the **last** element of the
  slice, so the maturity is *in* the schedule: `coupon_times[len - 1]` is the
  maturity and its payment is `coupon_rate + 1.0`. Both the doc comment and the
  inline parameter comment now say so, and the doctest pins the convention three
  ways (coupons-plus-principal decomposition, a one-element schedule pricing as a
  zero coupon bond, and a repeated maturity date being refused rather than paid
  twice).
* `get_coupon_times` had no doc comment at all.

### Added — enforcement, so this cannot quietly drift again

* `[lints.rust] missing_docs = "deny"` in `Cargo.toml`. Declared in the manifest
  rather than as `#![deny(missing_docs)]` in `src/lib.rs` so it covers every target
  cargo builds through that manifest (lib, doctests, benches), not just the one
  file carrying the attribute. `pub(crate)` internals are exempt — they are not the
  API.
* `tests/release_sync.rs`: README install snippet and every versioned `docs.rs`
  link checked against `package.version`; CHANGELOG required to carry a section for
  the released version **and** an empty `## [Unreleased]` for the next one; and
  mechanical doc-shape rules over `src/` for the patterns that caused the above —
  a `*_now` pricer summarised as pricing at a future date, a swaption summary
  naming the opposite side of the trade, the "volality" misspelling, and public
  functions whose summary line is too short to say what is being priced. Runs
  under plain `cargo test`.
* `Cargo.toml`: `keywords`, `categories`, and `rust-version = "1.85"` (the
  `edition = "2024"` floor — nothing in the crate needs more today).
* `package.exclude` for `documentation/` (R Sweave / LaTeX source and the
  rendered PDF): 188 KiB, most of the published tarball, read by nothing on the
  Rust side. The directory stays in the repository and the README says where.

### Changed — CI: pinned actions, publish failures that actually fail, a fmt/clippy gate

CI ran on floating refs and swallowed the one failure that mattered.

* Every third-party action is pinned to a full commit SHA with the version in a
  trailing comment (`actions/checkout` v4.4.0, `hecrj/setup-rust-action` v2.0.1,
  `taiki-e/install-action` v2.87.21, `actions/cache` v5.1.0,
  `benchmark-action/github-action-benchmark` v1.9.0,
  `coverallsapp/github-action` v2.3.8). `setup-rust-action` was on `@master` and
  the benchmark workflow was on `actions/checkout@v2` while the others were on `v4`
  — floating refs, so a stranger's push moved this repository's CI without a PR
  here. The pin table is in [README.md](README.md#ci).
* `cargo publish --token ... --allow-dirty || true` → `set -euo pipefail`, and the
  one tolerated case narrowed and stated in the workflow: the version in
  `Cargo.toml` is already on crates.io (a push that did not bump
  `package.version`), which emits a `::warning::` annotation rather than passing
  silently. An empty `CARGO_TOKEN` now fails explicitly; the tracked tree is
  asserted clean before publish, so the retained `--allow-dirty` cannot be hiding
  a modified tracked file.
* New `lint.yml`: `cargo fmt --all -- --check` and `cargo clippy ... -D warnings`
  as PR-blocking jobs (`fmt`, `clippy (stable)`,
  `clippy --all-targets (nightly)` — the nightly run adds the bench target, which
  needs `#![feature(test)]`). There was no style/lint job at all before, so
  formatting and lint drift were unbounded. The warnings that existed when the gate
  went up are fixed rather than allowed away; the only `allow`s are four
  `clippy::too_many_arguments` on the swaption pricers, whose arity mirrors the
  swap leg and changing which is a breaking API change in its own right.
* `cargo bench` is restricted to the default branch in both bench jobs (gate is
  `github.event.repository.default_branch`, not a hardcoded `master`), so an
  unreviewed branch can no longer drive `contents: write` against the published
  `gh-pages` trend or post a baseline-mismatch alert on an unrelated PR.
* Badges: the CI badge pointed at `workflows/Rust/badge.svg` and no workflow is
  named `Rust`; the coverage badge pointed at Codecov and nothing uploads there
  (`test.yml` uploads `lcov.info` to Coveralls). Both replaced with badges for what
  exists, each naming the default branch (`master`) explicitly.
* `permissions: contents: read` added to the workflows that need nothing more.

### Notes

* `html_root_url` is **not** set, deliberately. `#![html_root_url]` is no longer a
  recognised rustc attribute — adding it fails the build with "cannot find
  attribute" — and the doc root is supplied by whoever hosts the built docs
  (`docs.rs` sets its own), with cross-crate links resolved by rustdoc's
  `--extern-html-root-url` at build time. There is no version-pinned URL left to
  rot, which was the point.

## 0.9.0 — one curve object (breaking)

### The decision

`HullWhite<'a, T, U>` took **two** closures — a cumulative `yield_curve` and an
`forward_curve` — that are one thing. The forward is derivable from the yield
(`f(0,t) = d/dt y(t)` under this crate's convention), so nothing stopped a caller
from passing a pair that disagreed, and a model built on a mutually inconsistent
pair prices every instrument wrong with no diagnostic. The generic pair also leaked
into user code: every function that took a model had to name two curve closure
types it did not care about.

The chosen shape is **a trait, `YieldCurve`, with the cumulative yield as its one
required method and the forward derived from it** — plus a construction-time
consistency check for anybody who supplies the forward explicitly:

```rust
pub trait YieldCurve: Sync {
    fn zero_yield(&self, t: f64) -> f64;                     // required: y(t)
    fn discount(&self, t: f64) -> f64 { (-self.zero_yield(t)).exp() }   // P(0,t)
    fn forward(&self, t: f64) -> f64 { /* 5-point central difference */ }
}
```

Why this shape rather than the alternatives:

* **One required method** because the cumulative yield is the actual primitive:
  `P(0,t) = exp(-y(t))` and `f(0,t) = y'(t)`. A caller with only a yield curve
  gets a working model with nothing else to write, and there is then no second
  number for it to disagree with.
* **`forward` overridable** because a finite difference is not always wanted.
  Every closed-form curve in this crate — the test fixtures, a Vasicek calibration,
  a bootstrapped piecewise curve — knows its own derivative, and overriding costs
  one method and takes the finite difference out of the pricing path.
* **A consistency check anyway**, rather than "derive it and there is nothing to
  check", because the case that matters is the one where two numbers are supplied:
  a caller migrating from `init`, or a curve implemented from a calibration where
  the yield and the forward were fitted separately. The check cannot help someone
  who derived the forward; it exists precisely for the people who did not.
* **Dynamic dispatch over generics.** `HullWhite` holds a `&dyn YieldCurve`, so the
  type parameter goes away entirely. The curve evaluations are already two indirections
  deep inside a Jamshidian solve that spends ~1,400 evaluations per option; the
  loss against a monomorphised closure is not measurable against that, and the type
  a consumer writes is `HullWhite<'a>`.

### The convention, stated

`zero_yield` is **cumulative**: `y(t) = integral of f(0,.) over [0, t]`, so
`P(0,t) = exp(-y(t))` with no division by maturity. The consistency relation under
that convention is the plain `f = y'`.

This is worth spelling out because the *annualised* zero rate `z(t) = y(t)/t` — what
a quoting screen prints — satisfies a different relation to the same forward,
`f(t) = z(t) + t z'(t)`. A `YieldCurve` implemented against the annualised zero
is not "off by a factor of t" in some distant price: it fails the construction
check, by name, at the first probe.

### The check

`HullWhite::new` compares `curve.forward(t)` against the 5-point finite-difference
derivative of `curve.zero_yield(t)` at every time in `CURVE_PROBE_TIMES`
(`0.25, 0.5, 1.0, 2.0, 5.0, 10.0` years), and rejects the curve if any probe
disagrees by more than

```
tolerance = max(1e-6, 1e-4 * |forward(t)|)
```

i.e. a hundredth of a basis point, relaxing to one part in 10,000 of the forward
where the forward is large. The failure is
`HullWhiteError::InvalidInput` carrying the worst probe time, the supplied forward,
the derived one, the error and the tolerance. A curve that cannot be evaluated at
any probe time is rejected too — "nothing was checked" is not a pass.

A supplied forward that is **not finite** over a finite yield is always rejected, and
the relative term is deliberately not applied to it: `1e-4 * |inf|` is `inf`, and
`inf > inf` is false, so growing the tolerance by the magnitude of the supplied
forward would certify `forward(t) = infinity` as consistent and then carry that
infinity into `phi(t)`, `mu_r` and the Jamshidian bracket. For a non-finite
`forward` the tolerance stays at the `1e-6` absolute floor, so the comparison is
`inf > 1e-6` — true — and the curve is refused. `NaN` was already caught by this
rule; `±infinity` was the case that slipped, and it is the nastier one because a
`NaN` price announces itself and an infinite drift does not.

`0.0` is deliberately not a probe time: the derivative there is one-sided and the
least accurate, and the front of the curve already has its own gate — a curve whose
forward is not finite at `0` cannot supply `r(0)` and is reported by
`HullWhite::short_rate_now`, which is a different fact about a different thing.

Measured against the crate's own fixtures the check has six orders of headroom
(worst `2.0e-12` on `high_vol`, tolerance `5.1e-6`), so it fires on curve-shape
errors and not on rounding.

### Breaking changes

| Before 0.9 | 0.9 |
|---|---|
| `HullWhite<'a, T, U>` | `HullWhite<'a>` |
| `HullWhite::init(a, sigma, &yield_curve, &forward_curve)` | `HullWhite::new(a, sigma, &curve)` |
| `fn f<T: Fn(f64)->f64 + Sync, U: Fn(f64)->f64 + Sync>(m: &HullWhite<'_, T, U>)` | `fn f(m: &HullWhite)` |
| two closures, trusted | one curve, checked |

`HullWhite::init` is `#[deprecated(since = "0.9.0")]`. It stays compiling and
deprecated for the whole of `0.9.x`, and is **removed in `0.10.0`** — the breaking
type collapse ships with `0.9.0`, the removal trails it by one release so a caller
can adopt the new constructor without the two changing underneath at once.
**The deprecated call is not behaviour-preserving**: it wraps its two closures in
`from_yield_and_forward` and runs the same consistency check, so a pair that used
to price (wrongly) is now an error. Migrate off it rather than silencing the
deprecation.

### Migration

One closure — the cumulative yield. The forward is derived:

```rust
// 0.8
let y = |t: f64| 0.05 * t + 0.01 * t * t;
let f = |t: f64| 0.05 + 0.02 * t;
let m = HullWhite::init(a, sigma, &y, &f)?;

// 0.9 -- you already have everything the model needs in `y`
let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
let m = HullWhite::new(a, sigma, &curve)?;
```

Both halves in closed form — the literal migration, no finite difference:

```rust
// 0.9
let curve = hull_white::from_yield_and_forward(y, f);
let m = HullWhite::new(a, sigma, &curve)?;
```

Your own curve — a bootstrapped discount curve, a piecewise-linear forward, a
calibration object. Implement the trait; override `forward` if the derivative is
known and `discount` if the curve is stored as prices rather than yields:

```rust
struct MyCurve { /* ... */ }
impl YieldCurve for MyCurve {
    fn zero_yield(&self, t: f64) -> f64 { /* ... */ }
    fn forward(&self, t: f64) -> f64 { /* analytic derivative */ }
}
```

`HullWhite::curve()` returns the curve back (`&dyn YieldCurve`), so
`model.curve().discount(t)` is available wherever a price came from it.

### If the migration check rejects a pair that used to work

It is telling you the pair was inconsistent, which is to say the prices were
already wrong and nothing had said so. Three ways out, in order of preference:

1. Supply the actual derivative of the yield you passed. This is the fix.
2. Drop the forward and use `from_yield(y)` — the model derives it, and there is
   then no pair to be inconsistent.
3. Implement `YieldCurve` on your own type and make the two accessors agree; the
   error message names the probe time and both values, so the divergence is
   findable.

There is deliberately no "skip the check" flag. A model that will not build on an
inconsistent curve is the feature.

### Also in 0.9.0

* `pub mod curves` is public: `YieldCurve`, `from_yield`, `from_yield_and_forward`,
  `ForwardInconsistency`, `max_forward_inconsistency`, `validate_curve`,
  `CURVE_PROBE_TIMES`, `FORWARD_CONSISTENCY_TOLERANCE`,
  `FORWARD_CONSISTENCY_RELATIVE_TOLERANCE`, `forward_consistency_tolerance`.
* `HullWhite` implements `Debug` (model parameters; the curve is a trait object
  and is named rather than dumped).
* `test_support::hw_curves` / `Scenario::curves`, which returned a
  `(yield_closure, forward_closure)` tuple, are replaced by `HwCurve` /
  `hw_curve` / `Scenario::curve` returning the single object. Test-only and
  `#[doc(hidden)]`, but benches building on `--features test-support` need the
  rename.
* The test-suite fixtures are unchanged in value: every price the suite asserts is
  the same number it was before the refactor, rewired through the new constructor.
