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
minor bump rather than the **0.9.1** patch first sketched here. No price moves that
any caller can see: every number this crate returns is the number it was in 0.9.0.
One entry — the swap-annuity extraction — reassociates a sum, which moves prices in
their last bits (worst measured `|Δ| = 4.4e-16`, worst relative `3.7e-11` on a price
sitting within a rounding of its own zero, no sign flips). That is recorded by the tree
pins rather than smoothed over: `trees::tests::swaption_tree_at_t_is_bit_identical_to_pinned_values`
pins the new bit patterns and
`trees::tests::swaption_tree_at_t_is_within_its_documented_reassociation_slack_of_0_9_0`
bounds the distance from the 0.9.0 values at 16 ulps.

### Changed — identifiers that said the opposite of what they meant

Naming and plumbing; no price moves (checked by diffing a grid of ~1,850 prices across all
seven scenarios against the pre-change build — every number identical, the only differences
being the two error messages named below).

* `HullWhite::t_forward_bond_vol(t, t_m, t_f)` now takes `option_maturity` and
  `bond_maturity`. The old names were not merely vague: the doc described them as "the bond
  that matures at `t_m` is the asset, the option expires at `t_f`", which is the reverse of
  what the formula and every call site in the crate mean — `t_m` *was* the option maturity,
  `t_f` the bond maturity, and the function's own `t_f > t_m` check contradicts the sentence
  above it.  The doc now carries the Hull-White expression with each symbol defined (`T` =
  `option_maturity`, `T_b` = `bond_maturity`), and
  `model::tests::t_forward_bond_vol_is_the_integrated_forward_measure_variance` pins it by
  integrating the deliverable's forward-measure volatility
  `sigma * (B(s, T_b) - B(s, T))` over the option's life: the two dates are interchangeable
  in neither, so the roles are now checked in numbers rather than in prose.
* Invalid-input messages follow the rename — `t_forward_bond_vol(1.0, 2.0, nan)` now reports
  `bond_maturity = NaN is not finite` where it reported `t_f`, and the option leg reports
  `option_maturity`.  Callers matching on those strings are the only thing this can break.
* The affine bond maths is two crate-internal methods, `HullWhite::bond_b` (`B(t, T)`) and
  `HullWhite::bond_c` (`C(t, T)`), replacing three free functions that threaded the model's
  parameters by hand: `a_t(a, t_diff)`, `ct_t(a, sigma, t, t_m, curve)` and `at_t(a, t,
  t_m)` — the last a pure alias of the first with `t_diff = t_m - t` written into its own
  argument list.  Two names now, each documenting which half of `P(t, T) = exp(C - B * r)` it
  is, and `model::tests::the_affine_coefficients_rebuild_the_bond_price_and_its_duration`
  shows the reassembly is bit-identical to `bond_price_t` and that `B` really is that price's
  duration (`-dP/dr / P`, measured by finite difference).
* `gamma_edf` is a model method — `self.gamma_edf(t, option_maturity, delta)` — with its
  formula and its source links attached, and lives in `rates` with the instrument that uses
  it.  `edf_compute` and `compute_libor_rate` moved there too, and deliberately **stay** free
  functions: they read no model state at all (`a`, `sigma` and the curve are all spent
  upstream, into the bond prices and into `gamma`), so binding `&self` would advertise a
  dependency that is not in the arithmetic.
* The tree's diffusion callbacks are named for what they map.  What was `sigma_inv`, a closure
  whose argument is the lattice's pure-Brownian coordinate `w` and which *returns* the
  shifted state (`w -> sigma * w`), is `state_from_lattice_coord`; what was `sigma_prime`,
  i.e. `d sigma / d(state)`, is `d_sigma_d_state` (identically zero, a constant
  volatility).  The old names had the inverse-volatility callback's argument written as `y`,
  reading as though the engine handed over a rate that was then multiplied by sigma, when it
  is the `w -> y` map itself.
* `Sync` is spelled `Sync` in `HullWhite::init`'s closure bounds rather than
  `std::marker::Sync`.
* Crate docs: the "Key Concepts" `(0, t, T, TM)` block — itself the shape of the ambiguity
  above, since `T` and `TM` never say which maturity is whose — is replaced by the list of
  role-named time parameters (`t`, `option_maturity`, `bond_maturity`, `delta`) and the rule
  that textbook symbols appear only inside formulas, spelled against the name they stand for.
  The crate-layout table is updated for the moved items.
* `options.rs`: four `//volatility with maturity` end-of-line comments — which described
  neither the volatility nor the maturity — now say the volatility of the deliverable bond.

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

### Changed — the swap payment convention is a doc comment on one function

`forward_swap_rate_t` and `swap_price_t_init` each hand-rolled the same discounted
payment sum, over slightly different ranges, and `swap_price_t_init_raw` carried a
live `//open question, should num_swap_payments be num_swap_payments+1??` inside
pricing code that swaps, European swaptions and the American tree all depend on. The
sums are now one helper, and the convention it computes is written down on it.

* `HullWhite::annuity_t(r_t, t, start, n, delta)` is
  `delta * sum_{i=1..n} P(t, start + i * delta)` (crate-internal, both swap
  consumers call it, both hand-rolled loops deleted). Its doc comment states the
  convention: a swap's `n` payment dates are `start + delta, ..., start + n * delta`;
  **every one of them carries a coupon** (`K * delta` fixed, the reset Libor times
  `delta` floating); and **the principal rides on the last of those same `n` dates**,
  on top of that date's coupon — so `swap_maturity == start + n * delta` and the
  principal exchange adds no date of its own. Worked example in the comment:
  `start = 1.0, delta = 0.25, n = 4` → 1.25 / 1.50 / 1.75 / **2.00**, with the
  notional on 2.00.
* The open question is answered **no**, with the algebra rather than a shrug. The old
  body summed coupons over `1..n` and paid a combined `(1 + K*delta)` at the `n`-th
  date, which is why the `n`-th coupon looked missing — it is inside that
  `(1 + K*delta)`, paid with the principal. Extending the sum to `n` and folding the
  principal back out is the same number term by term:
  `P_start − Σ_{i<n} K·δ·P_i − (1 + K·δ)·P_n = P_start − K·A − P_n`. Reading the
  count as `n + 1` would not add a payment, it would add a *period*: the worked
  example would pay a fifth coupon, and the principal, at 2.25 — a quarter past the
  maturity.
* Both consumers now read the same `A`, which is what makes the price's zero *be* the
  rate: `forward_swap_rate_t = (P(t,start) − P(t,start+n·δ)) / A` is the fixed-equals-
  floating condition, written for `K`. Setting the legs equal and solving gives back
  exactly the formula the function computes — the two are counting the same dates.
* Numbers. The refactor changes the association of the fixed-leg sum (fold `P_i` over
  `1..=n`, scale the sum once by `K`, carry principal as its own term). Mathematically
  identical; in f64 it is a reordering. Measured over 16,800 swap prices — all seven
  named scenarios × 4 valuation times × 3 start offsets × `n ∈ {1, 4, 8, 20}` ×
  `delta ∈ {0.25, 0.5}` × 5 strikes around the forward × 5 short rates:
  worst `|Δ| = 4.4e-16`, worst relative `3.7e-11` (only on a one-period swap priced
  within a rounding of its own zero, where cancellation does the amplifying), **zero**
  sign flips, and 314 near-zero prices became *exactly* `0.0` where the old grouping
  left a `~1e-16` residue. `annuity_t` keeps the old fold order of the discount
  factors, so `forward_swap_rate_t` itself is bit-for-bit unchanged.
* `swap_price_t`'s off-schedule anchor is documented as a **known gap, deliberately
  not fixed here**. When the maturity is a whole number of periods past `t` the derived
  anchor is `t`, which satisfies the convention; otherwise the branch uses
  `swap_maturity − (n − 1) * delta`, which is the *next payment date* (`schedules`'
  own tests name that expression `next_exchange_date`), not the reset date one period
  before it. That branch therefore drops the first remaining payment and repays the
  principal one period past the stated maturity (`t = 0.5`, maturity `2.0`, `delta`
  `0.4` → coupons at 1.2 / 1.6 / 2.0 / **2.4** instead of 0.8 / 1.2 / 1.6 / 2.0).
  The payment *count* is right; the anchor is a period late. Correcting it moves
  prices by a whole period, not by rounding, so it is left for its own change instead
  of riding along inside a refactor whose other effects end in ulps.
* Tests. `swaps::tests::forward_swap_rate_is_the_exact_zero_of_the_swap_price_on_the_
  semiannual_schedule` — a second schedule, 8 payments at `delta = 0.5`, every named
  scenario × `t ∈ {0, 0.25, 0.5, 1.0}`, asserted `== 0.0` exactly (not to a
  tolerance), plus `±1 bp` of strike moving the price by `∓annuity_t * bp`: the
  annuity showing up as the swap's DV01 is a second, independent check that it sums
  the right cash flows.
  `swaps::tests::annuity_t_is_the_discounted_weight_of_the_payment_dates_it_documents`
  checks the helper against hand-written discount factors for the doc's four dates,
  against the one-payment case, against the `n + 1` additivity step, and against the
  leg-balance identity `K * A + P(maturity) == P(start)`.
  `swaps::tests::one_more_payment_is_a_different_swap_not_the_principal_being_counted`
  pins the resolved question: the `n + 1` misreading is a 7.0e-4 instrument change on
  that fixture — twelve orders of magnitude above the rounding this refactor moved —
  and is flat only at the 5-payment forward rate. Existing `test_swap` /
  `test_swap_init` pass unchanged.
* The `t = 0` tree goldens are re-pinned for the reassociation (they moved 4–10 ulps).
  `swaption_tree_at_t_is_bit_identical_to_pre_fix` is now
  `swaption_tree_at_t_is_bit_identical_to_pinned_values`, and the new
  `swaption_tree_at_t_is_within_its_documented_reassociation_slack_of_0_9_0` keeps
  the 0.9.0 numbers in the file and asserts the 16-ulp bound, so the pin stays exact
  against future changes while the size of this one stays checked rather than
  remembered.

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
