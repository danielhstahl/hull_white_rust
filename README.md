| [CI (build + test)][ci-link] | [Lint (fmt + clippy)][lint-link] | [Coverage (Coveralls)][cov-link] |
| :-------------------------: | :---------------------------: | :-----------------------------: |
| ![ci-badge]                 | ![lint-badge]                 | ![cov-badge]                    |

<!-- Both badges name the default branch (`master`) explicitly rather than relying on
     "whatever the default is".  The coverage badge is Coveralls, not Codecov:
     `test.yml` runs `cargo llvm-cov` and uploads `lcov.info` to Coveralls, and
     nothing in this repository ever uploaded to Codecov, so the Codecov badge had
     no producer.  The CI badge used the legacy `/workflows/Rust/badge.svg` form and
     no workflow here is named "Rust" (they are `rusttest`, `lint` and `bench`), so
     it was broken; `/actions/workflows/<file>/badge.svg` keys off the file that
     actually exists. -->
[ci-link]:   https://github.com/danielhstahl/hull_white_rust/actions/workflows/test.yml?query=branch%3Amaster
[ci-badge]:  https://github.com/danielhstahl/hull_white_rust/actions/workflows/test.yml/badge.svg?branch=master
[lint-link]: https://github.com/danielhstahl/hull_white_rust/actions/workflows/lint.yml?query=branch%3Amaster
[lint-badge]: https://github.com/danielhstahl/hull_white_rust/actions/workflows/lint.yml/badge.svg?branch=master
[cov-link]:  https://coveralls.io/github/danielhstahl/hull_white_rust?branch=master
[cov-badge]: https://coveralls.io/repos/github/danielhstahl/hull_white_rust/badge.svg?branch=master

## Hull White

This library implements functions that price fixed income products assuming that short rates follow a Hull-White process.

## Documentation

The [Documentation](./documentation) directory holds the model documentation for the various pricing
functions and assumptions: R Sweave / LaTeX source plus its rendered PDF.  It lives in the
repository only — `package.exclude` in `Cargo.toml` keeps it out of the published crate, because the
rendered PDF is most of the tarball and nothing on the Rust side reads it.

Library (API) documentation is at
[docs.rs](https://docs.rs/hull_white/0.9.0/hull_white/), built from `cargo doc` on publish.  To
read it locally:

```bash
cargo doc --open --no-deps
```

Every documented example in the crate is a doctest, so `cargo test` compiles and runs the API docs
along with the unit tests — prose that stops matching the code fails the build rather than drifting.

## Requirements

The documentation is written in [R Sweave](https://www.r-bloggers.com/getting-started-with-sweave-r-latex-eclipse-statet-texlipse/).  The application is written in [Rust](https://www.rust-lang.org/en-US/).

## Install

Add the following package to your Cargo.toml:

`hull_white = "0.9.0"`

## Doc and version bookkeeping

Three things are kept honest at compile time rather than in review:

* **`#![deny(missing_docs)]`**, declared as `[lints.rust] missing_docs = "deny"` in `Cargo.toml` so
  it covers every target cargo builds through this manifest.  Anything `pub` and reachable from the
  crate root without a doc comment is a build error.  `pub(crate)` items are exempt — they are not
  the API.
* **The MSRV** is `rust-version = "1.85"` in `Cargo.toml`: the `edition = "2024"` floor.  CI builds
  `stable` and `nightly`; raise the field only when the code actually needs more, and say so in
  `CHANGELOG.md`.
* **`html_root_url` is deliberately absent.**  `#![html_root_url]` is no longer a recognised rustc
  attribute — adding it fails the build with "cannot find attribute" — and the doc root is supplied
  by whatever hosts the built docs (`docs.rs` sets it itself), with inter-crate links resolved via
  rustdoc's `--extern-html-root-url` at build time.  There is nothing version-pinned to maintain
  here, which is the point: a hardcoded one in the sources is exactly the kind of number this crate
  used to let go stale.

The version string itself is read by cargo from one place, `package.version` in `Cargo.toml`, and
written by hand in two: the install snippet above, and the versioned docs.rs link above.  Those two
are the ones that rot quietly — they sat on `0.6.0` long after the crate left it — so
[`tests/release_sync.rs`](tests/release_sync.rs) reads `Cargo.toml` and asserts that the README
snippet matches `package.version`, that every `docs.rs/hull_white/<version>/` link in the README
names `package.version`, that `CHANGELOG.md` carries a section for the released version **and** an
`## [Unreleased]` section for the next one, and that the doc comments still describe what the code
does (`*_now` pricers are not summarised as pricing at a future date; a `receiver` function is not
summarised as the payer side; the misspelling that used to appear twice in `src/model.rs` is gone).
It runs under plain `cargo test`, so CI fails on a release that bumped the crate and forgot the
docs.  The per-release sequence is:

1. `CHANGELOG.md`: rename `## [Unreleased]` to `## <new version> - <date>`, add a fresh empty
   `## [Unreleased]` above it.
2. `Cargo.toml`: bump `package.version`.
3. `README.md`: the install snippet and the docs.rs link.  `cargo test` names each of these that is
   still wrong, so step 3 is the test going green rather than a memory exercise.

## The curve

A model is calibrated to **one** curve object, not to a pair of closures.  `YieldCurve`'s one
required method is the *cumulative* yield `y(t)` — the integral of the instantaneous forward,
so `P(0,t) = exp(-y(t))` — and the forward `f(0,t) = y'(t)` is derived from it:

```rust
use hull_white::{HullWhite, from_yield};

let curve = from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
let model = HullWhite::new(0.15, 0.02, &curve).unwrap();
let price = model.bond_price_now(2.0).unwrap();   // exp(-y(2))
```

Where the derivative is known in closed form, supply it and skip the finite difference; where it
is not, `from_yield` is the whole job.  A `struct` implementing the trait works just as well, and
is the right shape for a bootstrapped or piecewise curve.

Because the yield and the forward are one object, the relation between them is checkable, and the
constructor checks it: `HullWhite::new` compares `forward(t)` with the numerical derivative of
`zero_yield(t)` at probe times out to ten years and refuses a curve that disagrees by more than a
hundredth of a basis point (`1e-6`, relaxing to 1 part in 10,000 of the forward).  A model
calibrated to a yield curve and a *different* forward curve misprices everything quietly; now it
is an error at construction, naming the time, both values and the tolerance.  See
[`src/curves.rs`](src/curves.rs) and the [0.9.0 entry in the CHANGELOG](CHANGELOG.md).

Migrating from 0.8: `HullWhite<'a, T, U>` is `HullWhite<'a>`, and
`HullWhite::init(a, sigma, &yield_curve, &forward_curve)` is deprecated in 0.9 in favour of
`HullWhite::new(a, sigma, &from_yield_and_forward(yield_curve, forward_curve))` — the same two
closures, wrapped as one curve, with the consistency check the old call never ran. `init` stays
available (deprecated) through `0.9.x` and is removed in `0.10.0`; the full migration table, and
what to do if the check rejects a pair that used to "work", are in
[CHANGELOG.md](CHANGELOG.md#090--one-curve-object-breaking).

## Benchmarks

Benchmarks are at https://danielhstahl.github.io/hull_white_rust/dev/bench/.
## Tests

```bash
cargo test                     # the whole suite; the Monte-Carlo checks run at the light tier (~0.2s)
cargo test --features slow     # same checks at the heavy tier: ~20x the paths, ~4x finer grid (~18s)
cargo test --features slow -- --nocapture   # print the per-run Monte-Carlo error-budget report
```

The Monte-Carlo checks (caplet, floorlet, Eurodollar future) do not assert against a hard-coded
`epsilon`.  Each one asserts

```text
| mc - analytic |  <=  4 * standard_error  +  estimated discretisation residual
```

where the standard error is measured from the run itself (antithetic variates, so it is the spread
of the pair averages) and the residual is the Richardson estimate from running the same estimator
on two grids, `dt` and `dt / 2`.  The sample size and grid resolution are tiered on the `slow`
feature so the tight end of that curve is opt-in: a plain `cargo test` keeps a ~4s signal, and CI
runs `--features slow` as well (the nightly coverage run uses `--all-features`, which includes it).
See `src/mc.rs` for the harness, the step-size convention, and the measured numbers.

## CI

Four workflow files under [`.github/workflows`](.github/workflows):

| File | Workflow | Triggers | Runs |
| --- | --- | --- | --- |
| [`test.yml`](.github/workflows/test.yml) | `rusttest` | every push and PR, `stable` + `nightly` | build, `cargo test`, heavy MC tier (`--features slow`), nightly `cargo llvm-cov` → Coveralls, `cargo doc` |
| [`lint.yml`](.github/workflows/lint.yml) | `lint` | pushes to `master`/`main` and every PR | `cargo fmt --all -- --check`; `cargo clippy ... -D warnings` — on stable over lib + bins + tests, on nightly over `--all-targets`, which adds the benches |
| [`benchmark_ghpages.yml`](.github/workflows/benchmark_ghpages.yml) | `bench` | **default branch only** (+ `workflow_dispatch`) | `cargo bench` → the published trend chart on `gh-pages` |
| [`rust.yml`](.github/workflows/rust.yml) | `RustDeploy` | push to `master`/`main` | build, test, doc, `cargo publish` |

**Required checks.** A job that runs is not a job that blocks. In
Settings → Branches → *Require a pull request before merging* → *Required checks*, add
`fmt`, `clippy (stable)` and `clippy --all-targets (nightly)`. The lint workflow
produces a status worth requiring: green on a clean tree, red on any diff — no
`continue-on-error`, nothing allowed to be yellow.

**Why clippy runs twice.** `--all-targets` includes `benches/`, and those are
libtest `#[bench]` functions needing `#![feature(test)]` — nightly only, plus
`--features test-support` for the shared fixtures. So stable lints the crate that
ships and nightly lints everything including the benches; neither is a subset of
the other in what it can catch, because the lints themselves differ per toolchain.

**`cargo bench` runs on the default branch only.** The `fail-on-alert` comparison
is against the last *committed* baseline: on a feature branch it compares against
unrelated code, so the alert is noise, and `comment-on-alert` posts that noise on
the PR. The `gh-pages` variant force-pushes to a shared published branch, which is
not something an unreviewed branch should get to do with `contents: write`. Both
gate on `github.ref == format('refs/heads/{0}', github.event.repository.default_branch)`
rather than a hardcoded `master`, so the gate follows the default branch if it is
ever renamed instead of silently disabling.

**Running the jobs before pushing.** `bash .github/scripts/ci_local.sh` runs each
job's command locally and prints one PASS/FAIL/SKIP row per job — the branch-level
dry run, without needing a push. `--list` prints the job names, `--bench` actually
runs the benches (compile-only otherwise, since timings from shared hardware are
not the number the trend wants). A missing tool is a SKIP with a reason, never a
PASS. The deploy path is dry-runnable on its own:
`bash .github/scripts/publish.sh --dry-run` runs every check the real publish runs
and stops at `cargo publish --dry-run` — nothing is uploaded and no token needed.

**Action pins.** Every third-party action is pinned to a full commit SHA with the
version in a trailing comment. A tag like `@v4` is mutable — the owner can move it
under this repository — and `@master` is worse: each run takes whatever whoever
pushed last. That is a supply-chain exposure, so the pins are recorded here and are
meant to move deliberately, in a commit that says why:

| Action | Version | Commit |
| --- | --- | --- |
| `actions/checkout` | v4.4.0 | `11d5960a326750d5838078e36cf38b85af677262` |
| `hecrj/setup-rust-action` | v2.0.1 | `110f36749599534ca96628b82f52ae67e5d95a3c` |
| `taiki-e/install-action` | v2.87.21 | `4cef1412cce204788f482e778a0b9187f9626a29` |
| `actions/cache` | v5.1.0 | `caa296126883cff596d87d8935842f9db880ef25` |
| `benchmark-action/github-action-benchmark` | v1.9.0 | `fd31771ce86cc65eab85653da103f71ab1b4479c` |
| `coverallsapp/github-action` | v2.3.8 | `8d6379e14d29928660c4ba802d8e85393440b329` |

Each SHA was resolved from the tag it names (`git ls-remote` and the refs API) at
the time of pinning. To bump one, change the SHA and the version comment together.

**Publish failures are loud.** `cargo publish` used to end in `|| true`, which made
every failure invisible: expired token, rejected manifest, version already
published — the deploy job reported success either way. It now runs under
`set -euo pipefail`, and the single tolerated case is spelled out in the workflow:
the version in `Cargo.toml` is already on crates.io, meaning the push did not bump
`package.version`. That path emits a `::warning::` annotation on the run, so
"skipped, not shipped" is a visible state rather than the default one.
