| [Linux][lin-link] |  [Codecov][cov-link]  |
| :---------------: | :-------------------: |
| ![lin-badge]      | ![cov-badge]          |

[lin-badge]: https://github.com/danielhstahl/hull_white_rust/workflows/Rust/badge.svg
[lin-link]:  https://github.com/danielhstahl/hull_white_rust/actions
[cov-badge]: https://codecov.io/gh/danielhstahl/hull_white_rust/branch/master/graph/badge.svg
[cov-link]:  https://codecov.io/gh/danielhstahl/hull_white_rust

## Hull White

This library implements functions that price fixed income products assuming that short rates follow a Hull-White process.

## Documentation

The [Documentation](./documentation) holds the model documentation for the various pricing functions and assumptions.  Library (API) documentation is available at [docs.rs](https://docs.rs/hull_white/0.6.0/hull_white/)

## Requirements

The documentation is written in [R Sweave](https://www.r-bloggers.com/getting-started-with-sweave-r-latex-eclipse-statet-texlipse/).  The application is written in [Rust](https://www.rust-lang.org/en-US/).

## Install

Add the following package to your Cargo.toml:

`hull_white = "0.9.0"`

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
