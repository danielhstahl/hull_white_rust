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

`hull_white = "0.6.0"`

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
