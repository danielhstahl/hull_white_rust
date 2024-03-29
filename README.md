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
