//! `cargo bench` targets behind the published benchmark trend
//! (<https://danielhstahl.github.io/hull_white_rust/dev/bench/>).
//!
//! Two things the rest of the crate does not need and this one does:
//!
//! * `#![feature(test)]` — these are libtest `#[bench]` functions, so nightly only.
//!   That is why CI lints this target in the nightly job and not on stable
//!   (see `.github/workflows/lint.yml`).
//! * `--features test-support` — `benches/` is compiled as its own crate and cannot
//!   see `cfg(test)`, so the shared curve fixtures come from
//!   `hull_white::test_support`: `cargo bench --features test-support`.
//!
//! Every bench reads its inputs from those fixtures, so what is timed is the same
//! curve the unit tests assert against, not a copy of it.  This is a test-only
//! target: a consumer's `cargo build` never compiles it, but `cargo bench` and the
//! nightly lint job both do.
#![feature(test)]
extern crate test;

use hull_white::test_support::{BASELINE, FLAT_5PCT};
use test::Bencher;

#[bench]
fn bench_bond_t(bench: &mut Bencher) {
    let fixture = BASELINE;
    let curr_rate = fixture.curr_rate;
    let future_time = 0.0;
    let maturity = 1.5;
    let curve = fixture.curve();
    let hull_white = hull_white::HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();
    bench.iter(|| {
        hull_white
            .bond_price_t(curr_rate, future_time, maturity)
            .unwrap()
    })
}

#[bench]
fn bench_bond_now(bench: &mut Bencher) {
    let fixture = BASELINE;
    let maturity = 1.5;
    let curve = fixture.curve();
    let hull_white = hull_white::HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();
    bench.iter(|| hull_white.bond_price_now(maturity).unwrap())
}

#[bench]
fn bench_coupon_bond_t(bench: &mut Bencher) {
    let fixture = BASELINE;
    let curr_rate = fixture.curr_rate;
    let delta = fixture.delta;
    let future_time = 0.0;
    let coupon_rate = 0.05 * delta;
    let curve = fixture.curve();
    let hull_white = hull_white::HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();

    bench.iter(|| {
        let coupon_times = hull_white::get_coupon_times(5, future_time, delta).unwrap();
        hull_white
            .coupon_bond_price_t(curr_rate, future_time, &coupon_times, coupon_rate)
            .unwrap()
    })
}

#[bench]
fn bench_coupon_bond_now(bench: &mut Bencher) {
    let fixture = BASELINE;
    let delta = fixture.delta;
    let future_time = 0.0;
    let coupon_rate = 0.05 * delta;
    let curve = fixture.curve();
    let hull_white = hull_white::HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();

    bench.iter(|| {
        let coupon_times = hull_white::get_coupon_times(5, future_time, delta).unwrap();
        hull_white
            .coupon_bond_price_now(&coupon_times, coupon_rate)
            .unwrap()
    })
}

#[bench]
fn bench_swap_rate(bench: &mut Bencher) {
    let fixture = FLAT_5PCT;
    let curr_rate = fixture.curr_rate;
    let delta = fixture.delta;
    let future_time = 0.0;
    //Same instrument as the swaption benches below: a 5y swap (1.0 -> 6.0) struck
    //off a 1y fixing.  Here `swap_initiation` == `option_maturity`.
    let option_maturity = 1.0;
    let swap_tenor = 5.0;
    let num_swap_payments = (swap_tenor / delta) as usize; //20 quarterly payments
    let curve = fixture.curve();
    let hull_white = hull_white::HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();
    bench.iter(|| {
        hull_white
            .forward_swap_rate_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
            )
            .unwrap()
    })
}

#[bench]
fn bench_swaption_european(bench: &mut Bencher) {
    let fixture = FLAT_5PCT;
    let curr_rate = fixture.curr_rate;
    let delta = fixture.delta;
    let future_time = 0.0;
    //Instrument: 1y European option on a 5y swap paying quarterly (swap legs run
    //option_maturity .. option_maturity + swap_tenor = 1.0 .. 6.0), struck ATM-forward.
    //`european_payer_swaption_t` takes the *option* expiry in the third slot; the swap
    //tenor only reaches the pricer through `num_swap_payments`.
    let option_maturity = 1.0;
    let swap_tenor = 5.0;
    let num_swap_payments = (swap_tenor / delta) as usize; //20 quarterly payments
    let curve = fixture.curve();
    let hull_white = hull_white::HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();
    let swap_rate = hull_white
        .forward_swap_rate_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
        )
        .unwrap();
    //Sanity for the published trend: ATM-forward on this fixture prices around
    //1.47e-2 of notional (american ~1.52e-2), so non-trivial and finite.
    let expected = hull_white
        .european_payer_swaption_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
            swap_rate,
        )
        .unwrap();
    assert!(expected.is_finite() && expected > 0.0, "{expected}");
    bench.iter(|| {
        hull_white
            .european_payer_swaption_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
                swap_rate,
            )
            .unwrap()
    })
}

#[bench]
fn bench_swaption_american(bench: &mut Bencher) {
    let fixture = FLAT_5PCT;
    let curr_rate = fixture.curr_rate;
    let delta = fixture.delta;
    let future_time = 0.0;
    //Instrument: 1y American option on the same 5y quarterly swap (1.0 .. 6.0),
    //struck ATM-forward.  `american_payer_swaption_t` takes the *option* expiry in
    //the third slot; the swap tenor only reaches the pricer through `num_swap_payments`.
    let option_maturity = 1.0;
    let swap_tenor = 5.0;
    let num_swap_payments = (swap_tenor / delta) as usize; //20 quarterly payments
    let curve = fixture.curve();
    let hull_white = hull_white::HullWhite::new(fixture.a, fixture.sigma, &curve).unwrap();
    let swap_rate = hull_white
        .forward_swap_rate_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
        )
        .unwrap();
    bench.iter(|| {
        hull_white
            .american_payer_swaption_t(
                curr_rate,
                future_time,
                option_maturity,
                num_swap_payments,
                delta,
                swap_rate,
                200,
            )
            .unwrap()
    })
}
