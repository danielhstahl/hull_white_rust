#![feature(test)]
extern crate test;

use test::Bencher;

#[bench]
fn bench_bond_t(bench: &mut Bencher) {
    let curr_rate = 0.02;
    let sig: f64 = 0.02;
    let a: f64 = 0.3;
    let b = 0.04;
    let future_time = 0.0;
    let maturity = 1.5;
    let yield_curve = move |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = move |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let hull_white = hull_white::HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
    bench.iter(|| {
        hull_white
            .bond_price_t(curr_rate, future_time, maturity)
            .unwrap()
    })
}

#[bench]
fn bench_bond_now(bench: &mut Bencher) {
    let curr_rate = 0.02;
    let sig: f64 = 0.02;
    let a: f64 = 0.3;
    let b = 0.04;
    let maturity = 1.5;
    let yield_curve = move |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = move |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let hull_white = hull_white::HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
    bench.iter(|| hull_white.bond_price_now(maturity).unwrap())
}

#[bench]
fn bench_coupon_bond_t(bench: &mut Bencher) {
    let curr_rate = 0.02;
    let sig: f64 = 0.02;
    let a: f64 = 0.3;
    let delta = 0.25;
    let b = 0.04;
    let future_time = 0.0;
    let coupon_rate = 0.05 * delta;
    let yield_curve = move |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = move |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let hull_white = hull_white::HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();

    bench.iter(|| {
        let coupon_times = hull_white::get_coupon_times(5, future_time, delta).unwrap();
        hull_white
            .coupon_bond_price_t(curr_rate, future_time, &coupon_times, coupon_rate)
            .unwrap()
    })
}

#[bench]
fn bench_coupon_bond_now(bench: &mut Bencher) {
    let curr_rate = 0.02;
    let sig: f64 = 0.02;
    let a: f64 = 0.3;
    let delta = 0.25;
    let b = 0.04;
    let future_time = 0.0;
    let coupon_rate = 0.05 * delta;
    let yield_curve = move |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = move |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let hull_white = hull_white::HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();

    bench.iter(|| {
        let coupon_times = hull_white::get_coupon_times(5, future_time, delta).unwrap();
        hull_white
            .coupon_bond_price_now(&coupon_times, coupon_rate)
            .unwrap()
    })
}

#[bench]
fn bench_swap_rate(bench: &mut Bencher) {
    let curr_rate = 0.05;
    let sig: f64 = 0.01;
    let a: f64 = 0.05;
    let b = 0.05;
    let delta = 0.25;
    let future_time = 0.0;
    //Same instrument as the swaption benches below: a 5y swap (1.0 -> 6.0) struck
    //off a 1y fixing.  Here `swap_initiation` == `option_maturity`.
    let option_maturity = 1.0;
    let swap_tenor = 5.0;
    let num_swap_payments = (swap_tenor / delta) as usize; //20 quarterly payments
    let yield_curve = move |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = move |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let hull_white = hull_white::HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
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
    let curr_rate = 0.05;
    let sig: f64 = 0.01;
    let a: f64 = 0.05;
    let b = 0.05;
    let delta = 0.25;
    let future_time = 0.0;
    //Instrument: 1y European option on a 5y swap paying quarterly (swap legs run
    //option_maturity .. option_maturity + swap_tenor = 1.0 .. 6.0), struck ATM-forward.
    //`european_payer_swaption_t` takes the *option* expiry in the third slot; the swap
    //tenor only reaches the pricer through `num_swap_payments`.
    let option_maturity = 1.0;
    let swap_tenor = 5.0;
    let num_swap_payments = (swap_tenor / delta) as usize; //20 quarterly payments
    let yield_curve = move |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = move |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let hull_white = hull_white::HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
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
    let curr_rate = 0.05;
    let sig: f64 = 0.01;
    let a: f64 = 0.05;
    let b = 0.05;
    let delta = 0.25;
    let future_time = 0.0;
    //Instrument: 1y American option on the same 5y quarterly swap (1.0 .. 6.0),
    //struck ATM-forward.  `american_payer_swaption_t` takes the *option* expiry in
    //the third slot; the swap tenor only reaches the pricer through `num_swap_payments`.
    let option_maturity = 1.0;
    let swap_tenor = 5.0;
    let num_swap_payments = (swap_tenor / delta) as usize; //20 quarterly payments
    let yield_curve = move |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = move |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let hull_white = hull_white::HullWhite::init(a, sig, &yield_curve, &forward_curve).unwrap();
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
