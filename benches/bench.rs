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
    let hull_white = hull_white::HullWhite::init(a, sig, &yield_curve, &forward_curve);
    bench.iter(|| hull_white.bond_price_t(curr_rate, future_time, maturity))
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
    let hull_white = hull_white::HullWhite::init(a, sig, &yield_curve, &forward_curve);
    bench.iter(|| hull_white.bond_price_now(maturity))
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
    let hull_white = hull_white::HullWhite::init(a, sig, &yield_curve, &forward_curve);

    bench.iter(|| {
        let coupon_times = hull_white::get_coupon_times(5, future_time, delta);
        hull_white.coupon_bond_price_t(curr_rate, future_time, &coupon_times, coupon_rate)
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
    let hull_white = hull_white::HullWhite::init(a, sig, &yield_curve, &forward_curve);

    bench.iter(|| {
        let coupon_times = hull_white::get_coupon_times(5, future_time, delta);
        hull_white.coupon_bond_price_now(&coupon_times, coupon_rate)
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
    let num_swap_payments = 20;
    let option_maturity = 1.0;
    let yield_curve = move |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = move |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let hull_white = hull_white::HullWhite::init(a, sig, &yield_curve, &forward_curve);
    bench.iter(|| {
        hull_white.forward_swap_rate_t(
            curr_rate,
            future_time,
            option_maturity,
            num_swap_payments,
            delta,
        )
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
    let swap_tenor = 5.0;
    let option_maturity = 1.0;
    let num_swap_payments = 20;
    let yield_curve = move |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = move |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let hull_white = hull_white::HullWhite::init(a, sig, &yield_curve, &forward_curve);
    let swap_rate = hull_white.forward_swap_rate_t(
        curr_rate,
        future_time,
        option_maturity,
        num_swap_payments,
        delta,
    );
    bench.iter(|| {
        hull_white.european_payer_swaption_t(
            curr_rate,
            future_time,
            swap_tenor,
            num_swap_payments,
            delta,
            swap_rate,
        )
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
    let swap_tenor = 5.0;
    let option_maturity = 1.0;
    let num_swap_payments = 20;
    let yield_curve = move |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct = (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = move |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let hull_white = hull_white::HullWhite::init(a, sig, &yield_curve, &forward_curve);
    let swap_rate = hull_white.forward_swap_rate_t(
        curr_rate,
        future_time,
        option_maturity,
        num_swap_payments,
        delta,
    );
    bench.iter(|| {
        hull_white.american_payer_swaption_t(
            curr_rate,
            future_time,
            swap_tenor,
            num_swap_payments,
            delta,
            swap_rate,
            200,
        )
    })
}
