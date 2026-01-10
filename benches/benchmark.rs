use criterion::{black_box, criterion_group, criterion_main, Criterion};
use hull_white::{HullWhite, utils};

fn bench_bond_t(c: &mut Criterion) {
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
    let hull_white = HullWhite::new_panicking(a, sig, &yield_curve, &forward_curve);
    c.bench_function("bond_price_t", |b| {
        b.iter(|| hull_white.bond_price_t(black_box(curr_rate), black_box(future_time), black_box(maturity)))
    });
}

fn bench_bond_now(c: &mut Criterion) {
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
    let hull_white = HullWhite::new_panicking(a, sig, &yield_curve, &forward_curve);
    c.bench_function("bond_price_now", |b| {
        b.iter(|| hull_white.bond_price_now(black_box(maturity)))
    });
}

fn bench_coupon_bond_t(c: &mut Criterion) {
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
    let hull_white = HullWhite::new_panicking(a, sig, &yield_curve, &forward_curve);

    c.bench_function("coupon_bond_price_t", |b| {
        b.iter(|| {
            let coupon_times = utils::get_coupon_times(5, black_box(future_time), black_box(delta));
            hull_white.coupon_bond_price_t(
                black_box(curr_rate),
                black_box(future_time),
                &coupon_times,
                black_box(coupon_rate),
            )
        })
    });
}

fn bench_coupon_bond_now(c: &mut Criterion) {
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
    let hull_white = HullWhite::new_panicking(a, sig, &yield_curve, &forward_curve);

    c.bench_function("coupon_bond_price_now", |b| {
        b.iter(|| {
            let coupon_times = utils::get_coupon_times(5, black_box(future_time), black_box(delta));
            hull_white.coupon_bond_price_now(&coupon_times, black_box(coupon_rate))
        })
    });
}

fn bench_swap_rate(c: &mut Criterion) {
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
    let hull_white = HullWhite::new_panicking(a, sig, &yield_curve, &forward_curve);
    c.bench_function("swap_rate", |b| {
        b.iter(|| {
            hull_white.forward_swap_rate_t(
                black_box(curr_rate),
                black_box(future_time),
                black_box(option_maturity),
                black_box(num_swap_payments),
                black_box(delta),
            )
        })
    });
}

fn bench_swaption_european(c: &mut Criterion) {
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
    let hull_white = HullWhite::new_panicking(a, sig, &yield_curve, &forward_curve);
    let swap_rate = hull_white.forward_swap_rate_t(
        curr_rate,
        future_time,
        option_maturity,
        num_swap_payments,
        delta,
    );
    c.bench_function("european_swaption", |b| {
        b.iter(|| {
            let _ = hull_white.european_payer_swaption_t(
                black_box(curr_rate),
                black_box(future_time),
                black_box(swap_tenor),
                black_box(num_swap_payments),
                black_box(delta),
                black_box(swap_rate),
            );
        })
    });
}

criterion_group!(
    benches,
    bench_bond_t,
    bench_bond_now,
    bench_coupon_bond_t,
    bench_coupon_bond_now,
    bench_swap_rate,
    bench_swaption_european
);
criterion_main!(benches);