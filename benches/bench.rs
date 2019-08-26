
#[macro_use]
extern crate criterion;
use criterion::{Criterion, ParameterizedBenchmark, Benchmark};
fn bench_bond(c: &mut Criterion){
    let curr_rate = 0.02;
    let sig: f64 = 0.02;
    let a: f64 = 0.3;
    let b = 0.04;
    let future_time = 0.0;
    let maturity = 1.5;
    let yield_curve = move |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct =
            (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = move |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    c.bench("compare bonds", Benchmark::new("bond_t",move |b|{
        b.iter(||{
             hull_white::bond_price_t(
                curr_rate,
                a,
                sig,
                future_time,
                maturity,
                &yield_curve,
                &forward_curve,
            )
        })
    }).with_function("bond_now", move |b|{
        b.iter(||{
            hull_white::bond_price_now(maturity, &yield_curve)
        })
    }));
}

fn bench_coupon_bond(c: &mut Criterion){
    let curr_rate = 0.02;
    let sig: f64 = 0.02;
    let a: f64 = 0.3;
    let delta = 0.25;
    let b = 0.04;
    let future_time = 0.0;
    let coupon_rate = 0.05 * delta;
    let yield_curve = move |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct =
            (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = move |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    c.bench("compare coupon bonds", Benchmark::new("bond_t",move |b|{
        let coupon_times = hull_white::get_coupon_times(5, future_time, delta);
        let last_coupon=*coupon_times.last().unwrap() + delta;
        b.iter(||{
             hull_white::coupon_bond_price_t(
                curr_rate,
                a,
                sig,
                future_time,
                &coupon_times,
                last_coupon,
                coupon_rate,
                &yield_curve,
                &forward_curve,
            )
        })
    }).with_function("bond_now", move |b|{
        let coupon_times = hull_white::get_coupon_times(5, future_time, delta);
        let last_coupon=*coupon_times.last().unwrap() + delta;
        b.iter(||{
            hull_white::coupon_bond_price_now(
                &coupon_times,
                last_coupon,
                coupon_rate,
                &yield_curve,
            )
        })
    }));
}

fn bench_swap_rate(c: &mut Criterion){
    let curr_rate = 0.05;
    let sig: f64 = 0.01;
    let a: f64 = 0.05;
    let b = 0.05;
    let delta = 0.25;
    let future_time = 0.0;
    let swap_tenor = 5.0;
    let option_maturity = 1.0;
    let yield_curve = move |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct =
            (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = move |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    c.bench_function("swap_rate", move |b|{
        b.iter(||{
            hull_white::forward_swap_rate_t(
                curr_rate,
                a,
                sig,
                future_time,
                option_maturity,
                swap_tenor + option_maturity,
                delta,
                &yield_curve,
                &forward_curve,
            )
        })
    });
}

fn bench_swaption_bond(c: &mut Criterion){
    let curr_rate = 0.05;
    let sig: f64 = 0.01;
    let a: f64 = 0.05;
    let b = 0.05;
    let delta = 0.25;
    let future_time = 0.0;
    let swap_tenor = 5.0;
    let option_maturity = 1.0;
    let yield_curve = move |t: f64| {
        let at = (1.0 - (-a * t).exp()) / a;
        let ct =
            (b - sig.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sig * at).powi(2) / (4.0 * a);
        at * curr_rate - ct
    };
    let forward_curve = move |t: f64| {
        b + (-a * t).exp() * (curr_rate - b)
            - (sig.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
    };
    let swap_rate = hull_white::forward_swap_rate_t(
        curr_rate,
        a,
        sig,
        future_time,
        option_maturity,
        swap_tenor + option_maturity,
        delta,
        &yield_curve,
        &forward_curve,
    );
    c.bench("compare swaption tree", ParameterizedBenchmark::new("analytical",move |b, _i|{
        b.iter(||{
            hull_white::european_payer_swaption_t(
                curr_rate,
                a,
                sig,
                future_time,
                swap_tenor,
                option_maturity,
                delta,
                swap_rate,
                &yield_curve,
                &forward_curve,
            )
        })
    }, vec![50, 100, 200]).with_function("tree", move |b, i|{
        b.iter(||{
            hull_white::american_payer_swaption_t(
                curr_rate,
                a,
                sig,
                future_time,
                swap_tenor,
                option_maturity,
                delta,
                swap_rate,
                *i,
                &yield_curve,
                &forward_curve,
            )
        })
    }));
}


criterion_group!(benches, bench_bond, bench_coupon_bond, bench_swap_rate, bench_swaption_bond);
criterion_main!(benches);