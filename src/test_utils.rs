#[cfg(test)]
mod tests {
    use super::super::*;
    use crate::utils;
    use approx::*;

    #[test]
    fn test_get_num_payments_if_exact_integer() {
        let t: f64 = 0.5;
        let maturity: f64 = 2.0;
        let delta: f64 = 0.25;
        let (num_payments, is_exact) = utils::get_num_remaining_payments(t, maturity, delta);
        assert_eq!(num_payments, 6);
        assert!(is_exact);
        let next_exchange_date = maturity - (num_payments as f64 - 1.0) * delta;
        assert_abs_diff_eq!(next_exchange_date, t + delta, epsilon = 1e-7);
    }

    #[test]
    fn test_get_num_payments_not_even() {
        let t: f64 = 0.5;
        let maturity: f64 = 2.0;
        let delta: f64 = 0.4;
        let (num_payments, is_exact) = utils::get_num_remaining_payments(t, maturity, delta);
        assert_eq!(num_payments, 4);
        assert_eq!(is_exact, false);
        let next_exchange_date = maturity - (num_payments as f64 - 1.0) * delta;
        assert_abs_diff_eq!(next_exchange_date, 0.8, epsilon = 1e-7);
    }

    #[test]
    fn test_get_time_from_index() {
        let t: f64 = 0.5;
        let delta: f64 = 0.4;
        let index: usize = 3;
        let time = utils::get_time_from_index(index, t, delta);
        assert_abs_diff_eq!(time, 1.7, epsilon = 1e-7);
    }

    #[test]
    fn test_get_time_from_index_with_zero_index() {
        let t: f64 = 0.5;
        let delta: f64 = 0.4;
        let index: usize = 0;
        let time = utils::get_time_from_index(index, t, delta);
        assert_eq!(time, t);
    }

    #[test]
    fn test_get_coupon_times() {
        let num_payments: usize = 5;
        let t: f64 = 1.0;
        let delta: f64 = 0.25;
        let coupon_times = utils::get_coupon_times(num_payments, t, delta);
        let expected_coupon_times = vec![1.25, 1.5, 1.75, 2.0, 2.25];
        coupon_times
            .iter()
            .zip(expected_coupon_times.iter())
            .for_each(|(actual, expected)| assert_eq!(actual, expected))
    }

    #[test]
    fn test_get_coupon_times_no_payments() {
        let num_payments: usize = 0;
        let t: f64 = 1.0;
        let delta: f64 = 0.25;
        let coupon_times = utils::get_coupon_times(num_payments, t, delta);
        assert_eq!(coupon_times.len(), 0);
    }

    #[test]
    fn test_max_or_zero() {
        let v: f64 = 1.0;
        assert_eq!(utils::max_or_zero(v), 1.0);
        assert_eq!(utils::max_or_zero(-v), 0.0);
    }

    #[test]
    fn test_payoff_swaption() {
        let v: f64 = 1.0;
        assert_eq!(utils::payoff_swaption(true, v), 1.0);
        assert_eq!(utils::payoff_swaption(false, v), 0.0);
        assert_eq!(utils::payoff_swaption(true, -v), 0.0);
        assert_eq!(utils::payoff_swaption(false, -v), 1.0);
    }

    #[test]
    fn test_hull_white_creation() {
        let a: f64 = 0.1;
        let sigma: f64 = 0.01;
        let yield_curve = |t: f64| 0.05 * t;
        let forward_curve = |t: f64| t.ln();
        
        let hw = HullWhite::new_panicking(a, sigma, &yield_curve, &forward_curve);
        
        assert_eq!(hw.a, a);
        assert_eq!(hw.sigma, sigma);
    }

    #[test]
    fn test_bond_price_at_expiry() {
        let curr_rate: f64 = 0.02;
        let sigma: f64 = 0.02;
        let a: f64 = 0.3;
        let b: f64 = 0.04;
        let future_time: f64 = 0.5;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct = (b - sigma.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sigma * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sigma.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let hull_white = HullWhite::new_panicking(a, sigma, &yield_curve, &forward_curve);
        assert_abs_diff_eq!(
            hull_white.bond_price_t(curr_rate, future_time, future_time),
            1.0,
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_bond_now_same_as_t_when_t_is_zero() {
        let curr_rate: f64 = 0.02;
        let sigma: f64 = 0.02;
        let a: f64 = 0.3;
        let b: f64 = 0.04;
        let future_time: f64 = 0.0;
        let maturity: f64 = 1.5;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct = (b - sigma.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sigma * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sigma.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let hull_white = HullWhite::new_panicking(a, sigma, &yield_curve, &forward_curve);
        let bond_price_now = hull_white.bond_price_now(maturity);
        let bond_price_t = hull_white.bond_price_t(curr_rate, future_time, maturity);
        assert_abs_diff_eq!(bond_price_now, bond_price_t, epsilon = 1e-7);
    }

    #[test]
    fn test_coupon_bond_now_same_as_t_when_t_is_zero() {
        let curr_rate: f64 = 0.02;
        let sigma: f64 = 0.02;
        let a: f64 = 0.3;
        let b: f64 = 0.04;
        let delta: f64 = 0.25;
        let future_time: f64 = 0.0;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct = (b - sigma.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sigma * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sigma.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let coupon_times = utils::get_coupon_times(6, future_time, delta); // this was 5, but made six since last payment is now included
        let coupon_rate: f64 = 0.05 * delta;
        let hull_white = HullWhite::new_panicking(a, sigma, &yield_curve, &forward_curve);
        let bond_price_now = hull_white.coupon_bond_price_now(&coupon_times, coupon_rate);
        let bond_price_t = hull_white.coupon_bond_price_t(curr_rate, future_time, &coupon_times, coupon_rate);
        assert_abs_diff_eq!(bond_price_now, bond_price_t, epsilon = 1e-7);
    }

    #[test]
    fn test_swap() {
        let curr_rate: f64 = 0.02;
        let sigma: f64 = 0.02;
        let a: f64 = 0.3;
        let b: f64 = 0.04;
        let delta: f64 = 0.25;
        let future_time: f64 = 0.5;
        let swap_maturity: f64 = 5.5;
        let num_swap_payments: usize = 20;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct = (b - sigma.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sigma * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sigma.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let hull_white = HullWhite::new_panicking(a, sigma, &yield_curve, &forward_curve);
        assert_abs_diff_eq!(
            hull_white.swap_price_t(
                curr_rate,
                future_time,
                swap_maturity,
                delta,
                hull_white.swap_rate_t(curr_rate, future_time, num_swap_payments, delta)
            ),
            0.0,
            epsilon = 1e-9
        );
    }

    #[test]
    fn test_swap_init() {
        let curr_rate: f64 = 0.02;
        let sigma: f64 = 0.02;
        let a: f64 = 0.3;
        let b: f64 = 0.04;
        let delta: f64 = 0.25;
        let future_time: f64 = 0.5;
        let swap_maturity: f64 = 5.5;
        let num_swap_payments: usize = 20;
        let swap_rate: f64 = 0.03;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct = (b - sigma.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sigma * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sigma.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let hull_white = HullWhite::new_panicking(a, sigma, &yield_curve, &forward_curve);
        let sp_init = hull_white.swap_price_t_init(
            curr_rate,
            future_time,
            future_time,
            num_swap_payments,
            delta,
            swap_rate,
        );
        let sp = hull_white.swap_price_t(curr_rate, future_time, swap_maturity, delta, swap_rate);
        assert_eq!(sp_init, sp);
    }

    #[test]
    fn test_zero_coupon_reference() {
        // http://www.quantcalc.net/BondOption_Vasicek.html
        let curr_rate: f64 = 0.01;
        let sigma: f64 = 0.03;
        let a: f64 = 0.05;
        let b: f64 = 0.04;
        let strike: f64 = 0.96;
        let future_time: f64 = 0.0;
        let bond_maturity: f64 = 3.0;
        let option_maturity: f64 = 2.0;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct = (b - sigma.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sigma * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sigma.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let hull_white = HullWhite::new_panicking(a, sigma, &yield_curve, &forward_curve);
        let bond_call = hull_white.bond_call_t(
            curr_rate,
            future_time,
            option_maturity,
            bond_maturity,
            strike,
        );
        assert_abs_diff_eq!(bond_call, 0.033282, epsilon = 1e-4)
    }

    #[test]
    fn test_coupon_bond_call_with_zero_coupon() {
        let curr_rate: f64 = 0.01;
        let sigma: f64 = 0.03;
        let a: f64 = 0.05;
        let b: f64 = 0.04;
        let strike: f64 = 0.96;
        let future_time: f64 = 0.0;
        let bond_maturity: f64 = 3.0;
        let option_maturity: f64 = 2.0;
        let yield_curve = |t: f64| {
            let at = (1.0 - (-a * t).exp()) / a;
            let ct = (b - sigma.powi(2) / (2.0 * a.powi(2))) * (at - t) - (sigma * at).powi(2) / (4.0 * a);
            at * curr_rate - ct
        };
        let forward_curve = |t: f64| {
            b + (-a * t).exp() * (curr_rate - b)
                - (sigma.powi(2) / (2.0 * a.powi(2))) * (1.0 - (-a * t).exp()).powi(2)
        };
        let hull_white = HullWhite::new_panicking(a, sigma, &yield_curve, &forward_curve);
        let bond_call = hull_white.bond_call_t(
            curr_rate,
            future_time,
            option_maturity,
            bond_maturity,
            strike,
        );
        let coupon_rate: f64 = 0.0;
        let coupon_bond_call = hull_white
            .coupon_bond_call_t(
                curr_rate,
                future_time,
                option_maturity,
                &[bond_maturity],
                coupon_rate,
                strike,
            )
            .unwrap();

        assert_abs_diff_eq!(bond_call, coupon_bond_call, epsilon = 1e-4)
    }
}