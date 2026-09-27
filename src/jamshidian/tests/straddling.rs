//! Options whose coupon schedule straddles the option's expiry date.
//!
//! The convention under test is the one in the [module docs](crate::jamshidian): the underlying is
//! the bond *as it stands on the expiry date*, so payments made strictly before expiry are
//! dropped, a payment falling exactly on expiry is cash that folds into the strike, and only
//! payments strictly after expiry are decomposed.
//!
//! Every reference here is written from the schedule by way of the payoff, not from the
//! decomposition: [`deliverable_at_expiry`] rebuilds what the exercising holder has, then
//! [`super::payoff_integral`] integrates that payoff over the expiry-date distribution of the
//! rate, and [`tree_option_price`] rolls the same payoff back through a short-rate lattice.  The
//! two references share no code with each other or with the pricer; the lattice in
//! [`tree_option_price`] is the crate's shared tree engine ([`HullWhite::tree_price`]), which
//! knows nothing about bonds, coupons or Jamshidian -- it is the time-coordinate plumbing, and
//! making it one copy is what stops a fix having to be applied here as well as in `crate::trees`.

use crate::HullWhite;
use crate::test_support::{STEEP_CURVE, hw_curve};

use super::{FIXTURES, payoff_integral};

/// The expiry date every schedule in this module straddles, with the valuation at `t = 1.0`.
const U: f64 = 2.0;
const T: f64 = 1.0;
const R_T: f64 = 0.04;
const COUPON_RATE: f64 = 0.05;

/// Schedules that straddle `U`: payments before it, payments on it, payments after it.
const STRADDLING: [&[f64]; 5] = [
    &[1.25, 1.5, 2.25, 2.5, 3.0], //before and after, nothing on the expiry date
    &[1.25, 1.5, 2.0, 2.25, 2.5, 3.0], //before, on, and after
    &[1.9, 2.0, 2.05, 3.0],       //tightly around the expiry date
    &[2.0, 2.5, 3.0],             //cash on the expiry date plus a tail
    &[1.99, 2.0000001, 3.0],      //either side of it by a whisker
];

/// The same schedule with the pre-expiry payments removed -- the tail the holder actually gets,
/// for the schedules that carry nothing exactly on `U`.
fn tail_only(coupon_times: &[f64], expiry: f64) -> &[f64] {
    let first_kept = coupon_times.partition_point(|&time| time < expiry);
    &coupon_times[first_kept..]
}

/// What the exercising holder has on the expiry date, at expiry rate `rate`, read off the whole
/// schedule by hand: nothing before expiry, cash (`P(u, u) = 1`) on it, the discount bond on each
/// later payment.
fn deliverable_at_expiry(
    hull_white: &HullWhite,
    expiry: f64,
    coupon_times: &[f64],
    coupon_rate: f64,
    rate: f64,
) -> f64 {
    let last_index = coupon_times.len() - 1;
    coupon_times
        .iter()
        .enumerate()
        .map(|(index, &time)| {
            let weight = coupon_rate + if index == last_index { 1.0 } else { 0.0 };
            if time > expiry {
                weight * hull_white.bond_price_t(rate, expiry, time).unwrap()
            } else if time == expiry {
                weight
            } else {
                0.0
            }
        })
        .sum()
}

/// A short-rate tree reference, sharing nothing with the decomposition: the rate walks a
/// `binomial_tree` Black-Vasicek lattice from `t` to expiry, the deliverable is re-valued at each
/// expiry node, and the payoff is rolled back with each node's own discount factor.
///
/// The lattice is [`HullWhite::tree_price`], the engine the swaption trees run on, so the
/// tree-time-versus-absolute-time convention that `phi` is cached against is written once for the
/// whole crate (this function used to carry its own copy of it, which is how the same bug got paid
/// for three times).  What stays independent is the thing under test: the payoff, rebuilt from the
/// schedule by [`deliverable_at_expiry`] rather than decomposed.
#[allow(clippy::too_many_arguments)] //a pricer's instrument plus the tree's resolution
fn tree_option_price(
    hull_white: &HullWhite,
    r_t: f64,
    t: f64,
    option_maturity: f64,
    underlying_at_expiry: &dyn Fn(f64) -> f64,
    strike: f64,
    is_call: bool,
    num_steps: usize,
) -> f64 {
    hull_white.tree_price(
        r_t,
        t,
        option_maturity - t,
        num_steps,
        false,
        &|_tau: f64, rate: f64, _dt: f64, _j: usize| {
            let payoff = if is_call {
                underlying_at_expiry(rate) - strike
            } else {
                strike - underlying_at_expiry(rate)
            };
            if payoff > 0.0 { payoff } else { 0.0 }
        },
    )
}

fn assert_matches_payoff_integral(
    hull_white: &HullWhite,
    coupon_times: &[f64],
    coupon_rate: f64,
    strike: f64,
    is_call: bool,
) {
    let priced = if is_call {
        hull_white
            .coupon_bond_call_t(R_T, T, U, coupon_times, coupon_rate, strike)
            .unwrap()
    } else {
        hull_white
            .coupon_bond_put_t(R_T, T, U, coupon_times, coupon_rate, strike)
            .unwrap()
    };
    let deliverable =
        |rate: f64| deliverable_at_expiry(hull_white, U, coupon_times, coupon_rate, rate);
    let reference = payoff_integral(hull_white, R_T, T, U, &deliverable, strike, is_call);
    assert!(
        (priced - reference).abs() <= 1e-10f64.max(reference.abs() * 1e-9),
        "schedule {coupon_times:?} {} strike {strike}: {priced} vs {}",
        if is_call { "call" } else { "put" },
        reference,
    );
}

#[test]
fn a_straddling_schedule_matches_the_direct_payoff_integral() {
    //The headline check for the convention: the decomposition of the residual bond, against a
    //quadrature of the payoff the holder actually has.
    let strikes = [0.3f64, 0.7, 0.95, 1.0, 1.05, 1.3, 2.0];
    for s in FIXTURES.iter() {
        let curve = s.curve();
        let hull_white = HullWhite::new(s.a, s.sigma, &curve).unwrap();
        for schedule in STRADDLING.iter() {
            for &strike in strikes.iter() {
                for is_call in [true, false] {
                    assert_matches_payoff_integral(
                        &hull_white,
                        schedule,
                        COUPON_RATE,
                        strike,
                        is_call,
                    );
                }
            }
        }
    }
}

#[test]
fn a_straddling_schedule_with_a_negative_coupon_rate_matches_the_integral() {
    //A negative coupon is not nonsense (deep-discount and cross-currency legs pay one) and it makes
    //the split non-trivial: the dropped payments were *negative* cash, the payment on the expiry
    //date arrives as a strike increase, and the par still has to land on the schedule's final
    //payment rather than on whatever happens to be the residual's last one.
    let curve = STEEP_CURVE.curve();
    let hull_white = HullWhite::new(STEEP_CURVE.a, STEEP_CURVE.sigma, &curve).unwrap();
    for schedule in STRADDLING.iter() {
        for coupon_rate in [-0.02f64, -0.005, -0.0001] {
            for strike in [0.5f64, 0.9, 0.95, 1.0, 1.2] {
                for is_call in [true, false] {
                    assert_matches_payoff_integral(
                        &hull_white,
                        schedule,
                        coupon_rate,
                        strike,
                        is_call,
                    );
                }
            }
        }
    }
}

#[test]
fn dropping_the_pre_expiry_coupons_cannot_change_the_price() {
    //Everything strictly before expiry is worth nothing on the expiry date, so dropping it is a
    //no-op: the full schedule and the post-expiry tail are the same instrument.
    for s in FIXTURES.iter() {
        let curve = s.curve();
        let hull_white = HullWhite::new(s.a, s.sigma, &curve).unwrap();
        for schedule in STRADDLING.iter() {
            let tail = tail_only(schedule, U);
            if tail == *schedule {
                continue; //nothing was dropped; nothing to compare
            }
            for &strike in [0.7f64, 0.95, 1.0, 1.3].iter() {
                let full_call = hull_white
                    .coupon_bond_call_t(R_T, T, U, schedule, COUPON_RATE, strike)
                    .unwrap();
                let tail_call = hull_white
                    .coupon_bond_call_t(R_T, T, U, tail, COUPON_RATE, strike)
                    .unwrap();
                assert_eq!(
                    full_call, tail_call,
                    "call on {schedule:?} vs tail {tail:?}"
                );
                let full_put = hull_white
                    .coupon_bond_put_t(R_T, T, U, schedule, COUPON_RATE, strike)
                    .unwrap();
                let tail_put = hull_white
                    .coupon_bond_put_t(R_T, T, U, tail, COUPON_RATE, strike)
                    .unwrap();
                assert_eq!(full_put, tail_put, "put on {schedule:?} vs tail {tail:?}");
            }
        }
    }
}

#[test]
fn a_coupon_on_the_expiry_date_is_a_strike_reduction() {
    //A payment settled exactly on the expiry date is worth its weight there whatever the rate does,
    //so pricing the schedule that carries it must equal pricing the residual against a strike
    //shrunk by that cash: (cash + R - K)^+ == (R - (K - cash))^+.
    for s in FIXTURES.iter() {
        let curve = s.curve();
        let hull_white = HullWhite::new(s.a, s.sigma, &curve).unwrap();
        let residual = &[2.25, 2.5, 3.0][..];
        let carrying_cash = &[2.0, 2.25, 2.5, 3.0][..];
        for &strike in [0.7f64, 0.95, 1.0, 1.3].iter() {
            let folded_call = hull_white
                .coupon_bond_call_t(R_T, T, U, carrying_cash, COUPON_RATE, strike)
                .unwrap();
            let plain_call = hull_white
                .coupon_bond_call_t(R_T, T, U, residual, COUPON_RATE, strike - COUPON_RATE)
                .unwrap();
            assert_eq!(folded_call, plain_call, "call, cash folded into the strike");
            let folded_put = hull_white
                .coupon_bond_put_t(R_T, T, U, carrying_cash, COUPON_RATE, strike)
                .unwrap();
            let plain_put = hull_white
                .coupon_bond_put_t(R_T, T, U, residual, COUPON_RATE, strike - COUPON_RATE)
                .unwrap();
            assert_eq!(folded_put, plain_put, "put, cash folded into the strike");
        }
    }
}

#[test]
fn a_strike_at_or_below_the_expiry_cash_is_exercised_in_every_state() {
    //The cash settled on the expiry date is a floor on the deliverable that no rate can remove, so
    //at a strike at or below it the call is parity and the put is nothing, with nothing left to
    //solve for.  The quadrature checks the same thing from the payoff.
    for s in FIXTURES.iter() {
        let curve = s.curve();
        let hull_white = HullWhite::new(s.a, s.sigma, &curve).unwrap();
        let carrying_cash = &[2.0, 2.25, 2.5, 3.0][..];
        //The floor is one coupon here, and the deliverable really is above it at every rate the
        //option can reach: sample the whole state space the quadrature integrates over.
        let floor = COUPON_RATE;
        let (mean, sd) = super::expiry_rate_moments(&hull_white, R_T, T, U);
        for step in 0..=480 {
            let rate = mean + sd * (-12.0 + step as f64 * 0.05);
            assert!(
                deliverable_at_expiry(&hull_white, U, carrying_cash, COUPON_RATE, rate)
                    > floor - 1e-16,
                "deliverable dips below the expiry cash at rate {rate}"
            );
        }
        for strike in [0.0f64, 0.01, 0.04, 0.05] {
            let call = hull_white
                .coupon_bond_call_t(R_T, T, U, carrying_cash, COUPON_RATE, strike)
                .unwrap();
            let put = hull_white
                .coupon_bond_put_t(R_T, T, U, carrying_cash, COUPON_RATE, strike)
                .unwrap();
            assert_eq!(
                put, 0.0,
                "strike {strike}: a put on a deliverable the strike cannot beat"
            );
            assert_matches_payoff_integral(&hull_white, carrying_cash, COUPON_RATE, strike, true);
            //Call minus put is the deliverable minus the strike discounted.
            let deliverable_now = hull_white
                .coupon_bond_price_t(R_T, T, &[2.25, 2.5, 3.0], COUPON_RATE)
                .unwrap()
                + COUPON_RATE * hull_white.bond_price_t(R_T, T, U).unwrap();
            let parity = deliverable_now - strike * hull_white.bond_price_t(R_T, T, U).unwrap();
            assert!(
                (call - parity).abs() <= 1e-12f64.max(parity.abs() * 1e-11),
                "strike {strike}: call {call} vs parity {parity}"
            );
        }
    }
}

#[test]
fn a_straddling_schedule_agrees_with_a_short_rate_tree() {
    //A third route: the same payoff rolled back through a short-rate lattice, which discretises the
    //rate process instead of integrating its terminal law.  Tree error is O(steps^-1)-ish, so this
    //checks the convention rather than the last decimal: for this grid the worst |tree - analytic|
    //over 100..400 steps is ~2e-4 and it shrinks monotonically with the step count, so 1e-3 is
    //the tolerance and the tree's own convergence is asserted separately.
    let strikes = [0.7f64, 0.95, 1.0, 1.3];
    for s in FIXTURES.iter() {
        let curve = s.curve();
        let hull_white = HullWhite::new(s.a, s.sigma, &curve).unwrap();
        for schedule in STRADDLING.iter() {
            let deliverable =
                |rate: f64| deliverable_at_expiry(&hull_white, U, schedule, COUPON_RATE, rate);
            for &strike in strikes.iter() {
                for is_call in [true, false] {
                    let priced = if is_call {
                        hull_white
                            .coupon_bond_call_t(R_T, T, U, schedule, COUPON_RATE, strike)
                            .unwrap()
                    } else {
                        hull_white
                            .coupon_bond_put_t(R_T, T, U, schedule, COUPON_RATE, strike)
                            .unwrap()
                    };
                    let tree = tree_option_price(
                        &hull_white,
                        R_T,
                        T,
                        U,
                        &deliverable,
                        strike,
                        is_call,
                        400,
                    );
                    let finer = tree_option_price(
                        &hull_white,
                        R_T,
                        T,
                        U,
                        &deliverable,
                        strike,
                        is_call,
                        800,
                    );
                    assert!(
                        (priced - tree).abs() <= 1e-4f64.max(priced.abs() * 1e-2),
                        "schedule {schedule:?} {} strike {strike}: analytic {priced} vs tree {tree}",
                        if is_call { "call" } else { "put" },
                    );
                    assert!(
                        (finer - priced).abs() <= 1e-4f64.max(priced.abs() * 1e-2),
                        "finer tree disagrees with the analytic price: {schedule:?} {strike} \
                         {} {tree} -> {finer} vs {priced}",
                        if is_call { "call" } else { "put" },
                    );
                }
            }
        }
    }
}

#[test]
//The values below are captured at full 17-digit precision on purpose: they are the point of the
//test, and the assertion tolerance is far looser than the last digit anyway.
#[allow(clippy::excessive_precision)]
fn an_all_post_expiry_schedule_is_untouched_by_the_convention() {
    //The case that worked before the schedule split must keep working, to the bits.  These values
    //were produced by the pre-change pricers on the same inputs (option_maturity 1.5, schedule
    //1.75..2.5, all strictly after expiry), captured with `{:.17}` and reproduced here unchanged.
    let times = [1.75, 2.0, 2.25, 2.5];
    #[rustfmt::skip]
    let golden: [(f64, f64, f64, f64, f64, f64, f64); 12] = [
        // curr,    a,    b,  sigma,   strike,          call,             put
        (0.05, 0.05, 0.05, 0.01, 0.7, 0.44635401121770413, 0.0),
        (0.05, 0.05, 0.05, 0.01, 0.95, 0.20131903016022495, 0.0),
        (0.05, 0.05, 0.05, 0.01, 1.0, 0.15231203394872903, 0.0),
        (0.05, 0.05, 0.05, 0.01, 1.3, 0.0, 0.14172994332024591),
        (0.02, 0.20, 0.06, 0.03, 0.7, 0.44327125334607381, 0.0),
        (0.02, 0.20, 0.06, 0.03, 0.95, 0.19833583354304618, 0.0),
        (0.02, 0.20, 0.06, 0.03, 1.0, 0.14934874958244046, -0.00000000000000001),
        (0.02, 0.20, 0.06, 0.03, 1.3, 0.00000000000000369, 0.14457375418119650),
        (0.03, 0.40, 0.045, 0.005, 0.7, 0.44503883261650723, 0.0),
        (0.03, 0.40, 0.045, 0.005, 0.95, 0.20004642200790093, 0.0),
        (0.03, 0.40, 0.045, 0.005, 1.0, 0.15104793988617959, 0.0),
        (0.03, 0.40, 0.045, 0.005, 1.3, 0.0, 0.14294295284414782),
    ];
    for &(curr, a, b, sigma, strike, call_expected, put_expected) in golden.iter() {
        let curve = hw_curve(curr, a, b, sigma);
        let hull_white = HullWhite::new(a, sigma, &curve).unwrap();
        let call = hull_white
            .coupon_bond_call_t(R_T, 1.0, 1.5, &times, COUPON_RATE, strike)
            .unwrap();
        let put = hull_white
            .coupon_bond_put_t(R_T, 1.0, 1.5, &times, COUPON_RATE, strike)
            .unwrap();
        assert!(
            (call - call_expected).abs() <= 1e-15f64.max(call_expected.abs() * 1e-13),
            "all-post-expiry call changed: fixture {curr},{a},{b},{sigma} strike {strike}: {call} vs {call_expected}"
        );
        assert!(
            (put - put_expected).abs() <= 1e-15f64.max(put_expected.abs() * 1e-13),
            "all-post-expiry put changed: fixture {curr},{a},{b},{sigma} strike {strike}: {put} vs {put_expected}"
        );
    }
}
