extern crate nrfind;
extern crate black_scholes;

#[macro_use]
#[cfg(test)]
extern crate approx;

//nrfind::find_root(&obj_fn, &dfn, initial_guess, precision, iterations).unwrap()
/**
 * Note: the fundamental times 
 * here are (0, t, T, TM).  0 is 
 * current time (and is 
 * reflective of the current yield 
 * curve), t is some future time 
 * that we may want to price 
 * options at given the underlying 
 * at that time, T is an "initial" 
 * maturity and TM a "Final" 
 * maturity.  While it is natural 
 * to think of (0<t<T<TM), I only 
 * require 0<t and 0<T<TM.
 *  Note that ALL TIMES ARE WITH 
 * RESPECT TO 0!  
 * */

const PREC_1:f64=0.0000001;
const PREC_2:f64=0.0000001;
const R_INIT:f64=0.03;
const MAX_ITER:usize=50;

//tdiff=T-t
fn a_t(a:f64, t_diff:f64)->f64{
    (1.0-(-a*t_diff).exp())/a
}

//t is first future time
//T is second future time
fn at_t(a:f64, t:f64, T:f64)->f64{
    a_t(a, T-t)
}

fn ct_t(
    a:f64,
    sigma:f64,
    t:f64,
    T:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64 {
    let sqr=(-a*T).exp()-(-a*t).exp();
    yield_curve(t)-yield_curve(T)
    +forward_curve(t)*at_t(a, t, T)
    -(sigma*sqr).powi(2)
    *((2.0*a*t).exp()-1.0)/(4.0*a.powi(3))
}

fn t_forward_bond_vol(
    a:f64,
    sigma:f64,
    t:f64,
    T:f64,
    T_M:f64
)->f64{
    let exp_d=1.0-(-a*(T_M-T)).exp();
    let exp_t=1.0-(-2.0*a*(T-t)).exp();
    sigma*(exp_t/(2.0*a.powi(3)).sqrt())*exp_d
}

fn phi_t(
    a:f64,
    sigma:f64,
    T:f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let exp_T=1.0-(-a*T).exp();
    forward_curve(T)+(sigma*exp_T).powi(2)/(2.0*a.powi(2))
}

fn mu_r(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    T:f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    phi_t(a, sigma, T, forward_curve)+(r_t-phi_t(
        a, sigma, t, forward_curve
    ))*(-a*(T-t).exp())
}

fn variance_r(
    a:f64,
    sigma:f64,
    t:f64,
    T:f64
)->f64{
    sigma.powi(2)*(1.0-(-2.0*a*(T-t)).exp())/(2.0*a)
}

pub fn bond_price_t(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    bond_maturity:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    (
        -r_t*at_t(a, t, bond_maturity)+
        ct_t(a, sigma, t, bond_maturity, yield_curve, forward_curve)
    ).exp()
}

pub fn bond_price_now(
    bond_maturity:f64,
    yield_curve:&Fn(f64)->f64
)->f64{
    (-yield_curve(bond_maturity)).exp()
}

pub fn coupon_bond_price_t(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    coupon_times:&[f64], //needs last coupon time to be the largest
    coupon_rate:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let last_time=coupon_times.iter().fold(0.0 as f64, |running_max, coupon_time|{
        if *coupon_time>running_max {
            *coupon_time
        }
        else {
            running_max
        }
    });
    coupon_times.iter().fold(
        bond_price_t(r_t, a, sigma, t, last_time, yield_curve, forward_curve),
        |accum, coupon_time|{
            if *coupon_time>t {
                coupon_rate*bond_price_t(
                    r_t, a, sigma, t, *coupon_time, 
                    yield_curve, forward_curve
                )
            }
            else {
                0.0
            }
        }
    )
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
