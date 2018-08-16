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
const MAX_ITER:i32=50;

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

fn get_max_double(
    arr:&[f64]
)->f64{
    arr.iter().fold(0.0 as f64, |running_max, elem|{
        if *elem>running_max {
            *elem
        }
        else {
            running_max
        }
    })
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

fn bond_price_t_deriv(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    bond_maturity:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let at_t_c=at_t(a, t, bond_maturity);
    -(
        -r_t*at_t_c+
        ct_t(a, sigma, t, bond_maturity, yield_curve, forward_curve)
    ).exp()*at_t_c
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
    coupon_times:&[f64], 
    coupon_rate:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let last_time=get_max_double(coupon_times);
    coupon_times.iter().fold(
        bond_price_t(r_t, a, sigma, t, last_time, yield_curve, forward_curve),
        |accum, coupon_time|{
            if *coupon_time>t {
                accum+coupon_rate*bond_price_t(
                    r_t, a, sigma, t, *coupon_time, 
                    yield_curve, forward_curve
                )
            }
            else {
                accum
            }
        }
    )
}

fn coupon_bond_price_t_deriv(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    coupon_times:&[f64], 
    coupon_rate:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let last_time=get_max_double(coupon_times);
    coupon_times.iter().fold(
        bond_price_t_deriv(r_t, a, sigma, t, last_time, yield_curve, forward_curve),
        |accum, coupon_time|{
            if *coupon_time>t {
                accum+coupon_rate*bond_price_t_deriv(
                    r_t, a, sigma, t, *coupon_time, 
                    yield_curve, forward_curve
                )
            }
            else {
                accum
            }
        }
    )
}

pub fn coupon_bond_price_now(
    coupon_times:&[f64],
    coupon_rate:f64,
    yield_curve:&Fn(f64)->f64
)->f64{
    let last_time=get_max_double(coupon_times);
    coupon_times.iter().fold(
        bond_price_now(last_time, yield_curve),
        |accum, coupon_time|{
            accum+coupon_rate*bond_price_now(
                *coupon_time, 
                yield_curve
            )
        }
    )
}

pub fn bond_call_t(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    option_maturity:f64,
    bond_maturity:f64,
    strike:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    black_scholes::call_discount(
        bond_price_t(r_t, a, sigma, t, bond_maturity, yield_curve, forward_curve), //underlying
        strike,
        bond_price_t(r_t, a, sigma, t, option_maturity, yield_curve, forward_curve), //discount
        t_forward_bond_vol(a, sigma, t, option_maturity, bond_maturity) //volatility with maturity
    )
}

pub fn bond_call_now(
    a:f64,
    sigma:f64,
    option_maturity:f64,
    bond_maturity:f64,
    strike:f64,
    yield_curve:&Fn(f64)->f64
)->f64{
    black_scholes::call_discount(
        bond_price_now(bond_maturity, yield_curve), //underlying
        strike,
        bond_price_now(option_maturity, yield_curve), //discount
        t_forward_bond_vol(a, sigma, 0.0, option_maturity, bond_maturity) //volatility with maturity
    )
}

pub fn coupon_bond_call_t(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    option_maturity:f64,
    coupon_times:&[f64],
    coupon_rate:f64,
    strike:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let last_time=get_max_double(coupon_times);
    let fn_to_optimize=|r|{
        coupon_bond_price_t(r, a, sigma, t, coupon_times, coupon_rate, yield_curve, forward_curve)-strike
    };
    let fn_derv=|r|{
        coupon_bond_price_t_deriv(r, a, sigma, t, coupon_times, coupon_rate, yield_curve, forward_curve)
    };
    let r_optimal=nrfind::find_root(&fn_to_optimize, &fn_derv, R_INIT, PREC_1, MAX_ITER).expect("Requires convergence of optimal r");
    coupon_times.iter().fold(
        bond_call_t(
            r_t, a, sigma, t, option_maturity, 
            last_time, 
            bond_price_t(
                r_optimal,
                a, sigma, option_maturity,
                last_time,
                yield_curve,
                forward_curve
            ),
            yield_curve, forward_curve
        ),
        |accum, coupon_time|{
            accum+coupon_rate*bond_call_t(
                r_t, a, sigma, t, option_maturity, 
                *coupon_time, 
                bond_price_t(
                    r_optimal,
                    a, sigma, option_maturity,
                    *coupon_time,
                    yield_curve,
                    forward_curve
                ),
                yield_curve, forward_curve
            )
        }
    )
}

pub fn bond_put_t(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    option_maturity:f64,
    bond_maturity:f64,
    strike:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    black_scholes::put_discount(
        bond_price_t(r_t, a, sigma, t, bond_maturity, yield_curve, forward_curve), //underlying
        strike,
        bond_price_t(r_t, a, sigma, t, option_maturity, yield_curve, forward_curve), //discount
        t_forward_bond_vol(a, sigma, t, option_maturity, bond_maturity) //volatility with maturity
    )
}

pub fn bond_put_now(
    a:f64,
    sigma:f64,
    option_maturity:f64,
    bond_maturity:f64,
    strike:f64,
    yield_curve:&Fn(f64)->f64
)->f64{
    black_scholes::put_discount(
        bond_price_now(bond_maturity, yield_curve), //underlying
        strike,
        bond_price_now(option_maturity, yield_curve), //discount
        t_forward_bond_vol(a, sigma, 0.0, option_maturity, bond_maturity) //volatility with maturity
    )
}

pub fn coupon_bond_put_t(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    option_maturity:f64,
    coupon_times:&[f64],
    coupon_rate:f64,
    strike:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let last_time=get_max_double(coupon_times);
    let fn_to_optimize=|r|{
        coupon_bond_price_t(r, a, sigma, t, coupon_times, coupon_rate, yield_curve, forward_curve)-strike
    };
    let fn_derv=|r|{
        coupon_bond_price_t_deriv(r, a, sigma, t, coupon_times, coupon_rate, yield_curve, forward_curve)
    };
    let r_optimal=nrfind::find_root(&fn_to_optimize, &fn_derv, R_INIT, PREC_1, MAX_ITER).expect("Requires convergence of optimal r");
    coupon_times.iter().fold(
        bond_put_t(
            r_t, a, sigma, t, option_maturity, 
            last_time, 
            bond_price_t(
                r_optimal,
                a, sigma, option_maturity,
                last_time,
                yield_curve,
                forward_curve
            ),
            yield_curve, forward_curve
        ),
        |accum, coupon_time|{
            coupon_rate*bond_put_t(
                r_t, a, sigma, t, option_maturity, 
                *coupon_time, 
                bond_price_t(
                    r_optimal,
                    a, sigma, option_maturity,
                    *coupon_time,
                    yield_curve,
                    forward_curve
                ),
                yield_curve, forward_curve
            )
        }
    )
}

//equivalent to an option on a bond
//t=0
pub fn caplet_now(
    a:f64,
    sigma:f64,
    option_maturity:f64,
    delta:f64,//tenor of the simple yield
    strike:f64,
    yield_curve:&Fn(f64)->f64
)->f64{
    (strike*delta+1.0)*bond_put_now(
        a, sigma, option_maturity, 
        option_maturity+delta, 
        1.0/(delta*strike+1.0), 
        yield_curve
    )
}

pub fn euro_dollar_future(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    option_maturity:f64,
    delta:f64, //tenor of simple yield
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let exp_t=(-a*(option_maturity-t)).exp();
    let exp_d=(-a*delta).exp();
    let gamma=(sigma.powi(2)/a.powi(3))*(1.0-exp_d)*((1.0-exp_t)-exp_d*0.5*(1.0-exp_t.powi(2)));
    (
        (
            bond_price_t(
                r_t,
                a, sigma,
                t, option_maturity,
                yield_curve, 
                forward_curve
            )/bond_price_t(
                r_t,
                a, sigma,
                t, option_maturity+delta,
                yield_curve, 
                forward_curve
            )
        )*gamma.exp()-1.0
    )/delta
}

pub fn swap_rate_t(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    swap_maturity:f64,
    delta:f64, //tenor of simple yield
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let div_swap_maturity=(swap_maturity-t)/delta;
    let num_payments=div_swap_maturity.floor() as usize+1; //this should be an integer!  remember, T-t is the total swap length
    let denominator_swap=(1..(num_payments+1)).fold(0.0, |accum, curr|{
        accum+bond_price_t(r_t, a, sigma, t, t+delta*(curr as f64), yield_curve, forward_curve)
    });
    (1.0-bond_price_t(r_t, a, sigma, t, swap_maturity+delta, yield_curve, forward_curve))/denominator_swap
} 

pub fn forward_swap_rate_t(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    swap_initiation:f64,
    swap_maturity:f64,
    delta:f64, //tenor of simple yield
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let div_swap_maturity=(swap_maturity-t)/delta;
    let num_payments=div_swap_maturity.floor() as usize+1; //this should be an integer!  remember, T-t is the total swap length
    let denominator_swap=(1..(num_payments+1)).fold(0.0, |accum, curr|{
        accum+bond_price_t(r_t, a, sigma, t, t+delta*(curr as f64), yield_curve, forward_curve)
    });
    (bond_price_t(r_t, a, sigma, t, swap_initiation, yield_curve, forward_curve)-bond_price_t(r_t, a, sigma, t, swap_maturity+delta, yield_curve, forward_curve))/denominator_swap
} 

pub fn swap_price_t(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    swap_maturity:f64,
    delta:f64, //tenor of simple yield
    swap_rate:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let div_swap_maturity=(swap_maturity-t)/delta;
    let num_payments_fl=div_swap_maturity.floor()+1.0; //this should be an integer!  remember, T-t is the total swap length
    let num_payments=num_payments_fl as usize;
    let first_exchange_date=swap_maturity-(num_payments_fl-1.0)*delta;
    let sm_bond=(1..num_payments).fold(0.0, |accum, curr|{
        accum+bond_price_t(r_t, a, sigma, t, first_exchange_date+delta*(curr as f64), yield_curve, forward_curve)*swap_rate*delta
    });
    bond_price_t(r_t, a, sigma, t, first_exchange_date, yield_curve, forward_curve)-sm_bond-(1.0+swap_rate*delta)*bond_price_t(r_t, a, sigma, t, first_exchange_date+delta*num_payments_fl, yield_curve, forward_curve)

}


#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
