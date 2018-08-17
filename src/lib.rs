extern crate nrfind;
extern crate black_scholes;
extern crate binomial_tree;

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
const R_INIT:f64=0.03;
const MAX_ITER:i32=50;

//tdiff=T-t
fn a_t(a:f64, t_diff:f64)->f64{
    (1.0-(-a*t_diff).exp())/a
}

//t is first future time
//t_m is second future time
fn at_t(a:f64, t:f64, t_m:f64)->f64{
    a_t(a, t_m-t)
}

fn ct_t(
    a:f64,
    sigma:f64,
    t:f64,
    t_m:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64 {
    let sqr=(-a*t_m).exp()-(-a*t).exp();
    yield_curve(t)-yield_curve(t_m)
    +forward_curve(t)*at_t(a, t, t_m)
    -(sigma*sqr).powi(2)
    *((2.0*a*t).exp()-1.0)/(4.0*a.powi(3))
}

//t_forward_bond_vol(a, sigma, t, option_maturity, bond_maturity) 
fn t_forward_bond_vol(
    a:f64,
    sigma:f64,
    t:f64,
    t_m:f64,
    t_f:f64
)->f64{
    let exp_d=1.0-(-a*(t_f-t_m)).exp();
    let exp_t=1.0-(-2.0*a*(t_m-t)).exp();
    sigma*(exp_t/(2.0*a.powi(3))).sqrt()*exp_d
}

fn phi_t(
    a:f64,
    sigma:f64,
    t:f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let exp_t=1.0-(-a*t).exp();
    forward_curve(t)+(sigma*exp_t).powi(2)/(2.0*a.powi(2))
}
/*
fn mu_r(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    t_m:f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    phi_t(a, sigma, t_m, forward_curve)+(r_t-phi_t(
        a, sigma, t, forward_curve
    ))*(-a*(t_m-t).exp())
}

fn variance_r(
    a:f64,
    sigma:f64,
    t:f64,
    T:f64
)->f64{
    sigma.powi(2)*(1.0-(-2.0*a*(T-t)).exp())/(2.0*a)
}

fn get_max_float(
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
}*/

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

//used for newton's method
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

fn coupon_bond_generic_t(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    coupon_times:&[f64], //does not include the bond_maturity, but the function does check for that
    bond_maturity:f64,
    coupon_rate:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64,
    generic_fn:&Fn(f64, f64, f64, f64, f64, &Fn(f64)->f64, &Fn(f64)->f64)->f64
)->f64{
    let par_value=1.0; //without loss of generality
    let final_payment=generic_fn(
        r_t, a, sigma, t, 
        bond_maturity, 
        yield_curve, forward_curve
    )*(par_value+coupon_rate); //par of one+coupon rate
    coupon_times.iter().fold(
        final_payment,
        |accum, coupon_time|{
            if *coupon_time>t && *coupon_time<bond_maturity { //check for whether includes bond maturity
                accum+coupon_rate*generic_fn(
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
fn coupon_bond_generic_now(
    coupon_times:&[f64], //does not include the bond_maturity, but the function does check for that
    bond_maturity:f64,
    coupon_rate:f64,
    yield_curve:&Fn(f64)->f64,
    generic_fn:&Fn(f64, &Fn(f64)->f64)->f64
)->f64{
    let par_value=1.0; //without loss of generality
    let final_payment=generic_fn(
        bond_maturity, 
        yield_curve
    )*(par_value+coupon_rate); //par of one+coupon rate
    coupon_times.iter().fold(
        final_payment,
        |accum, coupon_time|{
            if *coupon_time<bond_maturity { //check for whether includes bond maturity
                accum+coupon_rate*generic_fn(
                    *coupon_time, 
                    yield_curve
                )
            }
            else {
                accum
            }
        }
    )
}
//par=1
pub fn coupon_bond_price_t(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    coupon_times:&[f64], //does not include the bond_maturity, but the function does check for that
    bond_maturity:f64,
    coupon_rate:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    coupon_bond_generic_t(
        r_t, 
        a, sigma, t, 
        coupon_times,
        bond_maturity,
        coupon_rate,
        yield_curve,
        forward_curve,
        &bond_price_t
    )
}

fn coupon_bond_price_t_deriv(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    coupon_times:&[f64], //does not include the bond_maturity, but the function does check for that
    bond_maturity:f64,
    coupon_rate:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    coupon_bond_generic_t(
        r_t, 
        a, sigma, t, 
        coupon_times,
        bond_maturity,
        coupon_rate,
        yield_curve,
        forward_curve,
        &bond_price_t_deriv
    )
}

pub fn coupon_bond_price_now(
    coupon_times:&[f64], //does not include the bond_maturity, but the function does check for that
    bond_maturity:f64,
    coupon_rate:f64, 
    yield_curve:&Fn(f64)->f64
)->f64{
    coupon_bond_generic_now(
        coupon_times, 
        bond_maturity,
        coupon_rate,
        yield_curve,
        &bond_price_now
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
    let t=0.0;//since "now"
    black_scholes::call_discount(
        bond_price_now(bond_maturity, yield_curve), //underlying
        strike,
        bond_price_now(option_maturity, yield_curve), //discount
        t_forward_bond_vol(a, sigma, t, option_maturity, bond_maturity) //volatility with maturity
    )
}
/*The price of a call option on coupon bond under Hull White...uses jamshidian's trick*/
fn coupon_bond_option_generic_t(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    option_maturity:f64,
    coupon_times:&[f64], //does not include final payment (bond maturity)
    bond_maturity:f64,
    coupon_rate:f64,
    strike:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64,
    generic_fn:&Fn(
        f64, f64, f64, 
        f64, f64, f64, 
        f64, &Fn(f64)->f64, &Fn(f64)->f64
    )->f64
)->f64{ 
    let par_value=1.0;
    let fn_to_optimize=|r|{
        coupon_bond_price_t(
            r, a, sigma, 
            option_maturity, 
            coupon_times, 
            bond_maturity,
            coupon_rate, 
            yield_curve, forward_curve
        )-strike
    };
    let fn_derv=|r|{
        coupon_bond_price_t_deriv(
            r, a, sigma, 
            option_maturity, 
            coupon_times, 
            bond_maturity,
            coupon_rate, 
            yield_curve, forward_curve
        )
    };
    let r_optimal=nrfind::find_root(
        &fn_to_optimize, &fn_derv, 
        R_INIT, PREC_1, MAX_ITER
    ).expect("Requires convergence of optimal r");
    let final_payment=generic_fn(
        r_t, a, sigma, t, option_maturity, 
        bond_maturity, 
        bond_price_t( 
            r_optimal,
            a, sigma, option_maturity,
            bond_maturity,
            yield_curve,
            forward_curve
        ),
        yield_curve, forward_curve
    )*(par_value+coupon_rate); //par of one+coupon rate
    coupon_times.iter().fold(
        final_payment,
        |accum, coupon_time|{
            if *coupon_time>t && *coupon_time<bond_maturity  {
                accum+coupon_rate*generic_fn(
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
            else {
                accum
            }
        }
    )
}

pub fn coupon_bond_call_t(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    option_maturity:f64,
    coupon_times:&[f64], //does not include final payment (bond maturity)
    bond_maturity:f64,
    coupon_rate:f64,
    strike:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    coupon_bond_option_generic_t(
        r_t, a, sigma,
        t, option_maturity,
        coupon_times,
        bond_maturity,
        coupon_rate, strike,
        yield_curve, 
        forward_curve,
        &bond_call_t
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
    coupon_times:&[f64], //these should be greater than the option maturity
    bond_maturity:f64,
    coupon_rate:f64,
    strike:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    coupon_bond_option_generic_t(
        r_t, a, sigma,
        t, option_maturity,
        coupon_times,
        bond_maturity,
        coupon_rate, strike,
        yield_curve, 
        forward_curve,
        &bond_put_t
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

fn get_num_payments(
    t:f64,
    maturity:f64,
    delta:f64
)->f64{
    (maturity-t)/delta+1.0
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
    let num_payments=get_num_payments(swap_initiation, swap_maturity, delta) as usize; //this should be an integer!  remember, T-t is the total swap length
    let denominator_swap=(1..(num_payments+1)).fold(0.0, |accum, curr|{
        accum+bond_price_t(r_t, a, sigma, t, swap_initiation+delta*(curr as f64), yield_curve, forward_curve)
    })*delta;
    (bond_price_t(r_t, a, sigma, t, swap_initiation, yield_curve, forward_curve)-bond_price_t(r_t, a, sigma, t, swap_maturity+delta, yield_curve, forward_curve))/denominator_swap
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
    forward_swap_rate_t(
        r_t, a, sigma,
        t, t, //swap is a forward swap starting at time "0" (t-t=0)
        swap_maturity,
        delta, 
        yield_curve,
        forward_curve
    )
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
    let num_payments_fl=get_num_payments(t, swap_maturity, delta); //this should be an integer!  remember, T-t is the total swap length
    let num_payments=num_payments_fl as usize;
    let first_exchange_date=swap_maturity-(num_payments_fl-1.0)*delta;
    let sm_bond=(1..num_payments).fold(0.0, |accum, curr|{
        accum+bond_price_t(r_t, a, sigma, t, first_exchange_date+delta*(curr as f64), yield_curve, forward_curve)*swap_rate*delta
    });
    bond_price_t(
        r_t, a, sigma, t, 
        first_exchange_date, 
        yield_curve, forward_curve
    )-sm_bond-(1.0+swap_rate*delta)*bond_price_t(
        r_t, a, sigma, t,
        first_exchange_date+delta*num_payments_fl, 
        yield_curve, forward_curve
    )
}

fn get_time_from_t_index(
    index:usize,
    t:f64,
    delta:f64
)->f64{
    t+(index as f64)*delta
}
fn get_coupon_times(
    num_payments:usize,
    t:f64,
    delta:f64
)->Vec<f64>{
    (1..(num_payments)).map(|index|{
        get_time_from_t_index(index, t, delta)
    }).collect()
}

pub fn european_payer_swaption_t(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    swap_tenor:f64,
    option_maturity:f64,
    delta:f64, //tenor of simple yield
    swap_rate:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let num_payments=get_num_payments(t, swap_tenor, delta) as usize;
    let coupon_times=get_coupon_times(num_payments, option_maturity, delta);
    let strike=1.0;
    coupon_bond_put_t(
        r_t, a, sigma, t, 
        option_maturity, 
        &coupon_times,
        get_time_from_t_index(num_payments, option_maturity, delta),
        swap_rate*delta, 
        strike, yield_curve, 
        forward_curve
    ) //swaption is equal to put on coupon bond with coupon=swaption swapRate*delta and strike 1.
}

pub fn european_receiver_swaption_t(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    swap_tenor:f64,
    option_maturity:f64,
    delta:f64, //tenor of simple yield
    swap_rate:f64,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let num_payments=get_num_payments(t, swap_tenor, delta) as usize;
    let coupon_times=get_coupon_times(num_payments, option_maturity, delta);
    let strike=1.0;
    coupon_bond_call_t(
        r_t, a, sigma, t, 
        option_maturity, 
        &coupon_times, 
        get_time_from_t_index(num_payments, option_maturity, delta),
        swap_rate*delta, 
        strike, yield_curve, 
        forward_curve
    ) //swaption is equal to call on coupon bond with coupon=swapRate*delta and strike 1.
}


fn max_or_zero(
    v:f64
)->f64{
    if v>0.0{
        v
    }
    else{
        0.0
    }
}
fn payoff_swaption(
    is_payer:bool,
    swp:f64
)->f64{
    match is_payer {
        true=>max_or_zero(swp),
        false=>max_or_zero(-swp)
    }
}

fn american_swaption(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    swap_tenor:f64,
    option_maturity:f64,
    delta:f64, //tenor of simple yield
    swap_rate:f64,
    is_payer:bool,
    num_steps:usize,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let alpha_div_sigma=|_t_step:f64, curr_val:f64, _dt:f64, _width:usize| -(a*curr_val)/sigma;
    let sigma_prime=|_t_step:f64, _curr_val:f64, _dt:f64, _j:usize| 0.0;
    let sigma_inv=|_t_step:f64, y:f64, _dt:f64, _j:usize| sigma*y;
    let t_of_option=option_maturity-t;
    let mut phi_cache:Vec<f64>=binomial_tree::get_all_t(t_of_option, num_steps).map(|t_a|phi_t(a, sigma, t_a, forward_curve)).collect();
    phi_cache.push(
        phi_t(a, sigma, t_of_option, forward_curve)
    );
    let payoff=|t_step:f64, curr_val:f64, _dt:f64, j:usize| {      
        let swp=swap_price_t(
            curr_val+phi_cache[j], 
            a, sigma, t_step, 
            swap_tenor+t_step, delta, 
            swap_rate, yield_curve, 
            forward_curve
        );
        payoff_swaption(is_payer, swp)

    };
    let discount=|_t_step:f64, curr_val:f64, dt:f64, j:usize| {
        (-(curr_val+phi_cache[j])*dt).exp()
    };
    binomial_tree::compute_price_american(
        &alpha_div_sigma,
        &sigma_prime,
        &sigma_inv,
        &payoff,
        &discount,
        (r_t-phi_t(a, sigma, t, forward_curve))/sigma, //initial "y"
        t_of_option,
        num_steps
    )
}

pub fn american_payer_swaption_t(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    swap_tenor:f64,
    option_maturity:f64,
    delta:f64, //tenor of simple yield
    swap_rate:f64,
    num_steps:usize,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let is_payer=true;
    american_swaption(
        r_t, a, sigma, t, swap_tenor,
        option_maturity, delta, 
        swap_rate, is_payer, num_steps, 
        yield_curve, forward_curve
    )
}

pub fn american_receiver_swaption_t(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    swap_tenor:f64,
    option_maturity:f64,
    delta:f64, //tenor of simple yield
    swap_rate:f64,
    num_steps:usize,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let is_payer=false;
    american_swaption(
        r_t, a, sigma, t, swap_tenor,
        option_maturity, delta, 
        swap_rate, is_payer, num_steps, 
        yield_curve, forward_curve
    )
}
#[cfg(test)]
fn european_swaption_tree(
    r_t:f64,
    a:f64,
    sigma:f64,
    t:f64,
    swap_tenor:f64,
    option_maturity:f64,
    delta:f64, //tenor of simple yield
    swap_rate:f64,
    is_payer:bool,
    num_steps:usize,
    yield_curve:&Fn(f64)->f64,
    forward_curve:&Fn(f64)->f64
)->f64{
    let alpha_div_sigma=|_t_step:f64, curr_val:f64, _dt:f64, _width:usize| -(a*curr_val)/sigma;
    let sigma_prime=|_t_step:f64, _curr_val:f64, _dt:f64, _j:usize| 0.0;
    let sigma_inv=|_t_step:f64, y:f64, _dt:f64, _j:usize| sigma*y;
    let mut phi_cache:Vec<f64>=binomial_tree::get_all_t(
        option_maturity-t, 
        num_steps
    ).map(|t_a|phi_t(a, sigma, t_a, forward_curve)).collect();
    phi_cache.push(
        phi_t(a, sigma, option_maturity-t, forward_curve)
    );
    let payoff=|t_step:f64, curr_val:f64, _dt:f64, j:usize| {
        let swp=swap_price_t(
            curr_val+phi_cache[j], 
            a, sigma, t_step, 
            swap_tenor+t_step, delta, 
            swap_rate, yield_curve, 
            forward_curve
        );
        payoff_swaption(is_payer, swp)
    };
    let discount=|_t_a:f64, curr_val:f64, dt:f64, j:usize| {
        (-(curr_val+phi_cache[j])*dt).exp()
    };
    binomial_tree::compute_price_raw(
        &alpha_div_sigma,
        &sigma_prime,
        &sigma_inv,
        &payoff,
        &discount,
        (r_t-phi_t(a, sigma, t, forward_curve))/sigma, //initial "y"
        option_maturity-t,
        num_steps,
        false
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_bond_now_same_as_t_when_t_is_zero(){
        let curr_rate=0.02;
        let sig:f64=0.02;
        let a:f64=0.3;
        let b=0.04;
        let future_time=0.0;
        let maturity=1.5;
        let yield_curve=|t:f64|{
            let at=(1.0-(-a*t).exp())/a;
            let ct=(b-sig.powi(2)/(2.0*a.powi(2)))*(at-t)-(sig*at).powi(2)/(4.0*a);
            at*curr_rate-ct
        };
        let forward_curve=|t:f64|{
            b+(-a*t).exp()*(curr_rate-b)-(sig.powi(2)/(2.0*a.powi(2)))*(1.0-(-a*t).exp()).powi(2)
        };
        let bond_price_now=bond_price_now(
            maturity, 
            &yield_curve
        );
        let bond_price_t=bond_price_t(
            curr_rate, a, sig, future_time, 
            maturity,
            &yield_curve, &forward_curve
        );
        assert_abs_diff_eq!(
            bond_price_now, 
            bond_price_t, 
            epsilon=0.0000001
        );
    }
    #[test]
    fn test_coupon_bond_now_same_as_t_when_t_is_zero(){
        let curr_rate=0.02;
        let sig:f64=0.02;
        let a:f64=0.3;
        let b=0.04;
        let delta=0.25;
        let future_time=0.0;
        let yield_curve=|t:f64|{
            let at=(1.0-(-a*t).exp())/a;
            let ct=(b-sig.powi(2)/(2.0*a.powi(2)))*(at-t)-(sig*at).powi(2)/(4.0*a);
            at*curr_rate-ct
        };
        let forward_curve=|t:f64|{
            b+(-a*t).exp()*(curr_rate-b)-(sig.powi(2)/(2.0*a.powi(2)))*(1.0-(-a*t).exp()).powi(2)
        };
        let coupon_times=get_coupon_times(5, future_time, delta);
        let coupon_rate=0.05*delta;
        let bond_price_now=coupon_bond_price_now(
            &coupon_times, 
            *coupon_times.last().unwrap()+delta, 
            coupon_rate, 
            &yield_curve
        );
        let bond_price_t=coupon_bond_price_t(
            curr_rate, a, sig, future_time, 
            &coupon_times, 
            *coupon_times.last().unwrap()+delta, 
            coupon_rate, 
            &yield_curve, &forward_curve
        );
        assert_abs_diff_eq!(
            bond_price_now, 
            bond_price_t, 
            epsilon=0.0000001
        );
    }
   
    #[test]
    fn test_bond_price() {
        let curr_rate=0.02;
        let sig:f64=0.02;
        let a:f64=0.3;
        let b=0.04;
        let future_time=0.5;
        let option_maturity=1.5;
        let yield_curve=|t:f64|{
            let at=(1.0-(-a*t).exp())/a;
            let ct=(b-sig.powi(2)/(2.0*a.powi(2)))*(at-t)-(sig*at).powi(2)/(4.0*a);
            at*curr_rate-ct
        };
        let forward_curve=|t:f64|{
            b+(-a*t).exp()*(curr_rate-b)-(sig.powi(2)/(2.0*a.powi(2)))*(1.0-(-a*t).exp()).powi(2)
        };
        assert_eq!(
            bond_price_t(curr_rate, a, sig, future_time, option_maturity, &yield_curve, &forward_curve), 
            bond_price_now(option_maturity-future_time, &yield_curve)
        );
    }
    #[test]
    fn test_bond_price_at_expiry() {
        let curr_rate=0.02;
        let sig:f64=0.02;
        let a:f64=0.3;
        let b=0.04;
        let future_time=0.5;
        let yield_curve=|t:f64|{
            let at=(1.0-(-a*t).exp())/a;
            let ct=(b-sig.powi(2)/(2.0*a.powi(2)))*(at-t)-(sig*at).powi(2)/(4.0*a);
            at*curr_rate-ct
        };
        let forward_curve=|t:f64|{
            b+(-a*t).exp()*(curr_rate-b)-(sig.powi(2)/(2.0*a.powi(2)))*(1.0-(-a*t).exp()).powi(2)
        };
        assert_eq!(
            bond_price_t(curr_rate, a, sig, future_time, future_time, &yield_curve, &forward_curve), 
            1.0
        );
    }
    #[test]
    fn test_swap(){
        let curr_rate=0.02;
        let sig:f64=0.02;
        let a:f64=0.3;
        let b=0.04;
        let delta=0.25;
        let future_time=0.5;
        let swap_maturity=5.5;
        let yield_curve=|t:f64|{
            let at=(1.0-(-a*t).exp())/a;
            let ct=(b-sig.powi(2)/(2.0*a.powi(2)))*(at-t)-(sig*at).powi(2)/(4.0*a);
            at*curr_rate-ct
        };
        let forward_curve=|t:f64|{
            b+(-a*t).exp()*(curr_rate-b)-(sig.powi(2)/(2.0*a.powi(2)))*(1.0-(-a*t).exp()).powi(2)
        };

        assert_abs_diff_eq!(
            swap_price_t(
                curr_rate, a, sig, future_time, 
                swap_maturity, delta, 
                swap_rate_t(
                    curr_rate, a, sig, future_time, 
                    swap_maturity, 
                    delta, 
                    &yield_curve, &forward_curve
                ),
                &yield_curve, &forward_curve
            ), 
            0.0,
            epsilon=0.000000001
        );
    }
    #[test]
    fn payer_swaption(){
        let curr_rate=0.05;
        let sig:f64=0.01;
        let a:f64=0.05;
        let b=0.05;
        let delta=0.25;
        let future_time=0.0;
        let swap_tenor=5.0;
        let option_maturity=1.0;
        let yield_curve=|t:f64|{
            let at=(1.0-(-a*t).exp())/a;
            let ct=(b-sig.powi(2)/(2.0*a.powi(2)))*(at-t)-(sig*at).powi(2)/(4.0*a);
            at*curr_rate-ct
        };
        let forward_curve=|t:f64|{
            b+(-a*t).exp()*(curr_rate-b)-(sig.powi(2)/(2.0*a.powi(2)))*(1.0-(-a*t).exp()).powi(2)
        };
        let swap_rate=forward_swap_rate_t(
            curr_rate, a, sig, 
            future_time, option_maturity, 
            swap_tenor+option_maturity, 
            delta, &yield_curve, 
            &forward_curve
        );
        let analytical=european_payer_swaption_t(
            curr_rate, a, sig, 
            future_time, swap_tenor, 
            option_maturity, delta, swap_rate,
            &yield_curve, &forward_curve
        );
        let is_payer=true;

        let tree=european_swaption_tree(
            curr_rate, a, sig, 
            future_time, swap_tenor, 
            option_maturity, delta, swap_rate,
            is_payer, 5000,
            &yield_curve, &forward_curve
        );
        assert_abs_diff_eq!(
            analytical,
            tree,
            epsilon=0.000001
        )
    }
    #[test]
    fn receiver_swaption(){
        let curr_rate=0.05;
        let sig:f64=0.01;
        let a:f64=0.05;
        let b=0.05;
        let delta=0.25;
        let future_time=0.0;
        let swap_tenor=5.0;
        let option_maturity=1.0;
        let yield_curve=|t:f64|{
            let at=(1.0-(-a*t).exp())/a;
            let ct=(b-sig.powi(2)/(2.0*a.powi(2)))*(at-t)-(sig*at).powi(2)/(4.0*a);
            at*curr_rate-ct
        };
        let forward_curve=|t:f64|{
            b+(-a*t).exp()*(curr_rate-b)-(sig.powi(2)/(2.0*a.powi(2)))*(1.0-(-a*t).exp()).powi(2)
        };
        let swap_rate=forward_swap_rate_t(
            curr_rate, a, sig, 
            future_time, option_maturity, 
            swap_tenor+option_maturity, 
            delta, &yield_curve, 
            &forward_curve
        );
        let analytical=european_receiver_swaption_t(
            curr_rate, a, sig, 
            future_time, swap_tenor, 
            option_maturity, delta, swap_rate,
            &yield_curve, &forward_curve
        );
        let is_payer=false;

        let tree=european_swaption_tree(
            curr_rate, a, sig, 
            future_time, swap_tenor, 
            option_maturity, delta, swap_rate,
            is_payer, 5000,
            &yield_curve, &forward_curve
        );
        assert_abs_diff_eq!(
            analytical,
            tree,
            epsilon=0.000001
        )
    }
    #[test]
    fn zero_coupon_reference(){ //http://www.quantcalc.net/BondOption_Vasicek.html
        let curr_rate=0.01;
        let sig:f64=0.03;
        let a=0.05;
        let b=0.04;
        let strike=0.96;
        let future_time=0.0;
        let bond_maturity=3.0;
        let option_maturity=2.0;
        let yield_curve=|t:f64|{
            let at=(1.0-(-a*t).exp())/a;
            let ct=(b-sig.powi(2)/(2.0*a.powi(2)))*(at-t)-(sig*at).powi(2)/(4.0*a);
            at*curr_rate-ct
        };
        let forward_curve=|t:f64|{
            b+(-a*t).exp()*(curr_rate-b)-(sig.powi(2)/(2.0*a.powi(2)))*(1.0-(-a*t).exp()).powi(2)
        };
        let bond_call=bond_call_t(
            curr_rate, a, sig, future_time, option_maturity, 
            bond_maturity, strike, 
            &yield_curve, &forward_curve
        );
        assert_abs_diff_eq!(
            bond_call,
            0.033282,
            epsilon=0.0001
        )

    }
    #[test]
    fn zero_coupon_to_coupon(){ //http://www.quantcalc.net/BondOption_Vasicek.html
        let curr_rate=0.01;
        let sig:f64=0.03;
        let a=0.05;
        let b=0.04;
        let strike=0.96;
        let future_time=0.0;
        let bond_maturity=3.0;
        let option_maturity=2.0;
        let yield_curve=|t:f64|{
            let at=(1.0-(-a*t).exp())/a;
            let ct=(b-sig.powi(2)/(2.0*a.powi(2)))*(at-t)-(sig*at).powi(2)/(4.0*a);
            at*curr_rate-ct
        };
        let forward_curve=|t:f64|{
            b+(-a*t).exp()*(curr_rate-b)-(sig.powi(2)/(2.0*a.powi(2)))*(1.0-(-a*t).exp()).powi(2)
        };
        let bond_call=bond_call_t(
            curr_rate, a, sig, future_time, option_maturity, 
            bond_maturity, strike, 
            &yield_curve, &forward_curve
        );
        let coupon_rate=0.0;
        let coupon_bond_call=coupon_bond_call_t(
            curr_rate, a, sig, future_time, option_maturity, 
            &[], bond_maturity, 
            coupon_rate, strike, 
            &yield_curve, &forward_curve
        );
        
        assert_abs_diff_eq!(
            bond_call,
            coupon_bond_call,
            epsilon=0.0001
        )

    }
}
