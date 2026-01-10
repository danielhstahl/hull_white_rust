//! Financial products in the Hull-White model

use crate::model::HullWhite;
use crate::Result;

/// Trait for pricing bonds
pub trait BondPricer<T, U> 
where
    T: Fn(f64) -> f64 + Sync,
    U: Fn(f64) -> f64 + Sync,
{
    /// Price a zero-coupon bond at time t
    fn price_zcb_t(&self, r_t: f64, t: f64, maturity: f64) -> f64;
    
    /// Price a zero-coupon bond at time 0
    fn price_zcb_now(&self, maturity: f64) -> f64;
    
    /// Price a coupon bond at time t
    fn price_coupon_bond_t(&self, r_t: f64, t: f64, coupon_times: &[f64], coupon_rate: f64) -> f64;
    
    /// Price a coupon bond at time 0
    fn price_coupon_bond_now(&self, coupon_times: &[f64], coupon_rate: f64) -> f64;
}

impl<'a, T, U> BondPricer<T, U> for HullWhite<'a, T, U>
where
    T: Fn(f64) -> f64 + Sync,
    U: Fn(f64) -> f64 + Sync,
{
    fn price_zcb_t(&self, r_t: f64, t: f64, maturity: f64) -> f64 {
        self.bond_price_t(r_t, t, maturity)
    }
    
    fn price_zcb_now(&self, maturity: f64) -> f64 {
        self.bond_price_now(maturity)
    }
    
    fn price_coupon_bond_t(&self, r_t: f64, t: f64, coupon_times: &[f64], coupon_rate: f64) -> f64 {
        self.coupon_bond_price_t(r_t, t, coupon_times, coupon_rate)
    }
    
    fn price_coupon_bond_now(&self, coupon_times: &[f64], coupon_rate: f64) -> f64 {
        self.coupon_bond_price_now(coupon_times, coupon_rate)
    }
}

/// Trait for pricing bond options
pub trait BondOptionPricer<T, U>
where
    T: Fn(f64) -> f64 + Sync,
    U: Fn(f64) -> f64 + Sync,
{
    /// Price a call option on a zero-coupon bond at time t
    fn price_zcb_call_t(&self, r_t: f64, t: f64, option_maturity: f64, bond_maturity: f64, strike: f64) -> f64;

    /// Price a call option on a zero-coupon bond at time 0
    fn price_zcb_call_now(&self, option_maturity: f64, bond_maturity: f64, strike: f64) -> f64;

    /// Price a put option on a zero-coupon bond at time t
    fn price_zcb_put_t(&self, r_t: f64, t: f64, option_maturity: f64, bond_maturity: f64, strike: f64) -> f64;

    /// Price a put option on a zero-coupon bond at time 0
    fn price_zcb_put_now(&self, option_maturity: f64, bond_maturity: f64, strike: f64) -> f64;

    /// Price a call option on a coupon bond at time t
    fn price_coupon_bond_call_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
    ) -> Result<f64>;

    /// Price a put option on a coupon bond at time t
    fn price_coupon_bond_put_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
    ) -> Result<f64>;
}

impl<'a, T, U> BondOptionPricer<T, U> for HullWhite<'a, T, U>
where
    T: Fn(f64) -> f64 + Sync,
    U: Fn(f64) -> f64 + Sync,
{
    fn price_zcb_call_t(&self, r_t: f64, t: f64, option_maturity: f64, bond_maturity: f64, strike: f64) -> f64 {
        self.bond_call_t(r_t, t, option_maturity, bond_maturity, strike)
    }

    fn price_zcb_call_now(&self, option_maturity: f64, bond_maturity: f64, strike: f64) -> f64 {
        self.bond_call_now(option_maturity, bond_maturity, strike)
    }

    fn price_zcb_put_t(&self, r_t: f64, t: f64, option_maturity: f64, bond_maturity: f64, strike: f64) -> f64 {
        self.bond_put_t(r_t, t, option_maturity, bond_maturity, strike)
    }

    fn price_zcb_put_now(&self, option_maturity: f64, bond_maturity: f64, strike: f64) -> f64 {
        self.bond_put_now(option_maturity, bond_maturity, strike)
    }

    fn price_coupon_bond_call_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
    ) -> Result<f64> {
        self.coupon_bond_call_t(r_t, t, option_maturity, coupon_times, coupon_rate, strike)
    }

    fn price_coupon_bond_put_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        coupon_times: &[f64],
        coupon_rate: f64,
        strike: f64,
    ) -> Result<f64> {
        self.coupon_bond_put_t(r_t, t, option_maturity, coupon_times, coupon_rate, strike)
    }
}

/// Trait for pricing interest rate derivatives
pub trait InterestRateDerivativePricer<T, U>
where
    T: Fn(f64) -> f64 + Sync,
    U: Fn(f64) -> f64 + Sync,
{
    /// Price a caplet at time 0
    fn price_caplet_now(&self, option_maturity: f64, delta: f64, strike: f64) -> f64;
    
    /// Price a caplet at time t
    fn price_caplet_t(&self, r_t: f64, t: f64, option_maturity: f64, delta: f64, strike: f64) -> f64;
    
    /// Price a Eurodollar futures contract at time t
    fn price_edf_t(&self, r_t: f64, t: f64, option_maturity: f64, delta: f64) -> f64;
    
    /// Price a Eurodollar futures contract at time 0
    fn price_edf_now(&self, option_maturity: f64, delta: f64) -> f64;
    
    /// Calculate forward LIBOR rate at time t
    fn forward_libor_rate_t(&self, r_t: f64, t: f64, maturity: f64, delta: f64) -> f64;
    
    /// Calculate forward LIBOR rate at time 0
    fn forward_libor_rate_now(&self, maturity: f64, delta: f64) -> f64;
    
    /// Calculate LIBOR rate at time t
    fn libor_rate_t(&self, r_t: f64, t: f64, delta: f64) -> f64;
}

impl<'a, T, U> InterestRateDerivativePricer<T, U> for HullWhite<'a, T, U>
where
    T: Fn(f64) -> f64 + Sync,
    U: Fn(f64) -> f64 + Sync,
{
    fn price_caplet_now(&self, option_maturity: f64, delta: f64, strike: f64) -> f64 {
        self.caplet_now(option_maturity, delta, strike)
    }
    
    fn price_caplet_t(&self, r_t: f64, t: f64, option_maturity: f64, delta: f64, strike: f64) -> f64 {
        self.caplet_t(r_t, t, option_maturity, delta, strike)
    }
    
    fn price_edf_t(&self, r_t: f64, t: f64, option_maturity: f64, delta: f64) -> f64 {
        self.euro_dollar_future_t(r_t, t, option_maturity, delta)
    }
    
    fn price_edf_now(&self, option_maturity: f64, delta: f64) -> f64 {
        self.euro_dollar_future_now(option_maturity, delta)
    }
    
    fn forward_libor_rate_t(&self, r_t: f64, t: f64, maturity: f64, delta: f64) -> f64 {
        self.forward_libor_rate_t(r_t, t, maturity, delta)
    }
    
    fn forward_libor_rate_now(&self, maturity: f64, delta: f64) -> f64 {
        self.forward_libor_rate_now(maturity, delta)
    }
    
    fn libor_rate_t(&self, r_t: f64, t: f64, delta: f64) -> f64 {
        self.libor_rate_t(r_t, t, delta)
    }
}

/// Trait for pricing swaps and swaptions
pub trait SwapPricer<T, U>
where
    T: Fn(f64) -> f64 + Sync,
    U: Fn(f64) -> f64 + Sync,
{
    /// Calculate forward swap rate at time t
    fn forward_swap_rate_t(
        &self,
        r_t: f64,
        t: f64,
        swap_initiation: f64,
        num_swap_payments: usize,
        delta: f64,
    ) -> f64;
    
    /// Calculate swap rate at time t
    fn swap_rate_t(&self, r_t: f64, t: f64, num_swap_payments: usize, delta: f64) -> f64;
    
    /// Price a swap at time t
    fn price_swap_t(
        &self,
        r_t: f64,
        t: f64,
        swap_maturity: f64,
        delta: f64,
        swap_rate: f64,
    ) -> f64;
    
    /// Price a swap at the start of the swap
    fn price_swap_t_init(
        &self,
        r_t: f64,
        t: f64,
        swap_start: f64,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
    ) -> f64;
    
    /// Price a European payer swaption at time t
    fn price_european_payer_swaption_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
    ) -> Result<f64>;
    
    /// Price a European receiver swaption at time t
    fn price_european_receiver_swaption_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
    ) -> Result<f64>;
}

impl<'a, T, U> SwapPricer<T, U> for HullWhite<'a, T, U>
where
    T: Fn(f64) -> f64 + Sync,
    U: Fn(f64) -> f64 + Sync,
{
    fn forward_swap_rate_t(
        &self,
        r_t: f64,
        t: f64,
        swap_initiation: f64,
        num_swap_payments: usize,
        delta: f64,
    ) -> f64 {
        self.forward_swap_rate_t(r_t, t, swap_initiation, num_swap_payments, delta)
    }
    
    fn swap_rate_t(&self, r_t: f64, t: f64, num_swap_payments: usize, delta: f64) -> f64 {
        self.swap_rate_t(r_t, t, num_swap_payments, delta)
    }
    
    fn price_swap_t(
        &self,
        r_t: f64,
        t: f64,
        swap_maturity: f64,
        delta: f64,
        swap_rate: f64,
    ) -> f64 {
        self.swap_price_t(r_t, t, swap_maturity, delta, swap_rate)
    }
    
    fn price_swap_t_init(
        &self,
        r_t: f64,
        t: f64,
        swap_start: f64,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
    ) -> f64 {
        self.swap_price_t_init(r_t, t, swap_start, num_swap_payments, delta, swap_rate)
    }
    
    fn price_european_payer_swaption_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
    ) -> Result<f64> {
        self.european_payer_swaption_t(r_t, t, option_maturity, num_swap_payments, delta, swap_rate)
    }
    
    fn price_european_receiver_swaption_t(
        &self,
        r_t: f64,
        t: f64,
        option_maturity: f64,
        num_swap_payments: usize,
        delta: f64,
        swap_rate: f64,
    ) -> Result<f64> {
        self.european_receiver_swaption_t(r_t, t, option_maturity, num_swap_payments, delta, swap_rate)
    }
}