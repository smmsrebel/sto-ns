#![cfg_attr(not(feature = "std"), no_std)]

mod math {
    #[cfg(not(feature = "std"))]
    pub use core::f64::consts::PI;
    #[cfg(feature = "std")]
    pub use std::f64::consts::PI;

    #[cfg(feature = "std")]
    #[inline(always)]
    pub fn exp(x: f64) -> f64 {
        x.exp()
    }
    #[cfg(not(feature = "std"))]
    #[inline(always)]
    pub fn exp(x: f64) -> f64 {
        libm::exp(x)
    }

    #[cfg(feature = "std")]
    #[inline(always)]
    pub fn abs(x: f64) -> f64 {
        x.abs()
    }
    #[cfg(not(feature = "std"))]
    #[inline(always)]
    pub fn abs(x: f64) -> f64 {
        libm::fabs(x)
    }

    #[cfg(feature = "std")]
    #[inline(always)]
    pub fn powi(x: f64, n: i32) -> f64 {
        x.powi(n)
    }
    #[cfg(not(feature = "std"))]
    #[inline(always)]
    pub fn powi(x: f64, n: i32) -> f64 {
        libm::pow(x, n as f64)
    }
}

use math::*;

const SINGULARITY_THRESHOLD: f64 = 1e-7;
const TAYLOR_THRESHOLD: f64 = 1e-5;

const MAX_FACT_N: usize = 171;
const MAX_AUX_SIZE: usize = 512;

static FACTORIALS: [f64; MAX_FACT_N] = compute_factorial_table();

const fn compute_factorial_table() -> [f64; MAX_FACT_N] {
    let mut table = [1.0; MAX_FACT_N];

    let mut val = 1.0;
    let mut n = 1;
    while n < MAX_FACT_N {
        val *= n as f64;
        table[n] = val;
        n += 1;
    }
    table
}

#[inline(always)]
fn fast_fact(n: usize) -> f64 {
    if n < MAX_FACT_N {
        unsafe { *FACTORIALS.get_unchecked(n) }
    } else {
        f64::INFINITY
    }
}

#[inline(always)]
fn fast_binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    fast_fact(n) / (fast_fact(k) * fast_fact(n - k))
}
