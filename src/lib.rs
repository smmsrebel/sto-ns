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
