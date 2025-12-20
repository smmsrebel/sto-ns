//! A lightweight, `no_std`, pure Rust mathematical kernel for the **exact** evaluation of
//! two-center Coulomb integrals over $ns$ Slater-Type Orbitals (STOs).
//!
//! # STO-ns: Slater-Type Orbitals (s-orbitals)
//!
//! This library implements the analytical expansion method using
//! ellipsoidal coordinates, specifically optimized for spherically symmetric (`l=0`) orbitals.
//! It serves as a high-performance primitive for semi-empirical methods (like QEq, EEM, ReaxFF)
//! and ab initio calculations involving s-orbitals.
//!
//! ## Usage
//!
//! ### 1. Structural API (Recommended)
//!
//! Abstracting parameters into orbital objects for clear, semantic interactions.
//!
//! ```
//! use sto_ns::NsOrbital;
//!
//! // Define orbital A: n=2, zeta=1.5
//! let orb_a = NsOrbital::new(2, 1.5).expect("Invalid parameters");
//!
//! // Define orbital B: n=1, zeta=1.0
//! let orb_b = NsOrbital::new(1, 1.0).expect("Invalid parameters");
//!
//! // Calculate the Coulomb integral J(A,B) at distance R=2.0
//! let integral = orb_a.repulsion(&orb_b, 2.0);
//!
//! assert!(integral > 0.0);
//! ```
//!
//! ### 2. Direct Functional API
//!
//! Low-level access to the mathematical kernel. Ideal for tight loops or mathematical
//! backends where struct overhead is unnecessary.
//!
//! ```
//! use sto_ns::sto_coulomb_integral;
//!
//! let r = 2.0;
//! // Calculate J(2s, 1s) with exponents 1.5 and 1.0
//! let result = sto_coulomb_integral(r, 2, 1.5, 1, 1.0);
//! ```
//!
//! ## Features
//!
//! - **Exact Arithmetic**: Uses analytical formulas, not approximations.
//! - **Numerical Stability**: Automatically handles singularities at short ranges.
//! - **High Performance**:
//!   - Zero heap allocation (stack-only).
//!   - Compile-time computed factorial tables.
//!   - **`O(N)`** recursive algorithms for auxiliary functions with minimal overhead.
//! - **Portable**: `no_std` compatible (via `libm`), runs on servers, WASM, and microcontrollers.

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

/// Represents a spherically symmetric **`ns`** Slater-Type Orbital (STO).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NsOrbital {
    /// Principal quantum number `n` (must be > 0).
    pub n: u8,
    /// Orbital exponent `ζ` (must be > 0).
    pub zeta: f64,
}

impl NsOrbital {
    /// Creates a new STO. Returns `None` if parameters are unphysical (`n=0` or `ζ ≤ 0`).
    pub fn new(n: u8, zeta: f64) -> Option<Self> {
        if n > 0 && zeta > 0.0 {
            Some(Self { n, zeta })
        } else {
            None
        }
    }

    /// Computes the exact Coulomb repulsion energy with another orbital.
    ///
    /// # Arguments
    ///
    /// * `other` - The other orbital.
    /// * `r_bohr` - The internuclear distance in Bohr.
    ///
    /// # Returns
    ///
    /// * Energy in Hartree.
    #[inline]
    pub fn repulsion(&self, other: &NsOrbital, r_bohr: f64) -> f64 {
        sto_coulomb_integral(r_bohr, self.n, self.zeta, other.n, other.zeta)
    }
}

/// Computes the exact two-center Coulomb integral **`J_AB`** (**Hartree**).
///
/// This is the low-level kernel. If you have [`NsOrbital`] structs, use [`NsOrbital::repulsion`] instead.
///
/// # Arguments
///
/// * `r` - Internuclear distance in Bohr.
/// * `n_a`, `zeta_a` - Parameters for orbital on center A.
/// * `n_b`, `zeta_b` - Parameters for orbital on center B.
pub fn sto_coulomb_integral(r: f64, n_a: u8, zeta_a: f64, n_b: u8, zeta_b: f64) -> f64 {
    if r < SINGULARITY_THRESHOLD {
        return calc_one_center_limit(n_a, zeta_a, n_b, zeta_b);
    }

    let alpha_a = 2.0 * zeta_a;
    let alpha_b = 2.0 * zeta_b;

    let term_v_a = calc_potential_v(r, n_a, alpha_a);

    let m_b = (2 * n_b) as usize;
    let mut term_overlap_sum = 0.0;
    let mut alpha_b_pow = 1.0;

    for k in 0..m_b {
        let c_k = (m_b - k) as f64 / (m_b as f64 * fast_fact(k));
        let prefactor = c_k * alpha_b_pow;

        let p_a = (2 * n_a as i32) - 2;
        let p_b = k as i32 - 1;

        let overlap = calc_overlap_ellipsoidal(r, p_a, alpha_a, p_b, alpha_b);
        term_overlap_sum += prefactor * overlap;

        alpha_b_pow *= alpha_b;
    }

    let norm_a_sq = calc_norm_sq_fast(n_a, zeta_a);

    term_v_a - (term_overlap_sum * norm_a_sq)
}

#[inline]
fn calc_potential_v(r: f64, n: u8, alpha: f64) -> f64 {
    let m = (2 * n) as usize;
    let ar = alpha * r;

    let mut sum_poly = 0.0;
    for k in (0..m).rev() {
        let c_k = (m - k) as f64 / (m as f64 * fast_fact(k));
        sum_poly = sum_poly * ar + c_k;
    }

    let damping = exp(-ar) * sum_poly;
    (1.0 - damping) / r
}

fn calc_one_center_limit(n_a: u8, zeta_a: f64, n_b: u8, zeta_b: f64) -> f64 {
    let alpha_a = 2.0 * zeta_a;
    let alpha_b = 2.0 * zeta_b;
    let norm_a_sq = calc_norm_sq_fast(n_a, zeta_a);

    let p_dens_a = (2 * n_a) as usize - 2;
    let term1_val = 4.0 * PI * fast_fact(p_dens_a + 1) / powi(alpha_a, (p_dens_a + 2) as i32);

    let m_b = (2 * n_b) as usize;
    let mut term2_sum = 0.0;
    let mut alpha_b_pow = 1.0;

    for k in 0..m_b {
        let c_k = (m_b - k) as f64 / (m_b as f64 * fast_fact(k));
        let p_total = p_dens_a + k + 1;
        let integral = fast_fact(p_total) / powi(alpha_a + alpha_b, (p_total + 1) as i32);
        term2_sum += c_k * alpha_b_pow * integral;
        alpha_b_pow *= alpha_b;
    }

    norm_a_sq * (term1_val - 4.0 * PI * term2_sum)
}

#[inline]
fn calc_overlap_ellipsoidal(r: f64, pa: i32, alpha_a: f64, pb: i32, alpha_b: f64) -> f64 {
    let rho = r / 2.0;
    let p_val = rho * (alpha_a + alpha_b);
    let q_val = rho * (alpha_a - alpha_b);

    let prefactor = 2.0 * PI * powi(rho, pa + pb + 3);
    let n_exp = (pa + 1) as usize;
    let m_exp = (pb + 1) as usize;

    let max_i = n_exp + m_exp;
    if max_i >= MAX_AUX_SIZE {
        return 0.0;
    }

    let mut a_arr = [0.0; MAX_AUX_SIZE];
    let mut b_arr = [0.0; MAX_AUX_SIZE];

    fill_aux_a(&mut a_arr, max_i, p_val);
    fill_aux_b(&mut b_arr, max_i, q_val);

    let mut sum = 0.0;

    for k in 0..=n_exp {
        let bin_n = fast_binomial(n_exp, k);
        for l in 0..=m_exp {
            let bin_m = fast_binomial(m_exp, l);
            let sign = if (m_exp - l) % 2 == 1 { -1.0 } else { 1.0 };

            let i = k + l;
            let j = (n_exp - k) + (m_exp - l);

            let a_val = unsafe { *a_arr.get_unchecked(i) };
            let b_val = unsafe { *b_arr.get_unchecked(j) };

            sum += sign * bin_n * bin_m * a_val * b_val;
        }
    }
    prefactor * sum
}

fn fill_aux_a(arr: &mut [f64], max_k: usize, p: f64) {
    let exp_neg_p = exp(-p);
    let inv_p = 1.0 / p;
    let mut val = exp_neg_p * inv_p;
    arr[0] = val;

    for k in 1..=max_k {
        val = (k as f64 * val + exp_neg_p) * inv_p;
        arr[k] = val;
    }
}

fn fill_aux_b(arr: &mut [f64], max_k: usize, q: f64) {
    if abs(q) < TAYLOR_THRESHOLD {
        for k in 0..=max_k {
            let mut sum = 0.0;
            let mut q_pow = 1.0;
            for i in 0..12 {
                if (k + i) % 2 == 0 {
                    let term = q_pow / fast_fact(i) * 2.0 / ((k + i + 1) as f64);
                    sum += term;
                }
                q_pow *= -q;
            }
            arr[k] = sum;
        }
    } else {
        let exp_q = exp(q);
        let exp_neg_q = exp(-q);
        let inv_q = 1.0 / q;

        let mut val = (exp_q - exp_neg_q) * inv_q;
        arr[0] = val;

        for k in 1..=max_k {
            let sign_term = if k % 2 == 0 { exp_q } else { -exp_q };
            val = (k as f64 * val + sign_term - exp_neg_q) * inv_q;
            arr[k] = val;
        }
    }
}

#[inline]
fn calc_norm_sq_fast(n: u8, zeta: f64) -> f64 {
    let m = (2 * n) as usize;
    let two_zeta = 2.0 * zeta;
    powi(two_zeta, (m + 1) as i32) / (fast_fact(m) * 4.0 * PI)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_long_range_limit() {
        let j = sto_coulomb_integral(100.0, 1, 1.0, 1, 1.0);
        assert_relative_eq!(j, 0.01, epsilon = 1e-10);
    }

    #[test]
    fn test_1s_1s_analytical() {
        let r: f64 = 2.0;
        let zeta: f64 = 1.0;
        let rho: f64 = r * zeta;
        let term = 1.0 + 1.375 * rho + 0.75 * rho.powi(2) + (1.0 / 6.0) * rho.powi(3);
        let expected = 1.0 / r - (-2.0 * rho).exp() / r * term;

        let calc = sto_coulomb_integral(r, 1, zeta, 1, zeta);
        assert_relative_eq!(calc, expected, epsilon = 1e-12);
    }

    #[test]
    fn test_symmetry() {
        let j1 = sto_coulomb_integral(2.5, 1, 1.2, 2, 0.8);
        let j2 = sto_coulomb_integral(2.5, 2, 0.8, 1, 1.2);
        assert_relative_eq!(j1, j2, epsilon = 1e-14);
    }

    #[test]
    fn test_stability_near_zero() {
        let j_limit = sto_coulomb_integral(1e-8, 2, 1.5, 2, 1.5);
        let j_close = sto_coulomb_integral(1e-5, 2, 1.5, 2, 1.5);

        assert_relative_eq!(j_limit, j_close, epsilon = 1e-5);
    }
}
