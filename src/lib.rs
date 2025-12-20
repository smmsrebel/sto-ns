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

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NsOrbital {
    pub n: u8,
    pub zeta: f64,
}

impl NsOrbital {
    pub fn new(n: u8, zeta: f64) -> Option<Self> {
        if n > 0 && zeta > 0.0 {
            Some(Self { n, zeta })
        } else {
            None
        }
    }

    #[inline]
    pub fn repulsion(&self, other: &NsOrbital, r_bohr: f64) -> f64 {
        sto_coulomb_integral(r_bohr, self.n, self.zeta, other.n, other.zeta)
    }
}

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
