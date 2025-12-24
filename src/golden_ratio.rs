//! Golden Ratio (φ) Constants and Mathematical Functions
//!
//! This module provides the fundamental golden ratio constants and
//! related mathematical functions that form the foundation of
//! Jean-Claude Perez's discoveries.

/// The Golden Ratio: φ = (1 + √5) / 2 ≈ 1.618033988749...
///
/// # Examples
/// ```
/// use rust_math::golden_ratio::PHI;
/// assert!((PHI - 1.618).abs() < 0.001);
/// ```
pub const PHI: f64 = 1.618_033_988_749_894_848_2;

/// φ² = φ + 1 ≈ 2.618033988749...
///
/// This is a key property: φ² = φ + 1
pub const PHI_SQUARED: f64 = 2.618_033_988_749_894_848_2;

/// 1/φ ≈ 0.618033988749...
///
/// Also equals φ - 1
pub const PHI_INVERSE: f64 = 1.0 / PHI;

/// √5 ≈ 2.236067977499...
///
/// Used in Binet's formula for Fibonacci numbers
pub const SQRT_5: f64 = 2.236_067_977_499_79;

/// Perez's Ratio: (3 - φ) / 2 ≈ 0.690983005625...
///
/// This is Jean-Claude Perez's key discovery in DNA:
/// the ratio (T + A) / (C + G) converges to this value.
///
/// # Mathematical Derivation
/// - From the golden ratio φ ≈ 1.618
/// - (3 - φ) / 2 = (3 - 1.618...) / 2 ≈ 0.691
/// - This ratio appears in codon populations across species
pub const PEREZ_RATIO: f64 = (3.0 - PHI) / 2.0;

/// Calculate the nth power of φ
///
/// # Arguments
/// * `n` - The exponent
///
/// # Examples
/// ```
/// use rust_math::golden_ratio::phi_pow;
/// assert!((phi_pow(2) - 2.618).abs() < 0.001);
/// ```
pub fn phi_pow(n: u32) -> f64 {
    PHI.powi(n as i32)
}

/// Calculate φ^n where n can be negative
///
/// For negative n: φ^(-n) = (-1)^n / φ^n
pub fn phi_pow_signed(n: i32) -> f64 {
    if n >= 0 {
        PHI.powi(n)
    } else {
        let abs_n = n.abs();
        let sign = if abs_n % 2 == 0 { 1.0 } else { -1.0 };
        sign / PHI.powi(abs_n)
    }
}

/// Check if a value is close to φ within tolerance
///
/// # Arguments
/// * `value` - The value to check
/// * `tolerance` - Maximum allowed difference (default: 1e-10)
pub fn is_golden_ratio(value: f64, tolerance: Option<f64>) -> bool {
    let tol = tolerance.unwrap_or(1e-10);
    (value - PHI).abs() < tol
}

/// Check if two values are in golden ratio proportion
///
/// Returns true if a/b ≈ φ or b/a ≈ φ
pub fn are_golden_proportion(a: f64, b: f64, tolerance: Option<f64>) -> bool {
    let tol = tolerance.unwrap_or(1e-6);
    if b == 0.0 {
        return false;
    }
    let ratio = a / b;
    (ratio - PHI).abs() < tol || (ratio - PHI_INVERSE).abs() < tol
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_phi_properties() {
        // φ² = φ + 1
        assert!((PHI_SQUARED - (PHI + 1.0)).abs() < 1e-10);

        // 1/φ = φ - 1
        assert!((PHI_INVERSE - (PHI - 1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_perez_ratio() {
        // (3-φ)/2 ≈ 0.691
        assert!((PEREZ_RATIO - 0.690983005625).abs() < 1e-9);
    }

    #[test]
    fn test_phi_pow() {
        assert!((phi_pow(0) - 1.0).abs() < 1e-10);
        assert!((phi_pow(1) - PHI).abs() < 1e-10);
        assert!((phi_pow(2) - PHI_SQUARED).abs() < 1e-10);
    }

    #[test]
    fn test_golden_proportion() {
        assert!(are_golden_proportion(PHI, 1.0, None));
        assert!(are_golden_proportion(1.0, PHI_INVERSE, None));
    }
}
