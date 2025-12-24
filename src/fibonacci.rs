//! Fibonacci and Lucas Number Sequences
//!
//! This module implements Fibonacci and Lucas number sequences using
//! Binet's formula, which relates these sequences to the golden ratio φ.
//!
//! # Mathematical Background
//! - Fibonacci: F(n) = (φ^n - (-φ)^(-n)) / √5
//! - Lucas: L(n) = φ^n + (-φ)^(-n)
//! - Relationship: L(n) = F(n-1) + F(n+1)

use crate::golden_ratio::{PHI, PHI_INVERSE, SQRT_5};

/// Generate the nth Fibonacci number using Binet's formula
///
/// # Arguments
/// * `n` - The index in the Fibonacci sequence (0-indexed)
///
/// # Returns
/// The nth Fibonacci number
///
/// # Examples
/// ```
/// use rust_math::fibonacci::fibonacci;
///
/// assert_eq!(fibonacci(0), 0);
/// assert_eq!(fibonacci(1), 1);
/// assert_eq!(fibonacci(10), 55);
/// ```
///
/// # Mathematical Formula
/// ```text
/// F(n) = (φ^n - (-φ)^(-n)) / √5
/// ```
pub fn fibonacci(n: u64) -> u64 {
    if n <= 1 {
        return n;
    }
    let n_f = n as f64;
    let fib_n = (PHI.powf(n_f) - (-PHI).powf(-n_f)) / SQRT_5;
    fib_n.round() as u64
}

/// Generate the first n Fibonacci numbers
///
/// # Arguments
/// * `n` - The number of Fibonacci numbers to generate
///
/// # Returns
/// A vector containing the first n Fibonacci numbers
///
/// # Examples
/// ```
/// use rust_math::fibonacci::fibonacci_sequence;
///
/// let seq = fibonacci_sequence(10);
/// assert_eq!(seq, vec![0, 1, 1, 2, 3, 5, 8, 13, 21, 34]);
/// ```
pub fn fibonacci_sequence(n: usize) -> Vec<u64> {
    (0..n).map(|i| fibonacci(i as u64)).collect()
}

/// Generate the nth Lucas number
///
/// # Arguments
/// * `n` - The index in the Lucas sequence (0-indexed)
///
/// # Returns
/// The nth Lucas number
///
/// # Examples
/// ```
/// use rust_math::fibonacci::lucas;
///
/// assert_eq!(lucas(0), 2);
/// assert_eq!(lucas(1), 1);
/// assert_eq!(lucas(10), 123);
/// ```
///
/// # Mathematical Formula
/// ```text
/// L(n) = φ^n + (-φ)^(-n)
/// ```
pub fn lucas(n: u64) -> u64 {
    if n == 0 {
        return 2;
    }
    if n == 1 {
        return 1;
    }
    let n_f = n as f64;
    let lucas_n = PHI.powf(n_f) + (-PHI).powf(-n_f);
    lucas_n.round() as u64
}

/// Generate the first n Lucas numbers
///
/// # Arguments
/// * `n` - The number of Lucas numbers to generate
///
/// # Returns
/// A vector containing the first n Lucas numbers
pub fn lucas_sequence(n: usize) -> Vec<u64> {
    (0..n).map(|i| lucas(i as u64)).collect()
}

/// Check if a number is a Fibonacci number
///
/// Uses the mathematical property that n is Fibonacci if and only if
/// one or both of (5*n² + 4) or (5*n² - 4) is a perfect square.
///
/// # Arguments
/// * `n` - The number to check
///
/// # Returns
/// true if n is a Fibonacci number, false otherwise
///
/// # Examples
/// ```
/// use rust_math::fibonacci::is_fibonacci;
///
/// assert!(is_fibonacci(0));
/// assert!(is_fibonacci(1));
/// assert!(is_fibonacci(8));
/// assert!(is_fibonacci(144));
/// assert!(!is_fibonacci(4));
/// assert!(!is_fibonacci(100));
/// ```
pub fn is_fibonacci(n: u64) -> bool {
    let n_sq = n as f64;
    let test1 = 5.0 * n_sq * n_sq + 4.0;
    let test2 = 5.0 * n_sq * n_sq - 4.0;
    is_perfect_square(test1) || is_perfect_square(test2)
}

/// Check if a number is a Lucas number
pub fn is_lucas(n: u64) -> bool {
    // Generate Lucas numbers up to n and check
    let mut l = 2;
    let mut l_next = 1;
    if n == 2 {
        return true;
    }
    while l_next <= n {
        if l_next == n {
            return true;
        }
        let temp = l_next;
        l_next = l + l_next;
        l = temp;
    }
    false
}

/// Calculate the ratio F(n+1)/F(n), which converges to φ
///
/// # Arguments
/// * `n` - The index (should be >= 1)
///
/// # Returns
/// The ratio F(n+1)/F(n), or None if n is 0
pub fn fibonacci_ratio(n: u64) -> Option<f64> {
    if n == 0 {
        return None;
    }
    let f_n = fibonacci(n);
    let f_n_plus_1 = fibonacci(n + 1);
    Some(f_n_plus_1 as f64 / f_n as f64)
}

/// Calculate the ratio L(n+1)/L(n), which also converges to φ
pub fn lucas_ratio(n: u64) -> Option<f64> {
    if n == 0 {
        return None;
    }
    let l_n = lucas(n);
    let l_n_plus_1 = lucas(n + 1);
    Some(l_n_plus_1 as f64 / l_n as f64)
}

/// Find the index of a Fibonacci number
///
/// # Arguments
/// * `n` - A Fibonacci number
///
/// # Returns
/// The index of n in the Fibonacci sequence, or None if n is not Fibonacci
pub fn fibonacci_index(n: u64) -> Option<u64> {
    if !is_fibonacci(n) {
        return None;
    }
    // Handle small Fibonacci numbers explicitly due to precision issues
    if n == 0 {
        return Some(0);
    }
    if n == 1 {
        return Some(1);  // Both F(1)=1 and F(2)=1, return the first
    }
    // Use the inverse of Binet's formula
    // F(n) ≈ φ^n / √5
    // n ≈ log(F(n) * √5) / log(φ)
    let index = ((n as f64) * SQRT_5).log(PHI).round() as u64;
    if fibonacci(index) == n {
        Some(index)
    } else {
        // Try adjacent indices
        for i in index.saturating_sub(2)..=index + 2 {
            if fibonacci(i) == n {
                return Some(i);
            }
        }
        None
    }
}

fn is_perfect_square(n: f64) -> bool {
    if n < 0.0 {
        return false;
    }
    let sqrt_n = n.sqrt();
    (sqrt_n - sqrt_n.round()).abs() < 1e-10
}

/// Verify the relationship L(n) = F(n-1) + F(n+1)
///
/// # Arguments
/// * `n` - The index (should be >= 1)
///
/// # Returns
/// true if the relationship holds for the given n
pub fn verify_lucas_fibonacci_relationship(n: u64) -> bool {
    if n == 0 {
        return false;
    }
    let l_n = lucas(n);
    let f_n_minus_1 = if n > 0 { fibonacci(n - 1) } else { 0 };
    let f_n_plus_1 = fibonacci(n + 1);
    l_n == f_n_minus_1 + f_n_plus_1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fibonacci_sequence() {
        assert_eq!(fibonacci(0), 0);
        assert_eq!(fibonacci(1), 1);
        assert_eq!(fibonacci(2), 1);
        assert_eq!(fibonacci(3), 2);
        assert_eq!(fibonacci(10), 55);
        assert_eq!(fibonacci(20), 6765);
    }

    #[test]
    fn test_lucas_sequence() {
        assert_eq!(lucas(0), 2);
        assert_eq!(lucas(1), 1);
        assert_eq!(lucas(2), 3);
        assert_eq!(lucas(3), 4);
        assert_eq!(lucas(10), 123);
    }

    #[test]
    fn test_fibonacci_sequence_array() {
        let seq = fibonacci_sequence(10);
        assert_eq!(seq, vec![0, 1, 1, 2, 3, 5, 8, 13, 21, 34]);
    }

    #[test]
    fn test_is_fibonacci() {
        assert!(is_fibonacci(0));
        assert!(is_fibonacci(1));
        assert!(is_fibonacci(2));
        assert!(is_fibonacci(3));
        assert!(is_fibonacci(5));
        assert!(is_fibonacci(8));
        assert!(is_fibonacci(13));
        assert!(is_fibonacci(144));
        assert!(!is_fibonacci(4));
        assert!(!is_fibonacci(6));
        assert!(!is_fibonacci(7));
        assert!(!is_fibonacci(100));
    }

    #[test]
    fn test_fibonacci_ratio_convergence() {
        // As n increases, F(n+1)/F(n) should approach φ
        let ratio_10 = fibonacci_ratio(10).unwrap();
        let ratio_20 = fibonacci_ratio(20).unwrap();

        let error_10 = (ratio_10 - PHI).abs();
        let error_20 = (ratio_20 - PHI).abs();

        // Higher n should give better approximation
        assert!(error_20 < error_10);
        assert!((ratio_20 - PHI).abs() < 0.0001);
    }

    #[test]
    fn test_lucas_fibonacci_relationship() {
        assert!(verify_lucas_fibonacci_relationship(1));
        assert!(verify_lucas_fibonacci_relationship(5));
        assert!(verify_lucas_fibonacci_relationship(10));
    }

    #[test]
    fn test_fibonacci_index() {
        assert_eq!(fibonacci_index(0), Some(0));
        assert_eq!(fibonacci_index(1), Some(1));
        assert_eq!(fibonacci_index(2), Some(3));
        assert_eq!(fibonacci_index(8), Some(6));
        assert_eq!(fibonacci_index(55), Some(10));
        assert_eq!(fibonacci_index(4), None);
    }
}
