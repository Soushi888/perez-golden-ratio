//! Pascal's Triangle and the Perez Hourglass Structure
//!
//! This module implements Pascal's triangle and the simplified Perez Hourglass
//! visualization. Note that the actual Perez Hourglass involves recursive
//! parity filtering and mod 2 operations beyond this basic visualization.
//!
//! # Mathematical Background
//! - Pascal's triangle: Binomial coefficients arranged in triangular form
//! - Each entry is the sum of the two entries above it
//! - Row n (starting from 0) contains C(n, k) for k = 0 to n
//!
//! # Perez Hourglass (Simplified)
//! The full Perez Hourglass involves:
//! - Recursive parity filtering (mod 2 operations)
//! - This creates Sierpiński-like fractal patterns
//! - Connects to OEIS A000975 (Lichtenberg sequence)
//!
//! This module provides the basic structure; the full fractal properties
//! require parity-filtered versions.

/// Generate Pascal's triangle
///
/// # Arguments
/// * `n` - The number of rows to generate (0-indexed, so n+1 total rows)
///
/// # Returns
/// A 2D vector where each inner vector is a row of Pascal's triangle
///
/// # Examples
/// ```
/// use rust_math::hourglass::pascal_triangle;
///
/// let triangle = pascal_triangle(4);
/// // Row 0: [1]
/// // Row 1: [1, 1]
/// // Row 2: [1, 2, 1]
/// // Row 3: [1, 3, 3, 1]
/// // Row 4: [1, 4, 6, 4, 1]
/// ```
pub fn pascal_triangle(n: usize) -> Vec<Vec<u64>> {
    let mut triangle = vec![vec![1u64]];

    for i in 1..=n {
        let prev = &triangle[i - 1];
        let mut row = vec![1u64];

        for j in 0..prev.len() - 1 {
            let sum = prev[j] + prev[j + 1];
            row.push(sum);
        }
        row.push(1);
        triangle.push(row);
    }
    triangle
}

/// Generate a single row of Pascal's triangle
///
/// # Arguments
/// * `n` - The row index (0-indexed)
///
/// # Returns
/// The nth row of Pascal's triangle as a vector
///
/// # Examples
/// ```
/// use rust_math::hourglass::pascal_row;
///
/// assert_eq!(pascal_row(0), vec![1]);
/// assert_eq!(pascal_row(1), vec![1, 1]);
/// assert_eq!(pascal_row(4), vec![1, 4, 6, 4, 1]);
/// ```
pub fn pascal_row(n: usize) -> Vec<u64> {
    if n == 0 {
        return vec![1];
    }

    let mut row = vec![0u64; n + 1];
    row[0] = 1;

    for i in 1..=n {
        for j in (1..=i).rev() {
            row[j] += row[j - 1];
        }
    }
    row
}

/// Calculate a specific binomial coefficient C(n, k)
///
/// Uses the multiplicative formula for efficiency
///
/// # Arguments
/// * `n` - The total number of items
/// * `k` - The number of items to choose
///
/// # Returns
/// The binomial coefficient C(n, k)
///
/// # Examples
/// ```
/// use rust_math::hourglass::binomial;
///
/// assert_eq!(binomial(5, 2), 10);
/// assert_eq!(binomial(10, 3), 120);
/// ```
pub fn binomial(n: u64, k: u64) -> u64 {
    if k > n {
        return 0;
    }
    if k == 0 || k == n {
        return 1;
    }
    if k > n - k {
        return binomial(n, n - k);
    }

    let mut result = 1u64;
    for i in 0..k {
        result = result.saturating_mul(n - i);
        result /= i + 1;
    }
    result
}

/// Generate Pascal's triangle modulo 2 (Sierpiński triangle)
///
/// This creates the famous Sierpiński triangle fractal pattern.
/// The actual Perez Hourglass uses parity filtering related to this.
///
/// # Arguments
/// * `n` - Number of rows to generate
///
/// # Returns
/// A 2D vector of 0s and 1s representing mod-2 Pascal triangle
pub fn pascal_mod_2(n: usize) -> Vec<Vec<u8>> {
    let mut triangle = vec![vec![1u8]];

    for i in 1..=n {
        let prev = &triangle[i - 1];
        let mut row = vec![1u8];

        for j in 0..prev.len() - 1 {
            // Mod 2 addition (XOR for binary)
            row.push((prev[j] + prev[j + 1]) % 2);
        }
        row.push(1);
        triangle.push(row);
    }
    triangle
}

/// Check if a position in Pascal's triangle contains an odd number
///
/// A position C(n, k) is odd if and only if (k & (n - k)) == 0
/// (using bitwise AND). This is a known result from Lucas's theorem.
///
/// # Arguments
/// * `n` - Row index
/// * `k` - Column index
///
/// # Returns
/// true if C(n, k) is odd, false otherwise
pub fn is_pascal_odd(n: u64, k: u64) -> bool {
    if k > n {
        return false;
    }
    // Using Lucas's theorem: C(n,k) is odd iff (k & (n-k)) == 0
    (k & (n - k)) == 0
}

/// Generate the Lichtenberg sequence (OEIS A000975)
///
/// The Lichtenberg sequence: 1, 2, 5, 10, 21, 42, 85, ...
/// Formula: a(n) = 2*a(n-1)     for odd n (indices 1, 3, 5, ...)
///          a(n) = 2*a(n-1) + 1 for even n (indices 2, 4, 6, ...)
///
/// This sequence appears in the Perez Hourglass structure.
///
/// # Arguments
/// * `n` - The number of terms to generate
///
/// # Returns
/// A vector containing the first n terms of the sequence
pub fn lichtenberg_sequence(n: usize) -> Vec<u64> {
    if n == 0 {
        return vec![];
    }

    let mut seq = vec![0u64; n];
    seq[0] = 1;

    for i in 1..n {
        if i % 2 == 1 {
            // Odd index: just multiply by 2
            seq[i] = 2 * seq[i - 1];
        } else {
            // Even index: multiply by 2 and add 1
            seq[i] = 2 * seq[i - 1] + 1;
        }
    }
    seq
}

/// Simplified Perez Hourglass (basic visualization)
///
/// Note: The actual Perez Hourglass involves recursive parity filtering
/// and infinite-dimensional folding beyond this basic structure.
///
/// This creates a mirrored structure where:
/// - Upper half: Pascal's triangle with positive values
/// - Lower half: Mirror image with negated values
///
/// # Arguments
/// * `n` - Size parameter for the hourglass
///
/// # Returns
/// A 2D vector representing the hourglass structure
///
/// # Examples
/// ```
/// use rust_math::hourglass::perez_hourglass_simple;
///
/// let hourglass = perez_hourglass_simple(3);
/// // Produces:
/// // [1]
/// // [1, 1]
/// // [1, 2, 1]
/// // [-2, -1]
/// // [-1]
pub fn perez_hourglass_simple(n: usize) -> Vec<Vec<i64>> {
    // For the hourglass, we want rows 0 to n-1 for the upper part,
    // then rows n-1 down to 0 negated for the lower part (excluding row 0).
    // This creates a mirrored structure without the middle row duplication.
    let upper = pascal_triangle(n.saturating_sub(1));
    let mut hourglass: Vec<Vec<i64>> = upper
        .iter()
        .map(|row| row.iter().map(|&x| x as i64).collect())
        .collect();

    // Partie négative (miroir) - mirror with negation
    let lower: Vec<Vec<i64>> = upper
        .iter()
        .rev()
        .skip(1)  // Skip the first row of reversed to avoid duplication
        .map(|row| row.iter().map(|&x| -(x as i64)).collect())
        .collect();

    hourglass.extend(lower);
    hourglass
}

/// Parity-filtered Perez Hourglass (closer to actual concept)
///
/// This applies mod 2 filtering to Pascal's triangle before
/// creating the hourglass structure, which is closer to Perez's
/// actual fractal construction.
///
/// # Arguments
/// * `n` - Size parameter for the hourglass
///
/// # Returns
/// A 2D vector representing the parity-filtered hourglass
pub fn perez_hourglass_parity(n: usize) -> Vec<Vec<i8>> {
    let mod2_triangle = pascal_mod_2(n);
    let mut hourglass: Vec<Vec<i8>> = mod2_triangle
        .iter()
        .map(|row| row.iter().map(|&x| x as i8).collect())
        .collect();

    // Mirror with negation
    let lower: Vec<Vec<i8>> = mod2_triangle
        .iter()
        .rev()
        .skip(1)
        .map(|row| row.iter().map(|&x| -(x as i8)).collect())
        .collect();

    hourglass.extend(lower);
    hourglass
}

/// Display Pascal's triangle in a formatted way
///
/// # Arguments
/// * `triangle` - The triangle to display
///
/// # Returns
/// A string representation of the triangle
pub fn format_pascal_triangle(triangle: &[Vec<u64>]) -> String {
    let mut result = String::new();

    for (i, row) in triangle.iter().enumerate() {
        // Add leading spaces for centering
        let indent = (triangle.len() - i) / 2;
        for _ in 0..indent {
            result.push(' ');
        }

        // Add row values
        for (j, val) in row.iter().enumerate() {
            if j > 0 {
                result.push(' ');
            }
            result.push_str(&val.to_string());
        }
        result.push('\n');
    }
    result
}

/// Calculate the sum of a row in Pascal's triangle
///
/// The sum of row n is 2^n
///
/// # Arguments
/// * `n` - The row index
///
/// # Returns
/// The sum of all elements in row n
pub fn pascal_row_sum(n: u64) -> u64 {
    1u64 << n  // 2^n using bit shift
}

/// Find the maximum value in a given row of Pascal's triangle
///
/// # Arguments
/// * `n` - The row index
///
/// # Returns
/// The maximum binomial coefficient C(n, k) for 0 ≤ k ≤ n
///
/// # Notes
/// For even n: max is at C(n, n/2)
/// For odd n: max is at C(n, (n-1)/2) or C(n, (n+1)/2)
pub fn pascal_row_max(n: u64) -> u64 {
    binomial(n, n / 2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pascal_triangle() {
        let triangle = pascal_triangle(4);
        assert_eq!(triangle[0], vec![1]);
        assert_eq!(triangle[1], vec![1, 1]);
        assert_eq!(triangle[2], vec![1, 2, 1]);
        assert_eq!(triangle[3], vec![1, 3, 3, 1]);
        assert_eq!(triangle[4], vec![1, 4, 6, 4, 1]);
    }

    #[test]
    fn test_pascal_row() {
        assert_eq!(pascal_row(0), vec![1]);
        assert_eq!(pascal_row(5), vec![1, 5, 10, 10, 5, 1]);
    }

    #[test]
    fn test_binomial() {
        assert_eq!(binomial(5, 2), 10);
        assert_eq!(binomial(10, 3), 120);
        assert_eq!(binomial(0, 0), 1);
        assert_eq!(binomial(5, 0), 1);
        assert_eq!(binomial(5, 5), 1);
    }

    #[test]
    fn test_pascal_mod_2() {
        let triangle = pascal_mod_2(4);
        assert_eq!(triangle[0], vec![1]);
        assert_eq!(triangle[1], vec![1, 1]);
        assert_eq!(triangle[2], vec![1, 0, 1]);
        assert_eq!(triangle[3], vec![1, 1, 1, 1]);
        assert_eq!(triangle[4], vec![1, 0, 0, 0, 1]);
    }

    #[test]
    fn test_is_pascal_odd() {
        assert!(is_pascal_odd(5, 0));  // 1 is odd
        assert!(is_pascal_odd(5, 1));  // 5 is odd
        assert!(!is_pascal_odd(5, 2)); // 10 is even
        assert!(!is_pascal_odd(5, 3)); // 10 is even
        assert!(!is_pascal_odd(4, 2)); // 6 is even
    }

    #[test]
    fn test_lichtenberg_sequence() {
        let seq = lichtenberg_sequence(8);
        assert_eq!(seq, vec![1, 2, 5, 10, 21, 42, 85, 170]);
    }

    #[test]
    fn test_pascal_row_sum() {
        assert_eq!(pascal_row_sum(0), 1);
        assert_eq!(pascal_row_sum(1), 2);
        assert_eq!(pascal_row_sum(5), 32);
    }

    #[test]
    fn test_pascal_row_max() {
        assert_eq!(pascal_row_max(4), 6);   // C(4,2) = 6
        assert_eq!(pascal_row_max(5), 10);  // C(5,2) = 10
    }

    #[test]
    fn test_perez_hourglass_simple() {
        let hourglass = perez_hourglass_simple(2);
        // Should have upper triangle plus mirrored lower
        assert_eq!(hourglass[0], vec![1]);
        assert_eq!(hourglass[1], vec![1, 1]);
        assert_eq!(hourglass[2], vec![-1]);
    }

    #[test]
    fn test_format_pascal_triangle() {
        let triangle = pascal_triangle(3);
        let formatted = format_pascal_triangle(&triangle);
        assert!(formatted.contains('1'));
        assert!(formatted.contains('\n'));
    }
}
