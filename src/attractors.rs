//! Fractal Attractor Analysis for Codon Populations
//!
//! This module implements the fractal attractor analysis discovered by
//! Jean-Claude Perez in 2010, where codon frequencies cluster around
//! two attractors related to the golden ratio φ.
//!
//! # Key Findings
//! - Codon populations cluster around **2 attractors**
//! - Attractor 1 ≈ 1/φ ≈ 0.618 (high frequency)
//! - Attractor 2 ≈ 1 - 1/φ ≈ 0.382 (low frequency)
//! - Ratio between attractors ≈ 2.0
//!
//! # References
//! Perez, J.C. (2010). "Codon Populations in Single-Stranded Whole Human
//! Genome DNA Are Fractal and Fine-Tuned by the Golden Ratio 1.618"

use crate::golden_ratio::{PHI, PHI_INVERSE};

/// Fractal attractor for codon populations
///
/// Based on Perez's 2010 discovery that codon frequencies
/// cluster around two attractors related to the golden ratio.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FractalAttractor {
    /// High-frequency attractor (usually ≈ 0.618 = 1/φ)
    pub attractor_1: f64,
    /// Low-frequency attractor (usually ≈ 0.382 = 1 - 1/φ)
    pub attractor_2: f64,
    /// Ratio between attractors (should be ≈ 2.0)
    pub ratio: f64,
}

impl FractalAttractor {
    /// Create a new fractal attractor struct
    pub fn new(attractor_1: f64, attractor_2: f64) -> Self {
        let ratio = if attractor_2 != 0.0 {
            attractor_1 / attractor_2
        } else {
            0.0
        };

        FractalAttractor {
            attractor_1,
            attractor_2,
            ratio,
        }
    }

    /// Analyze codon frequencies to detect fractal attractors
    ///
    /// This implements a simplified version of Perez's clustering method.
    /// The algorithm:
    /// 1. Find the median of all frequencies
    /// 2. Split frequencies into high and low groups
    /// 3. Calculate the mean of each group
    /// 4. These means are the attractors
    ///
    /// # Arguments
    /// * `frequencies` - A slice of codon frequency values
    ///
    /// # Returns
    /// A FractalAttractor struct with the detected attractors
    ///
    /// # Examples
    /// ```
    /// use rust_math::attractors::FractalAttractor;
    ///
    /// let freqs = vec![0.6, 0.62, 0.4, 0.38, 0.59, 0.41];
    /// let attractor = FractalAttractor::from_codon_frequencies(&freqs);
    /// ```
    pub fn from_codon_frequencies(frequencies: &[f64]) -> Self {
        if frequencies.is_empty() {
            return FractalAttractor::new(0.0, 0.0);
        }

        let mut sum_low = 0.0;
        let mut sum_high = 0.0;
        let mut count_low = 0usize;
        let mut count_high = 0usize;

        // Find median
        let median = {
            let mut sorted = frequencies.to_vec();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            sorted[sorted.len() / 2]
        };

        // Split and average
        for &freq in frequencies {
            if freq <= median {
                sum_low += freq;
                count_low += 1;
            } else {
                sum_high += freq;
                count_high += 1;
            }
        }

        let attractor_1 = if count_high > 0 {
            sum_high / count_high as f64
        } else {
            0.0
        };

        let attractor_2 = if count_low > 0 {
            sum_low / count_low as f64
        } else {
            0.0
        };

        FractalAttractor::new(attractor_1, attractor_2)
    }

    /// Alternative clustering using k-means approach
    ///
    /// This uses a simplified k-means with k=2, initialized
    /// at expected golden ratio positions.
    pub fn from_codon_frequencies_kmeans(frequencies: &[f64], iterations: usize) -> Self {
        if frequencies.is_empty() {
            return FractalAttractor::new(0.0, 0.0);
        }

        // Initialize at expected golden ratio positions
        let mut center_1 = PHI_INVERSE;  // ≈ 0.618
        let mut center_2 = 1.0 - PHI_INVERSE;  // ≈ 0.382

        for _ in 0..iterations {
            let mut cluster_1: Vec<f64> = Vec::new();
            let mut cluster_2: Vec<f64> = Vec::new();

            // Assign to nearest cluster
            for &freq in frequencies {
                let dist_1 = (freq - center_1).abs();
                let dist_2 = (freq - center_2).abs();

                if dist_1 < dist_2 {
                    cluster_1.push(freq);
                } else {
                    cluster_2.push(freq);
                }
            }

            // Update centers
            if !cluster_1.is_empty() {
                center_1 = cluster_1.iter().sum::<f64>() / cluster_1.len() as f64;
            }
            if !cluster_2.is_empty() {
                center_2 = cluster_2.iter().sum::<f64>() / cluster_2.len() as f64;
            }
        }

        FractalAttractor::new(center_1, center_2)
    }

    /// Check if attractors follow golden ratio patterns
    ///
    /// Returns true if:
    /// - Ratio is close to φ ≈ 1.618 (attractor1/attractor2)
    /// - Attractor 1 is close to 1/φ ≈ 0.618
    /// - Attractor 2 is close to 1 - 1/φ ≈ 0.382
    ///
    /// Note: The ratio of attractors equals φ because:
    /// (1/φ) / (1 - 1/φ) = (1/φ) / ((φ-1)/φ) = 1/(φ-1) = φ
    ///
    /// # Arguments
    /// * `tolerance` - Maximum allowed deviation (default: 0.1)
    pub fn validates_golden_ratio(&self, tolerance: Option<f64>) -> bool {
        let tol = tolerance.unwrap_or(0.1);

        // Ratio should equal PHI, not 2.0
        let ratio_close = (self.ratio - PHI).abs() < tol;
        let attractor1_close = (self.attractor_1 - PHI_INVERSE).abs() < tol;
        let attractor2_close = (self.attractor_2 - (1.0 - PHI_INVERSE)).abs() < tol;

        ratio_close && attractor1_close && attractor2_close
    }

    /// Calculate the "fractal dimension" based on attractor spacing
    ///
    /// This is a simplified metric based on Perez's work
    pub fn fractal_dimension(&self) -> f64 {
        if self.attractor_1 == 0.0 || self.attractor_2 == 0.0 {
            return 0.0;
        }
        (self.attractor_1 / self.attractor_2).log(PHI)
    }

    /// Calculate deviation from expected golden ratio values
    ///
    /// Returns a tuple of (attractor1_error, attractor2_error, ratio_error)
    pub fn deviation_from_golden_ratio(&self) -> (f64, f64, f64) {
        let expected_1 = PHI_INVERSE;
        let expected_2 = 1.0 - PHI_INVERSE;
        let expected_ratio = PHI;  // Ratio equals φ, not 2.0

        (
            (self.attractor_1 - expected_1).abs() / expected_1,
            (self.attractor_2 - expected_2).abs() / expected_2,
            (self.ratio - expected_ratio).abs() / expected_ratio,
        )
    }

    /// Get a summary of the attractor analysis
    pub fn summary(&self) -> String {
        format!(
            "FractalAttractor:\n  Attractor 1 (high): {:.6} (expected: {:.6})\n  Attractor 2 (low):  {:.6} (expected: {:.6})\n  Ratio: {:.6} (expected: {:.6} = φ)\n  Golden ratio validated: {}",
            self.attractor_1,
            PHI_INVERSE,
            self.attractor_2,
            1.0 - PHI_INVERSE,
            self.ratio,
            PHI,
            if self.validates_golden_ratio(None) { "YES" } else { "NO" }
        )
    }
}

/// Analyze a sequence and detect periodic patterns
///
/// This can reveal self-similar structures in DNA
pub fn analyze_periodic_pattern(sequence: &[f64], window_size: usize) -> Vec<f64> {
    if sequence.len() < window_size {
        return vec![];
    }

    let mut patterns = Vec::new();

    for i in 0..=(sequence.len() - window_size) {
        let window: Vec<f64> = sequence[i..i + window_size].to_vec();
        let mean: f64 = window.iter().sum::<f64>() / window.len() as f64;
        patterns.push(mean);
    }

    patterns
}

/// Calculate the autocorrelation of a sequence
///
/// Autocorrelation can reveal repeating patterns and periodicities
pub fn autocorrelation(sequence: &[f64], max_lag: usize) -> Vec<f64> {
    if sequence.is_empty() {
        return vec![];
    }

    let n = sequence.len();
    let mean: f64 = sequence.iter().sum::<f64>() / n as f64;
    let variance: f64 = sequence.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n as f64;

    if variance == 0.0 {
        return vec![1.0; max_lag + 1];
    }

    let mut result = Vec::with_capacity(max_lag + 1);

    for lag in 0..=max_lag {
        if lag >= n {
            result.push(0.0);
            continue;
        }

        let mut covariance = 0.0;
        for i in 0..(n - lag) {
            covariance += (sequence[i] - mean) * (sequence[i + lag] - mean);
        }
        covariance /= n as f64;

        let correlation = covariance / variance;
        result.push(correlation);
    }

    result
}

/// Test if a sequence exhibits fractal properties
///
/// Uses Hurst exponent estimation (simplified)
pub fn estimate_hurst_exponent(sequence: &[f64]) -> Option<f64> {
    if sequence.len() < 4 {
        return None;
    }

    // Simplified rescaled range analysis
    let n = sequence.len();
    let mean: f64 = sequence.iter().sum::<f64>() / n as f64;

    // Calculate cumulative deviations
    let mut cumulative_dev = Vec::with_capacity(n);
    let mut sum = 0.0;
    for &value in sequence {
        sum += value - mean;
        cumulative_dev.push(sum);
    }

    // Range
    let max_dev = *cumulative_dev.iter().max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)).unwrap_or(&0.0);
    let min_dev = *cumulative_dev.iter().min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)).unwrap_or(&0.0);
    let range = max_dev - min_dev;

    // Standard deviation
    let std_dev: f64 = sequence.iter().map(|&x| (x - mean).powi(2)).sum::<f64>().sqrt() / (n as f64).sqrt();

    if std_dev == 0.0 {
        return None;
    }

    // R/S ratio (simplified)
    let rs = range / std_dev;

    // Very rough Hurst estimate (for demonstration)
    // H ≈ log(R/S) / log(n)
    Some(rs.log((n as f64).log2()))
}

/// Simulate codon frequencies with golden ratio clustering
///
/// This generates synthetic codon frequencies that follow
/// Perez's discovered attractor patterns
pub fn simulate_codon_frequencies(
    n_codons: usize,
    high_attractor_ratio: f64,
    noise_level: f64,
) -> Vec<f64> {
    let n_high = (n_codons as f64 * high_attractor_ratio) as usize;
    let n_low = n_codons - n_high;

    let mut frequencies = Vec::with_capacity(n_codons);

    // High attractor codons
    for _ in 0..n_high {
        let base = PHI_INVERSE;
        let noise = (rand() - 0.5) * 2.0 * noise_level;
        frequencies.push((base + noise).clamp(0.0, 1.0));
    }

    // Low attractor codons
    for _ in 0..n_low {
        let base = 1.0 - PHI_INVERSE;
        let noise = (rand() - 0.5) * 2.0 * noise_level;
        frequencies.push((base + noise).clamp(0.0, 1.0));
    }

    frequencies
}

/// Simple random number generator [0, 1)
fn rand() -> f64 {
    use std::sync::atomic::{AtomicU64, Ordering};
    static SEED: AtomicU64 = AtomicU64::new(42);
    let s = SEED.fetch_add(1, Ordering::Relaxed);
    ((s as f64) % 1000.0) / 1000.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fractal_attractor_creation() {
        let attractor = FractalAttractor::new(0.618, 0.382);
        assert!((attractor.attractor_1 - 0.618).abs() < 0.001);
        assert!((attractor.attractor_2 - 0.382).abs() < 0.001);
        assert!((attractor.ratio - 1.618).abs() < 0.01);
    }

    #[test]
    fn test_from_codon_frequencies() {
        let freqs = vec![0.6, 0.62, 0.4, 0.38, 0.59, 0.41];
        let attractor = FractalAttractor::from_codon_frequencies(&freqs);

        // Should detect two clusters
        assert!(attractor.attractor_1 > 0.5);
        assert!(attractor.attractor_2 < 0.5);
    }

    #[test]
    fn test_validates_golden_ratio() {
        // Perfect golden ratio attractors
        let attractor = FractalAttractor::new(PHI_INVERSE, 1.0 - PHI_INVERSE);
        assert!(attractor.validates_golden_ratio(None));

        // Far from golden ratio
        let attractor_bad = FractalAttractor::new(0.9, 0.1);
        assert!(!attractor_bad.validates_golden_ratio(None));
    }

    #[test]
    fn test_simulate_codon_frequencies() {
        let freqs = simulate_codon_frequencies(64, 0.375, 0.05);

        // Should have the right number
        assert_eq!(freqs.len(), 64);

        // All values should be in valid range
        for &f in &freqs {
            assert!(f >= 0.0 && f <= 1.0);
        }
    }

    #[test]
    fn test_autocorrelation() {
        let sequence = vec![1.0, 2.0, 1.0, 2.0, 1.0, 2.0];
        let ac = autocorrelation(&sequence, 3);

        // Lag 0 should always be 1.0
        assert!((ac[0] - 1.0).abs() < 0.001);

        // Periodic pattern should show in autocorrelation
        assert!(ac.len() == 4);
    }

    #[test]
    fn test_fractal_dimension() {
        let attractor = FractalAttractor::new(PHI_INVERSE, 1.0 - PHI_INVERSE);
        let dim = attractor.fractal_dimension();

        // Should be positive
        assert!(dim > 0.0);
    }

    #[test]
    fn test_deviation_from_golden_ratio() {
        // Perfect match
        let attractor = FractalAttractor::new(PHI_INVERSE, 1.0 - PHI_INVERSE);
        let (dev1, dev2, dev_ratio) = attractor.deviation_from_golden_ratio();

        assert!(dev1 < 0.01);
        assert!(dev2 < 0.01);
        // Note: The ratio of attractors is PHI, not 2.0
        // attractor_1 / attractor_2 = (1/φ) / (1 - 1/φ) = φ
        let expected_ratio = PHI;
        let actual_ratio = attractor.ratio;
        assert!((actual_ratio - expected_ratio).abs() / expected_ratio < 0.01);
    }
}
