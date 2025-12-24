//! DNA Sequence Analysis Based on Jean-Claude Perez's Discoveries
//!
//! This module provides tools for analyzing DNA sequences using the
//! golden ratio patterns discovered by Jean-Claude Perez in the 1990s.
//!
//! # Key Discoveries
//! - The ratio (T+A)/(C+G) converges to Perez ratio ≈ 0.691
//! - Codon populations cluster around golden ratio attractors
//! - These patterns appear across species and genomic scales
//!
//! # References
//! - Perez, J.C. (1991). "Chaos DNA and neuro-computers: a golden link"
//! - Perez, J.C. (2010). "Codon Populations in Single-Stranded Whole Human
//!   Genome DNA Are Fractal and Fine-Tuned by the Golden Ratio 1.618"

use crate::golden_ratio::PEREZ_RATIO;

/// DNA nucleotides (bases)
///
/// The four nucleotides that form the building blocks of DNA:
/// - **A**denine: Purine, pairs with Thymine
/// - **T**hymine: Pyrimidine, pairs with Adenine
/// - **C**ytosine: Pyrimidine, pairs with Guanine
/// - **G**uanine: Purine, pairs with Cytosine
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Nucleotide {
    /// Adenine (A)
    Adenine,
    /// Thymine (T)
    Thymine,
    /// Cytosine (C)
    Cytosine,
    /// Guanine (G)
    Guanine,
}

impl Nucleotide {
    /// Create a nucleotide from a single character
    ///
    /// # Examples
    /// ```
    /// use rust_math::dna::Nucleotide;
    ///
    /// assert_eq!(Nucleotide::from_char('A'), Some(Nucleotide::Adenine));
    /// assert_eq!(Nucleotide::from_char('t'), Some(Nucleotide::Thymine));
    /// assert_eq!(Nucleotide::from_char('X'), None);
    /// ```
    pub fn from_char(c: char) -> Option<Self> {
        match c.to_ascii_uppercase() {
            'A' => Some(Nucleotide::Adenine),
            'T' => Some(Nucleotide::Thymine),
            'C' => Some(Nucleotide::Cytosine),
            'G' => Some(Nucleotide::Guanine),
            _ => None,
        }
    }

    /// Convert nucleotide to single character
    pub fn to_char(self) -> char {
        match self {
            Nucleotide::Adenine => 'A',
            Nucleotide::Thymine => 'T',
            Nucleotide::Cytosine => 'C',
            Nucleotide::Guanine => 'G',
        }
    }

    /// Check if nucleotide is a purine (A or G)
    ///
    /// Purines have a double-ring structure
    pub fn is_purine(&self) -> bool {
        matches!(self, Nucleotide::Adenine | Nucleotide::Guanine)
    }

    /// Check if nucleotide is a pyrimidine (C or T)
    ///
    /// Pyrimidines have a single-ring structure
    pub fn is_pyrimidine(&self) -> bool {
        !self.is_purine()
    }

    /// Check if nucleotide is amino (A or C)
    ///
    /// Amino groups in DNA structure
    pub fn is_amino(&self) -> bool {
        matches!(self, Nucleotide::Adenine | Nucleotide::Cytosine)
    }

    /// Check if nucleotide is keto (G or T)
    ///
    /// Keto groups in DNA structure
    pub fn is_keto(&self) -> bool {
        !self.is_amino()
    }

    /// Check if nucleotide forms strong hydrogen bonds (G or C)
    ///
    /// G-C bonds have 3 hydrogen bonds (stronger)
    pub fn is_strong(&self) -> bool {
        matches!(self, Nucleotide::Guanine | Nucleotide::Cytosine)
    }

    /// Check if nucleotide forms weak hydrogen bonds (A or T)
    ///
    /// A-T bonds have 2 hydrogen bonds (weaker)
    pub fn is_weak(&self) -> bool {
        !self.is_strong()
    }

    /// Get the complementary nucleotide (Watson-Crick base pairing)
    ///
    /// # Examples
    /// ```
    /// use rust_math::dna::Nucleotide;
    ///
    /// assert_eq!(Nucleotide::Adenine.complement(), Nucleotide::Thymine);
    /// assert_eq!(Nucleotide::Guanine.complement(), Nucleotide::Cytosine);
    /// ```
    pub fn complement(&self) -> Self {
        match self {
            Nucleotide::Adenine => Nucleotide::Thymine,
            Nucleotide::Thymine => Nucleotide::Adenine,
            Nucleotide::Cytosine => Nucleotide::Guanine,
            Nucleotide::Guanine => Nucleotide::Cytosine,
        }
    }
}

impl std::fmt::Display for Nucleotide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_char())
    }
}

/// A codon is a triplet of nucleotides
///
/// Codons are the basic units of the genetic code,
/// specifying amino acids during protein synthesis.
pub type Codon = [Nucleotide; 3];

/// Parse a DNA string into a sequence of nucleotides
///
/// # Arguments
/// * `sequence` - A string containing DNA characters (A, T, C, G)
///
/// # Returns
/// A vector of Nucleotide objects, or None if invalid characters are present
///
/// # Examples
/// ```
/// use rust_math::dna::parse_dna;
///
/// let dna = parse_dna("ATCG").unwrap();
/// assert_eq!(dna.len(), 4);
/// ```
pub fn parse_dna(sequence: &str) -> Option<Vec<Nucleotide>> {
    sequence
        .chars()
        .map(Nucleotide::from_char)
        .collect()
}

/// Convert a sequence of nucleotides to a string
pub fn dna_to_string(sequence: &[Nucleotide]) -> String {
    sequence.iter().map(|n| n.to_char()).collect()
}

/// DNA sequence analysis based on Perez's golden ratio discoveries
///
/// This struct contains the counts of each nucleotide and provides
/// methods for calculating the ratios discovered by Perez.
#[derive(Debug, Clone, PartialEq)]
pub struct DnaAnalysis {
    /// Count of Adenine (A)
    pub count_a: usize,
    /// Count of Thymine (T)
    pub count_t: usize,
    /// Count of Cytosine (C)
    pub count_c: usize,
    /// Count of Guanine (G)
    pub count_g: usize,
    /// Total count of nucleotides
    pub total: usize,
}

impl DnaAnalysis {
    /// Analyze a DNA sequence
    ///
    /// # Arguments
    /// * `sequence` - A slice of nucleotides
    ///
    /// # Returns
    /// A DnaAnalysis struct with counts and ratios
    ///
    /// # Examples
    /// ```
    /// use rust_math::dna::{DnaAnalysis, Nucleotide};
    ///
    /// let sequence = vec![
    ///     Nucleotide::Adenine, Nucleotide::Thymine,
    ///     Nucleotide::Cytosine, Nucleotide::Guanine,
    /// ];
    /// let analysis = DnaAnalysis::from_sequence(&sequence);
    /// ```
    pub fn from_sequence(sequence: &[Nucleotide]) -> Self {
        let mut count_a = 0;
        let mut count_t = 0;
        let mut count_c = 0;
        let mut count_g = 0;

        for &nuc in sequence {
            match nuc {
                Nucleotide::Adenine => count_a += 1,
                Nucleotide::Thymine => count_t += 1,
                Nucleotide::Cytosine => count_c += 1,
                Nucleotide::Guanine => count_g += 1,
            }
        }

        let total = count_a + count_t + count_c + count_g;

        DnaAnalysis {
            count_a,
            count_t,
            count_c,
            count_g,
            total,
        }
    }

    /// Create analysis from counts (constructor helper)
    pub fn from_counts(count_a: usize, count_t: usize, count_c: usize, count_g: usize) -> Self {
        let total = count_a + count_t + count_c + count_g;
        DnaAnalysis {
            count_a,
            count_t,
            count_c,
            count_g,
            total,
        }
    }

    /// Perez's key discovery: (T + A) / (C + G) ratio
    ///
    /// This ratio converges to PEREZ_RATIO = (3 - φ) / 2 ≈ 0.690983
    /// across long DNA sequences and many species.
    ///
    /// # Returns
    /// The actual (T+A)/(C+G) ratio
    pub fn perez_ratio(&self) -> f64 {
        let numerator = (self.count_t + self.count_a) as f64;
        let denominator = (self.count_c + self.count_g) as f64;
        if denominator == 0.0 {
            return 0.0;
        }
        numerator / denominator
    }

    /// Complementary ratio: (C + G) / (T + A)
    ///
    /// Should converge to 1 / PEREZ_RATIO ≈ 1.4472
    pub fn perez_ratio_complementary(&self) -> f64 {
        let ratio = self.perez_ratio();
        if ratio == 0.0 {
            return 0.0;
        }
        1.0 / ratio
    }

    /// Calculate how close the ratio is to the theoretical Perez ratio
    ///
    /// # Returns
    /// A value between 0 and 1, where 1 is perfect match
    pub fn perez_accuracy(&self) -> f64 {
        let actual = self.perez_ratio();
        let theoretical = PEREZ_RATIO;
        if theoretical == 0.0 {
            return 0.0;
        }
        1.0 - (actual - theoretical).abs() / theoretical
    }

    /// Calculate GC content (G + C) / Total
    ///
    /// GC content is an important genomic metric
    pub fn gc_content(&self) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        (self.count_g + self.count_c) as f64 / self.total as f64
    }

    /// Calculate AT content (A + T) / Total
    pub fn at_content(&self) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        (self.count_a + self.count_t) as f64 / self.total as f64
    }

    /// Check A/T ratio (Chargaff's first rule component)
    pub fn ratio_at(&self) -> Option<f64> {
        if self.count_t == 0 {
            return None;
        }
        Some(self.count_a as f64 / self.count_t as f64)
    }

    /// Check C/G ratio (Chargaff's first rule component)
    pub fn ratio_cg(&self) -> Option<f64> {
        if self.count_g == 0 {
            return None;
        }
        Some(self.count_c as f64 / self.count_g as f64)
    }

    /// Check if the sequence respects Chargaff's second rule
    ///
    /// Chargaff's second rule: For single-stranded DNA,
    /// A/T ≈ 1 and C/G ≈ 1 for long sequences.
    ///
    /// # Arguments
    /// * `tolerance` - Maximum allowed deviation from 1.0 (default: 0.05)
    pub fn respects_chargaff_second_rule(&self, tolerance: Option<f64>) -> bool {
        let tol = tolerance.unwrap_or(0.05);

        let ratio_at = match self.ratio_at() {
            Some(r) => r,
            None => return false,
        };

        let ratio_cg = match self.ratio_cg() {
            Some(r) => r,
            None => return false,
        };

        (ratio_at - 1.0).abs() < tol && (ratio_cg - 1.0).abs() < tol
    }

    /// Get nucleotide frequencies as percentages
    pub fn frequencies(&self) -> (f64, f64, f64, f64) {
        if self.total == 0 {
            return (0.0, 0.0, 0.0, 0.0);
        }
        let total_f = self.total as f64;
        (
            self.count_a as f64 / total_f * 100.0,
            self.count_t as f64 / total_f * 100.0,
            self.count_c as f64 / total_f * 100.0,
            self.count_g as f64 / total_f * 100.0,
        )
    }
}

/// Create simulated DNA with golden ratio proportions
///
/// This function creates a DNA sequence that follows Perez's discovered
/// proportions: (T+A)/(C+G) ≈ 0.691
///
/// # Arguments
/// * `length` - The desired length of the DNA sequence
/// * `seed` - Optional seed for reproducibility
pub fn create_golden_ratio_dna(length: usize, seed: Option<u64>) -> Vec<Nucleotide> {
    let mut dna = Vec::with_capacity(length);

    // Adjusted for Perez ratio: (T+A)/(C+G) = 0.691
    // Let C+G = 1, then T+A = 0.691, total = 1.691
    let p_ta = 0.691 / 1.691;  // ≈ 0.409
    let p_cg = 1.0 / 1.691;    // ≈ 0.591

    // Distribute within each group
    let prob_a = p_ta / 2.0;   // ≈ 0.204
    let prob_t = p_ta / 2.0;   // ≈ 0.204
    let prob_c = p_cg / 2.0;   // ≈ 0.296
    let prob_g = p_cg / 2.0;   // ≈ 0.296

    let mut rng_seed = seed.unwrap_or(42);

    for _ in 0..length {
        // Simple LCG for deterministic randomness
        rng_seed = rng_seed.wrapping_mul(1103515245).wrapping_add(12345);
        let r = (rng_seed as f64) / (u64::MAX as f64);

        let nuc = if r < prob_a {
            Nucleotide::Adenine
        } else if r < prob_a + prob_t {
            Nucleotide::Thymine
        } else if r < prob_a + prob_t + prob_c {
            Nucleotide::Cytosine
        } else {
            Nucleotide::Guanine
        };

        dna.push(nuc);
    }

    dna
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nucleotide_from_char() {
        assert_eq!(Nucleotide::from_char('A'), Some(Nucleotide::Adenine));
        assert_eq!(Nucleotide::from_char('t'), Some(Nucleotide::Thymine));
        assert_eq!(Nucleotide::from_char('X'), None);
    }

    #[test]
    fn test_nucleotide_complement() {
        assert_eq!(Nucleotide::Adenine.complement(), Nucleotide::Thymine);
        assert_eq!(Nucleotide::Thymine.complement(), Nucleotide::Adenine);
        assert_eq!(Nucleotide::Cytosine.complement(), Nucleotide::Guanine);
        assert_eq!(Nucleotide::Guanine.complement(), Nucleotide::Cytosine);
    }

    #[test]
    fn test_nucleotide_classifications() {
        assert!(Nucleotide::Adenine.is_purine());
        assert!(Nucleotide::Guanine.is_purine());
        assert!(Nucleotide::Cytosine.is_pyrimidine());
        assert!(Nucleotide::Thymine.is_pyrimidine());
        assert!(Nucleotide::Guanine.is_strong());
        assert!(Nucleotide::Cytosine.is_strong());
    }

    #[test]
    fn test_parse_dna() {
        assert!(parse_dna("ATCG").is_some());
        assert!(parse_dna("atcg").is_some());
        assert!(parse_dna("ATXG").is_none());
    }

    #[test]
    fn test_dna_analysis() {
        let sequence = vec![
            Nucleotide::Adenine, Nucleotide::Thymine,
            Nucleotide::Cytosine, Nucleotide::Guanine,
            Nucleotide::Adenine, Nucleotide::Thymine,
        ];
        let analysis = DnaAnalysis::from_sequence(&sequence);

        assert_eq!(analysis.count_a, 2);
        assert_eq!(analysis.count_t, 2);
        assert_eq!(analysis.count_c, 1);
        assert_eq!(analysis.count_g, 1);
        assert_eq!(analysis.total, 6);
    }

    #[test]
    fn test_perez_ratio() {
        // Create DNA with known proportions
        let dna = vec![
            Nucleotide::Thymine; 691  // 691 T
        ];
        let dna_extended: Vec<Nucleotide> = dna.into_iter()
            .chain(std::iter::repeat(Nucleotide::Cytosine).take(1000))
            .collect();

        let analysis = DnaAnalysis::from_sequence(&dna_extended);
        let ratio = analysis.perez_ratio();

        // Should be close to 0.691
        assert!((ratio - 0.691).abs() < 0.01);
    }

    #[test]
    fn test_gc_content() {
        let analysis = DnaAnalysis::from_counts(25, 25, 25, 25);
        assert_eq!(analysis.gc_content(), 0.5);
        assert_eq!(analysis.at_content(), 0.5);
    }

    #[test]
    fn test_chargaff_second_rule() {
        // Balanced DNA should respect Chargaff's rule
        let analysis = DnaAnalysis::from_counts(250, 250, 250, 250);
        assert!(analysis.respects_chargaff_second_rule(None));
    }

    #[test]
    fn test_dna_to_string_roundtrip() {
        let original = "ATCGATCG";
        let parsed = parse_dna(original).unwrap();
        let converted = dna_to_string(&parsed);
        assert_eq!(converted, original);
    }
}
