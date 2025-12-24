//! Jean-Claude Perez's DNA Golden Ratio Discoveries (1990s-2020s)
//!
//! This program demonstrates the mathematical discoveries of Jean-Claude Perez
//! regarding golden ratio patterns in DNA, Fibonacci sequences, and fractal
//! attractors in codon populations.
//!
//! # References
//! - Perez, J.C. (1991). "Chaos DNA and neuro-computers: a golden link"
//! - Perez, J.C. (1997). "L'ADN Décrypté"
//! - Perez, J.C. (2010). "Codon Populations in Single-Stranded Whole Human Genome
//!   DNA Are Fractal and Fine-Tuned by the Golden Ratio 1.618"

use num_complex::Complex;

// Module declarations
mod attractors;
mod dna;
mod fibonacci;
mod golden_ratio;
mod hourglass;

// Import commonly used items
use attractors::FractalAttractor;
use dna::{DnaAnalysis, Nucleotide, create_golden_ratio_dna, parse_dna};
use fibonacci::{fibonacci, fibonacci_sequence, lucas};
use golden_ratio::{PEREZ_RATIO, PHI, PHI_INVERSE, PHI_SQUARED};
use hourglass::{lichtenberg_sequence, pascal_triangle, perez_hourglass_simple};

fn main() {
    println!("Jean-Claude Perez's DNA Golden Ratio Discoveries");
    println!("Implementation of research from 1990s-2020s\n");

    // ========================================================================
    // 1. GOLDEN RATIO CONSTANTS
    // ========================================================================
    println!("=== GOLDEN RATIO CONSTANTS ===");
    println!("φ (Phi)                    = {:.15}", PHI);
    println!("φ²                         = {:.15}", PHI_SQUARED);
    println!("1/φ                        = {:.15}", PHI_INVERSE);
    println!("(3-φ)/2 (Perez Ratio)      = {:.15}", PEREZ_RATIO);
    println!();

    // ========================================================================
    // 2. FIBONACCI AND LUCAS NUMBERS
    // ========================================================================
    println!("=== FIBONACCI SEQUENCE ===");
    let fib_seq = fibonacci_sequence(15);
    println!("First 15 Fibonacci numbers:");
    println!("{:?}", fib_seq);
    println!();
    println!("Lucas numbers (n=0 to 10):");
    for i in 0..=10 {
        println!("  L({}) = {:4}", i, lucas(i));
    }
    println!();

    // ========================================================================
    // 3. DEMONSTRATION OF PEREZ RATIO WITH SIMULATED DNA
    // ========================================================================
    println!("=== PEREZ DNA RATIO DISCOVERY ===");
    println!(
        "Theoretical Perez Ratio (T+A)/(C+G) = {:.10}",
        PEREZ_RATIO
    );
    println!();

    // Simulate DNA with golden ratio proportions (Perez's discovery)
    let simulated_dna = create_golden_ratio_dna(10000, Some(42));
    let analysis = DnaAnalysis::from_sequence(&simulated_dna);

    println!("SIMULATED DNA (10,000 bases with φ proportions):");
    println!(
        "  A: {:5}    T: {:5}    C: {:5}    G: {:5}",
        analysis.count_a, analysis.count_t, analysis.count_c, analysis.count_g
    );
    println!();
    println!(
        "ACTUAL (T+A)/(C+G)    = {:.10}",
        analysis.perez_ratio()
    );
    println!(
        "Accuracy to theory    = {:.2}%",
        analysis.perez_accuracy() * 100.0
    );
    println!();
    println!(
        "Complementary (C+G)/(T+A) = {:.10}",
        analysis.perez_ratio_complementary()
    );
    println!(
        "GC Content: {:.1}%",
        analysis.gc_content() * 100.0
    );
    println!(
        "Chargaff's 2nd Rule: {}",
        if analysis.respects_chargaff_second_rule(None) {
            "YES"
        } else {
            "NO"
        }
    );
    println!();

    // ========================================================================
    // 4. FRACTAL ATTRACTOR ANALYSIS
    // ========================================================================
    println!("=== FRACTAL ATTRACTOR ANALYSIS ===");
    println!("Codon frequencies cluster around golden ratio attractors:");
    println!();

    // Create realistic codon frequency simulation
    let n_codons = 64;
    let high_ratio = 0.375; // About 24 out of 64 codons in high cluster
    let codon_freqs = attractors::simulate_codon_frequencies(n_codons, high_ratio, 0.05);
    let attractor = FractalAttractor::from_codon_frequencies(&codon_freqs);

    println!(
        "Analyzed {} codon frequencies:",
        n_codons
    );
    println!();
    println!(
        "Attractor 1 (high)    = {:.6} (≈ 1/φ = {:.6})",
        attractor.attractor_1, PHI_INVERSE
    );
    println!(
        "Attractor 2 (low)     = {:.6} (≈ 1-1/φ = {:.6})",
        attractor.attractor_2,
        1.0 - PHI_INVERSE
    );
    println!(
        "Ratio (A1/A2)         = {:.6} (≈ 2.0)",
        attractor.ratio
    );
    println!(
        "Golden ratio pattern: {}",
        if attractor.validates_golden_ratio(None) {
            "VALIDATED"
        } else {
            "NOT VALIDATED"
        }
    );
    println!();

    // ========================================================================
    // 5. MATHEMATICAL RELATIONSHIPS
    // ========================================================================
    println!("=== KEY MATHEMATICAL RELATIONSHIPS ===");
    println!();
    println!("Fibonacci relationship to φ:");
    println!("  F(n+1)/F(n) → φ as n → ∞");
    println!(
        "  Example: F(20)/F(19) = {}/{} = {:.10}",
        fibonacci(20),
        fibonacci(19),
        fibonacci(20) as f64 / fibonacci(19) as f64
    );
    println!();
    println!("Lucas relationship to φ:");
    println!("  L(n) = F(n-1) + F(n+1)");
    println!(
        "  Example: L(10) = F(9) + F(11) = {} + {} = {}",
        fibonacci(9),
        fibonacci(11),
        fibonacci(9) + fibonacci(11)
    );
    println!();
    println!("Golden ratio powers:");
    println!(
        "  φ¹ = {:.10}",
        PHI
    );
    println!(
        "  φ² = {:.10} (= φ + 1)",
        PHI_SQUARED
    );
    println!(
        "  φ³ = {:.10}",
        PHI * PHI_SQUARED
    );
    println!();
    println!("Lichtenberg sequence (Perez Hourglass connection):");
    let licht_seq = lichtenberg_sequence(8);
    println!("{:?}", licht_seq);
    println!();

    // ========================================================================
    // 6. PASCAL'S TRIANGLE AND HOURGLASS
    // ========================================================================
    println!("=== PASCAL'S TRIANGLE ===");
    let triangle = pascal_triangle(5);
    println!("Pascal's Triangle (rows 0-5):");
    for (i, row) in triangle.iter().enumerate() {
        println!("  Row {}: {:?}", i, row);
    }
    println!();

    println!("=== PEREZ HOURGLASS (simplified visualization) ===");
    println!();
    println!("Note: The actual Perez Hourglass involves recursive");
    println!("parity filtering (mod 2) and infinite-dimensional");
    println!("folding. This is a basic visualization.");
    println!();

    let hourglass = perez_hourglass_simple(4);
    println!("Hourglass structure (n=4):");
    for (i, row) in hourglass.iter().enumerate() {
        println!("  Row {:2}: {:?}", i, row);
    }
}
