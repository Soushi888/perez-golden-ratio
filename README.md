# Perez Golden Ratio 🧬

Rust implementation of Jean-Claude Perez's discoveries of golden ratio (φ) patterns in DNA codon populations and genetic code structure.

## Overview

This library demonstrates the mathematical relationships between the golden ratio φ ≈ 1.618 and biological systems, specifically:

- **Perez Ratio**: (T+A)/(C+G) converges to (3-φ)/2 ≈ 0.691 in DNA
- **Fractal Attractors**: Codon frequencies cluster around values related to φ
- **Fibonacci & Lucas**: Number sequences intrinsically connected to φ
- **Pascal's Triangle**: Shows golden ratio patterns and the Perez Hourglass structure

## Mathematical Background

### The Golden Ratio φ

```
φ = (1 + √5) / 2 ≈ 1.618033988749895
φ² = φ + 1 ≈ 2.618033988749895
1/φ ≈ 0.618033988749895
```

### Perez's Key Discovery

Jean-Claude Perez discovered that in single-stranded whole genome DNA:

```
(T + A) / (C + G) = (3 - φ) / 2 ≈ 0.690983005625
```

This ratio appears across species and genomic scales, suggesting a fundamental mathematical structure in genetic code.

### Fractal Attractors

Codon populations cluster around two attractors:

- **Attractor 1** ≈ 1/φ ≈ 0.618 (high frequency)
- **Attractor 2** ≈ 1 - 1/φ ≈ 0.382 (low frequency)
- **Ratio** = Attractor₁ / Attractor₂ = φ ≈ 1.618

## Features

- **DNA Analysis**: Calculate nucleotide ratios and validate golden ratio patterns
- **Fibonacci/Lucas Sequences**: Using Binet's formula with φ
- **Fractal Attractors**: Detect and validate golden ratio clustering
- **Pascal's Triangle**: Generate triangle, mod-2 patterns, and Perez Hourglass
- **Lichtenberg Sequence**: OEIS A000975 connected to hourglass structure

## Installation

```bash
git clone <repo-url>
cd perez-golden-ratio
cargo build --release
```

## Usage

### Run the Demo

```bash
cargo run
```

This displays:
- Golden ratio constants
- Fibonacci and Lucas sequences
- DNA analysis with simulated golden ratio proportions
- Fractal attractor detection
- Pascal's triangle and Perez Hourglass

### Library Usage

```rust
use rust_math::golden_ratio::{PHI, PHI_INVERSE, PEREZ_RATIO};
use rust_math::fibonacci::{fibonacci, lucas};
use rust_math::dna::{DnaAnalysis, Nucleotide};
use rust_math::attractors::FractalAttractor;

// DNA analysis
let dna = vec![
    Nucleotide::Adenine, Nucleotide::Thymine,
    Nucleotide::Cytosine, Nucleotide::Guanine,
];
let analysis = DnaAnalysis::from_sequence(&dna);
println!("Perez ratio: {}", analysis.perez_ratio());

// Fibonacci numbers
let fib_20 = fibonacci(20); // 6765

// Fractal attractors
let frequencies = vec![0.6, 0.62, 0.4, 0.38, 0.59, 0.41];
let attractor = FractalAttractor::from_codon_frequencies(&frequencies);
if attractor.validates_golden_ratio(None) {
    println!("Golden ratio pattern detected!");
}
```

### Run Tests

```bash
cargo test
```

All 37 tests validate the mathematical correctness of the implementation.

## Module Structure

```
src/
├── main.rs              # Demo program
├── golden_ratio.rs      # φ constants and functions
├── fibonacci.rs         # Fibonacci & Lucas sequences
├── dna.rs              # DNA nucleotide analysis
├── attractors.rs       # Fractal attractor detection
└── hourglass.rs        # Pascal's triangle & Perez Hourglass
```

## Mathematical Formulas

### Binet's Formula for Fibonacci

```
F(n) = (φ^n - (-φ)^(-n)) / √5
```

### Lucas Numbers

```
L(n) = φ^n + (-φ)^(-n)
```

### Relationship

```
L(n) = F(n-1) + F(n+1)
```

### Lichtenberg Sequence (OEIS A000975)

```
a(n) = 2 × a(n-1)     for odd n
a(n) = 2 × a(n-1) + 1 for even n
```

## References

### Jean-Claude Perez's Research

1. Perez, J.C. (1991). "Chaos DNA and neuro-computers: a golden link"
2. Perez, J.C. (1997). "L'ADN Décrypté"
3. Perez, J.C. (2010). "Codon Populations in Single-Stranded Whole Human Genome DNA Are Fractal and Fine-Tuned by the Golden Ratio 1.618"

### Key Concepts

- **Golden Ratio**: φ = (1 + √5) / 2
- **Chargaff's Second Rule**: A/T ≈ 1 and C/G ≈ 1 in single-stranded DNA
- **Lucas's Theorem**: C(n,k) is odd iff (k & (n-k)) == 0
- **Sierpiński Triangle**: Pascal's triangle modulo 2

## License

This is a mathematical demonstration project. Refer to original sources for scientific validation.

## Contributing

This implementation serves educational purposes. For scientific research, please refer to Jean-Claude Perez's original publications and peer-reviewed literature.

---

**Note**: This is a mathematical implementation based on published research. The golden ratio patterns in DNA remain an active area of scientific investigation.
