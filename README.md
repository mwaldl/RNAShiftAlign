# RNAShiftAlign

**Sankoff-Style RNA Bi-Alignments with Shift Detection**

RNAShiftAlign implements the simultaneous bi-alignment and RNA folding algorithm
for detecting incongruent evolution between sequence and secondary structure in
RNA molecules.

## Overview

This tool addresses the biological phenomenon where RNA sequences evolve with
independent selection pressures on sequence elements and secondary structure,
resulting in:
- Structurally analogous base pairs involving non-homologous sequence positions
- Homologous sequence positions that don't preserve base pairs

The algorithm computes coupled sequence and structure alignments (bi-alignments)
that can represent local shifts between sequence and structure conservation.

## Algorithm

Based on the paper "Incongruences Between Sequence and Secondary Structure
Alignments of Nucleic Acids" by Waldl et al., implementing:
- Extension of Sankoff algorithm for simultaneous alignment and folding
- 4-way alignment representation with shift scoring
- Heuristic optimizations for practical runtime (O(n²) with δ_max constraint)

## Development Status

🚧 **In active development** 🚧

This is a ground-up implementation based on LocARNA, developed in a test-driven manner.

## Building

Dependencies:
- LocARNA library (>= 2.0.0)
- ViennaRNA package (>= 2.5.1)
- C++14 compiler

Build instructions coming soon.

## Authors

Maria Waldl and collaborators
