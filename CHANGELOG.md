# Changelog

All notable changes to RNAShiftAlign will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned
- Shift penalty scoring function (w_s)
- Basic unpaired bi-alignment algorithm
- ShiftMatrixD implementation for structure matching
- Full Sankoff-style bi-alignment with structure
- Command-line interface
- Heuristic optimizations (δ_max, sparsity, envelope)

---

## [0.1.0] - 2025-11-06

### Added
- **Project Infrastructure**
  - Git repository with `main` and `develop` branches
  - Autotools build system (configure.ac, Makefile.am)
  - Integration with LocARNA library (>= 2.0.0)
  - C++14 compilation with proper dependency management
  - Catch2 testing framework for unit tests
  - Continuous testing with `make check`

- **ShiftMatrixM: 4D Dynamic Programming Matrix**
  - Template class for storing bi-alignment scores M(y1, y2, y3, y4)
  - Enforces shift constraint: |y1-y3| ≤ δ_max and |y2-y4| ≤ δ_max
  - Memory-efficient storage: O(n² · (2δ_max+1)²)
  - Comprehensive bounds and constraint checking
  - Clear exception messages for debugging
  - Full Doxygen-compatible documentation

- **Test Suite**
  - 7 test sections with 20+ individual test cases
  - Tests for basic operations (construction, resize, fill, transform)
  - Tests for index validation and bounds checking
  - Tests for maxshift constraint enforcement
  - Tests with various δ_max values (0, 1, 2)
  - All tests passing ✅

- **Documentation**
  - README.md with project overview
  - CONTRIBUTING.md with development guidelines
  - CHANGELOG.md (this file)
  - Inline code documentation (Doxygen-ready)

### Implementation Details

#### Build System
- Modern autotools configuration with `AC_CONFIG_HEADERS` and `subdir-objects`
- Proper LocARNA library integration via pkg-config
- Clean .gitignore for build artifacts and test binaries
- Separate test binary per component for modularity

#### ShiftMatrixM Design
- **Memory allocation**: Correct sizing with `(2*maxshift+1)²` for offset dimensions
- **Header guard**: Uses `RNASHIFTALIGN_` prefix to avoid namespace conflicts
- **Documentation**: Comprehensive Doxygen documentation for all methods
- **Error handling**: Descriptive exception messages with parameter context
- **Efficient shift calculation**: Direct conditional instead of min/max functions
- **Type safety**: Proper use of unsigned `size_type` without redundant checks
- **Modern C++**: Uses `std::fill()` and other standard algorithms
- **Clean includes**: Only necessary headers (`<stdexcept>`, `<string>`, `<tuple>`)

#### Test Infrastructure
- Proper Catch2 main implementation in dedicated test_main.cc
- One test binary per component for better organization
- Comprehensive test coverage (7 test sections, 20+ test cases)

### Technical Details

#### Memory Layout
ShiftMatrixM uses linearized 4D addressing:
```
index = y1 * bdim * a_offsets_dim * b_offsets_dim
      + y2 * a_offsets_dim * b_offsets_dim
      + (y3 - y1 + maxshift) * b_offsets_dim
      + (y4 - y2 + maxshift)
```

Where:
- `adim = len(sequence_a) + 1`
- `bdim = len(sequence_b) + 1`
- `a_offsets_dim = b_offsets_dim = 2 * maxshift + 1`

#### Constraint Enforcement
The matrix enforces two key constraints:
1. Sequence bounds: `0 ≤ y1, y3 < adim` and `0 ≤ y2, y4 < bdim`
2. Shift bounds: `|y1 - y3| ≤ maxshift` and `|y2 - y4| ≤ maxshift`

These correspond to the δ_max heuristic from Section 4 of the paper.

---

## Development Process

### Git Workflow
- **main**: Stable, release-ready code
- **develop**: Integration branch for active development
- **feature/**: Feature branches (branch from `develop`, merge back to `develop`)
- **bugfix/**: Bug fix branches (branch from `develop`, merge back to `develop`)
- **hotfix/**: Critical fixes (branch from `main`, merge to both `main` and `develop`)

### Testing Philosophy
Test-driven development approach:
1. Write tests first
2. Implement functionality to pass tests
3. Refactor with confidence that tests catch regressions

---

## References

- **Paper**: "Incongruences Between Sequence and Secondary Structure Alignments of Nucleic Acids"
  - Waldl et al. (draft)
- **LocARNA**: Base library for Sankoff-style RNA alignment
  - Used for: Scoring, RnaData, BasePairs, and other utilities
- **Development notes**: See planning_docs/ directory (not in repository)

---

## Contributors

- Maria Waldl <maria@bioinf.uni-leipzig.de>

---

[Unreleased]: https://github.com/mwaldl/RNAShiftAlign/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/mwaldl/RNAShiftAlign/releases/tag/v0.1.0
