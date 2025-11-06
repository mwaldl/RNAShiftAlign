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
  - Inline code documentation (Doxygen-ready)
  - CHANGELOG.md (this file)
  - Planning documents in separate directory

### Changed from RNAbialign

#### Build System Improvements
- **configure.ac**: Use modern `AC_CONFIG_HEADERS` instead of deprecated `AM_CONFIG_HEADER`
- **configure.ac**: Added `subdir-objects` to `AM_INIT_AUTOMAKE` for proper subdirectory handling
- **.gitignore**: Specific test binary names instead of wildcard `test_*` to avoid ignoring test source files
- **.gitignore**: Added `/RNAShiftAlign` path prefix to only ignore binary in root directory
- **Main binary**: Added version string and informative output messages

#### ShiftMatrixM Critical Fixes
- **🐛 CRITICAL BUG FIX**: Constructor memory allocation
  - **Before**: `mat_(adim * bdim * (maxshift * 2 + 1) * 2)` ❌
  - **After**: `mat_(adim * bdim * (maxshift * 2 + 1) * (maxshift * 2 + 1))` ✅
  - **Impact**: Was allocating 2× memory instead of (2δ+1)² for offset dimensions
  - **Result**: Would have caused memory corruption or access violations

#### ShiftMatrixM Code Quality Improvements
- **Header guard**: Changed from `LOCARNA_SHIFTMATRIX_M_HH` to `RNASHIFTALIGN_SHIFTMATRIX_M_HH` to prevent conflicts
- **Documentation**: Full Doxygen documentation for all methods and parameters
- **Error messages**:
  - Added context to exception messages (parameter names and constraint descriptions)
  - Better formatting with spaces for readability
  - Explicit constraint descriptions (e.g., "[y1-maxshift, y1+maxshift]")
- **Shift calculation**:
  - **Before**: `std::max(y1, y3) - std::min(y1, y3) > maxshift_`
  - **After**: `size_type shift_a = (y1 > y3) ? (y1 - y3) : (y3 - y1);`
  - **Benefit**: More efficient, clearer variable naming
- **Redundant checks removed**: Removed `y < 0` checks for `size_type` (unsigned type, always ≥ 0)
- **fill() method**: Use `std::fill()` instead of manual loop (more idiomatic, potentially better optimized)
- **Includes**: Only include headers that are actually used
  - Added: `<stdexcept>`, `<string>`, `<tuple>`
  - Removed: `<assert.h>` (not used)

#### Test Infrastructure
- **test_main.cc**: Proper Catch2 main implementation (was empty in RNAbialign)
- **Test organization**: One test binary per component instead of combined test binary
- **Test coverage**: Significantly expanded from 3 sections to 7 sections

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
  - Waldl et al. (draft in `../paper/shiftSankoff.pdf`)
- **Base implementation**: LocARNA (Sankoff-style RNA alignment)
  - Source: `../LocARNA/`
- **Previous development**: RNAbialign prototype
  - Source: `../RNAbialign/`

---

## Contributors

- Maria Waldl <maria@bioinf.uni-leipzig.de>

---

[Unreleased]: https://github.com/user/RNAShiftAlign/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/user/RNAShiftAlign/releases/tag/v0.1.0
