# Contributing to RNAShiftAlign

Thank you for your interest in contributing to RNAShiftAlign! This document provides guidelines for contributing to the project.

## Git Workflow

We use a simplified Git Flow branching model for development.

### Branch Structure

- **`main`** - Production-ready, stable code
  - Always deployable
  - Tagged with version numbers (v0.1.0, v0.2.0, etc.)
  - Only merge from `develop` when ready for release
  - Protected branch (may require pull requests in the future)

- **`develop`** - Integration branch for ongoing development
  - Default branch for development
  - Should always compile and pass tests
  - Features merge here first before going to `main`

- **`feature/*`** - Individual feature development
  - Branch from: `develop`
  - Merge back to: `develop`
  - Naming: `feature/descriptive-name` (e.g., `feature/shift-scoring`)
  - Delete after successful merge

- **`bugfix/*`** - Non-critical bug fixes
  - Branch from: `develop`
  - Merge back to: `develop`
  - Naming: `bugfix/issue-description`

- **`hotfix/*`** - Critical production fixes
  - Branch from: `main`
  - Merge back to: **both** `main` AND `develop`
  - Naming: `hotfix/critical-issue`
  - Used rarely, for urgent fixes

### Workflow Examples

#### Starting a New Feature

```bash
# Make sure develop is up to date
git checkout develop
git pull origin develop

# Create feature branch
git checkout -b feature/shift-penalty-calculator

# Work on feature, commit frequently
git add src/RNAShiftAlign/shift_scoring.hh
git commit -m "Implement shift penalty lookup table"

# Keep feature branch updated with develop
git fetch origin
git merge origin/develop

# When feature is complete and tested
git checkout develop
git merge --no-ff feature/shift-penalty-calculator
git branch -d feature/shift-penalty-calculator
git push origin develop
```

The `--no-ff` flag preserves feature branch history.

#### Preparing a Release

```bash
# Ensure develop is stable, all tests pass
git checkout develop
make clean && make check

# Update version numbers and CHANGELOG
# Edit configure.ac, CHANGELOG.md
git commit -m "Bump version to 0.2.0"

# Merge to main
git checkout main
git merge --no-ff develop

# Tag the release
git tag -a v0.2.0 -m "Release version 0.2.0: Shift scoring implementation"
git push origin main --tags

# Return to develop
git checkout develop
```

#### Hotfix for Critical Bug

```bash
# Branch from main
git checkout main
git checkout -b hotfix/memory-leak

# Fix the bug
# ... make changes ...
git commit -m "Fix critical memory leak in ShiftMatrixM"

# Merge to main
git checkout main
git merge --no-ff hotfix/memory-leak
git tag -a v0.1.1 -m "Hotfix v0.1.1: Fix memory leak"

# Also merge to develop
git checkout develop
git merge --no-ff hotfix/memory-leak

# Clean up
git branch -d hotfix/memory-leak
git push origin main develop --tags
```

## Coding Standards

### C++ Style
- **Standard**: C++14
- **Indentation**: 4 spaces (no tabs)
- **Naming**:
  - Classes: `PascalCase` (e.g., `ShiftMatrixM`)
  - Functions/methods: `snake_case` (e.g., `get_maxshift()`)
  - Variables: `snake_case_` with trailing underscore for members (e.g., `maxshift_`)
  - Constants: `UPPER_CASE`
- **Comments**: Doxygen-style documentation for all public methods

### Documentation
- Every public class and method must have Doxygen documentation
- Include `@brief`, `@param`, `@return`, `@throws` as appropriate
- Example:
  ```cpp
  /**
   * @brief Compute shift penalty for gap pattern mismatch
   *
   * @param c_U Gap pattern in sequence alignment U
   * @param c_V Gap pattern in structure alignment V
   * @param delta Shift penalty parameter
   *
   * @return Penalty score (0, Δ, or 2Δ)
   */
  ```

### Testing
- **Required**: Every new feature must have unit tests
- **Framework**: Catch2
- **Location**: `src/Tests/test_*.cc`
- **Coverage**: Aim for >90% code coverage
- **Test-driven development encouraged**: Write tests first, then implement

### Commit Messages
Follow conventional commit format:

```
<type>: <subject>

<body>

<footer>
```

**Types**:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `test`: Adding or updating tests
- `refactor`: Code refactoring
- `perf`: Performance improvements
- `build`: Build system changes
- `ci`: CI/CD changes

**Examples**:
```
feat: Implement shift penalty calculator

Add compute_shift_penalty() function with lookup table based on
equation 10 from the paper. Includes comprehensive tests for all
15x15 gap pattern combinations.

Closes #42
```

```
fix: Correct memory allocation in ShiftMatrixM constructor

Fixed critical bug where offset dimensions were multiplied by 2
instead of being squared (2*maxshift+1)^2. This was causing
memory corruption for maxshift > 1.

Fixes #38
```

## Pull Request Process

1. **Branch from `develop`** for new features/bugfixes
2. **Write tests** that cover your changes
3. **Ensure all tests pass**: `make check`
4. **Update documentation** if you've changed interfaces
5. **Update CHANGELOG.md** under `[Unreleased]` section
6. **Create pull request** with clear description
7. **Address review comments** if any
8. **Squash commits** if requested before merge

## Development Setup

### Requirements
- C++14 compiler (g++ or clang++)
- Autotools (autoconf, automake)
- LocARNA library (>= 2.0.0)
- ViennaRNA package (>= 2.5.1)
- Conda environment recommended

### Building from Source

```bash
# Clone repository
git clone <repository-url>
cd RNAShiftAlign

# Activate conda environment
conda activate shift

# Generate build system
export ACLOCAL_PATH=$CONDA_PREFIX/share/aclocal/
autoreconf -i

# Configure
mkdir _build && cd _build
../configure --prefix="$PWD/../_inst" \
  PKG_CONFIG_PATH="<path-to-locarna>/_inst/lib/pkgconfig"

# Build and test
make
make check

# Install (optional)
make install
```

## Questions?

If you have questions about contributing, please:
- Open an issue on GitHub
- Email: maria@bioinf.uni-leipzig.de

## Code of Conduct

Be respectful, constructive, and professional in all interactions.

---

Thank you for contributing to RNAShiftAlign!
