# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.6.0] - 2026-07-15

### Added

- Added the Antora documentation site, including a full algorithm catalog with one page per algorithm, an
  architecture overview, and installation instructions.
- Added Claude Code skills to automate common repository workflows (release management, documentation generation,
  code review, and more).

## [1.5.0] - 2026-03-03

### Fixed

- `KalmanFilter`: fixed an aliasing bug in the predict and correct steps where the internal `errorCovPre` and
  `gain` matrices were replaced by reference to an internal temporary instead of having their contents copied,
  which could cause previously computed matrices to be silently mutated on subsequent filter iterations.

## [1.4.0] - 2025-12-17

### Added

- `LevenbergMarquardtMultiDimensionFitter` and `LevenbergMarquardtMultiVariateFitter`: added `getReducedChisq()`
  to expose the reduced chi-square goodness-of-fit value (chi-square divided by its degrees of freedom).

### Fixed

- `LevenbergMarquardtMultiDimensionFitter` and `LevenbergMarquardtMultiVariateFitter`: `getP()` no longer computes
  an invalid chi-square CDF when the degrees of freedom are zero or negative, returning `1.0` in that case instead.

## [1.3.2] - 2025-09-22

### Changed

- Updated dependencies and improved test reliability. No public API or behavior changes.

## [1.3.1] - 2025-09-19

### Changed

- Updated Maven plugins and CI workflows, and improved test reliability. No public API or behavior changes.

## [1.3.0] - 2024-10-26

### Changed

- Migrated the build to Java 17 (source/target compatibility) and migrated the test suite from JUnit 4 to JUnit 5.
- Refactored source code to use local variable type inference (`var`) where applicable.

## [1.2.1] - 2023-11-18

### Fixed

- Corrected Javadoc comments across the interpolation package (`CurveInterpolator`, `KrigingInterpolator`) and
  fixed related test issues.

## [1.2.0] - 2023-11-18

### Added

- New `com.irurueta.numerical.integration` package: quadrature-based numerical integrators combining quadrature
  rules (trapezoidal, mid-point, exponential mid-point, upper/lower square-root mid-point, double-exponential)
  with integration drivers (plain quadrature, Simpson, Romberg), each with scalar and matrix-valued variants.
- New `com.irurueta.numerical.interpolation` package: 1D/2D interpolators (linear, cubic spline, bicubic spline,
  barycentric rational, polynomial, rational, Shepard, Kriging) and radial basis functions (Gaussian, multiquadric,
  inverse multiquadric, thin-plate).

### Changed

- Updated project dependencies.

## [1.1.0] - 2021-12-11

### Added

- Initial public release of the numerical utilities library, including:
  - Single-variable root finders (bisection, false position, secant, Newton-Raphson, safe Newton-Raphson, Brent,
    Ridder) and closed-form/general polynomial root estimators.
  - Single-variable and multi-variable optimizers (golden section, Brent, Powell, conjugate gradient, quasi-Newton,
    simplex).
  - Numerical derivative, gradient, and Jacobian estimators.
  - Polynomial arithmetic and least-squares/robust polynomial estimation.
  - A generic robust-estimation framework (RANSAC, LMedS, MSAC, PROSAC, PROMedS) with subset selection utilities.
  - Curve/surface fitting (Levenberg-Marquardt and linear/SVD-based fitters).
  - Signal processing utilities: `KalmanFilter`, 1D convolution, and measurement noise covariance estimation.

[Unreleased]: https://github.com/albertoirurueta/irurueta-numerical/compare/1.6.0...HEAD
[1.6.0]: https://github.com/albertoirurueta/irurueta-numerical/compare/1.5.0...1.6.0
[1.5.0]: https://github.com/albertoirurueta/irurueta-numerical/compare/1.4.0...1.5.0
[1.4.0]: https://github.com/albertoirurueta/irurueta-numerical/compare/1.3.2...1.4.0
[1.3.2]: https://github.com/albertoirurueta/irurueta-numerical/compare/1.3.1...1.3.2
[1.3.1]: https://github.com/albertoirurueta/irurueta-numerical/compare/1.3.0...1.3.1
[1.3.0]: https://github.com/albertoirurueta/irurueta-numerical/compare/1.2.1...1.3.0
[1.2.1]: https://github.com/albertoirurueta/irurueta-numerical/compare/1.2.0...1.2.1
[1.2.0]: https://github.com/albertoirurueta/irurueta-numerical/compare/1.1.0...1.2.0
[1.1.0]: https://github.com/albertoirurueta/irurueta-numerical/releases/tag/1.1.0
