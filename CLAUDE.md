# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

`irurueta-numerical` is a Java 21 Maven library providing numerical utilities: root finding, single/multi-variable
optimization, quadrature integration, interpolation, polynomial fitting/estimation, robust (RANSAC-family)
estimation, curve fitting, and Kalman-filter based signal processing. It is part of the `com.irurueta` family of
libraries and depends on `irurueta-algebra`, `irurueta-statistics`, and `irurueta-sorting` (also
`com.irurueta` artifacts, likely siblings checked out separately).

## Build and test commands

Build system is Maven (`pom.xml`), no wrapper script — use the system `mvn`.

```bash
# Full build (compile, test, package, jacoco report, javadoc/source jars)
mvn clean install

# Skip the extras profile (source/javadoc jars) — matches CI's fast path
mvn clean jacoco:prepare-agent install jacoco:report javadoc:jar source:jar -P '!extras'

# Run the whole test suite
mvn test

# Run a single test class
mvn test -Dtest=BrentSingleRootEstimatorTest

# Run a single test method
mvn test -Dtest=BrentSingleRootEstimatorTest#testEstimate

# Generate JaCoCo coverage report (target/site/jacoco/jacoco.csv, index.html)
mvn clean jacoco:prepare-agent test jacoco:report

# Static analysis (Checkstyle, PMD, SpotBugs)
mvn clean compile checkstyle:checkstyle pmd:pmd spotbugs:spotbugs

# Full Maven site (javadoc, surefire report, jacoco, checkstyle, pmd, jxr)
mvn site -Djacoco.skip -DskipTests -P '!extras'
```

Notes:
- Java source/target level is 21 (`maven.compiler.source/target`).
- Checkstyle rules live in `checkstyle.xml` at the repo root (tabs forbidden, 120-char line limit, mandatory
  package-info.java per package, Javadoc enforcement, etc.) — run it before committing non-trivial changes.
- Tests use JUnit Jupiter (5.x) via `maven-surefire-plugin`; there is no separate integration-test source set,
  `maven-failsafe-plugin` is configured but no `*IT.java` classes currently exist.
- A `validate`-phase Groovy step (`groovy-maven-plugin`) regenerates
  `src/main/resources/com/irurueta/numerical/build-info.properties` on every build (build timestamp, group/artifact/
  version, and CI commit/branch info pulled from Jenkins/Travis/GitLab env vars). This is expected, generated output —
  don't hand-edit it, and don't be surprised if it shows as a diff after a local build.
- CI (`.github/workflows/develop.yml`, `master.yml`) runs on `develop` (deploys SNAPSHOT) and `master` (release),
  reporting to SonarCloud and publishing the Maven site to GitHub Pages plus artifacts to Maven Central.

## Architecture

### Package map (`com.irurueta.numerical.*`)

- **root package** — shared abstractions used across the whole library: function evaluator listener interfaces
  (`SingleDimensionFunctionEvaluatorListener`, `MultiDimensionFunctionEvaluatorListener`,
  `MultiVariateFunctionEvaluatorListener`), numerical derivative/gradient/Jacobian estimators (plain, symmetric,
  Savitzky-Golay), polynomial evaluators, and the common exception hierarchy.
- **`roots`** — single-variable root finders (bisection, false position, secant, Newton-Raphson, safe Newton-Raphson,
  Brent, Ridder) plus closed-form polynomial root estimators for 1st/2nd/3rd degree and a general
  `LaguerrePolynomialRootsEstimator`.
- **`optimization`** — single-variable bracketing optimizers (golden section, Brent, derivative Brent) and
  multi-variable optimizers (Powell, conjugate gradient with/without derivatives, quasi-Newton, simplex/Nelder-Mead).
- **`integration`** — quadrature-based numerical integrators. Organized as a matrix of {quadrature rule} x
  {integration method}: quadrature rules (trapezoidal, mid-point, exponential mid-point, upper/lower square-root
  mid-point, double-exponential) each combined with an integration driver (plain quadrature, Simpson, Romberg).
  Every rule/integrator has a scalar variant and a `Matrix*` variant that integrates matrix-valued functions.
- **`interpolation`** — 1D/2D interpolators (linear, cubic spline, bicubic spline, barycentric rational, polynomial,
  rational, Shepard, Kriging) and radial basis functions (Gaussian, multiquadric, inverse multiquadric, thin-plate).
- **`polynomials`** — the core `Polynomial` value type (arithmetic, roots, derivatives, integrals).
- **`polynomials/estimators`** — fits polynomials to data. Split into a plain least-squares estimator
  (`LMSEPolynomialEstimator`, with a weighted variant) and a robust-estimator family
  (`RANSAC`/`LMedS`/`MSAC`/`PROSAC`/`PROMedS` `PolynomialRobustEstimator` subclasses) that reuses the generic
  robust-estimation framework from the `robust` package. Evaluations fed to these estimators are typed via
  `PolynomialEvaluation` subclasses (direct value, derivative, integral, integral-over-interval).
- **`robust`** — the generic, domain-agnostic robust-estimation framework (RANSAC, LMedS, MSAC, PROSAC, PROMedS)
  built around `RobustEstimator<T>` + a listener interface per algorithm; other packages (e.g.
  `polynomials/estimators`) specialize it rather than reimplementing sampling/consensus logic. Includes subset
  selection utilities (`SubsetSelector`, `FastRandomSubsetSelector`) shared by all these estimators.
- **`fitting`** — non-linear/linear curve/surface fitting: Levenberg-Marquardt fitters (single-dimension,
  multi-dimension, multi-variate) and linear fitters (plain and SVD-based), plus a simple straight-line fitter.
- **`signal/processing`** — `KalmanFilter`, 1D convolution (`Convolver1D`), and measurement noise covariance
  estimation.

### Common patterns to know before editing

- **Listener-based function evaluation**: algorithms don't take lambdas/functional callbacks in most APIs (this
  predates widespread lambda usage in the codebase); instead they take a `*FunctionEvaluatorListener` or
  `*EvaluatorListener` instance whose `evaluate(...)` method may throw a checked `EvaluationException`. When adding
  a new algorithm, follow this same listener pattern rather than introducing `java.util.function` interfaces.
- **Estimator lifecycle**: most estimator/fitter classes share a `NotReadyException` / `LockedException` /
  `isReady()` / `isLocked()` protocol — they refuse to run until required inputs are set, and refuse to mutate state
  while an estimation is in progress (some robust estimators use a `volatile boolean locked` guard for thread-safety
  during long-running estimations with progress callbacks).
- **Exception hierarchy is per-package**: each package defines its own root exception (e.g. `FittingException`,
  `IntegrationException`, `InterpolationException`, `PolynomialsException`, `RootEstimationException`,
  `OptimizationException`, `SubsetSelectorException`) rather than sharing one generic exception type — match this
  when adding new failure modes.
- **Naming mirrors the algorithm/method**: class names directly encode the numerical method plus what it operates on
  (e.g. `SimpsonInfinityMidPointQuadratureMatrixIntegrator`, `RombergTrapezoidalQuadratureIntegrator`). When adding a
  new combination (e.g. a new quadrature rule x integration method pair), follow the existing `{Method}{Rule}{Scalar|
  Matrix}Integrator` naming convention rather than inventing a new scheme.
- **Package-info + Javadoc are enforced**: Checkstyle's `JavadocPackage` module requires a `package-info.java` in
  every package, and Javadoc coverage (including private members, per the `maven-javadoc-plugin` `<show>private`
  config) is checked as part of the site build.
