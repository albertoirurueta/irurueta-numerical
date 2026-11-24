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

### Checkstyle baseline

`mvn checkstyle:checkstyle` now applies the repository's `checkstyle.xml` (the plugin is configured under
`<build><plugins>` in `pom.xml`, not only under `<reporting>`). Measured on 2026-11-24 at version 1.8.0-SNAPSHOT
(`src/main` only, before the fixes of this PR): 1010 violations in total.

| Check | Violations |
|-------|-----------|
| HiddenField | 443 |
| Indentation | 203 |
| MagicNumber | 101 |
| JavadocStyle | 87 |
| RedundantModifier | 67 |
| ImportOrder | 28 |
| ThrowsCount | 14 |
| JavadocVariable | 12 |
| DeclarationOrder | 11 |
| ParameterAssignment | 9 |
| ExplicitInitialization | 9 |
| JavadocPackage | 7 |
| TrailingComment | 4 |
| VisibilityModifier | 4 |
| ParameterNumber | 3 |
| NeedBraces | 2 |
| NewlineAtEndOfFile | 2 |
| WhitespaceAround | 1 |
| ConstantName | 1 |
| MultipleStringLiterals | 1 |
| JavadocMethod | 1 |

After Group 3 (`HiddenField` configured to ignore constructor parameters and setters, intentionally): `HiddenField`
drops from 443 to 113 and the total from 1010 to 680 (re-measured 2026-11-24, `src/main` only).

After Group 4 (mechanical fixes: `JavadocPackage`, `JavadocVariable`, `JavadocMethod`, `NewlineAtEndOfFile`, `NeedBraces`,
`WhitespaceAround`, `ConstantName`, `TrailingComment`, `ExplicitInitialization`, `ImportOrder`, `DeclarationOrder` are all
at 0): the total drops from 680 to 602 (re-measured 2026-11-24, `src/main` only). Remaining categories, left for
follow-up PRs: Indentation 203, HiddenField 113, MagicNumber 101, JavadocStyle 87, RedundantModifier 67, ThrowsCount 14,
ParameterAssignment 9, VisibilityModifier 4, ParameterNumber 3, MultipleStringLiterals 1.

### PMD baseline

`mvn pmd:pmd` (default rules, no custom ruleset, report-only: no `failOnViolation` gate and no `<build>` binding) needs
JDK 21: set `JAVA_HOME` to a JDK 21 (on macOS: `export JAVA_HOME=$(/usr/libexec/java_home -v 21)`). Measured on
2026-11-24 at version 1.8.0-SNAPSHOT (`src/main` only): 35 violations before the fixes of issue #40 (UselessParentheses
23, EmptyCatchBlock 11, UnnecessaryModifier 1), 0 after. Convention: intentionally empty catch blocks name the parameter
`ignored` (PMD's default `allowExceptionNameRegex` is `^(ignored|expected)$`) and keep an explanatory comment; no
suppressions are used.

### SpotBugs baseline

`mvn clean compile spotbugs:spotbugs` (report-only: no `check` goal bound, no `failOnError`; needs JDK 21, set
`JAVA_HOME=$(/usr/libexec/java_home -v 21)` on macOS) reported 264 findings on 2026-11-24 at version 1.8.0-SNAPSHOT
(`src/main` only): CT_CONSTRUCTOR_THROW 118, EI_EXPOSE_REP 75, EI_EXPOSE_REP2 53, REC_CATCH_EXCEPTION 6,
DE_MIGHT_IGNORE 5, FE_FLOATING_POINT_EQUALITY 4, IM_AVERAGE_COMPUTATION_COULD_OVERFLOW 2,
AT_STALE_THREAD_WRITE_OF_PRIMITIVE 1. Issue #41 fixed the real bug (`BaseInterpolator` uses `>>>` for the binary search
midpoint) and excluded the rest through `spotbugs-exclude.xml` (configured on the plugin in `pom.xml`), each entry with
a rationale comment. The design-level types (`CT_CONSTRUCTOR_THROW`, `EI_EXPOSE_REP`, `EI_EXPOSE_REP2`) are excluded
project-wide on purpose, so **new classes are not checked for these three bug types**: review constructor validation
and array/matrix sharing by hand when adding classes. Everything else (the Numerical Recipes exact comparisons, the
"failed subset yields no solution" catch blocks and `RobustEstimator`'s `AT_STALE_THREAD_WRITE_OF_PRIMITIVE`) is
scoped to the reviewed class *and* method, so new occurrences elsewhere are still reported. Robust estimators
and their listeners are single-threaded by contract (documented in their Javadoc), so no `volatile` fields are added
for configuration state. Result: 0 unfiltered findings. No `@SuppressFBWarnings` annotations are used.

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
