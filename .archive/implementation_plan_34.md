# Implementation Plan: Configurable relative and minimum step for derivative/gradient/Jacobian estimators (Issue #34)

## 1. Task summary

Source: GitHub issue #34
Base branch: develop

`JacobianEstimator` and six sibling estimators compute the finite-difference step as `h = EPS * |x|`, falling back to
`EPS` only when the product is exactly `0.0`. For small non-zero parameters (1e-5 … 1e-10) the step collapses
(e.g. 1e-18) and precision is lost; the step is also discontinuous at zero. The change adds two validated, configurable
properties to each estimator, keeping today's behavior by default:

- `relativeStep` — default `EPS`, must be `> 0` (else `IllegalArgumentException`).
- `minStep` — default `0.0`, must be `>= 0` (else `IllegalArgumentException`).
- Step rule: `h = Math.max(relativeStep * Math.abs(x), minStep); if (h == 0.0) h = relativeStep;` followed by the
  existing "actually applied step" correction (`h = xh - x`).

Choices made on the user's behalf (no user input was needed):
- `GradientEstimator` and `DerivativeEstimator` are the base classes of the Symmetric*/SavitzkyGolay* estimators
  (`extends`), so the two properties, validation and a single protected step helper live in the two base classes and
  the four subclasses reuse it instead of repeating the rule (they currently copy `EPS * Math.abs(..)` inline).
  `JacobianEstimator` has no parent and gets its own implementation.
- The "optional" work-array reuse (`fold`, `fh`) in `JacobianEstimator.jacobian(double[], Matrix)` is included, mirroring
  how `xh` is already reused.
- Release `1.7.0` (acceptance criterion 8) is **not** part of this PR's code: it is carried out afterwards from `develop`
  with `/iru-release` (the skill bumps `pom.xml`, `docs/antora.yml`, README mentions, CHANGELOG and what's-new). The
  tasks below only leave the repo release-ready (CHANGELOG `[Unreleased]` and what's-new entries).

Open clarifications: none

## 2. Acceptance criteria coverage

| # | Criterion | Covered by tasks |
|---|---|---|
| 1 | "`JacobianEstimator` exposes `relativeStep` and `minStep`, validated, with defaults `EPS` and `0.0`." | Task 1, Task 1.1, Task 1.2 |
| 2 | "With defaults, all existing tests pass unchanged and results are identical to the current implementation." | Task 1.3, Task 2.3, Task 3.3, Task 4 |
| 3 | "With `relativeStep = 1e-7` and `minStep = 1e-7`, the step for parameter `x` is `1e-7 * max(1, |x|)`." | Task 1.3, Task 2.3, Task 3.3 |
| 4 | "The six sibling estimators listed above expose the same two settings with the same semantics." | Task 2, Task 3 |
| 5 | "New tests cover defaults, validation and the small-parameter accuracy case." | Task 1.3, Task 2.3, Task 3.3 |
| 6 | "Javadoc, Antora docs, `CHANGELOG.md` and what's-new are updated." | Task 1.2, Task 2.1, Task 3.1, Task 4 |
| 7 | "Version `1.7.0` is released." | Task 5 (readiness check and hand-off to `/iru-release`; the release itself is performed after merge) |

## 3. Current code state

All production classes are in `src/main/java/com/irurueta/numerical/`; tests in `src/test/java/com/irurueta/numerical/`.

- `JacobianEstimator` — no superclass. `public static final double EPS = 1e-8`; `private final
  MultiVariateFunctionEvaluatorListener listener`; reusable `private double[] xh`. `jacobian(double[], Matrix)`
  allocates `fold` and `fh` on every call and computes `h = EPS * Math.abs(temp)`, `if (h == 0.0) h = EPS`, then
  `xh[j] = temp + h; h = xh[j] - temp;`. Single constructor `JacobianEstimator(listener)`.
- `GradientEstimator` — `EPS = 1e-8`, `protected final listener`, private `xh`; `gradient(double[], double[])` uses the
  same rule. Subclasses: `SymmetricGradientEstimator` (inline rule, `h1/h2/hh` correction after) and
  `SavitzkyGolayGradientEstimator` (inline rule `var h = EPS * Math.abs(temp)` with a "trick to reduce finite-precision
  error" comment).
- `DerivativeEstimator` — `EPS = 1e-8`, `protected final listener`; `derivative(double x)` uses the same rule.
  Subclasses: `SymmetricDerivativeEstimator` and `SavitzkyGolayDerivativeEstimator` (both inline rule at lines ~47 / ~65).
- Existing tests: `JacobianEstimatorTest`, `GradientEstimatorTest`, `SymmetricGradientEstimatorTest`,
  `SavitzkyGolayGradientEstimatorTest`, `DerivativeEstimatorTest`, `SymmetricDerivativeEstimatorTest`,
  `SavitzkyGolayDerivativeEstimatorTest`.
- Docs: Antora component in `docs/` (`antora.yml`, version `1.6.1`); `docs/modules/ROOT/pages/algorithms.adoc` is the
  algorithm catalog index and `docs/modules/ROOT/pages/algorithms/*.adoc` hold one page per algorithm; no page currently
  covers derivative/gradient/Jacobian estimation (the estimators are only mentioned in the Levenberg-Marquardt and
  conjugate-gradient pages). `CHANGELOG.md` has an empty `## [Unreleased]`; `docs/modules/ROOT/pages/whats-new.adoc`
  starts at `== 1.6.1`.
- Stack: Java 21, Maven, JUnit Jupiter 5, Checkstyle (120 cols, Javadoc on everything incl. private, no tabs).

## 4. Implementation steps

### Group 1 — Estimator changes and their tests (Parallelizable: yes — each task edits disjoint main and test files)

- [x] Task 1. Add `relativeStep` / `minStep` to `JacobianEstimator` _(java)_ — implemented, tests added; JacobianEstimatorTest (20 tests), 92.2% line / 86.4% branch coverage; PMD clean, remaining Checkstyle HiddenField on new ctor/setter params and SpotBugs CT_CONSTRUCTOR_THROW (validating ctor) left as inherent to the specified API.
  - [x] Task 1.1. In `JacobianEstimator.java` add fields `relativeStep = EPS`, `minStep = 0.0`; getters
    `getRelativeStep()`/`getMinStep()`; setters `setRelativeStep(double)` (throws `IllegalArgumentException` if
    `<= 0`) and `setMinStep(double)` (throws `IllegalArgumentException` if `< 0`; NaN rejected too); constructor
    overload `JacobianEstimator(listener, relativeStep, minStep)` that validates through the setters' rules; keep the
    existing constructor delegating with the defaults. In `jacobian(double[], Matrix)` use
    `h = Math.max(relativeStep * Math.abs(temp), minStep); if (h == 0.0) { h = relativeStep; }` and keep
    `xh[j] = temp + h; h = xh[j] - temp;`. Also reuse `fold`/`fh` as private fields reallocated only when
    `numvars` changes (same pattern as `xh`).
  - [x] Task 1.2. Javadoc: class doc and method docs describe the step rule, the defaults, and that an absolute minimum
    step assumes parameters scaled to order 1; keep public `EPS` documented as the default relative step.
  - [x] Task 1.3. Extend `JacobianEstimatorTest`: defaults (`EPS`, `0.0`); setters/constructor reject
    `relativeStep <= 0` and `minStep < 0`; defaults produce results bit-for-bit equal to a reference computation using
    the original rule; `relativeStep = minStep = 1e-7` yields step `1e-7 * max(1, |x|)` (verify via a listener that
    records the perturbation); small-parameter case (x in 1e-5…1e-10) where `minStep = 1e-7` is closer to the analytical
    Jacobian than defaults; a parameter at exactly `0.0` still gets a non-zero step with defaults.

- [x] Task 2. Add the two settings to the gradient estimator family _(java)_ — implemented, tests added; Gradient tests (13 each), coverage 97.6%/96.3%/90.4% line; ExplicitInitialization fixed; remaining HiddenField and CT_CONSTRUCTOR_THROW left as inherent to the specified API.
  - [x] Task 2.1. In `GradientEstimator.java` add the same two properties (fields, getters, validated setters,
    constructor `GradientEstimator(listener, relativeStep, minStep)`, existing constructor kept) and a
    `protected double computeStep(double x)` helper returning `Math.max(relativeStep * Math.abs(x), minStep)` with the
    `== 0.0 → relativeStep` fallback; use it in `gradient(double[], double[])`. Add overload constructors with the
    new parameters to `SymmetricGradientEstimator` and `SavitzkyGolayGradientEstimator` (delegating to `super`).
    Javadoc everything, including the order-1 scaling caveat.
  - [x] Task 2.2. In `SymmetricGradientEstimator.java` and `SavitzkyGolayGradientEstimator.java` replace the inline
    `EPS * Math.abs(temp)` / `h == 0.0` blocks with `computeStep(temp)`, keeping the existing `h1/h2/hh` and
    Savitzky-Golay logic untouched.
  - [x] Task 2.3. Extend `GradientEstimatorTest`, `SymmetricGradientEstimatorTest`, `SavitzkyGolayGradientEstimatorTest`
    with: defaults, validation (setters and new constructors), step `1e-7 * max(1, |x|)` with `1e-7/1e-7`, unchanged
    results with defaults, zero-valued parameter still perturbed, and a small-parameter (1e-5…1e-10) accuracy case
    where `minStep = 1e-7` beats the defaults. For Savitzky-Golay, pick a function/tolerance so the accuracy
    comparison is meaningful for its parabolic fit (it may be less accurate for large steps — assert relative
    closeness to the analytical gradient between the two settings).

- [x] Task 3. Add the two settings to the derivative estimator family _(java)_ — implemented, tests added; Derivative tests (11 each), coverage 100%/100%/94.4% line; remaining HiddenField and CT_CONSTRUCTOR_THROW left as inherent to the specified API.
  - [x] Task 3.1. In `DerivativeEstimator.java` add the same properties, validated setters, constructor
    `DerivativeEstimator(listener, relativeStep, minStep)` and `protected double computeStep(double x)`; use it in
    `derivative(double)`. Add the overload constructors to `SymmetricDerivativeEstimator` and
    `SavitzkyGolayDerivativeEstimator`. Javadoc everything. — DerivativeEstimator, SymmetricDerivativeEstimator, SavitzkyGolayDerivativeEstimator.
  - [x] Task 3.2. In `SymmetricDerivativeEstimator.java` and `SavitzkyGolayDerivativeEstimator.java` replace the inline
    step computation with `computeStep(x)`. — done.
  - [x] Task 3.3. Extend `DerivativeEstimatorTest`, `SymmetricDerivativeEstimatorTest`,
    `SavitzkyGolayDerivativeEstimatorTest` with the same case set as Task 2.3 for scalar `x`. — three test classes extended.

- Group 1 validation (once, after all tasks): delegate to the `iru-gate-runner` agent to run the seven touched test
  classes (`mvn test -Dtest='JacobianEstimatorTest,GradientEstimatorTest,SymmetricGradientEstimatorTest,SavitzkyGolayGradientEstimatorTest,DerivativeEstimatorTest,SymmetricDerivativeEstimatorTest,SavitzkyGolayDerivativeEstimatorTest'`),
  then license headers, Javadoc audit, coverage (>= 80% on the touched classes) and Checkstyle/PMD/SpotBugs.

### Group 2 — Documentation (Parallelizable: yes — single task; depends on Group 1 API being final)

- [x] Task 4. Update documentation
  - [x] Task 4.1. Add an Antora page `docs/modules/ROOT/pages/algorithms/numerical-derivatives.adoc` (follow the format of
    existing algorithm pages, e.g. `derivative-brent-optimization.adoc`) covering `DerivativeEstimator`,
    `SymmetricDerivativeEstimator`, `SavitzkyGolayDerivativeEstimator`, `GradientEstimator` and variants, and
    `JacobianEstimator`: step rule, `relativeStep`/`minStep` defaults, a code snippet, and guidance on when to set a
    minimum step (parameters small but non-zero, scaled to order 1, e.g. `relativeStep = minStep = 1e-7`). Link it from
    `docs/modules/ROOT/pages/algorithms.adoc` and add it to the nav file under `docs/modules/ROOT/` (`nav.adoc`, check
    its actual location) if one exists. Do not change the Antora component `version`.
  - [x] Task 4.2. Add under `## [Unreleased]` in `CHANGELOG.md` an `### Added` entry describing the new settings and the
    affected classes (#34).
  - [x] Task 4.3. Add a `== 1.7.0` section at the top of `docs/modules/ROOT/pages/whats-new.adoc` with the same entry
    (emoji style matches existing bullets, e.g. `✨`).
  - [x] Task 4.4. Run `/iru-update-docs` and the docs build (`/iru-build-docs`) to verify the new page renders.

### Group 3 — Release readiness (Parallelizable: yes — single task)

- [x] Task 5. Verify release readiness for `1.7.0`
  - [x] Task 5.1. Run the full test suite and static analysis (`mvn clean compile checkstyle:checkstyle pmd:pmd
    spotbugs:spotbugs`) and confirm no regressions and unchanged results with defaults. Passed: 1166 tests, 0 failures; quality warnings accepted as-is by user.
  - [x] Task 5.2. Record in the PR description that the release itself (version bump `1.7.0-SNAPSHOT` → `1.7.0`, moving
    `[Unreleased]` to `[1.7.0]`, Maven Central publish) is done afterwards from `develop` with `/iru-release`.
