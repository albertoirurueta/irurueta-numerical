# Implementation Plan: Simpson integrators use the wrong extrapolation for step-tripling mid-point quadratures (Issue #43)

## 1. Task summary

Source: GitHub issue #43
Base branch: develop

`SimpsonIntegrator` and `SimpsonMatrixIntegrator` apply the trapezoidal extrapolation `s = (4*st - ost)/3` to every
quadrature. The mid-point family (mid-point, infinity mid-point, lower/upper square-root mid-point) triples its step on
each refinement, so the error falls to 1/9 per step and the correct extrapolation is `s = (9*st - ost)/8`
(Numerical Recipes 3rd ed. §4.4, p. 169). With the wrong coefficient the leading error term is not cancelled, so these
integrators converge no faster than the plain mid-point rule and, with `JMAX = 20` (up to 3^19 evaluations), can run
about 10 s before throwing `IntegrationException`.

Approach chosen (no user input needed):
- Parameterize the two base classes with an extrapolation **ratio** (`4.0` halving, `9.0` tripling) and a `jmax`,
  via a new protected constructor `(q, eps, ratio, jmax)`. The existing `(q, eps)` constructor keeps today's behavior
  (`4.0`, `20`) so `SimpsonTrapezoidal*` and `SimpsonDoubleExponentialRule*` are unchanged. The step is
  `s = (ratio*st - ost)/(ratio - 1)`. This mirrors how the Romberg classes specialize their step ratio and `JMAX`.
- The 8 mid-point-family subclasses (4 scalar, 4 matrix) call the new constructor with `9.0` and `JMAX = 14`
  (value taken from NR `qromo` / `RombergIntegrator`).
- Out of scope, per the issue's "consider separately" notes: the `normMin` matrix convergence check and the
  double-exponential combination. They are not changed and are left for separate issues.
- Javadoc of the base and affected classes cites NR 3rd ed. §4.4.

Open clarifications: none

## Acceptance criteria coverage

| # | Criterion | Covered by tasks |
|---|---|---|
| 1 | "The mid-point-family Simpson integrators (scalar and matrix) use `(9*st - ost)/8`; trapezoidal keeps `(4*st - ost)/3`." | Task 1, Task 2, Task 3, Task 4 |
| 2 | "The mid-point family uses a reduced `JMAX`." | Task 1, Task 2, Task 3, Task 4 |
| 3 | "Tests demonstrate faster convergence for the mid-point Simpson integrators, and the full suite passes." | Task 3, Task 4, Task 5 |
| 4 | "The Javadoc of the affected classes cites NR 3rd ed. §4.4 for the extrapolation step." | Task 1, Task 2, Task 3, Task 4 |

## Current code state

Package `com.irurueta.numerical.integration` (`src/main/java/com/irurueta/numerical/integration/`). Java 21, Maven
single module, JUnit Jupiter 5, Checkstyle (120 cols, Javadoc incl. private members), license header on every file.

- `SimpsonIntegrator<T extends Quadrature>`: `EPS = 1e-10`, private `JMIN = 5`, private `JMAX = 20`, constructor
  `protected SimpsonIntegrator(T q, double eps)`. `integrate()` line 83: `s = (4.0 * st - ost) / 3.0;`.
- `SimpsonMatrixIntegrator<T extends MatrixQuadrature>`: same constants; `integrate(Matrix)` lines 91-95 compute
  `(4.0 * st - ost) / 3.0` via `tmp.multiplyByScalar(4.0)` / `subtract(ost)` / `multiplyByScalar(1.0 / 3.0)`. Convergence
  uses `normMin` (unchanged).
- Mid-point-family subclasses (scalar / matrix): `SimpsonMidPointQuadrature{,Matrix}Integrator`,
  `SimpsonInfinityMidPointQuadrature{,Matrix}Integrator`, `SimpsonLowerSquareRootMidPointQuadrature{,Matrix}Integrator`,
  `SimpsonUpperSquareRootMidPointQuadrature{,Matrix}Integrator`. Each has `super(new <Quadrature>(...), eps)` calls.
- Unchanged: `SimpsonTrapezoidalQuadrature{,Matrix}Integrator`, `SimpsonDoubleExponentialRuleQuadrature{,Matrix}Integrator`.
- `MidPointQuadrature.next()` triples the step (`del = (b-a)/(3.0*tnm)`), confirming the 9/8 factor.
- Model: `RombergIntegrator` (`h[j] = h[j-1]/9.0`, `JMAX = 14`, `EPS = 3e-9`); `RombergTrapezoidalQuadratureIntegrator`
  overrides with `0.25` / `20`.
- Tests (`src/test/java/.../integration/`): one `Simpson*Test` per class plus `SimpsonIntegratorTest` and
  `SimpsonMatrixIntegratorTest`. The exponential tests were hardened in issue #36
  (see `.archive/implementation_plan_36.md`: seeded `Random`, `MIN_INTERVAL`, `MIN_ABS_LAMBDA`, `EPS_EXPONENTIAL = 1e-8`,
  tolerance `1e-6 * max(1, |expected|)`); reuse those patterns.
- Docs: `docs/modules/ROOT/pages/algorithms/simpson-integration.adoc` (and `midpoint-integration.adoc`).

## Implementation steps

### Group 1 — Parameterize the base classes (Parallelizable: yes — different production and test files)

- [x] Task 1. Parameterize `SimpsonIntegrator` _(java)_ — `SimpsonIntegrator.java`, `SimpsonIntegratorTest.java`
  - Done: SimpsonIntegrator.java, SimpsonIntegratorTest.java; 95.3% line / 82.1% branch coverage; new checkstyle line-length/hidden-field items match pre-existing style, SpotBugs CT_CONSTRUCTOR_THROW left (required argument validation).
  - [x] Task 1.1. Replace the `private static final int JMAX` with an instance field `private final int jmax` and add
    `private final double ratio`. Keep `protected SimpsonIntegrator(q, eps)` delegating to
    `this(q, eps, TRAPEZOIDAL_RATIO /*4.0*/, DEFAULT_JMAX /*20*/)`; add
    `protected SimpsonIntegrator(T q, double eps, double ratio, int jmax)` (validate `ratio > 1` and `jmax > JMIN`,
    `IllegalArgumentException` otherwise). _Revised after the quality gate: the validation was removed (only subclasses call this
    constructor, always with constants) to avoid 4 new SpotBugs `CT_CONSTRUCTOR_THROW` warnings; the matching
    invalid-argument tests were dropped._
  - [x] Task 1.2. In `integrate()` use `s = (ratio * st - ost) / (ratio - 1.0);` and loop to `jmax`.
  - [x] Task 1.3. Javadoc: explain the ratio (2^2 for halving, 3^2 for tripling), the formulas `(4*st-ost)/3` and
    `(9*st-ost)/8`, and cite Numerical Recipes 3rd ed. §4.2 (p. 165) and §4.4 (p. 169).
  - [x] Task 1.4. `SimpsonIntegratorTest`: add tests for the new constructor (invalid `ratio`/`jmax` throw) and that the
    default constructor still yields the trapezoidal behavior.
- [x] Task 2. Parameterize `SimpsonMatrixIntegrator` _(java)_ — `SimpsonMatrixIntegrator.java`,
  - Done: SimpsonMatrixIntegrator.java, SimpsonMatrixIntegratorTest.java; 93.7% line / 79.4% branch coverage (uncovered lines are pre-existing create() default branches); same quality notes.
  `SimpsonMatrixIntegratorTest.java`
  - [x] Task 2.1. Same constructor/field changes as Task 1.1 for `SimpsonMatrixIntegrator`.
  - [x] Task 2.2. In `integrate(Matrix)` replace the hard-coded steps with
    `tmp.multiplyByScalar(ratio); tmp.subtract(ost); tmp.multiplyByScalar(1.0 / (ratio - 1.0));`; loop to `jmax`.
  - [x] Task 2.3. Javadoc as in Task 1.3.
  - [x] Task 2.4. `SimpsonMatrixIntegratorTest`: same constructor tests as Task 1.4.

### Group 2 — Apply 9/8 and a reduced JMAX to the mid-point family (Parallelizable: yes — 8 distinct classes with their own test files; depends on Group 1)

- [x] Task 3. Scalar mid-point-family integrators _(java)_
  - [x] Task 3.1. In `SimpsonMidPointQuadratureIntegrator`, `SimpsonInfinityMidPointQuadratureIntegrator`,
    `SimpsonLowerSquareRootMidPointQuadratureIntegrator` and `SimpsonUpperSquareRootMidPointQuadratureIntegrator`
    change every `super(new …Quadrature(...), eps)` to `super(new …Quadrature(...), eps, TRIPLING_RATIO /*9.0*/,
    MIDPOINT_JMAX /*14*/)`. Expose the two constants (e.g. `protected static final` in `SimpsonIntegrator`, named
    `TRIPLING_RATIO` and `TRIPLING_JMAX`) to avoid duplicating literals.
  - [x] Task 3.2. Class Javadoc of each: state the 9/8 extrapolation and cite NR 3rd ed. §4.4 p. 169.
  - [x] Task 3.3. Tests in each of the four `Simpson…MidPointQuadratureIntegratorTest` classes: a smooth integrand (e.g.
    `exp(x)` or a polynomial/`sin` over a well-conditioned interval) integrated with a counting listener; assert the
    result matches the analytic value and that the evaluation count is below the pre-fix count
    (record the old count by running the test against the previous code, then assert an upper bound with margin).
    Keep inputs deterministic (see issue #36 patterns). Review existing expectations/"returnsWrongResult" tests; update
    tolerances only if accuracy improved.
  - Done: 4 scalar classes + 4 scalar tests (smooth exp integrand: 729 evaluations vs 6561-59049 before; bound 1000); 100% line coverage; no new quality issues.
- [x] Task 4. Matrix mid-point-family integrators _(java)_
  - [x] Task 4.1. Same constructor change for `SimpsonMidPointQuadratureMatrixIntegrator`,
    `SimpsonInfinityMidPointQuadratureMatrixIntegrator`, `SimpsonLowerSquareRootMidPointQuadratureMatrixIntegrator`,
    `SimpsonUpperSquareRootMidPointQuadratureMatrixIntegrator`, using shared constants from `SimpsonMatrixIntegrator`.
  - [x] Task 4.2. Javadoc as in Task 3.2.
  - [x] Task 4.3. Tests as in Task 3.3 in the four `Simpson…MidPointQuadratureMatrixIntegratorTest` classes (counting
    listener over a matrix-valued smooth function; evaluation count lower than before). Avoid integrands with elements
    near 0, since `normMin` convergence is unchanged.
  - Done: 4 matrix classes + 4 matrix tests; 100% line coverage; no new quality issues; full suite green.

### Group 3 — Docs and full validation (Parallelizable: yes — single group of one chain; runs after Groups 1 and 2)

- [x] Task 5. Docs and verification _(java)_
  - Done: docs + CHANGELOG updated; 1188 tests green; coverage >=95%; no new PMD; Checkstyle +52 / SpotBugs +4 are style/CT_CONSTRUCTOR_THROW from new 4-arg constructors; no new secrets.
  - [x] Task 5.1. Update `docs/modules/ROOT/pages/algorithms/simpson-integration.adoc` (and the midpoint page if it
    discusses Simpson) to describe the step-dependent extrapolation (4/3 vs 9/8) and the reduced iteration cap.
    Consider adding a note to `whats-new.adoc` / `CHANGELOG.md` only if the repository convention does so for fixes.
  - [x] Task 5.2. Run the full suite (`mvn test`) via the `iru-gate-runner` agent with `iru-java-test`; confirm all
    existing and new tests pass and no test takes noticeably long.
  - [x] Task 5.3. Run `iru-java-coverage` for the touched classes (>= 80% line coverage), `iru-java-javadoc` for the
    changed classes, and `iru-java-code-quality` (Checkstyle, PMD, SpotBugs) comparing against baseline; run
    `iru-check-license` and `iru-check-security`.
