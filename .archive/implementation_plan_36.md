# Implementation Plan: Fix flaky SimpsonMidPointQuadratureMatrixIntegratorTest (Issue #36)

## 1. Task summary

Source: GitHub issue #36
Base branch: develop

`SimpsonMidPointQuadratureMatrixIntegratorTest.integrate_whenExponential1_returnsExpectedResult` failed CI run
45422386783 with an `IntegrationException` thrown from `SimpsonMatrixIntegrator.integrate` (the "too many steps" path,
`JMAX = 20`, about 10 s). The test draws `a`, `b` and `lambda` from an unseeded `UniformRandomizer`, so some draws are
degenerate (`b ≈ a`, `lambda ≈ 0`) and the default relative convergence tolerance `EPS = 1e-10` can't be reached in
double precision. The fix is **test-only**: make the inputs deterministic and well conditioned, and ask the integrator
for a tolerance consistent with the test's own `1e-6` accuracy check. The scalar twin
(`SimpsonMidPointQuadratureIntegratorTest`) has the same pattern and is hardened the same way.

Choices made on the user's behalf (no user input needed):
- Seeded `java.util.Random` (constant seed) passed to `UniformRandomizer`, plus a small fixed loop of trials per test, so
  the test is reproducible yet still covers varied inputs.
- Inputs: interval width at least `MIN_INTERVAL`, and `|lambda|` drawn from `[MIN_ABS_LAMBDA, MAX_LAMBDA]` with a random sign.
- Integrator built with the existing `eps` constructor using `EPS_EXPONENTIAL = 1e-8`.
- No production code changes. The Richardson-extrapolation coefficient in `SimpsonMatrixIntegrator`/`SimpsonIntegrator`
  (trapezoidal-style `(4*st - ost)/3` applied to a step-tripling midpoint rule) is out of scope and, if worth pursuing,
  belongs in a separate issue.

Open clarifications: none

## Acceptance criteria coverage

| # | Criterion | Covered by tasks |
|---|---|---|
| 1 | "`integrate_whenExponential1_returnsExpectedResult` uses deterministic (seeded) inputs and can't draw degenerate `lambda ≈ 0` or `b ≈ a` values." | Task 1 |
| 2 | "The test passes reliably: e.g. 1000 local runs or a repeated-trial loop with no `IntegrationException`." | Task 1, Task 3 |
| 3 | "The test runtime is well under 1 s." | Task 1, Task 3 |
| 4 | "The scalar `SimpsonMidPointQuadratureIntegratorTest` counterpart is reviewed and hardened the same way if needed." | Task 2 |
| 5 | "No production code changes." | Task 3.3 |

## Current code state

- `src/test/java/com/irurueta/numerical/integration/SimpsonMidPointQuadratureMatrixIntegratorTest.java`
  - Constants: `MIN_VALUE = -10.0`, `MAX_VALUE = 10.0`, `MIN_LAMBDA = -1.0`, `MAX_LAMBDA = 1.0`,
    `ABSOLUTE_ERROR_EXPONENTIAL = 1e-6`.
  - `integrate_whenExponential1_returnsExpectedResult` (about line 55): `new UniformRandomizer()`;
    `a = nextDouble(MIN_VALUE, MAX_VALUE)`, `b = nextDouble(a, MAX_VALUE)`, `lambda = nextDouble(MIN_LAMBDA, MAX_LAMBDA)`;
    `new SimpsonMidPointQuadratureMatrixIntegrator(a, b, listener)` (default `EPS = 1e-10`); compares with
    `1/lambda * (exp(lambda*b) - exp(lambda*a))` via `Matrix.equals(..., ABSOLUTE_ERROR_EXPONENTIAL)`.
- `src/test/java/com/irurueta/numerical/integration/SimpsonMidPointQuadratureIntegratorTest.java`
  - Same constants; `integrate_whenExponential_returnsExpectedResult` has the identical random pattern, built with
    `new SimpsonMidPointQuadratureIntegrator(a, b, point -> Math.exp(lambda * point))` and `assertEquals(..., 1e-6)`.
- Production (reference only, **not modified**): `SimpsonMatrixIntegrator` (`EPS = 1e-10`, `JMIN = 5`, `JMAX = 20`,
  throws `IntegrationException` when it doesn't converge), `MidPointMatrixQuadrature.next` (triples evaluations per step).
  The `eps` constructor already exists: `SimpsonMidPointQuadratureMatrixIntegrator(a, b, listener, eps)`; the scalar class is
  expected to offer the equivalent (verify in Task 2.1).
- Other tests in the repo use `new UniformRandomizer(new Random())`, so `UniformRandomizer(Random)` is available.
- Project: Java 21, Maven single module, JUnit Jupiter 5, Checkstyle (120-char lines, Javadoc rules), license headers on every file.

## Implementation steps

### Group 1 — Harden the exponential tests (Parallelizable: no — Task 2 reuses the helper approach and constants chosen in Task 1)

- [x] Task 1. Make `integrate_whenExponential1_returnsExpectedResult` deterministic and well conditioned _(java)_ — `SimpsonMidPointQuadratureMatrixIntegratorTest.java` (SEED=7_736_412_987L, TIMES=20, MIN_INTERVAL=0.5, MIN_ABS_LAMBDA=0.1, EPS_EXPONENTIAL=1e-8; MIN_LAMBDA removed). Deviation: with EPS 1e-8 (relative) the fixed 1e-6 absolute assert failed for results ~2e3, so the tolerance is now `ABSOLUTE_ERROR_EXPONENTIAL * max(1, |expected|)`. Tests pass (1170/1170 full suite), `SimpsonMidPointQuadratureMatrixIntegrator` 100% line coverage, no new Checkstyle/PMD/SpotBugs issues.
  - [x] Task 1.1. In `SimpsonMidPointQuadratureMatrixIntegratorTest.java` add constants: `SEED` (fixed long),
    `TIMES` (e.g. 20), `MIN_INTERVAL` (e.g. 0.5), `MIN_ABS_LAMBDA` (e.g. 0.1), `EPS_EXPONENTIAL = 1e-8`.
    Keep `ABSOLUTE_ERROR_EXPONENTIAL = 1e-6`. Remove `MIN_LAMBDA` only if nothing else uses it.
  - [x] Task 1.2. Rewrite the test body: `new UniformRandomizer(new Random(SEED))`, loop `TIMES` trials; per trial
    `a = nextDouble(MIN_VALUE, MAX_VALUE - MIN_INTERVAL)`, `b = nextDouble(a + MIN_INTERVAL, MAX_VALUE)`,
    `lambda = nextDouble(MIN_ABS_LAMBDA, MAX_LAMBDA) * (randomizer.nextBoolean() ? 1.0 : -1.0)`; build the integrator with
    `new SimpsonMidPointQuadratureMatrixIntegrator(a, b, listener, EPS_EXPONENTIAL)`; assert as today.
    Keep names and the `@Test` signature unchanged.
  - [x] Task 1.3. Check the test imports `java.util.Random` and that lines stay within 120 chars.

- [x] Task 2. Apply the same hardening to the scalar twin _(java)_ — `SimpsonMidPointQuadratureIntegratorTest.java`, same constants, eps constructor used, same relative-scaled tolerance as Task 1. Tests pass, `SimpsonMidPointQuadratureIntegrator` 100% line coverage, no new Checkstyle/PMD/SpotBugs issues.
  - [x] Task 2.1. In `SimpsonMidPointQuadratureIntegratorTest.java` confirm `SimpsonMidPointQuadratureIntegrator` has an
    `eps` constructor (check `SimpsonMidPointQuadratureIntegrator` / `SimpsonIntegrator`); if it doesn't, harden inputs
    only (seeded, well-conditioned) and note it in the progress comment. — confirmed: `(a, b, listener, eps)` constructor exists.
  - [x] Task 2.2. Rewrite `integrate_whenExponential_returnsExpectedResult` with the same seeded loop, input ranges and
    `EPS_EXPONENTIAL`, mirroring Task 1. — done.

### Group 2 — Validation (Parallelizable: yes — single task)

- [x] Task 3. Validate the fix _(java)_ — validated: both test classes pass (6/6 each, every test <= 22 ms); 1000-seed stress (20,000 integrations per test) with 0 `IntegrationException` and 0 tolerance failures (worst error/tolerance 0.00125), constants unchanged; only change was removing the blank line before `import java.util.Random;` in both test files; license headers compliant; full suite 1170 tests green; target classes 100% line coverage; no new Checkstyle issues vs HEAD (project checkstyle.xml run directly, since the Maven quality goals skip test sources); no src/main/java change.
  - [x] Task 3.1. Delegate to `iru-gate-runner` (`Skill({skill: "iru-java-test", args:
    "SimpsonMidPointQuadratureMatrixIntegratorTest,SimpsonMidPointQuadratureIntegratorTest"})`) and confirm both pass and each
    runs in well under 1 s. — both classes 6/6 pass; suites ~0.05-0.1 s, exponential tests 9 ms / 22 ms.
  - [x] Task 3.2. One-off stress check, not committed: temporarily run the new input ranges with many different seeds
    (about 1000 trials per test, e.g. via a throwaway JUnit/`jshell` snippet in the scratchpad) and confirm no
    `IntegrationException` and no tolerance failure. If `EPS_EXPONENTIAL = 1e-8` or the ranges don't hold, adjust the
    constants (not the asserted `1e-6` accuracy) and re-run. — 1000 seeds x 20 trials per test: 0 exceptions, 0 failures, worst error/tolerance 0.00125; constants unchanged.
  - [x] Task 3.3. Confirm `git diff --stat` touches only the two test files (plus `implementation_plan.md` housekeeping),
    with no change under `src/main`. — only the two test files changed (rest is known build/scan noise); no src/main/java diff.
  - [x] Task 3.4. Delegate to `iru-gate-runner` (`iru-java-code-quality`) over the two test files and fix any Checkstyle issue. — fixed new ImportOrder blank line before `java.util.Random` in both files; remaining hits pre-existing.

## Quality decision

The project-wide static-analysis baseline was unchanged by this work: 4241 issues already existed (Checkstyle 3942,
PMD 35, SpotBugs 264). The user chose to file follow-up issues: #39 (Checkstyle; the CLI run ignores `checkstyle.xml`),
#40 (PMD) and #41 (SpotBugs). The deviation in Task 1 (assert tolerance scaled by `max(1, |expected|)`) was reviewed and
accepted by the user. No docs update is needed: the change is test-only.
