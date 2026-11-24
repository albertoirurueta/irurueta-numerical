# Implementation Plan: Fix running covariance bias in MeasurementNoiseCovarianceEstimator (Issue #17)

## 1. Task summary

Source: GitHub issue #17

`MeasurementNoiseCovarianceEstimator.update(double[])` updates the running mean first and then uses the
**post-update** mean for both factors of the outer product folded into the running covariance. The exact incremental
(Welford/West) population covariance update is

```
Σ_{N+1} = (N·Σ_N + (x − μ_N)(x − μ_{N+1})ᵀ) / (N+1)
```

i.e. one factor uses the pre-update mean, the other the post-update mean. Since `x − μ_{N+1} = (x − μ_N)·N/(N+1)`,
the term equals `δδᵀ · N/(N+1)` with `δ = x − μ_N`, which is symmetric. The current code uses
`(x − μ_{N+1})(x − μ_{N+1})ᵀ = δδᵀ · (N/(N+1))²`, under-weighting each sample by `N/(N+1)`.
Check: for `[0, 2, 4]` the current code yields 5/3 while the exact population variance is 8/3.

Approach chosen: compute `sampleNoMean` from the **pre-update** mean `δ`, build `singleCov = δδᵀ`, scale it by
`sampleCount / nextCount` (equivalently, scale the final term), and only then update `sampleAverage`. This avoids
extra fields/arrays, keeps the result exactly symmetric, and keeps the existing population-covariance (divide by N)
semantics. Public API and the sample-average behavior are unchanged.

Additionally the existing `testUpdate` compares the estimator against `Σ xxᵀ / N` (the raw second moment, **not**
mean-removed), tolerated by `LARGE_ABSOLUTE_ERROR`; it only passes because the bias and the uncentered form roughly
cancel for near-zero-mean data. It must be rewritten to compare against the true centered batch covariance.

## 2. Current code state

- `src/main/java/com/irurueta/numerical/signal/processing/MeasurementNoiseCovarianceEstimator.java`, method
  `update(double[] sample)`: updates `sampleAverage` in a loop, then `ArrayUtils.subtract(sample, sampleAverage,
  sampleNoMean)`, copies it into `sampleMatrix`, transposes into `transposedSampleMatrix`, `singleCov =
  sampleMatrix * transposedSampleMatrix`, then `measurementNoiseCov = (measurementNoiseCov*N + singleCov)/(N+1)`.
- `src/test/java/com/irurueta/numerical/signal/processing/MeasurementNoiseCovarianceEstimatorTest.java`:
  `testConstructor`, `testUpdate` (file-based, compares to uncentered second moment, loose tolerance).
- `KalmanFilterTest` also references the estimator; it must still pass (re-run it).
- Docs: `docs/modules/ROOT/pages/algorithms/measurement-noise-covariance-estimation.adoc` documents the current
  formula as "exact", contains a NOTE about the discrepancy, and a mermaid flowchart showing the old order of
  steps. `docs/modules/ROOT/pages/algorithms.adoc` and `reference.adoc` only link to the page.
- `CHANGELOG.md` has an `## [Unreleased]` section with a `### Changed` list.

## 3. Implementation steps

### Group 1 — Fix and tests (Parallelizable: no — Task 2 tests assert the behavior Task 1 introduces and both edit the same class/test pair)

- [x] Task 1. Fix the covariance update in `MeasurementNoiseCovarianceEstimator.update` _(java)_ — implemented in MeasurementNoiseCovarianceEstimator.java; validated (no new quality issues).
  - [x] Task 1.1. Reorder so `sampleNoMean = sample − sampleAverage` is computed with the pre-update mean, before
    the `sampleAverage` loop. Then update `sampleAverage`.
  - [x] Task 1.2. Scale the single-sample covariance by `sampleCount / (double) nextCount` before folding it in
    (e.g. `singleCov.multiplyByScalar((double) sampleCount / nextCount)` after the multiply), keeping the existing
    `(cov*N + singleCov)/(N+1)` rescale. For `sampleCount == 0` the factor is 0, so the first sample contributes 0
    covariance (correct for a single sample).
  - [x] Task 1.3. Update the Javadoc of `update` / class to state that the result is the exact running population
    covariance (matches the batch population covariance of all samples so far); update the comment on
    `sampleNoMean` if needed.
- [x] Task 2. Update unit tests in `MeasurementNoiseCovarianceEstimatorTest` _(java)_ — 5 tests pass; estimator 88.9% line / 100% branch coverage; no new PMD/SpotBugs/Checkstyle issues.
  - [x] Task 2.1. Add `testUpdateKnownSamples`: samples `[0]`, `[2]`, `[4]` with `measureParams = 1`; assert the
    covariance after each update equals the batch population variance (0, 1, 8/3) and the mean (0, 1, 2) with a
    tight tolerance (`1e-12`).
  - [x] Task 2.2. Add a multivariate case (e.g. 2–3 dims, handful of fixed samples) comparing against a
    batch covariance computed in the test as `Σ (xᵢ−μ)(xᵢ−μ)ᵀ / N`, tight tolerance; also assert symmetry.
  - [x] Task 2.3. Fix `testUpdate` to compare against the **centered** batch covariance
    (`Σ (xᵢ−μ)(xᵢ−μ)ᵀ / N`) instead of `Σ xxᵀ / N`, and tighten the tolerance from `LARGE_ABSOLUTE_ERROR` to
    a small relative/absolute error that still holds for the file data.
  - [x] Task 2.4. Delegate running `MeasurementNoiseCovarianceEstimatorTest` and `KalmanFilterTest` to the
    `iru-gate-runner` agent (`iru-java-test` skill) rather than running them in the main conversation.

### Group 2 — Documentation and changelog (Parallelizable: yes)

- [x] Task 3. Update the Antora algorithm page
  `docs/modules/ROOT/pages/algorithms/measurement-noise-covariance-estimation.adoc`
  - [x] Task 3.1. Rewrite the covariance formula to
    `Σ_{N+1} = (N·Σ_N + (x−μ_N)(x−μ_{N+1})ᵀ)/(N+1)`, mention the equivalent `δδᵀ·N/(N+1)` form, and state the
    result equals the batch population covariance at every step.
  - [x] Task 3.2. Remove the NOTE about under-weighting/bias for small N, and the sentence recommending a large
    number of samples because of it (keep the generic TIP about collecting samples in a fixed state).
  - [x] Task 3.3. Update the mermaid flowchart: compute `delta = x − average` first, then `singleCov = delta *
    delta^T`, scale by `N/(N+1)`, update the average, then fold into the covariance.
  - [x] Task 3.4. Run `iru-java-javadoc`-style check not needed; build docs via `iru-build-docs` only if the user
    asks (docs-only text/diagram change). Not run, per caller instruction.
- [x] Task 4. Add a `### Fixed` entry under `## [Unreleased]` in `CHANGELOG.md`: "Fixed
  `MeasurementNoiseCovarianceEstimator` running covariance under-weighting new samples (bias for small sample
  counts), now exactly matching the batch population covariance (#17)."

### Verification (end of plan)

- Run the whole suite and static analysis (`mvn clean compile checkstyle:checkstyle pmd:pmd spotbugs:spotbugs`)
  via `iru-gate-runner`; confirm no new Checkstyle/PMD/SpotBugs findings and coverage of the class stays >= 80%.
