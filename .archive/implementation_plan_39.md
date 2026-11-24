# Implementation Plan: Checkstyle CLI run ignores checkstyle.xml; fix config and existing violations (Issue #39)

## 1. Task summary

Source: GitHub issue #39
Base branch: develop

`mvn checkstyle:checkstyle` (the command documented in `CLAUDE.md`) does not use the project's `checkstyle.xml`,
because `<configLocation>checkstyle.xml</configLocation>` is only set on the `maven-checkstyle-plugin` entry under
`<reporting>` in `pom.xml` (lines ~308-322). A CLI run therefore falls back to the plugin's default Sun checks (80-char
lines etc.) and reports 3942 mostly meaningless violations. The work: make CLI and site use the same rules, record the
real baseline, and reduce it.

Real baseline measured while planning (`mvn checkstyle:checkstyle -Dcheckstyle.config.location=checkstyle.xml`, `src/main`):
**1010 violations** — HiddenField 443, Indentation 203, MagicNumber 101, JavadocStyle 87, RedundantModifier 67,
ImportOrder 28, ThrowsCount 14, JavadocVariable 12, DeclarationOrder 11, ParameterAssignment 9, ExplicitInitialization 9,
JavadocPackage 7, TrailingComment 4, VisibilityModifier 4, ParameterNumber 3, NeedBraces 2, NewlineAtEndOfFile 2,
WhitespaceAround 1, ConstantName 1, MultipleStringLiterals 1, JavadocMethod 1. The final numbers are re-measured in Task 2.

Choices made on the user's behalf (challenge during review):
- The plugin is configured under `<build><plugins>` (not only `<pluginManagement>`) so a bare CLI invocation picks up
  the config. The `<reporting>` entry keeps working; it is left as is.
- No `failsOnError`/`failOnViolation` gate and no `validate`-phase execution binding: Checkstyle stays report-only, as today.
- `HiddenField` (443, the largest category) is accepted on purpose for constructors and setters via
  `ignoreConstructorParameter`/`ignoreSetter` in `checkstyle.xml`; this is the example the issue itself gives.
- Mechanical, low-risk categories are fixed in this PR (small ones only); large or judgment-heavy categories
  (Indentation, MagicNumber, JavadocStyle, RedundantModifier, ThrowsCount, the remainder of HiddenField) are documented and left
  to follow-up PRs, since the issue says to fix "in focused PRs". Production behavior must not change.
- The baseline is documented in `CLAUDE.md` (static analysis note) rather than a new docs page.

Open clarifications: none

## Acceptance criteria coverage

| # | Criterion | Covered by tasks |
|---|---|---|
| 1 | "`mvn checkstyle:checkstyle` uses `checkstyle.xml`, so `target/checkstyle-checker.xml` shows `max=120`." | Task 1 |
| 2 | "The real violation count against the project rules is documented." | Task 2 |
| 3 | "Violations are fixed or deliberately suppressed, and the count goes down." | Task 3, Task 4 |

## 2. Current code state

- `pom.xml`: `<build><plugins>` (from line 133) has no checkstyle plugin; `<reporting><plugins>` has
  `maven-checkstyle-plugin` 3.6.0 with `configLocation=checkstyle.xml` and a `checkstyle` report set.
- `checkstyle.xml` (repo root, 102 lines): `Checker` with `JavadocPackage`, `NewlineAtEndOfFile`, `FileTabCharacter`,
  `LineLength` (max 120) and a `TreeWalker` with ~50 checks, including a bare `<module name="HiddenField"/>`.
  No suppressions file exists.
- `CLAUDE.md` documents `mvn clean compile checkstyle:checkstyle pmd:pmd spotbugs:spotbugs` and the config rules.
- Stack: Java 21, Maven (no wrapper), JUnit Jupiter; `src/main` only is analyzed (tests are not in the Checkstyle scope).
- Reports land in `target/checkstyle-result.xml`, `target/checkstyle-checker.xml`.

## 3. Implementation steps

### Group 1 — Fix the configuration (Parallelizable: yes)

- [x] Task 1. Make CLI Checkstyle use `checkstyle.xml` _(java)_
  - Touched `pom.xml` only (config change, no Java code, so no new tests/coverage applicable; Checkstyle now uses checkstyle.xml with max=120).
  - [x] Task 1.1. In `pom.xml`, add `maven-checkstyle-plugin` 3.6.0 under `<build><plugins>` with
    `<configuration><configLocation>checkstyle.xml</configLocation></configuration>` (same version as the reporting entry; no
    `<executions>`, so no build phase is affected).
  - [x] Task 1.2. Verify: run `mvn clean compile checkstyle:checkstyle` and confirm `target/checkstyle-checker.xml` contains
    `<property name="max" value="120"/>` and no Sun default `max=80`. Run `mvn -q -o test-compile` to confirm the POM is valid.
    Delegate the Maven runs to the `iru-gate-runner` agent so logs stay out of the main context.

### Group 2 — Record the real baseline (Parallelizable: yes; depends on Group 1)

- [x] Task 2. Document the real violation count _(untagged, documentation)_
  - [x] Task 2.1. Re-run Checkstyle with the fixed config and tally violations per check from `target/checkstyle-result.xml`.
  - [x] Task 2.2. In `CLAUDE.md`, add a short "Checkstyle baseline" note under the build/test notes: the config now applies to CLI
    runs, the per-check table (before this PR's fixes), and the date/version measured.

### Group 3 — Deliberate suppressions in config (Parallelizable: yes; depends on Group 2)

- [x] Task 3. Accept `HiddenField` in constructors and setters _(java)_
  - [x] Task 3.1. In `checkstyle.xml`, change `<module name="HiddenField"/>` to include
    `ignoreConstructorParameter=true`, `ignoreSetter=true` and `setterCanReturnItsClass=true`, with a comment saying this is
    intentional (parameters named after fields are the project convention).
  - [x] Task 3.2. Re-run Checkstyle, record the new HiddenField count and total in the `CLAUDE.md` note from Task 2.2.
    - Done: `checkstyle.xml`, `CLAUDE.md`. HiddenField 443 -> 113, total 1010 -> 680. Config-only, no code/tests changed.

### Group 4 — Fix small mechanical categories (Parallelizable: no — tasks edit overlapping `src/main` files and share one Checkstyle run)

- [x] Task 4. Fix low-risk violations in `src/main` without behavior changes _(java)_
  - [x] Task 4.1. `JavadocPackage` (7): add missing `package-info.java` files (with Javadoc describing each package).
  - [x] Task 4.2. `JavadocVariable` (12) and `JavadocMethod` (1): add the missing Javadoc.
  - [x] Task 4.3. `NewlineAtEndOfFile` (2), `NeedBraces` (2), `WhitespaceAround` (1), `ConstantName` (1), `TrailingComment` (4),
    `ExplicitInitialization` (9): fix in place.
  - [x] Task 4.4. `ImportOrder` (28) and `DeclarationOrder` (11): reorder only; no logic changes.
  - [x] Task 4.5. Run the full test suite via `iru-gate-runner` (`mvn test`) and confirm everything still passes; re-run
    Checkstyle and update the final count in the `CLAUDE.md` note (Task 2.2), noting which categories remain for follow-up PRs.
  - [x] Task 4.6. Run `/iru-java-javadoc`-style audit on touched files and `iru-check-license` on any new files (license header
    must match existing files).
    - Done: 1186 tests pass, overall line coverage 93.22%; Checkstyle total 680 -> 602 (11 mechanical categories at 0); 7 new
      `package-info.java` files carry the Apache header (iru-check-license: all compliant). Remaining categories documented in `CLAUDE.md`.
