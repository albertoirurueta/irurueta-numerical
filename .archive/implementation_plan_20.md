# Implementation Plan: Add bibliography links (Issue #20)

## 1. Task summary

Issue [#20](https://github.com/albertoirurueta/irurueta-numerical/issues/20) asks to update the bibliography
section of the Antora documentation so that:
- book entries link to the book's official editorial/publisher website, and
- paper entries link to the official site where the paper can be downloaded or its publisher/DOI landing page.

The bibliography lives in a single file:
[docs/modules/ROOT/pages/reference.adoc](docs/modules/ROOT/pages/reference.adoc), under the `[#bibliography]` /
`== Bibliography` section (lines 52-127). This is a **documentation-only change** — no Java source is touched, so
no `Language/framework` tag applies to the tasks below and no `*-code-one-task` skill dispatch is needed; edit
the AsciiDoc file directly.

Of the 12 bibliography entries, 6 already carry a working link (`wikipedia-quadratic`, `wikipedia-cubic`,
`wikipedia-convolution`, `wikipedia-variance`, `welch1995`, and `chum2005`'s link to a third-party reference
implementation). The other 6 need a link added, or an existing bare DOI turned into a clickable one, to the
verified official source found below:

| Anchor | Type | Link to add |
|---|---|---|
| `press2007` | Book | `https://numerical.recipes` (the book's own official site — same link already used for this book in sibling repositories such as `irurueta-sorting`'s [reference.adoc](../irurueta-sorting/docs/modules/ROOT/pages/reference.adoc)) |
| `fischler1981` | Paper | `https://dl.acm.org/doi/10.1145/358669.358692` (ACM Digital Library) |
| `rousseeuw1984` | Paper | `https://www.tandfonline.com/doi/abs/10.1080/01621459.1984.10477105` (Taylor & Francis Online — replaces the current bare DOI text) |
| `torr2000` | Paper | `https://doi.org/10.1006/cviu.1999.0832` (DOI resolver to the Elsevier/ScienceDirect page) |
| `chum2005` | Paper | `https://doi.org/10.1109/CVPR.2005.221` (DOI resolver to the IEEE Xplore page) — added alongside, not replacing, the existing reference-implementation link |
| `irurueta2012` | Paper | `https://ieeexplore.ieee.org/document/6365452/` (IEEE Xplore) |
| `kalman1960` | Paper | `https://asmedigitalcollection.asme.org/fluidsengineering/article/82/1/35/397706/A-New-Approach-to-Linear-Filtering-and-Prediction` (ASME Digital Collection — replaces the current bare DOI text) |

All URLs above were verified via web search against the publisher/digital-library record for each work before
writing this plan, except `press2007`, which reuses the exact link already established for this same book in
this project's sibling repositories (e.g. `irurueta-sorting`) for consistency across the `com.irurueta` family.

Per the issue, this work must be implemented on a branch cut from `release_1.6.0` (already done — branch
`feature/20`) so it is published under the `1.6.0` Antora doc version; after merge it must also be propagated to
`master` and then `develop` — that propagation happens after this PR is reviewed/merged and is outside the scope
of the tasks below.

## 2. Current code state

- [docs/modules/ROOT/pages/reference.adoc](docs/modules/ROOT/pages/reference.adoc) is a single AsciiDoc page with
  three sections: "Generated Reports", "External Dashboards", and "Bibliography" (`[#bibliography]`, line 52
  onward).
- Bibliography entries follow a consistent AsciiDoc definition-list pattern:
  ```asciidoc
  [[anchor-id]]Author(s) (Year)::
  _Title_ or "Title." _Venue_, details. Optional inline link(s) using
  https://example.com[Link text] syntax. Prose describing which class(es)/xref:algorithms/...adoc[] the
  reference supports.
  ```
- Existing entries that already use an inline link (e.g. `welch1995` at line 116-121, `wikipedia-*` entries) are
  the style to match: the URL is wrapped as `https://...[descriptive link text]` and placed inline within the
  entry's descriptive prose (not appended as a separate line).
- Two entries currently expose a DOI as **plain text**, not a link: `rousseeuw1984` ("DOI: 10.1080/...") at line
  83, and `kalman1960` ("DOI: 10.1115/...") at line 113. These need to become clickable links to the DOI's actual
  landing page.
- `docs/antora.yml` confirms the docs component is versioned `1.6.0`, matching the current branch — no version
  bump needed for this change.

## 3. Implementation steps

1. [x] **Add/upgrade the missing bibliography links in `docs/modules/ROOT/pages/reference.adoc`**
   - [x] `press2007` (line 58-63): add an inline link to `https://numerical.recipes` on the book title, matching
     the link and style already used for this book in `irurueta-sorting`'s reference.adoc (`Available at:
     https://numerical.recipes[numerical.recipes].`), e.g. rendering as:
     `_https://numerical.recipes[Numerical Recipes: The Art of Scientific Computing]_, 3rd ed. New York:
     Cambridge University Press. ISBN 978-0-521-88068-8. ...`
   - [x] `fischler1981` (line 76-79): add an inline link to `https://dl.acm.org/doi/10.1145/358669.358692` on the
     paper title, e.g.:
     `"https://dl.acm.org/doi/10.1145/358669.358692[Random Sample Consensus: A Paradigm for Model Fitting with
     Applications to Image Analysis and Automated Cartography]." _Communications of the ACM_, 24(6), 381–395. ...`
   - [x] `rousseeuw1984` (line 81-84): replace the plain-text `DOI: 10.1080/01621459.1984.10477105` with a link to
     `https://www.tandfonline.com/doi/abs/10.1080/01621459.1984.10477105`, e.g.:
     `... 79(388), 871–880. https://www.tandfonline.com/doi/abs/10.1080/01621459.1984.10477105[DOI:
     10.1080/01621459.1984.10477105]. ...`
   - [x] `torr2000` (line 86-90): add an inline link to `https://doi.org/10.1006/cviu.1999.0832` on the paper
     title, following the same pattern as `fischler1981`.
   - [x] `chum2005` (line 92-98): add an inline link to `https://doi.org/10.1109/CVPR.2005.221` on the paper
     title itself (distinct from the existing Devernay reference-C-implementation link, which should stay as-is),
     e.g.:
     `"https://doi.org/10.1109/CVPR.2005.221[Matching with PROSAC – Progressive Sample Consensus]."
     _Proceedings of the IEEE Computer Society Conference on Computer Vision and Pattern Recognition (CVPR)_,
     220–226. ...` (keep the trailing sentence about the Devernay C implementation unchanged).
   - [x] `irurueta2012` (line 100-104): add an inline link to `https://ieeexplore.ieee.org/document/6365452/` on
     the paper title.
   - [x] `kalman1960` (line 111-114): replace the plain-text `DOI: 10.1115/1.3662552` with a link to
     `https://asmedigitalcollection.asme.org/fluidsengineering/article/82/1/35/397706/A-New-Approach-to-Linear-Filtering-and-Prediction`,
     following the same pattern as `rousseeuw1984`.
   - [x] Leave the 6 already-linked entries (`wikipedia-quadratic`, `wikipedia-cubic`, `wikipedia-convolution`,
     `wikipedia-variance`, `welch1995`) untouched.
   - Note: all 7 links added/upgraded in `docs/modules/ROOT/pages/reference.adoc` exactly per the table above; no
     Java source touched; no tests applicable (documentation-only change).

2. [x] **Build the Antora documentation locally to verify the page renders correctly**
   - [x] Confirm `docs/` has its Antora toolchain already set up (per `antora-setup` conventions used elsewhere in
     this repo) and run the existing local build command for this repo's docs site (check `docs/package.json` or
     equivalent scripts, e.g. `npm run build` from `docs/`, whichever is already configured).
   - [x] Open the built `reference.html` bibliography section and confirm every new/updated link renders as a
     clickable link with the expected target and readable link text, and that no AsciiDoc syntax leaked into the
     rendered output (e.g. no literal `https://...[...]` text).
   - Note: ran `npx antora antora-playbook.yml` from `docs/`; build succeeded, producing
     `docs/build/site/irurueta-numerical/1.6.0/reference.html`. Inspected the rendered bibliography section: all
     7 new/updated links (press2007, fischler1981, rousseeuw1984, torr2000, chum2005, irurueta2012, kalman1960)
     render as clickable `<a href="...">` anchors with the correct target and readable text; no leaked AsciiDoc
     syntax. The 6 untouched entries render unchanged.

3. [x] **Spot-check every new/changed URL resolves correctly**
   - [x] For each of the 6 web-searched URLs listed in section 1's table (all except `press2007`, whose
     `https://numerical.recipes` link is reused as-is from `irurueta-sorting` and needs no re-verification),
     confirm it loads the expected publisher/digital-library page for that exact paper (not a 404, paywall-only
     stub with wrong metadata, or unrelated redirect).
     Delegate this to a sub-agent to keep the main context clean:
     ```
     Agent({
       description: "Verify bibliography URLs for issue 20",
       prompt: "Fetch each of these URLs and confirm it resolves (HTTP 200, no redirect to an unrelated page) and
         that the page's title/metadata matches the paper/book named next to it: [list the 7 URL -> work pairs
         from the plan]. Report any URL that fails to resolve or doesn't match.",
       run_in_background: false
     })
     ```
   - [x] If any URL fails verification, find and substitute a corrected official URL for that entry before
     proceeding.
   - Note: sub-agent verified all 6 URLs (fischler1981, rousseeuw1984, torr2000, chum2005, irurueta2012,
     kalman1960). All match the cited paper's title/authors/journal. Two (torr2000's ScienceDirect DOI target and
     kalman1960's ASME page) sit behind a Cloudflare anti-bot gate rather than rendering directly, but the DOI/
     URL itself resolves and the target's identity was independently confirmed via search — no substitution
     needed.

No test suite or coverage step applies here — this is a documentation-only AsciiDoc change with no Java code, no
new classes/methods, and nothing for JUnit/JaCoCo to exercise.
