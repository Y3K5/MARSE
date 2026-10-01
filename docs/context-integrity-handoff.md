# Context and numerical-reference handoff

This contribution extracts generic verification and provenance patterns from a
separate educational viewer. It adds no biological model, fitted parameter set,
pathogen interaction, immune dynamics or clinical claim to MARSE. It changes no
production trajectories.

## Start here

1. Read `docs/context-exchange.md` and the synthetic context fixture.
2. Run `python -m pytest tests/test_context.py` for contract rejection checks.
3. With Node.js installed, run `node examples/numerical_reference/check.cjs`.
4. Read its JSON report: its engine is `pocket-numerical-reference`, not MARSE.

The Node reference is optional example code, not a Python runtime dependency and not
a provider wired into MARSE. It uses anonymous synthetic fields and dimensionless
length, time, resources and biomass. The protocol fixes the domain, endpoint times,
refinement grids and tolerances before executing checks. It is an included test
protocol, not a claimed independently registered biological study.

The checks cover resource-supported growth, positive integration under mortality and
boundary stress, conservation, reproducibility, stable identity-based seeding,
permutation invariance and fixed-domain grid/time refinement. Growth credit follows
resource allocation; reactions use a shared snapshot. The common allocation factor
is conservative and can underutilize a resource. Negative/nonfinite state throws
instead of silently clipping. The implementation is a rectangular synthetic reference.

Exact cosine diffusion is an analytical comparison. Time and coupled-grid checks
compare total biomass with a finer numerical reference. The coupled-grid endpoint is
weakly sensitive in this smooth setup and does not establish convergence of individual
fields, sharp interfaces or other geometries. A passing ledger alone cannot establish
that resource uptake supports growth.

Four numerical weaknesses were reproduced in the separate educational pocket kernels:
growth credit before resource allocation, accepted negative-state settings, sequential
reaction dependence and list-order seeding. Those observations are not evidence of
the same bugs in MARSE's Python engine. No MARSE defect is asserted or patched here.
The local viewer's prior trajectories, data, case studies and private research records
are excluded from this contribution.

## Review priorities

- Reuse MARSE's `SeedRegistry`, configuration validation and run manifests when
  mapping this example to production interfaces; preserve provider ownership.
- Compare the reference checks with existing Python numerical/invariance tests.
  Add only genuinely missing synthetic cases; establish failures independently before
  changing production algorithms.
- Review field units, scope-transfer rules and typed digests before a provider or
  viewer consumes the context format. The current validator checks metadata only.
- Extend a viewer to expose source provenance, scenario assumptions, unresolved values
  and independent status axes. The standalone viewer and Mac bundles are not included.
- Keep model calibration and held-out biological evaluation separate. This contribution
  supplies neither a calibration dataset nor an independently evaluated biological model.

Molecular-coordinate import, mediator mapping, structural-tool integration and native
access assessment remain unimplemented. This contract makes those omissions visible
without replacing them with plausible-looking coordinates or neutral scores.

## Verification

At handoff, on the branch's first base (main at 248c44e), the full suite passed 364
tests with Python 3.14, including the 40 context contract tests. The optional Node
reference passed its synthetic stress, invariance, conservation and refinement checks.
The pinned pre-commit checks passed, including formatting, secrets, repository
structure, file privacy and commit identity. The selected contribution files were also
screened for private workspace references; none were found.

After main was merged in (stages 0 to S1) and the validator moved to
`marse.evidence.context`, the full suite passed on Python 3.12 and 3.14: 860 tests,
with 11 expected failures for the known defects of `docs/validation.md`. A slow test
runs the Node reference, about 35 s with Node 22, and checks that the hashes it reports
are the ones the synthetic example binds
(`python -m pytest -m slow tests/test_examples.py`). Every pre-commit hook and the
repository guard passed on all files. These checks do not establish biological
calibration or validation.
