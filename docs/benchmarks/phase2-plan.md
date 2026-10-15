# Phase 2: benchmark harness

This page is the technical plan for a benchmark harness in this repository.
A case that passes has met the comparison it declares. That is not a claim
that a model is a calibrated description of a habitat, and it is not a new
biological result.

The comparison of one array with an expected array already lives in
`marse.validation.benchmarks` (`rmse`, `max_absolute`, `max_relative`). This
harness runs a case and writes a scored report. It calls that registry when
a case supplies arrays. It does not replace it.

## Scope

In this work:

- A YAML file for each case, committed with the repository.
- An adapter interface. The first adapter runs MARSE in-process through the
  engine registry. It does not start a shell and it does not start a container.
- Room for a later adapter to run another tool. A container image is an
  optional field on an adapter. No image is added with the harness. A case
  is valid with the field empty.
- Metrics read from the run directory the engines already write.
- A scored report next to that directory.
- Three commands: `marse bench list`, `marse bench run`, `marse bench report`.

Left where it already is:

- The version-1 ecosystem engine and `periodontal_pathogen_biofilm.json`.
  They stay runnable. `validation.md` records that engine as unverified.
- The IWA biofilm problems. `validation.md` lists V9 as reachable and V10 as
  proposed. This harness does not mark either as passed.
- Example configuration checksums. A benchmark run writes new files beside a
  run. It does not change `config_sha256`.

## Case file

One YAML document per case. Paths in the file are relative to the repository
root. The harness refuses a case before any run when a required field is
missing.

```yaml
id: example-case
title: What this one comparison checks
adapter: marse
experiment: examples/experiments/two_species_batch.json
command: run
metrics:
  - id: config-recorded
    kind: manifest_present
  - id: carbon-balance
    kind: balance_residual
    quantity: C
    tolerance: 1.0e-12
notes: >
  The sentence that says what this case does not claim.
```

`command` is `run`, `ecosystem`, or `none`. `none` does not start an
engine. An `ecosystem` case uses the unverified
engine, so the report cannot call it a quantitative pass. `notes` is required
and is copied into the report.

A metric that compares arrays adds `metric` (`rmse`, `max_absolute`, or
`max_relative`), `tolerance`, `unit`, `expected`, and `simulated`. Those two
arrays are scored with `compare_benchmark`.

## Adapters

An adapter has a name, a `run(case, output_dir)` method, and a `container`
value. `container` is either empty or an image name. The MARSE adapter leaves
it empty and calls the registry:

- `command: run` loads the experiment with `registry.load` and calls `run`.
- `command: ecosystem` loads `legacy_ecosystem_v1` directly, matching
  `marse ecosystem`. The registry would otherwise infer `legacy_batch` for
  that file, which is the behaviour `tests/test_periodontal_v2.py` already
  records.

The method returns the paths of `manifest.json`, `summary.json`, and
`report.html`. Those three files are what `publish_run` already writes.
`marse compare` reads the manifest. The harness does not parse the text of
`marse compare`; it uses the same manifest fields that command compares
(`config_sha256` and `outputs`).

A second adapter is not shipped until a case needs one. The interface is what
keeps that from rewriting the scorer.

## Metrics

Every metric has an `id` and a `kind`. A numeric comparison also has a
`tolerance`. The harness PR implements only kinds that read records the
engines already write.

| Kind | Passes when |
|---|---|
| `manifest_present` | `manifest.json` exists and `config_sha256` is 64 hex characters |
| `summary_present` | `summary.json` has `layout_version` 1 and the same `experiment_id` and `config_sha256` as the manifest |
| `report_present` | `report.html` exists and is the page `publish_run` wrote |
| `unverified_disclosed` | If the engine is unverified, `summary.json` says so. The case then cannot pass as a quantitative result |
| `replay_matches` | Replaying the manifest reproduces `outputs`. Optional, because some runs are slow. Same comparison as `marse compare` after a replay |
| `balance_residual` | The manifest's `outputs.balance` entry for `quantity` has `largest_relative_residual` at most `tolerance`. Version-2 engines write this field. A case that asks for it on a run without the field fails |
| `array_reference` | `compare_benchmark` score is within `tolerance` |
| `not_scored` | Never. The statement is stored and the case fails. Use it when the repository has no reference to score |

A case passes only when every metric passes. An unverified engine forces the
case result to not-passed, and the report says why. Metric rows are still
recorded.

## Reports

`marse bench run` writes two files beside the case's run directory:

- `bench.json` records the case id, adapter name, container value, engine,
  unverified flag, `config_sha256`, and each metric's score, tolerance, and
  pass bit, plus the paths of `manifest.json`, `summary.json`, and
  `report.html`.
- `bench.html` is a short page. It points at the run's `report.html`. It
  repeats the limit that page already states: a replay shows the run can be
  repeated, and it does not show that the model is right.

`marse bench report` reads `bench.json` files and prints a table. It does not
re-run. `marse bench list` reads case files and prints id, adapter, and
experiment. It does not run.

The command exits 0 when every selected case passed, and 1 when a case ran
and failed a metric. A case file that does not load exits 2, the same
convention as an invalid experiment.

## Milestones

These are pull requests, stacked on `cursor/phase1-finish-3075`. None of them
is a date.

1. **This page.** No harness code.
2. **Harness core.** Schema and loader, the adapter interface, the in-process
   MARSE adapter, the metrics above, `bench.json` and `bench.html`, and
   `marse bench list`, `run`, and `report`. Tests use a tiny batch fixture
   already covered by the suite. No new biological case.
3. **Periodontal column.** One case for
   `examples/networks/periodontal_pathogen_v2.json`. The balance residual is
   the check `tests/test_periodontal_v2.py` already makes. The version-1 file
   is recorded only as an unverified engine: its biomass is not scored
   against the column, because the two files are not the same experiment.
   The report repeats the limits in
   [periodontal-pathogen-study.md](../periodontal-pathogen-study.md). The
   version-1 file stays. Retiring it means a later change that stops treating
   its numbers as a result to cite, after the defects in `validation.md` are
   fixed or the file is withdrawn. This milestone does not delete it.
   The cases are `examples/benchmarks/periodontal/periodontal-v2.yml` and
   `periodontal-v1.yml`. What they score is
   [periodontal.md](periodontal.md).

## Acceptance

Harness pull request:

- A case missing `id`, `experiment`, `notes`, or a required `tolerance` is
  refused before a run.
- `marse bench list` prints the fixture cases and writes nothing.
- `marse bench run` on a tiny committed experiment writes `manifest.json`,
  `summary.json`, `report.html`, `bench.json`, and `bench.html`.
- A metric outside its tolerance fails the case and the command exits 1.
- `replay_matches` agrees with the manifest `outputs` that `marse compare`
  reads.
- The MARSE adapter's container field is empty.
- Example `config_sha256` values are unchanged.
- New modules are clean under the mypy baseline.
- `tests/test_benchmarks.py` still passes.

Periodontal pull request, when it is opened:

- The case runs the version-2 file and checks the balance residual.
- The version-1 file contributes its engine name and unverified flag only.
- The written report repeats the public limits of that study.
- `marse ecosystem` and the version-1 file still run.

Parameter values, a publication claim, and whether the version-1 engine is
deleted stay in the example files, `validation.md`, and the periodontal study.
This page does not decide them.
