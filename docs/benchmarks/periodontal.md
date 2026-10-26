# Periodontal benchmark

Two case files, scored by `marse bench`. A pass means the declared check
held. It is not a claim that either file is a pocket, a measurement, or a
calibrated community.

```bash
marse bench list examples/benchmarks/periodontal
marse bench run examples/benchmarks/periodontal -o runs/periodontal-bench
marse bench report runs/periodontal-bench
```

`marse bench run` on that directory exits 1. The version-2 case passes its
balance checks. The version-1 case does not pass, because its engine is
unverified. To run only the column:

```bash
marse bench run examples/benchmarks/periodontal -o runs/periodontal-bench \
  --case periodontal-v2
```

## What is scored

`examples/benchmarks/periodontal/periodontal-v2.yml` runs
`examples/networks/periodontal_pathogen_v2.json` with `marse run`. The
engine is `reactive_transport`. It is not marked unverified. The case
checks the same residual `tests/test_periodontal_v2.py` already checks:
carbon, nitrogen and electrons in `outputs.balance`, each
`largest_relative_residual` at most `1e-12`. It also checks that the run
wrote a manifest, a summary and a report, that the summary does not mark
the engine unverified, and that replaying the manifest reproduces
`outputs`. Replay shows the run can be repeated. It does not show that the
model is right.

`examples/benchmarks/periodontal/periodontal-v1.yml` runs
`examples/experiments/periodontal_pathogen_biofilm.json` with
`marse ecosystem`. The engine is `legacy_ecosystem_v1`. The case records
that name, that `summary.json` says the engine is unverified, and that the
manifest and report exist. `marse bench run` prints the same warning as
`marse ecosystem`. The case result is not a pass. No balance residual and
no biomass array is scored.

## What is not a comparison

The two files are not the same experiment, so version-1 biomass is not
scored against the column.

| | Version 2 | Version 1 |
|---|---|---|
| File | `periodontal_pathogen_v2.json` | `periodontal_pathogen_biofilm.json` |
| Engine | `reactive_transport` | `legacy_ecosystem_v1`, unverified |
| Geometry | one column, 32 voxels of 5 µm | 60 by 40 cells of 10 µm |
| Duration | 2 h | 24 h |
| Substrates | alanine, ammonium, carbon dioxide, acetate | oxygen, heme, peptides, carbon |
| Oxygen | not a component; the three organisms are obligate anaerobes | growth uses an oxygen term |

The names *T. denticola*, *T. forsythia* and *P. gingivalis* appear in both
files. That is the overlap. It is not a shared state to subtract.

## References already in this repository

No file here is a measured reference series for this community, so the
cases do not use `array_reference`.

- `tests/test_periodontal_v2.py` is the balance check the version-2 case
  repeats. The yields and rates in the column file are placeholders.
- [validation.md](../validation.md#the-two-dimensional-ecosystem-engine-is-not-yet-verified)
  lists the open defects of the version-1 engine. Until those are fixed,
  that engine supports no quantitative or comparative conclusion. The
  column does not use that engine. It does not close those defects, and it
  is not an IWA benchmark score.
- `examples/experiments/periodontal_pathogen_provenance.json` has status
  `exploratory-placeholder`. Its calibration status is `not-calibrated`.
  Its paths are fields of the version-1 file. It is not a set of measured
  values for the column.
- The benchmark table in `validation.md` does not contain this study.
  V9 and V10 there are not marked passed by these cases.

## Retiring version 1

This page does not remove `periodontal_pathogen_biofilm.json` or
`marse ecosystem`. Both still run. The version-1 numbers stay off the list
of results to cite, which is already what
[periodontal-pathogen-study.md](../periodontal-pathogen-study.md) says.

Retiring the version-1 file means a later change that stops offering its
numbers as something a reader might cite. That change is not this one. It
waits until the repository records one of the following, and says which:

1. The defects of the two-dimensional ecosystem engine in `validation.md`
   are fixed, and the tests in `tests/test_ecosystem_known_defects.py` are
   ordinary regressions rather than expected failures. A fixed engine still
   does not make the 60 by 40 map a measurement of the column, or a result
   to set beside it.
2. Or the file is withdrawn from the examples a reader is invited to run,
   and the study page says it is gone.

Until one of those is true, the file stays, the benchmark keeps recording
only the engine name and the unverified flag, and nothing in this case
deletes the engine.
