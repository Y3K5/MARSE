# Scores the repository does not have

Two `marse bench` cases record limits. Neither runs a simulation. The file
they name, `examples/benchmarks/limits/not-run.json`, is not an experiment.
`command: none` does not load an engine. The metric `not_scored` fails on
purpose, so `marse bench run` on this directory exits 1. That exit is the
result. It is not a biological failure and not a pass.

```bash
marse bench list examples/benchmarks/limits
marse bench run examples/benchmarks/limits -o runs/limits
```

## Peptidoglycan

`peptidoglycan.yml` states that this repository has no peptidoglycan
reaction and no measurement of one. Nothing in the case is a cell wall, a
fragment, or a yield. There is no reference array to compare.

## Antigen selection

`antigen-selection.yml` states the limit already written in
[specification.md](../specification.md): vaccine efficacy prediction and de
novo antigen design are deferred. The case does not rank an antigen, does
not read a structure score, and does not write a parameter.

A structure pLDDT, an abundance table, and these two cases are separate.
None of them fills in the score the other is missing.
