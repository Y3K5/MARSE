# MARSE

MARSE (Microbial Adaptability Resource Simulation Engine) is open-source
software for reproducible simulation of microbial populations and biofilms.
This site is the documentation in this repository.

```bash
marse run examples/experiments/two_species_batch.json -o runs/batch
marse replay runs/batch/manifest.json
marse ls runs
```

Build this site locally with `python -m pip install -e ".[docs]"` and then
`mkdocs build -f docs/mkdocs.yml --strict`. Continuous integration runs that
build. The site is not deployed from pull requests.

- [Architecture](architecture.md) — the package as it is.
- [Table import](ingest.md) — checked tables. Abundance is not biomass.
- [Validation](validation.md) — what has been checked, and what has not.
- [Periodontal study](periodontal-pathogen-study.md) — the version-2 column, and the unverified version-1 file. The [benchmark](benchmarks/periodontal.md) scores the column's balance and does not score version 1.
