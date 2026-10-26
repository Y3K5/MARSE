# Table import

`marse ingest` reads one CSV or TSV table, checks it, and writes a JSON feed.
The feed format is `marse.ingest/1`. A feed is not a run, and a check is not
a claim that the rows were measured.

```bash
marse ingest environment examples/ingest/synthetic_environment.csv -o feed.json
marse ingest abundance examples/ingest/synthetic_abundance.csv -o feed.json
marse ingest transcriptomics examples/ingest/synthetic_transcripts.csv -o feed.json
marse ingest metabolomics examples/ingest/synthetic_metabolite.csv -o feed.json
```

The files under `examples/ingest/` are synthetic. They are not a microbiome,
a transcriptome, or a metabolite panel from a sample.

## What is checked

| Kind | Required columns | A row fails when |
|---|---|---|
| `environment` | `name`, `value`, `unit` | the name is not `temperature_c` (unit `C`, −20 to 150) or `ph` (unit `pH`, 0 to 14) |
| `abundance` | `sample`, `taxon`, `abundance`, `unit` | the abundance is negative, or relative rows in one sample do not sum to 1 |
| `transcriptomics` | `sample`, `feature`, `value`, `unit` | the unit is not `count` or `tpm`, or the value is negative |
| `metabolomics` | `sample`, `feature`, `value`, `unit` | the unit is not `mmol_per_l`, `mol_per_m3`, or `arbitrary`, or the value is negative |

A missing column, a repeated key, an empty table, or a non-finite value
exits 2 and writes nothing. Counts may be zero.

## What can enter an experiment file

`experiment_fragment` is the only part of the feed that uses experiment
field names.

- An environment table copies `temperature_c` and `ph` into
  `environment`. Those are the two fields a batch experiment already has.
  Nothing else is copied.
- One metabolite, one sample, with unit `mmol_per_l` or `mol_per_m3`,
  copies that number to `substrate.initial_mm`. One mmol per litre is one
  mol per cubic metre, which is that field. Two features, two samples, or
  unit `arbitrary` stay in the table and the fragment is null.
- An abundance table and a transcriptomics table always leave the fragment
  null. A relative abundance is not a biomass, and a count or a TPM is not
  a rate.

`apply_fragment` copies those keys onto a dict. If the experiment already
has a different value for the same key, it refuses. It does not add
organisms, yields, or a domain.

Another importer can register with `marse.ingest.register`. Only the four
kinds above are registered by the package.
