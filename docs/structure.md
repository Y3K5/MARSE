# Structure scores

Two adapters, ESMFold and AlphaFold, turn a sequence into a checked record.
Each adapter calls a function the caller passes in. This package does not
download a model, does not import a folding library, and does not contact a
server. A caller that is not configured is refused.

pLDDT is required and must be one finite value in 0..100 per residue. PAE,
when the caller sends it, must be a square matrix of finite, non-negative
numbers in angstroms. The record reports the mean and minimum pLDDT, the
fraction at least 90, the fraction at least 70, and, when PAE is present,
its mean and maximum. AlphaFold publishes those pLDDT cuts as its own
confidence bins (Jumper and others, Nature, 2021,
[doi:10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2)).
Reporting a bin is not a claim that the fold is correct, and a missing PAE
is recorded rather than filled in.

The sequence must be the twenty standard residues. Ambiguity codes are
refused. The record written for a person stores the length, not the
sequence, so a feed file does not have to carry the chain.

## Path to a MARSE parameter

A structure record is not a parameter.

1. Keep the record: adapter name, length, and the scores above.
2. A later parameter file may name that record as its source. This package
   does not write that file.
3. Nothing here sets `maximum_per_h`, `yield_mol_per_mol`, a biomass
   formula, or an initial concentration.

A high pLDDT does not choose a substrate, an organism, or a rate.
