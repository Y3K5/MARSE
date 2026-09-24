# Genotype-to-capability mapping

MARSE supports an explicit, evidence-bounded mapping from declared loci and
alleles to capability parameter modifiers. A `GenotypeRule` can adjust a
capability's maximum rate or substrate half-saturation when its locus has the
specified allele.

This is intentionally not a DNA interpreter and does not infer traits from
arbitrary sequences. Every rule must name an existing capability and may carry
an evidence source. Missing alleles leave capabilities unchanged. Future work
can connect these records to strain manifests and uncertainty ensembles once
the mapping is constrained by measurements.
