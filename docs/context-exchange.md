# Evidence context exchange, version 1

`marse.evidence.context.validate_context_exchange` checks metadata before a viewer
or an external consumer joins it to a run. It is a standalone contract, not a
scientific provider or a format already consumed by MARSE's simulation loop.

The schema identifier is `marse.context-exchange/1`. Its collections are `evidence`,
`contexts`, `organisms`, `parameters` and `relationships`. The synthetic example is
`examples/context_exchange/synthetic.json`; it contains no biological measurements,
patient data, sequences or molecular coordinates.

## Records and boundaries

Every collection uses stable IDs and five separate status axes: implementation,
execution, gate, evidence and approval. A contract can pass while its biological
evidence remains unresolved. The validator checks consistency; it cannot verify the
truth of a source, approve publication or establish biological validity.

| Record | Required scope |
|---|---|
| Evidence | Source, host, observation site, method, observation, annotation digest, source digest or null, independence group and limits |
| Context | Host/site/disease context; separate biofilm age, phase, physiology and host response; representation/frame; engine/version/run/units/seed policy and typed digests |
| Organism | Identity, host, name, strain or explicit unresolved label, representation; sequence digest or null |
| Parameter | Context/host/site/state, compartment, units, value or explicit null, evidence class, method, uncertainty, calibration role, identifiability, transfer limits and next observation |
| Relationship | Context, participants, observation class, conditions, direction, alternatives, source links, limits and unresolved molecular state |

Source-file, configuration-file, output-file, annotation-text and normalized-sequence
SHA-256 digests have different types. An annotation checksum is not a downloaded
article checksum. A FASTA file checksum cannot satisfy a normalized-sequence identity
join. A declared sequence reference needs an accession. The validator does not fetch,
normalize or analyze sequences and does not authenticate source files itself.

An executed context needs an output digest; the unexecuted template has none. Its
source/protocol digests bind the included generic numerical reference files. They do
not imply that the template's context was simulated. Run/source byte comparisons
remain the caller's responsibility.

Host scope is checked across organisms, relationships and evidence. Measurements
from saliva cannot silently fill pocket-lumen parameters. Site and disease context
remain attached to every parameter. Unknown measurements must remain null. Measured
or inferred values need source links; unknowns are not demonstrated zeros.

Supported relationship classes include co-occurrence, spatial neighbors, functional
perturbation, resource exchange, host-mediated observations and unresolved relations.
More than two participants are allowed; they do not automatically produce pairwise
edges. This revision only accepts uncoupled relationship records, unresolved mediator
and structure links, and explicit false contact/native-access claims. Molecular
transforms are not implemented. Density fields do not establish cell contacts.

## Check a document

```python
import json
from pathlib import Path

from marse.evidence.context import validate_context_exchange

document = json.loads(Path("examples/context_exchange/synthetic.json").read_text())
summary = validate_context_exchange(document)
print(summary)
```

Run `python -m pytest tests/test_context.py`. Tests include wrong IDs, missing units,
cross-host and cross-compartment joins, duplicate participants, false structural
resolution, nonfinite values, incompatible digest types and independent statuses.

Future versions need an explicitly reviewed mapping to existing MARSE manifests and
provider records. The generic exchange must not introduce a second source of engine
identity or bypass the core's privacy-preserving provenance machinery.
