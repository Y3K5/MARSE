"""Structure adapters score a stand-in payload. They do not call a model."""

from __future__ import annotations

import pytest

from marse.structure import alphafold_adapter, esmfold_adapter
from marse.structure.record import StructureError


def _payload(sequence: str, *, pae: bool, low: float = 80.0) -> dict[str, object]:
    length = len(sequence)
    body: dict[str, object] = {"plddt": [low] * length}
    if pae:
        body["pae"] = [[1.5] * length for _ in range(length)]
    return body


def test_esmfold_and_alphafold_use_the_injected_caller():
    sequence = "ACDE"

    def esm(seq: str) -> dict[str, object]:
        assert seq == sequence
        return _payload(seq, pae=False, low=95.0)

    def alpha(seq: str) -> dict[str, object]:
        return _payload(seq, pae=True, low=60.0)

    esm_record = esmfold_adapter(esm).predict(sequence)
    alpha_record = alphafold_adapter(alpha).predict(sequence)
    assert esm_record.adapter == "esmfold"
    assert esm_record.qc.pae_present is False
    assert esm_record.qc.mean_pae is None
    assert esm_record.qc.mean_plddt == 95.0
    assert esm_record.qc.high_fraction == 1.0
    assert alpha_record.adapter == "alphafold"
    assert alpha_record.qc.pae_present is True
    assert alpha_record.qc.mean_pae == 1.5
    assert alpha_record.qc.confident_fraction == 0.0
    assert "growth rate" in esm_record.to_dict()["limits"]
    assert "sequence" not in esm_record.to_dict()


def test_an_unconfigured_adapter_does_not_call_out():
    with pytest.raises(StructureError, match="not configured"):
        esmfold_adapter(None).predict("ACD")
    with pytest.raises(StructureError, match="not configured"):
        alphafold_adapter().predict("ACD")


def test_scores_outside_range_are_refused():
    def bad_plddt(_sequence: str) -> dict[str, object]:
        return {"plddt": [120.0]}

    def bad_pae(_sequence: str) -> dict[str, object]:
        return {"plddt": [80.0, 80.0], "pae": [[0.0, 1.0]]}

    with pytest.raises(StructureError, match="plddt"):
        esmfold_adapter(bad_plddt).predict("A")
    with pytest.raises(StructureError, match="square"):
        alphafold_adapter(bad_pae).predict("AC")
    with pytest.raises(StructureError, match="amino acids"):
        esmfold_adapter(lambda _seq: {"plddt": [1.0]}).predict("BX")
