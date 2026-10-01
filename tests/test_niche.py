"""Tests for interpretable niche responses and scans."""

from pathlib import Path

import pytest

from marse.microbes.niche import (
    Capability,
    NicheError,
    NicheScan,
    NicheSpecies,
    evaluate_capability,
    load_niche_scan,
    run_niche_scan,
)


def capability() -> Capability:
    return Capability(
        "aerobic-respiration",
        1.0,
        "carbon",
        0.2,
        temperature_c=(5.0, 37.0, 47.0),
        ph=(5.0, 7.0, 9.0),
        oxygen_half_saturation=0.1,
        evidence_source="illustrative benchmark",
        evidence_confidence="C",
    )


def test_optimum_conditions_have_high_rate_and_no_limiting_factor():
    result = evaluate_capability(
        capability(),
        {"carbon": 100.0, "oxygen": 100.0, "temperature_c": 37.0, "ph": 7.0},
    )
    assert result.rate_per_h > 0.99
    assert result.limiting_factor in {"carbon", "oxygen", "temperature", "ph"}


def test_niche_limits_rate_and_reports_the_strongest_limit():
    result = evaluate_capability(
        capability(),
        {"carbon": 10.0, "oxygen": 0.0, "temperature_c": 37.0, "ph": 7.0},
    )
    assert result.rate_per_h == 0.0
    assert result.limiting_factor == "oxygen"


def test_missing_condition_is_explicit():
    with pytest.raises(NicheError, match="missing condition"):
        evaluate_capability(capability(), {"carbon": 1.0})


def test_niche_scan_covers_cartesian_product_and_exports(tmp_path: Path):
    scan = NicheScan(
        "scan",
        {
            "temperature_c": (20.0, 37.0),
            "ph": (6.0, 7.0),
            "oxygen": (1.0,),
            "carbon": (1.0,),
        },
        (NicheSpecies("aerobe", (capability(),)),),
    )
    result = run_niche_scan(scan)
    assert len(result.rows) == 4
    path = result.write_json(tmp_path / "scan.json")
    assert '"growth_rate_per_h"' in path.read_text()
    result.write_csv(tmp_path / "scan.csv")
    assert (tmp_path / "scan.csv").is_file()


def test_niche_scan_rejects_missing_required_axis():
    scan = NicheScan(
        "missing-axis",
        {"carbon": (1.0,)},
        (
            NicheSpecies(
                "aerobe",
                (Capability("respiration", 1.0, "carbon", 0.5, oxygen_half_saturation=0.1),),
            ),
        ),
    )

    with pytest.raises(NicheError, match="oxygen"):
        run_niche_scan(scan)


def test_json_niche_scan_loader(tmp_path: Path):
    path = tmp_path / "scan.json"
    path.write_text(
        '{"experiment_id":"scan","axes":{"temperature_c":[20,37],"ph":[7],'
        '"carbon":[1]},"species":[{"name":"aerobe","capabilities":[{"id":"resp",'
        '"maximum_rate_per_h":1,"substrate":"carbon","half_saturation":0.2,'
        '"temperature_c":[5,37,47],"ph":[5,7,9]}]}]}'
    )
    scan = load_niche_scan(path)
    assert scan.species[0].capabilities[0].id == "resp"
