"""The documented examples must keep running and keep telling the truth.

Examples that quietly break are worse than no examples, so each one is run as a
subprocess and its self-checks are asserted here.
"""

import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


def run_example(name: str) -> str:
    result = subprocess.run(
        [sys.executable, str(EXAMPLES / name)],
        capture_output=True,
        text=True,
        check=True,
        timeout=300,
    )
    assert result.stderr == ""
    return result.stdout


@pytest.mark.parametrize("name", ["batch_growth.py", "oxygen_penetration.py"])
def test_example_runs_cleanly(name):
    assert run_example(name).strip()


def test_batch_growth_reproduces_the_analytical_solution():
    output = run_example("batch_growth.py")
    # The simulated and predicted final biomass are printed to six decimals.
    simulated, predicted = (
        float(line.split(":")[1].split()[0])
        for line in output.splitlines()
        if "final biomass" in line
    )
    assert simulated == pytest.approx(predicted, rel=1e-6)

    drift_line = next(line for line in output.splitlines() if "drift" in line)
    assert float(drift_line.split(":")[1].split()[0]) < 1e-10

    # Every reported relative error against the closed form must be tiny.
    errors = [
        float(parts[3])
        for line in output.splitlines()
        if len(parts := line.split()) == 4 and "e-" in parts[3]
    ]
    assert len(errors) >= 5
    assert max(errors) < 1e-6


def test_oxygen_penetration_is_physically_consistent():
    output = run_example("oxygen_penetration.py")
    depth = float(
        next(line for line in output.splitlines() if "penetration depth" in line)
        .split("=")[-1]
        .replace("um", "")
    )
    # Tens of micrometres, as reported for dense colony biofilms.
    assert 10.0 < depth < 200.0
    assert "anoxic" in output
    # Deeper than the penetration depth there must be no oxygen left.
    assert "0.0000      0.0%" in output
