"""Write and read the scored report next to a run directory.

``bench.json`` is the record. ``bench.html`` points at the run's own
``report.html`` and repeats its limit: a replay shows that the run can be
repeated, and it does not show that the model is right. Paths stored here are
file names, not absolute paths.
"""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any

from marse.bench.cases import Case
from marse.bench.metrics import Evaluation


def document(
    case: Case,
    evaluation: Evaluation,
    *,
    engine: str,
    unverified: bool,
    container: str,
) -> dict[str, Any]:
    """The scored record. ``container`` is null when the adapter did not use one."""
    return {
        "adapter": case.adapter,
        "case_id": case.id,
        "config_sha256": evaluation.config_sha256,
        "container": container or None,
        "engine": engine,
        "manifest": "manifest.json",
        "metrics": [metric.to_dict() for metric in evaluation.metrics],
        "notes": case.notes,
        "passed": evaluation.passed,
        "reason": evaluation.reason,
        "report": "report.html",
        "summary": "summary.json",
        "title": case.title,
        "unverified": unverified,
    }


def write_report(output_dir: Path, payload: dict[str, Any]) -> None:
    """Write ``bench.json`` and ``bench.html`` into the case directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    (output_dir / "bench.json").write_text(text, encoding="utf-8")
    (output_dir / "bench.html").write_text(_html(payload), encoding="utf-8")


def read_reports(directory: Path) -> list[dict[str, Any]]:
    """Every ``bench.json`` under ``directory``, ordered by case id."""
    found: list[dict[str, Any]] = []
    if not directory.is_dir():
        return found
    paths = [directory / "bench.json"] if (directory / "bench.json").is_file() else []
    paths.extend(sorted(path for path in directory.glob("*/bench.json") if path not in paths))
    for path in paths:
        loaded: Any = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(loaded, dict):
            found.append(loaded)
    return sorted(found, key=lambda row: str(row.get("case_id", "")))


def _html(payload: dict[str, Any]) -> str:
    title = html.escape(str(payload["title"]))
    status = "pass" if payload["passed"] else "fail"
    reason = html.escape(str(payload["reason"]))
    notes = html.escape(str(payload["notes"]))
    rows = []
    for metric in payload["metrics"]:
        flag = "pass" if metric["passed"] else "fail"
        rows.append(
            "<tr>"
            f"<td>{html.escape(str(metric['id']))}</td>"
            f"<td>{html.escape(str(metric['kind']))}</td>"
            f"<td>{flag}</td>"
            f"<td>{html.escape(str(metric['detail']))}</td>"
            "</tr>"
        )
    body = "\n".join(rows)
    limit = "A replay can be repeated. It does not show that the model is right."
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
</head>
<body>
<h1>{title}</h1>
<p>{status}. {reason}</p>
<p>{notes}</p>
<p>{limit}</p>
<p><a href="report.html">Run report</a></p>
<table>
<thead><tr><th>Metric</th><th>Kind</th><th>Result</th><th>Detail</th></tr></thead>
<tbody>
{body}
</tbody>
</table>
</body>
</html>
"""
