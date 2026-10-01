"""Frames are recorded on demand and can stream to disk, and recording never changes a result.

The central property is the first test: however frames are recorded (every
step, every few steps, only the ends, or streamed to a store), the final
state is bit-identical. Frames observe a run; they must never steer it.
"""

import json
import tracemalloc
from pathlib import Path

import numpy as np
import pytest

from marse.cli import main
from marse.ecosystem import (
    EcosystemConfig,
    EcosystemError,
    NutrientConfig,
    SeedRegion,
    SpeciesConfig,
    run,
    write_viewer,
)
from marse.ecosystem.framestore import (
    EcosystemFrameSink,
    FrameStore,
    FrameStoreError,
    StoredFrames,
)
from marse.ecosystem.model import recorded_steps
from marse.ecosystem.viewer import _rounded
from marse.microbes.niche import Capability


def config(**overrides) -> EcosystemConfig:
    values = {
        "experiment_id": "frames-test",
        "width": 8,
        "height": 6,
        "cell_size_um": 10.0,
        "duration_h": 1.0,
        "timestep_h": 0.01,
        "seed": 5,
        "mutation_interval_h": 0.25,
        "nutrients": (
            NutrientConfig("oxygen", 1.0, 10.0, 1.0, ("top",)),
            NutrientConfig("carbon", 2.0, 5.0),
        ),
        "species": (
            SpeciesConfig(
                "aerobe",
                0.0,
                0.8,
                (0.2, 0.4),
                (1.0, 1.0),
                mutation_probability=0.1,
                seed_regions=(SeedRegion(0.3, 0.6, 0.3, 0.2),),
                spreading_per_h=20.0,
                capabilities=(Capability("carbon-use", 0.8, "carbon", 0.4),),
            ),
            SpeciesConfig("competitor", 0.05, 0.4, (0.3, 0.2), (1.0, 1.0)),
        ),
    }
    return EcosystemConfig(**(values | overrides))


def digest(result) -> str:
    return result.manifest.outputs["final_state_sha256"]


def stream(conf: EcosystemConfig, directory: Path, frame_every: int | None):
    capacity = len(recorded_steps(conf.steps, frame_every))
    sink = EcosystemFrameSink.create(directory, conf, capacity=capacity)
    try:
        result = run(conf, frame_every=frame_every, sink=sink)
    finally:
        sink.close()
    return result, StoredFrames(FrameStore.open(directory))


# --- recording never changes the result ---------------------------------------


def test_recording_frames_never_changes_the_result(tmp_path: Path):
    conf = config()
    reference = digest(run(conf))
    for frame_every in (3, 7, None):
        assert digest(run(conf, frame_every=frame_every)) == reference
    streamed, _ = stream(conf, tmp_path / "frames", 4)
    assert digest(streamed) == reference


@pytest.mark.parametrize(
    ("steps", "frame_every", "expected"),
    [
        (10, 1, tuple(range(11))),
        (10, 3, (0, 3, 6, 9, 10)),
        (10, 5, (0, 5, 10)),
        (10, None, (0, 10)),
        (10, 50, (0, 10)),
    ],
)
def test_recorded_steps_always_include_the_first_and_last(steps, frame_every, expected):
    assert recorded_steps(steps, frame_every) == expected


@pytest.mark.parametrize("bad", [0, -2, 1.5, True])
def test_a_meaningless_frame_interval_is_refused(bad):
    with pytest.raises(EcosystemError, match="frame_every"):
        run(config(), frame_every=bad)


def test_in_memory_frames_follow_the_requested_steps():
    conf = config()
    result = run(conf, frame_every=30)
    times = [frame.time_h for frame in result.frames]
    assert len(result.frames) == len(recorded_steps(conf.steps, 30)) == 5
    assert times[0] == 0.0
    assert times[-1] == pytest.approx(conf.duration_h)


# --- the frame store ------------------------------------------------------------


def test_a_streamed_run_keeps_only_its_first_and_latest_frame(tmp_path: Path):
    conf = config()
    result, stored = stream(conf, tmp_path / "frames", 10)
    assert len(stored) == len(recorded_steps(conf.steps, 10))
    assert len(result.frames) == 2
    assert result.frames[0].time_h == 0.0
    assert result.frames[-1].time_h == stored.store.times_h[-1]


def test_the_store_holds_what_the_run_recorded(tmp_path: Path):
    conf = config()
    in_memory = run(conf, frame_every=10).frames
    _, stored = stream(conf, tmp_path / "frames", 10)
    assert stored.store.steps == recorded_steps(conf.steps, 10)
    assert len(stored) == len(in_memory)
    for kept, frame in zip(stored, in_memory, strict=True):
        assert kept.time_h == frame.time_h
        # Stored in single precision: relative agreement to float32's resolution,
        # and values below float32's normal range flush towards zero.
        for field in ("biomass", "nutrients"):
            np.testing.assert_allclose(
                getattr(kept, field), getattr(frame, field), rtol=1e-6, atol=1e-30
            )
        for a, b in zip(kept.niche_rates, frame.niche_rates, strict=True):
            np.testing.assert_allclose(a, b, rtol=1e-6, atol=1e-30)
        np.testing.assert_array_equal(kept.mutations, frame.mutations)
        np.testing.assert_array_equal(kept.phenotype_indices, frame.phenotype_indices)
        for a, b in zip(kept.niche_limiting_factors, frame.niche_limiting_factors, strict=True):
            np.testing.assert_array_equal(a, b)


def test_the_store_index_describes_the_run_without_identifying_the_machine(tmp_path: Path):
    conf = config()
    result, _ = stream(conf, tmp_path / "frames", 25)
    text = (tmp_path / "frames" / "index.json").read_text("utf-8")
    index = json.loads(text)
    assert index["complete"] is True
    assert index["metadata"]["species"] == ["aerobe", "competitor"]
    assert index["metadata"]["config_sha256"] == result.manifest.config_sha256
    assert str(tmp_path) not in text
    assert str(Path.home()) not in text


def test_a_store_refuses_misuse(tmp_path: Path):
    store = FrameStore.create(
        tmp_path / "s", capacity=1, fields={"a": ((2,), "<f8"), "b": ((1,), "<u1")}
    )
    with pytest.raises(FrameStoreError, match="missing"):
        store.append(time_h=0.0, step=0, arrays={"a": [1.0, 2.0]})
    with pytest.raises(FrameStoreError, match="shape"):
        store.append(time_h=0.0, step=0, arrays={"a": [1.0], "b": [1]})
    store.append(time_h=0.0, step=0, arrays={"a": [1.0, 2.0], "b": [1]})
    with pytest.raises(FrameStoreError, match="full"):
        store.append(time_h=1.0, step=1, arrays={"a": [1.0, 2.0], "b": [1]})
    store.close()
    with pytest.raises(FrameStoreError, match="already holds"):
        FrameStore.create(tmp_path / "s", capacity=1, fields={"a": ((2,), "<f8")})
    with pytest.raises(FrameStoreError, match="reading"):
        FrameStore.open(tmp_path / "s").append(time_h=0.0, step=0, arrays={})
    with pytest.raises(FrameStoreError, match="no index"):
        FrameStore.open(tmp_path / "nowhere")


def test_overwriting_a_store_removes_only_its_own_files(tmp_path: Path):
    root = tmp_path / "s"
    FrameStore.create(root, capacity=1, fields={"old": ((1,), "<f8")}).close()
    (root / "notes.txt").write_text("keep me", encoding="utf-8")
    FrameStore.create(root, capacity=1, fields={"new": ((1,), "<f8")}, overwrite=True).close()
    assert sorted(p.name for p in root.iterdir()) == ["index.json", "new.npy", "notes.txt"]


def test_an_interrupted_store_says_so_and_stays_readable(tmp_path: Path):
    store = FrameStore.create(tmp_path / "s", capacity=3, fields={"a": ((1,), "<f8")})
    store.append(time_h=0.0, step=0, arrays={"a": [1.5]})
    store.close()
    reopened = FrameStore.open(tmp_path / "s")
    index = json.loads((tmp_path / "s" / "index.json").read_text("utf-8"))
    assert index["complete"] is False
    assert len(reopened) == 1
    assert reopened.read(0)["a"].tolist() == [1.5]


def test_streaming_keeps_memory_flat_as_runs_get_longer(tmp_path: Path):
    """The point of the store. The unstreamed case shows the test can tell the difference.

    A warm-up run comes first, so one-time allocations are not mistaken for
    growth. Without it the first run measured carries them, and the result
    would depend on test order. What still scales with a streamed run is the
    store's index (a time and a step per frame), which on this grid is small
    against the working arrays.
    """
    grid = {"width": 60, "height": 60, "timestep_h": 0.02}

    def peak(duration_h: float, streamed: bool, name: str) -> int:
        conf = config(duration_h=duration_h, **grid)
        tracemalloc.start()
        if streamed:
            stream(conf, tmp_path / name, 1)
        else:
            run(conf, frame_every=1)
        _, top = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        return top

    peak(0.2, streamed=True, name="warm-up")
    assert peak(4.0, streamed=True, name="long") / peak(1.0, streamed=True, name="short") < 1.2
    assert peak(4.0, streamed=False, name="") / peak(1.0, streamed=False, name="") > 2.0


# --- the viewer -----------------------------------------------------------------


def test_the_viewer_is_capped_self_contained_and_keeps_the_ends(tmp_path: Path):
    result = run(config())
    viewer = write_viewer(result, tmp_path / "viewer.html", max_frames=10)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["viewer.html"]
    html = viewer.read_text("utf-8")
    payload = json.loads(html.split("const data = ", 1)[1].split(";\nconst canvas", 1)[0])
    assert len(payload["frames"]) == 10
    assert payload["recorded_frames"] == len(result.frames)
    assert payload["frames"][0]["time_h"] == 0.0
    assert payload["frames"][-1]["time_h"] == result.frames[-1].time_h
    assert "10 of 101 recorded frames" in html
    # Totals come from the unrounded values.
    exact = result.frames[-1].to_dict(("aerobe", "competitor"), ("oxygen", "carbon"), 1.0)
    assert payload["frames"][-1]["statistics"] == exact["statistics"]


def test_the_viewer_reads_from_a_store(tmp_path: Path):
    result, stored = stream(config(), tmp_path / "frames", 5)
    viewer = write_viewer(result, tmp_path / "viewer.html", frames=stored)
    assert f"of {len(stored)} recorded frames" in viewer.read_text("utf-8")


def test_viewer_rounding_keeps_four_significant_digits_of_the_largest_value():
    values = np.array([123.456789, 0.00012345, 0.0])
    np.testing.assert_array_equal(_rounded(values), [123.5, 0.0, 0.0])
    np.testing.assert_array_equal(_rounded(np.array([0.0012345678])), [0.001235])
    np.testing.assert_array_equal(_rounded(np.zeros(3)), np.zeros(3))


# --- the command line -----------------------------------------------------------


def _experiment(tmp_path: Path) -> Path:
    path = tmp_path / "experiment.json"
    path.write_text(json.dumps(config().to_dict()), encoding="utf-8")
    return path


def test_cli_streams_frames_and_writes_no_frames_json(tmp_path: Path):
    output = tmp_path / "run"
    assert main(["ecosystem", str(_experiment(tmp_path)), "-o", str(output)]) == 0
    assert sorted(p.name for p in output.iterdir()) == ["frames", "manifest.json", "viewer.html"]
    index = json.loads((output / "frames" / "index.json").read_text("utf-8"))
    assert index["complete"] is True
    assert index["frame_count"] == config().steps + 1  # 100 steps: under the default cap


def test_cli_frame_interval_sets_how_many_frames_are_stored(tmp_path: Path, capsys):
    output = tmp_path / "run"
    args = ["ecosystem", str(_experiment(tmp_path)), "-o", str(output)]
    args += ["--frame-interval-h", "0.25"]
    assert main(args) == 0
    index = json.loads((output / "frames" / "index.json").read_text("utf-8"))
    assert index["steps"] == [0, 25, 50, 75, 100]
    assert "5 stored, one every 25 step(s)" in capsys.readouterr().out
    assert main([*args[:-1], "0"]) == 2


def test_cli_replay_with_output_writes_the_same_layout(tmp_path: Path):
    output = tmp_path / "run"
    assert main(["ecosystem", str(_experiment(tmp_path)), "-o", str(output)]) == 0
    again = tmp_path / "again"
    assert main(["replay", str(output / "manifest.json"), "-o", str(again)]) == 0
    assert sorted(p.name for p in again.iterdir()) == ["frames", "manifest.json", "viewer.html"]
