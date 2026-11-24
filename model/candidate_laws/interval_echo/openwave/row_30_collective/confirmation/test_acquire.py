"""Local release-integrity and scientific falsifier checks on exposed selectors."""

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from model.candidate_laws.interval_echo.openwave.row_30_collective.confirmation import (
    acquire as trial,
)
from model.candidate_laws.interval_echo.openwave.row_30_collective.confirmation import (
    spatial,
)


@pytest.fixture(scope="module", autouse=True)
def declaration_fixture(tmp_path_factory):
    # Production reads the frozen candidate.json named in the verified launcher.
    # The detector unit tests exercise its checksum with a private synthetic declaration.
    declaration = tmp_path_factory.mktemp("binding") / "candidate.json"
    declaration.write_text('{"scope":"synthetic detector unit test"}')
    with patch.object(trial, "CANDIDATE_FILE", declaration):
        yield


def synthetic_binding():
    revision = "c" * 40
    return dict(
        mode="pilot",
        anchor=None,
        source_commit=revision,
        candidate_pin=revision,
        candidate_sha256=trial.legacy.file_hash(trial.CANDIDATE_FILE),
        protocol_sha256=trial.legacy.file_hash(trial.HERE / "TASK.md"),
    )


@pytest.fixture(scope="module")
def baseline(tmp_path_factory):
    dest = tmp_path_factory.mktemp("ow30") / "release"
    assert trial.acquire(synthetic_binding(), trial.PILOT_SEEDS, dest) == "passed"
    return dest


def test_completed_release_replays_and_matches_exposed_pilot(baseline):
    assert trial.replay(synthetic_binding(), trial.PILOT_SEEDS, baseline) == "passed"
    summary = trial.pilot_summary(baseline)
    assert [board["raw_arrays"] for board in summary["boards"]] == [38, 54, 54]
    assert all(board["spatial"]["launch_label_rejected"] for board in summary["boards"])
    golden = trial.read_json(Path(__file__).with_name("pilot.json"))
    # A live replay verifies every array against its own retained inventory.
    # The golden inventory records the acquisition host's floating sampler.
    manifest = trial.read_json(baseline / "manifest.json")
    for live, reference, board in zip(
        summary["boards"], golden["boards"], manifest["boards"], strict=True
    ):
        inventory = json.dumps(
            board["counts"]["files"], sort_keys=True, separators=(",", ":")
        ).encode()
        assert live["inventory_sha256"] == sha256(inventory).hexdigest()
        assert len(reference["inventory_sha256"]) == 64
        live.pop("inventory_sha256")
        reference.pop("inventory_sha256")
        for reading in (live, reference):
            spatial_reading = reading["spatial"]
            assert spatial_reading["candidate_q"] <= spatial_reading["quadratic_cutoff"]
            assert spatial_reading["launch_label_q"] > spatial_reading["quadratic_cutoff"]
            # Residual sums describe the retained host's counts and floating fit.
            # Each retained replay and each scientific cutoff is checked above.
            for key in list(spatial_reading):
                if key.endswith("_q"):
                    spatial_reading.pop(key)
    trial.experiment.verify(summary, golden)


def test_missing_or_changed_count_and_no_substitution(baseline, tmp_path):
    binding = synthetic_binding()
    original = baseline / "records/0/clock_cos.npy"
    backup = tmp_path / "clock.npy"
    original.rename(backup)
    try:
        assert (
            trial.inspect_release(binding, trial.PILOT_SEEDS, baseline) == "incomplete"
        )
        with patch.object(
            trial.legacy, "run_bank", side_effect=AssertionError("unexpected retry")
        ):
            assert trial.acquire(binding, trial.PILOT_SEEDS, baseline) == "incomplete"
            assert trial.replay(binding, trial.PILOT_SEEDS, baseline) == "incomplete"
            assert (
                trial.acquire(binding, [927016, 927016, 927017], baseline) == "invalid"
            )
        arr = np.load(backup, allow_pickle=False)
        arr[0] = int(arr[0]) - 1
        np.save(original, arr, allow_pickle=False)
        assert trial.inspect_release(binding, trial.PILOT_SEEDS, baseline) == "invalid"
    finally:
        backup.replace(original)


@pytest.mark.parametrize(
    "data",
    [
        np.array([1, 3], dtype="<i8"),
        np.array([1, 3], dtype=">u8"),
        np.array([1, 2, 3], dtype="<u8"),
        np.array([1, 8], dtype="<u8"),
    ],
)
def test_malformed_integer_word(data, tmp_path):
    np.save(tmp_path / "reading.npy", data, allow_pickle=False)
    inventory = trial.audit_counts(tmp_path, {"reading": ((2,), 7)})
    assert set(inventory["invalid"]) == {"reading"}


def test_fully_recorded_propagating_altered_law_gets_rejected(tmp_path):
    dest = tmp_path / "altered"
    spectrum = trial.experiment.dual
    spectroscopy = trial.experiment.spectroscopy
    flight_started = False

    def start(*args, **kwargs):
        nonlocal flight_started
        flight_started = True
        return spectroscopy(*args, **kwargs)

    def altered(length, coupling):
        return spectrum(
            length,
            7.5 if flight_started and coupling == trial.experiment.K else coupling,
        )

    with (
        patch.object(trial.experiment, "spectroscopy", start),
        patch.object(trial.experiment, "dual", altered),
    ):
        assert trial.acquire(synthetic_binding(), trial.PILOT_SEEDS, dest) == "rejected"
    rows = trial.read_json(dest / "manifest.json")["boards"]
    assert [len(row["counts"]["files"]) for row in rows] == [38, 54, 54]
    assert all(
        row["native_gates"]["transport"] and not row["native_gates"]["continuum"]
        for row in rows
    )
    assert (
        trial.inspect_release(synthetic_binding(), trial.PILOT_SEEDS, dest)
        == "rejected"
    )


def test_zero_counts_are_incomplete(tmp_path):
    dest = tmp_path / "outage"
    with patch.object(
        trial.experiment,
        "spectroscopy",
        side_effect=ArithmeticError("interrupted acquisition"),
    ):
        assert (
            trial.acquire(synthetic_binding(), trial.PILOT_SEEDS, dest) == "incomplete"
        )
    rows = trial.read_json(dest / "manifest.json")["boards"]
    assert all(
        row["acquisition_error"]["error"] == "ArithmeticError"
        and not row["counts"]["files"]
        for row in rows
    )
    assert (
        trial.inspect_release(synthetic_binding(), trial.PILOT_SEEDS, dest)
        == "incomplete"
    )


def test_physically_ineligible_material_is_measured_rejection(tmp_path):
    # Feed the enrollment selector only actual single-dimer components of each
    # exposed board. Their addresses and geometry remain available in the receipt.
    grow = trial.experiment.grow

    def single(board, seed):
        return [
            (c, p)
            for c, p in grow(board, seed)
            if not c and len(p) == 1 and p[0][0] // 12 == p[0][1] // 12
        ]

    dest = tmp_path / "short"
    with patch.object(trial.experiment, "grow", single):
        assert trial.acquire(synthetic_binding(), trial.PILOT_SEEDS, dest) == "rejected"
        assert (
            trial.inspect_release(synthetic_binding(), trial.PILOT_SEEDS, dest)
            == "rejected"
        )
    rows = trial.read_json(dest / "manifest.json")["boards"]
    assert all(
        row["enrollment"]["ineligible"] == "enrolled component lies in a single cell"
        for row in rows
    )


def test_launch_label_alone_cannot_pass(baseline):
    raw = baseline / "records"
    physical = trial.read_json(raw / "result.json")["boards"][0]
    ref = spatial.score(physical, raw / "0")
    assert ref["status"] == "pass"
    with patch.object(
        spatial, "score", return_value=dict(ref, launch_label_rejected=False)
    ):
        evidence = trial.board_receipt(
            trial.enrollment(0, trial.PILOT_SEEDS[0]),
            physical,
            raw / "0",
            trial.measured_meter(),
        )
    assert evidence["status"] == "rejected"


def test_only_fixed_validation_domain_or_guarded_public_domain():
    binding = dict(synthetic_binding(), mode="validation")
    assert trial.valid_inputs(binding, trial.VALIDATION_SEEDS)
    assert not trial.valid_inputs(binding, trial.REPEAT_SEEDS)
    assert trial.valid_inputs(dict(binding, mode="repeat"), trial.REPEAT_SEEDS)
    assert not trial.valid_inputs(dict(binding, mode="repeat"), trial.VALIDATION_SEEDS)
    assert not trial.valid_inputs(binding, trial.PILOT_SEEDS)
    public = dict(
        mode="public",
        anchor="a" * 40,
        source_commit="b" * 40,
        candidate_pin="b" * 40,
        candidate_sha256="c" * 64,
        protocol_sha256="d" * 64,
    )
    assert not trial.valid_inputs(public, trial.PILOT_SEEDS)
    assert not trial.valid_inputs(public, [2**63] * 3)
    assert not trial.valid_inputs(
        dict(public, candidate_pin="f" * 40), [2**63 + i for i in range(3)]
    )
    assert len(set(trial.EXIT_CODES.values())) == 4


@pytest.mark.parametrize(
    "mode,selectors",
    [
        ("pilot", trial.PILOT_SEEDS),
        ("repeat", trial.REPEAT_SEEDS),
        ("validation", trial.VALIDATION_SEEDS),
    ],
)
def test_unauthenticated_or_mismatched_local_binding_is_rejected_without_selection(
    tmp_path, mode, selectors
):
    binding = dict(synthetic_binding(), mode=mode)
    dest = tmp_path / "forbidden"
    for change in (
        {"source_commit": None},
        {"candidate_pin": "d" * 40},
        {"protocol_sha256": None},
        {"candidate_sha256": "e" * 64},
        {"protocol_sha256": "f" * 64},
    ):
        assert trial.acquire(dict(binding, **change), selectors, dest) == "invalid"
    assert not dest.exists()


def test_every_wavelength_gate_uses_only_the_even_fit(baseline):
    raw = baseline / "records"
    board = deepcopy(trial.read_json(raw / "result.json")["boards"][0])
    launch = (
        np.pi
        * np.arange(1, len(board["observer"]["modes"]) + 1)
        / (len(board["route"]) + 2)
    )
    for row, k in zip(board["observer"]["modes"], launch, strict=True):
        row["k"] = float(k)
    for key in ("continuum", "lattice", "shape"):
        board["gates"][key] = False
    scored = trial.board_receipt(
        trial.enrollment(0, trial.PILOT_SEEDS[0]),
        board,
        raw / "0",
        trial.measured_meter(),
    )
    assert scored["status"] == "passed"
    assert (
        scored["spatial"]["full_recurrence_diagnostic_q"]
        > scored["spatial"]["quadratic_cutoff"]
    )
    d = scored["dispersion"]
    assert [row["k"] for row in d["modes"]] == scored["spatial"]["even_trained_k"]
    assert [row["full_recurrence_diagnostic_k"] for row in d["modes"]] == list(launch)
    assert all(
        scored["native_gates"][key] and not scored["legacy_recurrence_gates"][key]
        for key in ("continuum", "lattice", "shape")
    )
    wrong_ruler = dict(scored["spatial"], even_trained_k=[1.5] * len(launch))
    wrong = trial.even_hold_dispersion(board, wrong_ruler, trial.measured_meter())
    assert wrong["gates"] == {"lattice": False, "continuum": False, "shape": False}
    assert wrong["low_k_lines"] == 0 and wrong["low_k_max"] is None
