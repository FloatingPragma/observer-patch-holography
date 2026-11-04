"""False-green mutations for geometry, execution, pressure and interpretation."""
from copy import deepcopy
from pathlib import Path

import pytest

import verify_maxwell_eos as audit

LEAN_MATHLIB_LIB = (Path(__file__).resolve().parents[3] / "Lean" / ".lake"
                    / "packages" / "mathlib" / ".lake" / "build" / "lib")
LEAN_REPLAY_SKIP_REASON = (
    "Lean/.lake cache absent; the replay builds the Whitney Lean targets against "
    "Mathlib and is covered by the Lean CI lane"
)


def mathlib_oleans_built() -> bool:
    """True when the Lean workspace holds a compiled Mathlib."""
    if not LEAN_MATHLIB_LIB.is_dir():
        return False
    return next(LEAN_MATHLIB_LIB.glob("**/Mathlib.olean"), None) is not None


def lean_replay_available() -> bool:
    """True when the resolved Lean toolchain can replay without a cold Mathlib build."""
    return audit.original.lake_path() is not None and mathlib_oleans_built()


requires_lean_replay = pytest.mark.skipif(
    not lean_replay_available(), reason=LEAN_REPLAY_SKIP_REASON)


@pytest.fixture(scope="module")
def receipt():
    return audit.load(Path(__file__).parent / "runs/maxwell_eos_receipt.json")


@requires_lean_replay
def test_full_independent_replay(receipt):
    result = audit.verify(receipt)
    assert result["all_checks_passed"]
    assert result["source_original_fully_replayed"]


def native_claim(p):
    p["contract"]["native_repair_eos"] = True


def wrong_volume(p):
    p["forms"]["volume"] *= 2


def wrong_mass(p):
    p["forms"]["M1"][0][0] *= 1.05


def wrong_pressure(p):
    p["source_free"]["frames"][10]["mean_pressure"] *= 1.1


def wrong_energy(p):
    p["source_free"]["frames"][10]["energy"] *= 1.1


def forged_state(p):
    p["source_free"]["frames"][10]["pi"][0] += 0.01


def invented_conduction(p):
    p["source_free"]["step_checks"][0]["conduction_current_l2"] = 1


def field_only_promoted_to_total(p):
    p["authenticated_observer_fields"]["scope"] = "total equation of state of substrate"


def wrong_deformation(p):
    p["source_free"]["metric_probes"][0]["isotropic"][0]["plus"]["energy"] *= 1.01


def wrong_stress_sign(p):
    p["source_free"]["frames"][0]["stress_tensor"][2][2] *= -1


def zero_ratio(p):
    p["controls"]["zero_field"]["w_mean_stress"] = 1 / 3


def fake_observer_origin(p):
    p["authenticated_observer_fields"]["rows"][0]["decode_event_ids"][0] += 1


def copied_future_observer_state(p):
    p["authenticated_observer_fields"]["rows"][0]["A"][0] += 0.5


def thermalization_claim(p):
    p["classical_gibbs"]["thermalization_demonstrated"] = True


def gauge_mode_in_partition(p):
    p["classical_gibbs"]["physical_modes"] = 42


def thermal_covariance_mutation(p):
    p["classical_gibbs"]["scale_cases"][0]["rows"][0]["covariance_pi"][0][0] += 1


def thermal_partition_mutation(p):
    p["classical_gibbs"]["scale_cases"][0]["rows"][0]["partition_plus"] += 1


@pytest.mark.parametrize("mutate", [native_claim, wrong_volume, wrong_mass, wrong_pressure,
    wrong_energy, forged_state, invented_conduction, field_only_promoted_to_total,
    wrong_deformation, wrong_stress_sign, zero_ratio, fake_observer_origin,
    copied_future_observer_state, thermalization_claim, gauge_mode_in_partition,
    thermal_covariance_mutation, thermal_partition_mutation])
def test_semantic_mutations_fail(receipt, mutate):
    packet = deepcopy(receipt)
    mutate(packet)
    with pytest.raises(ValueError):
        audit.verify(packet, replay_original=False)


def test_lean_replay_guard_decision(monkeypatch, tmp_path):
    """The replay guard skips only when the toolchain or the compiled Mathlib is absent."""
    monkeypatch.setattr(audit.original, "lake_path", lambda: "/usr/bin/lake")
    monkeypatch.setitem(globals(), "LEAN_MATHLIB_LIB", tmp_path / "cold")
    assert lean_replay_available() is False
    built = tmp_path / "built" / "lean"
    built.mkdir(parents=True)
    (built / "Mathlib.olean").write_bytes(b"olean")
    monkeypatch.setitem(globals(), "LEAN_MATHLIB_LIB", built.parent)
    assert lean_replay_available() is True
    monkeypatch.setattr(audit.original, "lake_path", lambda: None)
    assert lean_replay_available() is False
