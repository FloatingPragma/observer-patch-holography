from __future__ import annotations

import json

import pytest

import stratum_verdicts as sv


def test_all_four_verdicts_build_with_expected_statuses() -> None:
    verdicts = sv.build_all()
    assert sorted(verdicts) == [642, 643, 645, 646]
    expected_pin_counts = {642: 1, 643: 3, 645: 2, 646: 3}
    for issue, verdict in verdicts.items():
        assert verdict["status"] == sv.VERDICTS[issue]
        boundary = verdict["comparison_boundary"]
        assert boundary["public_measurement_read"] is False
        assert boundary["comparison_permitted"] is False
        assert len(verdict["parent_pins"]) == expected_pin_counts[issue]
        assert verdict["reopen_condition"]
    assert verdicts[643]["frame_lock_disposition"]["clause"] == "FZ02-R03b"
    blocked = {
        row["combination"]: row["blocked_by"]
        for row in verdicts[646]["tested_combinations"]
    }
    assert blocked["M_W / M_Z"] == "no_pole_promotion"
    assert (
        blocked["any physical-unit member of the vector"]
        == "physical_units_not_evaluable"
    )


def test_dyonic_lattice_certificate_is_exact() -> None:
    lattice = sv.dyonic_line_lattice_certificate()
    assert lattice["maximal_isotropic_count"] == 12
    forms = lattice["untilted_global_forms"]
    assert sorted(forms) == [
        "unquotiented",
        "z2_quotient",
        "z3_quotient",
        "z6_quotient",
    ]
    assert all(row["present_among_maximal_isotropic"] for row in forms.values())
    assert forms["z6_quotient"]["lattice"] == [[0, m] for m in range(6)]
    assert forms["unquotiented"]["lattice"] == [[e, 0] for e in range(6)]
    assert lattice["allowed_classes"] == 6
    assert lattice["forbidden_classes"] == 30
    assert lattice["magnetic_sector"][
        "z6_quotient_pure_magnetic_classes_allowed"
    ] == [0, 1, 2, 3, 4, 5]


def test_maximal_isotropic_enumeration_independent() -> None:
    from itertools import combinations

    elements = [(e, m) for e in range(6) for m in range(6)]

    def pairing(a, b):
        return (a[0] * b[1] - b[0] * a[1]) % 6

    def closure(generators):
        group = {(0, 0)}
        frontier = list(generators)
        while frontier:
            g = frontier.pop()
            if g in group:
                continue
            group.add(g)
            for h in list(group):
                s = ((g[0] + h[0]) % 6, (g[1] + h[1]) % 6)
                if s not in group:
                    frontier.append(s)
        return frozenset(group)

    subgroups = {closure(pair) for pair in combinations(elements, 2)}
    subgroups |= {closure((g,)) for g in elements}
    maximal = []
    for group in subgroups:
        if len(group) != 6:
            continue
        if any(pairing(a, b) != 0 for a in group for b in group):
            continue
        if all(
            any(pairing(x, a) != 0 for a in group)
            for x in elements
            if x not in group
        ):
            maximal.append(group)
    assert len(maximal) == 12
    assert frozenset((0, m) for m in range(6)) in maximal
    assert frozenset((e, 0) for e in range(6)) in maximal
    assert frozenset((e, m) for e in (0, 2, 4) for m in (0, 3)) in maximal
    assert frozenset((e, m) for e in (0, 3) for m in (0, 2, 4)) in maximal


def test_dyonic_pairing_rejects_a_tilted_lattice() -> None:
    def pairing(e1: int, m1: int, e2: int, m2: int) -> int:
        return (e1 * m2 - e2 * m1) % 6

    tilted = [(m % 2, m) for m in range(6)]
    isotropic = all(
        pairing(*left, *right) == 0 for left in tilted for right in tilted
    )
    assert isotropic is False


def test_registry_slot_extraction_counts() -> None:
    registry = sv._load_registry()
    assert len(sv._slot_candidates(registry, "z6_charge_line_congruences")) == 5
    assert len(sv._slot_candidates(registry, "a5_angular_rules")) == 6
    assert len(sv._slot_candidates(registry, "wz_scale_free_response")) == 4
    assert len(sv._slot_candidates(registry, "observer_overlap_cross_spectra")) == 3


def test_committed_verdicts_are_byte_exact() -> None:
    for issue, verdict in sv.build_all().items():
        committed = (sv.RUNTIME / f"stratum_verdict_{issue}.json").read_bytes()
        assert committed == sv.canonical_json_bytes(verdict)


def test_registry_drift_fails_closed(tmp_path, monkeypatch: pytest.MonkeyPatch) -> None:
    tampered = json.loads(sv.REGISTRY_PATH.read_text(encoding="utf-8"))
    tampered["scoring_boundary"]["comparison_access_permitted"] = True
    path = tmp_path / "candidate_registry.json"
    path.write_text(json.dumps(tampered), encoding="utf-8")
    monkeypatch.setattr(sv, "REGISTRY_PATH", path)
    with pytest.raises(sv.VerdictError, match="scoring boundary drift"):
        sv.build_all()


def test_wz_promotion_drift_fails_closed(tmp_path, monkeypatch: pytest.MonkeyPatch) -> None:
    tampered = json.loads(sv.WZ_MANIFEST_PATH.read_text(encoding="utf-8"))
    tampered["promotion_allowed"] = True
    path = tmp_path / "INTEGRATION_MANIFEST.json"
    path.write_text(json.dumps(tampered), encoding="utf-8")
    monkeypatch.setattr(sv, "WZ_MANIFEST_PATH", path)
    registry = sv._load_registry()
    with pytest.raises(sv.VerdictError, match="promotion drift"):
        sv.build_646(registry)
