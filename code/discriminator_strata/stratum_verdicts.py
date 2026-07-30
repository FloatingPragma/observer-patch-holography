#!/usr/bin/env python3
"""Typed closure verdicts for the discriminator producer strata.

Issues #642, #643, #645, and #646 are producer strata inside the issue #647
invariant-mining registry. Their generation contributions are complete, the
registry is sealed, and each stratum's remaining outcome is a typed verdict
that consumes the registry and the pinned source receipts without opening
any comparison surface. This producer emits the four verdicts:

* issue #642: the diagonal Z6 packet converts to an exact conditional kill
  packet. The registry carries the descent congruence, the charge
  corollaries, the flux refinement, and the electric line-class exclusion;
  this verdict adds the exact dyonic line-lattice arithmetic on the Z6
  charge torus and types the physical line-operator attachment open.
* issue #643: the angular selection rules are registered as the source
  template, and the internal-versus-rotation-remnant question is not
  decidable on the declared source interface because no screen-to-sky map
  exists; the forced nonzero sky statistic is not constructible without it.
* issue #645: the normalized overlap cross-spectrum exists as registered
  source structure, and no map from the declared source interface to an
  interferometer readout exists, so the comparison contract is correctly
  not frozen.
* issue #646: the physical pole, width, residue, and asymmetry vector is
  not source-evaluable on the declared interface (no unit attachment and no
  pole promotion), and the surviving scale-free combinations are exactly
  the registered candidates.

Every verdict pins its inputs by exact bytes, fails closed on drift, reads
no public measurement, and permits no comparison.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
RUNTIME = HERE / "runtime"
REGISTRY_PATH = (
    REPO_ROOT / "code" / "invariant_mining" / "outputs" / "candidate_registry.json"
)
WZ_MANIFEST_PATH = (
    REPO_ROOT
    / "code"
    / "particles"
    / "calibration"
    / "wz_upstream_completion"
    / "INTEGRATION_MANIFEST.json"
)

VERDICTS = {
    642: "CONDITIONAL_KILL_PACKET__PHYSICAL_LINE_ATTACHMENT_OPEN",
    643: "NO_SCREEN_TO_SKY_MAP__SOURCE_TEMPLATE_REGISTERED",
    645: "NO_SOURCE_TO_READOUT_MAP__SOURCE_SPECTRUM_REGISTERED",
    646: "TESTED_VECTOR_NOT_SOURCE_EVALUABLE__SCALE_FREE_SURVIVORS_REGISTERED",
}


class VerdictError(ValueError):
    """A stratum verdict refused to build."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerdictError(message)


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")


def tagged_sha256(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def _pin(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    return {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "bytes": len(payload),
        "sha256": tagged_sha256(payload),
    }


def _load_registry() -> dict[str, Any]:
    registry = json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))
    require(
        registry.get("schema") == "oph.invariant_mining.candidate_registry.v1",
        "candidate registry schema drift",
    )
    require(
        registry.get("status") == "GENERATION_COMPLETE__SCORING_SEALED",
        "candidate registry status drift",
    )
    boundary = registry.get("scoring_boundary", {})
    require(
        boundary.get("physical_scoring_permitted") is False
        and boundary.get("comparison_access_permitted") is False,
        "candidate registry scoring boundary drift",
    )
    return registry


def _slot_candidates(registry: dict[str, Any], slot_id: str) -> list[dict[str, Any]]:
    rows = [
        candidate
        for candidate in registry["candidates"]
        if candidate["slot_id"] == slot_id
    ]
    require(bool(rows), f"no registered candidates for slot {slot_id}")
    return rows


def _candidate_digest(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "candidate_id": row["candidate_id"],
            "grammar_class": row["grammar_class"],
            "rank": row["rank"],
            "score": row["score"],
            "kill_rule": row["kill_rule"],
        }
        for row in sorted(rows, key=lambda row: row["rank"])
    ]


# ---------------------------------------------------------------------------
# Issue #642: exact dyonic line-lattice arithmetic on the Z6 charge torus
# ---------------------------------------------------------------------------


def dyonic_line_lattice_certificate() -> dict[str, Any]:
    """Exact line-class arithmetic for the quotient global form.

    Line classes live on the Z6 x Z6 charge torus with the Dirac pairing
    ``<(e1, m1), (e2, m2)> = e1*m2 - e2*m1 mod 6``. The quotient form's
    electric classes are the descended (zero) classes, so its line lattice
    is ``{(0, m)}``. The certificate checks by exhaustion that this lattice
    is isotropic, that it is maximal (every outside class fails mutual
    locality with some lattice element), and that exactly six of the
    thirty-six classes are allowed.
    """

    def pairing(e1: int, m1: int, e2: int, m2: int) -> int:
        return (e1 * m2 - e2 * m1) % 6

    lattice = [(0, m) for m in range(6)]
    isotropic = all(
        pairing(*left, *right) == 0 for left in lattice for right in lattice
    )
    outside_all_fail = all(
        any(pairing(e, m, *element) != 0 for element in lattice)
        for e in range(6)
        for m in range(6)
        if (e, m) not in lattice
    )
    allowed = [[e, m] for e in range(6) for m in range(6) if (e, m) in lattice]
    forbidden_count = 36 - len(allowed)
    return {
        "charge_torus": "Z6 x Z6",
        "dirac_pairing": "e1*m2 - e2*m1 mod 6",
        "lattice": allowed,
        "lattice_isotropic": isotropic,
        "lattice_maximal_by_exhaustion": outside_all_fail,
        "allowed_classes": len(allowed),
        "forbidden_classes": forbidden_count,
        "theta_dependence": (
            "the dyonic tilt of the lattice under a vacuum angle is a named "
            "open direction owned by the strong-CP issue; this certificate "
            "fixes the untilted normalization"
        ),
    }


def build_642(registry: dict[str, Any]) -> dict[str, Any]:
    rows = _slot_candidates(registry, "z6_charge_line_congruences")
    lattice = dyonic_line_lattice_certificate()
    require(
        lattice["lattice_isotropic"]
        and lattice["lattice_maximal_by_exhaustion"]
        and lattice["allowed_classes"] == 6,
        "dyonic lattice certificate drift",
    )
    return {
        "schema": "oph.discriminator_stratum_verdict.v1",
        "issue": 642,
        "status": VERDICTS[642],
        "registered_kill_rules": _candidate_digest(rows),
        "dyonic_line_lattice": lattice,
        "kill_rule_summary": [
            "a physical state violating the descent congruence",
            "a free fractionally charged color singlet",
            "a color triplet outside the descending charge classes",
            "an electric line class outside the descended lattice",
            "a dyonic line class outside the maximal isotropic lattice",
            "a flux and line refinement with index other than six",
        ],
        "physical_attachment": (
            "OPEN: the line-operator and flux observables have no laboratory "
            "attachment map on the declared source interface; every kill "
            "rule is conditional on the named z6 global-form branch"
        ),
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
        },
        "reopen_condition": (
            "a physical line-operator or flux attachment map, at which point "
            "the registered kill rules become scorable through the issue "
            "#639 protocol"
        ),
    }


def build_643(registry: dict[str, Any]) -> dict[str, Any]:
    rows = _slot_candidates(registry, "a5_angular_rules")
    return {
        "schema": "oph.discriminator_stratum_verdict.v1",
        "issue": 643,
        "status": VERDICTS[643],
        "registered_source_template": _candidate_digest(rows),
        "decision_boundary": (
            "the declared source interface carries the A5 action on the "
            "port coefficient module as internal algebra; deciding a "
            "physical rotation-remnant reading requires a screen-to-sky map "
            "that no registered source observable supplies, so the question "
            "is not decidable on the declared interface"
        ),
        "forced_sky_statistic": (
            "not constructible: a forced nonzero sky statistic requires the "
            "missing map, so no comparison contract is frozen"
        ),
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
        },
        "reopen_condition": (
            "a source-derived screen-to-sky template map; the registered "
            "angular candidates are the retained source template for that "
            "event"
        ),
    }


def build_645(registry: dict[str, Any]) -> dict[str, Any]:
    rows = _slot_candidates(registry, "observer_overlap_cross_spectra")
    return {
        "schema": "oph.discriminator_stratum_verdict.v1",
        "issue": 645,
        "status": VERDICTS[645],
        "registered_source_spectrum": _candidate_digest(rows),
        "feasibility_result": (
            "the normalized cross-spectrum exists as registered source "
            "structure: the unique flat twist sector, the conjugate-sector "
            "symmetry, and the ladder stability, all bound to the pinned "
            "classical realization receipt"
        ),
        "readout_boundary": (
            "no map from the declared source interface to an interferometer "
            "readout exists: the seam physicalization and laboratory "
            "attachment are open producer obligations, so the comparison "
            "contract is correctly not frozen"
        ),
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
        },
        "reopen_condition": (
            "a physical seam-readout attachment; the registered spectrum "
            "rows are the retained source side of that future contract"
        ),
    }


def build_646(registry: dict[str, Any]) -> dict[str, Any]:
    rows = _slot_candidates(registry, "wz_scale_free_response")
    manifest = json.loads(WZ_MANIFEST_PATH.read_text(encoding="utf-8"))
    require(
        manifest.get("promotion_allowed") is False,
        "wz integration manifest promotion drift",
    )
    return {
        "schema": "oph.discriminator_stratum_verdict.v1",
        "issue": 646,
        "status": VERDICTS[646],
        "tested_vector": [
            "pole",
            "width",
            "residue",
            "asymmetry",
        ],
        "non_evaluability": (
            "the physical vector is not source-evaluable on the declared "
            "interface: physical units are not evaluable on the serialized "
            "local domain, the upstream W/Z stack forbids pole promotion, "
            "and the clock, vacuum-scale, normalization, threshold, and "
            "scheme directions are unresolved nuisance directions in the "
            "frozen registry"
        ),
        "scale_free_survivors": _candidate_digest(rows),
        "survivor_boundary": (
            "the surviving scale-free combinations are exactly the "
            "registered candidates; the full source-to-pole campaign stays "
            "gated behind its activation conditions"
        ),
        "wz_manifest_status": manifest.get("scientific_status"),
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
        },
        "reopen_condition": (
            "a positive unit and pole producer chain, at which point the "
            "minimal surviving source inputs are the registered scale-free "
            "combinations"
        ),
    }


BUILDERS = {642: build_642, 643: build_643, 645: build_645, 646: build_646}


def build_all() -> dict[int, dict[str, Any]]:
    registry = _load_registry()
    pins = [_pin(REGISTRY_PATH), _pin(WZ_MANIFEST_PATH)]
    verdicts = {}
    for issue, builder in sorted(BUILDERS.items()):
        verdict = builder(registry)
        verdict["parent_pins"] = pins
        verdict["verdict_sha256"] = tagged_sha256(canonical_json_bytes(verdict))
        verdicts[issue] = verdict
    return verdicts


def write_runtime() -> list[Path]:
    RUNTIME.mkdir(exist_ok=True)
    paths = []
    for issue, verdict in build_all().items():
        path = RUNTIME / f"stratum_verdict_{issue}.json"
        path.write_bytes(canonical_json_bytes(verdict))
        paths.append(path)
    return paths


def verify_runtime() -> None:
    for issue, verdict in build_all().items():
        path = RUNTIME / f"stratum_verdict_{issue}.json"
        if path.read_bytes() != canonical_json_bytes(verdict):
            raise SystemExit(f"stratum verdict {issue} is stale")
        committed = json.loads(path.read_text(encoding="utf-8"))
        if committed["status"] != VERDICTS[issue]:
            raise SystemExit(f"stratum verdict {issue} status drift")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write:
        for path in write_runtime():
            print(path)
    if args.verify:
        verify_runtime()
        print("STRATUM_VERDICTS_VALID")
    if not args.write and not args.verify:
        print(
            canonical_json_bytes(
                {issue: verdict["status"] for issue, verdict in build_all().items()}
            ).decode("ascii"),
            end="",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
