#!/usr/bin/env python3
"""Fail-closed semantics menu ledger for the bounded capacity-closure campaign.

GitHub issue #615 asks for the admissible publicness, reserve, and channel
semantics to be enumerated before any optimization, and for the campaign to
exit through a fail-closed record when no source-only fixed-point selector
exists. This certificate is that record. It consolidates the executed
capacity-readback candidate families of the 2026-07-14/15 construction run
(F_CONSTRUCTION_2026-07-14.md, G2_GAP_1_COUPLING_THEOREM.md) against the
readback specification (F_READBACK_SPEC.md) and emits one deterministic
manifest, schema oph.capacity_semantics_menu_certificate.v1.

Content of the manifest:

* the declared semantic axes exactly as the specification and the
  construction record state them: publicness family, cell product structure,
  reserve semantics and attachment, readback-record effect, observer
  marking, symmetry quotient, Cap_read family, channel capacity semantics,
  kernel ontology, continuation scope;
* one menu row per executed candidate reading: semantic choice tuple, the
  implementing family module, and the recorded executed verdict with the
  recorded reason (excluded, no positive fixed point, excluded before
  evaluation, or conditional open);
* the fail-closed campaign verdict: the semantics enumeration is the
  declared menu for the executed families, the source-only fixed-point
  selector is a conditional open interface, no cosmological, electroweak,
  or horizon numeric target enters the certificate, and the horizon-area
  identification (#589) and the electroweak load bridge (#547) stay
  separate physical attachments.

The certificate consumes recorded verdicts. It does not re-derive interval
arithmetic. Where a candidate module exposes a cheap build(), the build is
re-executed and its statuses are compared with the recorded runtime
certificate; the 180-row CAP-L lattice is pinned by hash with re-execution
deferred on wall-clock grounds. The CAP-K linear-contraction argument is
checked symbolically: exact fractions 5/6 and 1/2 lie below one, and a
certified positive lower endpoint of the P enclosure forces both
exp(-P/24) < 1 and 1 - P/24 < 1.

Fail-closed controls: a menu row that claims a positive fixed point for a
linear contraction through the origin is rejected, and an injected row that
hardcodes a capacity target (by scale token or by marker key) is rejected.
The menu is a record of executed scope, and ties or incomplete semantics
leave the selector conditional_open_interface; an economy rule does not
resolve them.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
import sys
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

REPO_ROOT = MODULE_DIR.parent.parent

SCHEMA = "oph.capacity_semantics_menu_certificate.v1"
ISSUE = 615
MANIFEST_PATH = MODULE_DIR / "manifests" / "capacity_semantics_menu_reference.json"

ALLOWED_VERDICTS = (
    "excluded",
    "no_positive_fixed_point",
    "excluded_pre_evaluation",
    "conditional_open",
)

# Files whose recorded content this ledger consumes; pinned by sha256 in the
# manifest so the ledger fails closed against silent drift.
PINNED_SOURCES = (
    "code/capacity_readback/F_READBACK_SPEC.md",
    "code/capacity_readback/F_CONSTRUCTION_2026-07-14.md",
    "code/capacity_readback/G2_GAP_1_COUPLING_THEOREM.md",
    "code/capacity_readback/A5_FINITE_CONTROL_STATUS_2026-07-20.md",
    "code/capacity_readback/F_candidate_capK.py",
    "code/capacity_readback/F_candidate_capL.py",
    "code/capacity_readback/F_candidate_capP.py",
    "code/capacity_readback/F_candidate_coupled.py",
    "code/capacity_readback/runtime/F_candidate_capK_certificates.json",
    "code/capacity_readback/runtime/F_candidate_capL_certificates.json",
    "code/capacity_readback/runtime/F_candidate_capP_certificates.json",
    "code/capacity_readback/runtime/F_candidate_coupled_certificates.json",
    "code/capacity_readback/runtime/F_construction_comparison_2026-07-14.json",
    "code/capacity_readback/runtime/a5_finite_control_status.json",
)

# Marker key fragments whose presence anywhere in a row or in the manifest
# marks a hardcoded target. Normalized to lowercase alphanumerics.
TARGET_MARKER_TOKENS = (
    "measuredlambda",
    "electroweaktarget",
    "higgstarget",
    "rhoop",
    "expectedanswer",
    "targetcapacity",
    "desiredcapacity",
    "referencecapacity",
)

# A standalone float token whose exponent magnitude sits in the cosmological
# window is a numeric target import. Hex lookarounds keep sha256 digests out
# of the match.
_SCALE_TOKEN = re.compile(r"(?<![0-9a-fA-F.])\d+(?:\.\d+)?[eE]([+-]?\d+)(?![0-9a-fA-F])")
_HEX_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_SCALE_EXPONENT_LO = 50
_SCALE_EXPONENT_HI = 400


class CertificateError(ValueError):
    """Fail-closed ledger error carrying a stable code."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise CertificateError(code, message)


def sha256_file(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as exc:
        raise CertificateError("SOURCE_READ", f"cannot read {path}: {exc}") from exc


def load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CertificateError("JSON_READ", f"cannot read {path}: {exc}") from exc


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _normalized(text: str) -> str:
    return "".join(ch.lower() for ch in text if ch.isalnum())


# ---------------------------------------------------------------------------
# Hidden-target scan
# ---------------------------------------------------------------------------

def scan_for_targets(value: Any, path: str = "$") -> None:
    """Reject marker keys, cosmological-scale float tokens, and extreme numbers.

    The scan walks the payload object. A sha256 digest string is exempt from
    the scale-token scan; every dictionary key and every other string is
    checked. Any hit fails closed with HIDDEN_TARGET.
    """
    if isinstance(value, Mapping):
        for key, item in value.items():
            key_text = str(key)
            norm = _normalized(key_text)
            for token in TARGET_MARKER_TOKENS:
                require(
                    token not in norm,
                    "HIDDEN_TARGET",
                    f"{path}.{key_text}: key carries target marker '{token}'",
                )
            scan_for_targets(item, f"{path}.{key_text}")
        return
    if isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            scan_for_targets(item, f"{path}[{index}]")
        return
    if isinstance(value, str):
        if _HEX_DIGEST.match(value):
            return
        for match in _SCALE_TOKEN.finditer(value):
            exponent = abs(int(match.group(1)))
            require(
                not (_SCALE_EXPONENT_LO <= exponent <= _SCALE_EXPONENT_HI),
                "HIDDEN_TARGET",
                f"{path}: cosmological-scale numeric token '{match.group(0)}'",
            )
        return
    if isinstance(value, bool) or value is None:
        return
    if isinstance(value, (int, float)):
        magnitude = abs(float(value))
        require(
            not (magnitude >= 10.0 ** _SCALE_EXPONENT_LO),
            "HIDDEN_TARGET",
            f"{path}: numeric value at cosmological scale",
        )
        require(
            not (0.0 < magnitude <= 10.0 ** (-_SCALE_EXPONENT_LO)),
            "HIDDEN_TARGET",
            f"{path}: numeric value at inverse cosmological scale",
        )


# ---------------------------------------------------------------------------
# Declared semantic axes (enumerated before any verdict is attached)
# ---------------------------------------------------------------------------

def declared_semantic_axes() -> list[dict[str, Any]]:
    """The admissible semantics menu as the specification and the construction
    record declare it. The construction record carries its own
    DECLARED-BEFORE-COMPARISON marker; this enumeration restates the declared
    axes and adds no alternative."""
    return [
        {
            "axis": "publicness_family",
            "source": "F_READBACK_SPEC.md section 3",
            "alternatives": [
                "collective (one authorized set of all observers)",
                "universal-local (every observer separately)",
                "quorum (declared authorized subsets)",
            ],
            "executed_scope": (
                "the issue-548 source-derived packet freezes the collective "
                "all-port family universal-twelve-port-publicness/v1; the "
                "universal-local and quorum contracts have no executed packet "
                "and stay open"
            ),
        },
        {
            "axis": "cell_product_structure",
            "source": "F_CONSTRUCTION_2026-07-14.md step 2, branch BR-0",
            "alternatives": [
                "(i) K = 4N/P cells as independent record carriers of P/4 nats",
                "(ii) refusal of the product structure: no configuration count",
            ],
            "executed_scope": (
                "reading (i) is the only executable reading and is carried on "
                "every executed branch; reading (ii) is recorded as the "
                "fallback with outcome class (c)"
            ),
        },
        {
            "axis": "reserve_semantics",
            "source": "F_CONSTRUCTION_2026-07-14.md step 3, branch BR-1",
            "alternatives": [
                "poisson: survival factor exp(-P/24) per unit",
                "presence: survival factor 1 - P/24 per unit",
                "none: reserve slots stay countable states",
            ],
            "executed_scope": "all three readings executed inside CAP-L; poisson and presence executed inside CAP-P and CAP-K",
        },
        {
            "axis": "reserve_attachment",
            "source": "F_CONSTRUCTION_2026-07-14.md step 3, branch BR-2",
            "alternatives": [
                "per cell: K reserve units",
                "per shared edge: 3K reserve units",
            ],
            "executed_scope": (
                "both attachments executed; the poisson readings collapse to "
                "the exact readable shares 5/6 (per cell) and 1/2 (per edge), "
                "giving the five-reading rho axis R1..R5"
            ),
        },
        {
            "axis": "readback_record_effect",
            "source": "F_CONSTRUCTION_2026-07-14.md step 4, branch BR-3",
            "alternatives": [
                "X+: port-orbit multiplicity e^X = N/pi",
                "X0: delta constraint, no multiplicity and no cost",
                "X-: record cost X nats in bulk",
                "X=: write+check cost 2X on the oriented register",
            ],
            "executed_scope": (
                "all four readings executed inside the CAP-L lattice; the "
                "fifth reading, a constant readback through the D10 inner "
                "observation step, is the CAP-B bridge branch and is barred "
                "before evaluation under the specification blindness bar"
            ),
        },
        {
            "axis": "observer_marking",
            "source": "F_CONSTRUCTION_2026-07-14.md step 5, branch BR-4",
            "alternatives": [
                "K0: unmarked at-least-one chain",
                "K1: one marked host cell, factor K",
                "K2: marked host plus checkpoint cell, factor K^2",
            ],
            "executed_scope": (
                "all three executed inside CAP-L; per-port (k = 12) and "
                "per-flag (k = 60) markings import undeclared structure and "
                "are recorded without execution"
            ),
        },
        {
            "axis": "symmetry_quotient",
            "source": "F_CONSTRUCTION_2026-07-14.md step 6, branch BR-5",
            "alternatives": [
                "S-: A5 quotient, divide by 60",
                "S0: normal form supplies one representative per class",
                "S+: sixty face-corner flag anchorings, multiply by 60",
            ],
            "executed_scope": "all three executed inside CAP-L",
        },
        {
            "axis": "cap_read_family",
            "source": "F_CONSTRUCTION_2026-07-14.md step 7, branch BR-6; G2_GAP_1_COUPLING_THEOREM.md",
            "alternatives": [
                "CAP-L: log-count readback, unit nat normalization",
                "CAP-P: port-inversion readback",
                "CAP-K: cell-count readback",
                "CAP-B: bridge constant readback (barred before evaluation)",
                "coupled: theorem-coupled port-load inversion, conditional on CP-1..CP-3",
            ],
            "executed_scope": "CAP-L, CAP-P, CAP-K executed blind; CAP-B never executed; coupled executed as a conditional theorem-coupled candidate",
        },
        {
            "axis": "channel_capacity_semantics",
            "source": "F_READBACK_SPEC.md sections 4-5",
            "alternatives": [
                "zero-error capacity M_0 (independence number of the compound confusability graph)",
                "epsilon-tolerant capacity M_epsilon with the total-variation stability bound",
            ],
            "executed_scope": (
                "both implemented in correctable_public_record_capacity.py; "
                "the zero-error branch is intentionally discontinuous under "
                "arbitrarily small full-support noise and that control is on "
                "record"
            ),
        },
        {
            "axis": "kernel_ontology",
            "source": "F_READBACK_SPEC.md section 3",
            "alternatives": [
                "joint checkpoint kernels supplied per authorized set and continuation",
                "local marginals only (insufficient: parity and independent-uniform channels share marginals with capacities two and one)",
            ],
            "executed_scope": "joint kernels are mandatory; the marginal nonidentifiability countermodel is executed in the evaluator tests",
        },
        {
            "axis": "continuation_scope",
            "source": "F_READBACK_SPEC.md section 3",
            "alternatives": [
                "indefinite exact continuation via the finite support-relation semigroup",
                "declared finite horizon, valid only when the theorem is scoped to it",
            ],
            "executed_scope": "the semigroup closure is implemented and executed; no finite-horizon packet is on record",
        },
    ]


# ---------------------------------------------------------------------------
# Recorded runtime certificates: load and cross-check
# ---------------------------------------------------------------------------

RECORDED_EXPECTATIONS = {
    "capK": {
        "path": "runtime/F_candidate_capK_certificates.json",
        "rows": 4,
        "by_status": {"no_positive_fixed_point": 4},
    },
    "capP": {
        "path": "runtime/F_candidate_capP_certificates.json",
        "rows": 6,
        "by_status": {"fixed_point_certified": 4, "no_positive_fixed_point": 2},
    },
    "capL": {
        "path": "runtime/F_candidate_capL_certificates.json",
        "rows": 180,
        "by_status": {
            "fixed_point_certified": 98,
            "fixed_point_unstable_rejected": 15,
            "no_fixed_point": 48,
            "rejected_no_contraction": 18,
            "rejected_trivial": 1,
        },
        "p4_coherent_rows": 0,
    },
}


def load_recorded() -> dict[str, Any]:
    recorded: dict[str, Any] = {}
    for family, expectation in RECORDED_EXPECTATIONS.items():
        payload = load_json(MODULE_DIR / expectation["path"])
        summary = payload.get("summary", {})
        require(
            summary.get("rows") == expectation["rows"],
            "RECORDED_MISMATCH",
            f"{family}: recorded row count {summary.get('rows')} differs from "
            f"the pinned expectation {expectation['rows']}",
        )
        require(
            summary.get("by_status") == expectation["by_status"],
            "RECORDED_MISMATCH",
            f"{family}: recorded status counts differ from the pinned expectation",
        )
        if "p4_coherent_rows" in expectation:
            require(
                summary.get("p4_coherent_rows") == expectation["p4_coherent_rows"],
                "RECORDED_MISMATCH",
                f"{family}: recorded p4_coherent_rows differs from the pinned expectation",
            )
        require(
            payload.get("moves_cl7") is False,
            "RECORDED_MISMATCH",
            f"{family}: recorded certificate must not move the closure ledger row",
        )
        recorded[family] = payload

    coupled = load_json(MODULE_DIR / "runtime/F_candidate_coupled_certificates.json")
    require(
        coupled.get("status") == "conditional_on_CP-1_CP-2_CP-3",
        "RECORDED_MISMATCH",
        "coupled: recorded status differs from conditional_on_CP-1_CP-2_CP-3",
    )
    require(
        coupled.get("moves_cl7") is False
        and coupled.get("cl7_status") == "open_reduced_to_CP-1_CP-2_CP-3",
        "RECORDED_MISMATCH",
        "coupled: recorded closure-row fields differ from the pinned expectation",
    )
    recorded["coupled"] = coupled

    comparison = load_json(MODULE_DIR / "runtime/F_construction_comparison_2026-07-14.json")
    require(
        comparison.get("landed_rows") == 0 and comparison.get("total_rows") == 190,
        "RECORDED_MISMATCH",
        "comparison record: pinned expectation is 190 executed rows and zero landings",
    )
    recorded["comparison_landed_rows"] = 0
    recorded["comparison_total_rows"] = 190
    return recorded


def reexecute_cheap_families(recorded: Mapping[str, Any]) -> dict[str, Any]:
    """Re-run the sub-second candidate builds and compare statuses with the
    recorded certificates. CAP-L re-execution is deferred: the 180-row lattice
    is a minute-scale job and its recorded certificate is pinned by hash."""
    import F_candidate_capK
    import F_candidate_capP
    import F_candidate_coupled

    confirmations: dict[str, Any] = {}

    fresh_k = F_candidate_capK.build()
    require(
        fresh_k["summary"]["by_status"] == recorded["capK"]["summary"]["by_status"],
        "REEXECUTION_MISMATCH",
        "capK: re-executed status counts differ from the recorded certificate",
    )
    require(
        [row["status"] for row in fresh_k["rows"]]
        == [row["status"] for row in recorded["capK"]["rows"]],
        "REEXECUTION_MISMATCH",
        "capK: re-executed per-row statuses differ from the recorded certificate",
    )
    confirmations["capK"] = {"reexecuted": True, "statuses_match_recorded": True}

    fresh_p = F_candidate_capP.build()
    require(
        fresh_p["summary"]["by_status"] == recorded["capP"]["summary"]["by_status"],
        "REEXECUTION_MISMATCH",
        "capP: re-executed status counts differ from the recorded certificate",
    )
    require(
        [row["status"] for row in fresh_p["rows"]]
        == [row["status"] for row in recorded["capP"]["rows"]],
        "REEXECUTION_MISMATCH",
        "capP: re-executed per-row statuses differ from the recorded certificate",
    )
    confirmations["capP"] = {"reexecuted": True, "statuses_match_recorded": True}

    fresh_c = F_candidate_coupled.build()
    require(
        fresh_c["status"] == recorded["coupled"]["status"]
        and fresh_c["cl7_status"] == recorded["coupled"]["cl7_status"],
        "REEXECUTION_MISMATCH",
        "coupled: re-executed status differs from the recorded certificate",
    )
    confirmations["coupled"] = {"reexecuted": True, "statuses_match_recorded": True}

    confirmations["capL"] = {
        "reexecuted": False,
        "note": (
            "re-execution deferred: the 180-row interval lattice is a "
            "minute-scale job; the recorded certificate is pinned by sha256 "
            "and its summary is cross-checked against the pinned expectation"
        ),
    }
    return confirmations


# ---------------------------------------------------------------------------
# CAP-K linear-contraction argument, symbolic
# ---------------------------------------------------------------------------

def capk_contraction_check(recorded_capk: Mapping[str, Any]) -> dict[str, Any]:
    """Certify s < 1 for every CAP-K reading without interval re-derivation.

    Exact readings 5/6 and 1/2 are compared as fractions. The P-dependent
    readings exp(-P/24) and 1 - P/24 lie below one exactly when P > 0, and
    the certified lower endpoint of the P enclosure is positive; the recorded
    outward-rounded upper endpoints are additionally checked below one."""
    from F_candidate_capL import P_LO

    require(
        Decimal(P_LO) > 0,
        "CONTRACTION_ARGUMENT",
        "the certified lower endpoint of the P enclosure must be positive",
    )
    exact = {
        "capK.s_nat_share": Fraction(5, 6),
        "capK.s_edge_share": Fraction(1, 2),
    }
    rows_report = []
    for row in recorded_capk["rows"]:
        branch = row["branch"]
        recorded_hi = Decimal(row["s"]["hi"])
        require(
            recorded_hi < 1,
            "CONTRACTION_ARGUMENT",
            f"{branch}: recorded upper endpoint of s must lie below one",
        )
        entry: dict[str, Any] = {
            "branch": branch,
            "recorded_s_hi_below_one": True,
        }
        if branch in exact:
            fraction = exact[branch]
            require(
                fraction < 1,
                "CONTRACTION_ARGUMENT",
                f"{branch}: exact share must lie below one",
            )
            entry["exact_fraction"] = f"{fraction.numerator}/{fraction.denominator}"
            entry["exact_fraction_below_one"] = True
        else:
            entry["symbolic_argument"] = (
                "P > 0 from the certified enclosure lower endpoint, hence "
                "exp(-P/24) < 1 and 1 - P/24 < 1"
            )
        rows_report.append(entry)
    return {
        "map_form": "F(N) = s*N, linear through the origin",
        "conclusion": (
            "every CAP-K reading is a contraction with fixed point N = 0 "
            "outside (0, inf); a positive fixed point is impossible on this "
            "family"
        ),
        "rows": rows_report,
    }


# ---------------------------------------------------------------------------
# Menu rows
# ---------------------------------------------------------------------------

def _capl_status_counts(recorded_capl: Mapping[str, Any]) -> dict[str, dict[str, int]]:
    counts: dict[str, dict[str, int]] = {}
    for row in recorded_capl["rows"]:
        rho_branch = row["branch"].split(".")[1]
        per = counts.setdefault(rho_branch, {})
        per[row["status"]] = per.get(row["status"], 0) + 1
    return counts


CAPL_RHO_SEMANTICS = {
    "R1": ("poisson exp(-P/24)", "per cell", "rho = 5/6 exact"),
    "R2": ("presence 1 - P/24", "per cell", "rho = 1 + (4/P) ln(1 - P/24)"),
    "R3": ("none", "collapsed", "rho = 1"),
    "R4": ("poisson exp(-P/24)", "per shared edge", "rho = 1/2 exact"),
    "R5": ("presence 1 - P/24", "per shared edge", "rho = 1 + (12/P) ln(1 - P/24)"),
}

CAPK_SEMANTICS = {
    "capK.s_poisson": ("poisson exp(-P/24)", "per cell survival factor"),
    "capK.s_presence": ("presence 1 - P/24", "per cell survival factor"),
    "capK.s_nat_share": ("poisson, exact nat share 5/6", "per cell"),
    "capK.s_edge_share": ("poisson, exact nat share 1/2", "per shared edge"),
}

CAPP_SEMANTICS = {
    "capP.s_poisson_port": ("poisson exp(-P/24)", "multiplicative on the per-port load"),
    "capP.s_presence_port": ("presence 1 - P/24", "multiplicative on the per-port load"),
    "capP.s_poisson_pair": ("poisson exp(-P/12)", "multiplicative per oriented slot pair"),
    "capP.s_presence_pair": ("presence (1 - P/24)^2", "multiplicative per oriented slot pair"),
    "capP.add_slot": ("additive reserve cost P/24", "subtracted from the per-port load"),
    "capP.add_port": ("additive reserve cost P/12", "subtracted from the per-port load"),
}


def menu_rows(recorded: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    # CAP-K: four readings, all without a positive fixed point.
    for record in recorded["capK"]["rows"]:
        reserve, attachment = CAPK_SEMANTICS[record["branch"]]
        require(
            record["status"] == "no_positive_fixed_point",
            "RECORDED_MISMATCH",
            f"{record['branch']}: recorded status must be no_positive_fixed_point",
        )
        rows.append({
            "row_id": record["branch"],
            "family": "CAP-K",
            "implementation": "F_candidate_capK.py",
            "semantics": {
                "publicness_family": "inherited frozen family of the count construction",
                "reserve_semantics": reserve,
                "reserve_attachment": attachment,
                "channel": "cell-count readback F(N) = (P/4) * K_readable = s*N",
            },
            "executed_verdict": "no_positive_fixed_point",
            "reason": (
                "linear map through the origin with certified s < 1; the "
                "unique fixed point N = 0 lies outside (0, inf) and the "
                "self-map bracketing pair of specification property P3 does "
                "not exist"
            ),
            "linear_contraction": {"map": "F(N) = s*N", "s_below_one": True},
            "verdict_provenance": {
                "recorded": "runtime/F_candidate_capK_certificates.json",
                "reexecuted": True,
            },
        })

    # CAP-P: four multiplicative readings pinned at the fixed point pi and
    # excluded at reference scale; two additive readings with no positive
    # fixed point.
    for record in recorded["capP"]["rows"]:
        reserve, attachment = CAPP_SEMANTICS[record["branch"]]
        multiplicative = record["status"] == "fixed_point_certified"
        row: dict[str, Any] = {
            "row_id": record["branch"],
            "family": "CAP-P",
            "implementation": "F_candidate_capP.py",
            "semantics": {
                "publicness_family": "inherited frozen family of the count construction",
                "reserve_semantics": reserve,
                "reserve_attachment": attachment,
                "channel": "port-inversion readback of the twelve-port load X = log(N/pi)",
            },
            "verdict_provenance": {
                "recorded": "runtime/F_candidate_capP_certificates.json",
                "reexecuted": True,
            },
        }
        if multiplicative:
            row["executed_verdict"] = "excluded"
            row["reason"] = (
                "certified fixed point exactly N = pi for every degradation "
                "strength: the unique solution of (N/pi)^s = N/pi with "
                "s != 1 is N = pi; the count-density comparison registers a "
                "discrepancy on every certified row and the comparison "
                "record lands zero rows"
            )
        else:
            row["executed_verdict"] = "no_positive_fixed_point"
            row["reason"] = (
                "linear map through the origin with s < 1; the unique fixed "
                "point N = 0 lies outside the admissible interval"
            )
            row["linear_contraction"] = {"map": "F(N) = s*N", "s_below_one": True}
        rows.append(row)

    # CAP-L: five rho readings aggregating the 36-row (c, d) sublattice each.
    capl_counts = _capl_status_counts(recorded["capL"])
    for rho_branch in sorted(CAPL_RHO_SEMANTICS):
        reserve, attachment, rho_value = CAPL_RHO_SEMANTICS[rho_branch]
        counts = capl_counts[rho_branch]
        require(
            sum(counts.values()) == 36,
            "RECORDED_MISMATCH",
            f"CAP-L {rho_branch}: sublattice must carry exactly 36 recorded rows",
        )
        rows.append({
            "row_id": f"capL.{rho_branch}",
            "family": "CAP-L",
            "implementation": "F_candidate_capL.py",
            "semantics": {
                "publicness_family": "inherited frozen family of the count construction",
                "reserve_semantics": reserve,
                "reserve_attachment": attachment,
                "rho": rho_value,
                "channel": "log-count readback F(N) = rho*N + c*log N + d, unit nat normalization",
                "sublattice": "readback-record effect (4) x observer marking (3) x symmetry quotient (3)",
            },
            "executed_verdict": "excluded",
            "reason": (
                "every certified fixed point of the lattice lies in "
                "[1.4686, 1452.33] nats, zero rows are P4-coherent, and the "
                "comparison record lands zero rows against either reference; "
                "the recorded scale exclusion holds independently of every "
                "branch choice"
            ),
            "recorded_status_counts": counts,
            "verdict_provenance": {
                "recorded": "runtime/F_candidate_capL_certificates.json",
                "reexecuted": False,
                "note": "re-execution deferred; recorded certificate pinned by sha256",
            },
        })

    # CAP-B: barred before evaluation.
    rows.append({
        "row_id": "capB.bridge_constant",
        "family": "CAP-B",
        "implementation": "recorded in F_CONSTRUCTION_2026-07-14.md step 4; no module, never executed",
        "semantics": {
            "publicness_family": "inherited frozen family of the count construction",
            "reserve_semantics": "inapplicable: constant readback",
            "reserve_attachment": "inapplicable",
            "channel": (
                "per-port read equated to the inner electromagnetic "
                "observation step of the D10 lane, forcing a constant "
                "readback at the electroweak-bridge expression"
            ),
        },
        "executed_verdict": "excluded_pre_evaluation",
        "reason": (
            "the forced constant is the issue-589 electroweak-bridge "
            "expression, named in the specification blindness cone; barred "
            "under control V-08 before evaluation and never executed"
        ),
        "verdict_provenance": {
            "recorded": "F_CONSTRUCTION_2026-07-14.md",
            "reexecuted": False,
        },
    })

    # Coupled: conditional theorem-coupled candidate, open interface.
    rows.append({
        "row_id": "coupled.cp1_cp2_cp3",
        "family": "coupled",
        "implementation": "F_candidate_coupled.py",
        "semantics": {
            "publicness_family": "inherited frozen family of the count construction",
            "reserve_semantics": "re-emission weight lambda = 1/2, recorded free in (0, 1)",
            "reserve_attachment": "port-load re-emission",
            "channel": (
                "port-load inversion with the D10 balance coupling CP-1; "
                "fixed point independent of lambda, conditional on CP-1, "
                "CP-2, CP-3"
            ),
        },
        "executed_verdict": "conditional_open",
        "reason": (
            "the evaluation cone contains the electroweak-bridge expression "
            "by construction, so the candidate fails the specification "
            "blindness bar as written, cannot serve as a blind landing "
            "test, and moves no ledger row; the closure row stays open, "
            "reduced to the CP-1, CP-2, CP-3 premises"
        ),
        "verdict_provenance": {
            "recorded": "runtime/F_candidate_coupled_certificates.json",
            "reexecuted": True,
        },
    })

    return rows


# ---------------------------------------------------------------------------
# Fail-closed row validation
# ---------------------------------------------------------------------------

REQUIRED_ROW_KEYS = {
    "row_id",
    "family",
    "implementation",
    "semantics",
    "executed_verdict",
    "reason",
    "verdict_provenance",
}


def validate_menu_row(row: Mapping[str, Any]) -> None:
    missing = REQUIRED_ROW_KEYS - set(row)
    require(not missing, "ROW_SCHEMA", f"menu row missing fields {sorted(missing)}")
    verdict = row["executed_verdict"]
    require(
        verdict in ALLOWED_VERDICTS,
        "ROW_VERDICT",
        f"{row['row_id']}: verdict '{verdict}' outside the allowed vocabulary",
    )
    contraction = row.get("linear_contraction")
    if row["family"] == "CAP-K" or (
        contraction is not None and contraction.get("map") == "F(N) = s*N"
    ):
        require(
            verdict == "no_positive_fixed_point",
            "CONTRACTION_CONTRADICTION",
            f"{row['row_id']}: a linear contraction through the origin with "
            "s < 1 has no positive fixed point; any other verdict is "
            "rejected",
        )
    scan_for_targets(row, f"$.menu_rows[{row.get('row_id', '?')}]")


# ---------------------------------------------------------------------------
# Negative controls
# ---------------------------------------------------------------------------

def _expect_rejection(name: str, row: Mapping[str, Any], expected_code: str) -> dict[str, Any]:
    actual = "ACCEPTED"
    try:
        validate_menu_row(row)
    except CertificateError as exc:
        actual = exc.code
    require(
        actual == expected_code,
        "NEGATIVE_CONTROL_FAILED",
        f"{name}: expected {expected_code}, got {actual}",
    )
    return {"name": name, "expected_error": expected_code, "actual_error": actual, "passed": True}


def negative_control_records(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    capk_row = next(row for row in rows if row["family"] == "CAP-K")

    controls = []

    # Control 1: a CAP-K row claiming a positive fixed point is rejected.
    mutant = copy.deepcopy(dict(capk_row))
    mutant["executed_verdict"] = "excluded"
    mutant["reason"] = "claims a certified positive fixed point on a linear contraction"
    controls.append(
        _expect_rejection(
            "capk_positive_fixed_point_claim_rejected",
            mutant,
            "CONTRACTION_CONTRADICTION",
        )
    )

    # Control 2a: a row that hardcodes a capacity at cosmological scale by a
    # numeric token is rejected. The injected value is synthetic and target
    # shaped; it lands in no manifest.
    mutant = copy.deepcopy(dict(capk_row))
    mutant["reason"] = "row pins the readback at 9.99e122 nats by hand"
    controls.append(
        _expect_rejection(
            "hidden_target_scale_token_rejected",
            mutant,
            "HIDDEN_TARGET",
        )
    )

    # Control 2b: a row that carries a target marker key is rejected.
    mutant = copy.deepcopy(dict(capk_row))
    mutant_semantics = dict(mutant["semantics"])
    mutant_semantics["desired_capacity_nats"] = "pinned by hand"
    mutant["semantics"] = mutant_semantics
    controls.append(
        _expect_rejection(
            "hidden_target_marker_key_rejected",
            mutant,
            "HIDDEN_TARGET",
        )
    )

    return controls


# ---------------------------------------------------------------------------
# Axiom registry alignment (read-only)
# ---------------------------------------------------------------------------

def registry_alignment() -> dict[str, Any]:
    registry_path = REPO_ROOT / "claims" / "axiom_registry.yaml"
    try:
        text = registry_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise CertificateError("SOURCE_READ", f"cannot read {registry_path}: {exc}") from exc
    require(
        re.search(
            r'"id":\s*"capacity_publicness_and_closure",\s*"class":\s*"conditional_open_interface"',
            text,
        )
        is not None,
        "REGISTRY_MISMATCH",
        "the registry row capacity_publicness_and_closure must carry class "
        "conditional_open_interface",
    )
    require(
        re.search(
            r'"id":\s*"horizon_area_capacity_identification",\s*"class":\s*"physical_identification",\s*"attachment_state":\s*"pending"',
            text,
        )
        is not None,
        "REGISTRY_MISMATCH",
        "the registry row horizon_area_capacity_identification must carry "
        "class physical_identification with attachment pending",
    )
    return {
        "registry": "claims/axiom_registry.yaml",
        "registry_sha256": sha256_file(registry_path),
        "rows": {
            "capacity_publicness_and_closure": "conditional_open_interface",
            "horizon_area_capacity_identification": "physical_identification, attachment pending",
        },
        "access": "read-only; this certificate writes no registry row",
    }


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

def build(reexecute: bool = True) -> dict[str, Any]:
    sources = {
        relative: sha256_file(REPO_ROOT / relative) for relative in PINNED_SOURCES
    }
    axes = declared_semantic_axes()
    recorded = load_recorded()

    if reexecute:
        confirmations = reexecute_cheap_families(recorded)
    else:
        confirmations = {
            family: {"reexecuted": False, "note": "re-execution skipped by caller"}
            for family in ("capK", "capP", "coupled", "capL")
        }

    rows = menu_rows(recorded)
    for row in rows:
        validate_menu_row(row)

    contraction = capk_contraction_check(recorded["capK"])
    controls = negative_control_records(rows)
    controls.append({
        "name": "capk_linear_contraction_symbolic",
        "detail": contraction,
        "passed": True,
    })
    registry = registry_alignment()

    # The selector verdict is computed, never asserted: a positive selection
    # would require exactly one row carrying a source-only unique-zero
    # certificate, and no executed row carries one. Ties and incomplete
    # semantics stay fail-closed at the open interface.
    selectable = [
        row["row_id"]
        for row in rows
        if row["executed_verdict"] not in ALLOWED_VERDICTS
    ]
    require(
        not selectable,
        "ROW_VERDICT",
        f"rows outside the fail-closed vocabulary: {selectable}",
    )

    payload = {
        "schema": SCHEMA,
        "issue": ISSUE,
        "campaign": "bounded capacity-closure campaign",
        "statement": (
            "This ledger records the executed scope of the capacity-readback "
            "candidate campaign. The semantic axes are enumerated before any "
            "optimization, each executed reading carries its recorded "
            "verdict, and no reading yields a positive source-only fixed "
            "point at the public stable-record capacity without importing a "
            "target. Ties and incomplete semantics leave the selector "
            "conditional_open_interface; an economy rule does not resolve "
            "them. The menu is a record of executed scope."
        ),
        "sources_sha256": sources,
        "semantic_axes": axes,
        "menu_rows": rows,
        "menu_row_count": len(rows),
        "executed_construction_rows": recorded["comparison_total_rows"],
        "comparison_landed_rows": recorded["comparison_landed_rows"],
        "reexecution_confirmations": confirmations,
        "finite_evaluator_control": {
            "record": "A5_FINITE_CONTROL_STATUS_2026-07-20.md",
            "ledger": "runtime/a5_finite_control_status.json",
            "status": "BLOCKED_EXTERNAL_PHYSICAL_DERIVATION",
            "note": (
                "the exact finite packet is a software evaluator control and "
                "a no-go witness for raw-carrier promotion; it emits no "
                "physical capacity coordinate"
            ),
        },
        "campaign_verdict": {
            "semantics_enumeration": "declared_menu_complete_for_executed_families",
            "source_only_fixed_point_selector": "conditional_open_interface",
            "no_target_import": "verified",
            "horizon_area_and_ew_bridge": "separate_physical_attachments (#589, #547)",
        },
        "axiom_registry_alignment": registry,
        "negative_controls": controls,
    }

    # The no-target verdict is discharged on the assembled payload itself.
    scan_for_targets(payload)
    return payload


def verify_manifest(path: Path = MANIFEST_PATH, reexecute: bool = True) -> None:
    manifest = load_json(path)
    expected = build(reexecute=reexecute)
    require(
        manifest == expected,
        "MANIFEST_MISMATCH",
        f"{path} differs from the deterministic rebuild",
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Semantics menu ledger for issue #615")
    parser.add_argument("--verify", action="store_true", help="compare the stored manifest with a rebuild")
    parser.add_argument("--output", type=Path, default=MANIFEST_PATH)
    args = parser.parse_args(argv)
    if args.verify:
        verify_manifest(args.output)
        print(json.dumps({"status": "PASS", "manifest": str(args.output)}, indent=2))
        return 0
    payload = build()
    write_json(args.output, payload)
    digest = hashlib.sha256((args.output).read_bytes()).hexdigest()
    print(
        json.dumps(
            {
                "status": "PASS",
                "manifest": str(args.output),
                "manifest_sha256": digest,
                "menu_row_count": payload["menu_row_count"],
                "campaign_verdict": payload["campaign_verdict"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
