#!/usr/bin/env python3
"""Reduced-domain checks for the global uniqueness certificate.

Runs ``global_uniqueness_certificate.sweep_mode`` at low cutoffs and digits on
a narrow alpha domain around the source fixed point (fast) and asserts the
machinery:

- every piece certifies (no exceptional set) and the conclusion synthesizes
  to at-most-one on the domain;
- the emitted block carries the certificate schema keys;
- the run-analysis synthesis handles mixed and partial verdict patterns;
- the construction is deterministic.
"""

from __future__ import annotations

import global_uniqueness_certificate as guc
import interval_contraction_certificate as icc


REDUCED = dict(
    mp_dps=30,
    iv_dps=30,
    su2_cutoff=24,
    su3_cutoff=16,
    domain_lo="0.00728",
    domain_hi="0.00732",
    initial_pieces=4,
    depth_cap=6,
    eval_budget=200,
)

SCHEMA_KEYS = {
    "mode",
    "map_definition",
    "backend",
    "iv_dps",
    "point_dps",
    "su2_cutoff",
    "su3_cutoff",
    "edge_sum_tail_bounds",
    "subdivision",
    "verdict_counts",
    "worst_piece",
    "verdict_blocks",
    "exceptional_set",
    "exceptional_set_total_width",
    "conclusion",
}


def _sweep() -> dict:
    return guc.sweep_mode(icc.MODE_SOURCE, **REDUCED)


def test_reduced_domain_sweep_certifies() -> None:
    block = _sweep()
    assert SCHEMA_KEYS <= set(block.keys())

    counts = block["verdict_counts"]
    assert sum(counts.values()) == block["subdivision"]["total_pieces"]
    decisive = sum(counts.get(v, 0) for v in guc.DECISIVE_VERDICTS)
    assert decisive == block["subdivision"]["total_pieces"]
    assert block["exceptional_set"] == []

    conclusion = block["conclusion"]
    assert conclusion["at_most_one_on_domain"] is True
    assert conclusion["status"] == "certified"

    # Worst Lipschitz bound is recorded and below one on this all-contraction run.
    worst = block["worst_piece"]["max_sup_abs_gprime"]
    assert worst is not None
    assert float(worst) < 1.0

    # Tail bounds are declared included.
    assert block["edge_sum_tail_bounds"]["included"] is True


def test_sweep_is_deterministic() -> None:
    first = _sweep()
    second = _sweep()
    assert first == second


# ---------------------------------------------------------------------------
# Synthesis unit tests on fabricated piece patterns (no interval work).
# ---------------------------------------------------------------------------


def _piece(lo: float, hi: float, verdict: str) -> dict:
    return {"lo": lo, "hi": hi, "verdict": verdict}


def _no_signs(_x: float) -> None:
    return None


def test_synthesize_all_contraction() -> None:
    pieces = [
        _piece(0.0, 0.1, guc.VERDICT_CONTRACTION),
        _piece(0.1, 0.2, guc.VERDICT_CONTRACTION),
    ]
    out = guc.synthesize(pieces, _no_signs)
    assert out["at_most_one_on_domain"] is True
    assert out["status"] == "certified"


def test_synthesize_decreasing_family() -> None:
    pieces = [
        _piece(0.0, 0.1, guc.VERDICT_CONTRACTION),
        _piece(0.1, 0.2, guc.VERDICT_DECREASING),
    ]
    out = guc.synthesize(pieces, _no_signs)
    assert out["at_most_one_on_domain"] is True


def test_synthesize_partial_on_unresolved() -> None:
    pieces = [
        _piece(0.0, 0.1, guc.VERDICT_CONTRACTION),
        _piece(0.1, 0.2, guc.VERDICT_UNDECIDED),
    ]
    out = guc.synthesize(pieces, _no_signs)
    assert out["at_most_one_on_domain"] is False
    assert out["status"] == "partial"


def test_synthesize_single_candidate_run_with_root_free_gap() -> None:
    # Root-free gap splits the domain; only one monotone run can hold a root.
    pieces = [
        _piece(0.0, 0.1, guc.VERDICT_DECREASING),
        _piece(0.1, 0.2, guc.VERDICT_ROOT_FREE_NEG),
        _piece(0.2, 0.3, guc.VERDICT_DECREASING),
    ]

    def signs(x: float) -> str | None:
        # Residual positive at 0.0, negative from 0.1 on: root only in run one.
        return "+" if x == 0.0 else "-"

    out = guc.synthesize(pieces, signs)
    assert out["at_most_one_on_domain"] is True
    assert out["monotone_runs"] == 2
    assert out["candidate_runs"] == 1


def test_synthesize_two_candidate_runs_not_established() -> None:
    pieces = [
        _piece(0.0, 0.1, guc.VERDICT_DECREASING),
        _piece(0.1, 0.2, guc.VERDICT_ROOT_FREE_NEG),
        _piece(0.2, 0.3, guc.VERDICT_INCREASING),
    ]

    def signs(x: float) -> str | None:
        return {0.0: "+", 0.1: "-", 0.2: "-", 0.3: "+"}.get(x)

    out = guc.synthesize(pieces, signs)
    assert out["at_most_one_on_domain"] is False
    assert out["status"] == "not_established"
    assert out["candidate_runs"] == 2
