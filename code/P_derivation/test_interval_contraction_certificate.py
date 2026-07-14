#!/usr/bin/env python3
"""Reduced-precision checks for the interval contraction certificate.

Runs ``interval_contraction_certificate.certify_mode`` at low cutoffs and
digits (fast) and asserts the certificate invariants:

- the centered-form image ``K(I) = g(mid) + g'(I)(I - mid)`` is contained in
  the interior of the certified interval ``I``;
- the interval Lipschitz bound satisfies ``L < 1``;
- the interval-chain midpoint agrees with the repository Decimal chain
  (``paper_math.py``) to at least 30 significant digits;
- the emitted JSON carries the certificate schema keys;
- the construction is deterministic.
"""

from __future__ import annotations

from decimal import Decimal

import interval_contraction_certificate as icc


REDUCED = dict(
    mp_dps=35,
    iv_dps=35,
    su2_cutoff=24,
    su3_cutoff=16,
    half_width="0.000004",
    refine_passes=2,
)

SCHEMA_KEYS = {
    "mode",
    "map_definition",
    "backend",
    "iv_dps",
    "point_dps",
    "su2_cutoff",
    "su3_cutoff",
    "padding_policy",
    "edge_sum_tail_bounds",
    "uniqueness_interval_alpha",
    "banach",
    "refinement_passes",
    "certified_enclosure",
    "inner_root_certificates",
}

BANACH_KEYS = {
    "g_of_interval",
    "g_maps_interval_into_interior",
    "gprime_enclosure",
    "lipschitz_bound",
    "contraction",
    "existence",
    "uniqueness_in_interval",
}


def _certify(mode: str) -> dict:
    return icc.certify_mode(mode, **REDUCED)


def _assert_certificate_invariants(block: dict) -> None:
    assert SCHEMA_KEYS <= set(block.keys())
    assert BANACH_KEYS <= set(block["banach"].keys())

    interval = block["uniqueness_interval_alpha"]
    image = block["banach"]["g_of_interval"]
    lo_i, hi_i = Decimal(interval["lo"]), Decimal(interval["hi"])
    lo_k, hi_k = Decimal(image["lo"]), Decimal(image["hi"])
    # K(I) strictly inside I.
    assert lo_i < lo_k and hi_k < hi_i
    assert block["banach"]["g_maps_interval_into_interior"] is True

    # Contraction bound L < 1.
    lipschitz = Decimal(block["banach"]["lipschitz_bound"])
    assert Decimal(0) <= lipschitz < Decimal(1)
    assert block["banach"]["contraction"] is True
    assert block["banach"]["existence"] is True
    assert block["banach"]["uniqueness_in_interval"] is True

    # Certified enclosure inside the uniqueness interval, with finite width.
    enclosure = block["certified_enclosure"]["alpha"]
    lo_e, hi_e = Decimal(enclosure["lo"]), Decimal(enclosure["hi"])
    assert lo_i < lo_e <= hi_e < hi_i
    assert Decimal(enclosure["width"]) < Decimal(interval["half_width"])

    # Inner-root certificates carry the sign-definiteness verdicts.
    roots = block["inner_root_certificates"]
    assert roots["alpha_U_pixel_closure"]["endpoint_signs_verified"] is True
    assert roots["alpha_U_pixel_closure"]["R_u_sign_definite"] is True
    assert roots["m_Z_tree_closure"]["endpoint_signs_verified"] is True
    assert roots["m_Z_tree_closure"]["h_m_sign_definite"] is True

    # Tail bounds are declared included.
    assert block["edge_sum_tail_bounds"]["included"] is True


def test_source_mode_certificate() -> None:
    block = _certify(icc.MODE_SOURCE)
    _assert_certificate_invariants(block)
    alpha_inv = block["certified_enclosure"]["alpha_inv"]
    # The reduced-cutoff fixed point sits near the source value 136.99.
    assert Decimal("136.9") < Decimal(alpha_inv["lo"]) < Decimal("137.0")


def test_gauge_width_mode_certificate() -> None:
    block = _certify(icc.MODE_GAUGE_WIDTH)
    _assert_certificate_invariants(block)
    alpha_inv = block["certified_enclosure"]["alpha_inv"]
    # The reduced-cutoff mixed fixed point sits near 137.03.
    assert Decimal("137.0") < Decimal(alpha_inv["lo"]) < Decimal("137.1")


def test_midpoint_consistency_with_decimal_chain() -> None:
    """The mp midpoint chain matches paper_math's Decimal chain to >= 30 digits.

    The Decimal chain resolves its internal bisections to about
    ``2^-(precision+8)`` of the initial brackets, so precision 110 leaves a
    truncation floor near 1e-35; agreement at >= 30 digits confirms the two
    implementations encode the same map.
    """
    crosscheck = icc.decimal_crosscheck(
        mp_dps=60,
        su2_cutoff=6,
        su3_cutoff=4,
        decimal_precision=110,
        alpha_probe="0.0073",
    )
    assert crosscheck["agreement_digits"] >= 30


def test_certificate_is_deterministic() -> None:
    first = _certify(icc.MODE_SOURCE)
    second = _certify(icc.MODE_SOURCE)
    assert first == second
