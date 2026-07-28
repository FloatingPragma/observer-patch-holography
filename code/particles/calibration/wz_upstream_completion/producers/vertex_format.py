#!/usr/bin/env python3
"""Canonical vertex-record format for Workstream C.

Format only: this module defines the record container, the canonical
monomial coefficient encoding, and the vertex hash.  It contains no rule
content, no field list, and no coupling assignment, so both rule engines
may import it without sharing derivation logic.

A vertex record is:

    fields   sorted tuple of component-field labels (the multiset),
    coefficient   exact monomial: a Fraction prefactor and a sorted
                  tuple of (symbol, integer power) pairs,
    structure     a named tensor/momentum basis label from the frozen
                  vocabulary below.

The canonical hash is the SHA-256 of the canonical JSON of the record.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from typing import Any, Sequence

STRUCTURE_VOCABULARY = (
    "lorentz_metric_pair",          # g_{mu nu} type contraction
    "yang_mills_three_point",       # f^{abc} (p1-p2)^rho g^{mu nu} + cyclic
    "yang_mills_four_point",        # f f (gg - gg) double structure
    "scalar_gauge_gauge",           # g_{mu nu} S V V seagull or v-induced mass
    "scalar_scalar_gauge",          # (p1 - p2)^mu derivative coupling
    "scalar_scalar_gauge_gauge",    # g_{mu nu} seagull
    "scalar_potential",             # pure scalar, no tensor structure
    "fermion_vector_current",       # gamma^mu (chiral projector) current
    "fermion_scalar_yukawa",        # chiral scalar coupling
    "ghost_gauge_derivative",       # p^mu ghost-gauge coupling
    "fermion_bilinear_mass",        # v-induced fermion bilinear
    "scalar_bilinear_mass",         # v-induced scalar bilinear
    "vector_bilinear_mass",         # v-induced vector bilinear
)


def monomial(prefactor: Fraction | int, *powers: tuple[str, int]) -> dict[str, Any]:
    ordered = tuple(sorted((s, int(p)) for s, p in powers if p != 0))
    return {
        "prefactor": str(Fraction(prefactor)),
        "powers": [list(pair) for pair in ordered],
    }


def record(fields: Sequence[str], coefficient: dict[str, Any], structure: str) -> dict[str, Any]:
    if structure not in STRUCTURE_VOCABULARY:
        raise ValueError(f"unknown structure label {structure}")
    return {
        "fields": sorted(fields),
        "coefficient": coefficient,
        "structure": structure,
    }


def vertex_hash(entry: dict[str, Any]) -> str:
    payload = json.dumps(entry, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def table_digest(entries: Sequence[dict[str, Any]]) -> str:
    hashes = sorted(vertex_hash(entry) for entry in entries)
    payload = json.dumps(hashes, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()
