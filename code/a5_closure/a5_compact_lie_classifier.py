#!/usr/bin/env python3
"""Finite dimension certificate for the A5 inner-action compact-closure theorem.

External mathematical input (not proved by this script): compact real Lie
algebras are reductive, and the compact simple Lie algebras of dimension at
most 12 have dimensions 3 (su2), 8 (su3), and 10 (so5=sp2).
"""
from __future__ import annotations

import json
from functools import lru_cache

SIMPLE_DIMS = {3: "su(2)", 8: "su(3)", 10: "so(5)=sp(2)"}


def partitions(total: int, allowed: tuple[int, ...] = (3, 8, 10)) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []
    def rec(rem: int, start: int, acc: list[int]) -> None:
        if rem == 0:
            out.append(tuple(acc))
            return
        for i in range(start, len(allowed)):
            d = allowed[i]
            if d <= rem:
                rec(rem - d, i, acc + [d])
    rec(total, 0, [])
    return out


def payload() -> dict:
    p12 = partitions(12)
    p11 = partitions(11)
    assert p12 == [(3, 3, 3, 3)]
    assert p11 == [(3, 8)]
    return {
        "schema": "A5 inner-action compact closure dimension certificate v1",
        "hypotheses": [
            "g is a compact real Lie algebra of dimension 12",
            "as an A5-module, g = 1 + 3 + 3prime + 5, so dim(g^A5)=1",
            "A5 acts by inner automorphisms",
            "classification of compact simple Lie algebras through dimension 12",
        ],
        "dimension_partitions": {
            "semisimple_dimension_12": [list(x) for x in p12],
            "semisimple_dimension_11": [list(x) for x in p11],
        },
        "proof_spine": [
            "Inner automorphisms fix the center pointwise, hence dim Z(g) <= 1.",
            "If Z(g)=0, g is semisimple of dimension 12, hence su(2)^4.",
            "For an inner A5 action on each su(2) factor, its fixed space has dimension 0 (icosahedral action) or 3 (trivial action); four factors give a fixed dimension divisible by 3, contradicting dim(g^A5)=1.",
            "Therefore dim Z(g)=1 and dim [g,g]=11.",
            "The only compact semisimple dimension partition of 11 is 8+3, hence [g,g]=su(3)+su(2).",
            "Thus g=su(3)+su(2)+u(1). Under the hypotheses, the 8-dimensional ideal is 5 plus one triplet and the 3-dimensional ideal is the other triplet, up to the outer automorphism of A5.",
        ],
        "conclusion": "g is isomorphic to su(3) + su(2) + u(1), conditional on the hypotheses",
        "what_is_not_proved": [
            "the physical port-to-current map",
            "innerness of the screen A5 action from bare OPH repair axioms",
            "the kinetic form, locality, or normalization",
            "trace-balanced physical integration or axis-center descent",
            "a physical matter selection, hypercharge attachment, or global gauge quotient",
        ],
        "counterbranches_if_innerness_is_dropped_but_group_integrality_is_retained": [
            "u(1)^12 with a nontrivial A5 lattice action",
            "su(2)+su(2)+u(1)^6 with center module 1+5",
        ],
        "additional_real_module_counterbranch_only_if_group_level_integral_central_lattice_constraints_are_also_dropped": [
            "su(3)+u(1)^4 with a nonintegral real A5 center module 1+3",
        ],
    }


if __name__ == "__main__":
    print(json.dumps(payload(), indent=2))
