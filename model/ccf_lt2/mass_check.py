"""Exact actual-parent local mass compatibility with CCF-2 Dirac kinetics."""
from __future__ import annotations
from itertools import product
import argparse
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import sympy as sy

try:
    from . import core, dirac
except ImportError:
    from . import core
    from . import dirac


def _inverse_pair_map():
    """Native c expressed in exact canonical a; no PMNS matrix enters."""
    mapping = {}
    for _, left, right, _ in dirac.KINDS:
        for f, s in product(range(3), range(2)):
            e, p = core.index(0, s, left, f), core.index(0, 1-s, right, f)
            sigma = 1-2*s
            mapping[(e, False)] = core.scale(core.add(core.letter(e), core.scale(core.letter(p, True), -sigma)), 1/sy.sqrt(2))
            mapping[(p, False)] = core.scale(core.add(core.letter(p), core.scale(core.letter(e, True), sigma)), 1/sy.sqrt(2))
            mapping[(e, True)] = core.adjoint(mapping[(e, False)])
            mapping[(p, True)] = core.adjoint(mapping[(p, False)])
    return mapping


def _local_mass(local_poly):
    allowed = [core.index(0, s, c, f) for s, c, f in product(range(2), (9, 10, 11, 15), range(3))]
    mass = core.vacuum_restrict(local_poly, allowed)
    assert all(len(w) in (0, 2) for w in mass), "local lepton restriction is not purely quadratic"
    assert core.adjoint(mass) == mass, "local lepton mass is not Hermitian"
    charge = {9: 3, 11: -3, 10: 0, 15: 0}
    for word in mass:
        assert sum(charge[(i % 48)//3]*(1 if creation else -1) for i, creation in word) == 0
    return mass


def _evaluate(local_poly):
    mass = _local_mass(local_poly)
    canonical = core.substitute(mass, _inverse_pair_map())
    kinds = {}
    for name, left, right, q in dirac.KINDS:
        spinor = sum([dirac._spinor(0, left, right, f) for f in range(3)], [])
        basis = {next(iter(p)): (i, next(iter(p.values()))) for i, p in enumerate(spinor)}
        assert len(basis) == 12
        # All Nambu field entries have the same annihilation-charge change.
        coordinate_q = {9: 3, 11: -3, 10: 0, 15: 0}
        for word in basis:
            (i, creation), = word
            assert coordinate_q[(i % 48)//3]*(1 if creation else -1) == -q
        h = sy.zeros(12)
        for i, p in enumerate(spinor):
            commutator = core.add(core.multiply(p, canonical), core.scale(core.multiply(canonical, p), -1))
            for word, value in commutator.items():
                assert word in basis, "mass commutator leaves the complete Nambu field"
                j, scale = basis[word]
                h[i, j] += value/scale
        h = h.applyfunc(sy.expand)
        assert h == h.H, "Nambu mass is not Hermitian"
        for alpha in dirac.ALPHA:
            full_alpha = sy.kronecker_product(sy.eye(3), alpha)
            residual = (h*full_alpha+full_alpha*h).applyfunc(sy.expand)
            assert residual == sy.zeros(12), "actual mass does not anticommute with Dirac principal part"
        payload = json.dumps([[str(h[i, j]) for j in range(12)] for i in range(12)], separators=(",", ":"))
        kinds[name] = {"nambu_dimension": 12, "complete_nambu_closure": True,
                       "charge_pairing_exact": True, "hermitian": True,
                       "all_three_mass_dirac_anticommutators_exact_zero": True,
                       "mass_matrix_sha256": hashlib.sha256(payload.encode()).hexdigest(),
                       "nonzero_matrix_entries": len(h.todok())}
    return {"local_lepton_mass_monomials": len(mass),
            "exact_canonical_mass_monomials": len(canonical),
            "mass_coefficients_unchanged": True, "kinds": kinds}


def parent_mass_compatibility(local_poly):
    """Audit an imported matter local block, including negative control."""
    result = _evaluate(local_poly)
    # Native left-chiral density preserves charge, spin degeneracy and locality.
    # It changes the mass/kinetic relation rather than corrupting the instrument.
    control = core.add(local_poly, *[core.number(core.index(0, s, 11, 0)) for s in range(2)])
    _local_mass(control)  # Hermiticity and charge gates must still pass.
    try:
        _evaluate(control)
    except AssertionError as exc:
        assert str(exc) == "actual mass does not anticommute with Dirac principal part", str(exc)
    else:
        raise AssertionError("charge-preserving chiral-density control was not rejected")
    result["charge_preserving_spin_summed_density_control_rejected"] = True
    result["scope"] = "Actual local matter mass and Clifford principal part; not full interacting dispersion."
    return result


def verify():
    try:
        from . import source
    except ImportError:
        from . import source
    parent, _, _ = source.parent_source()
    local, _ = source.parent_parts(parent)
    return parent_mass_compatibility(local)

