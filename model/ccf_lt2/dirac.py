"""Signed-band charged/neutral Dirac terms on the common native CAR carrier.

The polynomial has 96 native modes per vertex. No positive-band projection,
lookup arrival law, or separately particle-number-conserving replacement is used.
"""
from __future__ import annotations

from collections import defaultdict
import argparse
import hashlib
import json
from pathlib import Path
import sys

import sympy as sy

try:
    from . import core
except ImportError:
    from . import core


I2 = sy.eye(2)
PAULI = (sy.Matrix([[0, 1], [1, 0]]),
         sy.Matrix([[0, -sy.I], [sy.I, 0]]), sy.diag(1, -1))
BETA = sy.diag(1, 1, -1, -1)
ALPHA = tuple(sy.kronecker_product(PAULI[0], p) for p in PAULI)
KINDS = (("charged", 11, 9, -3), ("neutral", 10, 15, 0))


def _clean_g(poly):
    return {k: c for k, v in poly.items() if (c := sy.expand(v)) != 0}


def _gadd(*polys):
    result = defaultdict(lambda: sy.Integer(0))
    for poly in polys:
        for key, coefficient in poly.items():
            result[key] += coefficient
    return _clean_g(result)


def _gscale(poly, coefficient):
    return _clean_g({key: sy.sympify(coefficient)*value for key, value in poly.items()})


def _dress(poly, shift):
    return {(word, tuple(shift)): value for word, value in poly.items()}


def _gadjoint(poly):
    result = {}
    for (word, shift), value in poly.items():
        result = _gadd(result, _dress(core.adjoint({word: value}), tuple(-x for x in shift)))
    return result


def pair_map(nvertices):
    """Exact a -> c substitution, including both neutral and charged pairs."""
    mapping = {}
    for v in range(nvertices):
        for _, left, right, _ in KINDS:
            for f in range(3):
                for s in range(2):
                    e = core.index(v, s, left, f)
                    p = core.index(v, 1-s, right, f)
                    sigma = 1-2*s
                    ae = core.scale(core.add(core.letter(e), core.scale(core.letter(p, True), sigma)), 1/sy.sqrt(2))
                    ap = core.scale(core.add(core.letter(p), core.scale(core.letter(e, True), -sigma)), 1/sy.sqrt(2))
                    mapping[(e, False)] = ae
                    mapping[(p, False)] = ap
                    mapping[(e, True)] = core.adjoint(ae)
                    mapping[(p, True)] = core.adjoint(ap)
    return mapping


def _substitute(poly, mapping):
    answer = {}
    for word, coefficient in poly.items():
        expanded = core.mono(coefficient=coefficient)
        for letter in word:
            expanded = core.multiply(expanded, mapping.get(letter, core.letter(*letter)))
        answer = core.add(answer, expanded)
    return answer


def to_native(poly, nvertices):
    mapping = pair_map(nvertices)
    answer = {}
    for (word, shift), coefficient in poly.items():
        answer = _gadd(answer, _dress(_substitute({word: coefficient}, mapping), shift))
    return answer


def _spinor(v, left, right, family):
    return [core.letter(core.index(v, s, left, family)) for s in range(2)] + [
        core.scale(core.letter(core.index(v, 1-s, right, family), True), 1-2*s)
        for s in range(2)]


def _bilinear(left, matrix, right):
    answer = {}
    for i in range(4):
        for j in range(4):
            if matrix[i, j] != 0:
                answer = core.add(answer, core.scale(core.multiply(core.adjoint(left[i]), right[j]), matrix[i, j]))
    return answer


def canonical_dirac_terms(nvertices, edges, directions, spacings=None, c=1, wilson_r=1):
    """Same full polynomial before the exact inverse canonical CAR map.

    Coordinates label canonical modes in this function only. ``dirac_terms``
    returns native coordinates for the common-source assembler.
    """
    if len(edges) != len(directions):
        raise ValueError("one direction is required for each edge")
    if spacings is None:
        spacings = [sy.Integer(1)]*len(edges)
    if len(spacings) != len(edges):
        raise ValueError("one spacing is required for each edge")
    c, wilson_r = sy.sympify(c), sy.sympify(wilson_r)
    if c.is_positive is not True or wilson_r.is_nonnegative is not True:
        raise ValueError("positive speed and nonnegative Wilson parameter required")
    answer = {}
    zero = (0,)*len(edges)
    for k, ((u, v), direction, raw_a) in enumerate(zip(edges, directions, spacings)):
        if u == v or not (0 <= u < nvertices and 0 <= v < nvertices):
            raise ValueError("edge endpoints must be distinct vertices in the carrier")
        n = tuple(map(sy.sympify, direction))
        if len(n) != 3 or sy.simplify(sum(x*x for x in n)-1) != 0:
            raise ValueError("edge direction must be an exact unit three-vector")
        a = sy.sympify(raw_a)
        if a.is_positive is not True:
            raise ValueError("spacing must be positive")
        alpha = sum((n[i]*ALPHA[i] for i in range(3)), sy.zeros(4))
        for _, left, right, q in KINDS:
            shift = tuple(-q if j == k else 0 for j in range(len(edges)))
            for f in range(3):
                x, y = _spinor(u, left, right, f), _spinor(v, left, right, f)
                derivative = _dress(core.scale(_bilinear(y, alpha, x), sy.I*c/(2*a)), shift)
                offdiagonal_wilson = _dress(core.scale(_bilinear(y, BETA, x), -wilson_r/(2*a)), shift)
                local_wilson = _dress(core.scale(core.add(_bilinear(x, BETA, x), _bilinear(y, BETA, y)), wilson_r/(2*a)), zero)
                answer = _gadd(answer, derivative, _gadjoint(derivative),
                               offdiagonal_wilson, _gadjoint(offdiagonal_wilson), local_wilson)
    return answer


def dirac_terms(nvertices, edges, directions, spacings=None, c=1, wilson_r=1):
    """Return the full native GaugedPoly for all lepton spins/families/signs."""
    return to_native(canonical_dirac_terms(nvertices, edges, directions, spacings, c, wilson_r), nvertices)


def clifford_checks():
    for i in range(3):
        assert ALPHA[i].H == ALPHA[i]
        assert ALPHA[i]*BETA+BETA*ALPHA[i] == sy.zeros(4)
        for j in range(3):
            assert ALPHA[i]*ALPHA[j]+ALPHA[j]*ALPHA[i] == 2*int(i == j)*sy.eye(4)
    assert BETA*BETA == sy.eye(4)
    sigmas = tuple(sy.kronecker_product(I2, p) for p in PAULI)
    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        assert ALPHA[i]*ALPHA[j] == sy.I*sigmas[k]
    assert sigmas[2]*ALPHA[0]-ALPHA[0]*sigmas[2] != sy.zeros(4)
    # A compatible mass is a hypothesis of this exact symbolic dispersion.
    p = sy.symbols("p0:3", real=True)
    m, a, c, r = sy.symbols("m a c r", positive=True)
    w = m+r*sum(1-sy.cos(a*x) for x in p)/a
    h = sum((c*ALPHA[i]*sy.sin(a*p[i])/a for i in range(3)), sy.zeros(4))+BETA*w
    target = (c*c*sum(sy.sin(a*x)**2 for x in p)/a**2+w*w)*sy.eye(4)
    assert (h*h-target).applyfunc(sy.expand) == sy.zeros(4)
    return {"clifford_and_beta_exact": True,
            "same_principal_speed_charged_and_neutral": "c",
            "compatible_mass_symbol_square_exact": True,
            "wilson_regulator_retained": True,
            "spin_sensitive_directional_derivative": True,
            "pauli_commutator_coefficient_exact": True,
            "actual_rounded_parent_neutral_mass_compatibility_not_assumed": True}


def charge_and_pair_checks(poly, edges):
    """Native local Gauss identity and actual canonical pair-current witness."""
    native_charge = {9: 3, 11: -3, 10: 0, 15: 0}
    for (word, shift), coefficient in poly.items():
        changes = [0]*(max(max(e) for e in edges)+1)
        for i, creation in word:
            v, rem = divmod(i, 96)
            coordinate = (rem % 48)//3
            changes[v] += native_charge[coordinate]*(1 if creation else -1)
        boundary = [0]*len(changes)
        for (u, v), n in zip(edges, shift):
            boundary[u] += n
            boundary[v] -= n
        assert changes == boundary
    assert _gadjoint(poly) == poly
    canonical = canonical_dirac_terms(2, [(0, 1)], [(0, 0, 1)])
    pair = {(w, sh): v for (w, sh), v in canonical.items()
            if len(w) == 2 and w[0][1] == w[1][1] and sh != (0,)}
    assert pair
    # [N_e+N_p,H] has signed degree +/-2 on these monomials.
    number_commutator = {(w, sh): (sum(1 if cre else -1 for _, cre in w))*v
                         for (w, sh), v in canonical.items()
                         if sum(1 if cre else -1 for _, cre in w)}
    assert number_commutator
    return {"native_hermiticity_and_local_gauss_exact": True,
            "charged_pair_current_monomials": len(pair),
            "separate_particle_number_commutator_nonzero": True,
            "canonical_total_occupation_commutator_terms": len(number_commutator)}


def _act_word(mask, word):
    sign = 1
    for i, creation in reversed(word):
        occupied = (mask >> i) & 1
        if occupied == creation:
            return None
        sign *= -1 if (mask & ((1 << i)-1)).bit_count() % 2 else 1
        mask ^= 1 << i
    return mask, sign


def edge_reduction(poly, cap=6):
    """Exact Pauli-blocked sector of the *full generated lepton edge term*.

    All inactive charged up modes are occupied, down antimodes are occupied;
    all neutral particle modes are occupied. The active charged up family0
    has its four occupations retained. This reduces the complete added
    term on a z edge. Parent source preservation is a separate assembler
    check, not inferred here.
    """
    active = [core.index(v, s, coord, 0) for v in range(2) for s, coord in ((0, 11), (1, 9))]
    frozen = 0
    for v in range(2):
        for f in range(3):
            for s in range(2):
                neutral = core.index(v, s, 10, f)
                frozen |= 1 << neutral
            for s, coordinate in ((0, 11), (0, 9)):
                i = core.index(v, s, coordinate, f)
                if i not in active:
                    frozen |= 1 << i
    active_mask = sum(1 << i for i in active)
    bases = []
    lookup = {}
    for bits in range(16):
        mask = frozen | sum(((bits >> j) & 1) << i for j, i in enumerate(active))
        for flux in range(-cap, cap+1):
            lookup[(mask, flux)] = len(bases)
            bases.append((mask, flux, bits))
    h = sy.zeros(len(bases))
    for col, (mask, flux, bits) in enumerate(bases):
        leakage = defaultdict(lambda: sy.Integer(0))
        for (word, shift), coefficient in poly.items():
            result = _act_word(mask, word)
            if result is None:
                continue
            target, sign = result
            next_flux = flux+shift[0]
            if not -cap <= next_flux <= cap:
                continue
            key = (target, next_flux)
            if key in lookup:
                h[lookup[key], col] += sign*coefficient
            else:
                leakage[key] += sign*coefficient
        assert all(sy.expand(v) == 0 for v in leakage.values()), "frozen spectator leakage"
    assert h == h.H
    charges = []
    for mask, flux, bits in bases:
        q = [0, 0]
        for v in range(2):
            for s in range(2):
                for f in range(3):
                    q[v] += 3*((mask >> core.index(v, s, 9, f)) & 1)
                    q[v] -= 3*((mask >> core.index(v, s, 11, f)) & 1)
        charges.append(q)
    r = -6
    gauss = [(flux-q[0]-r, -flux-q[1]) for (_, flux, _), q in zip(bases, charges)]
    for i, j in h.todok():
        assert gauss[i] == gauss[j]
    return h, bases, lookup, active, frozen, gauss


def finite_pair_witness():
    import numpy as np
    from scipy.linalg import expm
    edges, directions = [(0, 1)], [(0, 0, 1)]
    canonical = canonical_dirac_terms(2, edges, directions)
    kinetic, bases, lookup, active, frozen, gauss = edge_reduction(canonical)
    dimension = len(bases)
    mass_and_electric = sy.diag(*[sy.Integer(bits.bit_count())+sy.Rational(flux*flux, 16)
                                 for _, flux, bits in bases])
    h = kinetic+mass_and_electric
    pairmask = frozen | (1 << active[0]) | (1 << active[3])
    initial = lookup[(pairmask, -6)]
    vacuum = lookup[(frozen, -3)]
    assert gauss[initial] == (0, 0) and gauss[vacuum] == (0, 0)
    assert h[vacuum, initial] != 0
    control_poly = {(w, sh): v for (w, sh), v in canonical.items()
                    if sum(1 if cre else -1 for _, cre in w) == 0}
    control, *_ = edge_reduction(control_poly)
    control += mass_and_electric
    assert control[vacuum, initial] == 0
    vector = np.zeros(dimension, complex)
    vector[initial] = 1
    numeric = np.array(h.tolist(), dtype=complex)
    time = 0.25
    output = expm(-1j*time*numeric) @ vector
    unchanged_number = expm(-1j*time*np.array(control.tolist(), dtype=complex)) @ vector
    vacuum_probability = float(sum(abs(output[i])**2 for i, (_, _, bits) in enumerate(bases) if bits == 0))
    control_probability = float(sum(abs(unchanged_number[i])**2 for i, (_, _, bits) in enumerate(bases) if bits == 0))
    mismatch = float(sum(abs(output[i])**2 for i, g in enumerate(gauss) if g != (0, 0)))
    assert vacuum_probability > 1e-8
    assert control_probability < 1e-24
    assert mismatch < 1e-24
    assert abs(np.vdot(output, output).real-1) < 1e-11
    energy_error = float(abs(np.vdot(output, numeric @ output)-numeric[initial, initial]))
    assert energy_error < 1e-10
    work = float(np.dot(abs(output)**2, np.diag(numeric).real)-np.vdot(output, numeric @ output).real)
    return {"full_active_fock_rotor_dimension": dimension,
            "cap": 6, "time": time, "all_outcomes_retained": True,
            "spectators_reduce_full_generated_kinetic_exactly": True,
            "parent_mass_source_reduction_not_claimed_by_this_component": True,
            "pair_to_vacuum_matrix_element_exact": str(h[vacuum, initial]),
            "canonical_active_vacuum_probability": vacuum_probability,
            "number_conserving_control_vacuum_probability": control_probability,
            "gauss_mismatch_probability": mismatch,
            "source_energy_residual": energy_error,
            "terminal_occupation_flux_meter_work": work,
            "outcome_probabilities": [{"active_mask": bits, "flux": flux, "probability": float(abs(output[i])**2)}
                                      for i, (_, flux, bits) in enumerate(bases)]}


def verify():
    checks = clifford_checks()
    for name, nv, edges, directions in (
            ("edge", 2, [(0, 1)], [(0, 0, 1)]),
            ("triangle", 3, [(0, 1), (1, 2), (0, 2)], [(1, 0, 0), (0, 1, 0), (sy.sqrt(2)/2, sy.sqrt(2)/2, 0)])):
        polynomial = dirac_terms(nv, edges, directions)
        checks[name] = charge_and_pair_checks(polynomial, edges)
        checks[name]["native_polynomial_monomials"] = len(polynomial)
        checks[name]["all_native_mode_labels_retained_by_common_source"] = 96*nv
    # Local four-state CAR map: symbolic identities on one charged pair.
    mapping = pair_map(1)
    e, p = core.index(0, 0, 11, 0), core.index(0, 1, 9, 0)
    for i in (e, p):
        a = mapping[(i, False)]
        assert core.add(core.multiply(a, core.adjoint(a)), core.multiply(core.adjoint(a), a)) == core.mono()
    assert _substitute(core.add(core.scale(core.number(e), -3), core.scale(core.number(p), 3)), mapping) == core.add(core.scale(core.number(e), -3), core.scale(core.number(p), 3))
    checks["pair_map_car_and_charge_exact"] = True
    checks["pair_witness"] = finite_pair_witness()
    return checks

