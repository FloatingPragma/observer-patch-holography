"""Exact exterior calculation of the weak-epsilon odd composite.

The spinor and CAR structures are declared inputs. This is a portable
canonical-chart recomputation, not a continuum spin-statistics derivation.
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

from functools import lru_cache
from itertools import combinations, product
import json

import sympy as s


SOURCE_REFS = [
    "model/particle_primary_composite_spin/model.py",
    "model/particle_primary_composite_spin/PROOF.md",
    "model/native_quantum_gauss/algebra.py",
    "model/cubic_particle_channels/group.py",
]
PAULI = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))


def clean(a):
    return a.applyfunc(s.expand)


def basis_masks():
    # Four weak modes (spin, weak component), followed by two singlet spin modes.
    return tuple(sum(1 << j for j in pair)+(1 << k)
                 for pair in combinations(range(4), 2) for k in (4, 5))


def apply_word(mask, modes, *, fermionic=True):
    """Apply an ascending product of creators, rightmost factor first."""
    coefficient = 1
    for j in reversed(modes):
        if mask & (1 << j):
            return None, 0
        if fermionic:
            coefficient *= (-1)**((mask & ((1 << j)-1)).bit_count())
        mask |= 1 << j
    return mask, coefficient


def one_body_lift(h, masks):
    """Exterior derivative of h using literal Jordan-Wigner signs."""
    where = {m: i for i, m in enumerate(masks)}
    out = s.zeros(len(masks))
    for col, mask in enumerate(masks):
        for (i, j), value in h.todok().items():
            if not mask & (1 << j):
                continue
            removed = mask ^ (1 << j)
            if removed & (1 << i):
                continue
            sign = (-1)**((mask & ((1 << j)-1)).bit_count()
                          +(removed & ((1 << i)-1)).bit_count())
            target = removed | (1 << i)
            if target not in where:
                raise ValueError("the fixed-occupation subspace is not invariant")
            out[where[target], col] += value*sign
    return clean(out)


@lru_cache(maxsize=1)
def local_kernel():
    masks = basis_masks()
    gauges = tuple(one_body_lift(s.diag(s.kronecker_product(s.eye(2), p/2), s.zeros(2)), masks)
                   for p in PAULI)
    kernel = s.Matrix.hstack(*s.Matrix.vstack(*gauges).nullspace())
    gram = clean(kernel.H*kernel)
    left = gram.inv()*kernel.H
    if kernel.cols != 6 or any(g*kernel != s.zeros(12, 6) for g in gauges):
        raise ArithmeticError("the weak-epsilon gauge kernel differs")
    charge = one_body_lift(s.diag(-3, -3, -3, -3, 6, 6), masks)
    if charge != s.zeros(12):
        raise ArithmeticError("the native central charge is not neutral")
    ambient_spin = tuple(one_body_lift(s.diag(s.kronecker_product(p/2, s.eye(2)), p/2), masks)
                         for p in PAULI)
    spin = tuple(clean(left*j*kernel) for j in ambient_spin)
    for full, reduced in zip(ambient_spin, spin):
        if clean(full*kernel-kernel*reduced) != s.zeros(12, 6):
            raise ArithmeticError("spin leaves the gauge kernel")
    for a in range(3):
        if clean(spin[a]*spin[(a+1) % 3]-spin[(a+1) % 3]*spin[a]-s.I*spin[(a+2) % 3]) != s.zeros(6):
            raise ArithmeticError("spin commutator failed")
    casimir = clean(sum((j*j for j in spin), s.zeros(6)))
    p_half = clean((s.Rational(15, 4)*s.eye(6)-casimir)/3)
    p_three_halves = s.eye(6)-p_half
    for p in (p_half, p_three_halves):
        if clean(p*p-p) != s.zeros(6) or clean(p.H*gram-gram*p) != s.zeros(6):
            raise ArithmeticError("Casimir projector failed")
    return masks, kernel, gram, left, spin, casimir, p_half, p_three_halves


def quaternion_elements():
    rows = []
    for j in range(4):
        for sign in (-1, 1):
            rows.append(tuple(s.Integer(sign if i == j else 0) for i in range(4)))
    rows.extend(tuple(s.Rational(sign, 2) for sign in signs)
                for signs in product((-1, 1), repeat=4))
    for pair in combinations(range(4), 2):
        for signs in product((-1, 1), repeat=2):
            row = [s.S.Zero]*4
            for i, sign in zip(pair, signs):
                row[i] = sign*s.sqrt(2)/2
            rows.append(tuple(row))
    if len(set(rows)) != 48:
        raise ArithmeticError("binary cube group has incorrect size")
    return tuple(rows)


def quaternion_product(q, r):
    a, b, c, d = q
    e, f, g, h = r
    return tuple(map(s.expand, (a*e-b*f-c*g-d*h, a*f+b*e+c*h-d*g,
                                a*g-b*h+c*e+d*f, a*h+b*g-c*f+d*e)))


def spinor(q):
    return clean(q[0]*s.eye(2)-s.I*sum((q[a+1]*PAULI[a] for a in range(3)), s.zeros(2)))


def exterior_rotation(u):
    masks, kernel, gram, left, _, _, _, _ = local_kernel()
    one = s.diag(s.kronecker_product(u, s.eye(2)), u)
    occupations = [tuple(i for i in range(6) if mask & (1 << i)) for mask in masks]
    ambient = s.zeros(12)
    for i, target in enumerate(occupations):
        for j, source in enumerate(occupations):
            ambient[i, j] = s.expand(one.extract(target, source).det())
    reduced = clean(left*ambient*kernel)
    if clean(ambient*kernel-kernel*reduced) != s.zeros(12, 6):
        raise ArithmeticError("rotation leaves the gauge kernel")
    if clean(reduced.H*gram*reduced-gram) != s.zeros(6):
        raise ArithmeticError("rotation is not unitary in the physical Gram matrix")
    return reduced


def verify_turns(quarter_turn):
    """An actual quarter turn, four turns and eight turns; no phase quotient."""
    if clean(quarter_turn**4) != -s.eye(6) or clean(quarter_turn**8) != s.eye(6):
        raise ValueError("360-degree sign or 720-degree identity failed")
    return True


def _polynomial(column, site):
    masks, kernel, *_ = local_kernel()
    return tuple((tuple(i+6*site for i in range(6) if mask & (1 << i)), kernel[row, column])
                 for row, mask in enumerate(masks) if kernel[row, column])


def _act_polynomial(state, terms, *, fermionic=True):
    out = {}
    for mask, amplitude in state.items():
        for modes, coefficient in terms:
            target, phase = apply_word(mask, modes, fermionic=fermionic)
            if phase:
                out[target] = out.get(target, 0)+amplitude*coefficient*phase
    return {mask: value for mask, value in out.items() if value != 0}


def verify_exchange(*, fermionic=True, full_domain=True):
    """All 36 kernel-column pairs, on every state of the twelve-mode carrier."""
    polynomials = [[_polynomial(col, site) for col in range(6)] for site in range(2)]
    nonzero = 0
    for a, b in product(polynomials[0], polynomials[1]):
        for mask in range(4096) if full_domain else (0,):
            ab = _act_polynomial(_act_polynomial({mask: 1}, b, fermionic=fermionic), a, fermionic=fermionic)
            ba = _act_polynomial(_act_polynomial({mask: 1}, a, fermionic=fermionic), b, fermionic=fermionic)
            if ab != {key: -value for key, value in ba.items()}:
                raise ValueError("disjoint odd-composite exchange sign failed")
            nonzero += bool(ab)
    if nonzero == 0:
        raise ValueError("exchange check is vacuous")
    return nonzero


def _uniform_orbit_check(q):
    w, x, y, z = q
    r = clean(s.Matrix([[1-2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)],
                        [2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x)],
                        [2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x*x+y*y)]]))
    sites = tuple(product((-1, 1), repeat=3))
    mapped = [tuple(r*s.Matrix(v)) for v in sites]
    return set(mapped) == set(sites) and r.det() == 1


@lru_cache(maxsize=1)
def run():
    masks, kernel, gram, _, spin, casimir, p_half, p_three_halves = local_kernel()
    elements = quaternion_elements()
    element_set = set(elements)
    if any(quaternion_product(q, r) not in element_set for q, r in product(elements, repeat=2)):
        raise ArithmeticError("binary cube multiplication is not closed")
    characters = []
    for q in elements:
        u = spinor(q)
        d = exterior_rotation(u)
        t = s.trace(u)
        c_half, c_three_halves = s.expand(s.trace(p_half*d)), s.expand(s.trace(p_three_halves*d))
        if c_half != t or c_three_halves != s.expand(t**3-2*t):
            raise ArithmeticError("the actual kernel is not G1 plus H")
        if not _uniform_orbit_check(q):
            raise ArithmeticError("the uniform cube orbital factor is not invariant")
        characters.append({"quaternion": [str(x) for x in q], "G1": str(c_half), "H": str(c_three_halves)})
    quarter = exterior_rotation((s.eye(2)-s.I*PAULI[2])/s.sqrt(2))
    verify_turns(quarter)
    nonzero = verify_exchange()
    mutations = {}
    # Multiplication by exp(i*pi/4) changes the fourth power by -1.
    altered = (1+s.I)/s.sqrt(2)*quarter
    try:
        verify_turns(altered)
    except ValueError:
        mutations["central_sign_erasure_rejected"] = True
    try:
        verify_exchange(fermionic=False, full_domain=False)
    except ValueError:
        mutations["commuting_creator_replacement_rejected"] = True
    if len(mutations) != 2 or not all(mutations.values()):
        raise ArithmeticError("a claimed broken-input test did not fail")
    checks = {
        "exact_12_weight_ambient_space": len(masks) == 12,
        "six_dimensional_native_weak_singlet_kernel": kernel.cols == 6,
        "all_spin_commutators_and_kernel_invariance": True,
        "spin_half_projector_rank_two": p_half.rank() == 2,
        "spin_three_halves_projector_rank_four": p_three_halves.rank() == 4,
        "all_48_binary_cube_characters_G1_plus_H": True,
        "uniform_eight_site_orbit_invariant": True,
        "full_turn_minus_identity": True, "double_turn_identity": True,
        "all_36_composite_pairs_on_all_4096_Fock_states_anticommute": True,
    }
    if not all(checks.values()):
        raise ArithmeticError("a spin claim failed")
    return {
        "schema": "oph.openwave.spin.v1", "computational_pass": True,
        "criterion_rows": [7], "full_criterion_met": False,
        "checks": checks, "mutations": mutations, "source_refs": SOURCE_REFS,
        "source_restriction": "positive weak-pair/charge-six composite, native increment (0,2,0,0,1), canonical equivalent spin and weak charts",
        "native_central_charges": [-3, -3, -3, -3, 6, 6],
        "ambient_masks": list(masks), "kernel_basis": [[str(v) for v in row] for row in kernel.tolist()],
        "kernel_gram": [[str(v) for v in row] for row in gram.tolist()],
        "spin_generators": [[[str(v) for v in row] for row in j.tolist()] for j in spin],
        "Casimir_multiplicities": {str(k): int(v) for k, v in casimir.eigenvals().items()},
        "binary_cube_characters": characters,
        "exchange_basis_comparisons": 36*4096, "nonzero_exchange_comparisons": nonzero,
        "parameter_ancestry": {
            "declared": ["the primary source spinor action on each fermion", "exterior CAR quantization",
                         "native weak doublet and charge-six singlet representations",
                         "weak-pair/charge-six occupation sector", "uniform cube orbital preparation"],
            "empirical_calibrations_used": [],
        },
        "scope": {
            "licenses": "Exact double-cover sign and graded disjoint-support exchange of the computed finite gauge-singlet composite kernel, including its two-dimensional G1 component.",
            "does_not_license": "The CAR and spinor action are inputs. No classical emergence of anticommutation, relativistic spin-statistics theorem, physical electron assignment, source-Hamiltonian persistence, mass or laboratory interference experiment is supplied.",
            "full_primary_Hamiltonian_executed": False,
            "continuous_spatial_rotation_symmetry_derived": False,
            "same_site_free_CAR_claimed": False,
            "official_OpenWave_cell_assigned": False,
        },
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
