"""Source-bound phase nonidentifiability and exact cyclic-symmetry controls.

No experimental particle data enter this module. The examples distinguish a
complete classical record instrument from a chosen coherent record erasure.
Universal arguments, domains and nonclaims are stated in PROOF.md.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as s

from model.common_source import native_family_record_memory as memory
from model.common_source import native_gauge_object as geometry

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCES = (
    'model/common_source/native_family_record_memory.py',
    'model/common_source/native_repair_mass.py',
    'model/common_source/native_gauge_object.py',
    'model/particle_source_selection/selection.py',
    'model/particle_source_selection/PROOF.md',
)


def clean(m):
    return s.Matrix(m).applyfunc(s.simplify)


def cyclic_rotation():
    return s.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])


def _unit_phase(z):
    z = s.sympify(z)
    if z.has(s.Float) or z.free_symbols or s.simplify(z*s.conjugate(z)-1) != 0:
        raise ValueError('an exactly certified unit complex number is required')
    return z


def coherent_lift(z):
    """Complete diagnostic lift, not a native record implementation.

    Four raw outcomes retain weights q,p,p,1-q-2p. A selected coherent
    erasure of the first three records has Kraus operator C/sqrt(3).
    The fourth outcome is retained throughout.
    """
    z = _unit_phase(z)
    d = memory.architecture()
    p, q, R = d['p'], d['q'], cyclic_rotation()
    raw = (s.sqrt(q)*s.eye(3), s.sqrt(p)*z*R,
           s.sqrt(p)*s.conjugate(z)*R.T, s.sqrt(1-q-2*p)*s.eye(3))
    C = clean(sum(raw[:3], s.zeros(3)))
    # A unitary Fourier change of the first three environment outputs.
    w = -s.Rational(1, 2)+s.I*s.sqrt(3)/2
    erased = tuple(clean(sum((w**(j*k)*raw[k] for k in range(3)), s.zeros(3))/s.sqrt(3))
                   for j in range(3))+(raw[3],)
    return {'raw': raw, 'erased': erased, 'C': C, 'p': p, 'q': q}


def superoperator(K):
    """Full matrix-space channel; equality also preserves any reference."""
    return clean(s.kronecker_product(K, K.conjugate()))


def phase_lift_certificate():
    # Two exact unit phases, chosen as algebraic controls without particle data.
    phases = (s.Integer(1), (99+20*s.I)/101)
    a, b = (coherent_lift(z) for z in phases)
    identity = s.eye(3)
    for lift in (a, b):
        assert clean(sum((K.H*K for K in lift['raw']), s.zeros(3))) == identity
        assert clean(sum((K.H*K for K in lift['erased']), s.zeros(3))) == identity
        assert clean(sum((superoperator(K) for K in lift['raw']), s.zeros(9))) == clean(
            sum((superoperator(K) for K in lift['erased']), s.zeros(9)))
        C = lift['C']
        assert C == C.H
        # Sylvester's criterion gives strict positivity, without float gates.
        assert all(s.simplify(C[:i, :i].det()).is_positive is True for i in range(1, 4))
        assert s.simplify(s.trace(C*C)/s.trace(C)**2) == s.Rational(2, 3)
    # Outcome-by-outcome equality, stronger than equality of unread channels.
    assert all(superoperator(Ka) == superoperator(Kb)
               for Ka, Kb in zip(a['raw'], b['raw'], strict=True))
    x = s.Symbol('x')
    polys = [s.factor(lift['C'].charpoly(x).as_expr()) for lift in (a, b)]
    assert s.simplify(polys[0]-polys[1]) != 0
    assert s.factor(a['C'].charpoly(x).discriminant()) == 0
    assert s.factor(b['C'].charpoly(x).discriminant()) != 0
    # A real input vector witnesses different probabilities after coherent erasure.
    v = s.ones(3, 1)/s.sqrt(3)
    probability = [s.simplify((v.T*lift['erased'][0].H*lift['erased'][0]*v)[0])
                   for lift in (a, b)]
    assert s.simplify(probability[0]-probability[1]).is_positive is True
    return {
        'raw_complete_record_instruments_identical': True,
        'reference_safe': True,
        'all_erasure_outcomes_retained': True,
        'native_implementation_of_coherent_lift_claimed': False,
        'phases': [str(z) for z in phases],
        'positive_square_root_operators': True,
        'koide_invariant_both': '2/3',
        'characteristic_polynomials': [str(v) for v in polys],
        'distinct_spectra': True,
        'coherent_success_probabilities_uniform_input': [str(v) for v in probability],
        'universal_statement': 'independent Kraus branch phases cancel from every retained classical outcome map',
        'closure_boundary': 'classical outcome laws do not select a coherent erasure phase',
    }


def real_cyclic_twirl(X, R):
    """Exact twirl for a real symmetric matrix and real regular C3 action."""
    X, R = s.Matrix(X), s.Matrix(R)
    if X.shape != (3, 3) or R.shape != (3, 3):
        raise ValueError('three-dimensional family matrices required')
    if X != X.T or any(v.is_real is not True for v in X):
        raise ValueError('an exactly real symmetric kernel is required')
    if any(v.is_real is not True for v in R) or clean(R.T*R) != s.eye(3):
        raise ValueError('a real orthogonal cyclic action is required')
    if clean(R**3) != s.eye(3) or s.simplify(s.trace(R)) != 0:
        raise ValueError('a regular order-three action is required')
    twirl = clean(sum((R**j*X*(R.T)**j for j in range(3)), s.zeros(3))/3)
    fixed = clean((s.eye(3)+R+R**2)/3)
    singlet = s.simplify(s.trace(fixed*X))
    doublet = s.simplify((s.trace(X)-singlet)/2)
    assert clean(twirl-singlet*fixed-doublet*(s.eye(3)-fixed)) == s.zeros(3)
    return twirl, singlet, doublet


@lru_cache(maxsize=1)
def source_certificate():
    d = memory.architecture()
    p, q, joint, Y = d['p'], d['q'], d['coincidence'], d['Y']
    assert (p, q, joint) == (s.Rational(1, 12), s.Rational(1, 6), s.Rational(1, 72))
    assert s.simplify(2*s.sqrt(p/q)-s.sqrt(2)) == 0
    assert p+q-2*joint == s.Rational(2, 9)
    assert d['T']*d['stationary'] == d['stationary']
    native = geometry.geometry()
    B = native['frame']
    norm = s.simplify((B.T*B)[0, 0])
    assert clean(B.T*B) == norm*s.eye(3)
    rows = []
    for permutation in sorted(native['rotations']):
        if geometry.permutation_order(permutation) != 3:
            continue
        R = clean(B.T*geometry.permutation_matrix(permutation)*B/norm)
        J = clean(Y*R-R*Y)
        assert J != s.zeros(3)
        twirl, singlet, doublet = real_cyclic_twirl(Y, R)
        x = s.Symbol('x')
        assert s.simplify(twirl.charpoly(x).as_expr()-(x-singlet)*(x-doublet)**2) == 0
        rows.append({'permutation': list(permutation),
                     'commutator_norm_squared': str(s.simplify(s.trace(J.T*J))),
                     'twirl_eigenvalues': [str(singlet), str(doublet), str(doublet)]})
    assert len(rows) == 20
    x = s.Symbol('x')
    discriminant = s.factor(Y.charpoly(x).discriminant())
    assert discriminant != 0
    return {
        'event_factors': {'p': str(p), 'q': str(q), 'joint': str(joint)},
        'rate_identities': {'2sqrt(p/q)': 'sqrt(2)', 'p+q-2joint': '2/9'},
        'all_raw_labels_retained': len(d['records']),
        'native_kernel_discriminant': str(discriminant),
        'native_kernel_has_simple_spectrum': True,
        'regular_cyclic_actions_checked': len(rows),
        'native_kernel_commutes_with_any_tested_C3': False,
        'every_real_C3_twirl_has_repeated_eigenvalue': True,
        'cyclic_actions': rows,
        'scope': 'declared native family record memory and its twenty inherited order-three rotations',
    }


def receipt():
    # The six-variable identity checks the complete real symmetric matrix space.
    a, b, c, d, e, f = s.symbols('a b c d e f', real=True)
    _, singlet, doublet = real_cyclic_twirl(s.Matrix([[a, d, e], [d, b, f], [e, f, c]]), cyclic_rotation())
    return {
        'schema': 'oph.particle_source_selection.v1',
        'empirical_inputs': [],
        'measured_spectrum_claimed': False,
        'native_record_certificate': source_certificate(),
        'phase_lift_counterfamily': phase_lift_certificate(),
        'general_real_twirl': {'singlet': str(singlet), 'doublet': str(doublet), 'doublet_multiplicity': 2},
        'source_bindings': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    text = json.dumps(receipt(), indent=2, sort_keys=True)+'\n'
    path = HERE/'results.json'
    if args.check:
        if not path.is_file() or path.read_text() != text:
            raise SystemExit('source selection receipt differs from exact replay')
        print('source selection receipt: exact replay passed')
    else:
        path.write_text(text)
        print(path.relative_to(ROOT))


if __name__ == '__main__':
    main()
