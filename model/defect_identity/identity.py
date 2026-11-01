"""Exact bridge from the native ordered cover to the declared matter module.

This is a capacity theorem and an explicit permitted memory design.
The additional spectator fiber is specified. All finite algebra is exact.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import sympy as sp

from ..defect_protection import protection

ROOT = Path(__file__).resolve().parents[2]
UPSTREAM = ROOT / 'code' / 'a5_closure'
if str(UPSTREAM) not in sys.path:
    sys.path.insert(0, str(UPSTREAM))
import super_tannakian_matter_lift_certificate as exterior  # noqa: E402


def _scalar(value):
    return (sp.Rational(value.re.a) + sp.sqrt(5)*sp.Rational(value.re.b)
            + sp.I*(sp.Rational(value.im.a) + sp.sqrt(5)*sp.Rational(value.im.b)))


@lru_cache(maxsize=1)
def matter_basis():
    """Reuse the existing auxiliary CAR basis: even parity minus the vacuum.

    Lambda^2(C^5) plus Lambda^4(C^5), not a newly invented label table.
    The conjugate odd-minus-top representative is not physically excluded.
    """
    fock = exterior.FockSpace()
    return tuple(subset for subset in fock.subsets if len(subset) in (2, 4))


def _field(subset):
    color = sum(i < 3 for i in subset)
    weak = len(subset)-color
    return {(2, 0): 'u_c', (1, 1): 'Q', (0, 2): 'e_c',
            (2, 2): 'd_c', (3, 1): 'L'}[(color, weak)]


def matter_action(one_body):
    """Existing FockSpace.dgamma, restricted to the actual fifteen subsets."""
    if one_body.shape != (5, 5):
        raise ValueError('the declared carrier has five modes')
    exact = []
    for row in one_body.tolist():
        converted = []
        for entry in row:
            re, im = sp.expand(entry).as_real_imag()
            if not re.is_Rational or not im.is_Rational:
                raise ValueError('rational complex one-body entries required')
            converted.append(exterior.C5(exterior.F5(str(re)), exterior.F5(str(im))))
        exact.append(converted)
    fock = exterior.FockSpace()
    result = fock.dgamma(exact)
    indices = [fock.index[s] for s in matter_basis()]
    return sp.ImmutableMatrix([[_scalar(result[i][j]) for j in indices] for i in indices])


@lru_cache(maxsize=1)
def generators():
    """An exact skew-Hermitian basis of s(u(3)+u(2)), lifted by dGamma."""
    values = []
    for block in (tuple(range(3)), (3, 4)):
        for i, j in combinations(block, 2):
            real, imaginary = sp.zeros(5), sp.zeros(5)
            real[i, j], real[j, i] = 1, -1
            imaginary[i, j] = imaginary[j, i] = sp.I
            values.extend((real, imaginary))
        for i in block[:-1]:
            diagonal = sp.zeros(5)
            diagonal[i, i], diagonal[block[-1], block[-1]] = sp.I, -sp.I
            values.append(diagonal)
    values.append(sp.I*sp.diag(-2, -2, -2, 3, 3))  # q=6Y
    return tuple((sp.ImmutableMatrix(g), matter_action(g)) for g in values)


@lru_cache(maxsize=1)
def global_faithfulness_certificate():
    """Prove the exterior-dual identity for a generic matrix, not samples.

    H is the volume pairing Lambda^4 V -> V*. The polynomial identity
    g^T H Lambda^4(g) H^-1 = det(g) I holds without invertibility. On
    det(g)=1 it proves the actual group representation is the faithful dual.
    """
    def determinant(matrix):
        n = matrix.rows
        return sp.Add(*(sp.prod(matrix[i, p[i]] for i in range(n)) *
                        (-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1, n))
                        for p in permutations(range(n))))
    g = sp.Matrix(5, 5, sp.symbols('g0:25'))
    basis = [s for s in matter_basis() if len(s) == 4]
    H = sp.zeros(5)
    for column, subset in enumerate(basis):
        missing = next(i for i in range(5) if i not in subset)
        H[missing, column] = (-1)**missing
    fourth = sp.Matrix([[determinant(g.extract(a, b)) for b in basis] for a in basis])
    residual = g.T*H*fourth*H.T-determinant(g)*sp.eye(5)
    if H*H.T != sp.eye(5) or any(sp.expand(value) != 0 for value in residual):
        raise ArithmeticError('global exterior-dual polynomial identity failed')
    return {'identity': 'g^T H Lambda^4(g) H^-1 = det(g) I5',
            'generic_matrix_variables': 25, 'exact_polynomial_entries_checked': 25,
            'H_volume_pairing': [[int(x) for x in row] for row in H.tolist()],
            'conclusion_on_det_one': 'Lambda^4(g)=I implies g^T=I, hence g=I; the full matter action is group-faithful'}


def involution(a, b):
    """Conjugacy representative in S(U3 x U2), negative multiplicities a,b."""
    if (type(a) is not int or type(b) is not int or not 0 <= a <= 3
            or not 0 <= b <= 2 or (a+b) % 2):
        raise ValueError('involution multiplicities must obey det(A)det(B)=1')
    signs = tuple([-1]*a+[1]*(3-a)+[-1]*b+[1]*(2-b))
    weights = tuple(sp.prod(signs[i] for i in subset) for subset in matter_basis())
    blocks = {}
    for subset, sign in zip(matter_basis(), weights, strict=True):
        row = blocks.setdefault(_field(subset), {'plus': 0, 'minus': 0})
        row['plus' if sign == 1 else 'minus'] += 1
    return {'negative_color': a, 'negative_weak': b,
            'one_body_signs': list(signs), 'matter_signs': [int(s) for s in weights],
            'plus': weights.count(1), 'minus': weights.count(-1),
            'central': a in (0, 3) and b in (0, 2), 'blocks': blocks}


def involution_classes():
    return [involution(a, b) for a in range(4) for b in range(3) if (a+b) % 2 == 0]


def embedding(signs, multiplicity):
    """Return W=sqrt(2)V with W*W=2I and XW=WJ, or reject capacity.

    This alone is Z2 equivariance; full gauge equivariance additionally needs
    J central. The certificate constructs that stronger case separately.
    """
    if (type(multiplicity) is not int or multiplicity < 1
            or not signs or any(type(s) is not int or s not in (-1, 1) for s in signs)):
        raise ValueError('positive exact fiber dimension and Z2 signs required')
    if max(signs.count(1), signs.count(-1)) > multiplicity:
        raise ValueError('spectator fiber has insufficient capacity in a charge sector')
    matrix = sp.zeros(2*multiplicity, len(signs))
    positions = {1: 0, -1: 0}
    for column, sign in enumerate(signs):
        slot = positions[sign]
        positions[sign] += 1
        matrix[slot, column] = 1
        matrix[multiplicity+slot, column] = sign
    return sp.ImmutableMatrix(matrix)


def validate_embedding(matrix, signs, multiplicity):
    """Validate an arbitrary proposed isometric intertwiner, not just ours."""
    if matrix.shape != (2*multiplicity, len(signs)):
        raise ValueError('carrier dimensions do not match')
    identity = sp.eye(len(signs))
    flip = sp.kronecker_product(sp.Matrix([[0, 1], [1, 0]]), sp.eye(multiplicity))
    if matrix.H*matrix != 2*identity:
        raise ValueError('the proposed matter map loses norm or merges states')
    if flip*matrix != matrix*sp.diag(*signs):
        raise ValueError('the proposed map does not intertwine the native charge')


def native_intertwining(law, matrix, signs, multiplicity):
    """Check every source column, zero/abort slot, and matter basis vector.

    Amplitude 1/sqrt(2) is common to both complete branch maps and cancels.
    The physical ordering is (orientation, source, spectator). No native
    source or memory map is replaced, and no postselection enters the check.
    """
    validate_embedding(matrix, signs, multiplicity)
    columns = [tuple((r, slot, matrix[r*multiplicity+slot, column])
                     for r in (0, 1) for slot in range(multiplicity)
                     if matrix[r*multiplicity+slot, column])
               for column in range(len(signs))]
    total = 0
    rows = protection.branches(law)
    for row in rows:
        lifted = dict(protection.lifted_pairs(row))
        source_map = {i: (j, sign) for i, j, sign in
                      zip(row.inputs, row.outputs, row.signs, strict=True)}
        for source in range(66):
            for column, character in enumerate(signs):
                actual = {}
                for r, slot, coefficient in columns[column]:
                    destination = lifted.get(r*66+source)
                    if destination is None:
                        continue
                    actual[destination*multiplicity+slot] = coefficient
                expected = {}
                if source in source_map:
                    target, parity = source_map[source]
                    phase = character if parity == -1 else 1
                    for r, slot, coefficient in columns[column]:
                        expected[(r*66+target)*multiplicity+slot] = phase*coefficient
                if actual != expected:
                    raise ArithmeticError('native complete branch intertwiner failed')
                total += 1
    return {'raw_slots': len(rows), 'all_source_and_matter_columns': total,
            'zero_slots_retained': sum(not r.inputs for r in rows),
            'abort_slots_retained': sum(r.abort for r in rows),
            'complete_native_branch_intertwining': True}


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def negative_character_extension():
    """Keep the native character negative on one collective 15-label module.

    Q is a separate statistics factor, not an element of the internal G.
    The supplied fiber is the existing matter module and its gauge action
    is I_memory tensor rho(g). This does not yet construct a species Fock
    space: it is an exact same-carrier collective-module/two-token assignment.
    The full per-particle colored construction uses -SWAP and both characters.
    """
    signs, multiplicity = [-1]*15, 15
    W = embedding(signs, multiplicity)
    validate_embedding(W, signs, multiplicity)
    if W != sp.Matrix.vstack(sp.eye(15), -sp.eye(15)):
        raise ArithmeticError('uniform-fermion inclusion is not the full negative sector')
    Q = sp.kronecker_product(sp.Matrix([[0, 1], [1, 0]]), sp.eye(15))
    gs = generators()
    action = [sp.kronecker_product(sp.eye(2), g) for _, g in gs]
    for (_, g), extended in zip(gs, action, strict=True):
        if extended*W != W*g or extended*Q != Q*extended or extended.H != -extended:
            raise ArithmeticError('uniform-fermion gauge action failed')
    for i in range(12):
        for j in range(12):
            if action[i]*action[j]-action[j]*action[i] != sp.kronecker_product(
                    sp.eye(2), gs[i][1]*gs[j][1]-gs[j][1]*gs[i][1]):
                raise ArithmeticError('uniform-fermion gauge Lie bracket failed')
    return {
        'required_character_on_all_collective_labels': -1,
        'matter_exchange_signature': [0, 15],
        'necessary_capacity': 'm >= 15 when one collective fifteen-label module is required to lie entirely in the negative memory sector',
        'chosen_m': multiplicity, 'native_extended_dimension': 66*2*multiplicity,
        'source_times_matter_dimension': 66*15,
        'W_sqrt2_times_isometry': [[int(x) for x in row] for row in W.tolist()],
        'exact_gram_W_star_W': '2 I15', 'exact_exchange_intertwiner': 'X_memory W = -W',
        'full_group_action': 'I2 tensor rho(g)',
        'statistics_action': 'Z2 factor commuting with G; Q is not identified with an internal gauge element',
        'all_12_gauge_generators_intertwined': True, 'all_144_Lie_brackets_preserved': True,
        'positive_exchange_sector': 'an unused second copy of the same internal module; excluded by the chosen negative-sector encoding',
        'initialization': 'prepare the ordering qubit as |-> from |0> by H then Z and tensor the internal state; unitary design choice, not postselection',
        'full_carrier_retained': 'both fifteen-dimensional exchange sectors remain in C2 tensor C15; native branches and gauge operations preserve each sector',
        'complete_branch_checks': {law: native_intertwining(law, W, signs, multiplicity) for law in ('C', 'E')},
        'design_status': 'explicit permitted observer-memory design; no uniqueness principle is needed',
        'scope': 'one collective internal module on the native two-token configuration carrier; not a universal fermionic memory bound. Per-particle internal labels require the colored -SWAP construction, which uses both characters.',
    }


def certificate():
    contract = json.loads((UPSTREAM/'manifests/super_tannakian_matter_reference.json').read_text())
    charges = contract['exterior_matter_contract']['block_trace_charges']
    if (sp.Rational(charges['color_block']) != -sp.Rational(1, 3)
            or sp.Rational(charges['weak_block']) != sp.Rational(1, 2)):
        raise ArithmeticError('the upstream declared matter charge convention changed')
    classes = involution_classes()
    central = involution(0, 2)
    signs = central['matter_signs']
    multiplicity = max(signs.count(1), signs.count(-1))
    matrix = embedding(signs, multiplicity)
    validate_embedding(matrix, signs, multiplicity)
    charge = sp.diag(*signs)
    flip = sp.kronecker_product(sp.Matrix([[0, 1], [1, 0]]), sp.eye(multiplicity))
    gs = generators()
    if len(gs) != 12 or sp.Matrix.hstack(*(g.reshape(225, 1) for _, g in gs)).rank() != 12:
        raise ArithmeticError('the exact matter action lost faithfulness')
    compatible = []
    for row in classes:
        J = sp.diag(*row['matter_signs'])
        commutes = all(J*g == g*J for _, g in gs)
        if commutes != row['central']:
            raise ArithmeticError('group-central and representation-central involutions disagree')
        if row['minus']:
            if (row['plus'], row['minus']) != (7, 8):
                raise ArithmeticError('nontrivial involution signature changed')
            if commutes:
                compatible.append([row['negative_color'], row['negative_weak']])
    if compatible != [[0, 2]]:
        raise ArithmeticError('the nontrivial central involution classification changed')
    if any(row['minus'] == 15 for row in classes):
        raise ArithmeticError('unexpected internal gauge representative of uniform fermion parity')
    embedded = [matrix*g*matrix.H/2 for _, g in gs]
    for (_, g), extended in zip(gs, embedded, strict=True):
        if g.H != -g or charge*g != g*charge or extended*matrix != matrix*g:
            raise ArithmeticError('full commuting internal gauge action failed')
        if flip*extended != extended*flip:
            raise ArithmeticError('internal action fails native charge conservation')
    for i in range(12):
        for j in range(12):
            original = gs[i][1]*gs[j][1]-gs[j][1]*gs[i][1]
            lifted = embedded[i]*embedded[j]-embedded[j]*embedded[i]
            if lifted != matrix*original*matrix.H/2:
                raise ArithmeticError('the extended action fails the Lie bracket')
    # A group extension is V rho(g) V* plus the identity on the complement.
    complement = sp.eye(16)-matrix*matrix.H/2
    if complement.rank() != 1 or flip*complement != complement:
        raise ArithmeticError('the minimal extension complement is not one even line')
    if matrix*charge*matrix.H/2+complement != flip:
        raise ArithmeticError('the full group extension does not map J to the native flip')
    native = {law: protection.algebra_certificate(law) for law in ('C', 'E')}
    if any(x['commutant_complex_dimension'] != 2 for x in native.values()):
        raise ArithmeticError('the asserted native commutant premise changed')
    files = ['model/defect_identity/identity.py',
             'model/defect_protection/protection.py',
             'code/a5_closure/super_tannakian_matter_lift_certificate.py',
             'code/a5_closure/manifests/super_tannakian_matter_reference.json']
    return {
        'schema': 'oph.defect_identity.v1',
        'arithmetic': 'exact rational complex matrices; no numerical eigenspace thresholds',
        'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in files},
        'group_convention': {
            'group': 'G=S(U(3) x U(2)) isomorphic to (SU(3) x SU(2) x U(1))/Z6',
            'cover_map': '(A,B,z) -> (z^-2 A,z^3 B)',
            'kernel_generator': '(exp(2*pi*i/3) I3,-I2,exp(pi*i/3))',
            'hypercharge': 'q=6Y; q(one-body)=(-2,-2,-2,3,3)',
            'conditional': 'This is the faithful image of the specified matter design. Relating it to the physical global form remains an identification question.'},
        'matter': {'dimension': 15, 'basis_subsets': [list(s) for s in matter_basis()],
                   'representative': 'Lambda^2(C^5) plus Lambda^4(C^5); even-minus-vacuum',
                   'conjugate_representative_not_excluded': True,
                   'q_weights': [sum(-2 if i < 3 else 3 for i in s) for s in matter_basis()],
                   'generator_rank': 12,
                   'generator_sha256': _digest([[str(x) for x in g] for _, g in gs]),
                   'group_faithfulness': global_faithfulness_certificate()},
        'all_involution_conjugacy_classes': classes,
        'only_nonidentity_central_class': compatible[0],
        'centrality_reason': 'A commuting internal G action must commute with Q; if Q restricts to rho(J), faithful rho forces J central.',
        'internal_gauge_uniform_fermion_parity': {
            'exists': False,
            'proof': 'rho(g)=-I15 would imply rho(g^2)=I; global faithfulness gives g^2=I. All six involution classes have minus dimension 0 or 8, never 15.',
            'weak_parity_is_uniform_negative_character': False,
            'weak_parity_exchange_assignment': 'seven positive and eight negative native exchange characters; incompatible with an imposed uniformly negative collective fifteen-label encoding'},
        'native': {law: {'raw_maps_sha256': row['raw_maps_sha256'],
                          'commutant_complex_dimension': row['commutant_complex_dimension'],
                          'odd_closed_native_history': json.loads(json.dumps(row['odd_closed_native_history']))}
                   for law, row in native.items()},
        'unextended_obstruction': 'Native commutant C direct_sum C is abelian; no nontrivial su(3) or su(2) action commutes with every complete native branch.',
        'spectator_extension': {
            'scope': 'internal weak-grading control only; not a uniformly negative collective-label encoding or a full per-particle exchange construction',
            'supplied_additional_factor': 'C^m, untouched by every native branch',
            'commutant': 'M_m(C) direct_sum M_m(C)',
            'necessary_capacity_for_nonidentity_gauge_involution': 'm >= max(7,8) = 8',
            'minimal_feasible_m_for_this_declared_extension': multiplicity,
            'native_extended_dimension': 66*2*multiplicity,
            'source_times_matter_dimension': 66*15,
            'W_sqrt2_times_isometry': [[int(x) for x in row] for row in matrix.tolist()],
            'exact_gram_W_star_W': '2 I15', 'exact_intertwiner': 'X_memory W = W rho(I3,-I2)',
            'all_144_Lie_brackets_preserved': True,
            'group_extension': 'rho_ext(g)=(W rho(g) W*)/2 + (I16-W W*/2)',
            'complement': 'one even gauge-trivial line in this chosen minimal construction; no particle identity asserted',
            'fiber_and_gauge_action_status': 'explicit permitted design; no uniqueness or natural-selection principle is required',
            'implementation_status': 'exact algebraic extension of every native branch; execution and resource accounting for the added fiber remain to be supplied'},
        'complete_branch_checks': {law: native_intertwining(law, matrix, signs, multiplicity) for law in ('C', 'E')},
        'uniform_negative_character_extension': negative_character_extension(),
        'quantifiers': [
            'Both supplied native C/E ordered-cover instruments, with all zero and abort slots retained.',
            'Every finite history and adaptive choice of these branches; arbitrary untouched reference, by exact branch intertwining.',
            'Internal transformations commute with each fixed complete branch, not merely an averaged channel or relabeled family.',
            'The additional C^m fiber is supplied, and the retained 15D matter subspace is invariant under a faithful action of the declared G.',
            'For the internal-grading subtheorem only, the native deck involution is identified with the image of a nonidentity J in G; its mixed signature excludes a uniformly negative collective encoding.',
            'For the uniform-character design, Q=-I on one collective fifteen-label module; the chosen spectator has dimension fifteen. This is not the full two-particle internal space.'
        ],
        'not_discharged': ['execution, complete failures, resources and memory-noise control for the specified extra fiber and gauge operations',
                           'same-excitation localization, propagation, particle pole, mass, spin or statistics observables',
                           'empirical identification of the specified matter representation and global gauge image',
                           'actual observer memory protection against laws that fail to commute with X',
                           'full multiparticle fifteen-species Fock space and Lorentz spin-statistics connection',
                           'three generations or neutrino identity'],
        'physical_identification': False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    payload = certificate()
    path = Path(__file__).with_name('results.json')
    if args.check:
        if json.loads(path.read_text()) != payload:
            raise SystemExit('defect identity receipt drift')
        print('defect identity exact receipt OK')
    elif args.write:
        path.write_text(json.dumps(payload, indent=2)+'\n')
        print(f'wrote {path}')
    else:
        print(json.dumps(payload, indent=2))


if __name__ == '__main__':
    main()
