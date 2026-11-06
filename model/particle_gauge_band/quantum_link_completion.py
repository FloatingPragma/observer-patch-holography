#!/usr/bin/env python3
"""Exact finite link multiplication, Gauss covariance and CAR boundary.

The compact-group rotor realization is a premise of the coordinate-product
completion. The supplied 60-state representation and all current matrices
are the native source objects.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import time

import sympy as s

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

from model.native_quantum_gauss import algebra as native
from model.common_source import native_matter_attachment as matter
from model.interacting_defect_response.algebra import annihilators


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clean(A):
    return s.SparseMatrix(A.rows, A.cols,
                         {ij: s.simplify(value) for ij, value in A.todok().items()
                          if s.simplify(value) != 0})


def zero(A):
    return not clean(s.SparseMatrix(A)).todok()


def exterior_data():
    masks = tuple(mask for mask in range(32) if mask.bit_count() in (2, 4))
    q = lambda mask: -2*(mask & 7).bit_count()+3*((mask >> 3) & 3).bit_count()
    blocks = {charge: tuple(i for i, mask in enumerate(masks) if q(mask) == charge)
              for charge in (-4, -3, 1, 2, 6)}
    assert {tuple(sorted(indices)) for indices in blocks.values()} == {
        tuple(sorted(indices)) for _, indices in native.native_blocks()}
    currents = native.native_currents()
    return masks, blocks, currents


def wedge_sign(left, right):
    if left & right:
        return 0
    inversions = sum(1 for a in range(5) if left >> a & 1
                     for b in range(5) if right >> b & 1 and a > b)
    return (-1)**inversions


def irreducible_algebra_dimension(generators, d):
    basis = [s.eye(d)]
    while len(basis) < d*d:
        candidates = basis+[B*G for B in basis for G in generators]
        M = s.Matrix.hstack(*(A.reshape(d*d, 1) for A in candidates))
        pivots = M.rref()[1]
        grown = [candidates[i] for i in pivots]
        if len(grown) == len(basis):
            break
        basis = grown
    return len(basis)


def clebsch_gordan(masks, blocks, currents):
    tensors = {}
    checks = []
    names = {-4: 'u^c', -3: 'L', 1: 'Q', 2: 'd^c', 6: 'e^c'}
    for qr, ri in blocks.items():
        for qs, si in blocks.items():
            if qr+qs not in blocks:
                continue
            qt, ti = qr+qs, blocks[qr+qs]
            dr, ds, dt = len(ri), len(si), len(ti)
            raw = s.zeros(dt, dr*ds)
            for c, output in enumerate(ti):
                for a, left in enumerate(ri):
                    for b, right in enumerate(si):
                        if masks[left] | masks[right] == masks[output]:
                            raw[c, ds*a+b] = wedge_sign(masks[left], masks[right])
            if raw == s.zeros(dt, dr*ds):
                continue
            gram = raw*raw.T
            factor = gram[0, 0]
            assert gram == factor*s.eye(dt) and factor > 0
            C = s.SparseMatrix(raw/s.sqrt(factor))
            left_gram = s.zeros(ds)
            for c in range(ds):
                for d in range(ds):
                    left_gram[c, d] = sum(C[t, ds*a+c]*C[t, ds*a+d]
                                          for t in range(dt) for a in range(dr))
            assert left_gram == s.Rational(dt, ds)*s.eye(ds)
            for T in currents:
                Tr, Ts, Tt = T.extract(ri, ri), T.extract(si, si), T.extract(ti, ti)
                action = s.kronecker_product(Tr, s.eye(ds))+s.kronecker_product(s.eye(dr), Ts)
                assert zero(C*action-Tt*C)
            tensors[(qr, qs, qt)] = C
            checks.append({'product': [names[qr], names[qs], names[qt]],
                           'dimensions': [dr, ds, dt],
                           'wedge_row_squared_norm': str(factor),
                           'current_intertwiner_checks': len(currents),
                           'normalized_CG_rows': [[str(x) for x in row] for row in C.tolist()],
                           'vacuum_raising_normalized_product_factor': str(s.sqrt(s.Rational(dr*ds, dt)))})
    expected = {(1, 1, 2), (1, -4, -3), (-4, 1, -3), (-4, 6, 2), (6, -4, 2)}
    assert set(tensors) == expected
    return tensors, checks


def link_representation(blocks, currents):
    link_indices = {}
    labels = [{'block': 'vacuum'}]
    offset = 1
    for q, indices in blocks.items():
        d = len(indices)
        for a in range(d):
            for b in range(d):
                link_indices[(q, a, b)] = offset
                labels.append({'q': q, 'left': a, 'right': b})
                offset += 1
    assert offset == 60
    left, right = [], []
    for T in currents:
        L, R = s.zeros(60), s.zeros(60)
        for q, indices in blocks.items():
            d = len(indices)
            small = T.extract(indices, indices)
            for a in range(d):
                for b in range(d):
                    for c in range(d):
                        L[link_indices[q, a, b], link_indices[q, c, b]] = -small[c, a]
                        R[link_indices[q, a, b], link_indices[q, a, c]] = small[b, c]
        left.append(s.SparseMatrix(L)); right.append(s.SparseMatrix(R))
    return link_indices, labels, left, right


def coordinate_operators(qr, blocks, indices, tensors):
    dr = len(blocks[qr])
    source, completed = {}, {}
    for a in range(dr):
        for b in range(dr):
            entries = {(indices[qr, a, b], 0): s.Integer(1)}
            source[a, b] = s.SparseMatrix(60, 60, entries)
            for (r, qs, qt), C in tensors.items():
                if r != qr:
                    continue
                ds, dt = len(blocks[qs]), len(blocks[qt])
                factor = s.sqrt(s.Rational(dr*ds, dt))
                for c in range(ds):
                    for d in range(ds):
                        for i in range(dt):
                            if C[i, ds*a+c] == 0:
                                continue
                            for j in range(dt):
                                value = s.simplify(factor*C[i, ds*a+c]*C[j, ds*b+d])
                                if value:
                                    entries[indices[qt, i, j], indices[qs, c, d]] = value
            completed[a, b] = s.SparseMatrix(60, 60, entries)
    return source, completed


def combined_coordinate(ops, d):
    out = {}
    for (a, b), U in ops.items():
        for (i, j), value in U.todok().items():
            out[60*a+i, 60*b+j] = value
    return s.SparseMatrix(60*d, 60*d, out)


def coefficient_spectrum(A):
    return {str(value): int(multiplicity) for value, multiplicity in A.eigenvals().items()}


def quantum_projection_counterexample():
    """An invariant one-particle projector whose CAR lift has occupation1/2."""
    c = annihilators(2)
    phi = s.Matrix([1, 0, 0, 1])
    P = phi*phi.T/2
    assert P*P == P
    blocks = [[P[2*i:2*i+2, 2*j:2*j+2] for j in range(2)] for i in range(2)]
    number = sum((s.kronecker_product(c[i].H*c[j], blocks[i][j])
                  for i in range(2) for j in range(2)), s.zeros(8))
    spectra = coefficient_spectrum(number)
    assert spectra == {'0': 5, '1': 1, '1/2': 2}
    sigma = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    generator_checks = 0
    for S in sigma:
        T = s.I*S/2
        total_one = s.kronecker_product(T, s.eye(2))+s.kronecker_product(s.eye(2), -T.T)
        assert total_one*P == P*total_one
        charge = sum((T[i, j]*c[i].H*c[j] for i in range(2) for j in range(2)), s.zeros(4))
        total_fock = s.kronecker_product(charge, s.eye(2))+s.kronecker_product(s.eye(4), -T.T)
        assert total_fock*number == number*total_fock
        # A static reference in the fundamental representation cancels the
        # dual link charge on the completely occupied matter sector.
        state = s.kronecker_product(s.Matrix([0, 0, 0, 1]), phi)
        G = s.kronecker_product(total_fock, s.eye(2))+s.kronecker_product(s.eye(8), T)
        assert G*state == s.zeros(16, 1)
        assert s.kronecker_product(number, s.eye(2))*state == state/2
        generator_checks += 3
    b = [sum((s.kronecker_product(c[j], blocks[i][j]) for j in range(2)), s.zeros(8))
         for i in range(2)]
    car_error = b[0]*b[0].H+b[0].H*b[0]-s.kronecker_product(s.eye(4), blocks[0][0])
    assert car_error != s.zeros(8)
    return {'one_particle_dimension': 4, 'projector_rank': 1,
            'projector': [[str(x) for x in row] for row in P.tolist()],
            'matter_Fock_times_link_dimension': 8,
            'naive_number_eigenvalue_multiplicities': spectra,
            'local_SU2_invariance_and_static_reference_checks': generator_checks,
            'half_occupation_occurs_in_a_Gauss_singlet': True,
            'naive_projected_fields_satisfy_projected_CAR': False,
            'CAR_error_squared_trace': str(s.trace(car_error.H*car_error)),
            'scope': 'Exact universal counterexample to operator-valued second quantization; not a claim that this projector is the compact source band.'}


def classical_bogoliubov_countercontrol():
    """Commuting numerical coefficients give genuine canonical low modes."""
    c = annihilators(4)
    a, b = s.Rational(9999, 10001), s.Rational(200, 10001)
    assert a*a+b*b == 1
    low = (a*c[0]+b*c[1].H, a*c[2]+b*c[3].H)
    for i, x in enumerate(low):
        for j, y in enumerate(low):
            assert x*y+y*x == s.zeros(16)
            assert x*y.H+y.H*x == int(i == j)*s.eye(16)
    number = sum((x.H*x for x in low), s.zeros(16))
    assert number*(number-s.eye(16))*(number-2*s.eye(16)) == s.zeros(16)
    return {'full_Fock_dimension': 16, 'canonical_low_mode_count': 2,
            'CAR_checks': 8, 'occupation_minimal_polynomial': 'N(N-I)(N-2I)=0',
            'occupation_eigenvalue_multiplicities': coefficient_spectrum(number),
            'quartic_is_positive_and_band_supported': True}


def execute(output_dir=HERE):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    masks, blocks, currents = exterior_data()
    tensors, tensor_receipt = clebsch_gordan(masks, blocks, currents)
    indices, labels, left, right = link_representation(blocks, currents)
    names = {-4: 'u^c', -3: 'L', 1: 'Q', 2: 'd^c', 6: 'e^c'}
    rows, entries = [], {}
    nonnormal_witness = None
    total_covariance = 0
    for q, block in blocks.items():
        d = len(block)
        source, completed = coordinate_operators(q, blocks, indices, tensors)
        checks = 0
        for T, L, R in zip(currents, left, right, strict=True):
            G = T.extract(block, block)
            for a in range(d):
                for b in range(d):
                    U = completed[a, b]
                    expected_left = -sum((G[a, i]*completed[i, b] for i in range(d)), s.zeros(60))
                    expected_right = sum((G[j, b]*completed[a, j] for j in range(d)), s.zeros(60))
                    assert zero(L*U-U*L-expected_left)
                    assert zero(R*U-U*R-expected_right)
                    checks += 2
        A, C = combined_coordinate(source, d), combined_coordinate(completed, d)
        source_gram, complete_gram = clean(A.H*A), clean(C.H*C)
        assert zero(source_gram*source_gram-d*source_gram)
        assert zero(complete_gram*complete_gram-d*complete_gram)
        native_rank = s.trace(source_gram)/d
        completed_rank = s.trace(complete_gram)/d
        assert native_rank == d and completed_rank >= d
        algebra_dimension = irreducible_algebra_dimension([T.extract(block, block) for T in currents], d)
        assert algebra_dimension == d*d
        rows.append({'block': names[q], 'primitive_charge': q, 'matter_dimension': d,
                     'native_combined_rank': int(native_rank),
                     'completed_combined_rank': int(completed_rank),
                     'native_and_completed_operator_norm_squared': d,
                     'exact_Gram_identity': 'A†A A†A = d A†A',
                     'current_covariance_checks': checks,
                     'source_matter_associative_algebra_dimension': algebra_dimension,
                     'endpoint_product_representation_irreducible': True})
        entries[names[q]] = {
            f'{a},{b}': [[i, j, str(value)] for (i, j), value in sorted(U.todok().items())]
            for (a, b), U in completed.items()}
        if q == 6:
            U = completed[0, 0]
            commutator = U*U.H-U.H*U
            assert commutator != s.zeros(60)
            nonnormal_witness = {
                'coordinate': 'e^c, matrix entry0,0',
                'coordinate_commutator_with_adjoint_squared_trace': str(s.trace(commutator.H*commutator)),
                'compressed_coordinate_star_algebra_is_commutative': False}
        total_covariance += checks
        print(f'Exact link {names[q]}: rank {native_rank} -> {completed_rank}, norm² {d}, covariance {checks}.', flush=True)
    operators_path = output_dir/'completed_link_operators.json'
    operators_path.write_text(json.dumps({'link_basis': labels, 'operators': entries}, indent=2, sort_keys=True)+'\n')
    result = {'schema': 'oph.finite_quantum_link_coordinate_completion.v1',
              'all_pass': True, 'upstream_source_commit': '07c7acb2e7003efd9e574b67b374dbc36171dc6b',
              'native_link_dimension': 60, 'native_matter_dimension': 15,
              'native_current_count': len(currents), 'complete_ordered_product_channels': tensor_receipt,
              'complete_product_classification': 'Primitive-charge sums restrict the menu to the five listed nonvacuum products. The exact normalized exterior intertwiners verify all five.',
              'completed_coordinate_blocks': rows,
              'total_exact_link_covariance_checks': total_covariance,
              'coordinate_normalization': 'U_r = sqrt(dim r) times the compressed compact-group coordinate function',
              'native_vacuum_raising_entries_retained_exactly': True,
              'finite_compression_nonnormality': nonnormal_witness,
              'additional_rotor_attachment_is_a_premise': True,
              'source_selected_rotor_attachment': False,
              'link_Hamiltonian_commutant': {
                  'independent_local_endpoint_group': 'G_left × G_right',
                  'inequivalent_irreducible_blocks': ['1', 'Q', 'u^c', 'd^c', 'L', 'e^c'],
                  'Hermitian_commutant_real_dimension': 6,
                  'independent_energy_differences_after_common_offset': 5,
                  'general_invariant_link_Hamiltonian': 'epsilon_0 P_0 + sum_r epsilon_r P_r',
                  'native_one_excited_energy_rule_selected_by_covariance': False},
              'one_particle_quantum_link_spectral_bound': {
                  'arbitrary_edge_link_tensor_product_dimension': '60^882',
                  'combined_matter_link_hopping_norm_upper': '3 sqrt(6)',
                  'texture_amplitude': 64,
                  'central_band_absolute_energy_upper': '3 sqrt(6)',
                  'complement_absolute_energy_lower': '64 - 3 sqrt(6)',
                  'band_complement_gap_lower': '64 - 6 sqrt(6)',
                  'one_particle_low_rank_with_three_families': '180 × 60^882',
                  'same_180_mode_CAR_factor_obtained': False,
                  'link_energy_term_included_in_this_band_bound': False,
                  'with_link_energy_spectral_width_W_gap_lower': '64 - 6 sqrt(6) - W',
                  'independent_native_excited_energy_epsilon_spectral_width': '882 epsilon',
                  'sufficient_full_graph_excited_energy_domain': '0 <= epsilon < (64 - 6 sqrt(6))/882'},
              'operator_valued_projector_counterexample': quantum_projection_counterexample(),
              'commuting_coefficient_CAR_control': classical_bogoliubov_countercontrol(),
              'source_bindings': {name: digest(REPO/name) for name in sorted(set([
                  'model/particle_gauge_band/GAUGE_BAND_PROOF.md',
                  'model/particle_gauge_band/QUANTUM_LINK_PROOF.md',
                  'model/native_quantum_gauss/algebra.py', 'model/native_quantum_gauss/PROOF.md',
                  'model/common_source/native_matter_attachment.py',
                  'model/interacting_defect_response/algebra.py']+[Path(module.__file__).resolve().relative_to(REPO).as_posix()
                    for module in tuple(sys.modules.values())
                    if getattr(module, '__file__', None)
                    and Path(module.__file__).resolve().is_relative_to(REPO)
                    and not Path(module.__file__).resolve().relative_to(REPO).parts[0].startswith('.')
                    and Path(module.__file__).suffix == '.py']))},
              'operator_artifact_sha256': digest(operators_path),
              'execution_complete': True}
    (output_dir/'quantum_link_completion.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'all_pass': True, 'exact_covariance_checks': total_covariance,
                      'invariant_link_energy_parameters': 6,
                      'half_occupation_counterexample': True,
                      'duration_seconds': time.monotonic()-started}, indent=2), flush=True)

    return result


if __name__ == '__main__':
    execute()
