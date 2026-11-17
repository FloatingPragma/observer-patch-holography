#!/usr/bin/env python3
"""Exact spectral/CAR witness with the native 60-state raising link.

The two-site scalar texture is an explicit algebra witness, not the
343-site compact texture or a computed physical particle spectrum.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import sympy as s

from .quantum_link_completion import (HERE, REPO, annihilators, clean,
                                     digest, exterior_data,
                                     link_representation, zero)


def execute(output_dir=HERE):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    masks, blocks, currents = exterior_data()
    indices, labels, left, right = link_representation(blocks, currents)
    e = indices[6, 0, 0]
    U = s.SparseMatrix(60, 60, {(e, 0): 1})
    P0 = s.SparseMatrix(60, 60, {(0, 0): 1})
    Pe = s.SparseMatrix(60, 60, {(e, e): 1})
    eye = s.eye(60, cls=s.SparseMatrix)
    R = s.sqrt(4097)
    t, k = (1-64/R)/2, -s.Rational(1, 2)/R
    pb = [[eye-t*Pe, k*U], [k*U.H, t*P0]]
    P = s.SparseMatrix(s.BlockMatrix(pb))
    H = s.SparseMatrix(s.BlockMatrix([[s.zeros(60), U/2],
                                      [U.H/2, 64*eye]]))
    assert zero(P*P-P) and zero(P.H-P)
    assert zero(P*H-H*P) and s.trace(P) == 60
    # The only mixed two-dimensional block has energies (64+-sqrt(4097))/2.
    coupled = H.extract([e, 60], [e, 60])
    assert coupled == s.Matrix([[0, s.Rational(1, 2)],
                                 [s.Rational(1, 2), 64]])
    assert zero(P.extract([e, 60], [e, 60])-
                ((64+R)/2*s.eye(2)-coupled)/R)
    assert s.simplify(s.trace(P.extract(list(range(60)), list(range(60))))-(60-t)) == 0

    c = [s.SparseMatrix(x) for x in annihilators(2)]
    N = sum((s.kronecker_product(c[i].H*c[j], pb[i][j])
             for i in range(2) for j in range(2)), s.zeros(240, cls=s.SparseMatrix))
    N = s.SparseMatrix(N)
    assert zero(N.H-N)
    vacuum = list(range(60))
    one = list(range(60, 180))
    filled = list(range(180, 240))
    assert zero(N.extract(vacuum, vacuum))
    assert zero(N.extract(one, one)-P)
    assert zero(N.extract(filled, filled)-(eye+t*(P0-Pe)))
    psi = s.SparseMatrix(240, 1, {(180, 0): 1})
    assert zero(N*psi-(1+t)*psi)

    gauss = 0
    for T, L, A in zip(currents, left, right, strict=True):
        tr = T[blocks[6][0], blocks[6][0]]
        one_x = s.kronecker_product(s.diag(tr, 0), eye)+s.kronecker_product(s.eye(2), L)
        one_y = s.kronecker_product(s.diag(0, tr), eye)+s.kronecker_product(s.eye(2), A)
        assert zero(one_x*P-P*one_x) and zero(one_y*P-P*one_y)
        # One static dual e^c reference at each endpoint has scalar action -tr.
        fock_x = s.kronecker_product(tr*c[0].H*c[0]-tr*s.eye(4), eye)+s.kronecker_product(s.eye(4), L)
        fock_y = s.kronecker_product(tr*c[1].H*c[1]-tr*s.eye(4), eye)+s.kronecker_product(s.eye(4), A)
        assert zero(fock_x*N-N*fock_x) and zero(fock_y*N-N*fock_y)
        assert zero(fock_x*psi) and zero(fock_y*psi)
        gauss += 6
    result = {
        'schema': 'oph.native_link_spectral_CAR_witness.v1', 'all_pass': True,
        'upstream_source_commit': '07c7acb2e7003efd9e574b67b374dbc36171dc6b',
        'native_link_dimension': 60, 'native_link_block': 'e^c',
        'native_primitive_charge': 6, 'native_link_excited_index': e,
        'one_particle_dimension': 120, 'matter_Fock_times_link_dimension': 240,
        'site_energies': [0, 64], 'hopping': '1/2',
        'one_particle_spectrum': {'0': 59, '64': 59,
                                 '(64-sqrt(4097))/2': 1,
                                 '(64+sqrt(4097))/2': 1},
        'low_projector_rank': 60, 'low_complement_gap': '64',
        'projector_block_formula': [['I-t P_e', '-U/(2 sqrt(4097))'],
                                   ['-U^dagger/(2 sqrt(4097))', 't P_0']],
        't_exact': '(1-64/sqrt(4097))/2',
        't_decimal': str(s.N(t, 30)),
        'naive_number_spectrum': {'0': 120, '1': 118,
                                 '1+(1-64/sqrt(4097))/2': 1,
                                 '1-(1-64/sqrt(4097))/2': 1},
        'noninteger_occupation_Gauss_singlet': {
            'state': 'both matter sites occupied; blank link; one static dual e^c reference at each endpoint',
            'occupation': '1+(1-64/sqrt(4097))/2',
            'all_native_current_endpoint_checks': gauss,
            'includes_Gauss_reference_charges': True},
        'scope': {'native_60_state_raising_link_used_exactly': True,
                  'spectral_projector_of_explicit_two_site_Hamiltonian': True,
                  'noninteger_naive_Fock_occupation_proved_exactly': True,
                  'actual_compact_343_site_projector_computed': False,
                  'universal_obstruction_to_all_quantum_link_CAR_constructions': False,
                  'physical_particle_mass_calculation': False},
        'source_bindings': {name: digest(REPO/name) for name in (
            'model/particle_gauge_band/native_link_projector_witness.py',
            'model/particle_gauge_band/NATIVE_LINK_CAR_WITNESS.md',
            'model/particle_gauge_band/quantum_link_completion.py',
            'model/native_quantum_gauss/algebra.py',
            'model/interacting_defect_response/algebra.py')},
        'construction_dependency_sha256': digest(HERE/'quantum_link_completion.py'),
        'execution_complete': True}
    (output_dir/'native_link_projector_witness.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'all_pass': True, 'one_particle_dimension': 120,
                      'naive_number_has_noninteger_Gauss_singlet_eigenvalue': True,
                      't': str(s.N(t, 20)), 'exact_Gauss_checks': gauss,
                      'duration_seconds': time.monotonic()-started}, indent=2))

    return result


if __name__ == '__main__':
    execute()
