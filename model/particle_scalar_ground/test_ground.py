"""Complete-source cell maps, exact ground admission and observable records."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse
import sympy as s

from model.particle_scalar_ground import ground as g


def test_complete_local_signed_frames_and_car_pair_maps():
    for spin_sign in (1, -1):
        masks = []
        for kind in range(4):
            size = 3 if kind in (1, 2) else 1
            for coordinate in range(size):
                mask, phase = g.local_frame(kind, coordinate, spin_sign)
                masks.append(mask)
                for scalar in (1, 2):
                    first = g.source._create(mask, 2)
                    second = None if first is None else g.source._create(first[0], scalar-1)
                    actual = {} if second is None else {second[0]: phase*spin_sign*first[1]*second[1]}
                    possible = kind in (1, 2) and coordinate == (0 if kind == 1 else scalar)
                    if possible:
                        dest = scalar if kind == 1 else 0
                        target_mask, target_phase = g.local_frame(kind, dest, spin_sign)
                        assert actual == {target_mask: target_phase}
                    else:
                        assert actual == {}
        assert sorted(masks) == list(range(8))
    assert sum(np.prod([3 if k in (1, 2) else 1 for k in kinds])
               for kinds in product(range(4), repeat=6)) == 1 << 18


def test_signed_full_car_embedding_against_literal_source_columns():
    masks = list(range(64)) + [((1103515245*i+12345) % (1 << 18)) for i in range(96)]
    cache = {}
    for mask in masks:
        kinds, coordinate, phase = g.decoded_column(mask)
        assert g.original_column(kinds, coordinate) == (mask, phase)
        if kinds not in cache:
            states, operators = g.pair_operators(kinds)
            cache[kinds] = states, operators, {state: i for i, state in enumerate(states)}
        states, operators, indices = cache[kinds]
        expected = {}
        for scalar, operator in enumerate(operators, start=1):
            column = operator[:, indices[coordinate]].tocoo()
            for row, weight in zip(column.row, column.data, strict=True):
                target, target_phase = g.original_column(kinds, states[int(row)])
                expected[scalar, target] = int(weight)*phase*target_phase
        actual = g.source.vertex_action(s.diag(1, 2, 3), 0, mask)
        assert actual == expected
        # Literal adjoint columns retain their complete zero/occupied rules.
        for (scalar, destination), coefficient in actual.items():
            assert g.source.vertex_action(s.diag(1, 2, 3), scalar, destination).get((0, mask)) == coefficient


def test_inactive_cell_compression_has_the_exact_positive_increment():
    for inactive in (0, 3):
        kinds = (inactive, 2, 1, 2, 1, 2)
        old_states, _, old = g.gram(kinds)
        active = (1,)+kinds[1:]
        states, _, enlarged = g.gram(active)
        index = {state: i for i, state in enumerate(states)}
        image = [index[(0,)+state[1:]] for state in old_states]
        compression = enlarged[image][:, image]
        assert (compression-old-2*sparse.eye(len(old_states), dtype=np.int64, format='csr')).nnz == 0


def test_exact_perron_bounds_reject_bad_domain_and_admit_no_float_premise():
    matrix = sparse.csr_matrix([[2, 1], [1, 3]], dtype=np.int64)
    assert g.check_perron_vector(matrix, (1, 1)) == (F(3), F(4))
    for vector in ((0, 1), (1.0, 2), (True, 1), (1,)):
        with pytest.raises(ValueError):
            g.check_perron_vector(matrix, vector)
    with pytest.raises(ValueError):
        g.check_perron_vector(matrix.astype(float), (1, 1))
    with pytest.raises(ValueError):
        g.check_perron_vector(sparse.csr_matrix([[1, -1], [-1, 2]], dtype=np.int64), (1, 1))


def test_complete_active_sector_comparison_and_exact_algebraic_global_ground():
    data = g.all_active_blocks(); root = g.exact_root_and_gap()
    assert len(data['sectors']) == 27
    assert sum(row['represented_active_type_patterns'] for row in data['sectors']) == 64
    assert all(row['complete_source_sector_dimension'] == 729 for row in data['sectors'])
    assert data['winning']['block'].shape == (53, 53)
    assert data['winning']['lower'] > data['other_component_upper']
    assert root['global_ground_multiplicity'] == 1
    assert root['global_ground_interval'][1] < -7
    assert root['full_ground_gap_lower'] > F(469, 10000)
    assert root['full_source_dimension'] == 786432
    assert root['polynomial'] == (1, -244, 23737, -1171722, 30638880, -394771968, 1881169920)


def test_global_ground_actual_gauge_charge_and_exact_singlet_trial():
    numbers = g.quantum_numbers(); trial = g.singlet_trial()
    assert numbers['scalar_empty_support_NL_Ne_NF_twice_hypercharge'] == [(4, 2, 6, 0)]
    assert numbers['scalar_occupied_support_NL_Ne_Qstar_twice_hypercharge'] == [(5, 3, 6, 0)]
    weak = tuple((6*f+spin, 6*f+2+spin) for f in range(3) for spin in range(2))
    spin = tuple((6*f+2*internal, 6*f+2*internal+1) for f in range(3) for internal in range(3))
    assert g.one_body_action(trial['source_state'], weak) == {}
    assert g.one_body_action(trial['source_state'], spin) == {}
    assert all(mask.bit_count() == 6 for mask in trial['source_state'])
    assert all(mask.bit_count() == 8 for _, mask in trial['source_image'])
    assert len(trial['source_state']) == 53
    for mask in list(trial['source_state'])[:5]:
        assert g.scalar_image({mask: F(1)}) == g.source.vertex_action(s.diag(1, 2, 3), 0, mask)


def test_complete_preparation_bound_and_sector_choice_are_explicit():
    report = g.prepared_ground_witness()
    assert report['ground_infidelity_upper'] < F(26, 10**12)
    assert report['approximate_preparation_is_exactly_gauge_and_spin_singlet']
    assert not report['closed_source_vertex_from_empty_Qstar0_prepares_Qstar6']
    assert not report['compiled_primitive_tape_executed']
    # The full Householder map needs only the normalized difference vector;
    # its complement is retained. Check the universal two-vector algebra.
    overlap = s.Rational(3, 5)
    seed = s.Matrix([1, 0, 0]); target = s.Matrix([overlap, s.Rational(4, 5), 0])
    w = (seed-target)/s.sqrt(2*(1-overlap)); unitary = s.eye(3)-2*w*w.H
    assert unitary.H*unitary == s.eye(3)
    assert unitary*seed == target
    assert unitary*s.Matrix([0, 0, 1]) == s.Matrix([0, 0, 1])


def test_same_h_scalar_joint_law_on_exact_invariant_singular_pair():
    lam = s.symbols('lam', positive=True); delta = s.sqrt(1+4*lam)
    alpha, beta = (1+1/delta)/2, (1-1/delta)/2
    h = s.Matrix([[0, s.sqrt(lam)], [s.sqrt(lam), 1]])
    energy = (1-delta)/2
    # An unnormalized exact eigenvector avoids radicals of radicals.
    psi = s.Matrix([s.sqrt(lam), energy]); norm = (psi.H*psi)[0]
    assert (h*psi-energy*psi).applyfunc(s.simplify) == s.zeros(2, 1)
    assert s.simplify(energy**2/norm-beta) == 0
    assert s.simplify(alpha*beta-lam/delta**2) == 0
    t = s.symbols('t', real=True)
    transition = 4*lam/delta**2*s.sin(delta*t/2)**2
    probs = (alpha*(1-transition), alpha*transition, beta*transition, beta*(1-transition))
    assert s.simplify(sum(probs)-1) == 0
    assert s.trigsimp(s.sin(delta*t)**2/(4*s.sin(delta*t/2)**2)-s.cos(delta*t/2)**2) == 0


def test_record_enclosures_keep_full_domain_and_certified_preparation_error():
    report = g.complete_ground_records()
    assert report['pointer_rank_empty']+report['pointer_rank_occupied'] == 786432
    assert report['pure_preparation_trace_distance_upper'] < F(51, 10**7)
    assert report['example_complete_law_TV_error_upper'] < F(101, 10**7)
    for tick in report['exact_record_enclosures']:
        bounds = tick['unconditional_joint_record_intervals']
        assert set(bounds) == {'00', '01', '10', '11'}
        assert all(0 < F(lo) <= F(hi) < 1 for lo, hi in bounds.values())
        assert sum(F(row[0]) for row in bounds.values()) <= 1 <= sum(F(row[1]) for row in bounds.values())
    assert not report['record_gap_identified_as_relativistic_particle_mass']


def test_general_family_rotation_preserves_literal_pair_matrix():
    # Exact family SVD convention: left creators U and right creators Vbar.
    u = s.Matrix([[s.Rational(3, 5), s.Rational(4, 5), 0], [-s.Rational(4, 5), s.Rational(3, 5), 0], [0, 0, s.I]])
    v = s.Matrix([[0, s.I, 0], [1, 0, 0], [0, 0, 1]])
    sigma = s.diag(1, 2, 3); y = u*sigma*v.H
    assert u.H*u == v.H*v == s.eye(3)
    reconstructed = sum((sigma[j, j]*u[:, j]*v[:, j].conjugate().T for j in range(3)), s.zeros(3))
    assert reconstructed == y
    assert s.trace(y.H*y) == 14


def test_source_paths_resolve_and_no_measured_inputs_are_needed():
    root = Path(__file__).resolve().parents[2]
    assert all((root/path).is_file() for path in g.SOURCE_INPUTS)
    assert g.WEIGHTS == (1, 1, 2, 2, 3, 3)
