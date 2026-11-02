"""Independent exact finite controls for the centered response bridge."""
from functools import lru_cache

import numpy as np
import pytest
from scipy import sparse
import sympy as s

from ..field_stress.energy import GaussianInteger as GI
from ..cube_eta_protection import cube as source, repair
from . import dimer, particles


@pytest.fixture(scope='module')
def proved():
    return particles.certificate()


def test_complete_source_sector_current_and_record_counts(proved):
    selection = proved['selection']
    assert selection['exact_CAR_tensor_identities'] == 96
    assert selection['joint_singlet_isotypic_dimension'] == 588
    assert selection['charged_joint_irrep_multiplicity'] == 1344
    assert selection['charged_joint_isotypic_dimension'] == 5376
    response = proved['response_and_current']
    assert response['oriented_bond_currents_checked'] == 12
    assert response['local_charge_and_continuity_equations_checked'] == 8
    assert response['physical_site_conjugation_signs'] == [1, 1, -1, -1, 1, 1, -1, -1]
    complete = proved['complete_repair']
    assert complete['joint_input_dimension'] == 22880
    assert complete['native_raw_slots_retained'] == 120
    assert complete['outside_slots_retained'] == 10
    assert complete['nontrivial_native_defect_record_obtained'] is False


def test_no_gap_preparation_mass_or_full_reachability_is_inferred(proved):
    assert proved['selection']['every_isotypic_copy_reachable_claimed'] is False
    assert proved['selection']['numerical_eigenvalue_or_residual_used_for_admission'] is False
    assert proved['selection']['exact_ground_certificate_used'] is True
    assert proved['complete_repair']['new_native_preparation_or_detector_executed'] is False
    assert proved['new_gap_or_particle_mass_claimed'] is False
    assert proved['physical_species_identified'] is False
    mass = proved['onsite_mass_order']
    assert mass['actual_unique_ground_order_expectations'] == ['0']*4
    assert mass['interacting_gap_or_thermodynamic_symmetry_breaking_excluded'] is False


def test_bit_complement_uses_the_actual_fermionic_sign_on_every_column():
    complement, phase = particles.particle_hole_data()
    assert np.array_equal(complement[complement], np.arange(65536))
    for array in (complement, phase):
        with pytest.raises(ValueError): array.setflags(write=True)
    # Independently evaluate C d_j C on each allowed occupation column.
    for mode in range(16):
        eta = (-1)**sum(source.GEOMETRY.sites[mode//2])
        for mask in range(65536):
            if mask >> mode & 1: continue
            opposite = int(complement[mask])
            after = opposite^(1 << mode)
            left = int(phase[mask])*((-1)**(opposite&((1 << mode)-1)).bit_count())*int(phase[after])
            right = eta*((-1)**(mask&((1 << mode)-1)).bit_count())
            assert int(complement[after]) == mask^(1 << mode) and left == right


def test_public_array_metadata_cannot_modify_the_cached_particle_hole_map(proved):
    complement, phase = particles.particle_hole_data()
    complement.dtype = np.uint8
    phase.shape = (256, 256)
    fresh_complement, fresh_phase = particles.particle_hole_data()
    assert fresh_complement.shape == fresh_phase.shape == (65536,)
    assert fresh_complement.dtype == fresh_phase.dtype == np.dtype(np.int64)
    assert np.array_equal(fresh_complement[fresh_complement], np.arange(65536))
    assert particles.response_certificate() == proved['response_and_current']


@lru_cache(maxsize=1)
def _small_response():
    data = dimer.data(); fields = [dimer.annihilator(j) for j in range(4)]
    g = data['ground']; norm = data['ground_norm_squared']
    # A complex onsite spin-frame control, unrelated to the cube frame table.
    u = s.Matrix([[s.Rational(3, 5), 4*s.I/5], [4*s.I/5, s.Rational(3, 5)]])
    r = s.diag(s.eye(2), u)
    physical = [sum((r.H[i, j]*fields[j] for j in range(4)), s.zeros(16)) for i in range(4)]
    f = r.H*s.diag(1, 1, -1, -1)*r.conjugate()
    add = s.Matrix.hstack(*(a.H*g for a in physical))
    remove = s.Matrix.hstack(*(a*g for a in physical))
    return data, g, norm, f, add, remove


def test_complete_complex_matrix_response_not_just_equal_diagonal_weights():
    data, _, norm, f, add, remove = _small_response()
    # All N=1/3 energies have two distinct values, so moments0,1 determine
    # their projectors; higher moments supply independent complete controls.
    for power in range(4):
        polynomial = data['K']**power
        a = s.simplify(add.H*polynomial*add/norm)
        b = s.simplify(remove.H*polynomial*remove/norm)
        assert s.simplify(a-f*b*f.H) == s.zeros(4)
    a = s.simplify(add.H*data['K']*add/norm)
    b = s.simplify(remove.H*data['K']*remove/norm)
    naive = s.diag(1, 1, -1, -1)
    assert s.simplify(a-naive*b*naive) != s.zeros(4)
    assert s.simplify(a-f*b.T*f.H) != s.zeros(4)


def test_independent_local_casimirs_select_both_charge_branches():
    data = dimer.data(); g = data['ground']; fields = [dimer.annihilator(j) for j in range(4)]
    plus = sum((fields[2*x].T*fields[2*x+1] for x in range(2)), s.zeros(16))
    sz = sum((fields[2*x].T*fields[2*x]-fields[2*x+1].T*fields[2*x+1] for x in range(2)), s.zeros(16))/2
    ep = fields[0].T*fields[1].T-fields[2].T*fields[3].T
    ez = (data['N']-2*s.eye(16))/2
    spin2 = sz*sz+(plus*plus.T+plus.T*plus)/2
    eta2 = ez*ez+(ep*ep.T+ep.T*ep)/2
    assert s.simplify(spin2*g) == s.zeros(16, 1)
    assert s.simplify(eta2*g) == s.zeros(16, 1)
    for field in fields:
        for operator in (field, field.T):
            vector = operator*g
            for casimir in (spin2, eta2):
                assert s.simplify(casimir*vector-3*vector/4) == s.zeros(16, 1)


def test_wrong_charge_clock_is_rejected_before_response_claim(monkeypatch):
    monkeypatch.setattr(particles, '_centered_action', lambda: source.operators().hubbard_H)
    with pytest.raises(ArithmeticError, match='centered clock'):
        particles.response_certificate()


def test_probability_preserving_wrong_particle_hole_phase_is_rejected(monkeypatch):
    complement, phase = particles.particle_hole_data()
    wrong = phase.copy(); wrong[0] *= -1; wrong[-1] *= -1
    assert np.array_equal(wrong[complement]*wrong, np.ones(65536, dtype=np.int64))
    monkeypatch.setattr(particles, 'particle_hole_data', lambda: (complement, wrong))
    with pytest.raises(ArithmeticError, match='CAR map'):
        particles.response_certificate()


def test_repair_support_leakage_cannot_be_hidden_by_record_counts(monkeypatch):
    a, b = repair.exact_embedding()
    changed = a.tolil(); changed[int(source.weights(7)[0]), 0] = 1; changed = changed.tocsr()
    monkeypatch.setattr(repair, 'exact_embedding', lambda: (changed, b))
    with pytest.raises(ArithmeticError, match='N6 sector'):
        particles.repair_certificate()


def test_every_zero_raw_slot_is_required(monkeypatch):
    original = repair.matrices
    monkeypatch.setattr(repair, 'matrices', lambda law, seam: original(law, seam)[:-1])
    with pytest.raises(ArithmeticError, match='raw slot'):
        particles.repair_certificate()


@pytest.mark.parametrize('corruption', ['int64_min', 'masked'])
def test_zero_slot_overflow_or_missingness_cannot_forge_complete_repair(monkeypatch, corruption):
    original = repair.matrices
    def changed(law, seam):
        slots = list(original(law, seam))
        index = next(j for j, (_, matrix) in enumerate(slots) if not np.any(matrix))
        label, matrix = slots[index]
        matrix = matrix.copy()
        if corruption == 'int64_min': matrix[0, 0] = np.iinfo(np.int64).min
        else: matrix = np.ma.array(matrix, mask=True)
        slots[index] = label, matrix
        return tuple(slots)
    monkeypatch.setattr(repair, 'matrices', changed)
    with pytest.raises(ArithmeticError, match='partial injection'):
        particles.repair_certificate()


def test_narrow_integer_complete_effect_sum_cannot_wrap(monkeypatch):
    original = repair.matrices
    def changed(law, seam):
        slots = original(law, seam)
        if law == 'C': return slots
        stacked = np.zeros((384, 24), dtype=np.int8)
        stacked[:258, 0] = 1
        for j in range(1, 24): stacked[258+2*(j-1):260+2*(j-1), j] = 1
        matrices = tuple(stacked[24*j:24*(j+1)] for j in range(16))
        assert np.array_equal(sum(m.T@m for m in matrices), 2*np.eye(24))
        assert int((stacked.astype(np.int64).T@stacked.astype(np.int64))[0, 0]) == 258
        return tuple((label, matrix) for (label, _), matrix in zip(slots, matrices, strict=True))
    monkeypatch.setattr(repair, 'matrices', changed)
    with pytest.raises(ArithmeticError, match='partial injection'):
        particles.repair_certificate()


def test_actual_complete_repair_on_complex_charge_reference_columns():
    state = np.zeros((65536, 2), complex)
    state[source.weights(7)[[0, 17, 101]], 0] = [1+2j, -3j, 2]
    state[source.weights(9)[[7, 31, 401]], 1] = [2-1j, 3, 1j]
    state /= np.linalg.norm(state)
    count = 0
    for law in ('C', 'E'):
        for seam in repair.SEAMS:
            outputs = repair.apply(state, law, seam)
            assert len(outputs) == (9 if law == 'C' else 17)
            assert outputs[-1][0] == ('outside',)
            for _, vector in outputs[:-1]: assert np.array_equal(vector, np.zeros_like(state))
            assert np.array_equal(outputs[-1][1], state)
            count += len(outputs)
    assert count == 130


def test_zero_mass_order_is_not_zero_fluctuation_or_zero_charged_gap():
    data = dimer.data(); g = data['ground']; norm = data['ground_norm_squared']
    fields = [dimer.annihilator(j) for j in range(4)]
    n = [a.T*a for a in fields]
    scalar = n[0]+n[1]-n[2]-n[3]
    neel = n[0]-n[1]-n[2]+n[3]
    for order in (scalar, neel):
        assert s.simplify((g.T*order*g)[0]) == 0
        assert s.simplify((g.T*order*order*g)[0]/norm) > 0
    assert all(row['energy'] > 0 for row in data['charged_lines'])


@pytest.mark.parametrize('bad', [None, np.eye(2), GI.from_exact(sparse.eye(2))])
def test_partial_or_untyped_particle_hole_operator_refused(bad):
    with pytest.raises(ValueError): particles.conjugate(bad)
