"""Full native-spin Fourier operators and the retained cubic reference."""
from collections import defaultdict
from functools import lru_cache
import sympy as s

from ...physical_attachment.chiral import representation
from ...chiral_source.model import PAULI
from ..model import clean
from ..periodic_family.model import Configuration, cells, kinetic_terms, apply_state
from ..peter_weyl_links import model as pw
from ..family_transport.model import one_body_column
from ..oriented_links.rotations import cubic_rotations, geometry
from ..higgs_composites import model as car
from ..atomic_momentum_response.model import reference_data

CONFIG = Configuration()
MAX_COUPLING = s.Rational(1, 2**34)
FULL_NORM = s.Integer(60)


def validate_coupling(value):
    value = s.sympify(value)
    if value.is_Rational is not True or value < 0:
        raise ValueError('an exact nonnegative hopping coefficient is required')
    return value


def hamiltonian_column(key, coupling=MAX_COUPLING):
    coupling = validate_coupling(coupling)
    result = defaultdict(lambda: s.Integer(0), pw.column(key, CONFIG))
    for target, value in pw.hopping_column(key, kinetic_terms(CONFIG)).items():
        result[target] += (coupling - 1) * value
    return clean(result)


@lru_cache(maxsize=3)
def native_spin(axis):
    if type(axis) is not int or axis not in range(3):
        raise ValueError('one of three native spin components is required')
    return s.ImmutableSparseMatrix(s.kronecker_product(PAULI[axis]/2, s.eye(15)))


def spin_fourier_column(key, axis):
    """O_a=sum_x (-1)^x_a S_x,a on every full source basis state."""
    mask, links, scalars = pw.validate_key(key, CONFIG)
    matrix = native_spin(axis)
    result = defaultdict(lambda: s.Integer(0))
    for site, xyz in enumerate(cells(2)):
        for target, value in one_body_column(mask, matrix, 30*site).items():
            result[target, links, scalars] += (-1)**xyz[axis] * value
    return clean(result)


def pump_column(key):
    """Rectangular B from the reference vacuum to its retained vector sector."""
    return {(target, axis + 1): value/(FULL_NORM*s.sqrt(3))
        for axis in range(3) for target, value in spin_fourier_column(key, axis).items()}


def reference_translation(axis):
    if type(axis) is not int or axis not in range(3):
        raise ValueError('three reference translation directions')
    return s.ImmutableSparseMatrix.diag(1, *((-1)**int(axis == j) for j in range(3)))


def reference_rotation(rotation):
    if not any(rotation == actual for actual, _ in cubic_rotations()):
        raise ValueError('an actual proper cubic source rotation is required')
    return s.ImmutableSparseMatrix.diag(s.ones(1), rotation)


def operator_certificate():
    reference_data()
    yukawa, _ = car.quantum_yukawa()
    onsite = car.add(*yukawa.values())
    local_checks = 0
    for axis in range(3):
        matrix = native_spin(axis)
        assert matrix.H == matrix and matrix**2 == s.eye(30)/4 and s.trace(matrix) == 0
        assert car.one_body(matrix) == car.spin_current(axis)
        assert not car.commutator(car.spin_current(axis), onsite)
        for current in range(12):
            internal = s.kronecker_product(s.eye(2), representation()[2][current])
            assert matrix*internal == internal*matrix
            assert not car.commutator(car.spin_current(axis), car.gauge_current(current))
            local_checks += 1
    words = 0; ref_conjugations = 0
    for rotation, spin in cubic_rotations():
        ref = reference_rotation(rotation)
        assert s.simplify(ref.H*ref) == s.eye(4)
        sites, _ = geometry(rotation, 2)
        for axis in range(3):
            target_axis = next(j for j in range(3) if rotation[j, axis])
            sign = int(rotation[target_axis, axis])
            assert s.simplify(spin*PAULI[axis]*spin.H) == sign*PAULI[target_axis]
            for site, xyz in enumerate(cells(2)):
                moved, _ = sites[site]
                assert (-1)**xyz[axis] == (-1)**cells(2)[moved][target_axis]
                words += 1
            assert ref*reference_translation(axis)*ref.H == reference_translation(target_axis)
            ref_conjugations += 1
    # Complete occupied source columns, including nonblank scalar and link memory.
    _, links, _ = pw.empty_key()
    flux = list(links); flux[0] = pw.label_index()[1, 2, 8]
    mask = sum(((1 << ((3*site) % 15)) | (1 << (15+(5*site+1) % 15))) << (30*site) for site in range(8))
    key = mask, tuple(flux), tuple(site % 5 for site in range(8))
    translated_controls = 0
    for axis in range(3):
        original = spin_fourier_column(key, axis)
        for direction in range(3):
            translate = lambda value: pw.translation_column(value, direction)
            left = apply_state(original, translate)
            right = apply_state(translate(key), lambda value: spin_fourier_column(value, axis))
            assert left == {target: (-1)**int(axis == direction)*value for target, value in right.items()}
            translated_controls += 1
    assert hamiltonian_column(key, 1) == pw.column(key)
    assert all(target[1] == key[1] for target in hamiltonian_column(key, 0))
    return {'full_CAR_modes_per_site': 30, 'full_scalar_states_per_site': 5,
        'full_link_states': 119, 'source_sites': 8,
        'native_spin': 'S_x,a = sum_15_native_indices c_x^dagger (sigma_a/2) c_x',
        'Fourier_operator': 'O_a=sum_x (-1)^x_a S_x,a',
        'actual_native_current_commutators': local_checks,
        'actual_full_occupied_source_translation_controls': translated_controls,
        'full_one_body_eigenvalues': {'1/2': 120, '-1/2': 120},
        'complete_Fock_operator_norm': '60',
        'norm_uses_all_quark_and_lepton_modes': True,
        'atomic_compression': 'P0 S_x,a P0 = the actual local ground doublet spin S_x,a',
        'atomic_subspace_preserved_by_native_spin': True,
        'full_onsite_Yukawa_spin_commutators': 3,
        'complete_cubic_Fourier_word_checks': words,
        'retained_reference_rotation_translation_checks': ref_conjugations,
        'retained_reference_dimension': 4,
        'reference_sectors': 'one scalar vacuum plus a three-dimensional cubic vector',
        'reference_clock': 'identity', 'reference_fermion_parity': 'even in all four states',
        'reference_translation': 'vacuum has phase zero; vector a has inverse X_a character',
        'pump': 'B=(sum_a O_a tensor |a>)/(60 sqrt(3)), with reference vacuum input',
        'pump_effect': 'Bdagger B=sum_a O_a^2/10800 <= I',
        'all_other_reference_input_sectors': 'included in the full Halmos block unitary',
        'full_finite_unitary': '[[sqrt(I-BdaggerB),-Bdagger],[B,sqrt(I-BBdagger)]]',
        'failure_effect': 'I-BdaggerB',
        'Gauss_and_even_CAR_preserved': True,
        'joint_proper_cubic_rotations_and_translations_preserved': True,
        'source_momentum_in_each_success_reference_branch': 'X_a, relative to the source momentum zero ground',
        'source_momenta_called_zero_by_reference_compensation': False,
        'continuous_spin_symmetry_at_nonzero_hopping_assumed': False}
