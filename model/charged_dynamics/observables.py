"""Actual charge, current and poles of the unchanged joined four-site field.

The source product-dimer input is not the coupled ground. A dressed spectral
quartet is invariant under this Hamiltonian, but its logical left/right labels
carry only the computed block-charge contrast. No native repair protection,
topological charge, physical electron assignment or common SI clock follows.
"""
from __future__ import annotations

import math
from numbers import Real

import numpy as np

from . import chain
from ..full_fock_field.field import FieldGeometry, FullFock
from ..particle_response import dimer
from ..reference_observation.observation import frozen


ADMISSION = 3e-10
J_BARE = (1+1/math.sqrt(17))/64
LOCAL_RESIDUE = (1+1/math.sqrt(17))/4
TAU_X = np.kron(np.array([[0., 1.], [1., 0.]]), np.eye(2))
TAU_Y = np.kron(np.array([[0., -1j], [1j, 0.]]), np.eye(2))
TAU_Z = np.diag([1., 1., -1., -1.])


def _axis(axis):
    if type(axis) is not int or axis not in (0, 1, 2):
        raise ValueError('integer coordinate axis in0..2 required')
    return axis


def _kind(kind):
    if type(kind) is not str or kind not in ('hole', 'particle'):
        raise ValueError('hole or particle charged doublet required')
    return -1 if kind == 'hole' else 1


def _duration(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError('finite Hamiltonian model duration in[-512,512] required')
    value = float(value)
    if not math.isfinite(value) or abs(value) > 512:
        raise ValueError('finite Hamiltonian model duration in[-512,512] required')
    return value


def _array(value, shape):
    if np.ma.isMaskedArray(value):
        raise ArithmeticError('unmasked complete field arrays required')
    value = np.asarray(value)
    if value.shape != shape or value.dtype.kind not in 'iufc' or not np.isfinite(value).all():
        raise ArithmeticError('finite complete field array shape required')
    return value


def _admit(*residuals):
    for value in residuals:
        residual = float(np.linalg.norm(value))
        if not math.isfinite(residual) or residual > ADMISSION:
            raise ArithmeticError('same-action observable identity failed numerical admission')


def operators(axis=0):
    """Fresh full-carrier operators, independently bound to the original action."""
    axis = _axis(axis)
    data = chain.setup(axis)
    space = FullFock(8)
    sides = tuple(4 if j == axis else 1 for j in range(3))
    number = np.diag([address.bit_count() for address in range(256)])
    left = np.diag([(address & 15).bit_count() for address in range(256)])
    right = number-left
    canonical = space.hamiltonian(FieldGeometry(sides)).toarray()-number/2
    local = dimer.setup(FieldGeometry(tuple(2 if j == axis else 1 for j in range(3))))
    disconnected = np.kron(local.generator, np.eye(16))+np.kron(np.eye(16), local.generator)
    k = _array(data.K, (256, 256))
    _admit(k-canonical, _array(data.K0, (256, 256))-disconnected,
           _array(data.V, (256, 256))-(canonical-disconnected),
           _array(data.N, (256, 256))-number, _array(data.NA, (256, 256))-left,
           _array(data.NB, (256, 256))-right)
    vectors = _array(data.vectors, (256, 256))
    energies = _array(data.energies, (256,))
    sectors = _array(data.sectors, (256,))
    if (energies.dtype.kind not in 'iuf' or sectors.dtype.kind not in 'iu'
            or any(not 0 <= int(n) <= 8 for n in sectors)):
        raise ArithmeticError('literal complete conserved-number spectrum required')
    _admit(vectors.conj().T@vectors-np.eye(256),
           k@vectors-vectors*energies, number@vectors-vectors*sectors)
    ground = _array(data.ground, (256,))
    _admit(np.vdot(ground, ground)-1, k@ground-data.ground_energy*ground,
           data.ground_energy-float(np.min(energies)), (number-4*np.eye(256))@ground,
           _array(data.ground_projector, (256, 256))-np.outer(ground, ground.conj()))
    if np.sort(energies)[1]-data.ground_energy <= 1e-8:
        raise ArithmeticError('unique coupled half-filled ground required')
    conjugation = np.kron(dimer.charge_conjugation(local), dimer.charge_conjugation(local))
    qa, qb = left-2*np.eye(256), right-2*np.eye(256)
    current = -1j*(k@left-left@k)
    _admit(conjugation.conj().T@conjugation-np.eye(256),
           conjugation@k@conjugation.conj().T-k,
           conjugation@qa@conjugation.conj().T+qa,
           conjugation@current@conjugation.conj().T+current,
           1j*(k@right-right@k)-current)
    return dict(data=data, K=frozen(k), N=frozen(number), QA=frozen(qa), QB=frozen(qb),
                current=frozen(current), conjugation=frozen(conjugation), space=space,
                local_reference=frozen(local.reference))


def _bare(kind, axis):
    """Bind both localized input fibers to the predecessor's actual projectors."""
    _kind(kind); _axis(axis)
    local = dimer.setup(FieldGeometry(tuple(2 if j == axis else 1 for j in range(3))))
    candidate = dimer.candidate_projectors(local)[
        'charge_minus_one_hole' if kind == 'hole' else 'charge_plus_one_particle']
    columns = _array(chain.bare_columns(kind, axis), (256, 4))
    _admit(columns.conj().T@columns-np.eye(4),
           columns[:, :2]@columns[:, :2].conj().T-np.kron(local.reference, candidate),
           columns[:, 2:]@columns[:, 2:].conj().T-np.kron(candidate, local.reference))
    return columns


def projected_coefficient(kind='hole', axis=0):
    """Exact local residue determines the first-order physical joining coefficient."""
    charge = _kind(kind)
    ops = operators(axis)
    bare = _bare(kind, axis)
    _admit(bare.conj().T@bare-np.eye(4),
           (ops['N']-(4+charge)*np.eye(256))@bare)
    effective = bare.conj().T@ops['data'].V@bare
    qa = bare.conj().T@ops['QA']@bare
    current = bare.conj().T@ops['current']@bare
    expected_charge = charge*(np.eye(4)+TAU_Z)/2
    _admit(effective+J_BARE*TAU_X, qa-expected_charge,
           current-charge*J_BARE*TAU_Y,
           current+1j*(effective@qa-qa@effective))
    return dict(coefficient=J_BARE, local_CAR_residue=LOCAL_RESIDUE,
                generator=frozen(effective), charge=frozen(qa), current=frozen(current))


def _band(kind, axis, data):
    """Do not trust a mutable cached or caller-replaced band catalogue."""
    charge = _kind(kind)
    band = chain.band(kind, axis)
    indices = np.flatnonzero(np.asarray(data.sectors) == 4+charge)
    indices = indices[np.argsort(np.asarray(data.energies)[indices])]
    vectors = np.asarray(data.vectors)[:, indices[:4]]
    expected = vectors@vectors.conj().T
    columns = _array(band.columns, (256, 4))
    gap = float(data.energies[indices[4]]-data.energies[indices[3]])
    if not math.isfinite(float(band.gap_to_complement)) or gap <= 1e-8:
        raise ArithmeticError('isolated complete lowest charged quartet required')
    _admit(_array(band.projector, (256, 256))-expected,
           columns@columns.conj().T-expected, columns.conj().T@columns-np.eye(4),
           _array(band.energies, (4,))-np.asarray(data.energies)[indices[:4]],
           band.gap_to_complement-gap)
    return band


def dressed_observables(kind='hole', axis=0):
    """The polar-dressed quartet and its actual, generally weak, charge contrast."""
    charge = _kind(kind)
    ops = operators(axis)
    bare = _bare(kind, axis)
    band = _band(kind, axis, ops['data'])
    projector = _array(band.projector, (256, 256))
    columns = _array(chain.dressed_columns(kind, axis), (256, 4))
    _admit(projector-projector.conj().T, projector@projector-projector,
           ops['K']@projector-projector@ops['K'],
           columns.conj().T@columns-np.eye(4), columns@columns.conj().T-projector)
    overlap = bare.conj().T@projector@bare
    values, basis = np.linalg.eigh(overlap)
    if values[0] <= 1e-8:
        raise ArithmeticError('canonical polar dressing has no positive overlap margin')
    canonical = projector@bare@(basis@np.diag(values**-.5)@basis.conj().T)
    _admit(columns-canonical)
    effective = columns.conj().T@ops['K']@columns
    mean = float(np.trace(effective).real/4)
    hopping = -float(np.trace(effective[:2, 2:]).real/2)
    qa = columns.conj().T@ops['QA']@columns
    qb = columns.conj().T@ops['QB']@columns
    current = columns.conj().T@ops['current']@columns
    contrast = float(np.trace(TAU_Z@qa).real/(2*charge))
    if not 0 < hopping < 1 or not 0 < contrast < 1:
        raise ArithmeticError('resolved physical hopping and finite charge contrast required')
    _admit(effective-(mean*np.eye(4)-hopping*TAU_X),
           qa-charge*(np.eye(4)+contrast*TAU_Z)/2,
           qb-charge*(np.eye(4)-contrast*TAU_Z)/2,
           current-charge*hopping*contrast*TAU_Y,
           current+1j*(effective@qa-qa@effective))
    offband_read = (np.eye(256)-projector)@ops['QA']@columns
    return dict(columns=frozen(columns), projector=frozen(projector),
                generator=frozen(effective), QA=frozen(qa), QB=frozen(qb), current=frozen(current),
                mean_energy=mean, hopping=hopping, charge_contrast=contrast,
                transfer_time=math.pi/(2*hopping), bare_overlap_eigenvalues=frozen(values),
                charge_read_offband_operator_norm=float(np.linalg.norm(offband_read, 2)),
                physical_current_amplitude=hopping*contrast)


def initial_columns(kind='hole', preparation='bare', axis=0):
    _kind(kind); _axis(axis)
    if type(preparation) is not str or preparation not in ('bare', 'dressed'):
        raise ValueError('bare product-dimer or polar-dressed preparation required')
    columns = (_bare(kind, axis) if preparation == 'bare'
               else dressed_observables(kind, axis)['columns'])
    return frozen(_array(columns, (256, 4))[:, :2])


def history(kind='hole', preparation='bare', duration=4., axis=0):
    """All occupation-output Kraus rows for an unknown internal input doublet.

    The two columns are an isometry, not a classical average. Multiplication
    by any normalized2-by-reference coefficient matrix retains that reference.
    Display probabilities below specify the equally mixed internal input.
    """
    charge, duration = _kind(kind), _duration(duration)
    ops = operators(axis)
    initial = initial_columns(kind, preparation, axis)
    data = ops['data']
    state = data.vectors@(np.exp(-1j*duration*data.energies)[:, None]*(data.vectors.conj().T@initial))
    _admit(state.conj().T@state-np.eye(2), (ops['N']-(4+charge)*np.eye(256))@state)
    probabilities = np.abs(state)**2
    blocks = np.zeros((5, 5, 2))
    for address in range(256):
        blocks[(address & 15).bit_count(), (address >> 4).bit_count()] += probabilities[address]
    qa = state.conj().T@ops['QA']@state
    qb = state.conj().T@ops['QB']@state
    current = state.conj().T@ops['current']@state
    energy = state.conj().T@ops['K']@state
    _admit(qa+qb-charge*np.eye(2), energy-initial.conj().T@ops['K']@initial,
           probabilities.sum(axis=0)-np.ones(2))
    band = _band(kind, axis, data)
    band_weight = state.conj().T@band.projector@state
    _admit(band_weight-initial.conj().T@band.projector@initial)
    bare_destination = _bare(kind, axis)[:, 2:]
    compressed = dressed_observables(kind, axis)
    dressed_destination = compressed['columns'][:, 2:]
    destination = dressed_destination.conj().T@state
    bare_destination_amplitude = bare_destination.conj().T@state
    if preparation == 'dressed':
        angle = 2*compressed['hopping']*duration
        expected = charge*(1+compressed['charge_contrast']*math.cos(angle))/2
        _admit(qa-expected*np.eye(2),
               current-charge*compressed['physical_current_amplitude']*math.sin(angle)*np.eye(2),
               destination.conj().T@destination-math.sin(angle/2)**2*np.eye(2))
    # The exterior is a complete reduced density, not a mean-number proxy.
    tensor = state.reshape(16, 16, 2)
    right_density = np.einsum('bar,car->bc', tensor, tensor.conj())/2
    left_density = np.einsum('bar,bcr->ac', tensor, tensor.conj())/2
    reference = ops['local_reference']
    coupled = np.asarray(data.ground).reshape(16, 16)
    coupled_right = coupled@coupled.conj().T
    coupled_left = coupled.T@coupled.conj()
    return dict(state=frozen(state), initial=frozen(initial),
                occupation_probabilities=frozen(probabilities),
                occupation_kraus_rows=tuple(frozen(state[j:j+1]) for j in range(256)),
                block_number_probabilities=frozen(blocks), QA=frozen(qa), QB=frozen(qb),
                current=frozen(current), energy=frozen(energy), band_weight=frozen(band_weight),
                bare_destination_weight=frozen(bare_destination_amplitude.conj().T@bare_destination_amplitude),
                dressed_destination_weight=frozen(destination.conj().T@destination),
                left_marginal_to_product_ground=float(np.linalg.norm(left_density-reference, 'nuc')/2),
                right_exterior_to_product_ground=float(np.linalg.norm(right_density-reference, 'nuc')/2),
                left_marginal_to_coupled_ground_marginal=float(np.linalg.norm(left_density-coupled_left, 'nuc')/2),
                right_exterior_to_coupled_ground_marginal=float(np.linalg.norm(right_density-coupled_right, 'nuc')/2))


def charged_response(axis=0):
    """Complete unrenormalized addition/removal measures of the coupled ground."""
    ops = operators(axis)
    data, space = ops['data'], ops['space']
    frequencies = data.energies-data.ground_energy
    addition, removal = [], []
    band_addition, band_removal = [], []
    plus, minus = _band('particle', axis, data), _band('hole', axis, data)
    first_moments = []
    for mode in range(8):
        destroy = space.annihilator(mode).toarray()
        created, removed = destroy.conj().T@data.ground, destroy@data.ground
        add = np.abs(data.vectors.conj().T@created)**2
        rem = np.abs(data.vectors.conj().T@removed)**2
        _admit(add.sum()+rem.sum()-1,
               add[np.asarray(data.sectors) != 5].sum(), rem[np.asarray(data.sectors) != 3].sum())
        for vector, weights in ((created, add), (removed, rem)):
            for order in (0, 1, 2):
                exact = np.vdot(vector, np.linalg.matrix_power(ops['K']-data.ground_energy*np.eye(256), order)@vector).real
                _admit(float(weights@frequencies**order)-exact)
        addition.append(add); removal.append(rem)
        band_addition.append(float(np.vdot(created, plus.projector@created).real))
        band_removal.append(float(np.vdot(removed, minus.projector@removed).real))
        first_moments.append((float(add@frequencies), float(rem@frequencies)))
    addition, removal = np.asarray(addition), np.asarray(removal)
    _admit(np.asarray(band_addition)-band_removal, addition.sum(axis=1)-.5,
           removal.sum(axis=1)-.5)
    return dict(ground_energy=float(data.ground_energy), frequencies=frozen(frequencies),
                numbers=frozen(data.sectors), addition_weights=frozen(addition), removal_weights=frozen(removal),
                addition_band_weights=frozen(np.asarray(band_addition)), removal_band_weights=frozen(np.asarray(band_removal)),
                first_frequency_moments=frozen(np.asarray(first_moments)),
                charged_band_gaps=frozen(np.asarray(minus.energies)-data.ground_energy))


def _display(value):
    value = round(float(value), 10)
    return 0. if value == 0 else value


def certificate():
    for axis in range(3):
        for kind in ('hole', 'particle'):
            projected_coefficient(kind, axis)
            dressed_observables(kind, axis)
    response = charged_response()
    summaries = {}
    for kind in ('hole', 'particle'):
        compressed = dressed_observables(kind)
        times = (0., 4., compressed['transfer_time']/2, compressed['transfer_time'])
        rows = []
        for preparation in ('bare', 'dressed'):
            for duration in times:
                result = history(kind, preparation, duration)
                scalar = lambda key: _display(np.trace(result[key]).real/2)
                rows.append(dict(preparation=preparation, duration=_display(duration),
                    left_excess_number=scalar('QA'), right_excess_number=scalar('QB'),
                    outward_left_CAR_current=scalar('current'), invariant_band_weight=scalar('band_weight'),
                    bare_destination_weight=scalar('bare_destination_weight'),
                    dressed_destination_weight=scalar('dressed_destination_weight'),
                    all25_block_number_probabilities=[[_display(x) for x in row]
                        for row in result['block_number_probabilities'].mean(axis=2)],
                    left_marginal_to_product_ground=_display(result['left_marginal_to_product_ground']),
                    right_exterior_to_product_ground=_display(result['right_exterior_to_product_ground']),
                    left_marginal_to_coupled_ground_marginal=_display(result['left_marginal_to_coupled_ground_marginal']),
                    right_exterior_to_coupled_ground_marginal=_display(result['right_exterior_to_coupled_ground_marginal'])))
        summaries[kind] = dict(hopping=_display(compressed['hopping']),
            actual_block_charge_contrast=_display(compressed['charge_contrast']),
            current_amplitude=_display(compressed['physical_current_amplitude']),
            transfer_time=_display(compressed['transfer_time']),
            bare_overlap_eigenvalues=[_display(x) for x in compressed['bare_overlap_eigenvalues']],
            charge_read_offband_operator_norm=_display(compressed['charge_read_offband_operator_norm']), histories=rows)
    return dict(scope='UNCHANGED_UNIFORM_FOUR_SITE_CHARGED_FIELD_OBSERVABLES',
        axes_checked=3, complete_Fock_dimension=256, charged_internal_input_dimension=2,
        complete_occupation_records_per_history=256, arbitrary_untouched_reference_columns=True,
        bare_joining_coefficient='(1+1/sqrt(17))/64 = hopping(1/16) times local CAR residue',
        canonical_polar_dressing=True, physical_charge_current_continuity_checked=True,
        charge_conjugation_checked=True, charged_dynamics=summaries,
        coupled_ground_response=dict(ground_energy=_display(response['ground_energy']),
            charged_band_gaps=[_display(x) for x in response['charged_band_gaps']],
            rank_per_charged_band=4, local_modes=8, retained_final_states_per_mode=256,
            addition_total_weights=[_display(x) for x in response['addition_weights'].sum(axis=1)],
            removal_total_weights=[_display(x) for x in response['removal_weights'].sum(axis=1)],
            addition_band_weights=[_display(x) for x in response['addition_band_weights']],
            removal_band_weights=[_display(x) for x in response['removal_band_weights']],
            addition_complement_weights=[_display(x) for x in response['addition_weights'].sum(axis=1)-response['addition_band_weights']],
            removal_complement_weights=[_display(x) for x in response['removal_weights'].sum(axis=1)-response['removal_band_weights']],
            first_frequency_moments=[[_display(x) for x in row] for row in response['first_frequency_moments']],
            all_CAR_sum_rules_and_moments_through_order2_checked=True,
            band_residue_renormalization_applied=False),
        initial_product_dimer_state_is_stationary=False,
        dressed_band_implies_strict_localized_charge=False,
        physical_electron_identification=False, native_repair_or_topological_protection=False,
        numerical_observable_admission='Frobenius residual <= 3e-10',
        complete_unrounded_state_and_measure_arrays_available=True)


report = certificate
