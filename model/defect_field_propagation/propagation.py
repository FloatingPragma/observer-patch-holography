"""A supplied native-seam ring, complete CAR blocks, and literal source gates.

The finite spatial law is a symmetric product of the calibrated quadratic
terms. Its static sum and complete lattice dispersion are separately checked.
No native balancing transition is identified with a Weyl unitary.
"""
from __future__ import annotations

from functools import lru_cache
import math

import numpy as np
from scipy.linalg import expm

from ..common_source import calibrated_matter
from ..common_source.chiral_candidate import gaussian_source as source
from ..defect_field_benchmark import quasiparticles
from ..exchange_sector import statistics as native


CYCLE = (0, 1, 5, 4, 2, 3, 9, 8, 6, 7, 10, 11)
LOCAL_ADMISSION = 2e-10


def geometry():
    edges = tuple((CYCLE[i], CYCLE[(i+1) % 12]) for i in range(12))
    if set(CYCLE) != set(range(12)) or any(tuple(sorted(edge)) not in native.SEAMS for edge in edges):
        raise ArithmeticError('the supplied ring left the actual native seam graph')
    return {'cycle': list(CYCLE), 'oriented_seams': [list(e) for e in edges],
            'dimension': 1, 'coordinate_and_orientation_are_supplied': True,
            'native_balancing_is_not_Weyl_propagation': True}


def derivative(sites=12, spacing=1.):
    if type(sites) is not int or sites < 3 or not math.isfinite(spacing) or spacing <= 0:
        raise ValueError('at least three sites and positive finite spacing required')
    matrix = np.zeros((sites, sites), complex)
    for a in range(sites):
        b = (a+1) % sites
        matrix[a, b], matrix[b, a] = -0.5j/spacing, 0.5j/spacing
    return matrix


@lru_cache(maxsize=1)
def local_pairing():
    result = np.array(calibrated_matter.pairing_kernel(), complex)
    if result.shape != (90, 90) or not np.array_equal(result, -result.T):
        raise ArithmeticError('calibrated local pairing is not the retained90-mode CAR kernel')
    result.flags.writeable = False
    return result


@lru_cache(maxsize=1)
def components():
    delta = local_pairing()
    seen, result = set(), []
    for start in range(90):
        if start in seen:
            continue
        found, pending = {start}, [start]
        seen.add(start)
        while pending:
            for child in np.flatnonzero(delta[pending.pop()]):
                child = int(child)
                if child not in seen:
                    seen.add(child)
                    found.add(child)
                    pending.append(child)
        result.append(tuple(sorted(found)))
    if sorted(i for row in result for i in row) != list(range(90)):
        raise ArithmeticError('internal component partition omitted a mode')
    return tuple(result)


def block(indices, *, sites=12, spacing=1., speed=1.):
    if tuple(indices) not in components() or not math.isfinite(speed) or speed <= 0:
        raise ValueError('one complete calibrated pairing component and positive speed required')
    indices = np.array(indices, int)
    signs = np.where(indices < 45, 1., -1.)
    h = -speed*np.kron(derivative(sites, spacing), np.diag(signs))
    delta = np.kron(np.eye(sites), local_pairing()[np.ix_(indices, indices)])
    return h, delta


def nambu(h, delta):
    return np.block([[h, delta], [delta.conj().T, -h.T]])


def majorana_change(modes):
    change = np.zeros((2*modes, 2*modes), complex)
    for mode in range(modes):
        change[2*mode, mode] = change[2*mode, modes+mode] = 1
        change[2*mode+1, mode], change[2*mode+1, modes+mode] = -1j, 1j
    return change


def nambu_from_source(orthogonal):
    change = majorana_change(len(orthogonal)//2)
    return change.conj().T@orthogonal@change/2


def audit_local_catalog(catalog):
    """Compare every executed source word with its independent four-state H."""
    # Independent bit-string CAR matrices, including the occupied spectator.
    ladders = []
    for mode in (0, 1):
        matrix = np.zeros((4, 4), complex)
        for bits in range(4):
            if bits >> mode & 1:
                matrix[bits ^ (1 << mode), bits] = (-1)**((bits & ((1 << mode)-1)).bit_count())
        ladders.append(matrix)
    a, b = ladders
    for item in catalog.values():
        z = complex(*item['coefficient'])
        kind = item['kind']
        if kind == 'fswap':
            expected = np.array([[1, 0, 0, 0], [0, 0, 1, 0],
                                 [0, 1, 0, 0], [0, 0, 0, -1]], complex)
        else:
            if kind == 'onsite':
                H = z.real*a.conj().T@a
            elif kind == 'hopping':
                H = z*a.conj().T@b+z.conjugate()*b.conj().T@a
            elif kind == 'pairing':
                H = z*a.conj().T@b.conj().T+z.conjugate()*b@a
            else:
                raise ValueError('unknown local source generator')
            expected = expm(-1j*item['duration']*H)
        recorded = item['source_fock_unitary']
        actual = np.array(recorded['real'])+1j*np.array(recorded['imag'])
        if np.linalg.norm(actual-expected, 2) > LOCAL_ADMISSION:
            raise ArithmeticError('a full local source unitary differs from its CAR Hamiltonian')


def product_error_bound(h, delta, duration, steps):
    """Conservative full-Fock bound from equal Taylor coefficients to order2."""
    S = (sum(abs(h[i, i]) for i in range(len(h)))
         +sum(abs(h[i, j])+abs(delta[i, j])
              for i in range(len(h)) for j in range(i+1, len(h))))
    x = abs(duration)*float(S)/steps
    if x >= 4:
        return 2.
    term, tail = x*x*x/6, 0.
    for order in range(3, 200):
        tail += term
        term *= x/(order+1)
        if term < 1e-20:
            break
    return min(2., 2*steps*tail)


def dispersion():
    """All2160 Nambu eigenvalues, with every doubler and zero direction."""
    delta = local_pairing()
    signs = np.diag(np.r_[np.ones(45), -np.ones(45)])
    d0 = nambu(np.zeros_like(delta), delta)
    kinetic = np.block([[-signs, np.zeros_like(delta)], [np.zeros_like(delta), -signs]])
    if np.max(abs(kinetic@d0+d0@kinetic)) > 2e-11:
        raise ArithmeticError('same-action mass and Weyl kinetic factors do not anticommute')
    all_values, rows = [], []
    for m in range(12):
        k = 2*math.pi*m/12
        p = math.sin(k)
        matrix = d0+p*kinetic
        if np.max(abs(matrix@matrix-d0@d0-p*p*np.eye(180))) > 2e-10:
            raise ArithmeticError('complete finite lattice mass-dispersion identity failed')
        values = np.linalg.eigvalsh(matrix)
        all_values.append(values)
        rows.append({'momentum_index': m, 'momentum': k,
                     'finite_difference_momentum': p, 'nambu_directions': 180,
                     'massless_zero_directions': 12 if m in (0, 6) else 0})
    return {'momenta': rows, 'all_nambu_directions': 2160, 'physical_CAR_modes': 1080,
            'kinetic_zero_momentum_indices': [0, 6], 'zero_nambu_directions': 24,
            'all_doublers_retained': True,
            'identity': 'B(k)^2=B(0)^2+sin(k)^2 I; unit spacing and speed',
            'three_dimensional_extension': 'three commuting ring derivatives and Pauli factors give sum sin(k_a)^2 plus same mass square; eight nodes retained',
            'three_dimensional_source_execution_claimed': False}, np.stack(all_values)


def execute_block(indices, current=0, duration=1/64, steps=2, covariance=None):
    if type(current) is not int or current not in (-1, 0, 1):
        raise ValueError('signed native current required')
    h, delta = block(indices)
    effective = duration*(1+current/4)
    result = source.evolve(h, delta, effective, steps=steps, covariance=covariance)
    audit_local_catalog(result.local_catalog)
    accumulated = math.expm1(result.resources['completed_local_blocks']*math.log1p(LOCAL_ADMISSION))
    return result, {'effective_duration': effective,
        'clock_rule': 'entire action multiplied by 1+native_current/4',
        'sum_action_nambu_error': float(np.linalg.norm(nambu_from_source(result.orthogonal)
                                        -expm(-1j*effective*nambu(h, delta)), 2)),
        'executed_law': 'exact specified symmetric local product; not exact exponential of the sum',
        'ideal_product_to_static_sum_full_Fock_operator_bound': product_error_bound(h, delta, effective, steps),
        'source_local_operator_admission': LOCAL_ADMISSION,
        'source_product_full_Fock_operator_bound': accumulated,
        'source_product_normalized_input_TV_bound': accumulated+accumulated**2/2,
        'error_bound_scope': 'analytic composition conditional on numerical local admission; not outward rounding',
        'resources': result.resources,
        'diagnostics': result.diagnostics}


def purification(modes, occupied_mode=None):
    """Source-prepare faithful work occupations with all purifying CAR modes."""
    delta = np.zeros((2*modes, 2*modes), complex)
    for i in range(modes):
        angle = math.pi/3 if i == occupied_mode else math.pi/6
        delta[i, modes+i], delta[modes+i, i] = angle, -angle
    vacuum = source.vacuum_covariance(2*modes)
    result = source.evolve(np.zeros_like(delta), delta, 1., covariance=vacuum)
    return result


def number_matrix(covariance):
    modes = len(covariance)//2
    annihilator = np.zeros((modes, 2*modes), complex)
    for i in range(modes):
        annihilator[i, 2*i], annihilator[i, 2*i+1] = .5, .5j
    return annihilator.conj()@(np.eye(2*modes)-1j*covariance)@annihilator.T


def wavepacket():
    """A localized faithful excess and its complete source QND number records."""
    indices = min((row for row in components() if len(row) == 2),
                  key=lambda row: abs(local_pairing()[row[0], row[1]]))
    h, delta = block(indices)
    modes = len(h)
    arrays, records, profiles = {}, [], []
    for label, raised in [('reference', None), ('localized_excess', 0)]:
        prepared = purification(modes, raised)
        evolved, metadata = execute_block(indices, 1, duration=.4, steps=4,
                                          covariance=prepared.covariance)
        gamma = evolved.covariance
        profile = source.occupations(gamma)[:modes]
        profiles.append(profile)
        arrays[label+'_covariance_with_purifiers'] = gamma
        arrays[label+'_all_occupations'] = profile
        number = number_matrix(gamma)[:modes, :modes]
        left, right = 0, len(indices)
        current_value = -2*float(np.imag(h[left, right]*number[left, right]))
        # Rotate the same bond current to a difference of two occupations.
        probe_h = np.zeros_like(h)
        value = -h[left, right]/abs(h[left, right])
        probe_h[left, right], probe_h[right, left] = value, value.conjugate()
        probed = source.evolve(probe_h, np.zeros_like(delta), math.pi/4,
                              covariance=gamma)
        pointer_rows = []
        for site in (left, right):
            read = source.occupation_read(probed.covariance, site)
            pointer_rows.append({'mode': site,
                'probabilities': [branch.probability for branch in read.branches],
                'resources': read.resources})
            for branch in read.branches:
                if branch.covariance is not None:
                    arrays[f'{label}_current_mode{site}_outcome{branch.occupation}'] = branch.covariance
        observed = abs(h[left, right])*(pointer_rows[0]['probabilities'][1]
                                        -pointer_rows[1]['probabilities'][1])
        if abs(observed-current_value) > 2e-8:
            raise ArithmeticError('source bond-current pointer differs from complete covariance')
        records.append({'state': label, 'initial_work_occupations': '.25 everywhere, .75 at raised mode',
            'supplied_blank_CAR_modes_counted': 2*modes,
            'raised_mode': raised, 'preparation_resources': prepared.resources,
            'evolution': metadata, 'bond_current': current_value,
            'source_current_readout': observed, 'probe_resources': probed.resources,
            'pointer_records': pointer_rows})
    excess = profiles[1]-profiles[0]
    arrays['occupation_excess_profile'] = excess
    if max(abs(excess[len(indices):])) < 1e-5:
        raise ArithmeticError('the native-seam source wavepacket did not propagate')
    return {'internal_component': list(indices), 'local_mass': abs(local_pairing()[indices[0], indices[1]]),
            'retained_work_modes': modes, 'retained_purifying_modes': modes,
            'all_work_occupation_marginals_retained': True,
            'current_read_is_difference_of_separate_complete_pointer_experiments': True,
            'records': records}, arrays


def build():
    geometry_record = geometry()
    spectrum, values = dispersion()
    arrays = {'complete_nambu_eigenvalues': values}
    rows = []
    for i, indices in enumerate(components()):
        for current in (-1, 0, 1):
            result, record = execute_block(indices, current)
            arrays[f'component_{i}_current_{current}_complete_majorana'] = result.orthogonal
            # Local source blocks are the actual full-four-column matrices,
            # including their scalar phases; keep every distinct local word.
            record['local_source_catalog'] = result.local_catalog
            rows.append({'component': i, 'internal_indices': list(indices),
                         'current': current, **record})
    packet, packet_arrays = wavepacket()
    arrays.update(packet_arrays)
    report = {'schema': 'oph.native_seam_calibrated_propagation.v1',
              'geometry': geometry_record, 'spectrum': spectrum,
              'component_sizes': [len(row) for row in components()],
              'all_1080_CAR_modes_and_all_occupation_sectors_retained': True,
              'same_port_collisions_retained': True,
              'full_fock_dimension': '2^1080',
              'static_positive_quasiparticle_number_preserved_by_kinetic_claimed': False,
              'hard_core_colored_native_domain_invariant_claimed': False,
              'complete_native_instrument_extension_required_and_separate': True,
              'source_execution': rows, 'wavepacket': packet,
              'physical_three_plus_one_dimensional_identification_complete': False,
              'floating_execution_is_not_outward_certified': True}
    return report, arrays
