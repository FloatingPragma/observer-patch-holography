"""Exact source-instance attachment to existing finite energy Ward theorems."""
from fractions import Fraction as F
from functools import lru_cache

import numpy as np

from model.native_ordered_ground import source
from model.protected_defect_catalog import charge
from model.source_energy_effort.energy import expectation
from model.spin_modular_geometry import regional


def norm_bound(matrix, denominator=1):
    return F(max(sum(abs(int(x)) for x in row) for row in matrix), denominator)


def permutation(axis):
    return np.array([sum(((word >> x) & 1) << (x ^ (1 << axis)) for x in range(8))
                     for word in range(256)])


@lru_cache(maxsize=3, typed=True)
def operators(axis):
    if type(axis) is not int or not 0 <= axis < 3:
        raise ValueError('one source spatial axis required')
    data = regional.state()
    h = data['energy_integer']
    density = tuple(regional.integer_energy(x) for x in range(8))
    coordinates = tuple(2*((x >> axis) & 1)-1 for x in range(8))
    b = sum((s*e for s, e in zip(coordinates, density, strict=True)), np.zeros_like(h))
    c = h@b-b@h
    force = -(h@c-c@h)
    d0 = -(b@c-c@b)+256*h  # D = D0/2048 - E0 I
    currents = {(x, y): density[x]@density[y]-density[y]@density[x]
                for x in range(8) for y in range(x+1, 8)}
    weighted = sum(((coordinates[y]-coordinates[x])*j for (x, y), j in currents.items()),
                   np.zeros_like(h))
    if not np.array_equal(c, weighted):
        raise ArithmeticError('boost-polarization derivative differs from source energy flux')
    for x in range(8):
        outgoing = sum((currents[x, y] if x < y else -currents[y, x]
                        for y in range(8) if y != x), np.zeros_like(h))
        if np.any(h@density[x]-density[x]@h+outgoing):
            raise ArithmeticError('source instance fails the inherited finite continuity premise')
    if not np.any(force):
        raise ArithmeticError('source energy current unexpectedly became a conserved momentum')
    fields = []
    indices = np.arange(256)
    odd = source.cube()[1]
    for x in range(8):
        field = np.zeros_like(h)
        field[indices ^ (1 << x), indices] = 1 if x in odd else -1
        fields.append(field)
    h0 = h-sum(fields)
    b0 = b-sum(s*field for s, field in zip(coordinates, fields, strict=True))
    if not np.array_equal(h0@b0-b0@h0, c):
        raise ArithmeticError('bipartite selection-field energy-current cancellation failed')
    zero_field_force = -(h0@c-c@h0)
    if not np.any(zero_field_force):
        raise ArithmeticError('finite-box current force was incorrectly attributed only to selection field')
    q = np.zeros_like(h)
    for x in range(8): q[indices ^ (1 << x), indices] += 1
    if any(np.any(q@m-m@q) for m in (h, b, c)):
        raise ArithmeticError('source boost candidate lost the same conserved charge')
    # Actual selected-field source symmetry: coordinate reflection followed by
    # global pi spin rotation around Z. It reverses B and preserves H.
    p = permutation(axis)
    sign = np.array([(-1)**word.bit_count() for word in range(256)], dtype=np.int64)
    transform = lambda m: sign[:, None]*m[np.ix_(p, p)]*sign[None, :]
    if not np.array_equal(transform(h), h) or not np.array_equal(transform(b), -b):
        raise ArithmeticError('source symmetry does not certify zero ground boost mean')
    return {'axis': axis, 'H8': h, 'densities8': density, 'B16': b, 'P_imag128': c,
            'force1024': force, 'boost_defect_without_origin2048': d0,
            'zero_field_H8': h0, 'zero_field_force1024': zero_field_force,
            'currents_imag64': currents, 'twice_coordinates': coordinates,
            'twice_charge': q}


@lru_cache(maxsize=1)
def source_controls():
    rows = []
    data = regional.state()
    vector, epsilon = data['vector'], data['ground_vector_error']
    ground = data['full_source_certificate']['ground_energy_interval']
    elo, ehi = F(ground[0]), F(ground[1])
    for axis in range(3):
        op = operators(axis)
        h, b = op['H8'], op['B16']
        b2 = b@b
        second = expectation(b2, vector, 256)
        bnorm = norm_bound(b, 16)
        second_error = 2*bnorm**2*epsilon
        first = expectation(b@h@b, vector, 2048)
        # Bounds for <B(H-E0)B> use the actual source norm and its certified E0 interval.
        moment_lo = first-ehi*second-2*bnorm**2*(14+max(abs(elo), abs(ehi)))*epsilon
        moment_hi = first-elo*second+2*bnorm**2*(14+max(abs(elo), abs(ehi)))*epsilon
        if second-second_error <= F(24, 1000) or moment_lo <= F(65, 1000):
            raise ArithmeticError('selected-ground boost response lost strict source weight')
        d0 = op['boost_defect_without_origin2048']
        d_expect = expectation(d0, vector, 2048)
        defect_bound = norm_bound(d0, 2048)+max(abs(elo), abs(ehi))
        dlo, dhi = d_expect-ehi-2*defect_bound*epsilon, d_expect-elo+2*defect_bound*epsilon
        if dhi >= F(-13, 100):
            raise ArithmeticError('same-ground boost Ward defect lost its strict witness')
        # The finite trace obstruction applies to any Hermitian candidate,
        # not just this weighted-energy realization.
        trace_l_lower = F(int(np.trace(h)), 8)-256*ehi
        if trace_l_lower <= 0:
            raise ArithmeticError('positive source trace obstruction failed')
        rows.append({'axis': axis, 'boost_norm_upper': bnorm,
                     'energy_flux_norm_upper': norm_bound(op['P_imag128'], 128),
                     'energy_flux_force_norm_upper': norm_bound(op['force1024'], 1024),
                     'zero_field_energy_flux_force_norm_upper': norm_bound(op['zero_field_force1024'], 1024),
                     'selection_field_flux_torque_norm_upper': norm_bound(op['force1024']-op['zero_field_force1024'], 1024),
                     'H_B_and_H_J_energy_flux_are_same_operator': True,
                     'prepared_state_replaced_by_zero_field_ground': False,
                     'boost_defect_norm_upper': defect_bound,
                     'ground_boost_mean_exact': F(0),
                     'ground_boost_second_moment_lower': second-second_error,
                     'ground_boost_second_moment_upper': second+second_error,
                     'ground_boost_energy_moment_lower': moment_lo,
                     'ground_boost_energy_moment_upper': moment_hi,
                     'ground_boost_defect_lower': dlo, 'ground_boost_defect_upper': dhi,
                     'positive_trace_of_relative_H_lower': trace_l_lower,
                     'same_charge_commutators_exactly_zero': True,
                     'energy_flux_is_conserved_spatial_momentum': False})
    return rows


@lru_cache(maxsize=1)
def coarse_transport():
    """Instantiate block Ward sums and retain their exact subcell dipole."""
    h = operators(0)['H8']
    count, nonzero_dipoles = 0, 0
    for axis in range(3):
        op = operators(axis)
        densities, currents = op['densities8'], op['currents_imag64']
        for mask in range(8):
            groups = sorted(set(x & mask for x in range(8)))
            coarse_e = {g: sum((densities[x] for x in range(8) if x & mask == g), np.zeros_like(h))
                        for g in groups}
            coarse_j = {(a, b): sum((j if x & mask == a and y & mask == b else
                                     -j if x & mask == b and y & mask == a else np.zeros_like(h)
                                     for (x, y), j in currents.items()), np.zeros_like(h))
                        for a in groups for b in groups if a < b}
            for a in groups:
                outgoing = sum((coarse_j[a, b] if a < b else -coarse_j[b, a]
                                for b in groups if b != a), np.zeros_like(h))
                if np.any(h@coarse_e[a]-coarse_e[a]@h+outgoing):
                    raise ArithmeticError('instantiated coarse energy Ward identity failed')
            # Representatives are geometric cell means. Missing coordinate bits
            # have mean0; retained ones have the actual ±1/2 coordinate.
            coarse_b = sum(((2*((g >> axis) & 1)-1)*coarse_e[g]
                            for g in groups if mask >> axis & 1), np.zeros_like(h))
            dipole = op['B16']-coarse_b
            coarse_flux = h@coarse_b-coarse_b@h
            if not np.array_equal(op['P_imag128'], coarse_flux+h@dipole-dipole@h):
                raise ArithmeticError('coarse spatial readout omitted the subcell energy dipole')
            nonzero_dipoles += int(np.any(dipole))
            count += 1
    return {'axis_block_maps_checked': count,
            'maps_with_nonzero_subcell_energy_dipole': nonzero_dipoles,
            'matrix_valued_CoarseEnergyWard_premises_instantiated': True,
            'energy_polarization_and_flux_transport_with_dipole_exact': True,
            'discarding_subcell_dipoles_preserves_all_boost_observables': False}


@lru_cache(maxsize=1)
def forced_stress():
    """Three local forced columns for the actual energy-flux candidates."""
    results = []
    source_data = regional.state()
    vector, epsilon = source_data['vector'], source_data['ground_vector_error']
    # At site0 the Marshall diagonal is trivial. These exact physical A_q
    # matrices are the same coherent local creation used in the energy readout.
    creation = charge.operators()['twice_creation']
    for axis in range(3):
        op = operators(axis)
        h, e, c = op['H8'], op['densities8'], op['P_imag128']
        coordinates, currents = op['twice_coordinates'], op['currents_imag64']
        density = []
        for x in range(8):
            value = np.zeros_like(h)
            for y in range(8):
                if y == x: continue
                j = currents[x, y] if x < y else -currents[y, x]
                value += (coordinates[y]-coordinates[x])*j
            density.append(value)  # p_x=i*density_x/256
        if not np.array_equal(sum(density), 2*c):
            raise ArithmeticError('local energy-flux densities lost their total')
        stress = {(x, y): -(e[x]@density[y]-density[y]@e[x])
                             +(e[y]@density[x]-density[x]@e[y])
                  for x in range(8) for y in range(x+1, 8)}
        force = [-2*(e[x]@c-c@e[x]) for x in range(8)]
        delta = [-(h@p-p@h) for p in density]
        for x in range(8):
            outgoing = sum((stress[x, y] if x < y else -stress[y, x]
                            for y in range(8) if y != x), np.zeros_like(h))
            if not np.array_equal(delta[x]+outgoing, force[x]):
                raise ArithmeticError('actual local forced stress identity failed')
        if not np.array_equal(sum(force), 2*op['force1024']):
            raise ArithmeticError('total force differs from the actual current derivative')
        for mask in range(8):
            groups = sorted(set(x & mask for x in range(8)))
            for a in groups:
                inside = [x for x in range(8) if x & mask == a]
                coarse_delta = sum((delta[x] for x in inside), np.zeros_like(h))
                coarse_force = sum((force[x] for x in inside), np.zeros_like(h))
                crossing = np.zeros_like(h)
                for (x, y), value in stress.items():
                    if x in inside and y not in inside: crossing += value
                    if y in inside and x not in inside: crossing -= value
                if not np.array_equal(coarse_delta+crossing, coarse_force):
                    raise ArithmeticError('forced stress violates coarse_forced_balance')
        channels = []
        for q, a in sorted(creation.items()):
            observable = a.T@op['force1024']@a
            value = expectation(observable, vector, 4096)
            error = 2*norm_bound(observable, 4096)*epsilon
            channels.append({'charge': q,
                             'unconditioned_energy_flux_acceleration_lower': value-error,
                             'unconditioned_energy_flux_acceleration_upper': value+error})
        results.append({'axis': axis, 'local_balance_cells': 8,
                        'nonzero_antisymmetric_stress_pairs': sum(int(np.any(s)) for s in stress.values()),
                        'all_eight_coordinate_block_maps_preserve_force': True,
                        'total_force_nonzero_as_source_operator': True,
                        'closed_internal_stress_can_absorb_total_force': False,
                        'created_charge_channel_force': channels,
                        'candidate_densities_are_physical_momentum_densities': False})
    return results


def speed_normalization(speed):
    if type(speed) is not F or speed <= 0:
        raise ValueError('positive exact candidate clock/length conversion required')
    rows = []
    for source_row in source_controls():
        rows.append({'axis': source_row['axis'], 'candidate_speed': speed,
                     'normalized_boost_norm_upper': source_row['boost_norm_upper']/speed,
                     'normalized_energy_flux_norm_upper': source_row['energy_flux_norm_upper']/speed,
                     'ground_normalized_boost_defect_lower':
                         -2*source_row['ground_boost_energy_moment_upper']/speed**2,
                     'ground_normalized_boost_defect_upper':
                         -2*source_row['ground_boost_energy_moment_lower']/speed**2,
                     'physical_speed_selected_by_this_normalization': False})
    return rows


def finite_translations():
    """Actual charge-preserving even cube shifts, with finite Fourier logs."""
    op = operators(0)
    h, charge_matrix = op['H8'], op['twice_charge']
    data = regional.state()
    vector = tuple(map(int, data['vector']))
    norm = sum(v*v for v in vector)
    rows = []
    eye = np.eye(256, dtype=np.int64)
    z0 = np.diag(np.array([1-2*(word & 1) for word in range(256)], dtype=np.int64))
    flip0 = eye[np.arange(256)^1]
    for shift in (3, 5, 6):
        p = np.array([sum(((word >> x) & 1) << (x ^ shift) for x in range(8))
                      for word in range(256)])
        t = eye[:, p]
        twice_log_over_pi = eye-t
        if (not np.array_equal(t.T, t) or not np.array_equal(t@t, eye)
                or np.any(h@t-t@h) or np.any(charge_matrix@t-t@charge_matrix)
                or np.any(h@twice_log_over_pi-twice_log_over_pi@h)):
            raise ArithmeticError('actual source shift failed the finite polynomial momentum premises')
        for x in range(8):
            if not np.array_equal(t@op['densities8'][x]@t, op['densities8'][x ^ shift]):
                raise ArithmeticError('finite source shift does not transport its energy readout')
        for q in (-1, 1):
            z = np.diag(np.array([1-2*((word >> shift) & 1) for word in range(256)], dtype=np.int64))
            flip = eye[np.arange(256)^(1 << shift)]
            if not np.array_equal(t@(z0+q*flip0@z0)@t, z+q*flip@z):
                raise ArithmeticError('finite source shift does not transport coherent charge creation')
        difference = F(sum((vector[int(p[j])]-vector[j])**2 for j in range(256)), norm)
        if difference >= F(1, 16) or 2*data['ground_vector_error'] >= F(1, 4):
            raise ArithmeticError('selected ground source-shift eigenvalue not fixed to +1')
        flipped = []
        for axis in range(3):
            c = operators(axis)['P_imag128']
            if shift >> axis & 1:
                if not np.array_equal(t@c@t, -c):
                    raise ArithmeticError('energy-flux sign under actual source shift changed')
                flipped.append(axis)
        rows.append({'coordinate_xor_shift': shift, 'order': 2,
                     'source_energy_density_and_local_creation_transport_exact': True,
                     'same_H_B_charge_and_selected_ground_preserved': True,
                     'principal_log_over_pi': '(I-T)/2',
                     'principal_log_conserved_by_same_H_B': True,
                     'energy_flux_axes_reversed_by_shift': flipped,
                     'energy_flux_is_this_conserved_translation_log': False,
                     'finite_log_is_continuum_physical_momentum': False})
    return {'even_shift_controls': rows,
            'nonidentity_shifts_form_Z2_times_Z2': True,
            'three_independent_continuum_momenta_obtained': False}
