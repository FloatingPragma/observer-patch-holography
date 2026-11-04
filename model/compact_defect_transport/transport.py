"""Complete finite controls for the analytic projector and transport theorem.

Floating residuals check the implementation; the proof and rational gap bounds
establish the quantified statement. Local gate timing and drift cancellation
are specified control services.
"""
from fractions import Fraction as F

import numpy as np
from scipy import linalg as la

from model.compact_defect_pair import spectral
from . import bounds, source


def hamiltonian(state, radius=3, mass=64, half_separation=1):
    state = source.validate(state, half_separation)
    value = spectral.source(radius, mass, half_separation)
    h = value['h'].toarray()
    baseline = source.initial(half_separation)
    for index, point in enumerate(value['points']):
        if point not in state or state[point] == baseline[point]:
            continue
        direction = np.asarray(state[point], dtype=float)
        direction /= la.norm(direction)
        h[2*index:2*index+2, 2*index:2*index+2] = mass*sum(
            (direction[j]*spectral.SIGMA[j] for j in range(3)), np.zeros((2, 2), complex))
    return h, value


def band(h):
    eigenvalues = la.eigvalsh(h, driver='evr')
    indices = np.flatnonzero(abs(eigenvalues) < 2)
    if len(indices) != 4:
        raise ArithmeticError('four-dimensional complete compact band required')
    low, frame = la.eigh(h, subset_by_index=[int(indices[0]), int(indices[-1])], driver='evr')
    return low, frame, min(abs(eigenvalues[np.abs(eigenvalues) >= 2]))


def polar_transport(first, second):
    """Direct rotation, evaluated on span of both complete band frames."""
    first, second = np.asarray(first, complex), np.asarray(second, complex)
    if first.shape != second.shape or first.ndim != 2:
        raise ValueError('two equally ranked band isometries required')
    n, rank = first.shape
    if max(la.norm(first.conj().T@first-np.eye(rank)),
           la.norm(second.conj().T@second-np.eye(rank))) > 1e-9:
        raise ValueError('band frames must be isometries')
    q, r, _ = la.qr(np.column_stack((first, second)), mode='economic', pivoting=True)
    active = int(np.count_nonzero(abs(np.diag(r)) > 1e-12))
    q = q[:, :active]
    p0 = (q.conj().T@first)@(first.conj().T@q)
    p1 = (q.conj().T@second)@(second.conj().T@q)
    identity = np.eye(active)
    s = p1@p0+(identity-p1)@(identity-p0)
    if min(la.svdvals(s)) < 1e-8:
        raise ValueError('direct rotation is singular')
    small, _ = la.polar(s)
    return np.eye(n)+q@(small-identity)@q.conj().T


def snake_path(points):
    """All spatial-spin coordinates, adjacent onsite or on lattice neighbors."""
    radius = max(abs(x) for point in points for x in point)
    grid = list(range(-radius, radius+1))
    addresses = {p: i for i, p in enumerate(points)}
    path = []
    row = 0
    for i, x in enumerate(grid):
        for y in (grid if i % 2 == 0 else grid[::-1]):
            for z in (grid if row % 2 == 0 else grid[::-1]):
                path.extend((2*addresses[(x, y, z)], 2*addresses[(x, y, z)]+1))
            row += 1
    if len(path) != 2*len(points) or set(path) != set(range(len(path))):
        raise ValueError('complete cubic carrier required')
    for first, second in zip(path, path[1:]):
        if sum(abs(a-b) for a, b in zip(points[first//2], points[second//2])) > 1:
            raise ArithmeticError('nonlocal edge in physical control path')
    return np.asarray(path)


def local_compile(unitary, path):
    """Adjacent-row elimination and complete reverse replay of every gate."""
    n = len(path)
    if sorted(map(int, path)) != list(range(n)) or unitary.shape != (n, n):
        raise ValueError('complete permutation of the unitary coordinates required')
    target = unitary[np.ix_(path, path)]
    if la.norm(target.conj().T@target-np.eye(n)) > 1e-8:
        raise ValueError('unitary control required')
    matrix = target.copy()
    gates = []
    for column in range(n-1):
        for bottom in range(n-1, column, -1):
            a, b = matrix[bottom-1, column], matrix[bottom, column]
            if b == 0:
                continue
            magnitude = np.hypot(abs(a), abs(b))
            c, s = a/magnitude, b/magnitude
            first, second = matrix[bottom-1].copy(), matrix[bottom].copy()
            matrix[bottom-1] = c.conjugate()*first+s.conjugate()*second
            matrix[bottom] = -s*first+c*second
            gates.append((bottom, c, s))
    diagonal = np.diag(matrix)
    phases = diagonal/abs(diagonal)
    replay = np.diag(phases)
    for bottom, c, s in reversed(gates):
        first, second = replay[bottom-1].copy(), replay[bottom].copy()
        replay[bottom-1] = c*first-s.conjugate()*second
        replay[bottom] = s*first+c.conjugate()*second
    error = la.norm(replay-target)
    if error > 1e-8:
        raise ArithmeticError('complete local gate replay failed')
    return {'spatial_spin_coordinates': n, 'nearest_neighbor_two_mode_gates': len(gates),
            'two_mode_gate_count_upper': n*(n-1)//2, 'onsite_phase_gates': n,
            'complete_replay_Frobenius_error_below_1e_minus_8': bool(error < 1e-8),
            'onsite_or_neighbor_controls_only': True,
            'supplied_pulse_timing_and_drift_cancellation': True,
            'compiled_into_unchanged_native_CE_repairs': False}


def control():
    from model.common_source.chiral_candidate import action
    states, records, transition = source.continuing_protocol()
    sources = [hamiltonian(state) for state in states]
    bands = [band(h) for h, _ in sources]
    h0, native = sources[0]
    c = native['C'].toarray()
    path = snake_path(native['points'])
    exact = bounds.certificate()
    beta = float(F(exact['spectral_projection_difference_norm_upper']))
    epsilon = float(F(exact['any_two_admitted_source_difference_norm_upper']))
    rows, unitaries = [], {}
    current_norm = max(la.norm(np.asarray(g, complex)) for g in action.generators())
    for before, after in ((0, 1), (0, 2), (1, 2)):
        h0 = sources[before][0]
        h1 = sources[after][0]
        w0 = bands[before][1]
        p0 = w0@w0.conj().T
        pair0 = w0.conj().T@c@w0.conjugate()
        eigenvalues, raw_frame, high = bands[after]
        u = polar_transport(w0, raw_frame)
        unitaries[before, after], unitaries[after, before] = u, u.conj().T
        w1 = u@w0
        p1 = raw_frame@raw_frame.conj().T
        difference = la.norm(p1-p0, 2)
        delta = la.norm(h1-h0, 2)
        raw_leakage = la.norm(w0-raw_frame@(raw_frame.conj().T@w0), 2)**2
        controlled_leakage = la.norm(w1-raw_frame@(raw_frame.conj().T@w1), 2)**2
        pairing_error = la.norm(u@c@u.T-c)
        a0, a1 = w0.conj().T@h0@w0, w1.conj().T@h1@w1
        current_error = la.norm(w1-raw_frame@(raw_frame.conj().T@w1))*current_norm
        residuals = [la.norm(u.conj().T@u-np.eye(len(u))), pairing_error,
                     la.norm(w1.conj().T@c@w1.conjugate()-pair0), current_error,
                     la.norm(h1@w1-w1@a1)]
        if max(residuals) > 1e-9 or controlled_leakage > 1e-18:
            raise ArithmeticError('complete band, pair or native current transport failed')
        if difference > beta or delta > epsilon or raw_leakage <= 1e-20:
            raise ArithmeticError('resolvent bound or bare quench negative control failed')
        if high < 62 or max(abs(eigenvalues)) > float(F(exact['band_abs_energy_upper'])):
            raise ArithmeticError('the common compact spectral bounds failed')
        rows.append(dict(source_state=before, target_state=after,
                         reverse_is_adjoint_control=True,
                         source_operator_difference_norm=round(delta, 12),
                         projection_difference_norm=round(difference, 12),
                         bare_quench_leakage_probability=round(raw_leakage, 16),
                         controlled_band_leakage_below_1e_minus_18=True,
                         all_transport_residuals_below_1e_minus_9=True,
                         pair_tensor_preserved=True,
                         band_energy_operator_change_norm=round(la.norm(a1-a0, 2), 12),
                         numerical_band_energies=[round(float(e), 12) for e in eigenvalues],
                         local_gate_compilation=local_compile(u, path)))
    for before in range(len(states)):
        completeness = np.zeros(h0.shape, complex)
        for record in records:
            if record['source_state'] != before:
                continue
            after = record['target_state']
            u = np.eye(len(h0)) if before == after else unitaries[before, after]
            completeness += .25*u.conj().T@u
        if la.norm(completeness-np.eye(len(h0))) > 1e-9:
            raise ArithmeticError('complete retained instrument does not preserve trace')
    return {'baseline_numerical_band_energies':[round(float(e), 12) for e in bands[0][0]],
            'source_spatial_spin_coordinates': len(h0), 'native_internal_dimension':45,
            'native_current_count': len(action.generators()), 'localized_annihilators':180,
            'complete_instrument_trace_preservation_below_1e_minus_9':True,
            'all_continuing_records_retained': len(records),
            'nonidentity_transport_pairs': rows,
            'numerical_diagnostics_are_exact_spectral_certificates':False}
