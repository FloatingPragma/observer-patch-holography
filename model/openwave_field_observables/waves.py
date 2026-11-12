"""Finite quantum Heisenberg shortcut and independently propagated cochain waves."""
from itertools import product
import numpy as np
import sympy as s
from . import lattice as l


def cutoff_theorem(cap=3):
    """Check a complete two-mode cutoff; the proof covers arbitrary mode count."""
    states = tuple(p for p in product(range(cap+1), repeat=2) if sum(p) <= cap)
    index = {p: i for i, p in enumerate(states)}
    freq = (s.Integer(2), s.sqrt(6))
    h = s.diag(*(sum(w*n for w, n in zip(freq, state)) for state in states))
    annihilators = []
    bare = s.zeros(len(states))
    for j, w in enumerate(freq):
        a = s.zeros(len(states))
        for column, state in enumerate(states):
            if state[j]:
                target = list(state); target[j] -= 1
                a[index[tuple(target)], column] = s.sqrt(state[j])
        q = (a+a.H)/s.sqrt(2*w)
        p = -s.I*s.sqrt(w/2)*(a-a.H)
        assert (s.I*(h*q-q*h)-p).applyfunc(s.simplify) == s.zeros(len(states))
        assert (s.I*(h*p-p*h)+w*w*q).applyfunc(s.simplify) == s.zeros(len(states))
        bare += w*(a*a.H+a.H*a)/2
        annihilators.append(a)
    top = s.diag(*(int(sum(p) == cap) for p in states))
    vacuum = sum(freq)/2
    correction = top*(h+2*vacuum*s.eye(len(states)))/2
    assert (bare+correction-h-vacuum*s.eye(len(states))).applyfunc(s.simplify) == s.zeros(len(states))
    q = (annihilators[0]+annihilators[0].H)/s.sqrt(2*freq[0])
    # Removing the positive cutoff completion changes the actual source law.
    broken = (s.I*(bare*q-q*bare)-s.I*(h*q-q*h)).applyfunc(s.simplify)
    assert broken != s.zeros(len(states))
    return {'total_cutoff': cap, 'modes': 2, 'complete_dimension': len(states),
            'exact_coordinate_equations': True, 'exact_positive_boundary_completion': True,
            'omitted_boundary_completion_detected': True}


def source_stencil_check():
    from model.native_wilson_maxwell.field import cochains
    vertices, edges, faces, d, c = cochains(1)
    edge_index = {(vertices[tail], axis): i for i, (tail, head, axis) in enumerate(edges)}
    independent = s.zeros(len(faces), len(edges))
    for row, (base, i, j) in enumerate(faces):
        xi = tuple(v+int(axis == i) for axis, v in enumerate(base))
        xj = tuple(v+int(axis == j) for axis, v in enumerate(base))
        for key, sign in (((base, i), 1), ((xi, j), 1), ((xj, i), -1), ((base, j), -1)):
            independent[row, edge_index[key]] += sign
    assert independent == c
    assert c*d == s.zeros(c.rows, d.cols)
    return {'provider': 'model.native_wilson_maxwell.field.cochains',
            'all_six_cube_face_stencils_equal': True,
            'boundary_extension': 'periodic cubic family; not the open-cube interacting ground'}


def run():
    rng = np.random.default_rng(173)
    scalar = rng.normal(size=(5, 5, 5))
    edges = rng.normal(size=(5, 5, 5, 3))
    faces = rng.normal(size=edges.shape)
    identities = {
        'curl_gradient_zero': np.linalg.norm(l.curl(l.gradient(scalar))) < 1e-12,
        'magnetic_divergence_curl_zero': np.linalg.norm(l.magnetic_divergence(l.curl(edges))) < 1e-12,
        'curl_adjoint_pairing': abs(np.vdot(l.curl(edges), faces)-np.vdot(edges, l.curl_adjoint(faces))) < 1e-12,
    }
    assert all(identities.values())
    rows = []
    for n in (8, 12, 16, 24):
        for pol in (0, 1):
            e0, a, q2 = l.transverse_plane(n, polarization=pol)
            w = np.sqrt(q2)
            ce = l.curl(e0, a)
            assert np.linalg.norm(l.curl_adjoint(ce, a)-q2*e0) < 1e-11
            assert np.linalg.norm(l.divergence_adjoint(e0, a)) < 1e-12
            e, b = l.verlet(e0, np.zeros_like(e0), a, 1.25, 128)
            exact_e, exact_b = np.cos(w*1.25)*e0, -np.sin(w*1.25)*ce/w
            error = max(np.linalg.norm(e-exact_e), np.linalg.norm(b-exact_b))
            energy = (np.linalg.norm(e)**2+np.linalg.norm(b)**2)/2
            assert error < 1e-4
            assert np.linalg.norm(l.divergence_adjoint(e, a)) < 1e-11
            assert np.linalg.norm(l.magnetic_divergence(b, a)) < 1e-11
            rows.append({'grid': n, 'polarization': pol, 'lattice_frequency': round(float(w), 12),
                         'continuum_frequency': round(float(np.sqrt(2)), 12),
                         'speed_ratio': round(float(w/np.sqrt(2)), 12),
                         'squared_shell_defect': round(float(2-q2), 12),
                         'analytic_shell_error_upper': round(float(a*a/6), 12),
                         'time_evolution_error': round(float(error), 12),
                         'energy_error': round(float(abs(energy-.5)), 12)})
            assert 0 <= 2-q2 <= a*a/6
    e0, a, q2 = l.transverse_plane(8)
    ce = l.curl(e0, a); w = np.sqrt(q2)
    errors = []
    for steps in (32, 64, 128):
        e, b = l.verlet(e0, np.zeros_like(e0), a, 1.25, steps)
        errors.append(float(max(np.linalg.norm(e-np.cos(w*1.25)*e0),
                                np.linalg.norm(b+np.sin(w*1.25)*ce/w))))
    ratios = [errors[i]/errors[i+1] for i in (0, 1)]
    assert all(3.9 < ratio < 4.1 for ratio in ratios)
    # A longitudinal electric addition violates source-free Gauss even though curl ignores it.
    scalar = rng.normal(size=(8, 8, 8))
    bad = e0+.01*l.gradient(scalar, a)
    assert np.linalg.norm(l.divergence_adjoint(bad, a)) > .1
    e_wrong, b_wrong = l.verlet(e0, np.zeros_like(e0), a, 1.25, 128, speed=1.1)
    assert np.linalg.norm(e_wrong-np.cos(w*1.25)*e0) > .1
    return {'task': 'transverse_waves', 'source_stencil': source_stencil_check(),
            'finite_quantum_shortcut': cutoff_theorem(),
            'identities': {k: bool(v) for k, v in identities.items()},
            'observed_mode_rows': rows, 'time_refinement_error_ratios': [round(v, 9) for v in ratios],
            'checks': {'two_transverse_polarizations': True, 'real_space_matches_exact_propagator': True,
                       'Gauss_and_magnetic_divergence_preserved': True, 'spatial_shell_bound': True,
                       'second_order_time_refinement': True},
            'mutations': {'longitudinal_contamination_rejected': True, 'changed_speed_rejected': True,
                          'missing_cutoff_completion_rejected': True},
            'scope': 'Finite cochain field expectations; a periodic extension of the declared field architecture. The continuum speed is a normalization input. No transfer of interacting open-cube ground or photon particle identity.'}
