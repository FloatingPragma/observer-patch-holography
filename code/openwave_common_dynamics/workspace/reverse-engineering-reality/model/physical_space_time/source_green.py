"""Controlled occupation Green operator on the existing joined source screens.

Cell routing is a declared input. This module eliminates actual within-carrier
variables; it does not substitute a geometric mesh Laplacian for native repair.
"""
from __future__ import annotations

from fractions import Fraction as F

import numpy as np
import sympy as sp
from scipy import sparse
from scipy.sparse.linalg import eigsh, spsolve

from ..federated_observer import source
from . import construction as passage

SOURCE_PATHS = tuple(sorted(set((*source.SOURCE_PATHS, *passage.SOURCE_PATHS,
    'model/physical_space_time/source_green.py',
    'model/federated_observer/__init__.py',
    'model/anharmonic_transfer/__init__.py', 'model/anharmonic_transfer/environment.py',
    'model/common_source/protocol.py', 'model/common_source/smooth_protocol.py',
    'model/primitive_observer/__init__.py', 'model/relativistic_scalar/__init__.py',
    'model/relativistic_scalar/lattice.py', 'model/scalar_observations/__init__.py',
    'model/scalar_observations/experimental_readout/__init__.py',
    'model/scalar_observations/experimental_readout/canonical_source.py',
    'model/scalar_observations/experimental_readout/executable_source.py',
    'model/scalar_observations/lattice.py', 'model/scalar_observations/source.py'))))


def matrices(level):
    """Build exact integer blocks and check every source seam and matching."""
    graph = source.federation(level)
    _, local = passage.recovered_carrier()
    green = (local+sp.ones(12)/12).inv()-sp.ones(12)/12
    denominator = int(sp.ilcm(*(entry.q for entry in green)))
    green_integer = np.asarray(green*denominator, dtype=np.int64)
    cross = graph.inter_seams
    endpoints = [port for edge in cross for port in edge]
    if len(endpoints) != len(set(endpoints)):
        raise ArithmeticError('source cross seams are not a port matching')
    if any(a//12 == b//12 for a, b in cross):
        raise ArithmeticError('an intercarrier seam is internal')
    rows, columns, values = [], [], []
    for j, (a, b) in enumerate(cross):
        rows.extend((a, b)); columns.extend((j, j)); values.extend((1, -1))
    d = sparse.csc_matrix((values, (rows, columns)), shape=(graph.ports, len(cross)), dtype=np.int64)
    if (d.T@d-2*sparse.eye(len(cross), dtype=np.int64)).nnz:
        raise ArithmeticError('matching incidence norm certificate failed')
    # V repeats a cell load without the isometry factor 1/sqrt(12).
    v = sparse.csc_matrix((np.ones(graph.ports, dtype=np.int64),
        (np.arange(graph.ports), np.arange(graph.ports)//12)), shape=(graph.ports, graph.cells))
    r = (d.T@v).toarray()
    cell_l = r.T@r
    if not np.array_equal(np.diag(cell_l), np.full(graph.cells, 3)):
        raise ArithmeticError('the actual source cell graph is not degree three')
    g0_integer = sparse.kron(sparse.eye(graph.cells, dtype=np.int64), green_integer, format='csc')
    f_integer = (denominator*sparse.eye(len(cross), dtype=np.int64)+d.T@g0_integer@d).toarray()
    l = sparse.kron(sparse.eye(graph.cells, dtype=np.int64),
                    np.asarray(local, dtype=np.int64), format='csc')+d@d.T
    # Check the full port Laplacian against the original neighbor table.
    for port, neighbors in enumerate(graph.neighbours):
        row = l.getrow(port)
        expected = {port: len(neighbors), **{other: -1 for other in neighbors}}
        if dict(zip(row.indices, row.data, strict=True)) != expected:
            raise ArithmeticError('block elimination is not using the actual source graph')
    return graph, l, r, f_integer, denominator


def approximation_relative_bound(order):
    if type(order) is not int or not 1 <= order <= 32:
        raise ValueError('an exact approximation order from one through thirty-two is required')
    return F(1, 2**(order+1)-1)


def schur_series(f, r, order):
    """Numerical evaluation of the analytically bounded exact polynomial."""
    approximation_relative_bound(order)
    t = np.eye(len(f))-f/2
    term, result = r.astype(float), r.astype(float)/2
    for _ in range(order):
        term = t@term
        result += term/2
    return (r.T@result)/12


def _green(matrix):
    cells = len(matrix)
    constant = np.ones((cells, cells))/cells
    return np.linalg.solve(matrix+constant, np.eye(cells))-constant


def diagnostic(level):
    """Actual finite screens; floating evaluations are diagnostics, not proofs."""
    graph, l, r, integer_f, denominator = matrices(level)
    f = integer_f/denominator
    a = r.T@r/12
    s = (r.T@np.linalg.solve(f, r))/12
    exact_green, cell_green = _green(s), _green(a)
    gap = 5-np.sqrt(5)
    factor = (gap+2)/gap
    if (np.linalg.eigvalsh(exact_green-cell_green)[0] < -2e-9
            or np.linalg.eigvalsh(factor*cell_green-exact_green)[0] < -2e-9):
        raise ArithmeticError('evaluated source Green operators violate the analytic comparison')
    # An independent full-port Dirichlet solve audits a fixed cell contrast.
    readout = np.zeros(graph.cells); readout[0], readout[1] = 1/np.sqrt(2), -1/np.sqrt(2)
    port_readout = np.repeat(readout/np.sqrt(12), 12)
    solution = np.zeros(graph.ports)
    solution[1:] = spsolve(l[1:, 1:].astype(float), port_readout[1:])
    solution -= solution.mean()
    residual = float(np.linalg.norm(l@solution-port_readout))
    direct = float(port_readout@solution)
    schur = float(readout@exact_green@readout)
    if residual > 2e-10 or abs(direct-schur) > 2e-9:
        raise ArithmeticError('full native Green readout differs from carrier elimination')
    incident = tuple((a, b) for a, b in graph.seams if a//12 == 0 or b//12 == 0)
    if len(incident) != 33:
        raise ArithmeticError('the root-cell incident clock changed its seam count')
    local_green_pairing = sum((port_readout[a]-port_readout[b])*(solution[a]-solution[b])
                              for a, b in incident)
    incident_correction = (1-2*local_green_pairing)/(4*len(incident))
    attempt_correction = -1/(4*len(graph.seams))
    rows = []
    for order in (4, 8, 16):
        approximate = _green(schur_series(f, r, order))
        relative = float(np.linalg.norm(approximate-exact_green, 2)/np.linalg.norm(exact_green, 2))
        if relative > float(approximation_relative_bound(order))+2e-10:
            raise ArithmeticError('source Schur polynomial exceeded its relative allowance')
        rows.append({'order': order, 'relative_Green_error_bound': str(approximation_relative_bound(order)),
                     'relative_Green_error_diagnostic': relative})
    values, vectors = eigsh(l.astype(float), k=4, which='SM',
        v0=np.linspace(1., 2., graph.ports), tol=2e-11, maxiter=40*graph.ports)
    order = np.argsort(values); values, vectors = values[order], vectors[:, order]
    eigen_residual = float(np.linalg.norm(l@vectors-vectors*values))
    if eigen_residual > 2e-9 or abs(values[0]) > 1e-9:
        raise ArithmeticError('native low-mode diagnostic failed its residual check')
    modes = []
    for index in range(1, 4):
        vector = vectors[:, index].reshape(graph.cells, 12)
        within = float(np.sum((vector-vector.mean(axis=1, keepdims=True))**2))
        allowance = float(values[index]/gap)
        if within > allowance+2e-9:
            raise ArithmeticError('source low mode exceeds the internal-gap allowance')
        modes.append({'native_L_eigenvalue_diagnostic': float(values[index]),
            'within_cell_nonconstant_weight_diagnostic': within,
            'weight_bound_evaluated_at_diagnostic_eigenvalue': allowance})
    return {'source_level': level, 'carriers': graph.cells, 'ports': graph.ports,
        'intercarrier_seams': len(graph.inter_seams), 'all_intercarrier_endpoints_distinct': True,
        'exact_carrier_green_denominator': denominator,
        'exact_matching_and_all_source_laplacian_rows_checked': True,
        'cell_contrast_native_Green_diagnostic': direct,
        'cell_contrast_Schur_Green_diagnostic': schur,
        'same_history_clock_fields': {
            'cell_readout': '(cell0 total-cell1 total)/sqrt(24)',
            'canonical_total': graph.ports//3,
            'full_initial_occupation_and_seam_current_disclosure_supplied': True,
            'covariance_normalization': 'divide CLT covariance by 4*chi*clock_seam_count',
            'attempt_covariance_diagnostic': direct+attempt_correction,
            'root_cell_incident_covariance_diagnostic': float(direct+incident_correction),
            'root_cell_incident_clock_seams': len(incident),
            'incident_clock_correction_diagnostic': float(incident_correction),
            'readout_identified_across_different_screen_levels': False},
        'full_port_solve_residual_diagnostic': residual,
        'controlled_Schur_approximations': rows,
        'low_global_modes': modes, 'eigenpair_residual_diagnostic': eigen_residual,
        'continuum_target_or_physical_length_used': False}


def run():
    native_audit = source.local_audit()
    if not native_audit['diagonal_transition_is_half_identity_plus_swap']:
        raise ArithmeticError('source Green law lost the actual native exchange identity')
    return {'actual_native_two_endpoint_cases': len(native_audit['complete_local_cases']),
        'actual_source_levels': [diagnostic(level) for level in (0, 1, 2)],
        'exact_theorem': {
            'internal_nonconstant_gap': '5-sqrt(5)',
            'cross_matching_operator_norm_bound': '2',
            'cell_isometry': 'U u repeats u_c/sqrt(12) on each carrier',
            'coarse_energy': 'A=L_cell/12',
            'eliminated_energy': 'S=(1/12) R^T (I+D^T L_internal^+ D)^(-1) R',
            'source_Green': 'U^T L^+ U=S^+',
            'Loewner_Green_comparison': 'A^+ <= S^+ <= ((7-sqrt(5))/(5-sqrt(5))) A^+',
            'low_mode_within_carrier_weight_bound': '||Q psi||^2 <= omega/(5-sqrt(5)) for ||psi||=1',
            'series_order_p_relative_Green_error_bound': '1/(2^(p+1)-1), p>=1',
            'attempt_average_cell_field_covariance_after_normalization': 'S^+ - P_cell/(4M)',
            'same_native_endpoint_law': True,
            'matching_premise_required_for_factor_two': True},
        'acceptance': {'actual_source_eliminated_Green_operator_and_uniform_error_bound_proved': True,
            'source_target_manifold_Green_convergence_proved': False,
            'rank_three_internal_band_is_three_dimensional_bulk': False,
            'Lorentzian_evolution_or_physical_radial_coordinate_derived': False}}
