"""One fixed current-mediated auxiliary interaction, with every pole retained.

The trace metric of the actual six-dimensional current fixes one invariant
quadratic form. Its matter Casimir couples to the native family defect
projector through two auxiliary hopping bands. Both the coupling and hopping
are fixed to one before comparison. The complete 2250-dimensional operator
has an exact current-image reduction, with all complementary poles retained.

The auxiliary modes, trace metric, common unit and coherent interaction are
additional architecture choices. The bare Higgs coefficient is the proved
stationary averaged feedback kernel; simultaneous finite-memory driving is
not executed here. Numerical poles of this finite operator are not measured
particle masses, and no source-word compilation or detector history is claimed.
"""
from __future__ import annotations

from functools import lru_cache
from fractions import Fraction as F
import hashlib
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy import sparse
import sympy as s

from . import native_matter_attachment as matter
from . import native_repair_mass as feedback
from . import native_repair_mass_attachment as attachment
from . import native_code_transport as transport

ROOT = Path(__file__).resolve().parents[2]
COUPLING = s.Integer(1)
AUXILIARY_HOPPING = s.Integer(1)
NUMERICAL_ADMISSION = 2e-10
SCHEMA = {'dressing_charged_hamiltonians': ((3, 18, 18), 'float64'),
          'dressing_charged_poles': ((3, 18), 'float64'),
          'dressing_charged_residues': ((3, 18, 6, 6), 'float64'),
          'dressing_neutral_hamiltonian': ((9, 9), 'float64'),
          'dressing_neutral_poles': ((9,), 'float64'),
          'dressing_neutral_residues': ((5, 3, 3), 'float64'),
          'dressing_finite_covariance_errors': ((12,), 'float64'),
          'dressing_full_coupling': ((1080, 90), 'complex128'),
          'dressing_full_bare_response': ((90, 90), 'float64')}


def _zero(matrix):
    return all(s.expand(x) == 0 or s.simplify(x) == 0 for x in matrix)


@lru_cache(maxsize=1)
def current_data():
    """Actual port currents, an orthonormal trace frame and its adjoint action."""
    rep = matter.representation(); canonical = matter.canonical_currents()
    frame = tuple([g/s.sqrt(2) for g in canonical[:6]]
                  +[canonical[6]/s.sqrt(2), (canonical[6]+2*canonical[7])/s.sqrt(6)]
                  +[g/s.sqrt(2) for g in canonical[8:11]]+[canonical[11]/s.sqrt(3)])
    gram = s.Matrix(12, 12, lambda a, b: -s.trace(frame[a]*frame[b]))
    if not _zero(gram-s.eye(12)):
        raise ArithmeticError('source trace-frame orthonormality failed')
    coefficients = tuple(matter.native_current.admit_generator(g) for g in frame)
    lifted = tuple(matter.exterior_derivative(matter.spin_lift(g)) for g in frame)
    casimir = matter.clean(-sum((g*g for g in lifted), s.zeros(32)))
    ops = rep['operators']
    expected = 2*ops['color_C2']+ops['weak_C2']/2+ops['q']**2/12
    if not _zero(casimir-expected) or any(not _zero(casimir*g-g*casimir) for g in rep['full_currents']):
        raise ArithmeticError('the same-current trace Casimir identity failed')
    selected = list(rep['selected'])
    small = tuple(s.ImmutableSparseMatrix(g.extract(selected, selected)) for g in lifted)
    c15 = s.ImmutableMatrix(casimir.extract(selected, selected))
    if not c15.is_diagonal() or any(c <= 0 for c in c15.diagonal()):
        raise ArithmeticError('current-image normalization is singular on selected matter')
    adjoints = []
    for native, rho in zip(rep['native']['generators'], rep['matter_currents'], strict=True):
        columns = []
        for g in frame:
            bracket = native*g-g*native
            column = s.Matrix([-s.trace(a*bracket) for a in frame])
            column = matter.clean(column)
            if not _zero(bracket-sum((x*a for x, a in zip(column, frame, strict=True)), s.zeros(6))):
                raise ArithmeticError('an actual port bracket left the trace frame')
            columns.append(column)
        adjoint = s.ImmutableMatrix(s.Matrix.hstack(*columns))
        if not _zero(adjoint+adjoint.T) or any(s.im(x) != 0 for x in adjoint):
            raise ArithmeticError('the auxiliary adjoint action is not real orthogonal')
        for a in range(12):
            residual = rho*small[a]-small[a]*rho+sum(
                (adjoint[a, b]*small[b] for b in range(12)), s.zeros(15))
            if not _zero(residual):
                raise ArithmeticError('current-to-adjoint auxiliary coupling is not an intertwiner')
        adjoints.append(adjoint)
    return {'representation': rep, 'trace_frame': frame, 'port_coefficients': coefficients,
            'matter_frame': small, 'casimir32': casimir, 'casimir15': c15,
            'adjoints': tuple(adjoints)}


def finite_covariance_errors():
    """All twelve actual finite current actions, including auxiliary indices."""
    data = current_data()
    K = np.vstack([np.asarray(g, complex) for g in data['matter_frame']])
    errors = []
    for rho, adjoint in zip(data['representation']['matter_currents'], data['adjoints'], strict=True):
        U = expm(np.asarray(rho, complex)/8)
        A = expm(np.asarray(adjoint, float)/8)
        errors.append(float(np.linalg.norm(np.kron(A, U)@K-K@U, 2)))
    return np.asarray(errors)


@lru_cache(maxsize=1)
def architecture():
    source = feedback.architecture(); Y = source['Y']; Q = s.eye(3)-source['P']
    if not _zero(Q*Q-Q) or s.trace(Q) != 2:
        raise ArithmeticError('the fixed native defect projector changed')
    data = current_data(); rep = data['representation']; selected = list(rep['selected'])
    weak = rep['operators']['weak_C2'].extract(selected, selected)
    c = data['casimir15']; channels = []
    full_internal = s.zeros(15)
    for row in attachment.higgs_kernels(Y):
        kernel = row['internal'][1]; full_internal += kernel
        pairs = [(i, j) for i in range(15) for j in range(i+1, 15) if kernel[i, j] != 0]
        constants = set()
        for i, j in pairs:
            if weak[i, i] == 0: i, j = j, i
            if weak[i, i] != s.Rational(3, 4) or weak[j, j] != 0 or abs(kernel[i, j]) != 1:
                raise ArithmeticError('source Higgs line is not a normalized doublet-singlet pairing')
            constants.add((c[i, i], c[j, j]))
        if len(constants) != 1 or 2*len(pairs) != row['neutral_rank']:
            raise ArithmeticError('the current-derived sector constants are not uniform')
        left, right = constants.pop()
        channels.append({'name': row['channel'], 'left': left, 'right': right,
                         'multiplicity': row['neutral_rank']})
    channels.sort(key=lambda row: ('up', 'down', 'charged_lepton').index(row['name']))
    unpaired = [i for i in range(15) if all(full_internal[i, j] == 0 for j in range(15))]
    if len(unpaired) != 1 or sum(row['multiplicity'] for row in channels) != 14:
        raise ArithmeticError('complete neutral Higgs domain changed')
    return {'Y': Y, 'Q': s.ImmutableMatrix(Q), 'channels': tuple(channels),
            'neutral_casimir': c[unpaired[0], unpaired[0]]}


def charged_block(left, right):
    """An exact18-state current-image block; four auxiliary dark states remain.

    Chiral plus coordinates are(L3, bL3, aR3), minus(R3, aL3, bR3).
    Auxiliary a/b hopping is one. Its symmetric/antisymmetric combinations
    give the two signs used in the full2250-domain decomposition.
    """
    a = architecture(); D = s.zeros(9)
    D[:3, :3] = a['Y']; D[:3, 3:6] = s.sqrt(left)*a['Q']
    D[3:6, 3:6] = s.eye(3); D[6:9, :3] = s.sqrt(right)*a['Q']
    D[6:9, 6:9] = s.eye(3)
    if s.simplify(D.det()-a['Y'].det()) != 0:
        raise ArithmeticError('auxiliary hopping changed the exact determinant identity')
    return s.ImmutableMatrix(s.BlockMatrix([[s.zeros(9), D], [D.T, s.zeros(9)]]).as_explicit())


def neutral_block():
    a = architecture(); H = s.zeros(9); coupling = s.sqrt(a['neutral_casimir']/2)*a['Q']
    H[:3, 3:6] = coupling; H[3:6, :3] = coupling
    H[:3, 6:9] = coupling; H[6:9, :3] = coupling
    H[3:6, 3:6] = s.eye(3); H[6:9, 6:9] = -s.eye(3)
    return s.ImmutableMatrix(H)


@lru_cache(maxsize=1)
def full_factors():
    """All2250 input/output columns are specified by these unprojected factors."""
    a, data = architecture(), current_data()
    internal = sum((row['internal'][1] for row in attachment.higgs_kernels(a['Y'])), s.zeros(15))
    kernel = np.asarray(s.kronecker_product(internal, a['Y']), float)
    bare = np.block([[np.zeros((45, 45)), kernel], [kernel.T, np.zeros((45, 45))]])
    Q = np.asarray(a['Q'], float)
    coupling = []
    for rho in data['matter_frame']:
        rho = np.asarray(rho, complex)
        nambu = np.block([[rho, np.zeros((15, 15))], [np.zeros((15, 15)), rho.conj()]])
        coupling.append(np.kron(nambu, Q))
    K = np.vstack(coupling)
    c = np.asarray(data['casimir15'], float)
    expected = np.kron(np.kron(np.eye(2), c), Q)
    if np.max(abs(K.conj().T@K-expected)) > NUMERICAL_ADMISSION:
        raise ArithmeticError('complete current-family coupling lost its Gram identity')
    return bare, K


def full_sparse_operator():
    """The actual2250-dimensional Hermitian response, including dark states."""
    bare, K = full_factors(); V = sparse.csr_matrix(K/np.sqrt(2))
    eye = sparse.eye(1080, format='csr')
    return sparse.bmat([[sparse.csr_matrix(bare), V.conj().T, V.conj().T],
                        [V, eye, None], [V, None, -eye]], format='csr')


def shape_obstruction():
    """Exclude a large middle/small gap ratio for every nonnegative dressing.

    For the fixed positive Y, let alpha=Tr(PY), S=Y-YPY/alpha. Then S>=0,
    SP=0 and Y-S has rank one. With A_r(E)=P+a_r(E)Q, 0<a_r<=1,
    D(E)=A_L Y A_R is a rank-one matrix plus a_L*a_R*S. The singular-value
    min-max inequalities give sigma_2(D)<=a_L*a_R*Tr(S) and
    sigma_1(D)>=a_L*a_R*lambda_min(Y), in ascending three-value order.

    For a common Delta>0, a_r(E)=sqrt((Delta^2-E^2)/(Delta^2+c_r-E^2))
    decreases on0<=E<Delta for c_r>=0.
    Left/right contraction therefore makes every singular value nonincreasing.
    Any existing exact Feshbach root E_j=sigma_j(D(E_j)) is unique on its
    ordered branch. Conditional on two positive poles below Delta,
    E2<=sigma_2(D(E1))<32E1. Three below-threshold poles are not asserted
    for every Delta. This is an exact frequency-dependent
    argument, not the constant-Z approximation or a parameter scan.
    """
    a = architecture(); Y, P = a['Y'], s.eye(3)-a['Q']
    alpha = s.simplify(s.trace(P*Y)); S = (Y-Y*P*Y/alpha).applyfunc(s.simplify)
    if not _zero(S*P) or S.rank() != 2 or (Y-S).rank() != 1:
        raise ArithmeticError('the source rank-one Schur decomposition changed')
    # The fixed-source spectral certificate proves Y>0. Cauchy--Schwarz in
    # that positive form proves S>=0; retain exact principal-minor controls.
    import itertools
    for size in (1, 2, 3):
        for indices in itertools.combinations(range(3), size):
            if s.simplify(S.extract(indices, indices).det()).is_nonnegative is not True:
                raise ArithmeticError('a Schur-complement principal minor is negative')
    cert = attachment.spectral_certificate(Y)
    smallest = F(cert['positive_eigenvalue_intervals_exact'][0][0])
    sqrt_upper = F(cert['sqrt5_rational_interval'][1])
    trace = s.expand(s.simplify(s.trace(S)))
    rational, radical = trace.coeff(s.sqrt(5), 0), trace.coeff(s.sqrt(5), 1)
    if rational.is_Rational is not True or radical.is_Rational is not True or radical < 0:
        raise ArithmeticError('source Schur-trace rational enclosure changed')
    trace_upper = F(rational)+F(radical)*sqrt_upper
    bound = trace_upper/smallest
    if not bound < 32:
        raise ArithmeticError('the exact no-tuning source bound changed')
    return {'subject': 'two smallest positive exact finite-response poles below a common auxiliary threshold Delta',
            'fixed_bare_source_Y': True,
            'allowed_dressing_constants': 'arbitrary c_left>=0,c_right>=0,Delta>0 multiplying the same fixed Qf',
            'existence_premise': 'at least two positive exact poles satisfy0<E1<=E2<Delta',
            'three_poles_below_threshold_asserted_for_every_Delta': False,
            'subthreshold_poles_have_nonzero_light_component': True,
            'Schur_trace_exact': str(trace), 'source_smallest_eigenvalue_lower_exact': str(smallest),
            'Schur_trace_rational_upper_exact': str(trace_upper),
            'middle_over_smallest_ratio_upper_exact': str(bound),
            'middle_over_smallest_ratio_strictly_below_32': True,
            'proof': ['Y=YPY/Tr(PY)+S; S>=0, SP=0, rank(Y-S)=1',
                      'D(E)=A_left(E)Y A_right(E); sigma2/sigma1<=Tr(S)/lambda_min(Y)',
                      'A_r(E)=P+sqrt((Delta^2-E^2)/(Delta^2+c_r-E^2))*Q contracts monotonically',
                      'E_j=sigma_j(D(E_j)); E2<=sigma2(D(E1))<32E1'],
            'parameter_scan_or_low_energy_kinetic_approximation_used': False,
            'above_threshold_reassignment_excluded': False,
            'above_threshold_reassignment_removes_lower_charged_states': False,
            'other_family_projectors_or_interactions_excluded': False}


def _spectrum(H, light_rows):
    values, vectors = np.linalg.eigh(H)
    residual = float(np.linalg.norm(H@vectors-vectors*values, 2))
    if residual > NUMERICAL_ADMISSION:
        raise ArithmeticError('complete auxiliary eigenpair admission failed')
    projected = vectors[list(light_rows)]
    residues = np.einsum('ik,jk->kij', projected, projected)
    if np.max(abs(residues.sum(axis=0)-np.eye(len(light_rows)))) > NUMERICAL_ADMISSION:
        raise ArithmeticError('complete spectral residues do not resolve the light identity')
    return values, residues


def execute():
    a = architecture(); hs, poles, residues = [], [], []
    for row in a['channels']:
        H = np.asarray(charged_block(row['left'], row['right']), float)
        values, rs = _spectrum(H, tuple(range(3))+tuple(range(9, 12)))
        hs.append(H); poles.append(values); residues.append(rs)
    neutral = np.asarray(neutral_block(), float)
    values, rs = _spectrum(neutral, range(3))
    # Exact degeneracies of the neutral block:2,1,3,1,2. Aggregate complete
    # eigenspaces so arbitrary numerical bases inside them are not evidence.
    groups = ((0, 1), (2,), (3, 4, 5), (6,), (7, 8))
    arrays = {'dressing_charged_hamiltonians': np.asarray(hs),
              'dressing_charged_poles': np.asarray(poles),
              'dressing_charged_residues': np.asarray(residues),
              'dressing_neutral_hamiltonian': neutral, 'dressing_neutral_poles': values,
              'dressing_neutral_residues': np.asarray([rs[list(group)].sum(axis=0) for group in groups]),
              'dressing_finite_covariance_errors': finite_covariance_errors(),
              'dressing_full_bare_response': full_factors()[0].copy(),
              'dressing_full_coupling': full_factors()[1].copy()}
    validate_arrays(arrays)
    return _report(arrays), arrays


def validate_arrays(arrays):
    if set(arrays) != set(SCHEMA): raise ValueError('all auxiliary source arrays required')
    for name, (shape, dtype) in SCHEMA.items():
        value = arrays[name]
        if (type(value) is not np.ndarray or value.shape != shape or value.dtype != np.dtype(dtype)
                or not np.isfinite(value).all()):
            raise ValueError('invalid auxiliary source array '+name)
    if np.max(arrays['dressing_finite_covariance_errors']) > NUMERICAL_ADMISSION:
        raise ValueError('finite current covariance failed')
    bare, K = full_factors()
    if (np.max(abs(bare-arrays['dressing_full_bare_response'])) > NUMERICAL_ADMISSION
            or np.max(abs(K-arrays['dressing_full_coupling'])) > NUMERICAL_ADMISSION):
        raise ValueError('complete2250 source factors changed')
    for i, row in enumerate(architecture()['channels']):
        expected = np.asarray(charged_block(row['left'], row['right']), float)
        if np.max(abs(expected-arrays['dressing_charged_hamiltonians'][i])) > NUMERICAL_ADMISSION:
            raise ValueError('a primitive auxiliary operator changed')
        poles, residues = _spectrum(expected, tuple(range(3))+tuple(range(9, 12)))
        if (np.max(abs(poles-arrays['dressing_charged_poles'][i])) > NUMERICAL_ADMISSION
                or np.max(abs(residues-arrays['dressing_charged_residues'][i])) > NUMERICAL_ADMISSION):
            raise ValueError('complete charged poles or residues changed')
    neutral = np.asarray(neutral_block(), float)
    values, residues = _spectrum(neutral, range(3))
    groups = ((0, 1), (2,), (3, 4, 5), (6,), (7, 8))
    grouped = np.asarray([residues[list(g)].sum(axis=0) for g in groups])
    if (np.max(abs(neutral-arrays['dressing_neutral_hamiltonian'])) > NUMERICAL_ADMISSION
            or np.max(abs(values-arrays['dressing_neutral_poles'])) > NUMERICAL_ADMISSION
            or np.max(abs(grouped-arrays['dressing_neutral_residues'])) > NUMERICAL_ADMISSION):
        raise ValueError('complete neutral operator, poles or residues changed')


def _report(arrays):
    a = architecture(); rows = []
    for i, row in enumerate(a['channels']):
        poles = arrays['dressing_charged_poles'][i]
        rs = arrays['dressing_charged_residues'][i]
        positive = np.flatnonzero(poles > NUMERICAL_ADMISSION)
        light = np.flatnonzero((poles > NUMERICAL_ADMISSION)&(poles < 1-NUMERICAL_ADMISSION))
        dark = np.flatnonzero(abs(abs(poles)-1) < NUMERICAL_ADMISSION)
        if len(positive) != 9 or len(light) != 3 or len(dark) != 4:
            raise ArithmeticError('frozen light/auxiliary pole count changed')
        if np.max(abs(rs[dark])) > NUMERICAL_ADMISSION:
            raise ArithmeticError('declared auxiliary dark pole is visible on light inputs')
        energy = poles[light]
        intervals = [[float((x-NUMERICAL_ADMISSION)/(energy[0]+NUMERICAL_ADMISSION)),
                      float((x+NUMERICAL_ADMISSION)/(energy[0]-NUMERICAL_ADMISSION))] for x in energy[1:]]
        rows.append({'Higgs_channel': row['name'], 'left_Casimir_exact': str(row['left']),
            'right_Casimir_exact': str(row['right']), 'complete_Nambu_color_multiplicity': row['multiplicity'],
            'positive_poles': list(map(float, poles[positive])),
            'positive_pole_light_residue_traces': list(map(float, np.trace(rs[positive], axis1=1, axis2=2))),
            'light_window_declared_before_comparison': '0<E<auxiliary_hopping=1',
            'light_poles': list(map(float, energy)), 'light_ratio_numerical_intervals': intervals,
            'dark_poles_plus_minus_one_each_multiplicity': 2,
            'four_additional_positive_bright_auxiliary_poles_retained': True})
    return {'classification': 'fixed finite charged auxiliary response with complete poles',
        'common_coupling_exact': '1', 'common_auxiliary_hopping_exact': '1',
        'trace_metric': '-Tr_M6(D_a D_b); an explicitly chosen invariant form',
        'Casimir_exact': '2*color_C2+weak_C2/2+q^2/12',
        'native_current_auxiliary_intertwining_identities': 144,
        'finite_actual_port_current_actions_checked': 12,
        'family_coupling': 'Qf=I-P0 from the fixed native port seed',
        'complete_domain': {'light_Nambu_matter_family': 90, 'auxiliary_adjoint_indices': 12,
            'internal_domain': 'selected exterior M15; the other17 source states remain outside this declared response domain',
            'auxiliary_signs': 2, 'total_dimension': 2250,
            'current_image_reduction_dimension': 270, 'current_orthogonal_auxiliary_dimension': 1980,
            'family_dark_dimension_inside_reduction': 60, 'bright_dimension': 210,
            'all_dark_poles': [{'energy_exact': '-1', 'multiplicity': 1020},
                               {'energy_exact': '1', 'multiplicity': 1020}],
            'full_operator_zero_multiplicity': 6,
            'all2250_input_output_columns_specified_by_retained_factors': True,
            'auxiliary_representation': 'two copies of adjoint12 tensor Nambu(M15) tensor neutral family3',
            'dark_states_absent_from_the_complete_operator': False},
        'exact_reduction': {'current_map': 'K psi=(rho(E_a) psi)_a; K* K=C_R',
            'current_image_isometry': 'J=K C_R^(-1/2); J*J=I on all15 selected internal states',
            'auxiliary_current': 'Ad(E) on the current index tensor rho(E) on matter',
            'family_complement_retained': 'Qf rank2, I-Qf rank1; both auxiliary signs retained',
            'light_self_energy': 'Sigma_R(E)=C_R Qf E/(E^2-1)',
            'low_energy_kinetic_matrix': 'Z_R=I+C_R Qf',
            'self_energy_remainder_upper': 'C_R*e^3/(1-e^2) for |E|<=e<1',
            'kinetic_approximation_used_as_reported_poles': False},
        'charged_channels': rows,
        'exact_no_tuning_shape_bound': shape_obstruction(),
        'neutral_channel': {'full_Nambu_multiplicity': 2, 'Casimir_exact': str(a['neutral_casimir']),
            'all9_poles': list(map(float, arrays['dressing_neutral_poles'])),
            'grouped_pole_multiplicities': [2, 1, 3, 1, 2],
            'light_zero_modes_per_Nambu_copy': 3},
        'numerical_operator_admission': NUMERICAL_ADMISSION,
        'scope': {'auxiliary_rule_trace_form_and_common_unit_supplied': True,
            'bare_Y_uses_stationary_averaging_limit': True,
            'simultaneous_finite_memory_feedback_and_dressing_executed': False,
            'complete2250_operator_compiled_into_native_source_words': False,
            'poles_inferred_from_detector_records_here': False,
            'auxiliary_Weyl_kinetics_and_physical_spin_attached': False,
            'full_physical_spin_particle_hole_history_derived': False,
            'finite_bipartite_sign_symmetry_is_physical_particle_hole_theorem': False,
            'measured_mass_or_pole_values_used_as_couplings': False,
            'additional_auxiliary_particles_removed': False,
            'physical_interacting_particle_spectrum_derived': False,
            'numerical_intervals_are_outward_physical_certificates': False}}


def report():
    return execute()[0]


def report_from_arrays(arrays):
    """Bind numerical pole/ratio diagnostics to their own complete source arrays."""
    validate_arrays(arrays)
    return _report(arrays)


def source_bindings():
    paths = set(matter.source_bindings()) | set(transport.source_bindings()) | {
        'model/common_source/native_repair_mass.py', 'model/common_source/test_native_repair_mass.py',
        'model/common_source/native_repair_mass_attachment.py',
        'model/common_source/test_native_repair_mass_attachment.py',
        'model/common_source/native_mass_dressing.py', 'model/common_source/test_native_mass_dressing.py'}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}
