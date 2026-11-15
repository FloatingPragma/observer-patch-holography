"""Same-current Higgs attachment and held-out test of fixed repair feedback.

The common family coefficient is supplied by the separate native repair
producer. This module neither selects its record probability nor fits a
matrix to empirical masses. The physical Weyl kinetic and Higgs attachment
are declared architecture; their internal current identities are exact.
"""
from __future__ import annotations

from fractions import Fraction as F
import ast
from functools import lru_cache
import hashlib
import itertools
import json

import sympy as s

from tools import strict_json
from . import native_matter_attachment as matter
from . import native_gauge_object as objects, slow_geometry
from .particle_identification import family_attachment

ROOT = matter.ROOT
BENCHMARK = 'model/common_source/chiral_candidate/results.json'
E2 = s.Rational(1392, 22621)
DET = s.Rational(100224, 511709641)-s.Rational(418176, 12792741025)*s.sqrt(5)
ROOT_BOXES = ((2064153515, 2064153516), (63659825356, 63659825357),
              (934276021128, 934276021129))
ROOT_DENOMINATOR = 10**12
FAMILY_THEOREM = 'Lean/Screen/A5FamilyBand.lean'
PORT_ACTION_THEOREM = 'Lean/Screen/A5PortAction.lean'
FAMILY_MANIFEST = 'code/a5_closure/manifests/echosahedral_federation_reference.json'


def _zero(matrix):
    return all(s.expand(value) == 0 or s.simplify(value) == 0 for value in matrix)


def _family(value):
    if (not isinstance(value, s.MatrixBase) or value.shape != (3, 3)
            or any(x.has(s.Float) or x.free_symbols or x.is_real is not True
                   or x.is_finite is not True for x in value)
            or not _zero(value-value.T)):
        raise ValueError('the exact real symmetric three-family source coefficient is required')
    return s.ImmutableSparseMatrix(value)


def _lean_list(path, name):
    """Read only a literal integer-list definition from the pinned theorem."""
    text = (ROOT/path).read_text(encoding='utf-8')
    start = text.index('[', text.index(':=', text.index('def '+name+' '))+2)
    depth = 0
    for end in range(start, len(text)):
        depth += (text[end] == '[')-(text[end] == ']')
        if depth == 0:
            return ast.literal_eval(text[start:end+1])
    raise ValueError('the pinned Lean list is incomplete')


@lru_cache(maxsize=1)
def family_band_certificate():
    """Bind the neutral coefficient copy to the actual finite family theorem."""
    data = objects.geometry(); B = data['frame']; metric = s.expand((B.T*B)[0, 0])
    if not _zero(B.T*B-metric*s.eye(3)):
        raise ArithmeticError('the native coefficient frame is not isotropic')
    adjacency = s.zeros(12)
    for i, j in data['seams']: adjacency[i, j] = adjacency[j, i] = 1
    manifest = strict_json.load(ROOT/FAMILY_MANIFEST)['carrier']
    native_to_carrier = family_attachment._chart_map(data['seams'], slow_geometry.source().oriented_faces(), manifest)
    carrier_to_lean = _lean_list(FAMILY_THEOREM, 'carrierToLean')
    chart = tuple(carrier_to_lean[i] for i in native_to_carrier)
    inverse = tuple(chart.index(i) for i in range(12))
    if adjacency != s.Matrix(_lean_list(FAMILY_THEOREM, 'At')).extract(chart, chart):
        raise ArithmeticError('the native family chart does not bind to Lean incidence')
    X = s.Matrix(_lean_list(FAMILY_THEOREM, 'Xt')).extract(chart, chart)
    Y = s.Matrix(_lean_list(FAMILY_THEOREM, 'Yt')).extract(chart, chart)
    projectors = (s.ones(12)/12, (X+s.sqrt(5)*Y)/40, (X-s.sqrt(5)*Y)/40,
                  s.Matrix(_lean_list(FAMILY_THEOREM, 'M5t')).extract(chart, chart)/24)
    P = objects._clean(B*B.T/metric)
    if not _zero(adjacency*B-s.sqrt(5)*B) or not _zero(P-projectors[1]):
        raise ArithmeticError('the native family isometry differs from the exact Lean P3')
    if not _zero(sum(projectors, s.zeros(12))-s.eye(12)):
        raise ArithmeticError('the full12 coefficient bands do not resolve identity')
    for projector, rank, value in zip(projectors, (1, 3, 3, 5), (5, s.sqrt(5), -s.sqrt(5), -1), strict=True):
        if s.trace(projector) != rank or not _zero(projector*projector-projector) or not _zero(adjacency*projector-value*projector):
            raise ArithmeticError('an unselected native coefficient band changed')
    listed = set(map(tuple, _lean_list(PORT_ACTION_THEOREM, 'perms')))
    mapped = set(); kernels = 0
    for permutation in data['rotations']:
        mapped.add(tuple(chart[permutation[inverse[i]]] for i in range(12)))
        action = objects.permutation_matrix(permutation)
        rotation = objects._clean(B.T*action*B/metric)
        if not _zero(action*B-B*rotation) or not _zero(rotation.T*rotation-s.eye(3)):
            raise ArithmeticError('a native rotation does not intertwine the family coefficient isometry')
        kernels += _zero(rotation-s.eye(3))
    if mapped != listed or len(mapped) != 60 or kernels != 1:
        raise ArithmeticError('the actual faithful family rotation action differs from the Lean table')
    if not (1+s.sqrt(5)).is_positive or not (s.sqrt(5)-1).is_positive:
        raise ArithmeticError('the faithful in-window seam-cost order is not strict')
    V = B/s.sqrt(metric)
    return {'native_to_family_manifest_chart': list(native_to_carrier), 'native_to_Lean_chart': list(chart),
            'normalized_coefficient_isometry12_by3_exact': [[str(x) for x in row] for row in V.tolist()],
            'frame_metric_exact': str(metric), 'V_transpose_V_identity_exact': True,
            'adjacency_V_equals_sqrt5_V_exact': True, 'V_Vtranspose_equals_Lean_P3_exact': True,
            'actual_rotation_intertwining_count': 60, 'actual_family_rotation_kernel_order': kernels,
            'full_coefficient_band_dimensions_retained': [1, 3, 3, 5],
            'full_coefficient_adjacency_eigenvalues_exact': ['5', 'sqrt(5)', '-sqrt(5)', '-1'],
            'faithful_in_window_cost_order_exact': ['5-sqrt(5)', '6', '5+sqrt(5)'],
            'selected_band_tensor_M15_dimension': 45, 'unfiltered_coefficient_tensor_M15_dimension': 180,
            'separate_neutral_family_copy_is_the_colored_current_triplet': False,
            'coefficient_isometry_is_a_native_load_two_state_embedding': False,
            'single_band_physical_realization_and_operational_selection_rule_are_architecture': True,
            'physical_family_number_forced_by_gauge_current_alone': False,
            'finite_native_filter_history_for_all12_coefficient_inputs_executed_here': False,
            'theorem_connections': [
                'OPH.A5FamilyBand.band3_eigen', 'OPH.A5FamilyBand.band_equivariant',
                'OPH.A5FamilyBand.resolution_complete', 'OPH.A5FamilyBand.family_band_selected',
                'SMThreeFamilyPremiseData.band_realized', 'SMThreeFamilyPremiseData.band_faithful',
                'SMThreeFamilyPremiseData.band_window', 'SMThreeFamilyPremiseData.band_minimal']}


def spectral_certificate(family):
    """Isolate all three roots with rational arithmetic and a sqrt5 bracket."""
    family = _family(family)
    e2 = s.simplify((s.trace(family)**2-s.trace(family*family))/2)
    determinant = s.simplify(family.det())
    if not _zero(s.Matrix([s.trace(family)-1, e2-E2, determinant-DET])):
        raise ArithmeticError('the frozen source-feedback characteristic polynomial differs')
    lower = F(2236067977499789696409173668731, 10**30)
    upper = lower+F(1, 10**30)
    if not lower*lower < 5 < upper*upper:
        raise ArithmeticError('the rational square-root enclosure failed')

    def signs(x):
        rational = x**3-x*x+F(1392, 22621)*x-F(100224, 511709641)
        return (rational+F(418176, 12792741025)*lower,
                rational+F(418176, 12792741025)*upper)

    intervals = [tuple(F(x, ROOT_DENOMINATOR) for x in row) for row in ROOT_BOXES]
    for lo, hi in intervals:
        left, right = signs(lo), signs(hi)
        if not ((left[1] < 0 < right[0]) or (right[1] < 0 < left[0])):
            raise ArithmeticError('a frozen positive source root lost its sign bracket')
    if not 0 < intervals[0][0] < intervals[0][1] < intervals[1][0] < intervals[1][1] < intervals[2][0] < intervals[2][1] < 1:
        raise ArithmeticError('the three positive root brackets are not disjoint')
    ratios = [(lo/intervals[0][1], hi/intervals[0][0]) for lo, hi in intervals[1:]]
    return {'characteristic_polynomial': 'x^3-x^2+e2*x-determinant',
            'trace_exact': '1', 'e2_exact': str(E2), 'determinant_exact': str(DET),
            'sqrt5_rational_interval': list(map(str, (lower, upper))),
            'positive_eigenvalue_intervals_exact': [list(map(str, row)) for row in intervals],
            'ascending_ratios_after_one_exact': [list(map(str, row)) for row in ratios],
            'ascending_ratios_after_one_diagnostic': [[float(x) for x in row] for row in ratios],
            'three_sign_changes_exhaust_the_degree_three_polynomial': True,
            'empirical_inputs': []}


def higgs_kernels(family):
    """The same real symmetric Y on every actually derived Higgs line."""
    family = _family(family)
    _, mapping = matter.action_intertwiner()
    vertices = {tuple(map(str, row['origin'])): row for row in matter.invariant_higgs_lines()}
    rows = []
    for match in mapping:
        vertex = vertices[tuple(match['charge_origin'])]
        scale = s.sympify(match['same_current_tensor_scale_exact'])
        internal = tuple(matter.clean(scale*c) for c in vertex['components'])
        rows.append({'channel': match['existing_action_channel'], 'conjugate': vertex['conjugate'],
                     'internal': internal,
                     'components': tuple(s.ImmutableSparseMatrix(s.kronecker_product(c, family)) for c in internal),
                     'neutral_rank': vertex['neutral_rank']})
    return tuple(rows)


def gauge_attachment(family):
    family = _family(family)
    spectral_certificate(family)  # Positivity/rank claims below require this fixed source spectrum.
    r = matter.representation(); rows = higgs_kernels(family)
    checks = 0
    for row in rows:
        for current, lifted in zip(r['matter_currents'], r['lifted'], strict=True):
            h = -lifted[3:, 3:].T if row['conjugate'] else lifted[3:, 3:]
            for k in (0, 1):
                c = row['internal'][k]
                residual = current.T*c+c*current+sum((h[j, k]*row['internal'][j] for j in (0, 1)), s.zeros(15))
                if not _zero(residual):
                    raise ArithmeticError('an actual twelve-current Higgs tensor identity failed')
                # Exact Kronecker identity: currents are current15 tensor I3.
                if not _zero(s.kronecker_product(residual, family)):
                    raise ArithmeticError('the family lift changed a current identity')
                checks += 1
    internal = s.ImmutableSparseMatrix(sum((row['internal'][1] for row in rows), s.zeros(15)))
    square = matter.clean(internal.T*internal)
    if not square.is_diagonal() or sorted(square.diagonal()) != [0]+[1]*14:
        raise ArithmeticError('the common neutral internal matching changed')
    kernel = s.ImmutableSparseMatrix(s.kronecker_product(internal, family))
    if not _zero(kernel.T-kernel) or not _zero(kernel.T*kernel-s.kronecker_product(square, family*family)):
        raise ArithmeticError('the exact neutral mass-spectrum factorization failed')
    selected = list(r['selected'])
    electric = (r['operators']['T3']+r['operators']['q']/6).extract(selected, selected)
    charge = s.kronecker_product(electric, s.eye(3))
    if not _zero(charge.T*kernel+kernel*charge):
        raise ArithmeticError('the same-current electric charge does not preserve the pairing')
    pairing = s.kronecker_product(s.Matrix([[0, 1], [-1, 0]]), kernel)
    if not _zero(pairing.T+pairing):
        raise ArithmeticError('the separately chosen physical CAR pairing is not antisymmetric')
    serial = [[int(i), int(j), str(s.expand(x))] for (i, j), x in sorted(kernel.todok().items())]
    return {'current_count': 12, 'complete_Higgs_component_identities_exact': checks,
            'family_gauge_action': 'identity on a separately encoded neutral family factor',
            'physical_family_dimension_three_is_architecture': True,
            'all_three_channel_coefficients_are_the_identical_source_Y': True,
            'independent_species_strengths_or_empirical_yukawas_inserted': False,
            'neutral_Higgs_component_norm_in_model_units': '1',
            'physical_Higgs_scale_or_clock_calibration_derived': False,
            'canonical_family_kinetic_matrix': 'I3 on every internal and physical spin component',
            'same_internal_frame_preserves_all_three_Weyl_kinetic_axes': True,
            'neutral_K45_rank': 42, 'neutral_K45_zero_dimension': 3,
            'every_positive_family_mass_has_K45_singular_multiplicity': 14,
            'neutral_K45_squared_identity': 'diag(14 unit entries, one zero) tensor Y^2, in the actual internal ordering',
            'neutral_K45_sparse_sha256': hashlib.sha256(json.dumps(serial, separators=(',', ':')).encode()).hexdigest(),
            'physical_Delta90_antisymmetric_exact': True,
            'channels': [{'channel': row['channel'], 'same_Y': True,
                          'K45_singular_multiplicity_per_family': row['neutral_rank']} for row in rows],
            'equal_mass_degeneracy_across_up_down_and_charged_lepton_channels': True,
            'quark_charged_current_family_mixing': 'I3; identical left singular frames',
            'weak_CP_commutator_invariant_zero_exact': True,
            'three_minimal_left_Weyl_neutrinos_remain_massless': True,
            'quantum_gauge_Higgs_loops_or_confinement_computed': False,
            'laboratory_particle_poles_identified': False}


def whole_probability_shape_certificate(projectors):
    """A stationary-kernel exclusion independent of the event frequency.

    For either fivefold class the five distinct native port rays are
    equiangular. Cauchy--Binet and a nonnegative polynomial remainder bound
    determinant/e2^3 uniformly over all geometric-memory probabilities.
    """
    if len(projectors) != 5:
        raise ValueError('all five actual source orbit projectors are required')
    ps = tuple(_family(p) for p in projectors)
    for p in ps:
        if not _zero(p*p-p) or s.trace(p) != 1:
            raise ValueError('exact rank-one orthogonal source projectors required')
    for i, j in itertools.combinations(range(5), 2):
        if s.simplify(s.trace(ps[i]*ps[j])-s.Rational(1, 5)) != 0:
            raise ArithmeticError('the actual orbit is not the equiangular five-ray source')
    minimum = (10-2*s.sqrt(5))/25; maximum = (10+2*s.sqrt(5))/25
    triples = []
    for ids in itertools.combinations(range(5), 3):
        determinant = s.simplify(sum((ps[i] for i in ids), s.zeros(3)).det())
        if s.simplify(determinant-minimum) != 0 and s.simplify(determinant-maximum) != 0:
            raise ArithmeticError('a source triple Gram determinant left its exact two values')
        triples.append((ids, determinant))
    p = s.Symbol('p', real=True); S = sum(p**k for k in range(5))
    T = s.expand(sum(p**(sum(ids)-3) for ids, _ in triples))
    expected_remainder = p*(1-p)+p**3*(2-p+p**2)
    if s.expand(T-(1+p**2)**3-expected_remainder) != 0:
        raise ArithmeticError('the whole-probability nonnegative remainder identity failed')
    pairs = sum(p**(i+j) for i, j in itertools.combinations(range(5), 2))
    if s.expand(4*pairs-4*p*(1+p**2)*S) != 0:
        raise ArithmeticError('the exact second symmetric polynomial identity failed')
    lower = s.simplify(125*minimum/64)
    if lower != (25-5*s.sqrt(5))/32 or not (lower-s.Rational(1, 8)).is_positive:
        raise ArithmeticError('the uniform shape bound was not separated from the benchmark bound')
    return {'probability_domain': '0<p<=1; p=0 has rank one and two exact zero masses',
            'exclusion_subject': 'positive eigenvalue ratios of stationary averaged Yp with the declared canonical kinetic and Higgs attachment',
            'arbitrary_finite_step_CQ_channel_resonances_excluded': False,
            'p_equals_one_is_a_uniform_kernel_endpoint_not_a_reset_mixing_claim': True,
            'fivefold_classes_covered': ['72 degrees', '144 degrees'],
            'actual_source_pair_overlap_squared_exact': '1/5',
            'actual_source_triple_Gram_determinant_values_exact': [str(minimum), str(maximum)],
            'normalized_weights': 'p^k/S, k=0..4, S=1+p+p^2+p^3+p^4',
            'source_Y_trace_exact': '1',
            'scale_invariant_shape': 'det(Y)*tr(Y)^3/e2(Y)^3; source tr(Y)=1',
            'second_symmetric_polynomial': 'e2=4*p*(1+p^2)/(5*S)',
            'Cauchy_Binet_lower_bound': 'det(Y)>=((10-2*sqrt(5))/25)*p^3*T/S^3',
            'T_polynomial': str(T),
            'nonnegative_remainder_identity': 'T-(1+p^2)^3=p*(1-p)+p^3*(2-p+p^2)>=0',
            'determinant_over_e2_cubed_lower_exact': str(lower),
            'lower_bound_exceeds_one_eighth_exact': True,
            'source_inputs_or_probability_selected_by_empirical_comparison': False}


def heldout_comparison(spectrum, shape):
    """Read the retained benchmark only after all source calculations."""
    retained = strict_json.load(ROOT/BENCHMARK)['separate_empirical_comparison']
    windows = [tuple(F(x) for x in row) for row in retained['conservative_ratio_windows']]
    if windows != [(F(206), F(207)), (F(3470), F(3485))]:
        raise ArithmeticError('the retained broad held-out windows changed')
    bounds = [tuple(F(x) for x in row) for row in spectrum['ascending_ratios_after_one_exact']]
    passes = [max(a, c) <= min(b, d) for (a, b), (c, d) in zip(bounds, windows, strict=True)]
    (rlo, rhi), (klo, khi) = windows
    # For eigenvalue ratios (1,r,k), det/e2^3=r*k*(1+r+k)^3/(r+k+r*k)^3.
    upper = rhi*khi*(1+rhi+khi)**3/(rlo+klo+rlo*klo)**3
    if not upper < F(1, 8) or not shape['lower_bound_exceeds_one_eighth_exact']:
        raise ArithmeticError('the exact whole-law empirical exclusion failed')
    return {'comparison_only_source': retained['source'],
            'retained_broad_charged_lepton_ratio_windows': [list(map(str, row)) for row in windows],
            'windows_are_experimental_standard_uncertainties': False,
            'fixed_source_ratio_intersects_each_window': passes,
            'fixed_source_passes_both_windows': all(passes),
            'benchmark_determinant_over_e2_cubed_upper_exact': str(upper),
            'benchmark_shape_for_unnormalized_ratios': 'r*k*(1+r+k)^3/(r+k+r*k)^3',
            'benchmark_upper_less_than_one_eighth_exact': True,
            'all_probabilities_and_both_fivefold_classes_excluded_for_joint_windows': True,
            'other_interacting_repair_laws_excluded': False,
            'measured_values_used_to_generate_family_operator': False,
            'changing_common_Higgs_scale_or_clock_changes_ratios': False}


def report(family, *, projectors):
    spectrum = spectral_certificate(family)
    attachment = gauge_attachment(family)
    shape = whole_probability_shape_certificate(projectors)
    return {'source_spectrum': spectrum, 'same_current_Higgs_attachment': attachment,
            'native_family_band_theorem_connection': family_band_certificate(),
            'whole_probability_shape': shape,
            'separate_heldout_comparison': heldout_comparison(spectrum, shape)}


def source_bindings():
    names = set(matter.source_bindings()) | {BENCHMARK, FAMILY_THEOREM, PORT_ACTION_THEOREM, FAMILY_MANIFEST,
        'Lean/Screen/SMStructureAdequacySurface.lean',
        'model/common_source/particle_identification/family_attachment.py',
        'model/common_source/native_repair_mass.py', 'model/common_source/test_native_repair_mass.py',
        'model/common_source/native_repair_mass_attachment.py',
        'model/common_source/test_native_repair_mass_attachment.py'}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}
