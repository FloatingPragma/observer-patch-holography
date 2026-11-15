"""Exact native geometric loop phases and conditional C3 mass diagnostics.

Spin-coherent projectors on actual native port directions define a specified
geometric connection. Its phase is computed from their Bargmann products,
not supplied by a particle target. Koide balance and the identification
Y=C^2 with a canonical family kinetic form are separate supplied hypotheses.
No physical mass or complete native transport implementation follows merely
from this finite connection calculation.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math

import sympy as s

from tools import strict_json
from . import native_gauge_object as objects, slow_geometry

ROOT = objects.ROOT
BENCHMARK = 'model/common_source/chiral_candidate/results.json'
KOIDE_THEOREM = 'Lean/ObserverPatchHolography/KoideCirculant.lean'


def _zero(value):
    return all(s.expand(x) == 0 or s.simplify(x) == 0 for x in value)


def _scalar(value):
    return s.simplify(s.expand(value))


def _cyclic(face):
    at = face.index(min(face))
    return tuple(face[at:]+face[:at])


@lru_cache(maxsize=1)
def geometry():
    """Actual incidence, all twenty exact loops, and integral cycle coverage."""
    data = objects.geometry(); B = data['frame']
    squared = s.expand((B[0, :]*B[0, :].T)[0])
    directions = B/s.sqrt(squared)
    pauli = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    projectors = tuple((s.eye(2)+sum((directions[i, j]*pauli[j] for j in range(3)), s.zeros(2)))/2
                       for i in range(12))
    for projector in projectors:
        if not _zero(projector*projector-projector) or s.trace(projector) != 1:
            raise ArithmeticError('a native direction does not produce a spin-coherent projector')
    faces = []
    for source_face in slow_geometry.source().oriented_faces():
        face = list(source_face)
        volume = _scalar(directions[face, :].det())
        if volume.is_negative: face[1], face[2] = face[2], face[1]
        elif not volume.is_positive: raise ArithmeticError('native face orientation is not definite')
        faces.append(_cyclic(face))
    faces = tuple(sorted(faces)); edges = tuple(data['seams'])
    if len(faces) != len(set(faces)) or len(faces) != 20:
        raise ArithmeticError('all twenty native faces are required')
    q = s.sqrt(5)/5
    gram_volume_squared = (10+2*s.sqrt(5))/25
    real = (1+3*q)/4
    imag = s.sqrt(gram_volume_squared)/4
    norm_squared = ((1+q)/2)**3
    cosine = s.sqrt(10+2*s.sqrt(5))/4
    sine = (s.sqrt(5)-1)/4
    # Positive real/imaginary parts fix the root branches without numeric arg.
    if (not _zero(s.Matrix([real*real-norm_squared*cosine*cosine,
                            imag*imag-norm_squared*sine*sine]))
            or not (real-imag).is_positive or not sine.is_positive):
        raise ArithmeticError('the native Bargmann phase is not the stated first-quadrant root')
    phase = cosine+s.I*sine
    if _scalar(phase**5-s.I) != 0 or _scalar(phase**20-1) != 0:
        raise ArithmeticError('the exact native phase lattice changed')
    products = []
    for face in faces:
        n = [directions[i, :].T for i in face]
        dots = [_scalar((n[i].T*n[j])[0]) for i, j in ((0, 1), (1, 2), (2, 0))]
        volume = _scalar(s.Matrix.hstack(*n).det())
        if dots != [q, q, q] or _scalar(volume*volume-gram_volume_squared) != 0 or not volume.is_positive:
            raise ArithmeticError('a native face left the exact equilateral source geometry')
        actual = _scalar(s.trace(projectors[face[0]]*projectors[face[1]]*projectors[face[2]]))
        expected = (1+sum(dots)+s.I*volume)/4
        if _scalar(actual-expected) != 0 or _scalar(volume-4*imag) != 0:
            raise ArithmeticError('the actual projector product differs from the Pauli loop identity')
        products.append(actual)
    boundary1 = s.zeros(12, 30); boundary2 = s.zeros(30, 20)
    for j, (a, b) in enumerate(edges): boundary1[a, j], boundary1[b, j] = -1, 1
    for j, (a, b, c) in enumerate(faces):
        for x, y in ((a, b), (b, c), (c, a)):
            boundary2[edges.index(tuple(sorted((x, y)))), j] += 1 if x < y else -1
    if boundary1*boundary2 != s.zeros(12, 20) or boundary2*s.ones(20, 1) != s.zeros(30, 1):
        raise ArithmeticError('native face boundaries are not consistently oriented cycles')
    reached = {0}; tree = []
    while len(reached) != 12:
        crossing = next((j for j, (a, b) in enumerate(edges) if (a in reached) != (b in reached)), None)
        if crossing is None: raise ArithmeticError('native edge graph is not connected')
        tree.append(crossing); reached.update(edges[crossing])
    chords = tuple(j for j in range(30) if j not in tree)
    if boundary1.extract(range(12), tree).rank() != 11:
        raise ArithmeticError('the selected spanning tree has a nonzero integral cycle')
    reduced = boundary2.extract(chords, range(19))
    determinant = int(reduced.det())
    if abs(determinant) != 1:
        raise ArithmeticError('integral cycle coverage needs a unimodular reduced face matrix')
    inverse = reduced.inv()
    if any(x.is_Integer is not True for x in inverse):
        raise ArithmeticError('integer closed walks do not have integer face coefficients')
    # The selected face stabilizer is a regular C3 on the same native family band.
    face = faces[0]
    permutation = next(g for g in sorted(data['rotations'])
                       if all(g[face[i]] == face[(i+1) % 3] for i in range(3)))
    metric = s.expand((B.T*B)[0, 0])
    rotation = objects._clean(B.T*objects.permutation_matrix(permutation)*B/metric)
    if not _zero(rotation**3-s.eye(3)) or s.trace(rotation) != 0 or not _zero(rotation.T*rotation-s.eye(3)):
        raise ArithmeticError('the actual native family restriction is not regular C3')
    return {'directions': directions, 'projectors': projectors, 'faces': faces, 'edges': edges,
            'products': tuple(products), 'phase': phase, 'phase_norm_squared': norm_squared,
            'boundary1': boundary1, 'boundary2': boundary2, 'chords': chords, 'tree': tuple(tree),
            'reduced_determinant': determinant, 'reduced_inverse': inverse,
            'C3_permutation': permutation, 'C3_family_action': rotation}


def integral_face_decomposition(chain):
    """Every integral closed edge chain, without a bounded-walk enumeration."""
    data = geometry()
    if (type(chain) not in (list, tuple) or len(chain) != 30
            or any(type(x) is not int for x in chain)):
        raise ValueError('thirty original integer oriented-edge coefficients required')
    vector = s.Matrix(chain)
    if data['boundary1']*vector != s.zeros(12, 1):
        raise ValueError('the edge chain is not closed')
    coefficients = data['reduced_inverse']*vector.extract(data['chords'], [0])
    full = coefficients.col_join(s.zeros(1, 1))
    if data['boundary2']*full != vector or any(x.is_Integer is not True for x in full):
        raise ArithmeticError('integral native face decomposition failed')
    return tuple(int(x) for x in full)


def _sqrt_interval(bounds, digits=40):
    lo, hi = bounds
    if not 0 <= lo <= hi: raise ValueError('a nonnegative rational interval is required')
    scale = 10**digits
    lower = math.isqrt((lo.numerator*scale*scale)//lo.denominator)
    upper = math.isqrt((hi.numerator*scale*scale)//hi.denominator)+1
    result = F(lower, scale), F(upper, scale)
    if not result[0]**2 <= lo <= hi < result[1]**2:
        raise ArithmeticError('integer square-root enclosure failed')
    return result


def conditional_spectrum(k):
    """Certified roots for supplied balance and Y=C²; no measured inputs."""
    if type(k) is not int or k not in (0, 1, 2):
        raise ValueError('the three positive unordered phase classes are0,1,2')
    sqrt2 = _sqrt_interval((F(2), F(2))); sqrt5 = _sqrt_interval((F(5), F(5)))
    if k == 0: cosine = (F(1), F(1))
    elif k == 1: cosine = tuple(x/4 for x in _sqrt_interval((10+2*sqrt5[0], 10+2*sqrt5[1])))
    else: cosine = tuple((1+x)/4 for x in sqrt5)
    det = sqrt2[0]*cosine[0]/2-F(1, 2), sqrt2[1]*cosine[1]/2-F(1, 2)
    if k == 0:
        small = 1-sqrt2[1]/2, 1-sqrt2[0]/2
        roots = [small, small, (1+sqrt2[0], 1+sqrt2[1])]
    else:
        def sign(x):
            value = x**3-3*x*x+F(3, 2)*x
            if value-det[1] > 0: return 1
            if value-det[0] < 0: return -1
            raise ArithmeticError('root sign is unresolved at the declared rational precision')
        roots = []
        for left, right in ((F(0), F(1, 4)), (F(1, 3), F(1)), (F(2), F(3))):
            initial = sign(left)
            if initial == sign(right): raise ArithmeticError('conditional cubic root bracket failed')
            for _ in range(40):
                middle = (left+right)/2
                if sign(middle) == initial: left = middle
                else: right = middle
            roots.append((left, right))
    if roots[0][0] <= 0:
        raise ArithmeticError('the conditional square-root-mass carrier is not positive')
    ratios = [((row[0]/roots[0][1])**2, (row[1]/roots[0][0])**2) for row in roots[1:]]
    if k == 0: ratios[0] = F(1), F(1)
    return {'equal_link_phase': ('0' if k == 0 else 'pi/30' if k == 1 else 'pi/15'),
            'conditional_root_cubic': 'x^3-3*x^2+(3/2)*x-(-1/2+cos(3*delta)/sqrt(2))',
            'positive_root_intervals_exact': [list(map(str, row)) for row in roots],
            'conditional_mass_ratios_after_one_exact': [list(map(str, row)) for row in ratios],
            'conditional_mass_ratio_diagnostics': [[float(x) for x in row] for row in ratios],
            'signed_Koide_identity_is_physical_within_this_positive_mass_map': True}


def _heldout(rows):
    comparison = strict_json.load(ROOT/BENCHMARK)['separate_empirical_comparison']
    windows = [tuple(F(x) for x in row) for row in comparison['conservative_ratio_windows']]
    if windows != [(F(206), F(207)), (F(3470), F(3485))]:
        raise ArithmeticError('the retained broad lepton windows changed')
    result = []
    for row in rows:
        ratios = [tuple(F(x) for x in pair) for pair in row['conditional_mass_ratios_after_one_exact']]
        intersections = [max(a, c) <= min(b, d) for (a, b), (c, d) in zip(ratios, windows, strict=True)]
        result.append({'phase': row['equal_link_phase'], 'window_intersections': intersections,
                       'joint_pass': all(intersections)})
    return {'comparison_only_source': comparison['source'],
            'broad_lepton_ratio_windows': [list(map(str, row)) for row in windows],
            'windows_are_experimental_standard_uncertainties': False,
            'conditional_phase_classes': result,
            'complete_balanced_geometric_loop_class_excluded': not any(row['joint_pass'] for row in result),
            'continuously_chosen_or_refined_connection_phases_excluded': False,
            'phase_or_balance_adjusted_after_comparison': False}


def report():
    data = geometry()
    # Equal links differ by2pi/3 only by a permutation of the C3 eigenvalues.
    positive_residues = [k for k in range(20) if abs(k if k <= 10 else k-20)*12 <= 30]
    if positive_residues != [0, 1, 2, 18, 19]:
        raise ArithmeticError('positive balanced C3 chamber classification changed')
    rows = [conditional_spectrum(k) for k in (0, 1, 2)]
    return {'source_connection': {
                'native_port_directions': 12, 'actual_oriented_faces_checked': 20,
                'faces': [list(face) for face in data['faces']],
                'spin_half_projector': '(I+n_i dot sigma)/2',
                'triangle_Bargmann_identity': 'Tr(F_i F_j F_k)=(1+n_i.n_j+n_j.n_k+n_k.n_i+i*n_i.(n_j cross n_k))/4',
                'native_adjacent_dot_exact': '1/sqrt(5)',
                'native_oriented_triple_squared_exact': '(10+2*sqrt(5))/25',
                'spin_half_oriented_face_phase': 'pi/10',
                'spin_one_rule': 'symmetric square of the spin-half coherent line; Bargmann phase doubles',
                'spin_one_oriented_face_phase': 'pi/5',
                'normalized_face_phase_exact': str(data['phase']),
                'phase_fifth_power_exact': 'I', 'phase_order': 20,
                'phase_argument_strictly_between_zero_and_pi_over4': True,
                'source_phase_depends_on_measured_masses': False},
            'integral_closed_walk_certificate': {
                'native_edges': 30, 'connected_spanning_tree_edges': list(data['tree']),
                'non_tree_edge_indices': list(data['chords']),
                'face_basis_omits_last_face': True,
                'reduced_face_boundary_determinant_exact': data['reduced_determinant'],
                'integer_cycle_coefficients_from_integer_inverse': True,
                'all_integer_closed_edge_chains_covered': True,
                'face_decomposition_ambiguity': 'adding every oriented face changes total coefficient by20 and phase by2pi',
                'spin_half_closed_walk_phase_lattice': 'integer multiples of pi/10 modulo2pi',
                'spin_one_closed_walk_phase_lattice': 'integer multiples of pi/5 modulo2pi'},
            'actual_family_C3_attachment': {
                'native_face_rotation': list(data['C3_permutation']),
                'order_three_orthogonal_family_action_exact': [[str(x) for x in row] for row in data['C3_family_action'].tolist()],
                'trace': '0', 'restriction_is_regular_C3': True,
                'same_current_Higgs_family_substitution': 'Y=C^2 in the already verified generic complex-family identity',
                'physical_family_copy_is_neutral_under_internal_gauge_currents': True,
                'full_A5_geometric3_plus3prime_equals_family3_tensor_trivial_orientation2': False},
            'conditional_mass_contract': {
                'circulant': 'C=a*I+b*R+conjugate(b)*R^2, a>0',
                'supplied_balance': '|b/a|=1/sqrt(2)',
                'actual_native_GNS_orientation_attachment_executed_here': False,
                'supplied_mass_map': 'Y=C^2 with canonical family kinetics; a common Higgs/clock scale is separate',
                'equal_link_connection_gauge': 'delta=loop_phase/3 modulo2pi/3',
                'positive_chamber': '|delta|<=pi/12 modulo2pi/3',
                'positive_chamber_reason': 'For 0<=delta<=pi/3 the smallest root is 1+sqrt(2)*cos(delta+2pi/3); this root decreases to zero at delta=pi/12. Orientation reversal preserves the unordered spectrum.',
                'positive_spin_half_lattice_residues_mod20': positive_residues,
                'unordered_positive_phase_classes': ['0', 'pi/30', 'pi/15'],
                'fixed_single_face_rules': {'spin_half': 'pi/30', 'spin_one': 'pi/15'},
                'positivity_and_balance_select_a_unique_continuous_phase': False,
                'source_selected_physical_square_root_mass_map': False,
                'interacting_poles_or_physical_clock_derived': False,
                'complete_native_connection_transport_words_executed': False},
            'conditional_spectra': rows, 'separate_heldout_comparison': _heldout(rows),
            'finite_theorem': KOIDE_THEOREM,
            'primary_geometric_phase_sources': [
                'https://journals.aps.org/pra/abstract/10.1103/PhysRevA.41.42',
                'https://arxiv.org/abs/quant-ph/9907080']}


def source_bindings():
    from . import native_matter_attachment as matter
    names = set(objects.source_bindings()) | set(matter.source_bindings()) | {BENCHMARK, KOIDE_THEOREM,
        'model/common_source/native_mass_holonomy.py',
        'model/common_source/test_native_mass_holonomy.py',
        'model/common_source/native_matter_attachment.py',
        'model/common_source/test_native_matter_attachment.py'}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}
