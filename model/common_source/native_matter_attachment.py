"""A chosen spin/exterior repair register for the actual native current.

The twelve-port carrier, incidence and geometric rotations are unchanged.
Auxiliary encoded registers implement a declared spin lift and nonvacuum
exterior sector. Their currents, charges and invariant Higgs lines are
computed from native_sm_current, and their complete source records are read
without a species-name table. This constructs an internal matter interface;
physical Weyl statistics, propagation and a mass law have separate contracts.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from fractions import Fraction
import hashlib
import itertools
import json
import math
from numbers import Real
from pathlib import Path

import numpy as np
from scipy.linalg import expm
import sympy as s

from tools import strict_json
from . import native_sm_current as native_current
from .chiral_candidate import action as existing_action, gaussian_source as gaussian
from ..observer.source_scalar import compiler

ROOT = Path(__file__).resolve().parents[2]
ADMISSION = 2e-8
TRIALS = 262144
SEED = 20261027
SAMPLER_BOUNDARY_GUARD = 2**-32


def clean(matrix):
    return s.ImmutableSparseMatrix(matrix.applyfunc(lambda x: s.expand(s.radsimp(x))))


def zero(matrix):
    return all(s.expand(x) == 0 for x in matrix)


def entries(matrix):
    return [[str(x) for x in matrix.row(i)] for i in range(matrix.rows)]


def unit(n, i, j):
    result = s.zeros(n); result[i, j] = 1
    return result


def spin_lift(generator):
    """Declared Lie lift; the spinor is an auxiliary encoded representation."""
    if (not isinstance(generator, s.MatrixBase) or generator.shape != (6, 6)
            or any(x.has(s.Float) or x.free_symbols or x.is_finite is not True for x in generator)
            or not zero(generator+generator.H) or not zero(generator[:3, 3:])
            or not zero(generator[3:, :3]) or not zero(s.im(generator[3:, 3:]))):
        raise ValueError('the declared u3 plus real-so3 current domain is required')
    a, b = generator[:3, :3], generator[3:, 3:]
    w = (b[2, 1], b[0, 2], b[1, 0])
    pauli = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    weak = sum((-s.I*v*p/2 for v, p in zip(w, pauli, strict=True)), s.zeros(2))
    return clean(s.diag(a, weak-s.trace(a)*s.eye(2)/2))


@lru_cache(maxsize=1)
def fock_data():
    """All32 occupation masks, with an independently checked CAR convention."""
    creation = []
    for site in range(5):
        matrix = s.zeros(32)
        for mask in range(32):
            if not mask >> site & 1:
                matrix[mask | 1 << site, mask] = (-1)**((mask & ((1 << site)-1)).bit_count())
        creation.append(s.ImmutableSparseMatrix(matrix))
    for i, a in enumerate(creation):
        for j, b in enumerate(creation):
            if a*b+b*a != s.zeros(32) or a.H*b+b*a.H != int(i == j)*s.eye(32):
                raise ArithmeticError('auxiliary exterior CAR convention failed')
    return tuple(creation)


def exterior_derivative(generator):
    if (not isinstance(generator, s.MatrixBase) or generator.shape != (5, 5)
            or any(x.has(s.Float) or x.free_symbols or x.is_finite is not True for x in generator)):
        raise ValueError('an exact five-mode generator is required')
    c = fock_data()
    return clean(sum((generator[i, j]*c[i]*c[j].H for i in range(5) for j in range(5)
                      if generator[i, j] != 0), s.zeros(32)))


@lru_cache(maxsize=1)
def canonical_currents():
    """A basis requested through the actual native twelve-port control map."""
    result = []
    for i, j in itertools.combinations(range(3), 2):
        result.extend((unit(6, i, j)-unit(6, j, i), s.I*(unit(6, i, j)+unit(6, j, i))))
    result.extend((s.I*(unit(6, 0, 0)-unit(6, 1, 1)), s.I*(unit(6, 1, 1)-unit(6, 2, 2))))
    for vector in s.eye(3).columnspace():
        value = s.zeros(6); value[3:, 3:] = native_current._hat(vector)
        result.append(value)
    result.append(native_current.exact_architecture()['central'])
    for value in result:
        coefficients = native_current.admit_generator(value)
        if not zero(sum((a*b for a, b in zip(coefficients, native_current.exact_architecture()['generators'], strict=True)), s.zeros(6))-value):
            raise ArithmeticError('a requested current is not supplied by the native ports')
    return tuple(s.ImmutableMatrix(x) for x in result)


@lru_cache(maxsize=1)
def representation():
    native = native_current.exact_architecture()
    lifted = tuple(spin_lift(d) for d in native['generators'])
    currents = tuple(exterior_derivative(d) for d in lifted)
    selected = tuple(mask for mask in range(32) if mask.bit_count() in (2, 4))
    encoder = s.eye(32)[:, list(selected)]
    native_mutable = tuple(s.MutableSparseMatrix(g) for g in native['generators'])
    lifted_mutable = tuple(s.MutableSparseMatrix(g) for g in lifted)
    current_mutable = tuple(s.MutableSparseMatrix(g) for g in currents)
    for a, b in itertools.combinations(range(12), 2):
        bracket = clean(native_mutable[a]*native_mutable[b]-native_mutable[b]*native_mutable[a])
        lifted_bracket = clean(lifted_mutable[a]*lifted_mutable[b]-lifted_mutable[b]*lifted_mutable[a])
        if not zero(spin_lift(bracket)-lifted_bracket):
            raise ArithmeticError('same-current spin lift lost a bracket')
        derived = exterior_derivative(lifted_bracket)
        current_bracket = clean(current_mutable[a]*current_mutable[b]-current_mutable[b]*current_mutable[a])
        if not zero(derived-current_bracket):
            raise ArithmeticError('same-current exterior lift lost a bracket')
    center = clean(spin_lift(native['central'])/s.I)
    scale = s.ilcm(*(x.q for x in center.diagonal()))
    integers = [int(scale*x) for x in center.diagonal()]
    charge5 = clean(-scale*center/s.igcd(*integers))
    q = exterior_derivative(charge5)
    e = {(i, j): exterior_derivative(unit(5, i, j)) for i in range(3) for j in range(3)}
    nc = sum((e[i, i] for i in range(3)), s.zeros(32))
    c3 = clean(sum((e[i, j]*e[j, i]/2 for i in range(3) for j in range(3)), s.zeros(32))-nc*nc/6)
    pauli = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    weak = []
    for p in pauli:
        generator = s.zeros(5); generator[3:, 3:] = p/2
        weak.append(exterior_derivative(generator))
    c2 = clean(sum((x*x for x in weak), s.zeros(32)))
    operators = {'q': q, 'color_C2': c3, 'weak_C2': c2, 'T3': weak[2],
                 'H1': e[0, 0]-e[1, 1], 'H2': e[1, 1]-e[2, 2]}
    if any(value != s.diag(*value.diagonal()) for value in operators.values()):
        raise ArithmeticError('joint occupation record does not diagonalize the declared readout')
    matter = tuple(clean(encoder.T*g*encoder) for g in currents)
    if any(not zero(g*encoder-encoder*small) for g, small in zip(currents, matter, strict=True)):
        raise ArithmeticError('the selected nonvacuum even register is not invariant')
    return {'native': native, 'lifted': lifted, 'full_currents': currents, 'selected': selected,
            'encoder': encoder, 'matter_currents': matter, 'operators': operators,
            'primitive_charge5': charge5, 'center5': center}


def decoded_multiplets():
    """Group only source-computed joint charge/Casimir eigenrecords."""
    r = representation(); ops = r['operators']; groups = {}
    for mask in r['selected']:
        key = tuple(ops[name][mask, mask] for name in ('q', 'color_C2', 'weak_C2'))
        groups.setdefault(key, []).append(mask)
    rows = []
    for (charge, color, weak), masks in sorted(groups.items()):
        weights = {}
        for mask in masks:
            weight = (ops['H1'][mask, mask], ops['H2'][mask, mask], ops['T3'][mask, mask])
            weights[weight] = weights.get(weight, 0)+1
        rows.append({'q6Y': str(charge), 'Y': str(charge/6), 'color_C2': str(color),
                     'weak_C2': str(weak), 'dimension': len(masks), 'occupation_records': list(masks),
                     'weights': [{'H1': str(a), 'H2': str(b), 'T3': str(t),
                                  'electric_Q': str(t+charge/6), 'multiplicity': n}
                                 for (a, b, t), n in sorted(weights.items())]})
    return rows


@lru_cache(maxsize=1)
def invariant_higgs_lines():
    """Solve each allowed charge-origin tensor equation from the new currents."""
    r = representation(); indices = list(r['selected'])
    q = r['operators']['q'].extract(indices, indices)
    currents = r['matter_currents']; origins = sorted({tuple(sorted((q[i, i], q[j, j])))
        for i in range(15) for j in range(i, 15) if abs(q[i, i]+q[j, j]) == 3})
    result = []
    for origin in origins:
        pairs = [(i, j) for i in range(15) for j in range(i, 15)
                 if tuple(sorted((q[i, i], q[j, j]))) == origin]
        basis = []
        for i, j in pairs:
            c = unit(15, i, j)
            if i != j: c += unit(15, j, i)
            basis.append(c)
        conjugate = sum(origin) > 0
        higgs = tuple(-g[3:, 3:].T if conjugate else g[3:, 3:] for g in r['lifted'])
        equations = []
        for g, h in zip(currents, higgs, strict=True):
            action = s.Matrix.hstack(*(s.Matrix([(g.T*c+c*g)[i, j] for i, j in pairs]) for c in basis))
            equations.append(s.kronecker_product(s.eye(2), action)+s.kronecker_product(h.T, s.eye(len(pairs))))
        equation = clean(s.Matrix.vstack(*equations))
        kernel = equation.to_DM(extension=True).to_field().nullspace(divide_last=True).to_Matrix()
        if kernel.rows != 1:
            raise ArithmeticError('the source-current Higgs origin does not have one invariant line')
        coefficients = kernel.row(0)
        cs = [clean(sum((coefficients[k*len(pairs)+j]*c for j, c in enumerate(basis)), s.zeros(15))) for k in range(2)]
        normalization = next(x for x in cs[1] if x != 0)
        cs = tuple(clean(c/normalization) for c in cs)
        for g, h in zip(currents, higgs, strict=True):
            for k in range(2):
                if not zero(g.T*cs[k]+cs[k]*g+sum((h[j, k]*cs[j] for j in range(2)), s.zeros(15))):
                    raise ArithmeticError('the derived Higgs tensor failed complete current substitution')
        result.append({'origin': origin, 'conjugate': conjugate, 'components': cs,
                       'neutral_rank': cs[1].to_DM(extension=True).rank()})
    return tuple(result)


@lru_cache(maxsize=1)
def action_intertwiner():
    """An explicit same-current map into the existing kinetic/Higgs action."""
    r = representation(); old = existing_action._representation()
    new_subsets = tuple(tuple(i for i in range(5) if mask >> i & 1) for mask in r['selected'])
    p = s.zeros(15)
    for old_index, subset in enumerate(old[0]): p[old_index, new_subsets.index(subset)] = 1
    for port, current in enumerate(r['matter_currents']):
        if not zero(old[2][r['native']['native_to_fixture_port'][port]]*p-p*current):
            raise ArithmeticError('native matter and existing field currents have no stated intertwiner')
    rows = []
    for vertex in invariant_higgs_lines():
        matches = [(name, data) for name, data in old[4].items() if tuple(data['origin']) == vertex['origin']]
        if len(matches) != 1: raise ArithmeticError('derived Higgs origin has no unique existing field counterpart')
        name, data = matches[0]
        transformed = [clean(p.T*c*p) for c in data['components']]
        a, b = next((a, b) for a, b in zip(transformed[1], vertex['components'][1], strict=True) if b != 0)
        scale = s.simplify(a/b)
        if any(not zero(a-scale*b) for a, b in zip(transformed, vertex['components'], strict=True)):
            raise ArithmeticError('existing and freshly derived Higgs lines disagree')
        rows.append({'charge_origin': list(map(str, vertex['origin'])), 'existing_action_channel': name,
                     'same_current_tensor_scale_exact': str(scale), 'neutral_rank': vertex['neutral_rank']})
    q15 = r['operators']['q'].extract(list(r['selected']), list(r['selected']))
    em = r['operators']['T3'].extract(list(r['selected']), list(r['selected']))+q15/6
    if not zero(existing_action.electric_charge(families=False)*p-p*em):
        raise ArithmeticError('same-current electric charge dictionary disagrees')
    return p, rows


def field_intertwiner():
    """Separate physical CAR ordering: spin, internal matter mask, family."""
    p, _ = action_intertwiner()
    return s.ImmutableSparseMatrix(s.kronecker_product(s.eye(2), p, s.eye(3)))


def kinetic_attachment():
    """Generic complex-family action handoff, independent of any mass packet.

    Each of three physical family matrices has nine independent complex
    symbols. The weak-doublet to singlet orientation is recovered from the
    source weak Casimir, and the opposite block is transposed, not adjointed.
    Creation pairing conjugates this all-left action kernel. Thus the same
    map applies to off-diagonal CKM-compatible matrices and every calibrated
    or uncalibrated choice without executing a parameter producer.
    """
    p, rows = action_intertwiner()
    family = s.ImmutableSparseMatrix(s.kronecker_product(p, s.eye(3)))
    old_rep = existing_action._representation()
    native = representation(); selected = list(native['selected'])
    new_weak = native['operators']['weak_C2'].extract(selected, selected)
    old_weak = old_rep[3]['weak_C2']
    symbols = s.symbols('y0:27')
    vertices = {tuple(map(str, row['origin'])): row for row in invariant_higgs_lines()}

    def oriented_kernel(component, weak, y):
        directed = s.MutableSparseMatrix.zeros(45, 45)
        for (i, j), value in component.todok().items():
            if weak[i, i] == s.Rational(3, 4) and weak[j, j] == 0:
                for a, b in itertools.product(range(3), repeat=2):
                    directed[3*i+a, 3*j+b] = value*s.conjugate(y[a, b])
        return clean(directed+directed.T)

    old_kernels = [s.MutableSparseMatrix.zeros(45, 45) for _ in range(2)]
    new_kernels = [s.MutableSparseMatrix.zeros(45, 45) for _ in range(2)]
    for channel, row in enumerate(rows):
        vertex = vertices[tuple(row['charge_origin'])]
        y = s.Matrix(3, 3, symbols[9*channel:9*(channel+1)])
        existing = old_rep[4][row['existing_action_channel']]['components']
        for component in (0, 1):
            old = oriented_kernel(existing[component], old_weak, y)
            new = oriented_kernel(s.sympify(row['same_current_tensor_scale_exact'])*vertex['components'][component], new_weak, y)
            if not zero(family.T*old*family-new):
                raise ArithmeticError('generic complex family action lost its same-current intertwiner')
            old_kernels[component] += old
            new_kernels[component] += new
    physical = field_intertwiner()
    if physical.T*physical != s.eye(90):
        raise ArithmeticError('the separate physical CAR frame is not canonical')
    pauli = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    for sigma in pauli:
        kinetic = s.kronecker_product(sigma, s.eye(45))
        if kinetic*physical != physical*kinetic:
            raise ArithmeticError('the matter frame changed the supplied Weyl kinetic symbol')
    epsilon = s.Matrix([[0, 1], [-1, 0]])
    for old, new in zip(old_kernels, new_kernels, strict=True):
        old_pair = s.kronecker_product(epsilon, old.conjugate())
        new_pair = s.kronecker_product(epsilon, new.conjugate())
        if not zero(new_pair.T+new_pair) or not zero(physical.T*old_pair*physical-new_pair):
            raise ArithmeticError('the generic physical CAR pairing and matter frame disagree')
    return {'same_existing_uncalibrated_action_kernel_exact': True,
            'arbitrary_complex_family_action_intertwines_exactly': True,
            'independent_complex_family_symbols': len(symbols),
            'both_Higgs_components_and_all_offdiagonal_family_entries_checked': True,
            'all_left_directed_family_coefficient': 'Y.conjugate(); reverse block is transpose',
            'creation_pairing_coefficient': 'epsilon_spin tensor K.conjugate()',
            'Nambu_intertwiner': 'diag(P90,P90.conjugate())',
            'all_three_Weyl_principal_symbol_axes_intertwine_exactly': True,
            'physical_CAR_frame_is_unitary_and_pairing_antisymmetric': True,
            'native_to_existing_tensor_signs_retained': True,
            'existing_action_hierarchy_failure_removed_by_frame_change': False,
            'mass_parameters_refitted_in_this_handoff': False,
            'calibrated_or_uncalibrated_parameter_producer_executed_for_this_identity': False}


def relative_spectrum_freedom():
    """Exact spectrum countermodels inside the same native twelve-current law.

    Gauge currents act trivially on the chosen three-copy family factor.
    Therefore each of the three verified one-dimensional Higgs invariant
    lines is tensored with all of Mat3(C). With canonical common kinetic
    normalization, neutral-Higgs mass singular values are v*sigma(Y).
    These are action parameters; no loop-corrected laboratory pole is implied.
    """
    p, mapping = action_intertwiner()
    vertices = {tuple(map(str, row['origin'])): row for row in invariant_higgs_lines()}
    controls = (
        {'up': (1, 2, 3), 'down': (2, 3, 4), 'charged_lepton': (3, 4, 5)},
        {'up': (1, 4, 9), 'down': (2, 6, 12), 'charged_lepton': (3, 8, 15)})
    selected = list(representation()['selected'])
    charge = representation()['operators']['T3'].extract(selected, selected)
    charge += representation()['operators']['q'].extract(selected, selected)/6
    charge45 = s.kronecker_product(charge, s.eye(3))
    rows = []
    for index, parameters in enumerate(controls):
        kernel = s.MutableSparseMatrix.zeros(45, 45)
        expected = [s.Integer(0)]*3
        channels = []
        for row in mapping:
            channel = row['existing_action_channel']
            values = tuple(map(s.Integer, parameters[channel]))
            if not 0 < values[0] < values[1] < values[2]:
                raise ArithmeticError('the fixed comparison hierarchy is not positive and ordered')
            vertex = vertices[tuple(row['charge_origin'])]
            internal = s.sympify(row['same_current_tensor_scale_exact'])*vertex['components'][1]
            kernel += s.kronecker_product(internal, s.diag(*values))
            expected += [value**2 for value in values for _ in range(vertex['neutral_rank'])]
            channels.append({'channel': channel, 'positive_Yukawa_singular_values_exact': list(map(str, values)),
                             'relative_to_lightest_exact': [str(value/values[0]) for value in values],
                             'K45_singular_multiplicity_per_family': vertex['neutral_rank']})
        square = clean(kernel.H*kernel)
        if square != s.diag(*square.diagonal()) or sorted(square.diagonal()) != sorted(expected):
            raise ArithmeticError('complete neutral action spectrum differs from its Higgs singular values')
        if not zero(charge45.T*kernel+kernel*charge45):
            raise ArithmeticError('a hierarchy comparison violates the same electric-current action')
        rows.append({'name': 'A' if index == 0 else 'B', 'neutral_Higgs_norm': '1',
                     'channels': channels, 'complete_K45_mass_squared_multiset_exact': list(map(str, sorted(expected))),
                     'neutral_massless_family_components': 3,
                     'same_verified_twelve_current_Higgs_tensor_factorization': True,
                     'canonical_common_kinetic_form': 'I45'})
    first = {row['channel']: row['relative_to_lightest_exact'] for row in rows[0]['channels']}
    second = {row['channel']: row['relative_to_lightest_exact'] for row in rows[1]['channels']}
    if first == second:
        raise ArithmeticError('the exact relative-spectrum countermodels accidentally agree')
    return {'Higgs_invariant_line_count': len(mapping), 'family_factor_dimension': 3,
            'family_gauge_action': 'identity on the declared three-copy factor',
            'complex_Yukawa_coefficient_dimension_before_family_basis_quotient': 27,
            'dimension_reason': 'three invariant Higgs lines tensor Mat3(C)',
            'charged_rest_action_mass_formula': 'm_a,f=v*sigma_f(Y_a), a=up,down,charged_lepton',
            'common_positive_Higgs_norm_cancels_all_mass_ratios': True,
            'physical_mass_scale_or_ratios_follow_from_twelve_currents_alone': False,
            'same_twelve_currents_fix_additional_Yukawa_entries': False,
            'full_repair_dynamics_or_vacuum_could_supply_additional_constraints': True,
            'full_mass_actions_of_controls_are_identical': False,
            'exact_positive_hierarchy_controls': rows,
            'laboratory_mass_or_coupling_inputs': [],
            'loop_or_confinement_corrected_particle_poles_computed': False}


def global_cover_contract():
    """The chosen lift has the usual compatible central quotient.

    On the five-mode space, a central color/weak/hypercharge triple acts as
    diag(omega3^k exp(-2i alpha) I3, (-1)^ell exp(3i alpha) I2).
    Its identity conditions give alpha=n*pi/3, k=n mod3, ell=n mod2.
    The nonvacuum even module is faithful on SU5: trivial wedge2 action
    forces all five eigenvalues equal to +/-1, and determinant1 excludes-1.
    """
    q5 = representation()['primitive_charge5']
    if s.trace(q5) != 0:
        raise ArithmeticError('the chosen lift does not have determinant-one central action')
    kernels = []
    for k, ell, n in itertools.product(range(3), range(2), range(6)):
        phases = [s.exp(2*s.pi*s.I*s.Rational(k, 3))]*3+[s.Integer(-1)**ell]*2
        one = s.diag(*(s.simplify(s.expand_complex(value*s.exp(s.pi*s.I*q5[j, j]*n/3)))
                       for j, value in enumerate(phases)))
        selected = representation()['selected']
        actions = [s.simplify(s.prod(one[j, j] for j in range(5) if mask >> j & 1)) for mask in selected]
        if all(value == 1 for value in actions): kernels.append([k, ell, n])
    if kernels != sorted([[n % 3, n % 2, n] for n in range(6)]):
        raise ArithmeticError('the primitive charge and spin lift have incompatible global centers')
    return {'primitive_charge_period': '2*pi',
            'compatible_group_quotient': '(SU3 x SU2 x U1_q)/Z6',
            'central_kernel_triples_color_weak_qsixthroot': kernels,
            'weak_spin_lift_descends_to_original_SO3_without_a_cover': False,
            'global_spin_cover_and_central_period_are_architecture': True,
            'auxiliary_even_parity_identifies_physical_fermion_statistics': False}


@lru_cache(maxsize=1)
def primitive_coefficients():
    """Reconstruct executable controls from actual currents, not fixture rows."""
    targets = {'X': s.I*(unit(6, 0, 1)+unit(6, 1, 0)),
               'Y': unit(6, 0, 1)-unit(6, 1, 0),
               'Z': s.I*(unit(6, 0, 0)-unit(6, 1, 1)),
               'I': s.I*(unit(6, 0, 0)+unit(6, 1, 1)),
               'J': representation()['native']['central']}
    targets['H'] = (targets['X']+targets['Z'])/s.sqrt(2)
    return {kind: native_current.admit_generator(target) for kind, target in targets.items()}


@lru_cache(maxsize=4096)
def native_pulse(kind, angle):
    if (kind not in compiler.KINDS or isinstance(angle, (bool, np.bool_))
            or not isinstance(angle, Real) or not math.isfinite(float(angle))):
        raise ValueError('a finite declared native pulse is required')
    output = expm(float(angle)*primitive_generator(kind))
    output.flags.writeable = False
    return output


@lru_cache(maxsize=6)
def primitive_generator(kind):
    """One source-span solve per kind, independent of programmed pulse angle."""
    if kind not in compiler.KINDS: raise ValueError('a declared native primitive kind is required')
    coefficients = primitive_coefficients()[kind]
    exact = clean(sum((a*b for a, b in zip(coefficients, representation()['native']['generators'], strict=True)), s.zeros(6)))
    result = np.asarray(exact, complex)
    result.flags.writeable = False
    return result


def execute_word(word, qubits, state=None):
    """Full logical source columns; proper-code complement effects are checked.

    This executes the declared ideal continuous-current service. The separate
    native_sm_current finite64-row construction supplies implementation
    evidence for this architecture at its own tested requests. It does not
    certify every angle or composed primitive in this matter word and is not
    silently substituted for the ideal service.
    """
    if type(qubits) is not int or not 1 <= qubits <= 10:
        raise ValueError('one through ten declared encoded qubits are required')
    word = compiler._validate(word, qubits)
    out = np.eye(1 << qubits, dtype=complex) if state is None else np.array(state, complex, copy=True)
    vector = out.ndim == 1
    if out.shape[0] != 1 << qubits or not np.isfinite(out).all():
        raise ValueError('complete finite register/reference columns required')
    out = out.reshape(1 << qubits, -1)
    for event in word:
        full = native_pulse(event['kind'], event['angle'])
        d = 1 << len(event['targets'])
        if np.linalg.norm(full[d:, :d], 2) > 2e-12:
            raise ArithmeticError('the declared native proper-code primitive leaks')
        out = compiler._apply(out, full[:d, :d], event['targets'], qubits)
    return out[:, 0] if vector else out


def resources(word):
    """Schedule counts and proper-code check from this actual native backend."""
    word = compiler._validate(word)
    pair = sum(event['kind'] == 'J' for event in word)
    unique = {(event['kind'], event['angle']) for event in word}
    leakage = max((np.linalg.norm(native_pulse(kind, angle)[(4 if kind == 'J' else 2):,
        :(4 if kind == 'J' else 2)], 2) for kind, angle in unique), default=0.)
    return {'native_pulses': len(word), 'single_code_pulses': len(word)-pair, 'pair_code_pulses': pair,
            'zero_angle_requests': sum(event['angle'] == 0 for event in word),
            'identity_phase_pulses': sum(event['kind'] == 'I' for event in word),
            'proper_code_load_events': len(word), 'proper_code_export_events': len(word),
            'proper_code_service_events': 2*len(word),
            'integrated_processor_norm_upper': math.fsum(abs(event['angle']) for event in word),
            'maximum_work_qubits': max((max(event['targets']) for event in word), default=-1)+1,
            'maximum_numerical_valid_code_leakage': float(leakage),
            'service_failure_outcomes_retained': True, 'valid_code_failure_probability_in_exact_arithmetic': 0,
            'initialization_storage_readout_physical_costs_included': False}


def onebody_word(generator, duration):
    """Exact commuting-term source synthesis for the selected canonical pulses."""
    h = clean(s.I*generator)
    if h.shape != (5, 5) or not zero(h-h.H): raise ValueError('an anti-Hermitian five-mode current is required')
    pieces, word = [], []
    for i in range(5):
        if h[i, i] != 0:
            pieces.append(h[i, i]*unit(5, i, i))
            word += compiler.remap(gaussian.local_word('onsite', complex(h[i, i]), duration), [i])
        for j in range(i+1, 5):
            if h[i, j] == 0: continue
            pieces.append(h[i, j]*unit(5, i, j)+h[j, i]*unit(5, j, i))
            route = [(k-1, k) for k in range(j, i+1, -1)]
            for pair in route: word += compiler.remap(gaussian.local_word('fswap'), pair)
            word += compiler.remap(gaussian.local_word('hopping', complex(h[i, j]), duration), [i, i+1])
            for pair in reversed(route): word += compiler.remap(gaussian.local_word('fswap'), pair)
    if any(not zero(a*b-b*a) for a, b in itertools.combinations(pieces, 2)):
        raise ValueError('this exact pulse compiler requires commuting one-body terms')
    return tuple(compiler._validate(word, 5))


def current_execution():
    arrays, rows = {}, []
    requests = [(f'current_{i}', spin_lift(g), 0.125) for i, g in enumerate(canonical_currents())]
    requests += [('weak_2pi', spin_lift(canonical_currents()[10]), 2*math.pi),
                 ('weak_4pi', spin_lift(canonical_currents()[10]), 4*math.pi)]
    for name, generator, duration in requests:
        word = onebody_word(generator, duration)
        output = execute_word(word, 5)
        target = expm(duration*np.asarray(exterior_derivative(generator), complex))
        error = float(np.linalg.norm(output-target, 2))
        if not math.isfinite(error) or error > ADMISSION: raise ArithmeticError('complete source current action failed')
        arrays[name+'_source'] = output
        rows.append({'name': name, 'response_parameter': duration, 'sourceword': list(word),
                     'resources': resources(word), 'input_columns': 32, 'output_rows': 32,
                     'numerical_operator_admission': ADMISSION})
    weak_sign = np.diag([(-1)**((mask >> 3).bit_count()) for mask in range(32)])
    if np.linalg.norm(arrays['weak_2pi_source']-weak_sign, 2) > ADMISSION or np.linalg.norm(arrays['weak_4pi_source']-np.eye(32), 2) > ADMISSION:
        raise ArithmeticError('actual source spin-cover holonomy failed')
    return rows, arrays


def occupation_instrument():
    """All32 pointer outcomes, all32 matter rows and arbitrary input columns."""
    pointer_word = compiler.cnot(0, 1)
    pointer = execute_word(pointer_word, 2)
    local = (pointer[:2, :2], pointer[2:, :2])
    branches = {0: np.eye(32, dtype=complex)}
    for site in range(5):
        branches = {record | bit << site: compiler._apply(matrix, kraus, [site], 5)
                    for record, matrix in branches.items() for bit, kraus in enumerate(local)}
    kraus = np.stack([branches[record] for record in range(32)])
    if np.linalg.norm(sum((k.conj().T@k for k in kraus), np.zeros((32, 32), complex))-np.eye(32), 2) > ADMISSION:
        raise ArithmeticError('complete source occupation instrument failed')
    return kraus, pointer_word


def proper_code_exports():
    """All complement outcomes of the declared ideal native load/export law."""
    arrays = {}
    for dimension in (2, 4):
        kraus = np.asarray(compiler.transfer_kraus(dimension))
        if np.linalg.norm(sum((k.conj().T@k for k in kraus), np.zeros((6, 6), complex))-np.eye(6), 2) > ADMISSION:
            raise ArithmeticError('complete native proper-code export failed')
        arrays[f'proper_code_{dimension}_complete_export'] = kraus
    return arrays


def lifted_repair_contract():
    """The recorded scalar must not be exponentiated by exterior degree."""
    weights = []
    for law, count in (('C', 4), ('E', 8)):
        code = native_current.transport.PAIR_CODES[law]
        rows = [code.branch_weight(law, native_current.transport.SEAMS[law], coin, raw)
                for coin in (0, 1) for raw in range(count)]
        if sum(rows) != 1: raise ArithmeticError('complete native repair weights do not sum to one')
        weights.append(rows)
    joint = [a*b for a in weights[0] for b in weights[1]]
    if len(joint) != 128 or sum(joint) != 1: raise ArithmeticError('joint native repair record menu changed')
    naive = [sum(w**degree for w in joint) for degree in (2, 4)]
    if naive[0] == naive[1] or any(value == 1 for value in naive):
        raise ArithmeticError('the naive full-Kraus exterior lift countercontrol disappeared')
    return {'complete_native_joint_record_weights_exact': list(map(str, joint)),
            'chosen_auxiliary_repair_lift': 'Ktilde_rs=sqrt(w_rs)*rho_spin_exterior(g_lift); lifted control path and scalar record weight retained',
            'lifted_control_path_or_cover_element_is_required': True,
            'six_by_six_endpoint_alone_determines_spin_cover_holonomy': False,
            'chosen_lift_completeness_on_all32_auxiliary_inputs': True,
            'naive_exterior_of_full_Kraus_degree2_and4_total_weights': list(map(str, naive)),
            'raw_Kraus_exterior_lift_is_trace_preserving': False,
            'literal_128_branch_auxiliary_native_repair_history_executed_here': False,
            'this_is_an_explicit_ideal_repair_lift_contract': True}


def decode_record_counts(counts):
    """Inverse takes complete counts and the source-derived observable dictionary.

    No executed amplitudes, action masses or species-name table enter. The
    selected rank15 uniform preparation and pointer calibration are inputs.
    """
    if (type(counts) is not np.ndarray or counts.dtype != np.dtype('int64')
            or counts.shape != (32,) or np.any(counts < 0) or np.any(counts > TRIALS)
            or counts.sum() != TRIALS):
        raise ValueError('all32 original nonnegative int64 pointer counts are required')
    # Simultaneous Hoeffding bounds for the five disjoint internal sectors.
    radius = math.sqrt(math.log(10/1e-9)/(2*TRIALS))+2*ADMISSION
    decoded = []
    for row in decoded_multiplets():
        number = int(counts[row['occupation_records']].sum())
        interval = (max(0., number/TRIALS-radius), min(1., number/TRIALS+radius))
        dimensions = [d for d in range(16) if interval[0] <= d/15 <= interval[1]]
        if len(dimensions) != 1: raise ArithmeticError('record counts do not identify an internal multiplicity')
        decoded.append({key: value for key, value in row.items() if key != 'dimension'} | {
            'count': number, 'probability_interval': list(interval), 'decoded_dimension': dimensions[0]})
    return decoded


def sample_complete_pointer_law(probabilities):
    """Fixed raw draws against exact cumulative source probabilities.

    Probability-dependent binomial sampling can change its random-number
    consumption under last-bit probability changes. Here the raw PCG64 stream
    is independent of the law. Rational cumulative thresholds avoid rounding
    their summation, and a fixed draw-to-boundary guard rejects fragile replay.
    This is a synthetic record sampler, not a hardware error certificate.
    """
    if (type(probabilities) is not np.ndarray or probabilities.dtype != np.dtype('float64')
            or probabilities.shape != (32,) or not np.isfinite(probabilities).all()
            or np.any(probabilities < 0)):
        raise ValueError('all32 finite nonnegative float64 source probabilities are required')
    weights = [Fraction(float(p)) for p in probabilities]
    total = sum(weights)
    if abs(total-1) > Fraction(ADMISSION):
        raise ArithmeticError('complete source record mass differs from one')
    grid = 1 << 64
    raw = np.random.PCG64(SEED).random_raw(TRIALS)
    draws = np.sort(raw)
    cumulative = Fraction(0); previous = 0; counts = []
    margin = Fraction(1)
    for weight in weights:
        cumulative += weight/total
        scaled = grid*cumulative
        threshold = (scaled.numerator+scaled.denominator-1)//scaled.denominator
        number = (0 if threshold <= 0 else TRIALS if threshold >= grid else
                  int(np.searchsorted(draws, np.uint64(threshold), side='left')))
        counts.append(number-previous); previous = number
        if 0 < cumulative < 1:
            for index in (number-1, number):
                if 0 <= index < TRIALS:
                    margin = min(margin, abs(Fraction(int(draws[index]), grid)-cumulative))
    if margin <= Fraction(SAMPLER_BOUNDARY_GUARD):
        raise ArithmeticError('source CDF lies too close to a fixed sampler draw')
    # Each CDF grid error is <1/grid, hence total variation is <32/grid.
    # It fits the second ADMISSION allowance already used by the inverse;
    # the first covers normalization. Neither statistical allowance changes.
    grid_error = Fraction(32, grid)
    if grid_error >= Fraction(ADMISSION):
        raise ArithmeticError('sampler grid exceeds the existing numerical allowance')
    return np.asarray(counts, dtype=np.int64), {
        'method': 'fixed_PCG64_uint64_draws_exact_rational_source_CDF',
        'raw_draw_sha256_little_endian': hashlib.sha256(raw.astype('<u8').tobytes()).hexdigest(),
        'minimum_draw_to_source_CDF_distance': float(margin),
        'fixed_boundary_distance_guard': SAMPLER_BOUNDARY_GUARD,
        'grid_total_variation_upper': float(grid_error),
        'fixed_seed_or_draws_selected_by_acceptance': False,
        'boundary_guard_failure_action': 'abort_without_resampling',
        'statistical_coverage_conditioned_on_boundary_guard_pass': False,
        'source_probabilities_replaced_by_ideal_multiplicities': False,
        'hardware_or_outward_rounding_error_certificate': False}


def record_identification():
    r = representation(); kraus, pointer_word = occupation_instrument()
    initial = np.asarray(r['encoder'], complex)/math.sqrt(15)
    outputs = np.stack([k@initial for k in kraus])
    probabilities = np.sum(abs(outputs)**2, axis=(1, 2))
    total = float(probabilities.sum())
    if abs(total-1) > ADMISSION: raise ArithmeticError('complete source record mass differs from one')
    counts, sampling = sample_complete_pointer_law(probabilities)
    decoded = decode_record_counts(counts)
    selected = set(r['selected'])
    excluded = [mask for mask in range(32) if mask not in selected]
    return {'trials': TRIALS, 'sampler_seed': SEED, 'simultaneous_failure_probability_upper': 1e-9,
            'synthetic_record_sampler': sampling,
            'full_outcome_count': 32, 'pointer_sourceword_local': pointer_word,
            'selected_uniform_register_is_architecture': True,
            'excluded_outcomes_retained': excluded, 'excluded_outcome_count': int(counts[excluded].sum()),
            'all_outcomes_normalized_together_for_sampling': True,
            'normalization_total_variation_allowance': ADMISSION,
            'species_name_table_used_by_decoder': False,
            'decoded_internal_multiplets': decoded}, {
                'occupation_instrument': kraus, 'occupation_outputs': outputs,
                'occupation_counts': counts, 'occupation_complete_probabilities': probabilities}


def build():
    r = representation(); mapping, vertices = action_intertwiner()
    executions, arrays = current_execution()
    records, record_arrays = record_identification()
    selected = list(r['selected']); q = r['operators']['q'].extract(selected, selected)/6
    c3 = r['operators']['color_C2'].extract(selected, selected); c2 = r['operators']['weak_C2'].extract(selected, selected)
    anomaly = [s.trace(q), s.trace(q**3), s.trace(q*c3), s.trace(q*c2)]
    if any(value != 0 for value in anomaly): raise ArithmeticError('selected matter anomaly traces failed')
    report = {'schema': 'oph.native_matter_attachment.v1',
              'architecture': {'native_port_carrier_changed': False,
                   'auxiliary_encoded_five_mode_register': True,
                   'spin_half_lift_and_central_trace_redistribution_are_chosen': True,
                   'nonvacuum_even_exterior_selection_is_chosen': True,
                   'physical_Weyl_fermion_statistics_inferred_from_auxiliary_even_parity': False,
                   'ideal_continuous_native_current_service_executed': True,
                   'finite_64_row_native_pulse_leakage_replaced_by_zero': False,
                   'entire_matter_word_executed_on_finite_64_row_native_pulse_backend': False},
              'same_current_lift': {'native_current_count': 12, 'spin_bracket_checks': 66,
                   'exterior_bracket_checks': 66, 'full_Fock_dimension': 32,
                   'selected_matter_dimension': 15, 'primitive_integer_charge5': entries(r['primitive_charge5']),
                   'source_current_order_to_existing_action': list(r['native']['native_to_fixture_port']),
                   'all12_current_intertwiner_exact': True,
                   'source_weak_2pi_is_nontrivial_on_spin_lift': True,
                   'source_weak_4pi_is_identity_on_spin_lift': True},
              'derived_internal_multiplets': decoded_multiplets(),
              'source_current_executions': executions, 'record_identification': records,
              'same_current_Higgs_action': {'invariant_lines': vertices,
                   'all72_current_component_identities_exact': True,
                   'same_kinetic_action_intertwiner': 'I_spin tensor(P15 tensor I_family), with corresponding Nambu conjugate',
                   'internal_frame_change_alters_the_existing_mass_spectrum': False,
                   'mass_parameters_or_family_count_selected_by_the_current_algebra': False},
              'same_kinetic_and_mass_action': kinetic_attachment(),
              'relative_spectrum_from_same_current_constraints': relative_spectrum_freedom(),
              'anomaly_traces': list(map(str, anomaly)), 'empirical_mass_or_coupling_inputs': [],
              'global_cover_contract': global_cover_contract(),
              'lifted_repair_contract': lifted_repair_contract(),
              'physical_field_attachment': {
                   'separate_one_particle_space': 'C2_WeylSpin tensor M15 tensor C3_family',
                   'physical_CAR_mode_count': 90,
                   'physical_annihilator_index': 'spin*45+3*matter_mask_index+family',
                   'physical_CAR_is_a_separate_architecture_functor': True,
                   'physical_creation_exchange_sign': '-1 independently of auxiliary even parity',
                   'same_current_field_intertwiner_preserves_CAR_exactly': True,
                   'generation_count_three_is_chosen': True,
                   'physical_spatial_Spin_is_internal_weak_SU2': False},
              'scope': {'actual_native_current_to_executed_matter_record_bridge': True,
                   'physical_local_gauge_redundancy_derived': False,
                   'physical_Spin_and_statistics_selected_by_repair': False,
                   'quantum_gauge_Higgs_loops_executed': False,
                   'measured_particle_mass_spectrum_predicted': False}}
    arrays.update(record_arrays)
    arrays.update(proper_code_exports())
    arrays['matter_to_existing_action_intertwiner'] = np.asarray(mapping, float)
    arrays['native_matter_generators'] = np.asarray(r['matter_currents'], complex)
    return report, arrays


def source_bindings():
    paths = set(native_current.source_bindings()) | set(existing_action.SOURCE_FILES) | {
        'model/common_source/native_matter_attachment.py', 'model/common_source/test_native_matter_attachment.py',
        'model/common_source/chiral_candidate/action.py', 'model/common_source/chiral_candidate/gaussian_source.py',
        'model/observer/source_scalar/compiler.py', 'tools/strict_json.py'}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def _compare(a, b, path='receipt'):
    if type(a) is not type(b): raise ValueError(path+' type differs')
    if type(a) is dict:
        if set(a) != set(b): raise ValueError(path+' keys differ')
        for key in a: _compare(a[key], b[key], path+'.'+key)
    elif type(a) is list:
        if len(a) != len(b): raise ValueError(path+' length differs')
        for i, (x, y) in enumerate(zip(a, b, strict=True)): _compare(x, y, path+f'[{i}]')
    elif type(a) is float:
        if not math.isfinite(b) or abs(a-b) > 1e-12: raise ValueError(path+' finite value differs')
    elif a != b: raise ValueError(path+' differs')


def check(path):
    path = Path(path); pins = source_bindings(); report = strict_json.load(path)
    if type(report) is not dict or report.get('source_bindings') != pins: raise ValueError('matter source bindings differ')
    artifact = path.with_suffix('.npz')
    if report.get('artifact') != {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}:
        raise ValueError('matter artifact custody differs')
    fresh, arrays = build()
    if source_bindings() != pins: raise ValueError('matter sources changed during replay')
    _compare(fresh, {k: v for k, v in report.items() if k not in ('source_bindings', 'artifact')})
    with np.load(artifact, allow_pickle=False) as stored:
        if len(stored.files) != len(set(stored.files)) or set(stored.files) != set(arrays): raise ValueError('complete matter array keys differ')
        for key, value in arrays.items():
            old = stored[key]
            if old.shape != value.shape or old.dtype != value.dtype or not np.isfinite(old).all(): raise ValueError('matter array schema differs')
            if value.dtype.kind in 'iu':
                if not np.array_equal(old, value): raise ValueError('matter counts differ')
            elif np.max(abs(old-value), initial=0) > 3e-10: raise ValueError('matter source array differs')


def main():
    parser = argparse.ArgumentParser(description=__doc__); mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--output', type=Path); mode.add_argument('--check', type=Path); args = parser.parse_args()
    if args.check:
        check(args.check); print('Actual native current, exterior matter, source records and Higgs intertwiner: PASS'); return
    artifact = args.output.with_suffix('.npz')
    if args.output.exists() or artifact.exists(): parser.error('choose fresh receipt and artifact paths')
    pins = source_bindings(); report, arrays = build()
    if source_bindings() != pins: raise ValueError('matter sources changed during generation')
    with artifact.open('xb') as stream: np.savez_compressed(stream, **arrays)
    report['source_bindings'] = pins
    report['artifact'] = {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    with args.output.open('x') as stream: json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False); stream.write('\n')


if __name__ == '__main__': main()
