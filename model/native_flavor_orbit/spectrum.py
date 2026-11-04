"""Certified joint dimensionless spectra and quark mixing of source duties."""
from fractions import Fraction as F
from functools import lru_cache
import sympy as s

from . import algebra, intervals as iv


def eigenvalue_boxes(matrix, rank):
    """Deterministic rational bisection, certified by every LDL inertia."""
    result = []
    for index in range(3):
        if index < 3-rank:
            result.append(iv.I(0))
            continue
        lo, hi = F(0), F(1)
        for _ in range(55):
            middle = (lo+hi)/2
            if iv.inertia_below(matrix, middle) <= index:
                lo = middle
            else:
                hi = middle
        if (lo <= 0 or iv.inertia_below(matrix, lo) != index
                or iv.inertia_below(matrix, hi) != index+1):
            raise ArithmeticError('no certified simple positive eigenvalue bracket')
        result.append(iv.I(lo, hi))
    if any(a.hi >= b.lo for a, b in zip(result, result[1:])):
        raise ArithmeticError('the complete source eigenvalue brackets overlap')
    return result


def spectral_projector(matrix, eigenvalues, index):
    result = [[iv.I(int(i == j)) for j in range(3)] for i in range(3)]
    for j, other in enumerate(eigenvalues):
        if j != index:
            factor = iv.shifted(matrix, other)
            divisor = eigenvalues[index]-other
            result = [[x/divisor for x in row] for row in iv.matmul(result, factor)]
    return result


def mixing_boxes(up, down, up_values, down_values):
    left = [spectral_projector(up, up_values, i) for i in range(3)]
    right = [spectral_projector(down, down_values, i) for i in range(3)]
    result = []
    for a in left:
        row = []
        for b in right:
            product = iv.matmul(a, b)
            value = sum((product[i][i] for i in range(3)), iv.I(0))
            if value.hi < 0 or value.lo > 1:
                raise ArithmeticError('a spectral overlap missed the probability domain')
            row.append(iv.rounded_out(value, 9, probability=True))
        result.append(row)
    return result


@lru_cache(maxsize=1)
def catalogue():
    data = algebra.symbolic_memory()
    gates = algebra.source_gates()
    exact_gates = {}
    gate_boxes = {}
    for name, row in gates.items():
        gate_boxes[name] = {key: iv.expression(row[key]) for key in ('p', 'q')}
        exact_gates[name] = {'p': str(row['p']), 'q': str(row['q']),
                             'p_interval': iv.rounded_out(gate_boxes[name]['p'], 16).strings(),
                             'q_interval': iv.rounded_out(gate_boxes[name]['q'], 16).strings(),
                             'left_Casimir': str(row['left_casimir']),
                             'right_Casimir': str(row['right_casimir'])}
    if not gate_boxes['up']['q'].lo > gate_boxes['down']['q'].hi:
        raise ArithmeticError('the declared up and down current gates coincide')
    policies = []
    for policy in algebra.POLICIES:
        matrices, values = {}, {}
        rank = 3 if policy in ('always', 'C') else 2
        channels = []
        full_squared_rows = []
        for name in ('up', 'down', 'charged_lepton'):
            substitutions = {data['p']: gate_boxes[name]['p'], data['q']: gate_boxes[name]['q']}
            matrices[name] = iv.matrix(data['coefficients'][policy], substitutions)
            values[name] = eigenvalue_boxes(matrices[name], rank)
            multiplicity = 2 if name == 'charged_lepton' else 6
            mass_intervals = [iv.rounded_out(x, 12).strings() for x in values[name]]
            channels.append({'channel': name, 'rank': rank,
                             'dimensionless_singular_values': mass_intervals,
                             'K45_multiplicity_per_family': multiplicity,
                             'unconditional_trace': iv.rounded_out(
                                 sum((matrices[name][i][i] for i in range(3)), iv.I(0)), 12).strings()})
            for family, value in enumerate(values[name]):
                full_squared_rows.extend([{'channel': name, 'family': family+1,
                                           'squared_singular_value': iv.rounded_out(value**2, 12).strings()}]
                                         *multiplicity)
        full_squared_rows.extend([{'channel': 'minimal_neutrino', 'family': i+1,
                                   'squared_singular_value': ['0', '0']} for i in range(3)])
        if len(full_squared_rows) != 45:
            raise ArithmeticError('the complete action spectrum lost an internal row')
        mixing = mixing_boxes(matrices['up'], matrices['down'], values['up'], values['down'])
        if policy == 'always' and not any(mixing[i][j].lo > 0 for i in range(3)
                                          for j in range(3) if i != j):
            raise ArithmeticError('the certified source mixing is trivial')
        policies.append({'policy': policy, 'channels': channels,
                         'full_K45_squared_spectrum': full_squared_rows,
                         'full_K45_rank': 14*rank, 'full_K45_nullity': 45-14*rank,
                         'quark_mixing_modulus_squared': [[x.strings() for x in row] for row in mixing],
                         'charged_current_CP_commutator_trace_cube_exact': '0',
                         'mixing_convention': 'ascending source singular values; |U_up^dagger U_down|^2',
                         'units': 'stationary first-order coefficient per declared pulse step; Higgs norm1',
                         'physical_particle_poles_or_empirical_fit': False})
    return {'declared_angle_exact': str(s.pi/5), 'declared_reset_port': 2,
            'declared_alternative_angles': [str(s.pi/5), str(s.pi/3), str(s.pi/2)],
            'source_gates': exact_gates, 'policies': policies,
            'enclosure_method': 'rational Machin and alternating sine bounds; interval LDL inertia; spectral-projector overlaps',
            'charged_current_CP_invariant_is_basis_independent': True,
            'full_45_row_action_counts_include_all_zero_directions': True,
            'joint_action_uses_declared_static_mean_feedback_service': True,
            'joint_action_equals_finite_correlated_memory_channel': False,
            'independent_source_selection_or_physical_clock_supplied': False}
