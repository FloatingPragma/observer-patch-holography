"""Compute the fifteen internal weights before applying any particle names."""
from functools import lru_cache

import sympy as s

from model.common_source import native_matter_attachment as matter


@lru_cache(maxsize=1)
def source_components():
    r = matter.representation()
    ops = r['operators']
    rows = []
    for mask in r['selected']:
        q = ops['q'][mask, mask]
        t = ops['T3'][mask, mask]
        rows.append({'occupation_mask': mask,
                     'q6Y': str(q), 'Y': str(q/6), 'T3': str(t),
                     'Q': str(t+q/6),
                     'color_C2': str(ops['color_C2'][mask, mask]),
                     'weak_C2': str(ops['weak_C2'][mask, mask]),
                     'color_weight': [str(ops[k][mask, mask]) for k in ('H1', 'H2')]})
    return tuple(rows)


# These are comparison labels and exact representation weights, not source
# inputs, experimental masses, or derived spacetime chiralities.
TARGET_COMPONENTS = {
    'up_left': ('1', '1/2', '4/3', '3/4', '3'),
    'down_left': ('1', '-1/2', '4/3', '3/4', '3'),
    'up_right_conjugate': ('-4', '0', '4/3', '0', '3bar'),
    'down_right_conjugate': ('2', '0', '4/3', '0', '3bar'),
    'neutrino_left': ('-3', '1/2', '0', '3/4', '1'),
    'charged_lepton_left': ('-3', '-1/2', '0', '3/4', '1'),
    'charged_lepton_right_conjugate': ('6', '0', '0', '0', '1'),
}


def match_components(rows=None):
    """Require every source state once, including color conjugation weights."""
    rows = source_components() if rows is None else rows
    if len(rows) != 15 or len({r['occupation_mask'] for r in rows}) != 15:
        raise ValueError('the complete fifteen-state source is required')
    weights = {'3': {('1', '0'), ('-1', '1'), ('0', '-1')},
               '3bar': {('-1', '0'), ('1', '-1'), ('0', '1')},
               '1': {('0', '0')}}
    result, used = {}, []
    for name, (q, t, c, w, rep) in TARGET_COMPONENTS.items():
        selected = [r for r in rows if tuple(r[k] for k in
                    ('q6Y', 'T3', 'color_C2', 'weak_C2')) == (q, t, c, w)]
        if (len(selected) != len(weights[rep]) or
                {tuple(r['color_weight']) for r in selected} != weights[rep]):
            raise ValueError('source weight multiplicity or color conjugation differs: '+name)
        if any(s.Rational(r['Q']) != s.Rational(q)/6+s.Rational(t) or
               s.Rational(r['Y']) != s.Rational(q)/6 for r in selected):
            raise ValueError('electric or hypercharge readout differs')
        masks = sorted(r['occupation_mask'] for r in selected)
        used.extend(masks)
        result[name] = {'source_masks': masks, 'color_representation': rep,
                        'q6Y': q, 'T3': t, 'Q': str(s.Rational(q)/6+s.Rational(t)),
                        'comparison_convention': 'all-left-Weyl internal representation',
                        'physical_chirality_derived': False}
    if sorted(used) != sorted(r['occupation_mask'] for r in rows):
        raise ValueError('unmatched or multiply assigned source components')
    return result


def gauge_and_higgs():
    """Check the actual weak raising/lowering charge and neutral Higgs line."""
    r = matter.representation()
    indices = list(r['selected'])
    raw = matter.canonical_currents()
    weak5 = tuple(s.I*matter.spin_lift(g) for g in raw[8:11])
    weak = tuple(matter.exterior_derivative(g).extract(indices, indices) for g in weak5)
    y = r['operators']['q'].extract(indices, indices)/6
    q = y+weak[2]
    plus, minus = weak[0]+s.I*weak[1], weak[0]-s.I*weak[1]
    if (q*plus-plus*q != plus or q*minus-minus*q != -minus or
            plus.H != minus or plus.rank() != 4 or minus.rank() != 4):
        raise ArithmeticError('same-current charged weak directions differ')
    color = tuple(matter.exterior_derivative(s.I*matter.spin_lift(g)).extract(indices, indices)
                  for g in raw[:8])
    if any(q*g-g*q != s.zeros(15) for g in color):
        raise ArithmeticError('a color current has nonzero electric charge')
    if s.Matrix.hstack(*(s.Matrix(g).reshape(225, 1) for g in color)).rank() != 8:
        raise ArithmeticError('the source color algebra is not eight dimensional')
    h_y = r['primitive_charge5'][3:, 3:]/6
    h_q = h_y+weak5[2][3:, 3:]
    if h_y != s.eye(2)/2 or h_q != s.diag(1, 0):
        raise ArithmeticError('the declared Higgs doublet charge differs')
    return {'independent_color_current_directions': 8,
            'color_current_electric_charge': '0',
            'weak_charged_direction_charges': ['1', '-1'],
            'weak_raising_and_lowering_ranks_on_matter': [4, 4],
            'neutral_current_span': ['T3', 'Y'],
            'unbroken_generator_for_declared_neutral_Higgs_direction': 'Q=T3+Y',
            'Higgs_doublet_Y': '1/2', 'Higgs_component_Q': ['1', '0'],
            'Higgs_vacuum_selected_by_source': False,
            'weak_mixing_angle_selected_by_source': False,
            'adjoint_directions_identified_with_propagating_bosons': False}
