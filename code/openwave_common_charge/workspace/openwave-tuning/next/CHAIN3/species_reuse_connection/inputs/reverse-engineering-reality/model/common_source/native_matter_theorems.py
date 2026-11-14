"""Map the actual current's exterior components to the finite matter theorem.

The five-mode lift and component menu are declared. Within that menu,
nonempty chiral anomaly-free selection forces one of two conjugate sectors.
The calculation does not derive those selection criteria, exclude singlets
outside the menu, or fix how many copies occur in a physical field theory.
"""
from functools import lru_cache
import hashlib
from pathlib import Path
import re

import sympy as s

from . import native_matter_attachment as matter

ROOT = Path(__file__).resolve().parents[2]


def source_bindings():
    pins = dict(matter.source_bindings())
    for name in ('model/common_source/native_matter_theorems.py',
                 'model/common_source/test_native_matter_theorems.py',
                 'Lean/Screen/ExteriorComponentBridge.lean',
                 'Lean/Screen/ExteriorSelection.lean',
                 'Lean/Screen/QuantumMatterIntegration.lean'):
        pins[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return pins


def _lean_table(source, name, kind=int):
    match = re.search(r'\bdef '+re.escape(name)+r'\s*:.*?:=\s*!\[([^]]+)\]', source, re.S)
    if match is None:
        raise ValueError('finite theorem table missing: '+name)
    words = [word.strip() for word in match.group(1).split(',')]
    if kind is bool:
        if any(word not in ('true', 'false') for word in words):
            raise ValueError('invalid Boolean theorem table')
        return tuple(word == 'true' for word in words)
    return tuple(kind(word) for word in words)


def component_map():
    """Read every component from the full source action, including complements."""
    r = matter.representation(); ops = r['operators']
    groups = {}
    for mask in range(32):
        degree = ((mask & 7).bit_count(), (mask >> 3).bit_count())
        groups.setdefault(degree, []).append(mask)
    bridge = (ROOT/'Lean/Screen/ExteriorComponentBridge.lean').read_text(encoding='utf-8')
    match = re.search(r'def componentDegree.*?:=\s*!\[([^]]+)\]', bridge, re.S)
    if match is None: raise ValueError('finite component-degree table missing')
    degrees = tuple(tuple(map(int, pair)) for pair in re.findall(r'\((\d+),\s*(\d+)\)', match.group(1)))
    if len(degrees) != 10 or len(set(degrees)) != 10:
        raise ValueError('finite component-degree menu changed')
    if set(degrees) | {(0, 0), (3, 2)} != set(groups):
        raise ArithmeticError('the complete32 source basis does not match the theorem menu')
    rows = []
    color_cartan = ops['H1']+2*ops['H2']
    for index, (a, b) in enumerate(degrees):
        masks = groups[a, b]; E = s.eye(32)[:, masks]
        for g in r['full_currents']:
            if not matter.zero(g*E-E*(E.T*g*E)):
                raise ArithmeticError('a source current leaves its exterior component')
        charge = ops['q'][masks[0], masks[0]]
        if E.T*ops['q']*E != charge*s.eye(len(masks)):
            raise ArithmeticError('source component charge is not scalar')
        color_dim, weak_dim = s.binomial(3, a), s.binomial(2, b)
        if len(masks) != color_dim*weak_dim:
            raise ArithmeticError('source component multiplicity changed')
        for name, expected in [('color_C2', s.Rational(4, 3) if color_dim == 3 else 0),
                               ('weak_C2', s.Rational(3, 4) if weak_dim == 2 else 0)]:
            if E.T*ops[name]*E != expected*s.eye(len(masks)):
                raise ArithmeticError('source Casimir does not match the finite component')
        conjugate = degrees.index((3-a, 2-b))
        anomalies = (len(masks)*charge,
                     weak_dim*charge if color_dim == 3 else 0,
                     color_dim*charge if weak_dim == 2 else 0,
                     len(masks)*charge**3)
        rows.append({'index': index, 'degree': [a, b], 'occupation_masks': masks,
                     'dimension': len(masks), 'color_dimension': int(color_dim),
                     'weak_dimension': int(weak_dim), 'charge': int(charge),
                     'is_colored': bool(color_dim == 3), 'is_doublet': bool(weak_dim == 2),
                     'odd': bool((a+b) % 2), 'conjugate': conjugate,
                     'four_anomaly_forms': list(map(int, anomalies)),
                     'pure_color_cubic': int(-s.trace((E.T*color_cartan*E)**3)/6),
                     'weak_doublet_count': int(color_dim if weak_dim == 2 else 0)})
    source = (ROOT/'Lean/Screen/ExteriorSelection.lean').read_text(encoding='utf-8')
    for name, key, kind in [('colorDim', 'color_dimension', int), ('weakDim', 'weak_dimension', int),
                            ('charge', 'charge', int), ('odd', 'odd', bool), ('conj', 'conjugate', int),
                            ('isColored', 'is_colored', bool), ('isDoublet', 'is_doublet', bool)]:
        if tuple(row[key] for row in rows) != _lean_table(source, name, kind):
            raise ArithmeticError('actual source differs from finite theorem column: '+name)
    return rows


def admissible_masks(rows):
    """Exhaustive theorem premises on the supplied component menu."""
    if len(rows) != 10 or [r['index'] for r in rows] != list(range(10)):
        raise ValueError('the complete ten-component source menu is required')
    result = []
    for mask in range(1, 1 << 10):
        chosen = [row for row in rows if mask >> row['index'] & 1]
        if any(mask >> row['conjugate'] & 1 for row in chosen): continue
        if any(sum(row['four_anomaly_forms'][j] for row in chosen) for j in range(4)): continue
        result.append(mask)
    return result


@lru_cache(maxsize=1)
def run():
    rows = component_map(); masks = admissible_masks(rows)
    if masks != [227, 796]:
        raise ArithmeticError('source no longer realizes finite two-sector selection')
    sectors = []
    for mask in masks:
        selected = [row for row in rows if mask >> row['index'] & 1]
        occupations = sorted(i for row in selected for i in row['occupation_masks'])
        parity = {row['odd'] for row in selected}
        if len(parity) != 1 or len(occupations) != 15:
            raise ArithmeticError('source sector lost parity or rank15')
        color = sum(row['pure_color_cubic'] for row in selected)
        weak = sum(row['weak_doublet_count'] for row in selected)
        if color or weak % 2:
            raise ArithmeticError('additional source anomaly check failed')
        if mask == 796 and occupations != list(matter.representation()['selected']):
            raise ArithmeticError('selected matter register differs from theorem sector')
        sectors.append({'component_mask': mask, 'occupation_masks': occupations,
                        'dimension': len(occupations), 'odd_auxiliary_degree': parity.pop(),
                        'four_anomaly_forms': [0]*4, 'pure_color_cubic': color,
                        'weak_doublet_count': weak})
    return {'classification': 'conditional actual-current instance of the finite exterior selection',
            'components': rows, 'source_current_component_invariance_checks': 120,
            'exhaustive_component_subsets': 1024, 'admissible_sectors': sectors,
            'excluded_menu_singlet_masks_retained': [0, 31],
            'singlets_removed_by_anomaly_or_chirality_test': False,
            'copy_count_constraint_from_anomalies': 'none: A(N copies)=N A(one copy)=0',
            'three_families_derived': False,
            'supplied_selection_criteria': ['nonempty', 'chiral', 'four finite anomaly forms vanish'],
            'supplied_representation': 'same-current five-mode spin lift and non-vacuum/non-top exterior component menu',
            'physical_statistics_or_continuum_anomaly_derivation': False,
            'finite_theorems': ['Lean/Screen/ExteriorComponentBridge.lean',
                                'Lean/Screen/ExteriorSelection.lean',
                                'Lean/Screen/QuantumMatterIntegration.lean'],
            'empirical_mass_inputs': []}
