"""Exact full local-Fock colour census and uniform open-chain inequalities."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import comb, isqrt
import json
from pathlib import Path
import sympy as sy
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.openwave_spin_gradient_nuclei import source as native

HERE = Path(__file__).resolve().parent
INPUTS = json.loads((HERE / 'inputs.json').read_text())
SELECTED_KAPPA = F(INPUTS['selected_kappa'])
ELECTRON_GEV = F(INPUTS['electron_unit_GeV'])
QUADRUPOLE_FM2 = F(INPUTS['quadrupole_comparison_fm_squared'])
ILLUSTRATIVE_TENSION_GEV_PER_FM = F(INPUTS['illustrative_tension_GeV_per_fm'])


@lru_cache(None)
def local_colour_census():
    """Decompose all 2**24 native colourful occupations, without enumerating them."""
    subsets, _, generators, _ = representation()
    qids = [[subsets.index((c, w)) for c in range(3)] for w in (3, 4)]
    aids = [[subsets.index(tuple(c for c in range(3) if c != missing))
             for missing in range(3)],
            [subsets.index(tuple(c for c in range(5) if c != missing))
             for missing in range(3)]]
    frames = [ids for family in (qids, aids) for spin in range(2) for ids in family]
    colour = list(generators[:7]) + [(generators[6] + 2 * generators[7]) / sy.sqrt(3)]
    weights = []
    for ids in frames:
        matrices = [g.extract(ids, ids) for g in colour]
        assert sum((g * g for g in matrices), sy.zeros(3)) == sy.Rational(4, 3) * sy.eye(3)
        for i in ids:
            w = (2 * generators[6][i, i], 2 * generators[7][i, i])
            assert all(x.is_Integer for x in w)
            weights.append(tuple(map(int, w)))
    assert sorted(weights) == sorted([(1, 0), (-1, 1), (0, -1)] * 4
                                     + [(-1, 0), (1, -1), (0, 1)] * 4)
    character = {(0, 0, 0): 1}
    for u, v in weights:
        output = dict(character)
        for (n, a, b), multiplicity in character.items():
            key = n + 1, a + u, b + v
            output[key] = output.get(key, 0) + multiplicity
        character = output
    roots = ((2, -1), (-1, 2), (1, 1))
    denominator = [(sum(x for x, _ in group), sum(y for _, y in group), (-1) ** n)
                   for n in range(4) for group in combinations(roots, n)]
    rows, dimensions, central_labels = [], [], set()
    for n in range(25):
        dimension = 0
        for a, b in product(range(17), repeat=2):
            multiplicity = sum(sign * character.get((n, a + x, b + y), 0)
                               for x, y, sign in denominator)
            assert multiplicity >= 0
            if not multiplicity:
                continue
            dim = (a + 1) * (b + 1) * (a + b + 2) // 2
            casimir = F(a * a + b * b + a * b + 3 * a + 3 * b, 3)
            cubic = F((a - b) * (2 * a + b + 3) * (a + 2 * b + 3), 18)
            assert (n, casimir, cubic) not in central_labels
            central_labels.add((n, casimir, cubic))
            assert casimir <= 3 * n
            rows.append({'number': n, 'highest_weight': [a, b], 'multiplicity': multiplicity,
                         'irrep_dimension': dim, 'Casimir': str(casimir),
                         'cubic_Casimir': str(cubic),
                         'three_number_minus_Casimir': str(3 * n - casimir)})
            dimension += dim * multiplicity
        assert dimension == comb(24, n)
        dimensions.append(dimension)
    assert sum(dimensions) == 2 ** 24 and len(rows) == 547
    sharp = [r for r in rows if r['number'] and F(r['Casimir']) == 3 * r['number']]
    assert len(sharp) == 1 and sharp[0]['number'] == 8 and sharp[0]['highest_weight'] == [4, 4]
    return {'local_modes': 24, 'full_local_Fock_dimension': 2 ** 24,
            'number_sector_dimensions': dimensions, 'irreducible_rows': rows,
            'sharp_nonvacuum_equality_rows': sharp, 'operator_bound': 'C_x^2 <= 3 N_x',
            'native_one_body_Casimir_exact': '4/3',
            'domain': 'One colourful weak doublet and its conjugates, both spins; all local occupations and colour representations.'}


def sqrt_interval(value, digits=18):
    """Rational enclosure using integer arithmetic only."""
    value = F(value)
    assert value >= 0
    scale = 10 ** digits
    floor = isqrt(value.numerator * scale * scale // value.denominator)
    lower, upper = F(floor, scale), F(floor + 1, scale)
    assert lower * lower <= value <= upper * upper
    return lower, upper


def uniform_bound(kappa):
    p = native.parameters()
    kappa = F(kappa)
    assert kappa >= 0
    m0 = p['mbar'] - abs(p['delta']) / 2
    floor = m0 - 3 * max(F(0), kappa / 4 - p['U'])
    maximum = 4 * p['U'] + 4 * m0 / 3
    return {'parameters': p, 'kappa': kappa, 'number_floor': floor,
            'strict_positive_kappa_upper': maximum,
            'adjusted_local_colour_coefficient': p['U'] - kappa / 2}


def shared_length_attachment():
    """A rational common length; both comparison targets have explicit ancestry."""
    parent = json.loads((HERE.parent / 'openwave_spin_gradient_nuclei/results.json').read_text())
    qlo, qhi = map(F, parent['selected_nuclear_record']['quadrupole_graph_squared_interval'])
    lower = sqrt_interval(QUADRUPOLE_FM2 / qhi)[0] * 3 / (4 * ELECTRON_GEV)
    upper = sqrt_interval(QUADRUPOLE_FM2 / qlo)[1] * 3 / (4 * ELECTRON_GEV)
    a = 4 * SELECTED_KAPPA * ELECTRON_GEV / (3 * ILLUSTRATIVE_TENSION_GEV_PER_FM)
    predicted = (a * a * qlo, a * a * qhi)
    assert lower <= SELECTED_KAPPA <= upper
    assert predicted[0] <= QUADRUPOLE_FM2 <= predicted[1]
    bound = uniform_bound(upper)
    assert bound['number_floor'] > 0
    return {'kappa': SELECTED_KAPPA, 'kappa_required_interval': (lower, upper),
            'a_fm': a, 'quadrupole_graph_interval': (qlo, qhi),
            'quadrupole_fm2_interval': predicted,
            'tension_GeV_per_fm': ILLUSTRATIVE_TENSION_GEV_PER_FM,
            'electron_unit_GeV': ELECTRON_GEV,
            'entire_required_interval_number_floor': bound['number_floor'],
            'ancestry': 'The five strong calibrations are unchanged. The quadrupole comparison and illustrative 1 GeV/fm tension jointly select this retrospective common-length attachment; neither is an independent prediction.'}


def incidence_certificates():
    """An exact factorization, valid for every finite chain and each colour axis."""
    rows = []
    for length in range(2, 18):
        incidence = sy.zeros(length, length - 1)
        for x in range(length - 1):
            incidence[x, x] = 1
            incidence[x + 1, x] = -1
        plus = sy.zeros(length, length - 1)
        for x in range(length - 1):
            plus[x, x] = 1
            plus[x + 1, x] = 1
        assert 4 * sy.eye(length - 1) - incidence.T * incidence == plus.T * plus
        rows.append({'cells': length, 'edge_count': length - 1,
                     'four_identity_minus_incidence_Gram_equals_plus_Gram': True})
    return rows


def three_quark_control():
    rows = []
    for length in range(2, 18):
        minimum, maximum = None, None
        for positions in combinations(range(length + 2), 3):
            # Stars-and-bars ordering: all multisets of three sites.
            positions = tuple(x - i for i, x in enumerate(positions))
            n = [positions.count(x) for x in range(length)]
            span = max(positions) - min(positions)
            local = sum(F(2, 3) * x * (3 - x) for x in n)
            remainder = F(4, 3) * span - local / 2
            assert remainder >= 0
            minimum = remainder if minimum is None else min(minimum, remainder)
            maximum = remainder if maximum is None else max(maximum, remainder)
        rows.append({'cells': length, 'minimum_kappa_remainder': str(minimum),
                     'maximum_kappa_remainder': str(maximum)})
    # Three-site N=6 singlet: 3 tensor bar6 tensor 3 contains a singlet.
    # Both prefix charges are 3 or bar3; local Casimirs are4/3,10/3,4/3.
    assert F(8, 3) - (F(4, 3) + F(10, 3) + F(4, 3)) / 2 == -F(1, 3)
    return {'all_three_quark_multisets': rows,
            'six_quark_counterexample': {'local_numbers': [1, 4, 1],
                'local_highest_weights': [[1, 0], [0, 2], [1, 0]],
                'local_Casimir_sum': '6', 'prefix_Casimir_sum': '8/3',
                'electric_minus_half_local_Casimir': '-1/3',
                'scope': 'The colour-only kappa remainder can be negative; the complete number-plus-colour positivity theorem retains this sector.'}}
