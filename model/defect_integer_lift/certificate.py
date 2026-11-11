"""Assemble the integer-lift certificate: closed carriers, the exact gates on a reference run, kinematics of the
lift, the readings, and arithmetic on the imported run counts."""
from collections import Counter
from itertools import combinations_with_replacement, product
import hashlib
from pathlib import Path
import random

from . import model as m
from .census import census

HERE = Path(__file__).parent
SCHEMA = 'oph.defect_integer_lift.v1'
REFERENCE = dict(level=1, sweeps=200, cuts=(1, 2, 4, 16, 64, 512))
NAMES = {1: 'Q', -1: 'anti-Q', 4: 'u_R', -2: 'd_R', -4: 'anti-u_R', 2: 'anti-d_R', -3: 'L, anti-H', 3: 'anti-L, H',
         -6: 'e_R', 6: 'anti-e_R', 0: 'nu_R, G, W, B'}
COMPOSITIONS = ((2, 4), (1, 5), (3, 3), (2, 2, 4, 4), (4, 4, 4), (2, 2, 2), (2, 3, 3, 4), (1, 2, 3), (3, 4, 5))


def structure():
    rows, checks = {}, {}
    for level in (0, 1, 2):
        s = m.Federation(level).structure()
        rows[str(level)] = s
        ok = (s['euler_characteristic'] == 2 and s['every_directed_seam_once'] and s['port_degrees'] == [5]
              and s['every_seam_has_two_distinct_faces'] and s['internal_seam_faces_on_one_carrier']
              and s['glue_seam_faces_are_plaquettes'] and s['ports_glued_once'] and s['glue_seams_per_carrier'] == [3]
              and s['glue_seams_join_two_carriers'] and set(s['seams_around_a_plaquette']) <= {5, 6}
              and s['screen_euler_characteristic'] == 2)
        checks[f'closed_oriented_carriers_and_screen_level_{level}'] = ok
    return rows, checks


def reference_run():
    fed = m.Federation(REFERENCE['level'])
    per_sweep = []
    full = m.run(fed, REFERENCE['sweeps'],
                 after_sweep=lambda s, x, reg, cls, link, val: per_sweep.append(
                     all(m.gates(fed, dict(values=val, classes=cls, link=link, register=reg)).values())))
    end = m.gates(fed, full)
    twists = full['twists']
    cuts = {}
    for depth in REFERENCE['cuts']:
        cut = m.run(fed, REFERENCE['sweeps'], max_depth=depth)
        prefix = [t for t in twists if t[4] <= depth]
        link = [0] * fed.seams
        for _, _, seam, j, _, _, _ in prefix:
            link[seam] += j
        cuts[str(depth)] = dict(
            applied_twists=len(cut['twists']), gates=m.gates(fed, cut),
            twists_are_those_of_depth_at_most_the_cut=[t[:5] for t in cut['twists']] == [t[:5] for t in prefix],
            values_are_the_sum_over_those_twists=m.curl(fed, link) == cut['values'],
            defect_faces=sum(1 for k in cut['classes'] if k))
    # every applied twist: bounded current, pair sum, class rule
    fp, fm = fed.face_plus, fed.face_minus
    val = [0] * fed.faces
    births = moves = 0
    births_minimal = moves_kept = emptied_small = emptied_large = 0
    emptied_small_zero = emptied_large_nonzero = True
    for _, _, seam, j, _, a, b in twists:
        p, q = fp[seam], fm[seam]
        assert val[p] == a and val[q] == b
        na, nb = a + j, b - j
        if a == 0 and b == 0:
            births += 1
            births_minimal += na in m.MINIMAL[na % 6] and nb in m.MINIMAL[nb % 6] and na == -nb
        for before, after in ((a, na), (b, nb)):
            if before % 6 and after % 6 == 0:
                if abs(before) <= 2:
                    emptied_small += 1
                    emptied_small_zero &= after == 0
                elif abs(before) >= 4:
                    emptied_large += 1
                    emptied_large_nonzero &= after != 0
        if (a % 6 != 0) != (b % 6 != 0):
            src, dst = (a, nb) if a % 6 else (b, na)
            left = na if a % 6 else nb
            if left % 6 == 0:
                moves += 1
                moves_kept += src == dst
        val[p], val[q] = na, nb
    rng = random.Random(1729)
    link2 = [0] * fed.seams
    for _, _, seam, j, _, _, _ in twists:
        link2[seam] += j if rng.random() < 0.5 else j - 6 * (1 if j > 0 else -1)
    val2 = m.curl(fed, link2)
    null = dict(
        classes_equal=all(v % 6 == k for v, k in zip(val2, full['classes'])),
        carrier_totals_zero=all(sum(val2[20 * c:20 * c + 20]) == 0 for c in range(fed.carriers)),
        values_differ=val2 != full['values'],
        defect_faces_minimal=sum(1 for v in val2 if v % 6 and v in m.MINIMAL[v % 6]))
    digests = dict(loads=m.digest(full['loads']), register=m.digest(full['register']), classes=m.digest(full['classes']),
                   twists=m.twist_digest(twists))
    out = dict(
        level=fed.level, carriers=fed.carriers, seams=fed.seams, faces=fed.faces, sweeps=REFERENCE['sweeps'],
        applied_twists=len(twists), applied_by_abs_j=dict(sorted(Counter(str(abs(t[3])) for t in twists).items())),
        depth_at_end=[full['depth_min'], full['depth_max']], gates_at_end=end, sweeps_with_all_gates=sum(per_sweep),
        cuts=cuts, births=births, births_minimal=births_minimal, moves=moves, moves_with_value_kept=moves_kept,
        emptied_faces_of_size_at_most_2=emptied_small, emptied_faces_of_size_at_least_4=emptied_large,
        random_representative=null, digests=digests, census=m.census(fed, full))
    checks = dict(
        gates_at_the_end=all(end.values()),
        gates_after_every_sweep=len(per_sweep) == REFERENCE['sweeps'] and all(per_sweep),
        gates_on_every_cut=all(all(c['gates'].values()) for c in cuts.values()),
        cut_state_is_the_state_after_its_own_twists=all(
            c['twists_are_those_of_depth_at_most_the_cut'] and c['values_are_the_sum_over_those_twists'] for c in cuts.values()),
        current_at_most_three=full['max_abs_j'] <= 3,
        births_on_two_zero_faces_are_minimal=births > 0 and births == births_minimal,
        small_values_empty_to_zero=emptied_small > 0 and emptied_small_zero,
        large_values_leave_a_multiple_of_six=emptied_large_nonzero,
        random_representative_keeps_classes_and_stokes=null['classes_equal'] and null['carrier_totals_zero'] and null['values_differ'])
    return out, checks


def kinematics():
    cost, kappa = m.COST, m.KAPPA
    currents = Counter()
    for xa, xb, coin in product(range(6), range(6), (0, 1)):
        tot = xa + xb
        lo = tot // 2
        ni = tot - lo if coin else lo
        currents[(xa - xb, xa - ni)] += 1
    band = {(d, j) for (d, j) in currents if abs(d) <= 1}
    work = {str(abs(d)): sorted({j * (d - j) for (dd, j) in currents if dd == d and j}) for d in range(6)}
    births = sorted((d, j) for (d, j) in currents if j and 2 * cost[j % 6] <= kappa * j * (d - j))
    transport = {}
    for v in range(-12, 13):
        if v % 6:
            ends = sorted({v + j for j in range(-3, 4) if j and (v + j) % 6 == 0})
            transport[str(v)] = ends
    band_rule = {}
    for k1, k2 in product(range(6), repeat=2):
        acc = []
        for j in (1, -1):
            delta = j % 6
            n1, n2 = (k1 + delta) % 6, (k2 - delta) % 6
            if cost[n1] + cost[n2] - cost[k1] - cost[k2] <= 0:
                acc.append(dict(j=j, classes_after=[n1, n2]))
        band_rule[f'{k1},{k2}'] = acc
    checks = dict(
        K1_current_at_most_three=max(abs(j) for _, j in currents) == 3,
        K1_band_current_at_most_one=max(abs(j) for _, j in band) == 1,
        K1_released_potential_by_contrast=work == {'0': [], '1': [0], '2': [1], '3': [2], '4': [4], '5': [6]},
        K2_births_need_contrast_3_for_class_1_and_5=all(abs(d) == 3 for d, j in births if abs(j) == 1),
        K2_births_of_class_2_and_4_need_contrast_4_or_5=all(abs(d) in (4, 5) for d, j in births if abs(j) == 2),
        K2_births_of_class_3_need_contrast_5=all(abs(d) == 5 for d, j in births if abs(j) == 3),
        K2_birth_values_are_minimal=all(j in m.MINIMAL[j % 6] and -j in m.MINIMAL[-j % 6] for j in (1, 2, 3, -1, -2, -3)),
        K3_size_at_most_2_empties_to_zero=all(transport[str(v)] == [0] for v in (-2, -1, 1, 2)),
        K3_size_3_empties_to_zero_or_six=transport['3'] == [0, 6] and transport['-3'] == [-6, 0],
        K3_size_4_or_5_keeps_six=all(transport[str(v)] == [6] for v in (4, 5)) and all(transport[str(v)] == [-6] for v in (-4, -5)),
        K3_larger_values_keep_a_multiple_of_six=all(0 not in transport[str(v)] for v in range(-12, 13) if v % 6 and abs(v) >= 4),
        K4_band_class_1_moves_only=band_rule['1,0'] == [dict(j=-1, classes_after=[0, 1])]
        and band_rule['0,1'] == [dict(j=1, classes_after=[1, 0])],
        K4_band_class_5_moves_only=band_rule['5,0'] == [dict(j=1, classes_after=[0, 5])]
        and band_rule['0,5'] == [dict(j=-1, classes_after=[5, 0])],
        K5_band_classes_2_3_4_with_flat_neighbour_rest=all(band_rule[f'{k},0'] == [] and band_rule[f'0,{k}'] == [] for k in (2, 3, 4)),
        K5_band_no_birth=band_rule['0,0'] == [])
    return dict(currents_by_contrast={str(d): sorted({j for (dd, j) in currents if dd == d}) for d in range(-5, 6)},
                released_potential_by_contrast=work, births_accepted=[list(b) for b in births],
                value_after_emptying=transport, band_rule=band_rule), checks


def words_allowed(composition, values_of_class):
    """Words with every defect face at an allowed value of its class, total zero, all flat faces zero."""
    out = set()
    for combo in product(*[values_of_class[k] for k in composition]):
        if sum(combo) == 0:
            out.add(tuple(sorted(zip(composition, combo), key=lambda t: (t[0], -t[1]))))
    return sorted('(' + ','.join(m.fmt(v) for _, v in w) + ')' for w in out)


def readings():
    allowed = m.hypercharge_integers()
    by_class = {k: sorted((v for v in allowed if v % 6 == k), key=lambda v: (abs(v), -v)) for k in range(6)}
    table = {str(k): [dict(value=v, six_Y_of=NAMES[v]) for v in by_class[k]] for k in range(6)}
    words = {'{' + ','.join(map(str, c)) + '}': words_allowed(c, by_class) for c in COMPOSITIONS}
    minimal = {'{' + ','.join(map(str, c)) + '}': words_allowed(c, m.MINIMAL) for c in COMPOSITIONS}
    expected = {
        '{4,4,4}': ['(+4,-2,-2)'], '{2,2,2}': ['(+2,+2,-4)'], '{3,3}': ['(+3,-3)'], '{1,5}': ['(+1,-1)'],
        '{2,4}': ['(+2,-2)', '(-4,+4)'], '{1,2,3}': ['(+1,+2,-3)', '(+1,-4,+3)'], '{3,4,5}': ['(+3,-2,-1)', '(-3,+4,-1)'],
        '{2,2,4,4}': ['(+2,+2,-2,-2)', '(+2,-4,+4,-2)', '(-4,-4,+4,+4)'],
        '{2,3,3,4}': ['(+2,+3,-3,-2)', '(+2,-3,-3,+4)', '(-4,+3,+3,-2)', '(-4,+3,-3,+4)']}
    checks = dict(
        hypercharge_integers_of_the_32_patterns=allowed == [-6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6],
        only_5_is_excluded_up_to_6=sorted(set(range(-6, 7)) - set(allowed)) == [-5, 5],
        class_1_and_5_have_one_value=by_class[1] == [1] and by_class[5] == [-1],
        class_4_separates_u_R_from_d_R=by_class[4] == [-2, 4] and by_class[2] == [2, -4],
        left_handed_quarks_share_one_value=NAMES[1] == 'Q',
        words_with_total_zero_as_stated=all(sorted(words[c]) == sorted(expected[c]) for c in expected),
        minimal_values_cannot_fill_444_and_222=minimal['{4,4,4}'] == [] and minimal['{2,2,2}'] == [],
        minimal_values_fill_the_other_contents=all(minimal[c] for c in minimal if c not in ('{4,4,4}', '{2,2,2}')))
    return dict(hypercharge_integers=allowed, values_by_class=table, words_with_total_zero=words,
                minimal_words_with_total_zero=minimal,
                scope='the value separates u_R from d_R; u_L and d_L share +1; L and anti-H share -3; T3 is outside',
                statistics='every word is a pure flux content and a boson; a word names the flux seen by quanta of those hypercharges'), checks


INPUTS = ('model.py', 'census.py', 'certificate.py', 'federation.json', 'imported_census.json')


def input_hashes():
    return {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in INPUTS}


def certificate():
    struct, c1 = structure()
    ref, c2 = reference_run()
    kin, c3 = kinematics()
    read, c4 = readings()
    cen = census()
    c5 = dict(reference_digests_equal_the_public_engine=ref['digests'] == cen['reference_cross_check']['digests']
              and ref['applied_twists'] == cen['reference_cross_check']['applied_twists'])
    checks = {}
    for prefix, block in (('1_closed_carriers', c1), ('2_exact_gates', c2), ('3_kinematics', c3), ('4_readings', c4),
                          ('5_reference_run', c5), ('6_imported_census', cen['checks'])):
        for name, ok in block.items():
            checks[f'{prefix}.{name}'] = bool(ok)
    return dict(schema=SCHEMA, checks=checks, all_checks_pass=all(checks.values()), closed_carriers=struct,
                reference_run=ref, kinematics=kin, readings=read,
                imported_census=dict(source=cen['source'], reference_cross_check=cen['reference_cross_check'],
                                     derived=cen['derived'], pooled_mirror_pairs=cen['pooled_mirror_pairs'], runs=cen['runs']),
                input_sha256=input_hashes())


_RESULT = None


def result():
    global _RESULT
    if _RESULT is None:
        _RESULT = certificate()
    return _RESULT
