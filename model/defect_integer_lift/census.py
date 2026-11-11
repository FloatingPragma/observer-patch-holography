"""Imported run counts of the universe runs of the public simulator, levels 4 to 7.

The counts in imported_census.json were produced outside this package: the twist list of each run was regenerated
with the public engine, the integer lift was built from that list, and the classes on every stored cut and the
identities were compared with the stored run files. This module checks arithmetic on the counts only: every listed
word has total zero, every value sits in its class, shares and mirror differences.
"""
import json
import math
from pathlib import Path

from . import model as m

CENSUS = json.loads(Path(__file__).with_name('imported_census.json').read_text())


def parse_word(word):
    """'(+4,-2,-2|+6)' -> ([4, -2, -2], [6])."""
    body = word.strip('()')
    left, _, right = body.partition('|')
    values = [int(v) for v in left.split(',') if v]
    hidden = [int(v) for v in right.split(',') if v]
    return values, hidden


def composition(word):
    values, _ = parse_word(word)
    return '{' + ','.join(str(k) for k in sorted(v % 6 for v in values)) + '}'


def shares(run):
    """Shares at the end state of one run, as exact fractions written num/den, and as floats."""
    allowed = set(m.hypercharge_integers())
    out = {}
    for tier in ('internal', 'screen'):
        table = run['end'][tier]['values_by_class']
        defects = {int(v): n for k in '12345' for v, n in table[k].items()}
        hidden = {int(v): n for v, n in table['0'].items()}
        nd = sum(defects.values())
        row = dict(
            defect_faces=nd, hidden_faces=sum(hidden.values()),
            at_a_hypercharge_integer=sum(n for v, n in defects.items() if v in allowed),
            minimal=sum(n for v, n in defects.items() if v in m.MINIMAL[v % 6]),
            abs_4=sum(n for v, n in defects.items() if abs(v) == 4),
            abs_5=sum(n for v, n in defects.items() if abs(v) == 5),
            abs_7_or_more=sum(n for v, n in defects.items() if abs(v) >= 7),
            hidden_abs_6=sum(n for v, n in hidden.items() if abs(v) == 6),
            hidden_abs_12_or_more=sum(n for v, n in hidden.items() if abs(v) >= 12))
        row['share_at_a_hypercharge_integer'] = row['at_a_hypercharge_integer'] / nd
        row['share_minimal'] = row['minimal'] / nd
        row['hidden_per_flat_face'] = row['hidden_faces'] / run['end'][tier]['flat_faces']
        c4, c2 = table['4'], table['2']
        row['class_4_share_plus_4'] = c4.get('+4', 0) / sum(c4.values())
        row['class_4_share_minus_2'] = c4.get('-2', 0) / sum(c4.values())
        row['class_2_share_minus_4'] = c2.get('-4', 0) / sum(c2.values())
        row['class_2_share_plus_2'] = c2.get('+2', 0) / sum(c2.values())
        row['null_random_representative'] = run['end'][tier]['null_N2']
        mirror = {}
        for name, (ka, va, kb, vb) in {
                '4:+4 against 2:-4': ('4', '+4', '2', '-4'), '4:-2 against 2:+2': ('4', '-2', '2', '+2'),
                '3:+3 against 3:-3': ('3', '+3', '3', '-3'), '0:+6 against 0:-6': ('0', '+6', '0', '-6'),
                '4:-8 against 2:+8': ('4', '-8', '2', '+8')}.items():
            a, b = table[ka].get(va, 0), table[kb].get(vb, 0)
            mirror[name] = dict(a=a, b=b, standard_errors=round((a - b) / math.sqrt(a + b), 2) if a + b else None)
        row['mirror_pairs'] = mirror
        out[tier] = row
    words = {}
    for comp, strict, minimal in (('{4,4,4}', '(+4,-2,-2)', '(-2,-2,-2|+6)'), ('{2,2,2}', '(+2,+2,-4)', '(+2,+2,+2|-6)')):
        if comp in run['words_end']:
            w = run['words_end'][comp]
            a, b, n = w['words'].get(strict, 0), w['words'].get(minimal, 0), w['carriers']
            words[comp] = dict(carriers=n, word_with_one_value_of_size_4=strict, count=a, share=a / n,
                               binomial_error=math.sqrt(a / n * (1 - a / n) / n),
                               minimal_word_with_a_flat_6=minimal, minimal_count=b, minimal_share=b / n, other=n - a - b)
    out['words'] = words
    return out


def census():
    runs = CENSUS['runs']
    allowed = set(m.hypercharge_integers())
    checks = {}
    checks['every_gate_passes_at_every_level'] = all(
        r['all_gates_pass'] and all(g['passed'] for g in r['gates'].values()) for r in runs.values())
    checks['no_class_mismatch_and_no_nonzero_total'] = all(
        g.get('mismatches', 0) == 0 and g.get('nonzero_totals', 0) == 0 for r in runs.values() for g in r['gates'].values())
    checks['every_value_sits_in_its_class'] = all(
        int(v) % 6 == int(k) for r in runs.values() for tier in ('internal', 'screen')
        for k, row in r['end'][tier]['values_by_class'].items() for v in row)
    checks['every_listed_word_has_total_zero'] = all(
        sum(parse_word(w)[0]) + sum(parse_word(w)[1]) == 0
        for r in runs.values() for row in r['words_end'].values() for w in row['words'])
    checks['every_listed_word_has_its_composition'] = all(
        composition(w) == comp and all(h % 6 == 0 and h for h in parse_word(w)[1])
        for r in runs.values() for comp, row in r['words_end'].items() for w in row['words'])
    checks['internal_values_sum_to_zero'] = all(
        sum(int(v) * n for row in r['end']['internal']['values_by_class'].values() for v, n in row.items()) == 0
        and sum(int(v) * n for row in r['end']['screen']['values_by_class'].values() for v, n in row.items()) == 0
        for r in runs.values())
    checks['twists_keep_the_pair_sum_and_events_are_counted'] = all(
        sum(e['count'] for e in r['events'].values()) == r['applied_twists'] for r in runs.values())
    derived = {name: shares(r) for name, r in runs.items()}
    checks['run_is_more_concentrated_than_the_random_representative'] = all(
        d['internal']['share_at_a_hypercharge_integer'] > d['internal']['null_random_representative']['share_hypercharge_integers'] + 0.3
        and d['internal']['share_minimal'] > d['internal']['null_random_representative']['share_minimal'] + 0.5
        for d in derived.values())
    checks['free_current_sum_matches_classes_at_chance'] = all(
        abs(r['null_N1_end']['internal']['share_mod_6_equal_to_stored_class'] - 1 / 6) < 0.005 for r in runs.values())
    checks['first_cut_is_minimal_at_every_level'] = all(r['time']['first_cut_all_minimal'] for r in runs.values())
    checks['hypercharge_reading_is_not_exact'] = all(
        d['internal']['abs_7_or_more'] > 0 for d in derived.values())
    checks['minimal_reading_is_not_exact'] = all(d['internal']['abs_4'] > 0 for d in derived.values())
    checks['flat_faces_hold_multiples_of_six'] = all(d['internal']['hidden_faces'] > 0 for d in derived.values())
    pooled = {}
    for d in derived.values():
        for name, p in d['internal']['mirror_pairs'].items():
            a, b = pooled.get(name, (0, 0))
            pooled[name] = (a + p['a'], b + p['b'])
    mirror = {name: dict(a=a, b=b, standard_errors=round((a - b) / math.sqrt(a + b), 2)) for name, (a, b) in pooled.items()}
    checks['pooled_mirror_pairs_within_three_standard_errors'] = all(abs(p['standard_errors']) < 3.0 for p in mirror.values())
    checks['imported_counts_are_labelled'] = CENSUS['source']['status'].startswith('imported run counts')
    checks['reference_equals_public_engine'] = all(CENSUS['reference_cross_check']['public_engine_equal'].values())
    return dict(source=CENSUS['source'], reference_cross_check=CENSUS['reference_cross_check'], runs=runs,
                derived=derived, pooled_mirror_pairs=mirror, checks=checks)
