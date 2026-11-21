"""Imported run counts of the level-10 universe run of the public simulator.

These numbers are copied from the simulator's published evidence files; this package does not recompute them.
It checks only arithmetic on them: class neutrality of each carrier content and the difference of the
{4,4,4} and {2,2,2} counts in Poisson standard errors.
"""
import math

SOURCE = dict(
    repository='https://github.com/muellerberndt/oph-physics-sim',
    commit='476aff7219fe05592563cb0a1281aba4cffc9d96',
    run='level-10 universe run (L10, 300 sweeps, preparation exact06)',
    files={
        'oph_defects/viewer2/run_evidence.json': '74ac209556c652229e6139421e5042aba45883833204ec633b2f9b11847d2a21',
        'oph_defects/particle_map.json': 'deada9a482abde138e4ceec7aba6099fb68603656c57cad071250e9d4e6d083e',
    },
    status='imported run counts, not recomputed by this package')

# oph_defects/particle_map.json, defect_map.groups (carriers holding each content)
CARRIER_CONTENTS = {
    '{2,4}': 4991967, '{1,5}': 1041883, '{3,3}': 785555, '{2,2,4,4}': 139090,
    '{2,2,2}': 86788, '{4,4,4}': 86476, '{2,3,3,4}': 55690, '{1,2,3} and {3,4,5}': 22351}
CARRIERS_TOTAL = 20971520

# oph_defects/viewer2/run_evidence.json, runs[L10].screen
PAIRS = {
    '(1,5)': dict(first_cut_pair_births=124786, whole_run_pair_births=757258, whole_run_annihilations=609494),
    '(2,4)': dict(first_cut_pair_births=208860, whole_run_pair_births=649302, whole_run_annihilations=3585),
    '(3,3)': dict(first_cut_pair_births=41328, whole_run_pair_births=89909, whole_run_annihilations=17)}
IDENTITIES = {'1': 3053295, '2': 2486092, '3': 1607810, '4': 2481326, '5': 3050633}
ALIVE_AT_END = {'1': 3402, '2': 163775, '3': 163721, '4': 163082, '5': 3351}


def contents(label):
    return [tuple(int(x) for x in part.strip(' {}').split(',')) for part in label.split(' and ')]


def census():
    neutral = {label: all(sum(c) % 6 == 0 for c in contents(label)) for label in CARRIER_CONTENTS}

    def sigma(a, b):
        return round((a - b) / math.sqrt(a + b), 2)
    return dict(
        source=SOURCE, carriers_total=CARRIERS_TOTAL, carrier_contents=CARRIER_CONTENTS,
        carrier_contents_class_neutral=neutral, pairs=PAIRS, identities_ever=IDENTITIES,
        alive_at_end=ALIVE_AT_END,
        matter_antimatter_in_standard_errors={
            '{4,4,4} against {2,2,2}': sigma(CARRIER_CONTENTS['{4,4,4}'], CARRIER_CONTENTS['{2,2,2}'])})
