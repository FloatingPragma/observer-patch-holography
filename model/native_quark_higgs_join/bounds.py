"""Finite native quark/Higgs interference on the full interacting source."""
from fractions import Fraction as F
from functools import lru_cache
import sympy as s

from model.native_quark_response.response import positive_response
from model.native_quark_response.full import windows, down, encoded
from model.native_charge_resolution.bounds import data as native_data, frac, root, round_up
from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.resolvent import (
    field_data, polynomial, encoded as field_encoded,
)
from model.coherent_defect_wave.chiral_joint.family_cell.quotient import evaluate


def interval_mul(a, b):
    values = [x*y for x in a for y in b]
    return min(values), max(values)


def interval_square(a):
    return (F(0) if a[0] <= 0 <= a[1] else min(a[0]**2, a[1]**2),
            max(a[0]**2, a[1]**2))


def field_interval(value, ground):
    return tuple(frac(v) for v in evaluate(polynomial(value), ground))


@lru_cache(None)
def atomic(family):
    data = positive_response(family)
    _, _, field = field_data()
    e = field.unit
    m0, m1, m2 = data['moments'][:3]
    cross = m1-e*m0
    cnorm = m2-2*e*m1+e*e*m0
    ground = data['source_ground_interval']
    totals = {sign: (m0+cnorm/4+sign*cross)/4 for sign in (1, -1)}
    assert m0 == field.convert(s.Rational(1, 2))
    records = []
    for row in data['records']:
        d = tuple(frac(x) for x in row['difference_interval'])
        w = tuple(frac(x) for x in row['weight_interval'])
        fringe = {}
        for sign in (1, -1):
            factor = tuple(sorted((1+sign*d[0]/2, 1+sign*d[1]/2)))
            fringe[sign] = tuple(v/4 for v in interval_mul(w, interval_square(factor)))
        records.append(dict(energy=tuple(frac(x) for x in row['energy_interval']),
                            excitation=d, elementary=w, cross=interval_mul(w, d),
                            composite=interval_mul(w, interval_square(d)), fringe=fringe,
                            exactly_zero_residue=row['exactly_zero_residue']))
    for sign in (1, -1):
        interval = field_interval(totals[sign], ground)
        assert 0 < interval[0] <= interval[1] < 1
        assert sum(row['fringe'][sign][0] for row in records) <= interval[1]
        assert sum(row['fringe'][sign][1] for row in records) >= interval[0]
    return dict(records=records, totals=totals, cross=cross, composite_norm_squared=cnorm,
                source_ground=ground)


@lru_cache(None)
def response(family):
    source = native_data()
    inherited = windows(family)
    selected = inherited['records'][0]
    local = atomic(family)
    v = source['v']
    eta = F(inherited['exact_ground_trial_vector_error'])
    delta = F(inherited['actual_preparation_trace_distance_upper'])
    qlo = F(inherited['preparation_success_lower'])
    wl, wu = map(F, selected['exact_ground_input_probability'])
    lower_energy = F(selected['source_excitation_relative_actual_ground'][0])
    gamma = down(lower_energy*wl-2*v)
    assert gamma > 0
    window = tuple(map(F, selected['absolute_active_source_plus_reference_bin']))
    elo, ehi = source['e']
    atom_energies = [(row['energy'][0]+7*elo+4, row['energy'][1]+7*ehi+4)
                     for row in local['records']]
    inner = min(atom_energies[0][0]-window[0], window[1]-atom_energies[0][1])
    assert inner > v
    settings = []
    for sign in (1, -1):
        weight = local['records'][0]['fringe'][sign]
        leakage = F(0)
        for energy, row in zip(atom_energies[1:], local['records'][1:], strict=True):
            distance = max(window[0]-energy[1], energy[0]-window[1])
            assert distance > 0
            leakage += round_up(v*root(row['fringe'][sign][1], True)/distance)
        lo = max(F(0), root(weight[0]*(1-(v/inner)**2))-leakage-eta)**2
        hi = min(F(1), (root(weight[1], True)+leakage+eta)**2)
        hit = down(lo), round_up(hi)
        total = field_interval(local['totals'][sign], local['source_ground'])
        total = max(F(0), total[0]-eta), min(F(1), total[1]+eta)
        outside = max(F(0), total[0]-hit[1]), min(F(1), total[1]-hit[0])
        failure = 1-total[1], 1-total[0]
        prepared = {key: (max(F(0), val[0]-delta), min(F(1), val[1]+delta))
                    for key, val in [('pump_failure', failure), ('bin_hit', hit), ('outside_bin', outside)]}
        assert all(0 <= interval[0] <= interval[1] <= 1 for interval in prepared.values())
        assert sum(interval[0] for interval in prepared.values()) <= 1 <= sum(interval[1] for interval in prepared.values())
        settings.append(dict(sign=sign, pump='(A+sign*C/2)/2',
                             exact_atomic_pump_probability_coefficients=field_encoded(local['totals'][sign]),
                             exact_ground_hit_probability=encoded(hit),
                             exact_ground_total_pump_probability=encoded(total),
                             prepared_conditional_records={key: encoded(val) for key, val in prepared.items()},
                             original_attempt_joint_setting_records={key: encoded((qlo*val[0]/2, val[1]/2))
                                                                     for key, val in prepared.items()},
                             diagnostic_joint_setting_record_lowers={key: float(qlo*val[0]/2)
                                                                     for key, val in prepared.items()},
                             failed_pump_energy_hit_probability='0',
                             other_atomic_lines_leakage_amplitude_upper=str(round_up(leakage))))
    contrast = down(qlo*(gamma/4-delta))
    assert contrast > 0
    return dict(family=family, native_source= 'h=H0+2^-34 T on the complete eight-cell carrier',
                exact_ground_trial_vector_error=str(eta), actual_preparation_trace_distance_upper=str(delta),
                preparation_success_lower=str(qlo), full_hopping_norm=str(v),
                absolute_active_source_plus_reference_bin=list(selected['absolute_active_source_plus_reference_bin']),
                source_excitation_relative_actual_ground=list(selected['source_excitation_relative_actual_ground']),
                reference_success_cost='4', reference_failure_cost='0',
                exact_ground_elementary_bin_probability=list(selected['exact_ground_input_probability']),
                exact_ground_elementary_composite_cross_lower=str(gamma),
                exact_ground_normalized_composite_bin_probability_lower=str(down(gamma*gamma/(4*wu))),
                conditional_hit_contrast_lower=str(down(gamma/2)),
                fair_setting_joint_contrast_lower=str(down(gamma/4)),
                original_attempt_joint_contrast_lower=str(contrast),
                diagnostic_cross_lower=float(gamma), diagnostic_original_attempt_joint_contrast_lower=float(contrast),
                exact_atomic_cross_total_coefficients=field_encoded(local['cross']),
                exact_atomic_composite_norm_squared_coefficients=field_encoded(local['composite_norm_squared']),
                all_atomic_records=[dict(index=i, source_excitation=encoded(row['excitation']),
                                         elementary_weight=encoded(row['elementary']),
                                         elementary_composite_cross=encoded(row['cross']),
                                         composite_weight=encoded(row['composite']),
                                         fringe_weights={str(sign): encoded(row['fringe'][sign]) for sign in (1, -1)},
                                         exactly_zero_residue=row['exactly_zero_residue'])
                                    for i, row in enumerate(local['records'])],
                settings=settings,
                complete_successful_preparation_labels=8,
                coarse_preparation_failure_probability=['0', str(1-qlo)],
                complete_law=['1-q, with every physical preparation failure record retained',
                              'q/2 Tr(K_sign rho K_sign^dagger), energy pointer zero',
                              'q/2 Tr(P M_sign rho M_sign^dagger), energy pointer one',
                              'q/2 Tr((I-P) M_sign rho M_sign^dagger), energy pointer zero',
                              'zero for failed pump and energy pointer one, for each sign'],
                positive_cross_selector='an actual bin eigenspace maximizing positive real elementary/composite cross weight',
                selector_equals_maximum_elementary_residue_selector=False)


def calculate():
    result = {family: response(family) for family in ('uc', 'dc', 'Q')}
    assert atomic('uc')['cross'] == atomic('dc')['cross']
    assert atomic('uc')['composite_norm_squared'] == atomic('dc')['composite_norm_squared']
    return result


if __name__ == '__main__':
    for name, row in calculate().items():
        print(name, row['diagnostic_cross_lower'], row['diagnostic_original_attempt_joint_contrast_lower'], flush=True)
        for setting in row['settings']:
            print(setting['sign'], setting['diagnostic_joint_setting_record_lowers'], flush=True)
