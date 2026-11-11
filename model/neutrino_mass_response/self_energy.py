"""Exact normal Feshbach selfenergy of the complete native neutral measure."""
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sympy as s

from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import moments
from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.resolvent import (
    encoded, field_data, field_value, polynomial,
)
from model.coherent_defect_wave.chiral_joint.family_cell.quotient import evaluate, outward
from model.defect_scalar_response.response import ground_pair

ROOT = Path(__file__).resolve().parents[2]
PARENT = 'model/coherent_defect_wave/chiral_joint/ground_composite_response/results.json'


def value_record(value):
    interval = evaluate(polynomial(value), ground_pair()['ground_interval'])
    return {'coefficients_descending_in_source_ground_energy': encoded(value),
            'exact_outward_interval': outward(interval, 60),
            'numeric_midpoint': float(sum(interval)/2)}


@lru_cache(maxsize=1)
def calculate():
    parent = json.loads((ROOT/PARENT).read_text())
    for name, expected in parent['input_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != expected:
            raise ValueError('the actual neutral-ground source changed: '+name)
    complete = parent['actual_ground_composite_spectrum']
    reference = next(row for row in complete['recorded_settings'] if row['kind'] == 'remove' and row['probe_spin'] == 0)
    data, p, q, cross = moments('remove', 0)
    assert [[str(x) for x in row] for row in cross] == reference['exact_cross_moments_target_rows_source_columns']
    assert str(p.as_expr()) == reference['exact_source_minimal_polynomial']
    assert str(q.as_expr()) == reference['exact_target_minimal_polynomial']
    _, _, field = field_data()
    x, z = p.gen, s.Symbol('z')
    coefficients = []
    for i in range(p.degree()):
        coefficient = s.Poly(sum(p.nth(j)*x**(j-i-1) for j in range(i+1, p.degree()+1)), x)
        coefficients.append(field_value(coefficient)/field_value(p.diff()))
    power, rank = s.eye(data['source_dimension']), field.zero
    for coefficient in coefficients:
        rank += coefficient*field.convert(s.trace(power))
        power = data['source_H']*power
    assert rank == field.convert(2)
    raw = [sum((coefficient*field.convert(trace) for coefficient, trace in zip(coefficients, row, strict=True)),
               field.zero)/2 for row in cross]
    entries = {}
    for i, moment in enumerate(raw):
        for j in range(i+1, q.degree()+1):
            entries[j-i-1] = entries.get(j-i-1, field.zero)+field.convert(q.nth(j))*moment
    numerator = s.Poly.from_dict({(i,): value for i, value in entries.items()}, z, domain=field)
    denominator = s.Poly.from_dict({(i,): field.convert(q.nth(i)) for i in range(q.degree()+1)}, z, domain=field)
    assert s.gcd(numerator, denominator).degree() == 0
    assert numerator.degree() == 8 and denominator.degree() == 9
    mean = raw[1]/raw[0]
    variance = raw[2]/raw[0]-mean*mean
    quotient, remainder = denominator.mul_ground(raw[0]).div(numerator)
    assert quotient.rep.nth(1) == field.one and quotient.rep.nth(0) == -mean
    assert remainder.degree() == 7
    assert -remainder.rep.nth(7)/numerator.rep.nth(8) == variance
    assert evaluate(polynomial(raw[0]), ground_pair()['ground_interval'])[0] > 0
    assert evaluate(polynomial(variance), ground_pair()['ground_interval'])[0] > 0
    return {'source': 'actual all-three-unit local30CAR/5scalar physical-ground reference',
        'fresh_complete_target_dimension': data['target_dimension'],
        'fresh_complete_source_spin_dimension': data['source_dimension'],
        'source_ground_projector_rank': 2,
        'all_source_target_CAR_metric_maps_and_cross_moments_regenerated': True,
        'complete_original_addition_removal_records': complete,
        'unconditioned_J_norm_squared': value_record(raw[0]),
        'joint_success_probability_for_one_of_four_settings': value_record(raw[0]/16),
        'normalized_first_energy_transfer_moment': value_record(mean-field.unit),
        'normalized_second_central_energy_moment': value_record(variance),
        'complete_cyclic_dimension_of_the_actual_measure': 9,
        'cyclic_Feshbach_complement_dimension': 8,
        'full_target_complement_and_degeneracies_discarded': False,
        'raw_resolvent': 'R(w)=N(w)/D(w); its original total weight is retained separately',
        'normalized_gap_resolvent': 'G(z)=N(z+e)/(mu0 D(z+e))',
        'normal_selfenergy': 'Sigma(z)=-(division_remainder)(z+e)/N(z+e)',
        'Feshbach_identity': 'G(z)^(-1)=z-(mu1/mu0-e)-Sigma(z)',
        'selfenergy_denominator_degree': 8,
        'selfenergy_poles': 'eight real simple roots, one strictly between each adjacent pair of the nine visible neutral energy lines',
        'selfenergy_residues': 'nonnegative, from P H Q(z-QHQ)^(-1)Q H P in the positive-metric cyclic realization',
        'selfenergy_high_frequency_coefficient': value_record(variance),
        'exact_N_coefficients_descending_in_w': [encoded(numerator.rep.nth(i)) for i in range(numerator.degree(), -1, -1)],
        'exact_D': str(denominator.as_expr()),
        'exact_selfenergy_numerator_coefficients_descending_in_w': [encoded(-remainder.rep.nth(i)) for i in range(remainder.degree(), -1, -1)],
        'Majorana_or_DeltaL_two_selfenergy_generated': False,
        'extra_coupling_or_heavy_mass_parameter_supplied': False,
        'pole_count_is_generation_count': False,
        'physical_neutrino_rest_mass_or_GeV_scale': False}
