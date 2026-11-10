"""Exact ground-selected positive quark measures without numerical eigenvectors."""
from functools import lru_cache

import sympy as s

from model.coherent_defect_wave.chiral_joint.atomic_reference_control.ground_exchange import actual_spin_coefficients
from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.resolvent import (
    encoded, field_data, field_value, polynomial,
)
from model.coherent_defect_wave.chiral_joint.family_cell.quotient import add, divide, evaluate, mul, outward
from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import source, validate_annihilator

from .local import CHANNELS, transfer


@lru_cache(maxsize=None)
def positive_response(family):
    _, _, field = field_data()
    sb = source(1)
    h = sb['H']
    x = s.Symbol('x')
    p = s.Poly(h.charpoly(x).as_expr(), x).sqf_part()
    block = transfer(family)
    characteristic = s.Poly(block['H'].charpoly(x).as_expr(), x)
    q = characteristic.sqf_part()
    validate_annihilator(block['H'], q)
    validate_annihilator(h, p)
    coefficients = []
    powers = [s.eye(h.rows)]
    for i in range(p.degree()):
        c = s.Poly(sum(p.nth(j)*x**(j-i-1) for j in range(i+1, p.degree()+1)), x)
        coefficients.append(field_value(c)/field_value(p.diff()))
        if i+1 < p.degree():
            powers.append(h*powers[-1])
    trace = sum((c*field.convert(s.trace(power)) for c, power in zip(coefficients, powers)), field.zero)
    assert trace == field.one
    image = block['map']
    moments = []
    for _ in range(q.degree()):
        pullback = block['adjoint']*image
        moments.append(sum((c*field.convert(s.trace(power*pullback))
                            for c, power in zip(coefficients, powers)), field.zero))
        image = block['H']*image
    assert moments[0] == field.convert(s.Rational(1, 2))
    entries = {}
    for i, moment in enumerate(moments):
        for j in range(i+1, q.degree()+1):
            power = j-i-1
            entries[power] = entries.get(power, field.zero)+field.convert(q.nth(j))*moment
    z = s.Symbol('z')
    numerator = s.Poly.from_dict({(k,): value for k, value in entries.items()}, z, domain=field)
    denominator = s.Poly.from_dict({(k,): field.convert(q.nth(k)) for k in range(q.degree()+1)}, z, domain=field)
    common = s.gcd(numerator, denominator)
    reduced_numerator = numerator.exquo(common)
    reduced_denominator = denominator.exquo(common)
    # Invisible factors are determined in the source number field first.
    # Their rationality is an executed property, not a factor-selection rule.
    invisible = s.Poly(common.as_expr(), z, domain=s.QQ)
    support = s.Poly(reduced_denominator.as_expr(), z, domain=s.QQ)
    ground = tuple(map(s.Rational, actual_spin_coefficients()['ground_energy_isolating_interval']))
    records = []
    for interval, multiplicity in q.intervals(eps=s.Rational(1, 10**70)):
        assert multiplicity == 1
        interval = tuple(map(s.Rational, interval))
        zero = (invisible.eval(interval[0]) == 0 if interval[0] == interval[1]
                else invisible.count_roots(*interval) == 1)
        value = (s.Integer(0), s.Integer(0))
        if zero:
            weight = (s.Integer(0), s.Integer(0))
        else:
            for k in range(reduced_numerator.degree(), -1, -1):
                coefficient = evaluate(polynomial(reduced_numerator.rep.nth(k)), ground)
                value = add(mul(value, interval), coefficient)
            weight = divide(value, evaluate(support.diff(), interval))
            assert 0 < weight[0] <= weight[1] < s.Rational(1, 2)
        records.append(dict(energy_interval=interval, weight_interval=weight,
                            difference_interval=(interval[0]-ground[1], interval[1]-ground[0]),
                            exactly_zero_residue=zero))
    assert len(records) == q.degree()
    assert sum(r['weight_interval'][0] for r in records) < s.Rational(1, 2)
    assert sum(r['weight_interval'][1] for r in records) > s.Rational(1, 2)
    return dict(q=q, numerator=numerator, denominator=denominator, moments=tuple(moments),
                source_ground_interval=ground, records=records, block=block,
                invisible_factor=invisible, visible_support=support,
                full_characteristic=characteristic)


def certificate(family):
    data = positive_response(family)
    block = data['block']
    return dict(family=family, native_removed_mode=block['native_removed_mode'],
                full_weight_dimension=block['full_weight_dimension'],
                complete_target_dimension=block['H'].rows,
                complete_target_characteristic_factorization=str(s.factor(data['full_characteristic'].as_expr())),
                complete_target_annihilator=str(data['q'].as_expr()),
                local_source_ground_energy_interval=outward(data['source_ground_interval'], 60),
                exact_zero_residue_factor=str(data['invisible_factor'].as_expr()),
                visible_support_annihilator=str(data['visible_support'].as_expr()),
                exact_ground_moments=[encoded(value) for value in data['moments']],
                complete_positive_records=[dict(
                    source_energy_interval=outward(r['energy_interval'], 60),
                    source_excitation_interval=outward(r['difference_interval'], 60),
                    probability_interval=outward(r['weight_interval'], 80),
                    exactly_zero_residue=r['exactly_zero_residue']) for r in data['records']],
                pump_success='1/2', pump_failure='1/2', outside_support='0',
                reference_cost_for_success='4', reference_cost_for_failure='0')


if __name__ == '__main__':
    for family in CHANNELS:
        data = positive_response(family)
        print(family, data['block']['H'].rows, data['q'].degree(),
              'visible', data['visible_support'].degree(),
              'zero_factor', data['invisible_factor'].as_expr(), flush=True)
        for row in data['records']:
            print(float(sum(row['difference_interval'])/2),
                  float(sum(row['weight_interval'])/2), flush=True)
