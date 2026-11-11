"""The actual full-thirty-mode physical ground and complete scalar images."""
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as s

from model.defect_scalar_response.response import ground_pair, scalar_expectation
from model.coherent_defect_wave.chiral_joint.local_reference import physical
from model.coherent_defect_wave.chiral_joint.local_reference.reduction import (
    coordinate_kernel, lepton_operators, quark_module,
)
from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.model import weight_basis
from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.resolvent import encoded, polynomial
from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import reduced, validate_annihilator
from model.coherent_defect_wave.chiral_joint.family_cell.quotient import (
    add, divide, evaluate, mul, outward, sqrt_interval,
)
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory

ROOT = Path(__file__).resolve().parents[2]
QUARK_RECEIPT = 'model/native_quark_response/results.json'
LEPTON_RECEIPT = 'model/coherent_defect_wave/chiral_joint/interacting_charged_response/results.json'


def expectation(operator):
    return scalar_expectation(operator)['value']


def encoded_value(value):
    interval = evaluate(polynomial(value), ground_pair()['ground_interval'])
    return {'coefficients_descending_in_e': encoded(value),
            'interval': outward(interval, 60)}


@lru_cache(maxsize=1)
def scalar_population():
    source = ground_pair()['block']
    raw_metric = physical.metric_on(source['keys'])
    populations = []
    for level in range(5):
        raw = s.diag(*(int(key[2] == level) for key in source['keys']))
        compressed = source['gram'].inv()*source['basis'].H*raw_metric*raw*source['basis']
        populations.append(expectation(compressed))
    field = ground_pair()['field']
    assert sum(populations, field.zero) == field.one
    delta = ground_pair()['delta']
    n = (field.one-field.one/delta)/2
    assert populations == [field.one-n, n/4, n/4, n/4, n/4]
    e = expectation(source['H'])
    assert e == field.unit
    covariance = field.one-3*n/4
    normalization = 2-7*n/4-3*e/4
    return dict(populations=populations, scalar_occupation=n,
                scalar_covariance=covariance, scalar_normalization=normalization,
                interaction_energy=e-n)


@lru_cache(maxsize=None)
def scalar_target(hyper_sign):
    """Complete spin/weak highest-weight scalar insertion target."""
    if hyper_sign not in (-1, 1):
        raise ValueError('scalar hypercharge orientation required')
    keys = weight_basis(3, 1, 1, hyper_sign)
    q, lepton, scalar = quark_module(), lepton_operators(), higgs_memory()[1]
    weak = physical.tensor_map(keys, weight_basis(3, 1, 3, hyper_sign),
        quark=q['currents'][0]+s.I*q['currents'][1],
        lepton=lepton['currents'][0]+s.I*lepton['currents'][1],
        scalar=scalar[8]+s.I*scalar[9])
    spin = physical.tensor_map(keys, weight_basis(3, 3, 1, hyper_sign),
        lepton=lepton['spin'][0]+s.I*lepton['spin'][1])
    basis, _ = coordinate_kernel(weak.col_join(spin))
    expected = (len(keys)-len(weight_basis(3, 1, 3, hyper_sign))
                -len(weight_basis(3, 3, 1, hyper_sign))
                +len(weight_basis(3, 3, 3, hyper_sign)))
    assert basis.cols == expected
    return reduced(keys, basis)


@lru_cache(maxsize=None)
def scalar_measure(hyper_sign):
    source = ground_pair()['block']
    target = scalar_target(hyper_sign)
    fields = higgs_memory()[0]
    operator = fields[0].H if hyper_sign == 1 else fields[1]
    raw = physical.tensor_map(source['keys'], target['keys'], scalar=operator)
    image = raw*source['basis']
    matrix = target['left']*image.extract(target['coordinates'], range(source['basis'].cols))
    assert target['basis']*matrix == image
    lower_keys = weight_basis(3, 1, -1, hyper_sign)
    q, lepton, currents = quark_module(), lepton_operators(), higgs_memory()[1]
    lowering = physical.tensor_map(target['keys'], lower_keys,
        quark=q['currents'][0]-s.I*q['currents'][1],
        lepton=lepton['currents'][0]-s.I*lepton['currents'][1],
        scalar=currents[8]-s.I*currents[9])
    lower_basis = lowering*target['basis']
    assert lower_basis.H*physical.metric_on(lower_keys)*lower_basis == target['gram']
    other_field = fields[1].H if hyper_sign == 1 else -fields[0]
    other_image = physical.tensor_map(source['keys'], lower_keys, scalar=other_field)*source['basis']
    assert lowering*image == other_image
    assert physical.hamiltonian_on(lower_keys)*lower_basis == lower_basis*target['H']
    adjoint = source['gram'].inv()*matrix.H*target['gram']
    x = s.Symbol('x')
    annihilator = s.Poly(target['H'].charpoly(x).as_expr(), x).sqf_part()
    validate_annihilator(target['H'], annihilator)
    zeroth = expectation(adjoint*matrix)
    first = expectation(adjoint*target['H']*matrix)-ground_pair()['field'].unit*zeroth
    assert zeroth == scalar_population()['scalar_covariance']
    return dict(target=target, map=matrix, adjoint=adjoint,
                zeroth=zeroth, first=first, annihilator=annihilator)


@lru_cache(maxsize=1)
def scalar_readout():
    positive, negative = scalar_measure(1), scalar_measure(-1)
    population = scalar_population()
    assert positive['zeroth'] == negative['zeroth']
    assert positive['first'] == negative['first']
    assert positive['first']+negative['first'] == population['scalar_normalization']
    field = ground_pair()['field']
    c = population['scalar_covariance']
    no_operation = 3*population['scalar_occupation']/4
    assert c+no_operation == field.one
    phi, _, number, _ = higgs_memory()
    effects = sum((p.H*p+p*p.H for p in phi), s.zeros(5))/4
    assert effects+3*number/4 == s.eye(5)
    data = ground_pair()
    assert data['n']**2 == data['n']
    assert data['n']*data['v']+data['v']*data['n'] == data['v']
    canonical = sqrt_interval(evaluate(polynomial(population['scalar_normalization']),
                                       data['ground_interval']))
    return {'scalar_population': [encoded_value(x) for x in population['populations']],
            'scalar_covariance_each_component': encoded_value(c),
            'scalar_retarded_first_moment_each_component': encoded_value(population['scalar_normalization']),
            'canonical_contact_modulus_interval': outward(canonical, 60),
            'scalar_retarded_inverse_z_term': '0',
            'complete_scalar_orientation_targets': [{
                'hypercharge_orientation': sign,
                'dimension': scalar_measure(sign)['target']['H'].rows,
                'full_weight_dimension': len(scalar_measure(sign)['target']['keys']),
                'full_target_annihilator': str(scalar_measure(sign)['annihilator'].as_expr()),
                'zeroth_spectral_moment': encoded_value(scalar_measure(sign)['zeroth']),
                'first_energy_transfer_moment': encoded_value(scalar_measure(sign)['first']),
                'complete_image_reconstruction': True} for sign in (1, -1)],
            'other_doublet_component_metric_and_H_intertwiner_checked': True,
            'complete_scalar_branch_probabilities': {
                'each_of_four_component_orientation_records': encoded_value(c/4),
                'no_operation': encoded_value(no_operation)},
            'literal_five_level_Kraus_completeness': True,
            'all_energy_projectors_included_in_moment_sum': True,
            'complete_energy_records_per_component_orientation': scalar_spectrum(),
            'no_operation_complete_support': [
                {'energy': 'e', 'probability': encoded_value(3*population['scalar_occupation']**2/4)},
                {'energy': '1-e', 'probability': encoded_value(3*population['scalar_occupation']*(1-population['scalar_occupation'])/4)},
                {'energy': 'all other energies', 'probability': encoded_value(field.zero)}],
            'individual_energy_residues_tabulated': True}


@lru_cache(maxsize=1)
def scalar_spectrum():
    data = scalar_measure(1)
    opposite = scalar_measure(-1)
    field = ground_pair()['field']
    ground = ground_pair()['ground_interval']
    q = data['annihilator']
    image = data['map']
    opposite_image = opposite['map']
    moments = []
    for _ in range(q.degree()):
        moment = expectation(data['adjoint']*image)
        assert moment == expectation(opposite['adjoint']*opposite_image)
        moments.append(moment)
        image = data['target']['H']*image
        opposite_image = opposite['target']['H']*opposite_image
    z = s.Symbol('z')
    entries = {}
    for i, moment in enumerate(moments):
        for j in range(i+1, q.degree()+1):
            entries[j-i-1] = entries.get(j-i-1, field.zero)+field.convert(q.nth(j))*moment
    numerator = s.Poly.from_dict({(i,): value for i, value in entries.items()}, z, domain=field)
    denominator = s.Poly(q.as_expr().subs(q.gen, z), z, domain=field)
    assert q == opposite['annihilator']
    assert numerator.rep.nth(q.degree()-1) == data['zeroth']
    common = s.gcd(numerator, denominator)
    visible_n = numerator.exquo(common)
    visible_d = s.Poly(denominator.exquo(common).as_expr(), z, domain=s.QQ)
    invisible = s.Poly(common.as_expr(), z, domain=s.QQ)
    records = []
    for bounds, multiplicity in q.intervals(eps=s.Rational(1, 10**75)):
        assert multiplicity == 1
        bounds = tuple(map(s.Rational, bounds))
        zero = (invisible.eval(bounds[0]) == 0 if bounds[0] == bounds[1]
                else invisible.count_roots(*bounds) == 1)
        if zero:
            weight = (s.S.Zero, s.S.Zero)
        else:
            value = (s.S.Zero, s.S.Zero)
            for power in range(visible_n.degree(), -1, -1):
                value = add(mul(value, bounds), evaluate(polynomial(visible_n.rep.nth(power)), ground))
            weight = divide(value, evaluate(visible_d.diff(), bounds))
            assert weight[0] > 0
        records.append({'target_energy_interval': outward(bounds, 60),
                        'source_energy_transfer_interval': outward((bounds[0]-ground[1], bounds[1]-ground[0]), 60),
                        'probe_probability_interval': outward((weight[0]/4, weight[1]/4), 70),
                        'exact_zero_residue': zero})
    return {'records': records, 'full_target_root_count': q.degree(),
            'exact_resolvent_numerator_coefficients_descending_in_z':
                [encoded(numerator.rep.nth(i)) for i in range(numerator.degree(), -1, -1)],
            'exact_resolvent_denominator': str(denominator.as_expr()),
            'exact_ground_moments_coefficients_descending_in_e': [encoded(value) for value in moments],
            'zero_residue_factor': str(invisible.as_expr()),
            'complete_original_weight': encoded_value(data['zeroth']/4),
            'coefficient_scaling_of_each_actual_phi_probe': '1/2',
            'same_complete_measure_for_both_doublet_components_and_orientations': True,
            'full_source_charged_gaps_assumed_positive': False}


@lru_cache(maxsize=1)
def removal_moments():
    field = ground_pair()['field']
    result = {}
    for path, sector, key, divisor in (
        (QUARK_RECEIPT, 'exact_local_response', 'exact_ground_moments', 1),
        (LEPTON_RECEIPT, 'exact_local_channels', 'positive_total_raw_moments_coefficients_in_E0', 2)):
        receipt = json.loads((ROOT/path).read_text())
        for name, digest in receipt['input_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
        for family, row in receipt[sector].items():
            moments = []
            for coefficients in row[key][:2]:
                value = field.zero
                for coefficient in coefficients:
                    value = value*field.unit+field.convert(s.Rational(coefficient))
                moments.append(value/divisor)
            assert moments[0] == field.convert(s.Rational(1, 2))
            result[family] = moments[1]-field.unit*moments[0]
    return result
