from copy import deepcopy
import pytest
from model.weyl_geometric_modular import packet


@pytest.fixture(scope='module')
def receipt():
    return packet.build()


def test_current_source_and_explicit_dependency_custody(receipt):
    packet.validate(receipt)
    assert receipt['new_checks_pass'] is True
    assert receipt['all_pass'] == receipt['dependency_custody_pass']
    assert receipt['custody']['prior_full_source_replay'] is False
    assert receipt['custody']['historical_theorem_basis']['source_mismatches_at_that_commit'] == 0
    prior = packet.strict_json.load(packet.ROOT / packet.PRIOR)
    assert receipt['custody']['prior_binding_count'] == len(prior['source_bindings'])


@pytest.mark.parametrize('flag', [
    'dynamical_Einstein_equation_derived', 'full_G1_through_G6_closed',
    'edge_area_density_or_Newton_constant_derived', 'interacting_or_massive_sector_supplied',
    'absolute_type_III_entropy_assumed_finite', 'common_lattice_collar_chart_limit_proved',
])
def test_reject_missing_premise_promotions(receipt, flag):
    bad = deepcopy(receipt); bad['theorem'][flag] = True
    with pytest.raises(ArithmeticError): packet.validate(bad)


def test_reject_deleted_valley_wrong_clock_or_spinor_correction(receipt):
    for path, value in [
        (('algebra', 'source', 'Weyl_fields'), 7),
        (('algebra', 'required_gradient_coefficient'), '0'),
        (('geometry', 'diamond', 'modular_energy_in_model_clock'), '2*pi*integral_B w_R*energy_density'),
    ]:
        bad = deepcopy(receipt); node = bad['certificate']
        for key in path[:-1]: node = node[key]
        node[path[-1]] = value
        with pytest.raises(ArithmeticError): packet.validate(bad)


def test_reject_silent_custody_relabeling_and_source_mutation(receipt):
    bad = deepcopy(receipt); bad['custody']['prior_full_source_replay'] = True
    with pytest.raises(ArithmeticError): packet.validate(bad)
    bad = deepcopy(receipt); bad['source_bindings'][next(iter(bad['source_bindings']))] = '0'*64
    with pytest.raises(ArithmeticError): packet.validate(bad)
    bad = deepcopy(receipt); bad['new_checks_pass'] = 1
    with pytest.raises(ArithmeticError): packet.validate(bad)


@pytest.mark.parametrize('location', ['theorem', 'certificate', 'custody'])
def test_scientific_and_custody_booleans_do_not_alias_integer_zero(receipt, location):
    bad = deepcopy(receipt)
    if location == 'theorem':
        bad['theorem']['dynamical_Einstein_equation_derived'] = 0
    elif location == 'certificate':
        bad['certificate']['algebra']['stress']['Einstein_equation_used'] = 0
    else:
        bad['custody']['prior_full_source_replay'] = 0
    with pytest.raises(ArithmeticError): packet.validate(bad)
