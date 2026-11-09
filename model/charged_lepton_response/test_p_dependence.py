from copy import deepcopy
from fractions import Fraction as F
import json

import mpmath as mp
import pytest

from model.native_flavor_orbit.intervals import I
from model.mcpr_native_memory_bridge import construction
from . import p_dependence as pd


@pytest.fixture(scope='module')
def retained():
    return pd.certificate()


def test_complete_authenticated_branch_set_and_scope(retained):
    assert [r['branch'] for r in retained['branches']] == [
        'forward', 'measured_endpoint', 'mssm_beta4_proxy', 'sm_color_beta7']
    assert retained['authenticated_source_pin_counts'] == [241, 237]
    assert all(r['admitted_complete_memory_inverse_branches'] == 2
               for r in retained['branches'])
    assert retained['physical_interacting_mass_predictions'] == []
    assert retained['cold_root_or_native_memory_producers_replayed_here'] is False
    assert retained['optional_MeV_unit_calibration_is_source_derived'] is False
    for row in retained['branches']:
        a, b = map(F, row['optional_electron_unit_mass_intervals_MeV'][2])
        lo, hi = map(F, row['source']['charged_action_ratios'][2])
        assert a <= lo*F('0.51099895069') <= hi*F('0.51099895069') <= b
    json.dumps(retained, allow_nan=False)


@pytest.mark.parametrize('branch', construction.BRANCHES)
def test_source_formulas_against_independent_high_precision_circulant(branch):
    # Independent diagonalization of the three circulant Fourier eigenvalues.
    # This uses no root solver, native memory calculation or empirical reference.
    p, u = construction.source(branch)
    value = pd.evaluate(p, u)
    with mp.workdps(80):
        p, u = mp.mpf(p), mp.mpf(u)
        chi = u*u/(512*(1-u/77))
        delta = mp.mpf(2)/9+(p-(1+mp.sqrt(5))/2)**2/(21*mp.pi*(1-u/27-u*u/135))
        rho = mp.sqrt(2)*mp.exp(-chi)
        roots = sorted(1+rho*mp.cos(delta+2*mp.pi*j/3) for j in range(3))
        ratios = [(v/roots[0])**2 for v in roots]
        for exact, box in zip(ratios, value['charged_action_ratios']):
            lo = mp.mpf(box.lo.numerator)/box.lo.denominator
            hi = mp.mpf(box.hi.numerator)/box.hi.denominator
            assert lo <= exact <= hi
        q = sum(v*v for v in roots)/sum(roots)**2
        box = value['positive_root_Koide_ratio']
        assert mp.mpf(box.lo.numerator)/box.lo.denominator <= q <= mp.mpf(box.hi.numerator)/box.hi.denominator


def test_both_closure_root_mass_boxes_are_reused_without_resolving(retained):
    for row in retained['branches'][2:]:
        assert 'closure_alpha_inverse' in row and 'v_over_declared_E_star' in row
        assert 'unique root within' in row['ancestry']
    outcomes = [(r['comparison']['muon_electron_ratio']['source_interval_inside_reference_one_sigma'],
                 r['comparison']['tau_electron_ratio']['source_interval_inside_reference_one_sigma'])
                for r in retained['branches']]
    assert outcomes == [(False, True), (True, True), (False, True), (False, False)]


def test_reference_mutation_cannot_change_source_coordinates_or_ratios():
    value = pd.evaluate(*construction.source('forward'))
    before = pd._json(value)
    reference = pd.strict_json.load(pd.ROOT/pd.REFERENCE)
    changed = deepcopy(reference)
    changed['tau_electron_ratio']['value'] = '4000'
    assert pd.comparison(value, reference) != pd.comparison(value, changed)
    assert pd._json(value) == before
    assert pd._json(pd.evaluate(*construction.source('forward'))) == before


def test_declared_law_has_capacity_dependence_at_fixed_coupling():
    first, second = pd.evaluate('1.63', '.04'), pd.evaluate('1.631', '.04')
    assert first['rho'] == second['rho']
    assert first['delta'].hi < second['delta'].lo
    assert first['charged_action_ratios'][1].hi < second['charged_action_ratios'][1].lo


def test_forged_receipt_is_not_admitted(tmp_path):
    path = tmp_path/'forged.json'
    actual = pd.strict_json.load(pd.bridge_owner.OUTPUT)
    actual['source_bindings'].pop(next(iter(actual['source_bindings'])))
    path.write_text(json.dumps(actual))
    with pytest.raises(ValueError, match='incomplete'):
        pd.authenticate(path, pd.bridge_owner.source_bindings())


@pytest.mark.parametrize('mutation', ['branch', 'count', 'records', 'domain'])
def test_inverse_branch_scope_is_consumed_from_owner(mutation):
    row = pd.strict_json.load(pd.closure_owner.OUTPUT)['closures'][0]['memory_admissibility']
    assert pd._inverse_count(row) == 2
    if mutation == 'branch': row['inverse_branches'].pop()
    elif mutation == 'count': row['quadratic_root_count'] = 1
    elif mutation == 'records': row['all_complete_record_probabilities_retained'] = False
    else: row['root_domain'] = 'center only'
    with pytest.raises(ValueError, match='root-box'):
        pd._inverse_count(row)


@pytest.mark.parametrize('p,u', [(True,'.04'),(1.63,'.04'),('1','.04'),('2','.04'),('1.63','0'),('1.63','1')])
def test_invalid_source_domain_rejected(p,u):
    with pytest.raises(ValueError): pd.evaluate(p,u)


def test_input_box_not_mutated():
    p = I(F(163,100)); u = I(F(1,25))
    value = pd.evaluate(p,u)
    assert value['P'] == p and value['u'] == u
    assert value['P'] is not p and value['P'].lo is not p.lo
