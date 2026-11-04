import copy
from fractions import Fraction

import numpy as np
import mpmath as mp
import pytest
from scipy.linalg import expm

from . import construction as c
from . import certificate as cert
from . import packet


@pytest.fixture(scope='module')
def default():
    return c.design()


def test_exact_generic_algebra():
    assert c.exact_identities()['all_three_inverse_weight_identities_exact']


@pytest.mark.parametrize('branch',c.BRANCHES)
def test_admission_is_exact_and_exhaustive(branch):
    proof = cert.admissibility(branch)
    assert (proof['cubic_root_count'],proof['inverse_root_count']) == (3,2)
    assert proof['absolute_mass_bound']['operator_norm_mass_matrix_over_v_strict_upper_bound'] == '1'


@pytest.mark.parametrize('branch',c.BRANCHES)
@pytest.mark.parametrize('inverse',(0,1))
def test_full_record_realization(branch,inverse):
    row = c.design(branch,inverse)
    assert packet.independent_check(row)['all_checks_pass']
    law = c.direct_memory(row)
    assert law['retained_old_inputs'] == 12
    assert law['retained_records_per_input'] == 512
    assert np.count_nonzero(law['stationary'] > 1e-14) == 4


def test_source_law_actually_changes_with_pixel_branch():
    a,b = (c.design(branch) for branch in c.BRANCHES)
    assert a['p'] != b['p'] and a['q'] != b['q']
    assert not np.array_equal(a['mass_matrix_over_v'],b['mass_matrix_over_v'])


@pytest.mark.parametrize('key,value',[
    ('p',.2),('q',.3),('p',0),('p',float('nan')),
    ('gain_over_v',float('nan')),('gain_over_v',1),
    ('theta_left',0),('theta_right',0),('P','1'),
])
def test_independent_checker_rejects_changed_design(default,key,value):
    row = copy.deepcopy(default); row[key] = value
    with pytest.raises((ValueError,np.linalg.LinAlgError)):
        packet.independent_check(row)


def test_no_postselection_or_duty_renormalization(default):
    bad = copy.deepcopy(default)
    bad['coefficient_matrix'] /= bad['p']
    with pytest.raises(ValueError):
        packet.independent_check(bad)


def test_pixel_shape_enclosures_contain_recomputed_values():
    with mp.workdps(100):
        for branch in c.BRANCHES:
            masses = c.target(branch)['masses_over_v']
            total = sum(masses)
            for value,box in zip(masses,cert.normalized_mass_boxes(branch)):
                lo = mp.mpf(box.lo.numerator)/box.lo.denominator
                hi = mp.mpf(box.hi.numerator)/box.hi.denominator
                assert lo <= value/total <= hi


@pytest.mark.parametrize('steps',(8,32,128))
@pytest.mark.parametrize('stationary',(False,True))
def test_correlated_evolution_obeys_averaging_bound(default,steps,stationary):
    # Rescale the entire declared source gain to test a visible, dimensionless
    # finite-time evolution; the general theorem admits every nonnegative g.
    row = copy.deepcopy(default); row['gain_over_v'] = .4
    law = c.direct_memory(row)
    ket = np.array([1,2j,3,1j,0,2],complex); ket /= np.linalg.norm(ket)
    rho = np.outer(ket,ket.conj())
    initial = np.zeros((12,6,6),complex)
    if stationary:
        initial = law['stationary'][:,None,None]*rho
    else:
        initial[7] = rho
    actual = c.evolve(row,initial,steps).sum(axis=0)
    H = row['gain_over_v']*np.kron([[0,1],[1,0]],row['coefficient_matrix'])
    U = expm(-1j*H)
    error = np.linalg.norm(actual-U@rho@U.conj().T,ord='nuc')
    bound = c.averaging_bound(steps,p=row['p'],gain=row['gain_over_v'],stationary_product=stationary)
    assert error <= bound
    assert abs(np.trace(actual)-1) < 2e-12
    assert np.linalg.eigvalsh(actual).min() > -2e-12


@pytest.mark.parametrize('steps,time,p,gain',[(0,1,.1,1),(2,-1,.1,1),(2,1,1,1),(2,1,.1,-1)])
def test_bad_averaging_domain(steps,time,p,gain):
    with pytest.raises(ValueError):
        c.averaging_bound(steps,time,p=p,gain=gain)


def test_receipt_matches_fresh_scientific_calculation():
    value = packet.strict_json.load(packet.OUTPUT)
    packet.validate(value)


def test_receipt_scope_mutation_rejected(monkeypatch):
    value = packet.strict_json.load(packet.OUTPUT)
    changed = copy.deepcopy(value)
    changed['scope']['physical_interacting_rest_mass_or_clock_derived'] = True
    monkeypatch.setattr(packet,'report',lambda:value)
    with pytest.raises(ValueError):
        packet.validate(changed)


def test_receipt_pin_omission_rejected():
    value = packet.strict_json.load(packet.OUTPUT)
    del value['source_bindings'][c.MCPR]
    with pytest.raises(ValueError):
        packet.validate(value)
