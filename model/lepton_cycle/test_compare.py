"""Independent numerical and adversarial checks of the comparison contract."""
import copy
import json
from fractions import Fraction as F

import mpmath as mp
import numpy as np
import pytest

from . import compare


def test_frozen_shape_by_independent_double_precision_formula():
    result=compare.architecture()
    e=(1+np.sqrt(2)*np.cos(2/9+np.arange(3)*2*np.pi/3))**2
    assert float(result['muon_electron_ratio']) == pytest.approx(e[2]/e[1],rel=2e-13)
    assert float(result['tau_electron_ratio']) == pytest.approx(e[0]/e[1],rel=2e-13)
    assert result['physical_pole_mass_map_proved'] is False


@pytest.mark.parametrize('key', ['muon_electron_ratio','electron_MeV','tau_MeV'])
def test_targets_cannot_change_architecture_or_implementation(key):
    inputs=json.loads((compare.HERE/'experimental_inputs.json').read_text())
    before=compare.comparison(inputs)
    inputs[key]['value']=str(float(inputs[key]['value'])*1.01)
    after=compare.comparison(inputs)
    for field in ('architecture','implementation','preparation_sensitivity','exposure'):
        assert before[field]==after[field]
    assert before['muon']!=after['muon'] or before['tau']!=after['tau']


def test_reference_and_missing_covariance_are_not_theory_errors():
    r=compare.comparison()
    assert 442 < float(r['muon']['reference_only_standardized_discrepancy']) < 443
    assert 9.83 < float(r['muon']['relative_residual_ppm']) < 9.84
    assert 0.6 < float(r['tau']['reference_only_standardized_discrepancy_range'][0]) < 0.7
    assert r['uncertainty_contract']['joint_chi_squared'] is None
    assert r['uncertainty_contract']['physical_theory_standard_uncertainty'] is None
    assert not r['exposure']['phase_displacement_applied']


def test_covariance_extremes_and_boxes():
    r=compare.comparison(); inp=r['experimental_inputs']
    e=float(inp['electron_MeV']['value']); ue=float(inp['electron_MeV']['standard_uncertainty'])
    t=float(inp['tau_MeV']['value']); ut=float(inp['tau_MeV']['standard_uncertainty'])
    grad=np.array([1/e,-t/e**2])
    endpoints=[np.sqrt(grad@np.array([[ut**2,rho*ut*ue],[rho*ut*ue,ue**2]])@grad) for rho in (1,-1)]
    assert list(map(float,r['tau']['reference_ratio_linearized_sd_range'])) == pytest.approx(endpoints)
    lo,hi=map(float,r['tau']['conditional_one_standard_uncertainty_input_box_ratio'])
    for a in (t-ut,t+ut):
        for b in (e-ue,e+ue): assert lo-1e-11 <= a/b <= hi+1e-11


def test_negative_tau_residual_has_ordered_discrepancy_range():
    inputs=json.loads((compare.HERE/'experimental_inputs.json').read_text())
    inputs['tau_MeV']['value']='1800'
    lo,hi=map(float,compare.comparison(inputs)['tau']['reference_only_standardized_discrepancy_range'])
    assert lo < hi < 0


def test_no_mass_error_inferred_from_channel_error():
    r=compare.finite_design()
    for kind, coefficient in (('weak',F(14,81)),('mass',F(349,1296))):
        d=r['designs'][kind]; n=d['rounds']
        assert coefficient/n <= compare.ACCURACY/2 < coefficient/(n-1)
        assert F(d['certified_total_error']) <= compare.ACCURACY
    assert 'physical_mass_error' not in r


def test_effect_ratio_bound_and_denominator_rejection():
    with mp.workdps(50):
        eta=mp.mpf('1e-10'); row=compare.preparation_sensitivity(str(eta))
        e=list(map(mp.mpf,compare.architecture()['effect_eigenvalues_k0_k1_k2']))
        lo,hi=map(mp.mpf,row['muon_electron_rate_ratio_interval'])
        # Extremal permitted eigenvalue perturbations, independent of producer box implementation.
        assert abs(lo-(e[2]-eta)/(e[1]+eta)) < mp.mpf('1e-20')
        assert abs(hi-(e[2]+eta)/(e[1]-eta)) < mp.mpf('1e-20')
        assert lo < e[2]/e[1] < hi
    with pytest.raises(ValueError): compare.preparation_sensitivity('0.1')
    with pytest.raises(ValueError): compare.ratio_box(1,1,0,1)


@pytest.mark.parametrize('bad',['nan','inf','-1','0'])
def test_invalid_reference_rejected(bad):
    inputs=json.loads((compare.HERE/'experimental_inputs.json').read_text())
    inputs['tau_MeV']['value']=bad
    with pytest.raises(ValueError): compare.comparison(inputs)


def test_receipt_and_source_mutation_guard():
    assert json.loads((compare.HERE/'physical_comparison.json').read_text())==compare.report()
    bad=copy.deepcopy(compare.report()); bad['architecture']['theta']='0.222'
    assert bad!=compare.report()
