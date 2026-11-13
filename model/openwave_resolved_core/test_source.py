from fractions import Fraction as F
import json
import numpy as np
from scipy.linalg import eigh
import pytest
from . import source as s

@pytest.fixture(scope='module')
def member():return s.member(*s.charged_masses()[0])

def test_larger_box_independently_checks_ground_and_resolved_core(member):
    o=s.polynomial_member(28,s.CHARGE,member['mass'],regulator=s.REGULATOR)
    n=len(o['basis'])//2;e,u=eigh(o['H_minus_mass'][:n,:n])
    lo,hi=map(F,member['certificate']['full_ground_energy_minus_bare_mass_interval_exact'])
    assert float(lo)-1e-13<=e[0]<=float(hi)+1e-13
    x=np.array([s.COORDINATES[z[1]][0]for z in o['basis'][:n]])
    p=float(abs(u[:,0])**2@x)
    assert float(member['p_interval'][0])<=p<=float(member['p_interval'][1])
    assert 2*member['p_interval'][0]>F(1,2)

def test_complete_instrument_retains_failure_and_terminal_work(member):
    r=s.probe_records(member)
    assert F(r['suppression_at_three_quarter_pi_interval_exact'][1])<F(1,2)
    assert F(r['suppression_at_seven_eighth_pi_interval_exact'][0])>F(1,2)
    for x in r['records']:
        total=sum(sum(sum(z)for z in spin)for spin in x['complete_joint_position_flux_probabilities'])
        assert abs(total+.5-1)<1e-14
        final=x['final_mean_source_energy_with_failure_electron_units']
        work=r['original_preparation_mean_work_electron_units']+x['original_input_probe_mean_work_electron_units']
        assert abs(final-float(s.MAGNETIC)-work)<1e-13
        assert abs(x['original_input_probe_mean_work_electron_units']-x['original_input_phase_pointer_mean_work_electron_units']-x['original_input_terminal_position_flux_mean_work_electron_units'])<1e-13

def test_responses_enclose_independent_dense_spectral_sum(member):
    r=s.response_records(member)
    for row in r['certified_infinite_responses']:
        lo,hi=map(F,row['response_interval_exact'])
        assert float(lo)<=row['response_diagnostic']<=float(hi)
        assert F(row['actual_error_upper_exact'])<=s.RESPONSE_BUDGET

def test_changed_ground_witness_is_rejected(tmp_path,monkeypatch):
    r=json.loads(s.WITNESS.read_text());name=s.charged_masses()[0][0]
    r[name]['ground'][0][0]+=2**40
    p=tmp_path/'bad.json';p.write_text(json.dumps(r));monkeypatch.setattr(s,'WITNESS',p)
    with pytest.raises(AssertionError):s.member(*s.charged_masses()[0])


def test_static_charge_derivative_and_core_use_same_source(member):
    r=s.scalar_probe_coupling(member)
    for row in r['complete_interventions']:
        hf=list(map(F,row['Hellmann_Feynman_total_charge_response_interval_exact']))
        finite=list(map(F,row['complete_source_energy_central_difference_interval_exact']))
        error=F(row['finite_difference_to_derivative_error_upper_exact'])
        assert finite[0]-error<=hf[0]<=hf[1]<=finite[1]+error
        for key in ('minus_intervention','plus_intervention'):
            record=row[key]
            total=sum(F(x)for spin in record['both_spin_original_input_position_flux_probabilities_exact']for x in spin)
            assert total+F(record['preparation_failure_probability_exact'])==1
    assert not r['static_effective_charge_is_universal_renormalized_charge']


def test_opposite_charge_reflects_every_flux_and_probe_sign():
    r=s.opposite_charge_control(F(1))
    assert r['full_cap_dimension']==8*(2*s.CAP+1)
    assert r['opposite_charge_has_identical_sector_energy_and_radius']
    assert not r['whole_Standard_Model_charge_conjugation_symmetry_claimed']


def test_complete_vacuum_witness_replays_without_an_eigensolver(monkeypatch):
    from model.openwave_transverse_completion import source as parent
    r=s.frozen_vacuum()
    assert r['complete_cap_certificate']['dimension']==2*s.CAP+1
    assert len(r['complete_cap_certificate']['all_eigenvalue_centres_exact'])==2*s.CAP+1
    n=np.arange(-32,33,dtype=np.int64)
    A=np.diag(4*n*n+4)+np.diag(np.full(64,-2),1)+np.diag(np.full(64,-2),-1)
    larger_ground=eigh(A.astype(float)/32,eigvals_only=True)[0]
    lo,hi=map(F,r['full_infinite_rotor_vacuum_interval_exact'])
    assert float(lo)<larger_ground<float(hi)
    with monkeypatch.context() as changed:
        changed.setattr(s,'ELECTRIC',2*s.ELECTRIC)
        with pytest.raises(AssertionError):s.frozen_vacuum()
    def forbidden(*args,**kwargs):
        raise AssertionError('a replay must verify the fixed full frame')
    monkeypatch.setattr(parent,'eigh',forbidden)
    assert s.frozen_vacuum()==r


def test_corrupted_complete_vacuum_frame_is_rejected(tmp_path,monkeypatch):
    r=json.loads(s.VACUUM_WITNESS.read_text())
    for row in r['integer_frame']:row[0]=0
    p=tmp_path/'bad-vacuum.json';p.write_text(json.dumps(r));monkeypatch.setattr(s,'VACUUM_WITNESS',p)
    with pytest.raises(AssertionError):s.frozen_vacuum()
