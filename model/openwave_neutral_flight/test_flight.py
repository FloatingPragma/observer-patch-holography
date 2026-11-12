import numpy as np
import pytest
from .flight import data,experiment,receiver

def test_complete_native_aperture_and_momentum_reference():
 d=data();K=d['K'];F=d['F'];assert np.max(abs(K.conj().T@K+F.conj().T@F-np.eye(len(K))))<2e-14
 assert d['native_Ward_identities']==12
 with pytest.raises(ValueError):data(4)

def test_sharp_recoil_removes_spatial_interference():
 d=data();sharp=[experiment(d,1,3,x,'sharp')['total_received_probability']for x in range(5)];local=[experiment(d,1,3,x)['total_received_probability']for x in range(5)]
 assert np.ptp(sharp)<2e-15 and np.ptp(local)>1e-5

def test_every_failure_energy_and_mass_record_retained():
 d=data();r=experiment(d,1,7,2,'partial',True)
 assert len(r['complete_success_records_by_retained_mass_tag'])==3
 assert abs(r['complete_probability_sum']-1)<1e-13
 assert r['free_source_energy_error']<1e-14
 assert abs(r['source_preparation_mean_energy_work'])>1e-5
 assert abs(r['detector_mean_energy_work'])>1e-7

def test_no_port_phase_has_no_wrong_family():
 d=data(5,True)
 for a in range(3):
  p=experiment(d,a,9,2)['received_by_charged_family'];assert max(p[b]for b in range(3)if b!=a)<1e-25
