import numpy as np
import pytest
from scipy.linalg import expm
from . import spectrum as s


def test_exact_native_minimum_and_full_charge_classification():
    d=s.minimum_theorem()
    assert d['charged_colourless_Dirac_families']==3
    assert d['massive_charged_Dirac_channels_all_families']==21
    assert d['full_native_symbol']['complete_Nambu_dimension']==180


def test_complete_native_spectrum_and_unchanged_control():
    out=s.complete_census();control=s.complete_census(wilson=0)
    assert out['Nambu_rows']==11520
    assert out['minimum_count_per_charged_family']==[1,1,1]
    assert control['minimum_count_per_charged_family']==[8,8,8]
    assert len(out['complete_momentum_records'])==64
    assert all(len(r['all_180_Nambu_energies'])==180 for r in out['complete_momentum_records'])


def test_native_full_evolution_and_continuum_bound():
    # Compare actual full180 evolution with the independently assembled
    # source-defined four-coordinate frame, at off-grid momenta.
    frames=s.native_frames()['frames'];angles=np.array([.07,-.03,.05])
    for a in (1.,.5,.25):
        H,p,w=s.symbol(a*angles,a)
        for k in range(3):
            U=np.asarray(frames[k],complex);h=s.dirac_hamiltonian(p,s.masses()[k]+w)
            assert np.linalg.norm(expm(-.3j*H)@U-U@expm(-.3j*h))<2e-12
            continuum=s.dirac_hamiltonian(angles,s.masses()[k])
            bound=s.continuum_bound(s.masses()[k],a,float(np.linalg.norm(angles)),.3)
            assert np.linalg.norm(expm(-.3j*h)-expm(-.3j*continuum),2)<=bound['complete_record_TV_upper']+2e-14


def test_complete_arrivals_retain_original_probability():
    for k in (0,1,2):
        for t in (0.,.7):
            r=s.complete_arrival_records(k,t)
            assert np.shape(r['complete_256_arrival_probabilities'])==(64,4)
            assert abs(r['original_input_total']-1)<2e-12
            assert r['creation_failure_probability']==.5


def test_nonzero_Wilson_changes_only_off_rest_energies():
    for k,m in enumerate(s.masses()):
        assert np.linalg.eigvalsh(s.dirac_hamiltonian([0,0,0],m))[-1]==pytest.approx(m)
        # A former corner is a strict maximum along its pi coordinate.
        x=np.array([np.pi,0,0]);d=.001
        energy=lambda y: np.sqrt(np.sin(y)@np.sin(y)+(m+np.sum(1-np.cos(y)))**2)
        assert energy(x)>energy(x+[d,0,0]) and energy(x)>energy(x-[d,0,0])
