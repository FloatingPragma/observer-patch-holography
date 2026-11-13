import numpy as np
from .source import source_packet,dispersion,boost_checks


def test_literal_native_frame_charge_and_complete_probe():
    d=source_packet();T=np.array(d['native_positive_frame_real'])+1j*np.array(d['native_positive_frame_imag'])
    assert T.shape==(96,48) and np.linalg.norm(T.conj().T@T-np.eye(48))<3e-15
    q=np.asarray(d['charges'])
    assert np.linalg.norm(np.r_[q,-q][:,None]*T-T*q[None,:])==0


def test_full_brillouin_speed_and_continuum_error_without_mode_cut():
    rng=np.random.default_rng(88)
    for a in (.125,.5,1.):
        p=rng.uniform(-np.pi/a,np.pi/a,(500,3));m=.7
        E=dispersion(m,p,a);continuum=np.sqrt(m*m+np.sum(p*p,axis=1))
        v=np.sin(a*p)/(a*E[:,None])
        assert np.max(np.linalg.norm(v,axis=1))<=1
        error=a*a*np.sum(p*p,axis=1)**2/(24*m)
        assert np.all(E<=continuum+1e-13) and np.all(continuum-E<=error+1e-13)
        assert np.all(E>m)
        assert dispersion(m,np.zeros(3),a)==m


def test_boost_defect_quadratic_refinement_and_retained_source_frame():
    rows=boost_checks()
    assert all(r['full_regulator_frame_retained'] for r in rows)
    assert all(rows[i+1]['boost_energy_shell_defect_electron_units_diagnostic']<.26*rows[i]['boost_energy_shell_defect_electron_units_diagnostic'] for i in range(3))
    assert all(rows[i+1]['boost_energy_shell_defect_upper_diagnostic']==rows[i]['boost_energy_shell_defect_upper_diagnostic']/4 for i in range(3))
