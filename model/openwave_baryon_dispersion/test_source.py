from fractions import Fraction as F
import numpy as np
import sympy as sy
import pytest
from scipy.linalg import eigh
from .source import matrix,repair_matrix,phase_record,repaired
from .native import certificate


def test_full_native_frame_and_channel_symmetries():
    r=certificate()
    assert r['complete_light_triquark_shape_channels']==100
    assert r['spin_and_isospin_commutators_exact_zero']
    assert r['actual_neutron_five_state_intertwiner']
    assert r['complete_rest_matrix_entries_from_exact_CAR_enumeration']
    assert F(r['complete_rest_matrix_positive_floor_source_exact'])>0


def test_inverse_rule_is_positive_and_rest_cancellation_exact():
    for proton in(False,True):
        A,G=matrix(F(1),F(0),proton);J,norm,zeta=repair_matrix(proton)
        assert A*G.inv()*J==sy.Rational(zeta)*G
        assert np.linalg.eigvalsh(np.asarray(J,float))[0]>0
        assert phase_record(F(1),F(0),True,proton)==phase_record(F(1),F(0),False,proton)
    with pytest.raises(ValueError):matrix(F(1),F(1))


def test_realification_retains_every_complex_compact_band():
    A,G=matrix(F(3,5),F(4,5));J,_,_=repair_matrix();J=sy.diag(J,J)
    source=A+sy.Rational(F(4,5))*J
    numeric=eigh(np.asarray(source,float),np.asarray(G,float),eigvals_only=True)
    row=phase_record(F(3,5),F(4,5),True)
    for value,pair in zip(numeric[::2],row['all_five_compact_source_intervals']):
        lo,hi=map(F,pair);assert float(lo)-1e-9<=value<=float(hi)+1e-9
    assert np.max(abs(numeric[::2]-numeric[1::2]))<1e-9


def test_constructive_kinetic_repair_retains_all_complements():
    r=repaired();lo,hi=map(F,r['repaired_inertial_mass_GeV_interval'])
    assert lo<F('.93956542194')<hi
    assert F(r['all_momentum_P_upper'])<F(r['full_complement_floor'])
    assert F(r['available_ground_bandwidth_GeV_interval'][0])>F('.14')
    assert not r['other_full_Fock_spectra_transferred']
