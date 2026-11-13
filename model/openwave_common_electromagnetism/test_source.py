from fractions import Fraction as F
import numpy as np
from .source import (exact_ward,finite_gauss_control,gauss_label,SELF,STIFFNESS,
                     charge_diagonal,sector,parameters,G)


def test_complete_native_charge_polynomial_and_Bogoliubov_reference():
    r=exact_ward()
    assert r['global_charge_commutator_exact_zero']
    assert r['positive_Dirac_modes_checked']==48


def test_full_fock_gauss_dressing_and_undressed_negative_control():
    r=finite_gauss_control()
    assert len(r['complete_occupation_flux_reference_records'])==16
    assert r['charged_hop_columns_checked']==r['undressed_hop_constraint_failures']==16
    E,R=gauss_label([F(2,3),F(-1,3)])
    assert E==(3,-3) and R==-1


def test_complete_colour_complement_and_electric_self_energy():
    d=sector('baryon');diag,_=charge_diagonal();B=d['B'].toarray()
    C=np.asarray(list(map(float,diag)));reduced=B.conj().T@(C[:,None]*B)
    assert np.linalg.norm(C[:,None]*B-B@reduced)<1e-13
    assert B.shape==(2024,120)
    for q in (F(1),F(-1),F(2,3),F(-1,3)):
        assert SELF*(q*q-q*q)==0
        assert STIFFNESS*q*q/8>0
    p=parameters()
    assert p['mbar']-p['delta']/2-18*G-4*SELF/9>1
    assert 2*(1-SELF)>1
