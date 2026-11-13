from fractions import Fraction as F
from itertools import product
import numpy as np
from scipy.linalg import eigh
from .source import run,coefficient,electrostatic
from model.openwave_light_nuclear_curve.source import projected_source,lift,UNIT
from model.openwave_light_nuclear_curve.certificate import INPUTS


def test_required_mirror_split_is_source_exact():
    r=run()['exact_isospin_obstruction']
    assert F(r['additional_charge_breaking_gap_needed_GeV'])==F('.00076380695')
    assert 2*F(r['minimum_charge_breaking_operator_norm_GeV'])==F('.00076380695')


def test_full_native_charge_negative_bound():
    g,_=coefficient()
    for row in product(range(4),repeat=6):
        N=sum(row)
        if N>9:continue
        charges=[F(2*row[2*x]-row[2*x+1],3)for x in range(3)]
        value=g*sum(charges[x]*charges[y]/(y-x)for x in range(3)for y in range(x+1,3))
        assert value>=-2*g*N*N/9


def test_independent_dense_charge_blocks_and_complete_records():
    g,_=coefficient();B=F(INPUTS['three_body_coefficient'])
    for L,Z in((3,1),(3,2),(4,2)):
        d=projected_source(L,B,True);Q=[lift(np.diag([1,0,1,0]),x,L)for x in range(L)]
        V=sum(Q[x]@Q[y]/(y-x)for x in range(L)for y in range(x+1,L))*float(g)
        ids=np.flatnonzero(np.diag(d['Iz'])==Z-L/2)
        e=eigh((d['H']+V)[ids][:,ids],eigvals_only=True)[0]*float(UNIT)
        r=electrostatic(L,Z);lo,hi=map(F,r['ground_GeV_interval'])
        assert float(lo)<=e<=float(hi)+1e-11
        assert F(r['complete_preparation_success_exact'])+F(r['complete_preparation_failure_exact'])==1
        assert r['complete_low_magnetic_dimension']==len(ids)


def test_failed_held_out_comparisons_are_retained():
    rows=run()['complete_finite_spectra']
    assert rows[0]['target_in_interval']
    assert not rows[1]['target_in_interval'] and not rows[2]['target_in_interval']
