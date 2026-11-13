from fractions import Fraction as F
import numpy as np
import pytest
from scipy.linalg import eigh
from model.openwave_uniform_pole_bands import source as s


def test_entire_circle_cover_and_positive_complete_gaps():
    result=s.run()
    assert result['circle_cover']['every_real_momentum_covered']
    for member in result['members']:
        assert len(member['records'])==64
        assert F(member['all_momentum_gap_lower_source'])>F(127,1000)
        assert F(member['all_momentum_high_probability_upper'])<F(27,10000000)
        for sign in(1,-1):
            rows=[r for r in member['records']if r['half_angle_chart_sign']==sign]
            assert rows[0]['half_angle_tangent_interval'][0]=='0'
            assert rows[-1]['half_angle_tangent_interval'][1]=='1'
            assert all(a['half_angle_tangent_interval'][1]==b['half_angle_tangent_interval'][0]for a,b in zip(rows,rows[1:]))


def test_nonmidpoint_exact_source_spectra_inside_uniform_cells():
    for proton in(False,True):
        for index,sign in((0,1),(12,1),(0,-1),(31,-1)):
            row=s.cell(index,sign,proton)
            for offset in(F(0),F(2,7),F(1)):
                c,z=s.halfangle((index+offset)/s.CELLS,sign);A,G=s.exact_symbol(c,z,proton)
                eig=eigh(np.array(A,complex),np.array(G,float),eigvals_only=True)
                interval=list(map(F,row['all_three_full_source_GeV_intervals'][0]));unit=s.member(proton)['unit']
                assert float(interval[0])<=eig[0]*float(unit)<=float(interval[1])
                assert eig[1]-eig[0]>=float(F(row['complete_ground_isolation_gap_lower_source']))


def test_light_high_frobenius_bound_and_invalid_chart_rejection():
    for proton in(False,True):
        d=s.member(proton)
        for t in(F(0),F(1,3),F(1)):
            c,z=s.halfangle(t,-1);A,G=s.exact_symbol(c,z,proton)
            low=[float(G[i,i])**-.5 for i in s.LOW];high=[float(G[i,i])**-.5 for i in s.HIGH]
            block=np.diag(low)@np.array(A.extract(s.LOW,s.HIGH),complex)@np.diag(high)
            assert np.linalg.norm(block,2)**2<=float(d['coupling_square'])
    with pytest.raises(ValueError):s.halfangle(F(2))
    with pytest.raises(ValueError):s.halfangle(F(0),0)
