from fractions import Fraction as F
from pathlib import Path
import numpy as np
import sympy as s
from scipy.linalg import expm
import pytest
HERE=Path(__file__).resolve().parent
from . import locality as l


def test_actual_full_edge_CAR_current_and_charge_identities():
    dim=256;c=[]
    for mode in range(8):
        entries={}
        for word in range(dim):
            if word&(1<<mode):
                entries[word^(1<<mode),word]=(-1)**((word&((1<<mode)-1)).bit_count())
        c.append(s.SparseMatrix(dim,dim,entries))
    alpha,beta=l.matrices();I=s.eye(dim)
    tail=sum((x.H*x for x in c[:4]),s.zeros(dim))-2*I
    head=sum((x.H*x for x in c[4:]),s.zeros(dim))-2*I
    for a in alpha:
        B=-(beta+s.I*a)/2
        Fop=sum((B.H[i,j]*c[4+i].H*c[j] for i in range(4) for j in range(4)),s.zeros(dim))
        J=s.I*(Fop-Fop.H)
        assert tail*Fop-Fop*tail==-Fop
        assert head*Fop-Fop*head==Fop
        h=s.zeros(8);h[:4,4:]=-s.I*B;h[4:,:4]=s.I*B.H
        for k in range(8):
            assert c[k]*J-J*c[k]==sum((h[k,j]*c[j] for j in range(8)),s.zeros(dim))
        assert J**5-5*J**3+4*J==s.zeros(dim)


def test_rational_majorant_encloses_independent_matrix_exponential():
    c=l.field.cube();initial=np.zeros(32)
    initial[12:24]=np.abs(np.array(c['C'][3,:],float).ravel())/64
    M=np.array(l.comparison_matrix(),float)
    for tau in (F(1,64),F(1,16),F(1,4),F(1)):
        out=l.commutator_bound(tau)
        actual=(expm(float(tau)*M)@initial)[12]/64
        lo=float(F(out['retained_nonnegative_series_lower']))
        hi=float(F(out['commutator_norm_upper']))
        assert lo-2e-18<=actual<=hi+2e-18
    coeff=l.commutator_bound(F(1,64))['first_nine_unnormalized_Taylor_coefficients']
    assert coeff[:4]==['0','0','1/4096','0']


def test_physical_gauss_fiber_and_current_geometry():
    assert l.source_algebra()['all_three_Wilson_hopping_spin_matrices_are_rank_two_partial_isometries']
    c=l.field.cube();D,C,G,PT=(c[k] for k in ('D','C','G','PT'))
    L=s.Matrix.hstack(*D.columnspace())
    assert D*G*D.T*L==L
    assert C*L==s.zeros(C.rows,L.cols)
    assert D.T*PT==s.zeros(8,12)
    assert C[3,0]==0


def test_complete_prepared_influence_and_all_outcomes():
    d=l.scientific();b=d['complete_influence_exact']
    assert F(b['complete_preparation_flags_and_Bob_bit_TV_upper'])<F(182,10**10)
    assert F(b['conditional_prepared_code_complement_weight_upper'])<F(6,10**14)
    assert d['same_finite_H_for_preparation_kick_wait_and_detection']
    assert d['initial_state_comparison_error']=='0'
    assert d['photon_cap']==2**64
    assert d['apparatus_resources']['number_of_duration_one_driven_pulses']==6
    assert 'all M^3' in d['complete_failures']
    assert not d['spatial_or_spacing_refinement']


def test_positive_path_comparison_and_parameter_domain():
    free=l.comparison_matrix(F(0));coupled=l.comparison_matrix(F(1,64))
    assert all(a<=b for ra,rb in zip(free,coupled) for a,b in zip(ra,rb))
    assert l.commutator_bound(F(0))['commutator_norm_upper']=='0'
    with pytest.raises(ValueError):l.commutator_bound(F(-1))
