"""Independent controls for routing, complete events and spectral bounds."""
from fractions import Fraction as F
import numpy as np
import pytest
from . import source as s


def test_cubic_routing_complete_native_columns_and_padding():
    for L in (2,3,4):
        g=s.graph(L)
        assert(g.cells,g.ports,len(g.seams),len(g.inter_seams))==(L**3,12*L**3,36*L**3,6*L**3)
        assert all(len(set(v))==6 for v in g.neighbours)
        c=s.native_columns(g);assert c['complete_population_matrix_entries']==g.ports**2
        u=s.complete_event_control(g);assert u['all_source_event_basis_columns_checked']==g.ports*g.events
        assert u['additional_padded_W_minus_one_modes']==g.ports*u['event_padding_dimension']
        assert u['full_padded_source_event_dimension']==g.ports*g.events+u['additional_padded_W_minus_one_modes']
        z=s.walk.full_arc_spectrum(g)
        assert z['unprepared_complement_plus_one_multiplicity']>0
        assert sum(z[k]for k in('full_W_plus_one_multiplicity','full_W_minus_one_multiplicity','dispersive_dimension'))==g.ports*g.events
    with pytest.raises(ValueError):s.graph(1)


def test_exact_hessian_and_fixed_anisotropic_control():
    j=s.exact_hessian();assert j['quadratic_coefficient']==[['1/8','0','0'],['0','1/8','0'],['0','0','1/8']]
    c=s.exact_hessian(True);assert c['quadratic_coefficient']==[['4/33','0','0'],['0','1/8','0'],['0','0','1/8']]
    assert len(s.symmetry_certificate())==2
    for L in(2,3):assert s.bloch_control(s.graph(L))['complete_source_basis_columns']==12*L**3


def test_complete_coherent_event_and_reference_effect_identity():
    g=s.graph(2);N,D=g.ports,g.events
    def swap(x):
        y=x.copy()
        for e,(a,b)in enumerate(g.seams):y[a,2*e+1],y[b,2*e+1]=x[b,2*e+1],x[a,2*e+1]
        return y
    a=np.array([[((i*7+3*j)%11)-5 for j in range(2)]for i in range(N)],float)
    b=np.array([[((i*3+5*j)%7)-3 for j in range(2)]for i in range(N)],float)
    expand=lambda v:np.repeat(v[:,None,:],D,axis=1)/np.sqrt(D)
    psi=expand(a)+swap(expand(b));actual=swap(2*psi.mean(axis=1,keepdims=True)-psi)
    Tb=b.copy()
    for x,y in g.seams:
        d=(b[x]-b[y])/D;Tb[x]-=d;Tb[y]+=d
    predicted=expand(-b)+swap(expand(a+2*Tb))
    assert np.max(abs(actual-predicted))<2e-14
    # Every terminal source record includes the two reference cross terms.
    for j in range(2):
        for k in range(2):
            product=b[:,j]*b[:,k];Tproduct=product.copy()
            for x,y in g.seams:
                d=(product[x]-product[y])/D;Tproduct[x]-=d;Tproduct[y]+=d
            expected=a[:,j]*a[:,k]+Tproduct+a[:,j]*Tb[:,k]+Tb[:,j]*a[:,k]
            assert np.max(abs((psi[:,:,j]*psi[:,:,k]).sum(axis=1)-expected))<1e-12
    assert np.max(abs(np.einsum('xei,xej->ij',actual,actual)-np.einsum('xei,xej->ij',psi,psi)))<1e-9


def test_packet_support_is_exact_and_axis_moments_agree():
    g=s.graph(4);a=[int(i<12)for i in range(g.ports)];b=[0]*g.ports
    for n in range(4):
        nums=s.walk.port_numerators(g,a,b)
        for x,value in enumerate(nums):
            c=g.cell_coordinates[x//12]
            if sum(min(v,g.size-v)for v in c)>n:assert value==0
        a,b=s.walk.cyclic_step(g,a,b)
    assert all(r['all_three_periodic_axis_second_moments_equal_exact']for r in s.cell_packet(s.graph(2))['records'])


def test_rational_intervals_contain_independent_high_precision_eigenvalues():
    import mpmath as mp
    mp.mp.dps=75
    pi=s.pi_interval();assert mp.mpf(pi[0].numerator)/pi[0].denominator<mp.pi<mp.mpf(pi[1].numerator)/pi[1].denominator
    for angles in (('1/10','0','0'),('1/7','-1/9','1/11'),('1/2','0','0')):
        q=[F(x)for x in angles];cert=s.low_branch_certificate(tuple(s.I(x)for x in q),'independent')
        H=mp.eye(12)*6
        for a,b in s.slow_geometry.source().seams():H[a,b]-=1;H[b,a]-=1
        for axis,row in enumerate(s.links()):
            phase=mp.exp(1j*mp.mpf(q[axis].numerator)/q[axis].denominator)
            for a,b in row:H[a,b]-=phase;H[b,a]-=mp.conj(phase)
        values,_=mp.eighe(H);lo,hi=map(F,cert['ground_eigenvalue_interval_exact'])
        assert mp.mpf(lo.numerator)/lo.denominator<=values[0]<=mp.mpf(hi.numerator)/hi.denominator
        beta=4-max(abs(x)for x in q);assert values[1]>=mp.mpf(beta.numerator)/beta.denominator
        assert F(cert['uniform_internal_preparation_ground_band_weight_lower_exact'])>0
    with pytest.raises(AssertionError):s.low_branch_certificate((s.I(1),s.I(0),s.I(0)),'outside domain')


def test_refinement_uses_declared_clock_and_retains_all_branches():
    rows=s.refinement_certificates();errors=[]
    for row in rows:
        lo,hi=map(F,row['continuum_frequency_squared_interval_exact']);targetlo,targethi=map(F,row['continuum_limit_frequency_squared_interval_exact'])
        assert lo<hi<targetlo and row['full_graph_enumerated']is False
        errors.append(targethi-lo)
    assert errors[2]<errors[1]/3 and errors[1]<errors[0]/3
    # Dropping the cell spacing from the clock would collapse this branch.
    assert float(F(rows[-1]['continuum_frequency_squared_interval_exact'][1]))/64**2<0.003
