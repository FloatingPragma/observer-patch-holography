from fractions import Fraction as F
import numpy as np
from .source import (energy, cycle, phase_records, covariance_components,
                     load_data, inverse_apply, chi_square_interval,
                     sqrt_interval, hubble_squared, decimal, comparison)


def test_exact_pressure_counterfamily_and_capacity_work():
    for k in range(1,7):
        v0,v1=F(1),F(3,2)
        assert -(energy(v1,k)-energy(v0,k))/(v1-v0) == -F(1,k)
        assert energy(v1,k,volume_coupled=False)==energy(v0,k,volume_coupled=False)
        c=cycle(k)
        assert sum(map(F,c['work_strokes']))==0
        assert F(c['spurious_work_if_capacity_strokes_omitted'])>0
        assert sum(F(r['probability']) for r in c['complete_appended_reference_records'])==1
        assert len(c['complete_appended_reference_records'])==2**(k+1)


def test_complete_native_phase_probabilities_and_work():
    # An independent dense controlled-phase circuit, including register records.
    for k in range(1,5):
        h=np.diag([float(energy(F(1),k)),float(energy(F(3,2),k))])
        plus=np.array([1,1])/np.sqrt(2)
        for phase,c in ((F(0),1),(F(1,2),0),(F(1),-1),(F(2),1)):
            r=phase_records(k,phase,c)
            state=np.exp(-1j*np.diag(h)*float(F(r['time_over_pi']))*np.pi)*plus
            for row in r['records']:
                bra=np.array([1,row['X']])/np.sqrt(2)
                assert abs(abs(bra@state)**2/2**k-float(F(row['probability'])))<1e-14
            assert abs(float(F(r['preparation_work']))-(plus@h@plus-h[0,0]))<1e-14
            assert len(r['records'])==2**(k+1)


def test_woodbury_matches_full_matrix_and_offdiagonal_matters():
    data=load_data('HzTable_MM_BC03.dat');d,u=covariance_components(data)
    uf=np.array(u,dtype=float);c=np.diag(np.array(d,dtype=float))+uf.T@uf
    r=[F(i-7,3) for i in range(len(data))]
    exact=np.array(inverse_apply(d,u,r),dtype=float)
    assert np.max(np.abs(exact-np.linalg.solve(c,np.array(r,dtype=float))))<1e-15
    chi=sum(a*b for a,b in zip(r,inverse_apply(d,u,r)))
    diagonal=sum(a*a/b for a,b in zip(r,d))
    assert abs(float(chi-diagonal))>.01


def test_rational_enclosures_cover_perturbed_vectors_and_round_outwards():
    for x in (F(0),F(2),F(9,4),F(1,10**24)):
        a,b=sqrt_interval(x)
        assert a*a<=x<=b*b
    for x in (F(-1,3),F(1,3),F(-1,10**15)):
        assert F(decimal(x))<=x<=F(decimal(x,upper=True))
    data=load_data('HzTable_MM_BC03.dat');d,u=covariance_components(data)
    rs=[(F(i,10)-F(1,100),F(i,10)+F(1,100)) for i in range(len(d))]
    lo,hi=chi_square_interval(d,u,rs)
    for parity in (0,1):
        r=[x[(i+parity)%2] for i,x in enumerate(rs)]
        assert lo<=sum(a*b for a,b in zip(r,inverse_apply(d,u,r)))<=hi


def test_all_rows_and_controls_retained_without_parameter_fitting():
    for name in ('HzTable_MM_BC03.dat','HzTable_MM_M11.dat'):
        r=comparison(name)
        assert r['rows']==15 and r['fitted_chronometer_parameters']==0
        assert len(r['full_covariance_rational'])==15
        assert r['systematics_clamped_at_upper_endpoint_rows']==[14]
        assert F(r['chi_square'][0])<=F(r['chi_square'][1])
    z=F(1)
    assert hubble_squared(z,radiation=F('.0001'))>hubble_squared(z)
